"""Aligning a group of frames along the slit, end to end with the real object-trace detection (DY-1612)."""

from __future__ import annotations

import gc
import weakref

import numpy as np
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit, align_frames_along_slit, frame_group
from tests.factories import pipeline_settings
from tests.factories.slit_alignment import (
    NIR_ORDERS,
    PIXEL_ARCSEC,
    synthetic_alignment_map,
    synthetic_nod_frame,
    zero_slit_samples,
)

pytestmark = pytest.mark.unit

A_POSITION = 2.0
B_POSITION = -2.0
# THE ACCEPTANCE LIMIT: A MEASURED SHIFT IS WITHIN A TENTH OF A PIXEL OF THE TRUTH
TOLERANCE_ARCSEC = 0.1 * PIXEL_ARCSEC

DETECT_CONTINUUM_SETTINGS = {
    "order-sample-count": 150,
    "slice-length": 40,
    "slice-width": 3,
    "peak-sigma-limit": 1.5,
    "disp-axis-deg": 5,
    "order-deg": 3,
    "poly-fitting-residual-clipping-sigma": 4,
    "poly-clipping-iteration-limit": 5,
}


def _use_synthetic_trace_samples(monkeypatch: pytest.MonkeyPatch, scene) -> None:
    """Sample the object trace along the slit-centre pixels of the synthetic map and skip the quicklook plot."""
    monkeypatch.setattr(detect_continuum, "create_pixel_arrays", lambda self: (zero_slit_samples(scene), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)


def _align(log, tmp_path, scene, frames, partners, *, perOrder: bool):
    return align_frames_along_slit(
        log=log,
        frames=frames,
        partnerFrames=partners,
        twoDMapPath=scene.path,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
        perOrder=perOrder,
        frameNames=[f"A{index}.fits" for index in range(len(frames))],
        startNightDate="2024-01-02",
    )


def test_measured_vis_shifts_match_the_known_shift_of_every_order_within_a_tenth_of_a_pixel(
    log, tmp_path, monkeypatch
) -> None:
    # ARRANGE: THREE A FRAMES WHOSE OBJECT MOVES BY A DIFFERENT KNOWN AMOUNT IN EACH ORDER. THE FIRST IS THE REFERENCE
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    moves = [
        {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0},
        {1: 0.3, 2: -0.2, 3: 0.4, 4: 0.1},
        {1: -0.3, 2: 0.2, 3: -0.1, 4: 0.3},
    ]
    frames = [
        synthetic_nod_frame(scene, objectSlit={o: A_POSITION + m for o, m in move.items()}, mjd=60000.0 + 0.01 * index)
        for index, move in enumerate(moves)
    ]
    partners = [
        synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005 + 0.01 * index, offset=-3.0)
        for index in range(3)
    ]

    # ACT
    result = _align(log, tmp_path, scene, frames, partners, perOrder=True)

    # ASSERT
    measured = result.shifts.pivot(index="frame", columns="order", values="shift_arcsec")
    for index, move in enumerate(moves):
        for order, expected in move.items():
            assert measured.loc[f"A{index}.fits", order] == pytest.approx(expected, abs=TOLERANCE_ARCSEC)
    assert result.fallbackCount == 0
    assert len(result.frames) == 3


def _centroid(scene, frame, order: int) -> float:
    """Return the flux-weighted mean slit coordinate of an object-only frame in the middle of an order."""
    middle = (slice(None), slice(200, 400)) if scene.isTransposed else (slice(200, 400), slice(None))
    flux = (scene.orderMap[middle] == order) * frame.data[middle]
    return float((np.nan_to_num(scene.slit[middle]) * flux).sum() / flux.sum())


def test_aligned_frames_put_the_object_of_every_order_where_the_reference_frame_has_it(
    log, tmp_path, monkeypatch
) -> None:
    # ARRANGE
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    moves = [{1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, {1: 0.3, 2: -0.2, 3: 0.4, 4: 0.1}]
    frames = [
        synthetic_nod_frame(
            scene, objectSlit={o: A_POSITION + m for o, m in move.items()}, mjd=60000.0 + 0.01 * index, hasSky=False
        )
        for index, move in enumerate(moves)
    ]
    partners = [
        synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005 + 0.01 * index, offset=-3.0, hasSky=False)
        for index in range(2)
    ]

    # ACT
    result = _align(log, tmp_path, scene, frames, partners, perOrder=True)

    # ASSERT
    for aligned in result.frames:
        for order in scene.orders:
            assert _centroid(scene, aligned, order) == pytest.approx(A_POSITION, abs=TOLERANCE_ARCSEC)
    assert frames[1].header is not result.frames[1].header


def test_an_order_with_no_detectable_trace_falls_back_to_the_median_of_the_other_orders_and_warns(
    log, tmp_path, monkeypatch
) -> None:
    # ARRANGE: THE SECOND FRAME CARRIES NO OBJECT IN ORDER 3. THE OTHER THREE ORDERS MOVE BY 0.3, -0.2 AND 0.1
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    reference = synthetic_nod_frame(scene, objectSlit=A_POSITION, mjd=60000.0, hasSky=False)
    moved = synthetic_nod_frame(
        scene, objectSlit={1: A_POSITION + 0.3, 2: A_POSITION - 0.2, 4: A_POSITION + 0.1}, mjd=60000.01, hasSky=False
    )
    partners = [
        synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005 + 0.01 * index, offset=-3.0, hasSky=False)
        for index in range(2)
    ]

    # ACT
    result = _align(log, tmp_path, scene, [reference, moved], partners, perOrder=True)

    # ASSERT
    fallback = result.shifts.loc[result.shifts["is_fallback"]]
    assert fallback[["frame", "order"]].to_dict("records") == [{"frame": "A1.fits", "order": 3}]
    assert fallback["shift_arcsec"].iloc[0] == pytest.approx(0.1, abs=TOLERANCE_ARCSEC)
    assert result.fallbackCount == 1
    warnings = [message for level, message in log.messages if level == "warning" and "A1.fits" in message]
    assert any("order 3" in message for message in warnings)


def test_nir_frames_get_one_shift_each_within_a_tenth_of_a_pixel(log, tmp_path, monkeypatch) -> None:
    # ARRANGE: NIR HAS THE DISPERSION ALONG THE OTHER AXIS. THE ORDERS OF A FRAME MOVE BY SLIGHTLY DIFFERENT AMOUNTS,
    # AND THE FRAME SHIFT IS THE MEDIAN OF THEM
    scene = synthetic_alignment_map(tmp_path, orders=NIR_ORDERS, isTransposed=True)
    _use_synthetic_trace_samples(monkeypatch, scene)
    reference = synthetic_nod_frame(scene, objectSlit=A_POSITION, mjd=60000.0, arm="NIR", hasSky=False)
    moved = synthetic_nod_frame(
        scene,
        objectSlit={order: A_POSITION + 0.4 + 0.02 * (order - 13) for order in NIR_ORDERS},
        mjd=60000.01,
        arm="NIR",
        hasSky=False,
    )
    partners = [
        synthetic_nod_frame(
            scene, objectSlit=B_POSITION, mjd=60000.005 + 0.01 * index, offset=-3.0, arm="NIR", hasSky=False
        )
        for index in range(2)
    ]

    # ACT
    result = _align(log, tmp_path, scene, [reference, moved], partners, perOrder=False)

    # ASSERT
    assert result.shifts["order"].tolist() == [-1, -1]
    assert result.shifts["shift_arcsec"].tolist() == pytest.approx([0.0, 0.41], abs=TOLERANCE_ARCSEC)
    assert result.fallbackCount == 0


def test_a_single_frame_is_returned_unchanged_with_a_zero_shift(log, tmp_path) -> None:
    # ARRANGE
    scene = synthetic_alignment_map(tmp_path)
    frame = synthetic_nod_frame(scene, objectSlit=A_POSITION, mjd=60000.0)
    partner = synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005, offset=-3.0)

    # ACT
    result = _align(log, tmp_path, scene, [frame], [partner], perOrder=True)

    # ASSERT
    assert result.frames == (frame,)
    assert result.shifts["shift_arcsec"].tolist() == [0.0]
    assert result.fallbackCount == 0


def test_two_groups_are_measured_before_either_is_resampled_and_originals_are_released(
    log, tmp_path, monkeypatch
) -> None:
    # ARRANGE: AN A GROUP AND A B GROUP THAT ARE EACH OTHER'S PARTNERS. THE SECOND FRAME OF EACH MOVES BY A KNOWN AMOUNT
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    framesA = [
        synthetic_nod_frame(scene, objectSlit=A_POSITION + move, mjd=60000.0 + 0.01 * index, hasSky=False)
        for index, move in enumerate((0.0, 0.3))
    ]
    framesB = [
        synthetic_nod_frame(
            scene, objectSlit=B_POSITION + move, mjd=60000.005 + 0.01 * index, offset=-3.0, hasSky=False
        )
        for index, move in enumerate((0.0, -0.2))
    ]
    originals = [weakref.ref(frame) for frame in framesA + framesB]
    groups = [
        frame_group(framesA, framesB, ["A0.fits", "A1.fits"]),
        frame_group(framesB, framesA, ["B0.fits", "B1.fits"]),
    ]

    # ACT
    resultA, resultB = align_frame_groups_along_slit(
        log=log,
        groups=groups,
        twoDMapPath=scene.path,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
        perOrder=True,
        startNightDate="2024-01-02",
    )
    gc.collect()

    # ASSERT: THE PARTNERS WERE UNSHIFTED WHEN SUBTRACTED, SO BOTH GROUPS' SHIFTS ARE THE KNOWN ONES
    assert resultA.shifts["shift_arcsec"].tolist() == pytest.approx([0.0] * 4 + [0.3] * 4, abs=TOLERANCE_ARCSEC)
    assert resultB.shifts["shift_arcsec"].tolist() == pytest.approx([0.0] * 4 + [-0.2] * 4, abs=TOLERANCE_ARCSEC)
    # THE LISTS HOLD THE ALIGNED FRAMES, AND NOTHING KEEPS AN ORIGINAL ALIVE
    assert framesA == list(resultA.frames)
    assert framesB == list(resultB.frames)
    assert [reference() for reference in originals] == [None] * 4
