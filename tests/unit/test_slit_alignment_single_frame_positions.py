"""Measuring the trace position of a group that holds one frame, so a nod pair can be combined (DY-1674)."""

from __future__ import annotations

import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit, frame_group
from tests.factories import pipeline_settings
from tests.factories.slit_alignment import synthetic_alignment_map, synthetic_nod_frame
from tests.unit.test_slit_alignment_frames import (
    A_POSITION,
    B_POSITION,
    DETECT_CONTINUUM_SETTINGS,
    TOLERANCE_ARCSEC,
    _use_synthetic_trace_samples,
)

pytestmark = pytest.mark.unit


def _align_pair(log, tmp_path, scene, frameA, frameB, **options):
    return align_frame_groups_along_slit(
        log=log,
        groups=[frame_group([frameA], [frameB], ["A0.fits"]), frame_group([frameB], [frameA], ["B0.fits"])],
        twoDMapPath=scene.path,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
        perOrder=True,
        startNightDate="2024-01-02",
        **options,
    )


@pytest.fixture
def pair(tmp_path, monkeypatch):
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    frameA = synthetic_nod_frame(scene, objectSlit=A_POSITION, mjd=60000.0, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005, offset=-3.0, hasSky=False)
    return scene, frameA, frameB


def test_single_frame_groups_report_the_measured_position_of_each_order_when_asked(log, tmp_path, pair) -> None:
    # ARRANGE
    scene, frameA, frameB = pair

    # ACT
    resultA, resultB = _align_pair(log, tmp_path, scene, frameA, frameB, measureSingleFrames=True)

    # ASSERT: THE POSITIONS ARE MEASURED, THE FRAMES ARE NOT RESAMPLED AND NOTHING COUNTS AS A FALLBACK
    assert resultA.shifts["position_arcsec"].tolist() == pytest.approx([A_POSITION] * 4, abs=TOLERANCE_ARCSEC)
    assert resultB.shifts["position_arcsec"].tolist() == pytest.approx([B_POSITION] * 4, abs=TOLERANCE_ARCSEC)
    assert resultA.shifts["shift_arcsec"].tolist() == [0.0] * 4
    assert resultA.shifts["is_reference"].all()
    assert resultA.frames == (frameA,)
    assert resultB.frames == (frameB,)
    assert resultA.fallbackCount == resultB.fallbackCount == 0
    assert resultA.mapImages is not None


def test_single_frame_groups_are_not_measured_by_default(log, tmp_path, pair) -> None:
    # ARRANGE
    scene, frameA, frameB = pair

    # ACT
    resultA, resultB = _align_pair(log, tmp_path, scene, frameA, frameB)

    # ASSERT
    assert resultA.shifts["position_arcsec"].isna().all()
    assert resultA.mapImages is None
    assert resultB.mapImages is None


def test_an_unreliable_single_frame_trace_is_not_counted_as_a_fallback(log, tmp_path, pair, monkeypatch) -> None:
    # ARRANGE: NO TRACE IS FOUND IN ANY ORDER
    scene, frameA, frameB = pair
    monkeypatch.setattr(detect_continuum, "sample_trace", lambda self: (None, 0.0))

    # ACT
    resultA, _ = _align_pair(log, tmp_path, scene, frameA, frameB, measureSingleFrames=True)

    # ASSERT
    assert resultA.shifts["position_arcsec"].isna().all()
    assert resultA.fallbackCount == 0


def test_a_group_of_several_frames_is_measured_while_a_single_frame_group_is_not_by_default(
    log, tmp_path, monkeypatch
) -> None:
    # ARRANGE: TWO A FRAMES AND ONE B FRAME
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    framesA = [
        synthetic_nod_frame(scene, objectSlit=A_POSITION, mjd=60000.0 + 0.01 * index, hasSky=False)
        for index in range(2)
    ]
    frameB = synthetic_nod_frame(scene, objectSlit=B_POSITION, mjd=60000.005, offset=-3.0, hasSky=False)

    # ACT
    resultA, resultB = align_frame_groups_along_slit(
        log=log,
        groups=[
            frame_group(framesA, [frameB, frameB], ["A0.fits", "A1.fits"]),
            frame_group([frameB], framesA[:1], ["B0.fits"]),
        ],
        twoDMapPath=scene.path,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
        perOrder=True,
        startNightDate="2024-01-02",
    )

    # ASSERT
    assert resultA.shifts["position_arcsec"].notna().all()
    assert resultA.mapImages is not None
    assert resultB.shifts["position_arcsec"].isna().all()
    assert resultB.mapImages is None
