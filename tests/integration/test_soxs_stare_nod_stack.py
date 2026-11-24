"""The stare path of a nodding SOF aligns and stacks the A and B frames each on their own (DY-2194, DY-2290).

The frames of one nod position are aligned to each other along the slit. The B frames are not shifted onto the A frames,
so the object stays where it is in each stack. The slit alignment and the stacks run for real on synthetic frames. The
detrending and the extraction tail are replaced by recorders.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.slit_alignment import WHOLE_FRAME, build_two_d_map_images, slit_alignment_result
from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import position_stack, reduction_inputs
from tests.factories import pipeline_settings, qc_table
from tests.factories.slit_alignment import FOOTPRINT_HALF, PIXEL_ARCSEC, synthetic_alignment_map, synthetic_nod_frame
from tests.integration.test_observing_recipe_orchestration import _empty_products
from tests.unit.test_slit_alignment_frames import DETECT_CONTINUUM_SETTINGS, _use_synthetic_trace_samples

pytestmark = pytest.mark.integration

# THE SLIT POSITION OF THE OBJECT IS MINUS THE CUMULATIVE OFFSET, AS MEASURED ON REAL NOD FRAMES
OFFSET_A = 2.1
OFFSET_B = -2.9
SLIT_A = -OFFSET_A
SLIT_B = -OFFSET_B
MIDDLE_ROWS = slice(200, 400)
SOF_NAME = "20251212T045200_VIS_1X1_1_NOD_OBJ_SOXS"
# A SHIFT LONGER THAN THE SYNTHETIC SLIT, AND THE DISTANCE FROM A COVERAGE EDGE AT WHICH A PIXEL IS SAFELY ON ONE SIDE
BEYOND_SLIT_LENGTH = 3 * FOOTPRINT_HALF
EDGE_MARGIN_ARCSEC = 0.4
# HOW FAR THE SECOND FRAME OF A NOD POSITION IS FROM THE FIRST, ALONG THE SLIT
JITTER_ARCSEC = 1.0
# THE LARGEST SHIFT THE ALIGNMENT OF A FRAME TO ITS OWN POSITION NEEDS, AND THE NOD THROW IT MUST NOT APPLY
ALIGNMENT_TOLERANCE_ARCSEC = 0.15
NOD_THROW = OFFSET_A - OFFSET_B


def _frames(
    scene: Any, *, countPerNod: int = 2, offsetsA: tuple[float, ...] | None = None, jitterB: float = 0.0
) -> list[Any]:
    """Return A and B frames in time order: A, B, A, B, ... The second B frame's object is ``jitterB`` arcsec higher."""
    frames = []
    for index in range(countPerNod):
        frames.append(synthetic_nod_frame(scene, objectSlit=SLIT_A, mjd=60000.0 + 0.02 * index, offset=OFFSET_A))
        frames.append(
            synthetic_nod_frame(
                scene,
                objectSlit=SLIT_B + (jitterB if index == 1 else 0.0),
                mjd=60000.01 + 0.02 * index,
                offset=OFFSET_B,
            )
        )
    return frames


def _recipe(log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, scene: Any) -> soxs_stare:
    """Return a stare recipe in nod mode, with its calibration inputs, detrending and stacking recorded."""
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.settings = pipeline_settings(tmpPath)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.imageType = "OBJECT"
    recipe.reduceNodAsStare = True
    recipe.verbose = False
    recipe.debug = False
    recipe.recipeName = "soxs-stare"
    recipe.recipeSettings = {
        "use_flat": True,
        "detect-continuum": DETECT_CONTINUUM_SETTINGS,
        "stacked-clipping-sigma": 15,
        "stacked-clipping-iterations": 3,
    }
    recipe.sofName = SOF_NAME
    recipe.productNameMarker = "_STAREMODE"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.products = _empty_products()
    recipe.dateObs = "2024-01-02T03:04:05.678"
    recipe.objectFrameNames = [f"frame-{index}.fits" for index in range(4)]
    _use_synthetic_trace_samples(monkeypatch, scene)
    return recipe


def _inputs(scene: Any, *, masterFlat: Any = "FLAT") -> tuple:
    """The calibration inputs: a bias, no dark, a flat or none, and the order, dispersion and 2D map paths."""
    return reduction_inputs("BIAS", False, masterFlat, "order.fits", "disp.fits", scene.path, False)


def _stack_positions(
    recipe: soxs_stare, scene: Any, frames: list[Any], *, masterFlat: Any = False
) -> tuple[position_stack, position_stack]:
    return recipe._stack_nod_positions(*recipe._split_nod_frames(frames), _inputs(scene, masterFlat=masterFlat))


def _record_detrend(recipe: Any, monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []

    def detrend(**kwargs: Any) -> Any:
        calls.append(kwargs)
        return kwargs["inputFrame"]

    monkeypatch.setattr(recipe, "detrend", detrend)
    return calls


def _record_stack(recipe: Any, monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []
    realStack = recipe.clip_and_stack

    def stack(**kwargs: Any) -> Any:
        result = realStack(**kwargs)
        # THE RECIPE EMPTIES THE LIST AFTER THE STACK, SO KEEP A COPY OF THE FRAMES IT WAS GIVEN
        calls.append({**kwargs, "frames": list(kwargs["frames"]), "result": result})
        return result

    monkeypatch.setattr(recipe, "clip_and_stack", stack)
    return calls


def _record_shifts(monkeypatch: pytest.MonkeyPatch) -> list[dict[int, float]]:
    """Record the order shifts of every frame the slit alignment resamples, and resample it as usual."""
    slitModule = import_module("soxspipe.commonutils.slit_alignment")
    realShift = slitModule.shift_frame_along_slit
    calls: list[dict[int, float]] = []

    def shift(frame: Any, orderShifts: dict[int, float], mapImages: Any) -> Any:
        calls.append(dict(orderShifts))
        return realShift(frame, orderShifts, mapImages)

    monkeypatch.setattr(slitModule, "shift_frame_along_slit", shift)
    return calls


def _peak(scene: Any, frame: Any, order: int) -> tuple[float, float]:
    """Return the slit coordinate and the counts of the brightest median column profile of an order, middle rows."""
    isOrder = (scene.orderMap == order)[MIDDLE_ROWS]
    columns = np.flatnonzero(isOrder.all(axis=0))
    slit = np.median(np.where(isOrder, scene.slit[MIDDLE_ROWS], np.nan)[:, columns], axis=0)
    data = np.where(frame.mask[MIDDLE_ROWS], np.nan, frame.data[MIDDLE_ROWS]) if frame.mask is not None else frame.data
    counts = np.nanmedian(np.where(isOrder, data, np.nan)[:, columns], axis=0)
    return float(slit[np.nanargmax(counts)]), float(np.nanmax(counts))


def _peak_slit(scene: Any, frame: Any, order: int) -> float:
    return _peak(scene, frame, order)[0]


def _middle(scene: Any, isPixel: np.ndarray) -> np.ndarray:
    """Return ``isPixel`` limited to the pixels of the orders in the middle rows."""
    limited = np.zeros(scene.shape, dtype=bool)
    limited[MIDDLE_ROWS] = (isPixel & np.isfinite(scene.orderMap))[MIDDLE_ROWS]
    return limited


@pytest.fixture
def scene(tmp_path: Path) -> Any:
    return synthetic_alignment_map(tmp_path)


def test_each_frame_is_detrended_without_the_flat_then_the_aligned_frames_are_flat_fielded(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    detrends = _record_detrend(recipe, monkeypatch)
    frames = _frames(scene)

    # ACT
    _stack_positions(recipe, scene, frames, masterFlat="FLAT")

    # ASSERT: FOUR FRAMES WITH THE BIAS AND THE SCATTERED-LIGHT TABLE AND NO FLAT, THEN FOUR WITH THE FLAT ONLY
    assert len(detrends) == 8
    for call in detrends[:4]:
        assert (call["master_bias"], call["dark"], call["master_flat"], call["order_table"]) == (
            "BIAS",
            False,
            False,
            "order.fits",
        )
    for call in detrends[4:]:
        assert (call["master_bias"], call["dark"], call["master_flat"], call["order_table"]) == (
            False,
            False,
            "FLAT",
            False,
        )


def test_the_dark_is_what_the_nir_frames_are_detrended_with(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    detrends = _record_detrend(recipe, monkeypatch)

    # ACT
    recipe._stack_nod_positions(
        *recipe._split_nod_frames(_frames(scene)),
        reduction_inputs(False, "DARK", False, "order.fits", "disp.fits", scene.path, False),
    )

    # ASSERT: NO FLAT, SO THE FOUR DETRENDS ARE THE ONLY ONES
    assert [(call["master_bias"], call["dark"]) for call in detrends] == [(False, "DARK")] * 4


def test_each_nod_position_stack_records_the_file_names_of_its_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT: THE FRAMES ARE READ A, B, A, B
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene))

    # ASSERT
    assert stackA.frameNames == ("frame-0.fits", "frame-2.fits")
    assert stackB.frameNames == ("frame-1.fits", "frame-3.fits")


def test_each_nod_position_is_stacked_on_its_own_and_the_b_object_stays_where_it_is(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)
    frames = _frames(scene)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, frames)

    # ASSERT: TWO STACKS OF TWO FRAMES EACH, STACKED AS STARE FRAMES THAT KEEP THEIR MASKS
    assert [call["recipe"] for call in stacks] == ["soxs_stare", "soxs_stare"]
    assert [call["ignore_input_masks"] for call in stacks] == [False, False]
    assert [call["post_stack_clipping"] for call in stacks] == [False, False]
    assert [len(call["frames"]) for call in stacks] == [2, 2]
    assert (stackA.position, stackB.position) == ("A", "B")
    # THE B OBJECT IS NOT MOVED, AND THE B STACK KEEPS THE FULL HEIGHT OF ONE FRAME
    for order in scene.orders:
        slitA, heightA = _peak(scene, stackA.notFlattened, order)
        slitB, heightB = _peak(scene, stackB.notFlattened, order)
        assert slitA == pytest.approx(SLIT_A, abs=2 * PIXEL_ARCSEC)
        assert slitB == pytest.approx(SLIT_B, abs=2 * PIXEL_ARCSEC)
        assert heightA == pytest.approx(_peak(scene, frames[0], order)[1], rel=0.1)
        assert heightB == pytest.approx(_peak(scene, frames[1], order)[1], rel=0.1)


def test_each_position_keeps_the_raw_frames_of_its_own_nod_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    frames = _frames(scene)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, frames)

    # ASSERT
    assert stackA.objectFrames == [frames[0], frames[2]]
    assert stackB.objectFrames == [frames[1], frames[3]]


def test_the_frames_are_resampled_only_by_their_alignment_to_their_own_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: THE SECOND B FRAME SITS 1 ARCSEC HIGHER ALONG THE SLIT
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    shifts = _record_shifts(monkeypatch)

    # ACT
    _stack_positions(recipe, scene, _frames(scene, jitterB=JITTER_ARCSEC))

    # ASSERT: EACH FRAME OF AN ALIGNED POSITION IS RESAMPLED ONCE. THE LARGEST SHIFT IS THE JITTER, NOT THE NOD THROW
    assert len(shifts) == 4
    largest = [max(abs(shift) for shift in frameShifts.values()) for frameShifts in shifts]
    assert largest[0] == pytest.approx(0.0, abs=ALIGNMENT_TOLERANCE_ARCSEC)
    assert largest[1] == pytest.approx(0.0, abs=ALIGNMENT_TOLERANCE_ARCSEC)
    assert largest[2] == pytest.approx(0.0, abs=ALIGNMENT_TOLERANCE_ARCSEC)
    assert largest[3] == pytest.approx(JITTER_ARCSEC, abs=ALIGNMENT_TOLERANCE_ARCSEC)
    assert max(largest) < NOD_THROW / 2


def test_a_position_with_no_jitter_has_nothing_uncovered(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene))

    # ASSERT: THE ALIGNMENT SHIFTS ARE TOO SMALL TO UNCOVER A PIXEL AWAY FROM THE SLIT ENDS
    inside = _middle(scene, np.abs(scene.slit) < FOOTPRINT_HALF - EDGE_MARGIN_ARCSEC)
    assert stackA.uncovered.dtype == bool
    assert not stackA.uncovered[inside].any()
    assert not stackB.uncovered[inside].any()


def test_the_jitter_of_one_position_leaves_the_slit_end_of_that_position_only_uncovered(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: THE SECOND B FRAME IS ALIGNED BY +1 ARCSEC, SO ITS TOP ARCSEC TAKES ITS VALUE FROM OFF THE SLIT
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    edge = FOOTPRINT_HALF - JITTER_ARCSEC
    isUncoveredBand = _middle(scene, scene.slit > edge + EDGE_MARGIN_ARCSEC)
    isCoveredBand = _middle(scene, scene.slit < edge - EDGE_MARGIN_ARCSEC)
    assert isUncoveredBand.any() and isCoveredBand.any()

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene, jitterB=JITTER_ARCSEC))

    # ASSERT: THE A FRAMES COVER THE WHOLE SLIT, AND THE B FRAMES DO NOT
    assert stackB.uncovered[isUncoveredBand].all()
    assert not stackB.uncovered[isCoveredBand].any()
    assert not stackA.uncovered[isUncoveredBand].any()
    assert not stackB.uncovered[~np.isfinite(scene.orderMap)].any()


def test_a_position_of_one_frame_needs_no_alignment_and_has_nothing_uncovered(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    recipe.objectFrameNames = ["a.fits", "b.fits"]
    _record_detrend(recipe, monkeypatch)
    shifts = _record_shifts(monkeypatch)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene, countPerNod=1))

    # ASSERT
    assert shifts == []
    assert not stackA.uncovered.any()
    assert not stackB.uncovered.any()
    assert stackA.uncovered.shape == scene.shape


def _order_shift_table(names: list[str], shiftsOfFrames: list[dict[int, float]]) -> pd.DataFrame:
    """Return the per-order shift table of a group whose frames were aligned by the given order shifts."""
    return pd.DataFrame(
        [
            {
                "frame": name,
                "mjd": float(index),
                "order": order,
                "position_arcsec": np.nan,
                "shift_arcsec": shift,
                "is_fallback": False,
                "is_reference": index == 0,
            }
            for index, (name, shifts) in enumerate(zip(names, shiftsOfFrames, strict=True))
            for order, shift in shifts.items()
        ]
    )


def _stub_alignment_with_shifts(
    monkeypatch: pytest.MonkeyPatch,
    aligned: tuple[list[Any], list[Any]],
    tables: tuple[pd.DataFrame, pd.DataFrame],
    mapImages: Any,
) -> None:
    """Replace the slit alignment with one that returns the given frames and shift tables."""
    results = [
        slit_alignment_result(tuple(frames), table, 0, mapImages) for frames, table in zip(aligned, tables, strict=True)
    ]
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.slit_alignment"), "align_frame_groups_along_slit", lambda **_: results
    )


def _stub_order_shifts(
    monkeypatch: pytest.MonkeyPatch, recipe: soxs_stare, scene: Any, shiftsB: list[dict[int, float]]
) -> None:
    """Make the slit alignment report no shift for the A frames and ``shiftsB`` for the B frames."""
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene))
    noShift = dict.fromkeys(scene.orders, 0.0)
    _stub_alignment_with_shifts(
        monkeypatch,
        (framesA, framesB),
        (_order_shift_table(namesA, [noShift, noShift]), _order_shift_table(namesB, shiftsB)),
        build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap),
    )


def test_an_order_with_no_common_slit_range_in_a_position_is_warned_about_and_left_out_whole(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: THE SECOND B FRAME IS ALIGNED BY MORE THAN THE SLIT LENGTH IN ORDER 1 ONLY
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    noShift = dict.fromkeys(scene.orders, 0.0)
    _stub_order_shifts(monkeypatch, recipe, scene, [noShift, {**noShift, 1: BEYOND_SLIT_LENGTH}])

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene))

    # ASSERT
    assert stackB.uncovered[scene.orderMap == 1].all()
    assert not stackB.uncovered[scene.orderMap == 2].all()
    assert not stackA.uncovered.any()
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("no slit range" in message and "order 1" in message and "B frames" in message for message in warnings)
    assert not any("order 2" in message for message in warnings)


def test_a_position_whose_alignment_leaves_no_common_slit_range_in_any_order_is_refused(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    noShift = dict.fromkeys(scene.orders, 0.0)
    _stub_order_shifts(monkeypatch, recipe, scene, [noShift, dict.fromkeys(scene.orders, BEYOND_SLIT_LENGTH)])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no slit range that every B frame covers"):
        _stack_positions(recipe, scene, _frames(scene))


def test_the_unflattened_stack_is_the_same_stack_when_no_flat_is_used(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene), masterFlat=False)

    # ASSERT: ONE STACK FOR EACH POSITION
    assert stackA.notFlattened is stackA.flattened
    assert stackB.notFlattened is stackB.flattened
    assert len(stacks) == 2


def test_with_a_flat_both_the_unflattened_and_the_flat_fielded_frames_of_each_position_are_stacked(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    detrends = _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene), masterFlat="FLAT")

    # ASSERT: A'S NOT-FLATTENED STACK, A'S FLAT-FIELDED STACK, THEN B'S TWO. A FLAT-FIELD IS A FLAT-ONLY DETREND
    assert len(stacks) == 4
    assert stackA.notFlattened is stacks[0]["result"]
    assert stackA.flattened is stacks[1]["result"]
    assert stackB.notFlattened is stacks[2]["result"]
    assert stackB.flattened is stacks[3]["result"]
    assert [len(call["frames"]) for call in stacks] == [2, 2, 2, 2]
    assert [call["master_flat"] for call in detrends[4:]] == ["FLAT"] * 4


def test_each_stack_header_comes_from_the_frames_of_its_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    offsetKey = recipe.kw("NOD_CUMULATIVE_OFFSETY")

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene))

    # ASSERT
    assert stackA.notFlattened.header[offsetKey] == OFFSET_A
    assert stackB.notFlattened.header[offsetKey] == OFFSET_B


def test_the_aligned_frames_come_from_the_alignment_results_not_from_a_change_to_the_input_lists(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene))
    alignedA, alignedB = list(reversed(framesA)), list(reversed(framesB))
    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    _stub_alignment(monkeypatch, (alignedA, alignedB), (namesA, namesB), mapImages, ([0.0, 0.7], [0.0, -0.4]))

    # ACT
    returnedA, returnedB, alignmentShifts, returnedMap = recipe._align_nod_frames(
        framesA, framesB, namesA, namesB, "map.fits", "disp.fits"
    )

    # ASSERT: THE ALIGNMENT SHIFT OF EVERY ORDER OF EACH FRAME COMES BACK FOR THE COVERAGE
    assert returnedA == alignedA
    assert returnedB == alignedB
    assert returnedMap is mapImages
    shiftsA, shiftsB = alignmentShifts
    assert [set(frame.values()) for frame in shiftsA] == [{0.0}, {0.7}]
    assert [set(frame.values()) for frame in shiftsB] == [{0.0}, {-0.4}]
    assert set(shiftsA[0]) == set(scene.orders)


def test_positions_the_alignment_read_no_map_for_have_nothing_uncovered(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: NO POSITION NEEDED ALIGNING, SO THE ALIGNMENT READ NO 2D MAP
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    recipe.objectFrameNames = ["a.fits", "b.fits"]
    _record_detrend(recipe, monkeypatch)
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene, countPerNod=1))
    _stub_alignment(monkeypatch, (framesA, framesB), (namesA, namesB), None)

    # ACT
    stackA, stackB = _stack_positions(recipe, scene, _frames(scene, countPerNod=1))

    # ASSERT
    assert not stackA.uncovered.any()
    assert stackB.uncovered.shape == scene.shape


def test_unbalanced_nodding_frames_are_refused_before_any_frame_is_detrended(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    detrends = _record_detrend(recipe, monkeypatch)

    # ACT / ASSERT
    recipe.objectFrameNames = recipe.objectFrameNames[:3]
    with pytest.raises(Exception, match="Found 2 A frames and 1 B frames"):
        recipe._split_nod_frames(_frames(scene)[:3])

    assert detrends == []
    assert any(level == "error" for level, _ in log.messages)


def test_the_frame_names_come_from_the_raw_frame_names_of_the_run(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    recipe.objectFrameNames = ["a1.fits", "b1.fits", "a2.fits", "b2.fits"]

    # ACT
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene))

    # ASSERT
    assert (namesA, namesB) == (["a1.fits", "a2.fits"], ["b1.fits", "b2.fits"])
    assert len(framesA) == len(framesB) == 2


def test_the_stacked_positions_reach_the_stare_tail_without_a_second_detrend(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    tailCalls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(recipe, "_read_reduction_inputs", lambda kw, arm: _inputs(scene))
    monkeypatch.setattr(recipe, "_stack_nod_positions", lambda *args: ("STACK A", "STACK B"))
    monkeypatch.setattr(
        recipe, "_reduce_nod_positions", lambda *args: tailCalls.append(args) or ("product.fits", False)
    )
    frames = _frames(scene)

    # ACT
    result = recipe._reduce_object_frames(frames, rawFrames=["raw.fits"])

    # ASSERT
    assert result == ("product.fits", False)
    assert tailCalls == [(["STACK A", "STACK B"], frames, ["raw.fits"], _inputs(scene))]


def test_a_single_frame_is_reduced_without_the_nod_stack(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    tailCalls: list[tuple[Any, ...]] = []
    stackCalls: list[dict[str, Any]] = []
    monkeypatch.setattr(recipe, "_stack_nod_positions", lambda *args: pytest.fail("a single frame needs no alignment"))
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **kwargs: stackCalls.append(kwargs) or "STACK")
    monkeypatch.setattr(
        recipe, "_reduce_stacked_object", lambda *args: tailCalls.append(args) or ("product.fits", False)
    )
    frame = _frames(scene)[0]

    # ACT
    recipe._reduce_object_frames([frame], rawFrames=["a.fits"])

    # ASSERT
    assert len(stackCalls) == 1
    assert tailCalls == [("STACK", [frame], ["a.fits"])]


def test_a_stare_recipe_outside_nod_mode_stacks_its_frames_as_before(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    recipe.reduceNodAsStare = False
    stackCalls: list[dict[str, Any]] = []
    tailCalls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(recipe, "_stack_nod_positions", lambda *args: pytest.fail("not in nod mode"))
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **kwargs: stackCalls.append(kwargs) or "STACK")
    monkeypatch.setattr(
        recipe, "_reduce_stacked_object", lambda *args: tailCalls.append(args) or ("product.fits", False)
    )
    frames = _frames(scene)

    # ACT
    recipe._reduce_object_frames(frames)

    # ASSERT
    assert stackCalls == [
        {"frames": frames, "recipe": "soxs_stare", "ignore_input_masks": True, "post_stack_clipping": False}
    ]
    assert tailCalls == [("STACK", frames, False)]


def test_a_repeated_background_product_row_is_kept_once(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: EACH FRAME'S SCATTERED-LIGHT FIT ADDS A ROW FOR THE ONE BACKGROUND PLOT THE LAST FRAME LEAVES
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    row = dict.fromkeys(recipe.products.columns, "x")

    def detrend(**kwargs: Any) -> Any:
        recipe.products = pd.concat([recipe.products, pd.DataFrame([row])], ignore_index=True)
        return kwargs["inputFrame"]

    monkeypatch.setattr(recipe, "detrend", detrend)

    # ACT
    _stack_positions(recipe, scene, _frames(scene))

    # ASSERT
    assert len(recipe.products) == 1


def test_frames_and_names_of_different_numbers_are_refused(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    recipe.objectFrameNames = ["a1.fits", "b1.fits"]

    # ACT / ASSERT
    with pytest.raises(ValueError, match="4 object frames but 2 raw file names"):
        recipe._split_nod_frames(_frames(scene))


def _shift_table(names: list[str], shifts: list[float]) -> pd.DataFrame:
    """Return the whole-frame shift table of a group whose frames were aligned by ``shifts`` arcsec."""
    return pd.DataFrame(
        {
            "frame": names,
            "mjd": np.arange(len(names), dtype=float),
            "order": WHOLE_FRAME,
            "position_arcsec": np.nan,
            "shift_arcsec": shifts,
            "is_fallback": False,
            "is_reference": [index == 0 for index in range(len(names))],
        }
    )


def _stub_alignment(
    monkeypatch: pytest.MonkeyPatch,
    aligned: tuple[list[Any], list[Any]],
    names: tuple[list[str], list[str]],
    mapImages: Any,
    shifts: tuple[list[float], list[float]] | None = None,
) -> None:
    """Replace the slit alignment with one that leaves the caller's lists alone and returns other frames."""
    shifts = shifts or ([0.0] * len(names[0]), [0.0] * len(names[1]))
    results = [
        slit_alignment_result(tuple(frames), _shift_table(frameNames, frameShifts), 0, mapImages)
        for frames, frameNames, frameShifts in zip(aligned, names, shifts, strict=True)
    ]
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.slit_alignment"), "align_frame_groups_along_slit", lambda **_: results
    )
