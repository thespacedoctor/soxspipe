"""The stacked stare-mode path of a nodding SOF aligns the A and B frames along the slit before stacking (DY-2194).

The slit alignment, the nominal-throw fallback and the stack run for real on synthetic frames. The detrending and the
extraction tail are replaced by recorders.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.slit_alignment import build_two_d_map_images
from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import reduction_inputs
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


def _frames(scene: Any, *, countPerNod: int = 2, offsetsA: tuple[float, ...] | None = None) -> list[Any]:
    """Return A and B frames in time order: A, B, A, B, ..."""
    frames = []
    for index in range(countPerNod):
        frames.append(synthetic_nod_frame(scene, objectSlit=SLIT_A, mjd=60000.0 + 0.02 * index, offset=OFFSET_A))
        frames.append(synthetic_nod_frame(scene, objectSlit=SLIT_B, mjd=60000.01 + 0.02 * index, offset=OFFSET_B))
    return frames


def _recipe(log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, scene: Any):
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
        calls.append({**kwargs, "result": result})
        return result

    monkeypatch.setattr(recipe, "clip_and_stack", stack)
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


def _patch_plan(monkeypatch: pytest.MonkeyPatch, scene: Any, orderShifts: dict[int, float]) -> None:
    """Make the slit alignment return a plan that shifts the B frames by ``orderShifts``."""
    from soxspipe.commonutils.nod_shift_and_add import shift_and_add_plan

    mapImages = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.nod_shift_and_add"),
        "build_shift_and_add_plan",
        lambda *_: shift_and_add_plan(orderShifts, mapImages),
    )


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


def test_each_frame_is_detrended_without_the_flat_then_the_aligned_frames_are_flat_fielded(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    detrends = _record_detrend(recipe, monkeypatch)
    frames = _frames(scene)

    # ACT
    recipe._stack_nod_frames(*recipe._split_nod_frames(frames), _inputs(scene))

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
    recipe._stack_nod_frames(
        *recipe._split_nod_frames(_frames(scene)),
        reduction_inputs(False, "DARK", False, "order.fits", "disp.fits", scene.path, False),
    )

    # ASSERT: NO FLAT, SO THE FOUR DETRENDS ARE THE ONLY ONES
    assert [(call["master_bias"], call["dark"]) for call in detrends] == [(False, "DARK")] * 4


def test_the_stack_puts_the_b_object_on_the_a_object_and_honours_the_masks(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)

    # ACT
    stacked, _, _ = recipe._stack_nod_frames(
        *recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False)
    )

    # ASSERT
    assert len(stacks) == 1
    assert stacks[0]["recipe"] == "soxs_stare"
    assert stacks[0]["ignore_input_masks"] is False
    assert stacks[0]["post_stack_clipping"] is False
    assert len(stacks[0]["frames"]) == 4
    for order in scene.orders:
        assert _peak_slit(scene, stacked, order) == pytest.approx(SLIT_A, abs=2 * PIXEL_ARCSEC)

    # THE B FRAMES ARE MASKED WHERE THEIR SOURCE FELL OFF THE SLIT, AND THE A FRAMES ARE NOT
    isOffSlit = (scene.orderMap == scene.orders[0]) & (scene.slit + (OFFSET_A - OFFSET_B) > FOOTPRINT_HALF + 0.4)
    isOffSlit[:200] = isOffSlit[400:] = False
    stackedFrames = stacks[0]["frames"]
    assert [bool(frame.mask[isOffSlit].all()) for frame in stackedFrames] == [False, True, False, True]
    # THE STACK KEEPS THOSE PIXELS, AS THE A FRAMES COVER THEM
    assert not stacked.mask[isOffSlit].any()


def test_the_stack_reports_the_slit_range_the_shifted_b_frames_do_not_cover(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT
    _, _, uncovered = recipe._stack_nod_frames(
        *recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False)
    )

    # ASSERT: THE B FRAMES ARE SHIFTED BY 5.0 ARCSEC, SO THE SLIT BEYOND FOOTPRINT_HALF - 5.0 HOLDS ONLY A FRAMES
    edge = FOOTPRINT_HALF - (OFFSET_A - OFFSET_B)
    isInside = np.isfinite(scene.orderMap)
    middle = np.zeros(scene.shape, dtype=bool)
    middle[MIDDLE_ROWS] = True
    assert uncovered.dtype == bool
    assert uncovered[isInside & middle & (scene.slit > edge + 0.4)].all()
    assert not uncovered[isInside & middle & (scene.slit < edge - 0.4)].any()
    assert not uncovered[~isInside].any()


def test_an_order_with_no_common_slit_range_is_warned_about_and_left_out_whole(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: THE B FRAMES ARE SHIFTED BY MORE THAN THE SLIT LENGTH IN ORDER 1 ONLY
    _patch_plan(monkeypatch, scene, {1: 3 * FOOTPRINT_HALF, 2: 4.6, 3: 4.6, 4: 4.6})
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT
    _, _, uncovered = recipe._stack_nod_frames(
        *recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False)
    )

    # ASSERT
    assert uncovered[scene.orderMap == 1].all()
    assert not uncovered[scene.orderMap == 2].all()
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("no slit range" in message and "order 1" in message for message in warnings)
    assert not any("order 2" in message for message in warnings)


def test_a_shift_that_leaves_no_common_slit_range_in_any_order_is_refused(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    _patch_plan(monkeypatch, scene, dict.fromkeys(scene.orders, 3 * FOOTPRINT_HALF))
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no slit range that every frame covers"):
        recipe._stack_nod_frames(*recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False))


def test_the_unflattened_stack_is_the_same_stack_when_no_flat_is_used(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)

    # ACT
    notFlattened, flattened, _ = recipe._stack_nod_frames(
        *recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False)
    )

    # ASSERT
    assert notFlattened is flattened
    assert len(stacks) == 1


def test_with_a_flat_both_the_unflattened_and_the_flat_fielded_frames_are_stacked(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    stacks = _record_stack(recipe, monkeypatch)

    # ACT
    notFlattened, flattened, _ = recipe._stack_nod_frames(*recipe._split_nod_frames(_frames(scene)), _inputs(scene))

    # ASSERT
    assert len(stacks) == 2
    assert notFlattened is stacks[0]["result"]
    assert flattened is stacks[1]["result"]


def test_a_measured_position_is_used_and_recorded_as_not_nominal(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)

    # ACT
    recipe._stack_nod_frames(*recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False))

    # ASSERT
    row = recipe.qc.set_index("qc_name").loc["NOD STARE NOMINAL THROW"]
    assert row["qc_value"] == 0
    assert not [message for level, message in log.messages if level == "warning" and "nominal" in message]


def test_without_a_measured_position_the_nominal_throw_shifts_the_b_frames_and_is_recorded(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE: THE OBJECT POSITION CANNOT BE MEASURED AT BOTH NOD POSITIONS
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.nod_shift_and_add"), "build_shift_and_add_plan", lambda *_: None
    )
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    _record_detrend(recipe, monkeypatch)
    frames = _frames(scene)

    # ACT
    stacked, _, _ = recipe._stack_nod_frames(*recipe._split_nod_frames(frames), _inputs(scene, masterFlat=False))

    # ASSERT: THE B OBJECT LANDS ON THE A OBJECT, SO THE STACK KEEPS THE FULL HEIGHT OF ONE FRAME (HALF IF IT DID NOT).
    # THE WARNING SAYS THE NOMINAL THROW WAS USED, AND THE QC SAYS SO
    for order in scene.orders:
        slit, height = _peak(scene, stacked, order)
        assert slit == pytest.approx(SLIT_A, abs=2 * PIXEL_ARCSEC)
        assert height == pytest.approx(_peak(scene, frames[0], order)[1], rel=0.1)
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("nominal nod throw" in message and "+5.000" in message for message in warnings)
    row = recipe.qc.set_index("qc_name").loc["NOD STARE NOMINAL THROW"]
    assert row["qc_value"] == 1


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


def test_the_stacked_frames_reach_the_stare_tail_without_a_second_detrend(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    tailCalls: list[tuple[Any, ...]] = []
    monkeypatch.setattr(recipe, "_read_reduction_inputs", lambda kw, arm: _inputs(scene))
    monkeypatch.setattr(recipe, "_stack_nod_frames", lambda *args: ("NOTFLAT", "FLAT", "UNCOVERED"))
    monkeypatch.setattr(
        recipe, "_reduce_stacked_object", lambda *args: tailCalls.append(args) or ("product.fits", False)
    )
    frames = _frames(scene)

    # ACT
    result = recipe._reduce_object_frames(frames, rawFrames=["raw.fits"])

    # ASSERT
    assert result == ("product.fits", False)
    assert tailCalls == [("NOTFLAT", frames, ["raw.fits"], _inputs(scene), "FLAT", "UNCOVERED")]


def test_a_single_frame_is_reduced_without_the_nod_stack(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    tailCalls: list[tuple[Any, ...]] = []
    stackCalls: list[dict[str, Any]] = []
    monkeypatch.setattr(recipe, "_stack_nod_frames", lambda *args: pytest.fail("a single frame needs no alignment"))
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
    monkeypatch.setattr(recipe, "_stack_nod_frames", lambda *args: pytest.fail("not in nod mode"))
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
    recipe._stack_nod_frames(*recipe._split_nod_frames(_frames(scene)), _inputs(scene, masterFlat=False))

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


def _stub_alignment(monkeypatch: pytest.MonkeyPatch, alignedA: list[Any], alignedB: list[Any], mapImages: Any) -> None:
    """Replace the slit alignment with one that leaves the caller's lists alone and returns other frames."""
    from soxspipe.commonutils.slit_alignment import slit_alignment_result

    results = [
        slit_alignment_result(tuple(alignedA), pd.DataFrame(), 0, mapImages),
        slit_alignment_result(tuple(alignedB), pd.DataFrame(), 0, mapImages),
    ]
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.slit_alignment"), "align_frame_groups_along_slit", lambda **_: results
    )


def test_the_aligned_frames_come_from_the_alignment_results_not_from_a_change_to_the_input_lists(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene))
    alignedA, alignedB = list(reversed(framesA)), list(reversed(framesB))
    _stub_alignment(
        monkeypatch, alignedA, alignedB, build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    )
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.nod_shift_and_add"), "build_shift_and_add_plan", lambda *_: "PLAN"
    )

    # ACT
    returnedA, returnedB, plan = recipe._align_nod_frames(framesA, framesB, namesA, namesB, "map.fits", "disp.fits")

    # ASSERT
    assert returnedA == alignedA
    assert returnedB == alignedB
    assert plan == "PLAN"


def test_a_nominal_throw_with_no_2d_map_from_the_alignment_is_refused_with_a_clear_error(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, scene: Any
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, monkeypatch, scene)
    framesA, framesB, namesA, namesB = recipe._split_nod_frames(_frames(scene))
    _stub_alignment(monkeypatch, framesA, framesB, None)
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.nod_shift_and_add"), "build_shift_and_add_plan", lambda *_: None
    )

    # ACT / ASSERT
    with pytest.raises(ValueError, match="read no 2D map"):
        recipe._align_nod_frames(framesA, framesB, namesA, namesB, "map.fits", "disp.fits")
