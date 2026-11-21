"""A nodding SOF reduced in stare mode reads its frames, names its products and stamps its headers apart (DY-2194).

The stacked path itself is in `test_soxs_stare_nod_stack.py`. These tests run the recipe's `produce_product` with the
collaborators replaced by recorders, as the other stare characterization tests do, and assert what reaches them.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.nddata import CCDData

from soxspipe.recipes import soxs_stare
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.nod_stare_stacking import reduction_inputs
from tests.factories import synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import _route, _write_prepared_frame
from tests.integration.test_soxs_stare_reduction_characterization import _patch_stare, _stare_recipe

pytestmark = pytest.mark.integration

NOD_TECHNIQUE = "ECHELLE,SLIT,NODDING"
SOF_NAME = "20251212T045200_VIS_1X1_1_NOD_OBJ_SOXS"
MARKER = "_STAREMODE"
SKY_SEEDS = (530, 531, 532)
# THE COLUMNS AT THE HIGH END OF THE SLIT THAT THE TAIL TESTS TREAT AS UNCOVERED
UNCOVERED_COLUMNS = 3


def _nod_stare_recipe(
    log: Any,
    tmpPath: Path,
    *,
    frameType: str = "OBJECT",
    frameCount: int = 1,
    recipeSettings: dict[str, object] | None = None,
) -> soxs_stare:
    """Return a stare recipe in nod mode, routed to nodding frames of one type."""
    paths = [
        str(
            _write_prepared_frame(
                tmpPath / f"nod-object-{index}_pre.fits",
                seed=600 + index,
                headerOverrides={"DPR_TYPE": frameType, "DPR_TECH": NOD_TECHNIQUE},
            )
        )
        for index in range(frameCount)
    ]
    recipe = _stare_recipe(
        log,
        tmpPath,
        extraRoutes={_route(DPR_TYPE=frameType, DPR_TECH=NOD_TECHNIQUE): paths},
        recipeSettings=recipeSettings or {"use_flat": False, "sky-subtraction": {"subtract_sky": True}},
    )
    recipe.reduceNodAsStare = True
    recipe.productNameMarker = MARKER
    recipe.sofName = SOF_NAME
    recipe.filenameTemplate = f"{SOF_NAME}{MARKER}.fits"
    return recipe


def _sky_result() -> tuple[CCDData, CCDData, CCDData]:
    return tuple(synthetic_ccd(seed=seed, prepared=True) for seed in SKY_SEEDS)


def test_stare_mode_reads_the_nodding_frame_and_reduces_it_as_a_stare_frame(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.readStareFrames is True
    assert len(captured["stack"]) == 1
    assert len(captured["stack"][0]["frames"]) == 1
    assert len(captured["sky"]) == 1
    assert len(captured["extractor"]) == 1


def test_the_sky_extraction_and_sky_products_are_named_with_the_marker(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["sky"][0]["sofName"] == SOF_NAME + MARKER
    assert captured["extractor"][0]["sofName"] == SOF_NAME + MARKER
    assert [call["filename"] for call in captured["write"]] == [
        f"{SOF_NAME}{MARKER}_SKYSUB.fits",
        f"{SOF_NAME}{MARKER}_SKYMODEL.fits",
        f"{SOF_NAME}{MARKER}_SKYSUB_RESIDUALS.fits",
    ]


def test_the_recipe_still_reduces_under_the_stare_recipe_name(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["sky"][0]["recipeName"] == "soxs-stare"
    assert captured["extractor"][0]["recipeName"] == "soxs-stare"


def test_the_flux_calibration_and_the_response_curve_are_named_with_the_marker(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    recipe.dateObs = "2024-01-02T03:04:05"
    commonutils = import_module("soxspipe.commonutils")
    stareModule = import_module("soxspipe.recipes.soxs_stare")
    airmass = import_module("soxspipe.commonutils.sequence_airmass")
    received: dict[str, dict[str, object]] = {}

    class FakeFluxCalibration:
        def __init__(self, **kwargs: object) -> None:
            received["flux"] = kwargs

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            return "fluxcal.fits", pd.DataFrame()

    class FakeResponse:
        def __init__(self, **kwargs: object) -> None:
            received["response"] = kwargs

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
            return recipe.qc, recipe.products, False

    class FakeExtractor:
        def __init__(self, **kwargs: object) -> None:
            pass

        def extract(self) -> tuple[pd.DataFrame, pd.DataFrame, None, None, str]:
            return recipe.qc, recipe.products, None, None, "extraction.fits"

    class FakeDetectorLookup:
        def __init__(self, **kwargs: object) -> None:
            pass

        def get(self, arm: str) -> dict[str, str]:
            return {"extinction": "extinction.txt"}

    frame = synthetic_ccd(seed=1, prepared=True)
    monkeypatch.setattr(commonutils, "flux_calibration", FakeFluxCalibration)
    monkeypatch.setattr(commonutils, "response_function", FakeResponse)
    monkeypatch.setattr(commonutils, "horne_extraction", FakeExtractor)
    monkeypatch.setattr(stareModule, "detector_lookup", FakeDetectorLookup)
    monkeypatch.setattr(stareModule, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(airmass, "require_equal_exposure_times", lambda *_: None)
    monkeypatch.setattr(airmass, "mean_airmass", lambda *_: 1.2)
    recipe.startNightDate = "2024-01-02"
    recipe.scienceHeaders = []

    # ACT
    recipe._flux_calibrate_spectrum("resp.fits", pd.DataFrame(), frame, [frame])
    recipe._generate_response_curve(frame, False, "map.fits", "disp.fits", "ext.fits", {}, objectFrames=[frame])

    # ASSERT
    assert received["flux"]["sofName"] == SOF_NAME + MARKER
    assert received["response"]["sofName"] == SOF_NAME + MARKER


def test_a_recipe_without_the_marker_names_everything_by_the_sof_name_as_before(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path, recipeSettings={"use_flat": False, "sky-subtraction": {"subtract_sky": True}})
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["sky"][0]["sofName"] == "OBJECT_VIS"
    assert captured["extractor"][0]["sofName"] == "OBJECT_VIS"
    assert captured["write"][0]["filename"] == "OBJECT_VIS_SKYSUB.fits"


def test_each_frame_of_a_no_stack_run_is_named_frame_then_marker_and_keyed_by_the_bare_frame_name(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(
        log,
        tmp_path,
        frameCount=2,
        recipeSettings={"use_flat": False, "sky-subtraction": {"subtract_sky": True}, "stack_before_extraction": False},
    )
    recipe.settings = {**recipe.settings, "no-stack": True}
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())
    reportedKeys: list[str] = []
    monkeypatch.setattr(recipe, "report_output", lambda: reportedKeys.append(recipe.sofName) or recipe.qc)

    # ACT
    recipe.produce_product()

    # ASSERT: PRODUCT NAMES CARRY THE FRAME NUMBER THEN THE MARKER. THE DATABASE KEYS CARRY NO MARKER
    assert [call["sofName"] for call in captured["sky"]] == [f"{SOF_NAME}_F1{MARKER}", f"{SOF_NAME}_F2{MARKER}"]
    assert [call["filename"] for call in captured["write"]][::3] == [
        f"{SOF_NAME}_F1{MARKER}_SKYSUB.fits",
        f"{SOF_NAME}_F2{MARKER}_SKYSUB.fits",
    ]
    assert reportedKeys == [f"{SOF_NAME}_F1", f"{SOF_NAME}_F2"]
    assert recipe.sofName == SOF_NAME
    assert recipe.filenameTemplate == f"{SOF_NAME}{MARKER}.fits"


def test_the_header_of_a_stare_mode_product_records_the_technique_and_the_mode(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    recipe.kw = import_module("soxspipe.commonutils").keyword_lookup(log=log, settings=recipe.settings).get
    monkeypatch.setattr(base_recipe, "update_fits_keywords", lambda self, frame, rawFrames=False: None)
    frame = synthetic_ccd(seed=2, prepared=True)
    frame.header["ESO PRO TECH"] = "IMAGE"

    # ACT
    recipe.update_fits_keywords(frame=frame)

    # ASSERT
    assert frame.header["ESO PRO TECH"] == "ECHELLE,SLIT,STARE"
    assert frame.header["ESO PRO STARE MODE"] is True


def test_the_header_of_a_normal_stare_product_gets_no_stare_mode_keyword(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _stare_recipe(log, tmp_path)
    recipe.kw = import_module("soxspipe.commonutils").keyword_lookup(log=log, settings=recipe.settings).get
    monkeypatch.setattr(base_recipe, "update_fits_keywords", lambda self, frame, rawFrames=False: None)
    frame = synthetic_ccd(seed=2, prepared=True)

    # ACT
    recipe.update_fits_keywords(frame=frame)

    # ASSERT
    assert "ESO PRO STARE MODE" not in frame.header
    assert "ESO PRO TECH" not in frame.header


def test_the_tail_uses_a_stack_that_was_detrended_before_stacking_without_detrending_it_again(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _nod_stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=_sky_result())
    notFlattened = synthetic_ccd(seed=3, prepared=True)
    flattened = synthetic_ccd(seed=4, prepared=True)
    inputs = reduction_inputs(False, False, False, "order.fits", "disp.fits", "map.fits", False)
    recipe.scienceHeaders = []
    recipe.subtractSky = True

    # ACT
    recipe._reduce_stacked_object(notFlattened, [notFlattened], False, inputs, flattened)

    # ASSERT
    assert captured["detrend"] == []
    assert captured["sky"][0]["objectFrame"] is flattened
    assert recipe.dateObs == notFlattened.header["DATE-OBS"]


def _tail_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, skyResult: tuple[Any, ...] | None = None
) -> tuple[soxs_stare, dict[str, list[dict[str, object]]], CCDData]:
    """Return a nod-mode stare recipe ready to run the tail, with its collaborators recorded, and the stacked frame."""
    recipe = _nod_stare_recipe(log, tmp_path)
    captured = _patch_stare(recipe, monkeypatch, tmp_path, skyResult=skyResult or _sky_result())
    recipe.scienceHeaders = []
    recipe.subtractSky = True
    recipe.slitIlluminationOffsets = None
    return recipe, captured, synthetic_ccd(seed=3, prepared=True)


def _tail_inputs() -> tuple[Any, ...]:
    return reduction_inputs(False, False, False, "order.fits", "disp.fits", "map.fits", False)


def _uncovered_columns(shape: tuple[int, ...]) -> np.ndarray:
    """Return a mask over the last three columns, which the stack does not cover in every frame."""
    uncovered = np.zeros(shape, dtype=bool)
    uncovered[:, -UNCOVERED_COLUMNS:] = True
    return uncovered


def test_the_pixels_not_every_frame_covers_are_left_out_of_the_sky_fit_without_joining_the_exclude_mask(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)
    uncovered = _uncovered_columns(stack.data.shape)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack, uncovered)

    # ASSERT: THE SKY FIT GETS THE TWO MASKS APART, AS ONLY THE EXCLUDE MASK SWITCHES OFF THE SLIT-EDGE CLIP
    np.testing.assert_array_equal(captured["sky"][0]["uncoveredMask"], uncovered)
    assert captured["sky"][0]["excludeMask"] is None


def test_the_pixels_not_every_frame_covers_are_flagged_for_the_extraction_when_no_sky_is_subtracted(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: WITH A SKY MODEL THE SKY SUBTRACTION FLAGS THE EXCLUDED PIXELS; WITHOUT ONE THE RECIPE DOES
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)
    recipe.subtractSky = False
    uncovered = _uncovered_columns(stack.data.shape)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack, uncovered)

    # ASSERT
    assert captured["extractor"][0]["skySubtractedFrame"].mask[uncovered].all()
    assert captured["extractor"][0]["unflattenedFrame"].mask[uncovered].all()


def test_the_pixels_not_every_frame_covers_are_flagged_for_the_extraction_when_the_sky_fit_gives_no_model(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE SKY SUBTRACTION RETURNS NO MODEL, SO THE RECIPE FLAGS THE FRAME ITSELF
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch, skyResult=(None, None, None))
    uncovered = _uncovered_columns(stack.data.shape)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack, uncovered)

    # ASSERT
    assert captured["extractor"][0]["skySubtractedFrame"].mask[uncovered].all()


def test_the_slit_illumination_exclude_mask_reaches_the_sky_fit_beside_the_uncovered_mask(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)
    illumination = np.zeros(stack.data.shape, dtype=bool)
    illumination[:2, :] = True
    uncovered = _uncovered_columns(stack.data.shape)
    assert not illumination[:, -UNCOVERED_COLUMNS:].all()
    monkeypatch.setattr(recipe, "_slit_illumination_exclude_mask", lambda *args: illumination)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack, uncovered)

    # ASSERT
    np.testing.assert_array_equal(captured["sky"][0]["excludeMask"], illumination)
    np.testing.assert_array_equal(captured["sky"][0]["uncoveredMask"], uncovered)


def test_the_no_sky_extraction_is_flagged_with_both_masks(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)
    recipe.subtractSky = False
    illumination = np.zeros(stack.data.shape, dtype=bool)
    illumination[:2, :] = True
    uncovered = _uncovered_columns(stack.data.shape)
    monkeypatch.setattr(recipe, "_slit_illumination_exclude_mask", lambda *args: illumination)
    expected = np.asarray(stack.mask, dtype=bool) | illumination | uncovered

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack, uncovered)

    # ASSERT
    np.testing.assert_array_equal(captured["extractor"][0]["skySubtractedFrame"].mask, expected)


def test_a_stack_with_no_uncovered_mask_gives_the_sky_fit_the_exclude_mask_as_it_was(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A GENUINE STARE FRAME NEVER GETS AN UNCOVERED MASK
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)
    illumination = np.zeros(stack.data.shape, dtype=bool)
    illumination[:2, :] = True
    monkeypatch.setattr(recipe, "_slit_illumination_exclude_mask", lambda *args: illumination)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack)

    # ASSERT
    assert captured["sky"][0]["excludeMask"] is illumination
    assert captured["sky"][0]["uncoveredMask"] is None


def test_a_stack_with_neither_mask_gives_the_sky_fit_no_mask(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, captured, stack = _tail_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe._reduce_stacked_object(stack, [stack], False, _tail_inputs(), stack)

    # ASSERT
    assert captured["sky"][0]["excludeMask"] is None
    assert captured["sky"][0]["uncoveredMask"] is None
