"""A nodding SOF reduced in stare mode reads its frames, names its products and stamps its headers apart (DY-2194).

The stacked path itself is in `test_soxs_stare_nod_stack.py`. These tests run the recipe's `produce_product` with the
collaborators replaced by recorders, as the other stare characterization tests do, and assert what reaches them.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy.nddata import CCDData

from soxspipe.recipes import soxs_stare
from soxspipe.recipes.base_recipe import base_recipe
from tests.factories import synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import _route, _write_prepared_frame
from tests.integration.test_soxs_stare_reduction_characterization import _patch_stare, _stare_recipe

pytestmark = pytest.mark.integration

NOD_TECHNIQUE = "ECHELLE,SLIT,NODDING"
SOF_NAME = "20251212T045200_VIS_1X1_1_NOD_OBJ_SOXS"
MARKER = "_STAREMODE"
SKY_SEEDS = (530, 531, 532)


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
    inputs = (False, False, False, "order.fits", "disp.fits", "map.fits", False)
    recipe.scienceHeaders = []
    recipe.subtractSky = True

    # ACT
    recipe._reduce_stacked_object(notFlattened, [notFlattened], False, inputs, flattened)

    # ASSERT
    assert captured["detrend"] == []
    assert captured["sky"][0]["objectFrame"] is flattened
    assert recipe.dateObs == notFlattened.header["DATE-OBS"]
