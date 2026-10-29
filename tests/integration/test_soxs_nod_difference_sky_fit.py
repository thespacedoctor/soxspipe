"""The sky fit on the nod A-B difference frame, behind the `difference_sky_fit` setting (DY-2057)."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils.difference_sky_fit import DIFFERENCE_VARIANCE_FACTOR, MIN_MODEL_SIGNIFICANCE
from soxspipe.recipes import soxs_nod, soxs_offset
from tests.integration.test_soxs_nod_cycle_characterization import (
    _ab_frames,
    _cycle_recipe,
    _patch_detrend,
    _patch_extraction_collaborators,
)

pytestmark = pytest.mark.integration

# THE COMBINED FRAME IS THE MEAN OF TWO INDEPENDENT DIFFERENCES, SO ITS VARIANCE IS THAT OF ONE FRAME
COMBINED_VARIANCE_FACTOR = 1.0
STARE_SKY_SETTINGS = {
    "bspline_order": 3,
    "subtract_sky": False,
    "sky_model_qc_plot": True,
    "object_wing_correction": True,
}


def _fit_recipe(log: Any, tmp_path: Path, *, recipeName: str = "soxs-nod", **recipeSettings: Any) -> soxs_nod:
    recipe = _cycle_recipe(log, tmp_path, recipeName=recipeName)
    recipe.arm = "NIR"
    recipe.readNoise = None
    recipe.recipeSettings = {**recipe.recipeSettings, **recipeSettings}
    recipe.settings = {"instrument": "soxs", "soxs-stare": {"nir": {"sky-subtraction": dict(STARE_SKY_SETTINGS)}}}
    return recipe


def _patch_sky_fit(monkeypatch: pytest.MonkeyPatch, *, returnsModel: bool = True) -> list[dict[str, Any]]:
    """Replace `subtract_sky` with a recorder whose sky-subtracted frame is a copy of the frame it was given."""
    commonutils = import_module("soxspipe.commonutils")
    calls: list[dict[str, Any]] = []

    class RecordingSkyFit:
        def __init__(self, **kwargs: Any) -> None:
            calls.append(kwargs)
            self.kwargs = kwargs

        def subtract(self) -> tuple[Any, Any, Any, Any, Any]:
            frame = self.kwargs["objectFrame"]
            if not returnsModel:
                return None, None, None, self.kwargs["qcTable"], self.kwargs["productsTable"]
            fitted = frame.copy()
            calls[-1]["fitted"] = fitted
            return frame.copy(), fitted, frame.copy(), self.kwargs["qcTable"], self.kwargs["productsTable"]

    monkeypatch.setattr(commonutils, "subtract_sky", RecordingSkyFit)
    return calls


def _run_cycle(recipe: soxs_nod, *, notFlattened: bool = False, masterFlat: Any = False) -> None:
    aFrame, bFrame = _ab_frames()
    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        notFlattened=notFlattened,
        masterFlat=masterFlat,
    )


def test_with_the_setting_absent_no_sky_is_fitted_and_the_differences_are_extracted_as_before(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)

    # ACT
    _run_cycle(recipe)

    # ASSERT
    aFrame, bFrame = _ab_frames()
    assert fitCalls == []
    assert len(extractCalls) == 2
    pd.testing.assert_frame_equal(
        pd.DataFrame(extractCalls[0]["skySubtractedFrame"].data), pd.DataFrame(aFrame.data - bFrame.data)
    )


def test_with_the_setting_false_no_sky_is_fitted(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=False)
    _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)

    # ACT
    _run_cycle(recipe)

    # ASSERT
    assert fitCalls == []


def test_a_nod_cycle_fits_the_sky_of_each_difference_with_the_frame_it_subtracted_as_the_noise_sky(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)
    aFrame, bFrame = _ab_frames()

    # ACT
    _run_cycle(recipe)

    # ASSERT
    assert len(fitCalls) == 2
    assert fitCalls[0]["differenceFit"].noiseSkyFrame.data.tolist() == bFrame.data.tolist()
    assert fitCalls[1]["differenceFit"].noiseSkyFrame.data.tolist() == aFrame.data.tolist()
    assert [call["differenceFit"].varianceFactor for call in fitCalls] == [DIFFERENCE_VARIANCE_FACTOR] * 2
    assert extractCalls[0]["skySubtractedFrame"] is fitCalls[0]["fitted"]
    assert extractCalls[1]["skySubtractedFrame"] is fitCalls[1]["fitted"]


def test_a_combined_frame_is_fitted_once_with_the_variance_of_one_frame(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)
    aFrame, bFrame = _ab_frames()
    combined, combinedUnflattened, skyFrame = aFrame.copy(), aFrame.copy(), bFrame.copy()
    combinedFrames = (combined, combinedUnflattened, skyFrame)

    # ACT
    recipe._extract_cycle_spectra(
        aFrame, bFrame, aFrame, bFrame, aFrame, bFrame, False, 1, False, combinedFrames=combinedFrames
    )

    # ASSERT
    assert len(fitCalls) == 1
    assert fitCalls[0]["objectFrame"] is combined
    assert fitCalls[0]["differenceFit"].noiseSkyFrame is skyFrame
    assert fitCalls[0]["differenceFit"].varianceFactor == COMBINED_VARIANCE_FACTOR
    assert extractCalls[0]["skySubtractedFrame"] is fitCalls[0]["fitted"]


def test_the_flattened_pass_fits_the_flat_fielded_frame_and_keeps_the_trace_frame_as_it_was(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)
    aFrame, bFrame = _ab_frames()
    combined, combinedUnflattened = aFrame.copy(), aFrame.copy()
    masterFlat = bFrame.copy()

    # ACT
    recipe._extract_cycle_spectra(
        aFrame,
        bFrame,
        aFrame,
        bFrame,
        aFrame,
        bFrame,
        False,
        1,
        False,
        combinedFrames=(combined, combinedUnflattened, bFrame),
        masterFlat=masterFlat,
    )

    # ASSERT
    assert fitCalls[0]["masterFlat"] is masterFlat
    assert extractCalls[0]["unflattenedFrame"] is combinedUnflattened


def test_the_not_flattened_pass_fits_the_unflattened_frame_without_a_flat_and_extracts_the_fit(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)
    aFrame, bFrame = _ab_frames()
    combined, combinedUnflattened = aFrame.copy(), aFrame.copy()

    # ACT
    recipe._extract_cycle_spectra(
        aFrame,
        bFrame,
        aFrame,
        bFrame,
        aFrame,
        bFrame,
        False,
        1,
        True,
        combinedFrames=(combined, combinedUnflattened, bFrame),
        masterFlat=bFrame.copy(),
    )

    # ASSERT
    assert fitCalls[0]["objectFrame"] is combinedUnflattened
    assert fitCalls[0]["masterFlat"] is None
    assert extractCalls[0]["unflattenedFrame"] is fitCalls[0]["fitted"]
    assert extractCalls[0]["skySubtractedFrame"] is fitCalls[0]["fitted"]


def test_the_fit_uses_the_stare_sky_settings_of_the_arm_with_sky_subtraction_on_and_its_plot_off(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)

    # ACT
    _run_cycle(recipe)

    # ASSERT
    skySettings = fitCalls[0]["recipeSettings"]["sky-subtraction"]
    assert skySettings["bspline_order"] == 3
    assert skySettings["subtract_sky"] is True
    assert skySettings["sky_model_qc_plot"] is False
    assert recipe.settings["soxs-stare"]["nir"]["sky-subtraction"] == STARE_SKY_SETTINGS


def test_the_minimum_significance_setting_reaches_the_fit_and_defaults_to_the_shipped_value(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    defaultRecipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    tunedRecipe = _fit_recipe(log, tmp_path, difference_sky_fit=True, difference_sky_fit_min_significance=5.0)
    for recipe in (defaultRecipe, tunedRecipe):
        _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)

    # ACT
    _run_cycle(defaultRecipe)
    _run_cycle(tunedRecipe)

    # ASSERT
    assert fitCalls[0]["differenceFit"].minSignificance == MIN_MODEL_SIGNIFICANCE
    assert fitCalls[2]["differenceFit"].minSignificance == 5.0


def test_a_fit_that_makes_no_model_leaves_the_difference_to_be_extracted_as_it_was(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    _patch_sky_fit(monkeypatch, returnsModel=False)
    aFrame, bFrame = _ab_frames()

    # ACT
    _run_cycle(recipe)

    # ASSERT
    assert len(extractCalls) == 2
    assert extractCalls[0]["skySubtractedFrame"].data.tolist() == (aFrame.data - bFrame.data).tolist()
    assert any(level == "warning" for level, _ in log.messages)


def test_an_offset_run_fits_its_one_difference_through_the_same_code_path(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, recipeName="soxs-offset", difference_sky_fit=True)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)

    # ACT
    _run_cycle(recipe)

    # ASSERT
    assert soxs_offset._extract_cycle_spectra is soxs_nod._extract_cycle_spectra
    assert soxs_offset.process_single_ab_nodding_cycle is soxs_nod.process_single_ab_nodding_cycle
    assert len(extractCalls) == 1
    assert len(fitCalls) == 1


def test_a_detrended_cycle_passes_its_master_flat_to_the_flattened_fit(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _fit_recipe(log, tmp_path, difference_sky_fit=True)
    _patch_extraction_collaborators(monkeypatch, recipe)
    _patch_detrend(monkeypatch, recipe)
    fitCalls = _patch_sky_fit(monkeypatch)
    masterFlat = _ab_frames()[0]

    # ACT
    _run_cycle(recipe, masterFlat=masterFlat)

    # ASSERT
    assert [call["masterFlat"] for call in fitCalls] == [masterFlat, masterFlat]


@pytest.mark.parametrize("recipeName", ["soxs-nod", "soxs-nod-std"])
def test_the_shipped_settings_have_the_difference_sky_fit_off(recipeName: str) -> None:
    # ARRANGE
    import yaml

    import soxspipe

    settingsPath = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"

    # ACT
    with settingsPath.open() as stream:
        recipeSettings = yaml.safe_load(stream)[recipeName]

    # ASSERT
    assert recipeSettings["difference_sky_fit"] is False
    assert recipeSettings["difference_sky_fit_min_significance"] == MIN_MODEL_SIGNIFICANCE
