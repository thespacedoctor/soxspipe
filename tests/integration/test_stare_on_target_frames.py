"""The stare recipe counts the frames that reach each of its spectra as on target (DY-2356).

A stacked stare spectrum counts every frame of the stack. With `stack_before_extraction` off, the spectrum of each
frame counts that frame only. The recipe builds the list when it reads the frames, before it edits their headers, and
flux calibration and the response keep dividing by the exposure time of one frame.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_stare
from tests.factories import response_product_fits, synthetic_ccd
from tests.factories.on_target_frames import raw_science_header
from tests.integration.test_observing_recipe_orchestration import _route, _write_prepared_frame
from tests.integration.test_soxs_stare_reduction_characterization import (
    STARE_TECHNIQUE,
    _patch_stare,
    _stare_recipe,
)

pytestmark = pytest.mark.integration

EXPOSURE_SECONDS = 600.0


def _frame_paths(tmpPath: Path, count: int, dprType: str = "OBJECT") -> list[str]:
    """Write ``count`` prepared stare frames whose ARCFILE differs from the file name, and return their paths."""
    paths = []
    for index in range(count):
        header = raw_science_header(
            "VIS", mjd=61000.0 + index / 10, arcfile=f"SOXS.{index}.fits", exptime=EXPOSURE_SECONDS
        )
        overrides = {key: header[key] for key in header}
        overrides.update({"DPR_TYPE": dprType, "DPR_TECH": STARE_TECHNIQUE})
        paths.append(
            str(_write_prepared_frame(tmpPath / f"stare-{index}_pre.fits", seed=520 + index, headerOverrides=overrides))
        )
    return paths


def _listed_stare(
    log: Any,
    tmpPath: Path,
    count: int = 2,
    *,
    dprType: str = "OBJECT",
    extraRoutes: dict | None = None,
    **settings: Any,
) -> soxs_stare:
    paths = _frame_paths(tmpPath, count, dprType)
    # THE OBJECT ROUTE COMES FIRST, SO IT IS EMPTIED WHEN THE FRAMES ARE OF ANOTHER TYPE
    routes = {
        _route(DPR_TYPE="OBJECT", DPR_TECH=STARE_TECHNIQUE): [],
        _route(DPR_TYPE=dprType, DPR_TECH=STARE_TECHNIQUE): paths,
        **(extraRoutes or {}),
    }
    recipe = _stare_recipe(log, tmpPath, extraRoutes=routes)
    recipe.recipeSettings = {**recipe.recipeSettings, **settings}
    return recipe


def _names(entries: Any) -> list[str]:
    return [entry.name for entry in entries]


def test_every_frame_of_a_stacked_stare_is_on_target(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    recipe = _listed_stare(log, tmp_path, count=3)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["extractor"]) == 1
    entries = captured["extractor"][0]["onTargetFrames"]
    assert _names(entries) == ["SOXS.0.fits", "SOXS.1.fits", "SOXS.2.fits"]
    assert [entry.exposureSeconds for entry in entries] == [EXPOSURE_SECONDS] * 3


def test_the_list_is_built_before_the_recipe_edits_the_frame_headers(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE HEADER SCRUB OF A PRODUCT WRITE POPS ARCFILE
    recipe = _listed_stare(log, tmp_path, count=1)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    def scrubbing_stack(**kwargs: Any) -> Any:
        frame = kwargs["frames"][0]
        frame.header.pop("ARCFILE")
        return frame

    monkeypatch.setattr(recipe, "clip_and_stack", scrubbing_stack)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _names(captured["extractor"][0]["onTargetFrames"]) == ["SOXS.0.fits"]


def test_with_stacking_off_the_spectrum_of_each_frame_counts_that_frame_only(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _listed_stare(log, tmp_path, count=2, **{"stack_before_extraction": False})
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [_names(call["onTargetFrames"]) for call in captured["extractor"]] == [["SOXS.0.fits"], ["SOXS.1.fits"]]
    assert _names(recipe.onTargetFrames) == ["SOXS.0.fits", "SOXS.1.fits"]


def test_a_standard_response_uses_the_exposure_time_of_one_frame_not_the_sum(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _listed_stare(log, tmp_path, count=2, dprType="STD,FLUX")
    recipe.generateReponseCurve = True
    _patch_stare(recipe, monkeypatch, tmp_path)
    responseArguments: list[dict[str, Any]] = []

    class FakeResponse:
        def __init__(self, **kwargs: Any) -> None:
            responseArguments.append(kwargs)

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, bool]:
            return recipe.qc, recipe.products, False

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "response_function", FakeResponse)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert responseArguments[0]["exptime"] == EXPOSURE_SECONDS


def test_flux_calibration_divides_by_one_frame_while_its_header_carries_the_summed_exposure(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    routes = {_route(PRO_CATG="RESP_TAB_VIS"): [str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))]}
    recipe = _listed_stare(log, tmp_path, count=2, extraRoutes=routes)
    # THE DETRENDED STACK CARRIES THE HEADER OF ITS FIRST FRAME, SO ITS EXPTIME IS THAT OF ONE FRAME
    detrended = synthetic_ccd(seed=511, prepared=True, headerOverrides={"EXPTIME": EXPOSURE_SECONDS})
    _patch_stare(recipe, monkeypatch, tmp_path, detrended=detrended)
    calibratorArguments: list[dict[str, Any]] = []

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: Any) -> None:
            calibratorArguments.append(kwargs)

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            return "fluxcal.fits", pd.DataFrame()

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "flux_calibration", FakeFluxCalibrator)
    monkeypatch.setattr(import_module("soxspipe.recipes.soxs_stare"), "get_calibrations_path", lambda **_: "/cal")
    monkeypatch.setattr(
        import_module("soxspipe.recipes.soxs_stare"),
        "detector_lookup",
        lambda **_: type("Lookup", (), {"get": staticmethod(lambda arm: {"extinction": "extinction.dat"})})(),
    )
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *_, **__: None)

    # ACT
    recipe.produce_product()

    # ASSERT
    arguments = calibratorArguments[0]
    assert arguments["exptime"] == EXPOSURE_SECONDS
    assert arguments["header"]["EXPTIME"] == 2 * EXPOSURE_SECONDS
    assert arguments["header"]["NCOMBINE"] == 2
    assert (arguments["header"]["PROV1"], arguments["header"]["PROV2"]) == ("SOXS.0", "SOXS.1")
