"""The merged spectrum carries the wavelength-calibration and resolution cards, in a real FITS file (DY-2358).

`SPEC_RES` comes from the nominal table by slit, and `LAMRMS`, `LAMNLIN` and `SPEC_ERR` from the `WAVE RMS` and
`WAVE NLINES` QCs of the `DISP_TAB` that the recipe already holds. The nod and offset recipes stack their extractions
and flux calibrate the stack, and the stare recipe flux calibrates its extraction. Each writes the cards, and each hands
the same cards to the flux calibration, which rewrites the file.
"""

from __future__ import annotations

import math
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.phase3_check import FAILURE, check_phase3
from soxspipe.recipes import soxs_nod, soxs_offset
from tests.factories import response_product_fits, synthetic_ccd
from tests.factories.disp_tab import disp_tab_with_qcs
from tests.integration.test_merged_spectrum_header import (
    NOD_NAMES,
    OFFSET_NAMES,
    _capture_flux_calibration,
    _frame,
    _recipe,
    _spectrum,
)
from tests.integration.test_observing_recipe_orchestration import _route
from tests.integration.test_soxs_stare_reduction_characterization import _patch_stare
from tests.integration.test_stare_on_target_frames import EXPOSURE_SECONDS, _listed_stare

pytestmark = pytest.mark.integration

WAVE_RMS_NM = 0.0036
WAVE_NLINES = 121
SLIT_CARD = "ESO INS VISE NAME"
WAVELENGTH_CARDS = ("SPEC_RES", "LAMRMS", "LAMNLIN", "SPEC_ERR")


def _nominal_spec_res(log: Any, slit: str) -> int:
    return detector_lookup(log=log, settings={"instrument": "soxs"}).get("VIS")["spec-res"][slit]


def _wave_qcs() -> dict[str, Any]:
    return {"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES}


def _with_slit_and_map(recipe: Any, slit: str | None, dispMap: str | None) -> Any:
    """Give the recipe the slit of its science frames and the `DISP_TAB` it reduced them with."""
    headerOverrides = {key: _frame(0).header[key] for key in _frame(0).header}
    if slit is not None:
        headerOverrides[SLIT_CARD] = slit
    recipe.masterHeaderFrame = synthetic_ccd(prepared=True, headerOverrides=headerOverrides)
    recipe.dispMap = dispMap
    return recipe


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


@pytest.mark.parametrize(
    ("recipeClass", "names"), [(soxs_nod, NOD_NAMES), (soxs_offset, OFFSET_NAMES)], ids=["nod", "offset"]
)
def test_the_stacked_spectrum_file_has_the_resolution_and_wavelength_accuracy_cards(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, recipeClass: type, names: list[str]
) -> None:
    # ARRANGE
    recipe = _with_slit_and_map(
        _recipe(recipeClass, tmp_path, log, monkeypatch, names), "SLIT1.5", disp_tab_with_qcs(tmp_path, **_wave_qcs())
    )

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])

    # ASSERT: THE CARDS ARE IN THE FILE, AND THE CHECKER HAS NOTHING TO SAY ABOUT THEM
    header = fits.getheader(filePath)
    assert header["SPEC_RES"] == _nominal_spec_res(log, "SLIT1.5")
    assert header["LAMRMS"] == WAVE_RMS_NM
    assert header["LAMNLIN"] == WAVE_NLINES
    assert header["SPEC_ERR"] == pytest.approx(WAVE_RMS_NM / math.sqrt(WAVE_NLINES), rel=1e-12)
    findings = check_phase3(filePath)
    assert [f.message for f in findings if any(card in f.message for card in WAVELENGTH_CARDS)] == []


@pytest.mark.parametrize(
    ("recipeClass", "names", "fluxCalibrate"),
    [
        (
            soxs_nod,
            NOD_NAMES,
            lambda recipe, path, stacked: recipe._flux_calibrate_nod_stack(
                path, stacked, [_frame(1), _frame(3)], [_frame(2), _frame(4)]
            ),
        ),
        (
            soxs_offset,
            OFFSET_NAMES,
            lambda recipe, path, stacked: recipe._flux_calibrate_stack(path, stacked, [_frame(1), _frame(3)]),
        ),
    ],
    ids=["nod", "offset"],
)
def test_the_header_given_to_the_flux_calibration_has_the_same_cards(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, recipeClass: type, names: list[str], fluxCalibrate: Any
) -> None:
    # ARRANGE
    recipe = _with_slit_and_map(
        _recipe(recipeClass, tmp_path, log, monkeypatch, names), "SLIT1.0", disp_tab_with_qcs(tmp_path, **_wave_qcs())
    )
    captured = _capture_flux_calibration(monkeypatch, f"soxspipe.recipes.{recipeClass.__name__}")
    stacked, _ = recipe.stack_extractions([_spectrum()])
    responsePath = str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))

    # ACT
    fluxCalibrate(recipe, responsePath, stacked)

    # ASSERT
    header = captured[0]["header"]
    assert header["SPEC_RES"] == _nominal_spec_res(log, "SLIT1.0")
    assert (header["LAMRMS"], header["LAMNLIN"]) == (WAVE_RMS_NM, WAVE_NLINES)


def test_a_dispersion_map_with_no_wave_qcs_leaves_out_the_accuracy_cards_with_a_warning(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A MAP FROM BEFORE THE QCs, WITH THE PIXEL RESIDUAL QCs THAT IT HAS ALWAYS HAD
    recipe = _with_slit_and_map(
        _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES),
        "SLIT1.0",
        disp_tab_with_qcs(tmp_path, **{"X RES MEDIAN": 0.1}),
    )

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])

    # ASSERT: THE RESOLUTION STAYS, THE ACCURACY GOES, AND THE FILE STILL PASSES THE CHECKER'S MANDATORY RULES
    header = fits.getheader(filePath)
    assert header["SPEC_RES"] == _nominal_spec_res(log, "SLIT1.0")
    assert [card for card in ("LAMRMS", "LAMNLIN", "SPEC_ERR") if card in header] == []
    assert any("WAVE RMS" in message for message in _warnings(log))
    assert [f.message for f in check_phase3(filePath) if f.severity == FAILURE and "SPEC_RES" in f.message] == []


@pytest.mark.parametrize("slit", ["BLANK", "UNDEFINED", None])
def test_a_slit_with_no_resolution_leaves_out_spec_res_and_the_checker_fails_the_arm(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, slit: str | None
) -> None:
    # ARRANGE
    recipe = _with_slit_and_map(
        _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES), slit, disp_tab_with_qcs(tmp_path, **_wave_qcs())
    )

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])

    # ASSERT
    assert "SPEC_RES" not in fits.getheader(filePath)
    assert any("SPEC_RES" in message and "arm fails" in message for message in _warnings(log))
    failures = [f.message for f in check_phase3(filePath) if f.severity == FAILURE]
    assert any("SPEC_RES" in message and "missing" in message for message in failures)


def test_a_recipe_with_no_dispersion_map_still_writes_the_resolution(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _with_slit_and_map(_recipe(soxs_offset, tmp_path, log, monkeypatch, OFFSET_NAMES), "SLIT1.0", None)

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])

    # ASSERT
    header = fits.getheader(filePath)
    assert "SPEC_RES" in header
    assert "LAMRMS" not in header


def test_the_stare_flux_calibration_header_has_the_cards_of_the_map_the_extraction_used(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE STARE RECIPE FINDS ITS DISP_TAB IN THE SOF, AND THE SLIT IS ON THE DETRENDED FRAME
    dispMap = disp_tab_with_qcs(tmp_path, **_wave_qcs())
    routes = {
        _route(PRO_CATG="RESP_TAB_VIS"): [str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))],
        _route(PRO_CATG="DISP_TAB_VIS"): [dispMap],
    }
    recipe = _listed_stare(log, tmp_path, count=2, extraRoutes=routes)
    detrended = synthetic_ccd(
        seed=511, prepared=True, headerOverrides={"EXPTIME": EXPOSURE_SECONDS, SLIT_CARD: "SLIT1.5"}
    )
    _patch_stare(recipe, monkeypatch, tmp_path, detrended=detrended)
    arguments: list[dict[str, Any]] = []

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: Any) -> None:
            arguments.append(kwargs)

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
    stamped = arguments[0]["header"]
    assert stamped["SPEC_RES"] == _nominal_spec_res(log, "SLIT1.5")
    assert (stamped["LAMRMS"], stamped["LAMNLIN"]) == (WAVE_RMS_NM, WAVE_NLINES)
