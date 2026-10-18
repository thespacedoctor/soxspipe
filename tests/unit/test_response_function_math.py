"""Contracts for the response-function recipe built on the spline fit."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import response_fit
from soxspipe.commonutils.response_function import _ResponseFitConvergenceError, response_function

pytestmark = pytest.mark.unit
responseModule = importlib.import_module("soxspipe.commonutils.response_function")
toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
phase3Module = importlib.import_module("soxspipe.commonutils.phase3")


def test_response_get_rejects_an_unknown_standard_star(log: object) -> None:
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame({"WAVE": [500.0], "FLUX_COUNTS": [10.0]})
    response.stdExtractionNotFlatDF = pd.DataFrame({"WAVE": [500.0], "FLUX_DENSITY_COUNTS": [10.0]})
    response.stdAbsFluxDF = pd.DataFrame({"WAVE": [500.0], "KNOWN_STAR": [1.0]})
    response.std_objName = "UNKNOWN_STAR"

    with pytest.raises(
        LookupError,
        match="Standard star UNKNOWN_STAR not found.*KNOWN_STAR",
    ):
        response.get()


def test_the_convergence_error_is_shared_with_the_fit_module() -> None:
    assert _ResponseFitConvergenceError is response_fit._ResponseFitConvergenceError


def _vis_response(log: object, monkeypatch: pytest.MonkeyPatch, absFlux: float = 2.0) -> response_function:
    """Return a VIS response_function over a flat synthetic standard, with disk writers captured."""
    wavelengths = np.linspace(400.0, 800.0, 4000)
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame({"WAVE": wavelengths, "FLUX_COUNTS": np.full(4000, 10.0)})
    response.stdExtractionNotFlatDF = pd.DataFrame({"WAVE": wavelengths, "FLUX_DENSITY_COUNTS": np.full(4000, 10.0)})
    response.stdAbsFluxDF = pd.DataFrame({"WAVE": wavelengths, "SYNTHETIC_STAR": np.full(4000, absFlux)})
    response.std_objName = "SYNTHETIC_STAR"
    response.arm = "VIS"
    response.instrument = "soxs"
    response.texp = 10.0
    response.recipeSettings = {"vis": {}}
    response.calibrationRootPath = "/calibration"
    response.detectorParams = {"extinction": "extinction.fits"}
    response.airmass = 1.0
    response.qc = pd.DataFrame()
    response.products = pd.DataFrame()
    response.recipeName = "soxs-response"
    response.sofName = "SYNTHETIC"
    response.qcDir = "/qc"
    response.productDir = "/products"
    response.header = {}
    response.kw = lambda key: key
    response.settings = {"instrument": "soxs"}
    response.dateObs = "2024-01-01T00:00:00"
    response.orderJoins = []
    monkeypatch.setattr(
        responseModule,
        "extinction_correction_factor",
        lambda wavelengths, *args: np.ones_like(wavelengths),
    )
    monkeypatch.setattr(
        toolkitModule,
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    monkeypatch.setattr(phase3Module, "write_fits_table_to_disk", lambda **kwargs: None)
    return response


def test_response_get_writes_the_spline_curve_and_legacy_coefficients(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    response = _vis_response(log, monkeypatch)
    captured: dict[str, object] = {}
    response.write_response_function_to_file = lambda **kwargs: captured.update({"responseWrite": kwargs})
    response.plot_response_curve = lambda **kwargs: captured.update({"plot": kwargs})

    # ACT
    qc, products, recipeError = response.get()

    # ASSERT: A FLAT RATIO OF 2E17 GIVES A FLAT SPLINE AND AN ORDER-5 LEGACY POLYNOMIAL THAT IS FLAT TOO
    write = captured["responseWrite"]
    assert recipeError is False
    assert write["polyOrder"] == 5
    assert len(write["responseFuncCoeffs"]) == 6
    np.testing.assert_allclose(write["response"], 2.0e17, rtol=1e-6)
    np.testing.assert_allclose(np.polyval(write["responseFuncCoeffs"], 600.0), 2.0e17, rtol=1e-6)
    np.testing.assert_array_equal(write["wavelength"], response.stdExtractionDF["WAVE"].values)
    assert products["product_label"].tolist() == ["EFFICIENCY"]
    assert captured["plot"]["qcWarnings"] == []


def test_response_get_adds_the_response_fit_qcs_before_writing_the_product(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    response = _vis_response(log, monkeypatch)
    seenAtWrite: dict[str, list[str]] = {}
    response.write_response_function_to_file = lambda **kwargs: seenAtWrite.update(
        {"names": response.qc["qc_name"].tolist()}
    )
    response.plot_response_curve = lambda **kwargs: None

    # ACT
    response.get()

    # ASSERT
    assert seenAtWrite["names"][:2] == ["RESP BIN CLIP FRAC", "RESP MAX GAP"]
    assert "RESP SELFCAL RES 350-370" not in seenAtWrite["names"]
    assert "RESP SELFCAL RES 500-600" in seenAtWrite["names"]
    selfCal = response.qc.loc[response.qc["qc_name"] == "RESP SELFCAL RES 500-600", "qc_value"].iloc[0]
    assert abs(selfCal) < 1e-3


def test_response_get_warns_and_passes_the_warnings_to_the_plot_when_qc_is_out_of_tolerance(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: A TOLERANCE OF ZERO ON THE MAX GAP FORCES A BREACH WITHOUT CHANGING THE FIT
    response = _vis_response(log, monkeypatch)
    response.recipeSettings = {"vis": {"max_gap_nm_max": 1}}
    captured: dict[str, object] = {}
    response.write_response_function_to_file = lambda **kwargs: None
    response.plot_response_curve = lambda **kwargs: captured.update(kwargs)

    # ACT
    _, _, recipeError = response.get()

    # ASSERT: THE RECIPE STILL SUCCEEDS
    assert recipeError is False
    assert len(captured["qcWarnings"]) == 1
    assert any(level == "warning" and "gap" in message for level, message in log.messages)


def test_response_get_raises_when_the_response_has_no_finite_samples(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: ZERO EXTRACTED FLUX EVERYWHERE
    response = _vis_response(log, monkeypatch)
    response.stdExtractionDF["FLUX_COUNTS"] = 0.0

    # ACT AND ASSERT
    with pytest.raises(Exception, match="did not converge"):
        response.get()
