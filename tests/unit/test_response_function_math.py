"""Analytic contracts for response-function polynomial fitting."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.response_function import (
    _fit_response_polynomial,
    _ResponseFitConvergenceError,
    response_function,
)

pytestmark = pytest.mark.unit
responseModule = importlib.import_module("soxspipe.commonutils.response_function")
toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
phase3Module = importlib.import_module("soxspipe.commonutils.phase3")


def _legacy_fit_response_polynomial(
    wavelength: np.ndarray,
    rawResponse: np.ndarray,
    polynomialOrder: int,
    maxIterations: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fittedWavelength = wavelength.copy()
    fittedResponse = rawResponse.copy()
    iteration = 0
    deletedPointCount = 1
    while iteration < maxIterations and deletedPointCount > 0:
        sampledWavelength = np.array(
            [
                np.median(fittedWavelength[max(0, index - 5) : index + 6])
                for index in range(0, len(fittedWavelength), 100)
            ]
        )
        sampledResponse = np.array(
            [
                np.median(fittedResponse[max(0, index - 5) : index + 6])
                for index in range(0, len(fittedResponse), 100)
            ]
        )
        coefficients = np.polyfit(
            sampledWavelength,
            sampledResponse,
            deg=polynomialOrder,
        )
        modelResponse = np.polyval(coefficients, fittedWavelength)
        deletedPoints = [
            index
            for index, (responseValue, modelValue) in enumerate(
                zip(fittedResponse, modelResponse)
            )
            if responseValue < 0
            or abs(abs(responseValue) - abs(modelValue)) / abs(responseValue) > 0.2
        ]
        fittedWavelength = np.delete(fittedWavelength, deletedPoints)
        fittedResponse = np.delete(fittedResponse, deletedPoints)
        deletedPointCount = len(deletedPoints)
        iteration += 1
    return coefficients, fittedWavelength, fittedResponse


def test_response_fit_recovers_a_linear_curve_and_rejects_outliers() -> None:
    wavelengths = np.linspace(400.0, 800.0, 201)
    rawResponse = 0.02 * wavelengths + 3.0
    rawResponse[50] = -1.0
    rawResponse[150] *= 4.0

    coefficients, fittedWavelengths, fittedResponse = _fit_response_polynomial(
        wavelengths,
        rawResponse,
        polynomialOrder=1,
        maxIterations=5,
    )

    np.testing.assert_allclose(coefficients, [0.02, 3.0], rtol=1e-12, atol=1e-12)
    assert len(fittedWavelengths) == 199
    assert wavelengths[50] not in fittedWavelengths
    assert wavelengths[150] not in fittedWavelengths
    np.testing.assert_allclose(
        fittedResponse,
        np.polyval(coefficients, fittedWavelengths),
        rtol=1e-12,
        atol=1e-12,
    )


def test_response_fit_matches_the_legacy_multi_iteration_algorithm() -> None:
    wavelengths = np.linspace(400.0, 800.0, 301)
    rawResponse = 0.02 * wavelengths + 3.0
    rawResponse[[50, 150]] = [-1.0, 80.0]

    expected = _legacy_fit_response_polynomial(
        wavelengths,
        rawResponse,
        polynomialOrder=1,
        maxIterations=5,
    )
    actual = _fit_response_polynomial(
        wavelengths,
        rawResponse,
        polynomialOrder=1,
        maxIterations=5,
    )

    np.testing.assert_allclose(actual[0], expected[0], rtol=0, atol=0)
    np.testing.assert_array_equal(actual[1], expected[1])
    np.testing.assert_array_equal(actual[2], expected[2])


def test_response_fit_excludes_requested_wavelength_regions() -> None:
    wavelengths = np.linspace(700.0, 1000.0, 301)
    rawResponse = 0.01 * wavelengths + 1.0

    coefficients, fittedWavelengths, _ = _fit_response_polynomial(
        wavelengths,
        rawResponse,
        polynomialOrder=1,
        maxIterations=2,
        excludedRegions=[(770.0, 850.0), (930.0, 970.0)],
    )

    assert not np.any((fittedWavelengths >= 770.0) & (fittedWavelengths <= 850.0))
    assert not np.any((fittedWavelengths >= 930.0) & (fittedWavelengths <= 970.0))
    np.testing.assert_allclose(coefficients, [0.01, 1.0], rtol=1e-12, atol=1e-12)


def test_response_fit_smooths_a_constant_nir_response() -> None:
    wavelengths = np.linspace(900.0, 1100.0, 201)
    rawResponse = np.full(wavelengths.shape, 7.5)

    coefficients, fittedWavelengths, fittedResponse = _fit_response_polynomial(
        wavelengths,
        rawResponse,
        polynomialOrder=0,
        maxIterations=2,
        smoothingSigma=5,
    )

    np.testing.assert_allclose(coefficients, [7.5], rtol=1e-14, atol=1e-14)
    np.testing.assert_array_equal(fittedWavelengths, wavelengths)
    np.testing.assert_allclose(fittedResponse, rawResponse, rtol=1e-14, atol=1e-14)


def test_response_fit_preserves_native_preprocessing_and_zero_iteration_errors() -> (
    None
):
    wavelengths = np.array([500.0, 600.0])
    rawResponse = np.array([5.0, 6.0])

    with pytest.raises(ValueError, match="too many values to unpack"):
        _fit_response_polynomial(
            wavelengths,
            rawResponse,
            polynomialOrder=1,
            maxIterations=1,
            excludedRegions=[(400.0, 450.0, 500.0)],
        )

    with pytest.raises(UnboundLocalError):
        _fit_response_polynomial(
            wavelengths,
            rawResponse,
            polynomialOrder=1,
            maxIterations=0,
        )


def test_response_fit_marks_only_iterative_convergence_failures() -> None:
    with pytest.raises(_ResponseFitConvergenceError) as excInfo:
        _fit_response_polynomial(
            np.array([]),
            np.array([]),
            polynomialOrder=1,
            maxIterations=1,
        )

    assert isinstance(excInfo.value.originalError, TypeError)


def test_response_get_rejects_an_unknown_standard_star(log: object) -> None:
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame({"WAVE": [500.0], "FLUX_COUNTS": [10.0]})
    response.stdExtractionNotFlatDF = pd.DataFrame(
        {"WAVE": [500.0], "FLUX_DENSITY_COUNTS": [10.0]}
    )
    response.stdAbsFluxDF = pd.DataFrame({"WAVE": [500.0], "KNOWN_STAR": [1.0]})
    response.std_objName = "UNKNOWN_STAR"

    with pytest.raises(
        LookupError,
        match="Standard star UNKNOWN_STAR not found.*KNOWN_STAR",
    ):
        response.get()


def test_response_get_writes_response_and_efficiency_products(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wavelengths = np.linspace(500.0, 520.0, 21)
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame(
        {"WAVE": wavelengths, "FLUX_COUNTS": np.full(21, 10.0)}
    )
    response.stdExtractionNotFlatDF = pd.DataFrame(
        {"WAVE": wavelengths, "FLUX_DENSITY_COUNTS": np.full(21, 10.0)}
    )
    response.stdAbsFluxDF = pd.DataFrame(
        {"WAVE": wavelengths, "SYNTHETIC_STAR": np.full(21, 2.0)}
    )
    response.std_objName = "SYNTHETIC_STAR"
    response.arm = "VIS"
    response.instrument = "soxs"
    response.texp = 10.0
    response.recipeSettings = {"vis": {"poly_order": 0, "max_iteration": 2}}
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
    captured: dict[str, object] = {}
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
    monkeypatch.setattr(
        phase3Module,
        "write_fits_table_to_disk",
        lambda **kwargs: captured.update({"efficiencyWrite": kwargs}),
    )
    response.write_response_function_to_file = lambda **kwargs: captured.update(
        {"responseWrite": kwargs}
    )
    response.plot_response_curve = lambda **kwargs: captured.update({"plot": kwargs})

    qc, products, recipeError = response.get()

    assert recipeError is False
    assert qc.empty
    assert products["product_label"].tolist() == ["EFFICIENCY"]
    np.testing.assert_allclose(
        captured["responseWrite"]["responseFuncCoeffs"],
        [2.0e17],
        rtol=1e-12,
    )
    assert captured["responseWrite"]["polyOrder"] == 0
    assert captured["efficiencyWrite"]["filePath"] == "/products/SYNTHETIC_EFFICIENCY.fits"
    assert captured["plot"]["stdEfficiencyEstimate"].shape == (21,)
