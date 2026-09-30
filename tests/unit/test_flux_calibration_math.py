"""Analytic contracts for extracted-spectrum flux calibration."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.flux_calibration import _calculate_flux_calibration

pytestmark = pytest.mark.unit


def test_flux_calibration_applies_exposure_extinction_and_response() -> None:
    wavelengths = np.array([400.0, 500.0, 600.0])
    counts = np.array([20.0, 40.0, 60.0])
    extinctionFactors = np.array([1.0, 1.1, 1.2])
    responseCoefficients = np.array([0.01, 2.0])

    result = _calculate_flux_calibration(
        wavelengths,
        counts,
        exposureTime=10.0,
        responseCoefficients=responseCoefficients,
        extinctionFactors=extinctionFactors,
    )

    expected = (
        counts
        / 10.0
        * extinctionFactors
        * np.polyval(responseCoefficients, wavelengths)
        * 1e-17
    )
    np.testing.assert_allclose(result, expected, rtol=1e-14, atol=0)


def test_flux_calibration_preserves_nan_samples() -> None:
    result = _calculate_flux_calibration(
        np.array([500.0, 600.0]),
        np.array([10.0, np.nan]),
        exposureTime=2.0,
        responseCoefficients=np.array([1.0]),
        extinctionFactors=1.0,
    )

    np.testing.assert_allclose(result[0], 5e-17, rtol=1e-14, atol=0)
    assert np.isnan(result[1])
