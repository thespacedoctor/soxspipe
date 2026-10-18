"""Analytic contracts for extracted-spectrum flux calibration."""

from __future__ import annotations

import warnings

import numpy as np
import pytest

from soxspipe.commonutils.flux_calibration import _calculate_flux_calibration, _response_factors_from_table

pytestmark = pytest.mark.unit


def test_flux_calibration_applies_exposure_extinction_and_response() -> None:
    counts = np.array([20.0, 40.0, 60.0])
    extinctionFactors = np.array([1.0, 1.1, 1.2])
    responseFactors = np.array([6.0, 7.0, 8.0])

    result = _calculate_flux_calibration(
        counts,
        exposureTime=10.0,
        responseFactors=responseFactors,
        extinctionFactors=extinctionFactors,
    )

    expected = counts / 10.0 * extinctionFactors * responseFactors * 1e-17
    np.testing.assert_allclose(result, expected, rtol=1e-14, atol=0)


def test_flux_calibration_preserves_nan_samples() -> None:
    result = _calculate_flux_calibration(
        np.array([10.0, np.nan]),
        exposureTime=2.0,
        responseFactors=np.array([1.0, 1.0]),
        extinctionFactors=1.0,
    )

    np.testing.assert_allclose(result[0], 5e-17, rtol=1e-14, atol=0)
    assert np.isnan(result[1])


def test_response_factors_interpolate_the_table_in_log_space() -> None:
    # ARRANGE
    tableWave = np.array([500.0, 600.0])
    tableResponse = np.array([100.0, 400.0])

    # ACT
    factors = _response_factors_from_table(np.array([500.0, 550.0, 600.0]), tableWave, tableResponse)

    # ASSERT
    np.testing.assert_allclose(factors, [100.0, 200.0, 400.0], rtol=1e-12)


def test_response_factors_are_nan_outside_the_table_range() -> None:
    # ACT
    factors = _response_factors_from_table(
        np.array([450.0, 550.0, 650.0]), np.array([500.0, 600.0]), np.array([1.0, 1.0])
    )

    # ASSERT
    assert np.isnan(factors[0])
    assert factors[1] == pytest.approx(1.0)
    assert np.isnan(factors[2])


@pytest.mark.parametrize("badResponse", [0.0, -3.0, np.nan, np.inf])
def test_response_factors_interpolate_over_the_good_rows_when_a_table_row_is_unusable(badResponse: float) -> None:
    # ARRANGE
    tableWave = np.array([500.0, 550.0, 600.0])
    tableResponse = np.array([100.0, badResponse, 400.0])

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        factors = _response_factors_from_table(np.array([550.0]), tableWave, tableResponse)

    # ASSERT
    np.testing.assert_allclose(factors, [200.0], rtol=1e-12)


def test_response_factors_ignore_a_table_row_with_a_non_finite_wavelength() -> None:
    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        factors = _response_factors_from_table(
            np.array([550.0]), np.array([500.0, np.nan, 600.0]), np.array([100.0, 1.0, 400.0])
        )

    # ASSERT
    np.testing.assert_allclose(factors, [200.0], rtol=1e-12)


def test_response_factors_sort_an_unordered_table_by_wavelength() -> None:
    # ACT
    factors = _response_factors_from_table(np.array([550.0]), np.array([600.0, 500.0]), np.array([400.0, 100.0]))

    # ASSERT
    np.testing.assert_allclose(factors, [200.0], rtol=1e-12)


def test_response_factors_raise_when_no_table_row_is_usable() -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="no usable rows"):
        _response_factors_from_table(np.array([550.0]), np.array([500.0, 600.0]), np.array([0.0, np.nan]))
