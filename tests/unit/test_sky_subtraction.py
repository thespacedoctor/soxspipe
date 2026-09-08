"""Analytic residual contracts for sky subtraction."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit


def _subtractor(log: object) -> subtract_sky:
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "residual_floor_percentile": 50,
            "noise_sigma": 0,
        }
    }
    return subtractor


def test_calculate_residuals_evaluates_polynomial_and_propagates_nan(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    pixels = pd.DataFrame(
        {
            "order_pow_0": [1.0, 1.0, 1.0],
            "wavelength_pow_0": [1.0, 1.0, 1.0],
            "wavelength_pow_1": [2.0, 3.0, 4.0],
            "slit_position_pow_0": [1.0, 1.0, 1.0],
            "sky_subtracted_flux": [4.0, 8.0, np.nan],
        }
    )

    mean, std, median, result = subtractor.calculate_residuals(
        pixels.copy(),
        fluxcoeff=[1.0, 2.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
    )

    np.testing.assert_allclose(
        result["fit_sky_subtracted_flux"],
        [5.0, 7.0, 9.0],
        rtol=1e-12,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        result["residuals_sky_subtracted_flux"],
        [1.0, -1.0, np.nan],
        rtol=1e-12,
        atol=1e-12,
        equal_nan=True,
    )
    assert mean == pytest.approx(0.0)
    assert std == pytest.approx(1.0)
    assert np.isnan(median)


def test_residual_floor_normalises_errors_and_flags_local_excess(
    log: object,
) -> None:
    subtractor = _subtractor(log)
    wavelength = np.linspace(500.0, 501.0, 32)
    model = np.full(32, 8.0)
    pixels = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": np.full(32, 10.0),
            "error": np.ones(32),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": np.zeros(32),
            "flux_windowed_long_median": np.full(32, 10.0),
            "sky_residual_rolling_average": np.full(32, 2.0),
        }
    )
    pixels.loc[5, "error"] = 0.0
    pixels.loc[10, "sky_residual_rolling_average"] = 3.0
    spline = splrep(wavelength, model, k=1)

    result, residualFloor = subtractor.determine_residual_floor(
        pixels.copy(),
        spline,
        iteration=1,
    )

    assert residualFloor == 5
    assert result.loc[5, "sky_residuals"] == 1.0
    assert np.isfinite(result["sky_residuals"]).all()
    assert result.loc[10, "flagged_sky_line"] == "line"
    assert result.loc[0, "flagged_sky_line"] == False
    assert not result["flagged_noisy_region"].any()
