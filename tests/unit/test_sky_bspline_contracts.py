"""Analytic B-spline sky-model contracts."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit


def _subtractor(log: Any, *, arm: str) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.bspline_order = 1
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_fitting_residual_clipping_sigma": 3,
            "bspline_iteration_limit": 1,
            "min_points_per_knot": 4,
            "noise_sigma": 0,
            "residual_floor_percentile": 50,
            "starting_points_per_knot": 100,
        }
    }
    return subtractor


def _constant_sky_pixels(*, include_nan_flux: bool = False) -> pd.DataFrame:
    wavelength = np.linspace(500.0, 510.0, 64)
    flux = np.full(wavelength.size, 12.0)
    if include_nan_flux:
        flux[-1] = np.nan
    residualStd = np.ones(wavelength.size)
    residualStd[1] = np.nan
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": wavelength,
            "flux": flux,
            "error": np.ones(wavelength.size),
            "residual_windowed_std": residualStd,
            "flagged_all_clipped": False,
            "flagged_noisy_region": False,
            "residual_windowed_long_median": np.zeros(wavelength.size),
            "flux_windowed_long_median": np.full(wavelength.size, 12.0),
        }
    )


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_fit_bspline_curve_models_constant_sky_for_both_weighting_paths(
    log: Any,
    arm: str,
) -> None:
    """A flat analytic sky remains flat after the public fitting workflow."""
    subtractor = _subtractor(log, arm=arm)

    modelled, spline, knots, fluxErrorRatio, residualFloor = (
        subtractor.fit_bspline_curve_to_sky(_constant_sky_pixels())
    )

    assert spline[2] == 1
    assert knots.size == 0
    assert residualFloor == 5
    np.testing.assert_allclose(modelled["sky_model"], 12.0)
    np.testing.assert_allclose(modelled["sky_subtracted_flux"], 0.0)
    np.testing.assert_allclose(modelled["residual_windowed_std"], 1.0)
    np.testing.assert_allclose(fluxErrorRatio, 1.0)


def test_fit_bspline_curve_excludes_nan_flux_from_model_quality_metrics(
    log: Any,
) -> None:
    """Invalid flux is clipped while the remaining analytic sky is modelled."""
    subtractor = _subtractor(log, arm="VIS")

    modelled, _, _, fluxErrorRatio, _ = subtractor.fit_bspline_curve_to_sky(
        _constant_sky_pixels(include_nan_flux=True)
    )

    assert modelled.iloc[-1]["flagged_all_clipped"]
    assert modelled.iloc[-1]["sky_model"] == pytest.approx(12.0)
    assert np.isnan(modelled.iloc[-1]["sky_subtracted_flux"])
    assert fluxErrorRatio.shape == (63,)
    np.testing.assert_allclose(fluxErrorRatio, 1.0)


def test_cross_dispersion_normaliser_retains_flux_and_adds_unity_ratio(
    log: Any,
) -> None:
    """Slit-illumination fitting leaves the legacy unity correction in place."""
    subtractor = _subtractor(log, arm="VIS")
    subtractor.recipeSettings["sky-subtraction"]["slit_illumination_order"] = 1
    slitPositions = np.linspace(-1.0, 1.0, 32)
    pixels = pd.DataFrame(
        {
            "order": 10,
            "slit_position": slitPositions,
            "sky_subtracted_flux": 2.0 + slitPositions,
            "residual_windowed_std": np.ones(slitPositions.size),
            "flagged_all_clipped": False,
        }
    )

    corrected = subtractor.cross_dispersion_flux_normaliser(pixels.copy())

    np.testing.assert_allclose(corrected["sky_subtracted_flux"], 2.0 + slitPositions)
    np.testing.assert_allclose(corrected["slit_normalisation_ratio"], 1.0)


def test_adjust_tilt_returns_wavelength_sorted_corrected_order(log: Any) -> None:
    """Tilt optimization preserves all pixels and leaves an ordered wavelength axis."""
    subtractor = _subtractor(log, arm="VIS")
    subtractor.recipeSettings["sky-subtraction"]["slit_illumination_order"] = 1
    subtractor.qcPlotOrder = 10
    wavelengths = np.linspace(500.0, 510.0, 80)
    slitPositions = np.linspace(-1.0, 1.0, 80)
    spline = splrep(wavelengths, 0.5 * wavelengths + 1.0, k=1)
    pixels = pd.DataFrame(
        {
            "order": 10,
            "wavelength": wavelengths[::-1],
            "slit_position": slitPositions,
            "pixelScale": np.ones(wavelengths.size),
            "flux": (0.5 * wavelengths[::-1]) + 1.0,
            "residual_windowed_std": np.ones(wavelengths.size),
            "flagged_all_clipped": False,
        }
    )

    corrected = subtractor.adjust_tilt(pixels.copy(), spline)

    assert len(corrected) == len(pixels)
    assert corrected["wavelength"].is_monotonic_increasing
    assert np.isfinite(corrected["wavelength"]).all()
