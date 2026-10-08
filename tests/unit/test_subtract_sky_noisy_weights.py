"""Fit weights of the noisy-region pixels in the B-spline sky fit (DY-1253)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit.test_subtract_sky_bspline_characterization import _subtractor

pytestmark = pytest.mark.unit

PIXEL_COUNT = 3000
TRUE_SKY = 3.0
NOISE = 4.0
NOISY_BLOCK = slice(1000, 2000)
SEEDS = range(8)
# THE BLOCK MEAN IS NOISY BECAUSE THE BLOCK IS DOWN-WEIGHTED BY DESIGN, SO IT IS AVERAGED OVER SEEDS. THE
# UNBOUNDED REWEIGHT GAVE ~0.5 AND THE BOUNDED ONE GIVES ~3.3 AGAINST A TRUE SKY OF 3
MIN_RECOVERED_SKY = 2.0


def _low_flux_order(seed: int) -> pd.DataFrame:
    """A flat 3 e- sky with 4 e- noise whose middle third is flagged as a noisy region."""
    rng = np.random.default_rng(seed)
    longMedian = np.zeros(PIXEL_COUNT)
    longMedian[NOISY_BLOCK] = 50.0
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": np.linspace(500.0, 510.0, PIXEL_COUNT),
            "flux": TRUE_SKY + rng.normal(0.0, NOISE, PIXEL_COUNT),
            "error": np.full(PIXEL_COUNT, NOISE),
            "residual_windowed_std": np.full(PIXEL_COUNT, NOISE),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": longMedian,
            "flux_windowed_long_median": np.full(PIXEL_COUNT, TRUE_SKY),
        }
    )


def test_noisy_low_flux_region_is_not_pulled_toward_zero(log: Any) -> None:
    """Noisy-region weights stay bounded where flux crosses zero, so the sky keeps its level (DY-1253)."""
    # ARRANGE
    blockMeans = []

    # ACT
    for seed in SEEDS:
        subtractor = _subtractor(log, iterationLimit=8, pointsPerKnot=1000, noiseSigma=3)
        subtractor.recipeSettings["sky-subtraction"]["min_points_per_knot"] = 155
        modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(_low_flux_order(seed))
        assert int(modelled["flagged_noisy_region"].sum()) == 1000
        blockMeans.append(modelled["sky_model"].iloc[NOISY_BLOCK].mean())

    # ASSERT
    assert np.mean(blockMeans) > MIN_RECOVERED_SKY


def _scale(flux: list[float], error: list[float], windowedStd: list[float]) -> np.ndarray:
    pixels = pd.DataFrame({"flux": flux, "error": error, "residual_windowed_std": windowedStd})
    return subtract_sky._noisy_region_flux_scale(pixels)


def test_flux_scale_is_the_absolute_flux_when_flux_exceeds_the_noise() -> None:
    """Bright noisy pixels keep the 1/abs(flux) down-weighting of DY-695."""
    scale = _scale([100.0, -80.0], [4.0, 4.0], [5.0, 5.0])

    assert scale.tolist() == [100.0, 80.0]


def test_flux_scale_is_the_pixel_error_when_flux_is_below_the_noise() -> None:
    """Flux near or at zero gives a finite scale equal to the pixel noise."""
    scale = _scale([0.0, -1.0, 2.0], [4.0, 4.0, 4.0], [9.0, 9.0, 9.0])

    assert scale.tolist() == [4.0, 4.0, 4.0]


@pytest.mark.parametrize("badError", [0.0, np.nan, np.inf, -1.0])
def test_flux_scale_falls_back_to_the_windowed_std_when_the_error_is_unusable(badError: float) -> None:
    """A zero or non-finite error (no dark subtraction) must not give a zero or infinite scale."""
    scale = _scale([0.0], [badError], [6.0])

    assert scale.tolist() == [6.0]


def test_residual_floor_leaves_infinite_values_in_other_columns_alone(log: Any) -> None:
    """Only sky_residuals has its infinities reset to 1; an infinite fit weight is not rewritten (DY-1253)."""
    # ARRANGE
    subtractor = _subtractor(log)
    wavelength = np.linspace(500.0, 501.0, 32)
    pixels = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": np.full(32, 10.0),
            "error": np.ones(32),
            "weights": np.ones(32),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": np.zeros(32),
            "flux_windowed_long_median": np.full(32, 10.0),
            "sky_residual_rolling_average": np.full(32, 2.0),
        }
    )
    pixels.loc[5, "error"] = 0.0
    pixels.loc[7, "weights"] = np.inf
    spline = splrep(wavelength, np.full(32, 8.0), k=1)

    # ACT
    result, _ = subtractor.determine_residual_floor(pixels.copy(), spline, iteration=1)

    # ASSERT
    assert result.loc[7, "weights"] == np.inf
    assert result.loc[5, "sky_residuals"] == 1.0
