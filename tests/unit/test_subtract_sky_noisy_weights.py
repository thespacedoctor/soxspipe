"""A noisy low-flux region keeps its sky level in the B-spline sky fit (DY-1253).

The 1/abs(flux) noisy-region reweight this file once pinned was removed with the inverse-noise
weights of DY-1282; the region-level test stays as the regression guard.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

from tests.unit.test_subtract_sky_bspline_characterization import _subtractor

pytestmark = pytest.mark.unit

PIXEL_COUNT = 3000
TRUE_SKY = 3.0
NOISE = 4.0
NOISY_BLOCK = slice(1000, 2000)
SEEDS = range(8)
# THE BLOCK MEAN IS AVERAGED OVER SEEDS. THE UNBOUNDED 1/abs(flux) REWEIGHT OF DY-1253 GAVE ~0.5 AGAINST A TRUE SKY OF 3
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
            "flux_percentile_smoothed": np.full(PIXEL_COUNT, TRUE_SKY),
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
