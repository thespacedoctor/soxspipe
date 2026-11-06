"""The VIS pre-fit clip must not walk down into the noise of a sky it cannot find an object in (DY-2081)."""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

subtract_sky_module = importlib.import_module("soxspipe.commonutils.subtract_sky")

PIXEL_COUNT = 60000
SKY_LEVEL = 400.0
WINDOW_SIZE = 35
START_SIGMA = 3.0
# THE PIXEL NOISE VARIES BY THIS FACTOR EITHER WAY AROUND THE SKY NOISE, AS THE FLAT-FIELD AND THE SLIT POSITION VARY IT
# IN A REAL ORDER. A WAVELENGTH WINDOW THEN MIXES PIXELS OF UNLIKE NOISE
NOISE_SPREAD = 4.0
# GAUSSIAN NOISE ABOVE THE 30TH PERCENTILE PLUS 3 STANDARD DEVIATIONS IS 0.66%, AND THE UNLIKE NOISE ABOUT DOUBLES THAT
# FIRST PASS. THE CLIP MUST STAY NEAR THAT, NOT WALK TO THE 5% THE RETRY ASKS FOR
MAX_NOISE_ONLY_SHARE = 0.03
COSMIC_RAY_SHARE = 0.001


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.ron = 3.3
    subtractor.debug = False
    return subtractor


def _sky_order(
    *, seed: int, noiseSpread: float = NOISE_SPREAD, objectShare: float = 0.0, objectSigma: float = 0.0
) -> pd.DataFrame:
    """Gaussian sky noise with cosmic rays, and optionally an object lifting a random `objectShare` of the pixels.

    A wavelength window of a real order mixes pixels of every slit position, so the object pixels are interspersed.
    """
    rng = np.random.default_rng(seed)
    noise = np.sqrt(SKY_LEVEL)
    pixelNoise = noise * np.exp(rng.uniform(-np.log(noiseSpread), np.log(noiseSpread), PIXEL_COUNT))
    flux = SKY_LEVEL + rng.normal(0.0, 1.0, PIXEL_COUNT) * pixelNoise
    isCosmic = rng.random(PIXEL_COUNT) < COSMIC_RAY_SHARE
    flux[isCosmic] += 50.0 * noise
    isObject = rng.random(PIXEL_COUNT) < objectShare
    flux[isObject] += objectSigma * noise
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": 500.0 + 0.01 * np.arange(PIXEL_COUNT),
            "flux": flux,
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
        }
    )


def _final_sigma(clipped: pd.DataFrame) -> float:
    """Return the sigma limit of the last clipping pass, recovered from the stored upper limit and std."""
    ratio = clipped["flux_minus_smoothed_residual_upper_limit"] / clipped["flux_minus_smoothed_residual_std"]
    return float(ratio.dropna().iloc[0])


def _clip(subtractor: subtract_sky, pixels: pd.DataFrame) -> pd.DataFrame:
    return subtractor.rolling_window_clipping(
        pixels, windowSize=WINDOW_SIZE, sigma_clip_limit=START_SIGMA, max_iterations=15
    )


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_a_vis_order_of_pure_sky_noise_clips_only_the_noise_tail_and_the_cosmic_rays(log: Any, seed: int) -> None:
    """Nothing in the order is an object, so the clip must not lower its limit until 5% of the pixels are gone."""
    subtractor = _subtractor(log)

    clipped = _clip(subtractor, _sky_order(seed=seed))

    assert clipped["flagged_object_clipped"].mean() < MAX_NOISE_ONLY_SHARE
    assert _final_sigma(clipped) == pytest.approx(START_SIGMA)


def test_a_vis_order_an_object_hides_from_the_first_pass_still_lowers_the_limit(log: Any) -> None:
    """Object light over 30% of the pixels inflates the rolling scatter, so the first pass clips next to nothing."""
    subtractor = _subtractor(log)
    pixels = _sky_order(seed=2, noiseSpread=1.0, objectShare=0.3, objectSigma=20.0)

    clipped = _clip(subtractor, pixels)

    assert _final_sigma(clipped) < START_SIGMA


def test_a_clip_at_the_noise_level_below_the_sigma_floor_neither_retries_nor_warns(log: Any) -> None:
    """Under 5% of all pixels are clipped, but that is the noise share of the few judged pixels: no retry or warning."""
    subtractor = _subtractor(log)
    pixels = _sky_order(seed=7, noiseSpread=1.0)
    isEdge = np.random.default_rng(7).random(PIXEL_COUNT) < 0.95
    pixels.loc[isEdge, "flagged_edge_clipped"] = True
    pixels.loc[isEdge, "flagged_all_clipped"] = True
    belowFloorSigma = 1.4

    clipped = subtractor.rolling_window_clipping(
        pixels, windowSize=WINDOW_SIZE, sigma_clip_limit=belowFloorSigma, max_iterations=15
    )

    assert _final_sigma(clipped) == pytest.approx(belowFloorSigma)
    assert [message for message in log.messages if message[0] == "warning"] == []


@pytest.mark.parametrize(
    ("sigmaLimit", "quantile", "expected"),
    [(3.0, 0.3, 0.0066), (2.7, 0.15, 0.0485), (1.5, 0.1, 0.4)],
)
def test_the_gaussian_clip_share_is_the_tail_above_the_percentile_plus_the_sigma_limit(
    sigmaLimit: float, quantile: float, expected: float
) -> None:
    """A Gaussian pixel is clipped when it lies above the `quantile` percentile plus `sigmaLimit` deviations."""
    share = subtract_sky_module._gaussian_clip_share(sigmaLimit, quantile)

    assert share == pytest.approx(expected, rel=0.05)
