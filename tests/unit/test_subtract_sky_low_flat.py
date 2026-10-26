"""Low-flat pixels no longer dominate the sky fit or the noisy-region test of VIS stare (DY-1517).

At each end of a VIS order the master flat falls to about 0.001. After flat division those pixels are garbage with
a huge uncertainty, and a model of the noise in raw electrons ignores that.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from soxspipe.commonutils.subtract_sky import LOW_FLAT_FRACTION, subtract_sky
from tests.unit.test_subtract_sky_bspline_characterization import _subtractor
from tests.unit.test_subtract_sky_inverse_noise_weights import RON, _production_like_subtractor

pytestmark = pytest.mark.unit

PIXEL_COUNT = 6000
SKY = 100.0
# THE FLAT FALLS LINEARLY FROM 1 TO THIS VALUE OVER THE LAST FRACTION OF THE ORDER
FLAT_FRACTION_OF_ORDER = 0.2
NOISE_UNIT_END_FLAT = 0.05
# THE RESIDUAL OF EVERY PIXEL IS THIS MULTIPLE OF ITS OWN NOISE: A SMALL SYSTEMATIC OFFSET THAT AN ABSOLUTE
# THRESHOLD IN ELECTRONS READS AS A NOISY REGION WHERE THE FLAT IS LOW
NOISE_UNIT_OFFSET = 0.3
NOISE_SIGMA = 5
# A STRETCH OF THE FULL-FLAT REGION WHOSE RESIDUAL IS FAR ABOVE ITS NOISE, AND THE PART OF IT THAT THE ROLLING
# WINDOW OF 500 PIXELS FILLS ENTIRELY
NOISY_STRETCH = slice(1000, 2500)
NOISY_STRETCH_MIDDLE = slice(1300, 2200)
NOISY_STRETCH_SIGMA = 20.0


def _flat_profile(pixelCount: int, endFlat: float, fractionOfOrder: float) -> np.ndarray:
    """A flat of 1 that falls linearly to ``endFlat`` over the last ``fractionOfOrder`` of the order."""
    flat = np.ones(pixelCount)
    rampStart = int(pixelCount * (1.0 - fractionOfOrder))
    flat[rampStart:] = np.linspace(1.0, endFlat, pixelCount - rampStart)
    return flat


def _noise_unit_order() -> pd.DataFrame:
    """A constant sky whose flat-fielded residuals are a fixed multiple of their own noise."""
    rng = np.random.default_rng(3)
    flat = _flat_profile(PIXEL_COUNT, NOISE_UNIT_END_FLAT, FLAT_FRACTION_OF_ORDER)
    sigma = np.sqrt(SKY * flat + RON**2) / flat
    residual = (NOISE_UNIT_OFFSET + 0.05 * rng.normal(size=PIXEL_COUNT)) * sigma
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": np.linspace(500.0, 510.0, PIXEL_COUNT),
            "slit_position": np.zeros(PIXEL_COUNT),
            "flat": flat,
            "flux": SKY + residual,
            "error": sigma,
            "flux_minus_smoothed_residual": residual,
            "flux_percentile_smoothed": np.full(PIXEL_COUNT, SKY),
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "sky_residual_rolling_average": np.zeros(PIXEL_COUNT),
        }
    )


def _noise_unit_subtractor(log: Any) -> subtract_sky:
    subtractor = _subtractor(log, noiseSigma=NOISE_SIGMA)
    subtractor.recipeSettings["sky-subtraction"].update(
        {"percentile_rolling_window_size": 35, "noise_rolling_window_size": 500}
    )
    return subtractor


def test_the_long_residual_median_is_in_noise_units(log: Any) -> None:
    # ARRANGE
    subtractor = _noise_unit_subtractor(log)

    # ACT
    [clipped] = subtractor.clip_object_slit_positions([_noise_unit_order()])

    # ASSERT: EVERY PIXEL SITS AT ABOUT 0.3 NOISE, WHATEVER ITS FLAT
    longMedian = clipped["residual_windowed_long_median"].dropna()
    np.testing.assert_allclose(longMedian, NOISE_UNIT_OFFSET, atol=0.05)


def test_a_low_flat_end_with_proportionally_larger_noise_is_not_flagged_as_a_noisy_region(log: Any) -> None:
    # ARRANGE
    subtractor = _noise_unit_subtractor(log)
    [pixels] = subtractor.clip_object_slit_positions([_noise_unit_order()])
    spline = scipy.interpolate.splrep(pixels["wavelength"], np.full(PIXEL_COUNT, SKY), k=1)

    # ACT
    flagged, _ = subtractor.determine_residual_floor(pixels, spline, iteration=1)

    # ASSERT
    assert not flagged["flagged_noisy_region"].any()


def test_the_long_residual_median_without_a_flat_uses_the_sky_and_read_noise(log: Any) -> None:
    # ARRANGE
    subtractor = _noise_unit_subtractor(log)
    pixels = _noise_unit_order().drop(columns="flat")
    pixels["flux_minus_smoothed_residual"] = 2.0 * np.sqrt(SKY + RON**2)

    # ACT
    [clipped] = subtractor.clip_object_slit_positions([pixels])

    # ASSERT
    np.testing.assert_allclose(clipped["residual_windowed_long_median"].dropna(), 2.0)


def test_a_region_whose_residual_is_many_times_its_noise_is_still_flagged_as_a_noisy_region(log: Any) -> None:
    # ARRANGE: A STRETCH OF THE FULL-FLAT REGION, LONGER THAN THE ROLLING WINDOW, SITS AT 20 TIMES ITS NOISE
    subtractor = _noise_unit_subtractor(log)
    order = _noise_unit_order()
    isNoisy = np.zeros(PIXEL_COUNT, dtype=bool)
    isNoisy[NOISY_STRETCH] = True
    order["flux_minus_smoothed_residual"] = np.where(
        isNoisy, NOISY_STRETCH_SIGMA * order["error"], order["flux_minus_smoothed_residual"]
    )
    [pixels] = subtractor.clip_object_slit_positions([order])
    spline = scipy.interpolate.splrep(pixels["wavelength"], np.full(PIXEL_COUNT, SKY), k=1)

    # ACT
    flagged, _ = subtractor.determine_residual_floor(pixels, spline, iteration=1)

    # ASSERT: THE MIDDLE OF THE STRETCH IS FLAGGED AND THE LOW-FLAT END IS NOT
    isFlagged = flagged["flagged_noisy_region"].to_numpy()
    assert isFlagged[NOISY_STRETCH_MIDDLE].all()
    assert not isFlagged[int(PIXEL_COUNT * (1.0 - FLAT_FRACTION_OF_ORDER)) :].any()


# AN ORDER OF POISSON SKY WHOSE LAST PIXELS ARE OUT OF BAND: A FLAT OF 0.001 AND A SMALL NEGATIVE BIAS IN
# ELECTRONS THAT THE FLAT DIVISION TURNS INTO -1000
OUT_OF_BAND_FRACTION = 0.06
OUT_OF_BAND_FLAT = 0.001
OUT_OF_BAND_BIAS_E = -1.0
# THE FLAT FALLS FROM 1 TO ITS OUT-OF-BAND VALUE OVER THIS MANY PIXELS
EDGE_RAMP_PIXELS = 30
GOOD_REGION_END_NM = 509.0
MAX_GOOD_REGION_ERROR_SIGMA = 0.06


def _band_edge_order(seed: int) -> pd.DataFrame:
    """A sloped sky of 100 e- whose last 6% is out of band: the flat is 0.001 and the raw frame is biased low."""
    rng = np.random.default_rng(seed)
    wavelength = np.linspace(500.0, 510.0, PIXEL_COUNT)
    flat = np.ones(PIXEL_COUNT)
    bandEdge = int(PIXEL_COUNT * (1.0 - OUT_OF_BAND_FRACTION))
    flat[bandEdge - EDGE_RAMP_PIXELS : bandEdge] = np.linspace(1.0, OUT_OF_BAND_FLAT, EDGE_RAMP_PIXELS)
    flat[bandEdge:] = OUT_OF_BAND_FLAT
    trueSky = SKY + 5.0 * (wavelength - 505.0)
    isOutOfBand = flat < 0.5
    raw = rng.poisson(trueSky * flat) + rng.normal(0.0, RON, PIXEL_COUNT) + np.where(isOutOfBand, OUT_OF_BAND_BIAS_E, 0)
    flux = raw / flat
    return pd.DataFrame(
        {
            "order": 10,
            "x": np.arange(PIXEL_COUNT),
            "y": 0,
            "wavelength": wavelength,
            "slit_position": rng.uniform(-5.0, 5.0, PIXEL_COUNT),
            "flat": flat,
            "true_sky": trueSky,
            "flux": flux,
            "error": np.sqrt(np.clip(raw, 0.0, None) + RON**2) / flat,
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "flagged_bspline_clipped": False,
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
        }
    )


def _fit_order(log: Any, pixels: pd.DataFrame) -> pd.DataFrame:
    subtractor = _production_like_subtractor(log)
    subtractor.objectFrame = SimpleNamespace(data=np.zeros((1, PIXEL_COUNT)))
    lowFlat = subtractor._low_flat_pixels(pixels)
    pixels = subtractor._drop_excluded_pixels(pixels, lowFlat).copy()
    pixels = subtractor.rolling_window_clipping(pixels, windowSize=35, sigma_clip_limit=3, max_iterations=15)
    pixels = subtractor.clip_object_slit_positions([pixels], aggressive=False)[0]
    modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(pixels)
    return modelled


def test_out_of_band_order_ends_do_not_disturb_the_sky_model_of_the_good_region(log: Any) -> None:
    # ARRANGE
    errors = []

    # ACT
    for seed in range(3):
        modelled = _fit_order(log, _band_edge_order(seed))
        isGood = modelled["wavelength"] < GOOD_REGION_END_NM
        errors.append(((modelled["sky_model"] - modelled["true_sky"]) / np.sqrt(modelled["true_sky"] + RON**2))[isGood])

    # ASSERT: THE MEAN ERROR IS ABOUT 0.04 SIGMA WITH THE FIX AND ABOUT 0.085 SIGMA WITHOUT IT
    assert np.abs(np.concatenate(errors)).mean() < MAX_GOOD_REGION_ERROR_SIGMA


def test_the_out_of_band_pixels_are_left_out_of_the_fit(log: Any) -> None:
    # ACT
    modelled = _fit_order(log, _band_edge_order(0))

    # ASSERT
    assert modelled["flat"].min() >= LOW_FLAT_FRACTION
    assert len(modelled) < PIXEL_COUNT
