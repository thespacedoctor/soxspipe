"""Inverse-noise weights of the B-spline sky fit (DY-1282)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit.test_subtract_sky_bspline_characterization import _skyline_order, _subtractor

pytestmark = pytest.mark.unit

RON = 3.3
PIXEL_COUNT = 20000
CONTINUUM = 10.0
SKYLINES = [(501.3, 300.0), (502.0, 2000.0), (503.1, 80.0), (504.5, 1500.0), (506.2, 600.0), (507.3, 2500.0)]
LINE_WIDTH_NM = 0.03
# A SLIT-DEPENDENT WAVELENGTH ERROR OF UP TO ±0.01 NM, SO A 1D SKY MODEL CANNOT FIT EVERY PIXEL EXACTLY
# AND THE WEIGHTS DECIDE WHERE THE MISFIT GOES. WITHOUT IT, ANY WEIGHTS GIVE AN UNBIASED FIT
SLIT_WAVELENGTH_ERROR_NM = 0.01
# THE LINE-PIXEL MEAN OF ONE SEED SCATTERS BY ABOUT ±0.2σ, SO THE BIAS IS AVERAGED OVER SEEDS
SEEDS = range(6)
# THE OLD FLUX-BOOSTED VIS WEIGHTS GAVE +0.83σ OVER LINE PIXELS ON THIS SKY; INVERSE-NOISE WEIGHTS GIVE +0.02σ
MAX_LINE_BIAS_SIGMA = 0.1
LINE_PERCENTILE = 95
# SLIT POSITIONS AT WHICH THE TRUE SKY IS AVERAGED TO GIVE THE BEST SKY A 1D MODEL CAN REACH
SLIT_AVERAGE_SAMPLES = 41


def _production_like_subtractor(log: Any) -> subtract_sky:
    """The soxs-stare VIS defaults for every setting the clipping and the fit read."""
    subtractor = _subtractor(log, arm="VIS", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.ron = RON
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "bspline_fitting_residual_clipping_sigma": 10,
            "min_points_per_knot": 15,
            "noise_rolling_window_size": 2000,
            "percentile_rolling_window_size": 35,
            "residual_floor_percentile": 80,
        }
    )
    return subtractor


def _true_sky(wavelength: np.ndarray, wavelengthError: np.ndarray) -> np.ndarray:
    sky = CONTINUUM + 0.5 * (wavelength - 505.0)
    for centre, amplitude in SKYLINES:
        sky = sky + amplitude * np.exp(-0.5 * ((wavelength + wavelengthError - centre) / LINE_WIDTH_NM) ** 2)
    return sky


def _poisson_sky_order(seed: int) -> tuple[pd.DataFrame, np.ndarray]:
    """An oversampled order of Poisson sky plus read noise, and the slit-averaged true sky of every pixel.

    Judging the model against each pixel's own slit-shifted sky would be biased: the brightest
    pixels are those shifted onto a line peak, which no 1D model reaches.
    """
    rng = np.random.default_rng(seed)
    wavelength = np.sort(rng.uniform(500.0, 510.0, PIXEL_COUNT))
    slitPosition = rng.uniform(-5.0, 5.0, PIXEL_COUNT)
    sky = _true_sky(wavelength, SLIT_WAVELENGTH_ERROR_NM * slitPosition / 5.0)
    slitGrid = np.linspace(-1.0, 1.0, SLIT_AVERAGE_SAMPLES) * SLIT_WAVELENGTH_ERROR_NM
    slitAveragedSky = _true_sky(wavelength[:, None], slitGrid[None, :]).mean(axis=1)
    flux = rng.poisson(sky) + rng.normal(0.0, RON, PIXEL_COUNT)
    pixels = pd.DataFrame(
        {
            "order": 10,
            "wavelength": wavelength,
            "slit_position": slitPosition,
            "flux": flux,
            "error": np.sqrt(np.clip(flux, 0.0, None) + RON**2),
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "flagged_bspline_clipped": False,
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
        }
    )
    return pixels, slitAveragedSky


def _fit_through_the_production_clipping(subtractor: subtract_sky, pixels: pd.DataFrame) -> pd.DataFrame:
    """Build the percentile-sky and noise columns as production does, then fit every pixel.

    The object clip is undone before the fit: it rejects only high pixels (DY-1281),
    which would bias the line peaks low whatever the weights.
    """
    pixels = subtractor.rolling_window_clipping(pixels, windowSize=35, sigma_clip_limit=3, max_iterations=15)
    pixels = subtractor.clip_object_slit_positions([pixels], aggressive=False)[0]
    pixels["flagged_object_clipped"] = False
    pixels["flagged_all_clipped"] = False
    modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(pixels)
    return modelled


def test_the_mean_sky_model_over_skyline_pixels_is_unbiased(log: Any) -> None:
    """VIS once boosted the weights by the measured flux, which biased line peaks high (DY-1282)."""
    # ARRANGE
    lineBiases = []

    # ACT
    for seed in SEEDS:
        pixels, sky = _poisson_sky_order(seed)
        pixels["true_sky"] = sky
        modelled = _fit_through_the_production_clipping(_production_like_subtractor(log), pixels)
        normalisedError = (modelled["sky_model"] - modelled["true_sky"]) / np.sqrt(modelled["true_sky"] + RON**2)
        isLine = modelled["true_sky"] > np.percentile(modelled["true_sky"], LINE_PERCENTILE)
        lineBiases.append(normalisedError[isLine].mean())

    # ASSERT
    assert abs(np.mean(lineBiases)) < MAX_LINE_BIAS_SIGMA


def test_inverse_noise_weights_follow_the_model_sky_and_floor_negative_sky_at_the_read_noise() -> None:
    # ARRANGE
    skyModel = np.array([-5.0, 0.0, 16.0, np.nan, -np.inf])

    # ACT
    weights = subtract_sky._inverse_noise_weights(skyModel, 3.0)

    # ASSERT
    np.testing.assert_allclose(weights, [1 / 3.0, 1 / 3.0, 1 / 5.0, 1 / 3.0, 1 / 3.0])


def test_an_infinite_model_sky_gets_zero_weight() -> None:
    weights = subtract_sky._inverse_noise_weights(np.array([np.inf]), 3.0)

    assert weights.tolist() == [0.0]


def test_nan_values_are_filled_by_interpolation_and_an_all_nan_array_becomes_zeros() -> None:
    # ACT
    filled = subtract_sky._fill_nans_by_interpolation(np.array([np.nan, 2.0, np.nan, 6.0, np.nan]))
    allNan = subtract_sky._fill_nans_by_interpolation(np.full(3, np.nan))

    # ASSERT
    np.testing.assert_allclose(filled, [2.0, 2.0, 4.0, 6.0, 6.0])
    np.testing.assert_allclose(allNan, 0.0)


@pytest.mark.parametrize("headerValue", [None, "0", 0.0, -1.0, np.nan, np.inf, "junk"])
def test_an_unusable_header_read_noise_falls_back_to_the_detector_default(log: Any, headerValue: Any) -> None:
    ron = subtract_sky._read_noise(headerValue, 3.8, log)

    assert ron == 3.8


@pytest.mark.parametrize("headerValue", [3.3, "3.3"])
def test_a_usable_header_read_noise_is_used(log: Any, headerValue: Any) -> None:
    ron = subtract_sky._read_noise(headerValue, 3.8, log)

    assert ron == 3.3
    assert not [message for level, message in log.messages if level == "warning"]


def test_an_unusable_header_read_noise_is_logged(log: Any) -> None:
    subtract_sky._read_noise("junk", 3.8, log)

    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "junk" in warnings[0]


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_each_fit_is_weighted_by_the_noise_of_the_previous_model_not_the_measured_flux(
    log: Any, arm: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The first fit uses the rolling-percentile sky; every later fit uses the fit before it."""
    # ARRANGE
    realSplrep = scipy.interpolate.splrep
    calls: list[tuple[np.ndarray, Any]] = []

    def recording_splrep(x: Any, y: Any, **kwargs: Any) -> Any:
        result = realSplrep(x, y, **kwargs)
        calls.append((np.asarray(kwargs["w"], dtype=float).copy(), result[0]))
        return result

    monkeypatch.setattr(scipy.interpolate, "splrep", recording_splrep)
    subtractor = _subtractor(log, arm=arm)
    pixels = _skyline_order()
    percentileSky = pixels["flux_percentile_smoothed"].to_numpy()
    wavelength = pixels["wavelength"].to_numpy()

    # ACT
    subtractor.fit_bspline_curve_to_sky(pixels)

    # ASSERT
    assert len(calls) > 2
    expectedFirst = 1.0 / np.sqrt(np.clip(percentileSky, 0.0, None) + subtractor.ron**2)
    # THE FIRST AND LAST SAMPLES CARRY THE ORDER-END ANCHOR WEIGHT
    np.testing.assert_allclose(calls[0][0][1:-1], expectedFirst[1:-1])
    for (_, previousSpline), (weights, _) in zip(calls[:-1], calls[1:], strict=True):
        previousModel = scipy.interpolate.splev(wavelength, previousSpline)
        expected = 1.0 / np.sqrt(np.clip(previousModel, 0.0, None) + subtractor.ron**2)
        np.testing.assert_allclose(weights[1:-1], expected[1:-1])
