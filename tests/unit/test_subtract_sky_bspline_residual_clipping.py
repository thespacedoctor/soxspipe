"""Two-sided residual clipping of the B-spline sky fit (DY-1281)."""

from __future__ import annotations

import importlib
from typing import Any

import numpy as np
import pandas as pd
import pytest
import scipy.interpolate

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit.test_subtract_sky_bspline_characterization import _subtractor

pytestmark = pytest.mark.unit

# THE PACKAGE RE-EXPORTS THE CLASS UNDER THE MODULE'S NAME, SO THE MODULE IS FETCHED BY ITS DOTTED PATH
subtract_sky_module = importlib.import_module("soxspipe.commonutils.subtract_sky")

RON = 3.3
PIXEL_COUNT = 20000
CONTINUUM = 40.0
SKYLINES = [(501.3, 300.0), (502.0, 2000.0), (503.1, 80.0), (504.5, 1500.0), (506.2, 600.0), (507.3, 2500.0)]
LINE_WIDTH_NM = 0.03
# THE SOXS DEFAULT FOR bspline_fitting_residual_clipping_sigma
CLIPPING_SIGMA = 10
OUTLIER_COUNT = 40
# A COSMIC-RAY RESIDUAL AND A COLD PIXEL, IN UNITS OF THE PIXEL'S NOISE; BOTH SIT WELL OUTSIDE THE CLIP
COSMIC_RAY_SIGMA = 200.0
COLD_PIXEL_SIGMA = -40.0
MAX_MODEL_DIFFERENCE_SIGMA = 0.1
MODEL_DIFFERENCE_PERCENTILE = 99
# THE HALF-WIDTH OF THE WAVELENGTH RANGE AROUND AN OUTLIER THAT IT WOULD PULL
NEIGHBOURHOOD_NM = 0.02
SEEDS = range(6)


def _production_like_subtractor(log: Any) -> subtract_sky:
    """The soxs-stare VIS defaults for every setting the fit reads."""
    subtractor = _subtractor(log, arm="VIS", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.ron = RON
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "bspline_fitting_residual_clipping_sigma": CLIPPING_SIGMA,
            "min_points_per_knot": 155,
            "residual_floor_percentile": 90,
        }
    )
    return subtractor


def _true_sky(wavelength: np.ndarray) -> np.ndarray:
    sky = CONTINUUM + 0.5 * (wavelength - 505.0)
    for centre, amplitude in SKYLINES:
        sky = sky + amplitude * np.exp(-0.5 * ((wavelength - centre) / LINE_WIDTH_NM) ** 2)
    return sky


def _sky_order(flux: np.ndarray, wavelength: np.ndarray, sky: np.ndarray) -> pd.DataFrame:
    """An order of sky pixels with the columns the earlier clipping steps leave for the fit."""
    noise = np.sqrt(sky + RON**2)
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": wavelength,
            "flux": flux,
            "error": noise,
            "residual_windowed_std": noise.copy(),
            "flux_percentile_smoothed": sky,
            "residual_windowed_long_median": np.zeros(PIXEL_COUNT),
            "flux_windowed_long_median": np.full(PIXEL_COUNT, CONTINUUM),
            "flagged_all_clipped": False,
            "flagged_bspline_clipped": False,
        }
    )


def _clean_and_contaminated_orders(seed: int) -> tuple[pd.DataFrame, pd.DataFrame, np.ndarray]:
    """The same Poisson sky with and without injected cosmic rays and cold pixels, and the injected rows."""
    rng = np.random.default_rng(seed)
    wavelength = np.sort(rng.uniform(500.0, 510.0, PIXEL_COUNT))
    sky = _true_sky(wavelength)
    cleanFlux = rng.poisson(sky) + rng.normal(0.0, RON, PIXEL_COUNT)
    # KEEP THE OUTLIERS AWAY FROM THE ORDER ENDS, WHOSE SAMPLES ARE REPLACED BY THE END ANCHORS
    outlierRows = rng.choice(np.arange(100, PIXEL_COUNT - 100), size=OUTLIER_COUNT, replace=False)
    noise = np.sqrt(sky + RON**2)
    shiftSigma = np.where(np.arange(OUTLIER_COUNT) % 2 == 0, COSMIC_RAY_SIGMA, COLD_PIXEL_SIGMA)
    contaminatedFlux = cleanFlux.copy()
    contaminatedFlux[outlierRows] += shiftSigma * noise[outlierRows]
    return (
        _sky_order(cleanFlux, wavelength, sky),
        _sky_order(contaminatedFlux, wavelength, sky),
        outlierRows,
    )


def test_cosmic_rays_and_cold_pixels_do_not_pull_the_sky_model(log: Any) -> None:
    """The fit to contaminated data matches the fit to the same data without the outliers (DY-1281).

    Without the clip, the outliers pull the model by 3 to 6σ. With it, a handful of pixels on bright-line
    flanks still differ by up to 0.25σ: knots sit at the mean wavelength of a group of pixels, so removing
    40 pixels moves a knot by about 0.001 nm. That is knot placement, not outlier pull, so the test reads
    the 99th percentile of every pixel and the mean around each outlier rather than the maximum.
    """
    # ARRANGE
    worstPercentile = []
    worstNearOutlierMean = []

    # ACT
    for seed in SEEDS:
        clean, contaminated, outlierRows = _clean_and_contaminated_orders(seed)
        wavelength = clean["wavelength"].to_numpy()
        noise = np.sqrt(_true_sky(wavelength) + RON**2)
        cleanModel, _, _, _, _ = _production_like_subtractor(log).fit_bspline_curve_to_sky(clean)
        contaminatedModel, _, _, _, _ = _production_like_subtractor(log).fit_bspline_curve_to_sky(contaminated)
        difference = (contaminatedModel["sky_model"].to_numpy() - cleanModel["sky_model"].to_numpy()) / noise
        isNearOutlier = (np.abs(wavelength[:, None] - wavelength[outlierRows][None, :]) < NEIGHBOURHOOD_NM).any(axis=1)
        worstPercentile.append(np.percentile(np.abs(difference), MODEL_DIFFERENCE_PERCENTILE))
        worstNearOutlierMean.append(abs(difference[isNearOutlier].mean()))

    # ASSERT
    assert max(worstPercentile) < MAX_MODEL_DIFFERENCE_SIGMA
    assert max(worstNearOutlierMean) < MAX_MODEL_DIFFERENCE_SIGMA


def test_injected_cosmic_rays_and_cold_pixels_are_flagged_as_bspline_clipped(log: Any) -> None:
    # ARRANGE
    _, contaminated, outlierRows = _clean_and_contaminated_orders(seed=3)
    contaminated["row"] = np.arange(PIXEL_COUNT)

    # ACT
    modelled, _, _, _, _ = _production_like_subtractor(log).fit_bspline_curve_to_sky(contaminated)

    # ASSERT
    clippedRows = set(modelled.loc[modelled["flagged_bspline_clipped"], "row"])
    allClippedRows = set(modelled.loc[modelled["flagged_all_clipped"], "row"])
    assert clippedRows == set(outlierRows.tolist())
    assert clippedRows <= allClippedRows


def test_a_clean_sky_with_bright_skylines_has_no_pixel_clipped(log: Any) -> None:
    """Early fits do not resolve the skylines; the clip must leave those lines to knot insertion."""
    # ARRANGE
    clean, _, _ = _clean_and_contaminated_orders(seed=3)

    # ACT
    modelled, _, _, _, _ = _production_like_subtractor(log).fit_bspline_curve_to_sky(clean)

    # ASSERT
    assert int(modelled["flagged_bspline_clipped"].sum()) == 0
    assert int(modelled["flagged_all_clipped"].sum()) == 0


def test_the_clip_limit_is_the_bspline_fitting_residual_clipping_sigma_setting(log: Any) -> None:
    """A limit above the injected outliers clips none of them."""
    # ARRANGE
    _, contaminated, _ = _clean_and_contaminated_orders(seed=3)
    subtractor = _production_like_subtractor(log)
    subtractor.recipeSettings["sky-subtraction"]["bspline_fitting_residual_clipping_sigma"] = 2 * COSMIC_RAY_SIGMA

    # ACT
    modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(contaminated)

    # ASSERT
    assert int(modelled["flagged_bspline_clipped"].sum()) == 0


def test_the_final_clip_and_refit_stops_after_the_maximum_number_of_passes(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A limit inside the noise clips new pixels on every pass, so only the pass limit stops the refits."""
    # ARRANGE
    realSplrep = scipy.interpolate.splrep
    callsPerLimit = {}

    def counting_splrep(*args: Any, **kwargs: Any) -> Any:
        calls.append(1)
        return realSplrep(*args, **kwargs)

    monkeypatch.setattr(scipy.interpolate, "splrep", counting_splrep)

    # ACT
    for passLimit in (2, 3):
        calls: list[int] = []
        monkeypatch.setattr(subtract_sky_module, "BSPLINE_CLIP_MAX_PASSES", passLimit)
        subtractor = _production_like_subtractor(log)
        subtractor.recipeSettings["sky-subtraction"]["bspline_fitting_residual_clipping_sigma"] = 1
        _, contaminated, _ = _clean_and_contaminated_orders(seed=3)
        subtractor.fit_bspline_curve_to_sky(contaminated)
        callsPerLimit[passLimit] = len(calls)

    # ASSERT
    assert callsPerLimit[3] - callsPerLimit[2] == 1
