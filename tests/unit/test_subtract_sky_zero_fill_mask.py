"""Zero-filled sky-subtraction pixels are flagged, and the zero clip of the sky model is pinned (DY-1284)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

FRAME_SHAPE = (3, 3)
# DETECTOR (ROW, COLUMN) POSITIONS THE SKY MODEL REACHES, AND THE ONE IT REACHES WITH A NAN RESULT
MODELLED_PIXELS = [(0, 0), (1, 1)]
NAN_MODEL_PIXEL = (2, 2)


def _subtractor(log: Any, *, uncertainty: np.ndarray) -> subtract_sky:
    """Build a subtractor whose fit stage returns the map pixels with fixed sky values."""
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = False
    subtractor.stopSubtraction = False
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    subtractor.objectFrame = CCDData(
        np.full(FRAME_SHAPE, 10.0),
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(uncertainty),
    )
    # THE LAST MAP PIXEL IS IN AN ORDER BUT THE FIT GAVE IT NO MODEL
    subtractor.mapDF = pd.DataFrame(
        {
            "order": [10, 10, 10],
            "x": [0, 1, 2],
            "y": [0, 1, 2],
            "sky_model": [4.0, 5.0, np.nan],
            "sky_subtracted_flux": [6.0, 5.0, np.nan],
            "error": [2.0, 1.0, 1.0],
        }
    )
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "bspline_order": 3,
            "clip-slit-edge-fraction": 0.1,
            "aggressive_object_masking": False,
            "sky_model_qc_plot": False,
        }
    }
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.copy()
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = lambda frame: (
        frame.copy(),
        (None, None, 3),
        np.array([501.0]),
        np.array([1.0, 2.0]),
        0.5,
    )
    return subtractor


def _modelled_pixel_mask() -> np.ndarray:
    modelled = np.zeros(FRAME_SHAPE, dtype=bool)
    for row, column in MODELLED_PIXELS:
        modelled[row, column] = True
    return modelled


def test_every_zero_filled_pixel_carries_a_qual_flag_in_the_sky_subtracted_frame(log: Any) -> None:
    subtractor = _subtractor(log, uncertainty=np.ones(FRAME_SHAPE))

    _, subtracted, _, _, _ = subtractor.subtract()

    modelled = _modelled_pixel_mask()
    filled = subtracted.data == 0
    assert np.count_nonzero(filled) == 7
    assert subtracted.mask[filled].all()
    assert not subtracted.mask[modelled].any()


def test_every_zero_filled_pixel_carries_a_qual_flag_in_the_sky_model(log: Any) -> None:
    subtractor = _subtractor(log, uncertainty=np.ones(FRAME_SHAPE))

    model, _, _, _, _ = subtractor.subtract()

    modelled = _modelled_pixel_mask()
    assert model.mask[model.data == 0].all()
    assert not model.mask[modelled].any()


def test_a_pixel_inside_the_map_with_no_model_is_flagged_and_zero_filled(log: Any) -> None:
    subtractor = _subtractor(log, uncertainty=np.ones(FRAME_SHAPE))

    model, subtracted, _, _, _ = subtractor.subtract()

    assert model.data[NAN_MODEL_PIXEL] == 0
    assert subtracted.data[NAN_MODEL_PIXEL] == 0
    assert model.mask[NAN_MODEL_PIXEL]
    assert subtracted.mask[NAN_MODEL_PIXEL]


def test_no_output_uncertainty_is_zero_unless_the_pixel_is_masked(log: Any) -> None:
    uncertainty = np.ones(FRAME_SHAPE)
    uncertainty[0, 2] = np.nan
    uncertainty[NAN_MODEL_PIXEL] = np.nan
    subtractor = _subtractor(log, uncertainty=uncertainty)

    model, subtracted, _, _, _ = subtractor.subtract()

    for frame in (model, subtracted):
        assert (frame.uncertainty.array[~frame.mask] != 0).all()


def test_a_nan_uncertainty_is_left_nan_and_masked(log: Any) -> None:
    uncertainty = np.ones(FRAME_SHAPE)
    uncertainty[0, 2] = np.nan
    subtractor = _subtractor(log, uncertainty=uncertainty)

    model, subtracted, _, _, _ = subtractor.subtract()

    for frame in (model, subtracted):
        assert np.isnan(frame.uncertainty.array[0, 2])
        assert frame.mask[0, 2]


def test_an_existing_object_frame_mask_survives_the_flagging(log: Any) -> None:
    subtractor = _subtractor(log, uncertainty=np.ones(FRAME_SHAPE))
    subtractor.objectFrame.mask[1, 1] = True

    model, subtracted, _, _, _ = subtractor.subtract()

    assert subtracted.mask[1, 1]
    assert model.mask[1, 1]
    assert not subtracted.mask[0, 0]


def _flat_sky_pixels(flux: np.ndarray) -> pd.DataFrame:
    wavelength = np.linspace(500.0, 510.0, flux.size)
    return pd.DataFrame(
        {
            "order": 10,
            "wavelength": wavelength,
            "flux": flux,
            "error": np.ones(wavelength.size),
            "residual_windowed_std": np.ones(wavelength.size),
            "flagged_all_clipped": False,
            "flagged_noisy_region": False,
            "residual_windowed_long_median": np.zeros(wavelength.size),
            "flux_windowed_long_median": np.zeros(wavelength.size),
        }
    )


def _fitter(log: Any) -> subtract_sky:
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "NIR"
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


def test_the_zero_clip_of_the_sky_model_reports_how_many_pixels_it_clipped_per_order(log: Any) -> None:
    pixelCount = 64
    subtractor = _fitter(log)

    subtractor.fit_bspline_curve_to_sky(_flat_sky_pixels(np.full(pixelCount, -5.0)))

    clipMessages = [message for _, message in log.messages if "clipped" in message and "zero" in message]
    assert len(clipMessages) == 1
    assert f"{pixelCount} of {pixelCount}" in clipMessages[0]
    assert "order 10" in clipMessages[0]


def test_a_sky_model_that_is_never_negative_reports_no_zero_clip(log: Any) -> None:
    subtractor = _fitter(log)

    subtractor.fit_bspline_curve_to_sky(_flat_sky_pixels(np.full(64, 12.0)))

    assert not [message for _, message in log.messages if "zero" in message]


@pytest.mark.xfail(
    strict=True,
    reason="the sky model is clipped at zero, so a sky that is 0 in truth gets a positive mean; DY-1282 removes it",
)
def test_a_sky_that_is_zero_in_truth_gives_a_model_with_mean_zero(log: Any) -> None:
    # THIS SEED GIVES A FIT THAT CROSSES ZERO (UNCLIPPED MEAN 0.015, CLIPPED MEAN 0.234)
    noise = np.random.default_rng(10).normal(0.0, 1.0, 64)
    noise -= noise.mean()
    subtractor = _fitter(log)

    modelled, _, _, _, _ = subtractor.fit_bspline_curve_to_sky(_flat_sky_pixels(noise))

    assert modelled["sky_model"].mean() == pytest.approx(0.0, abs=0.05)
