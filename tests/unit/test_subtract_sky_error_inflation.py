"""The sky-subtracted ERRS is inflated by the local scatter of the sky-fit residuals (DY-1973)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.subtract_sky import (
    ERROR_INFLATION_MIN_SAMPLES,
    ERROR_INFLATION_WINDOW,
    MAD_TO_SIGMA,
    _sky_residual_scatter_factor,
    subtract_sky,
)

pytestmark = pytest.mark.unit

INFLATION_KEY = "sky_error_inflation"
# MARKS A SETTINGS FILE THAT DOES NOT CARRY THE KEY AT ALL (NONE AND FALSE ARE VALUES A FILE CAN HOLD)
_MISSING = object()
# THE SAMPLES OF THE SYNTHETIC ORDERS: SIX WINDOWS, SO A REGION OF ONE SCATTER FILLS WHOLE WINDOWS
ORDER_SAMPLES = 6 * ERROR_INFLATION_WINDOW
# THE SCATTER (IN UNITS OF THE STATED ERROR) OF THE SYNTHETIC SKYLINE RESIDUAL
LINE_SCATTER = 6.0
# THE TOLERANCE ON THE MEAN OF A SCATTER MEASURED FROM MANY WINDOWS OF NOISE
SCATTER_TOLERANCE = 0.2
# THE LARGEST FACTOR THE NOISE OF THE MAD OF ONE WINDOW (ABOUT 14% FOR 100 SAMPLES) LIFTS A QUIET ORDER TO
NOISE_FACTOR_LIMIT = 1.5
STATED_ERROR = 2.0
FRAME_COLUMNS = 4000
# THE AMPLITUDE (IN E-) OF THE SMOOTH OBJECT-WING LIGHT THE SKY FIT LEFT IN THE FLUX: 6 STATED ERRORS
WING_AMPLITUDE = 6.0 * STATED_ERROR


def _wavelength(samples: int) -> np.ndarray:
    return np.linspace(1500.0, 1600.0, samples)


def _noise(scatter: np.ndarray | float, samples: int, seed: int = 7) -> np.ndarray:
    return np.random.default_rng(seed).normal(0.0, 1.0, samples) * scatter


def test_scale_factor_is_one_when_the_residuals_are_smaller_than_their_errors() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    residual = _noise(0.5, ORDER_SAMPLES)

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    assert (factor == 1.0).all()


def test_scale_factor_stays_close_to_one_for_residuals_that_are_pure_noise() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    residual = _noise(1.0, ORDER_SAMPLES)

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    assert factor.min() >= 1.0
    assert factor.max() < NOISE_FACTOR_LIMIT
    assert factor.mean() < 1.0 + SCATTER_TOLERANCE / 2


def test_scale_factor_follows_the_local_scatter_of_a_skyline_residual() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    isLine = np.zeros(ORDER_SAMPLES, dtype=bool)
    isLine[2 * ERROR_INFLATION_WINDOW : 4 * ERROR_INFLATION_WINDOW] = True
    residual = _noise(np.where(isLine, LINE_SCATTER, 1.0), ORDER_SAMPLES)

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    lineCentre = 3 * ERROR_INFLATION_WINDOW
    lineCore = factor[lineCentre - ERROR_INFLATION_WINDOW // 2 : lineCentre + ERROR_INFLATION_WINDOW // 2]
    assert lineCore.mean() == pytest.approx(LINE_SCATTER, rel=SCATTER_TOLERANCE)
    assert factor[:ERROR_INFLATION_WINDOW].max() < NOISE_FACTOR_LIMIT
    assert factor[-ERROR_INFLATION_WINDOW:].max() < NOISE_FACTOR_LIMIT


def test_scale_factor_is_not_raised_by_a_few_outliers() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    residual = _noise(1.0, ORDER_SAMPLES)
    residual[:: ERROR_INFLATION_WINDOW // 5] = 100.0

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    assert factor.max() < NOISE_FACTOR_LIMIT


def test_scale_factor_follows_the_row_order_of_the_input_not_the_wavelength_order() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    isLine = (wavelength > 1500.0 + 100.0 * 2 / 6) & (wavelength < 1500.0 + 100.0 * 4 / 6)
    residual = _noise(np.where(isLine, LINE_SCATTER, 1.0), ORDER_SAMPLES)
    shuffle = np.random.default_rng(3).permutation(ORDER_SAMPLES)

    # ACT
    inOrder = _sky_residual_scatter_factor(wavelength, residual)
    shuffled = _sky_residual_scatter_factor(wavelength[shuffle], residual[shuffle])

    # ASSERT
    np.testing.assert_allclose(shuffled, inOrder[shuffle])


def test_non_finite_residuals_are_left_out_and_every_pixel_still_gets_a_finite_factor() -> None:
    # ARRANGE
    wavelength = _wavelength(ORDER_SAMPLES)
    residual = _noise(LINE_SCATTER, ORDER_SAMPLES)
    residual[::7] = np.nan
    residual[3::11] = np.inf

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    assert factor.shape == wavelength.shape
    assert np.isfinite(factor).all()
    assert factor.mean() == pytest.approx(LINE_SCATTER, rel=SCATTER_TOLERANCE)


@pytest.mark.parametrize("samples", [0, 1, ERROR_INFLATION_MIN_SAMPLES - 1])
def test_scale_factor_is_one_when_there_are_too_few_residuals_to_measure_scatter(samples: int) -> None:
    # ARRANGE
    wavelength = _wavelength(50)
    residual = np.full(50, np.nan)
    residual[:samples] = LINE_SCATTER * np.array([1.0, -1.0] * 25)[:samples]

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    np.testing.assert_array_equal(factor, np.ones(50))


def test_an_order_shorter_than_the_window_is_measured_with_one_window_of_its_own_length() -> None:
    # ARRANGE
    samples = ERROR_INFLATION_WINDOW // 2
    wavelength = _wavelength(samples)
    residual = _noise(LINE_SCATTER, samples)

    # ACT
    factor = _sky_residual_scatter_factor(wavelength, residual)

    # ASSERT
    expected = MAD_TO_SIGMA * np.median(np.abs(residual - np.median(residual)))
    np.testing.assert_allclose(factor, expected)


def _order_frame(scatter: np.ndarray) -> pd.DataFrame:
    """One fitted order whose residuals (in stated errors) have the given scatter."""
    skyModel = np.full(FRAME_COLUMNS, 5.0)
    flux = skyModel + STATED_ERROR * _noise(scatter, FRAME_COLUMNS)
    return pd.DataFrame(
        {
            "order": 10,
            "x": np.arange(FRAME_COLUMNS),
            "y": 0,
            "wavelength": _wavelength(FRAME_COLUMNS),
            "flux": flux,
            "error": STATED_ERROR,
            "sky_model": skyModel,
            "sky_subtracted_flux": flux - skyModel,
            "flagged_all_clipped": False,
        }
    )


def _stub_fit_stages(subtractor: subtract_sky) -> None:
    """Make the clip and fit stages hand the map pixels straight through."""
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.copy()
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = lambda frame: (
        frame.copy(),
        (None, None, 3),
        np.array([501.0]),
        np.array([1.0, 2.0]),
        0.5,
    )


def _order_subtractor(log: Any, *, arm: str, inflationSetting: Any, scatter: np.ndarray) -> subtract_sky:
    """Build a subtractor whose fit returns one order whose residuals (in stated errors) have the given scatter."""
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = False
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    subtractor.objectFrame = CCDData(
        np.full((1, FRAME_COLUMNS), 10.0),
        unit=u.electron,
        mask=np.zeros((1, FRAME_COLUMNS), dtype=bool),
        uncertainty=StdDevUncertainty(np.full((1, FRAME_COLUMNS), STATED_ERROR)),
    )
    subtractor.mapDF = _order_frame(scatter)
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()
    skySettings = {
        "bspline_order": 3,
        "clip-slit-edge-fraction": 0.1,
        "aggressive_object_masking": False,
        "object_wing_correction": False,
        "sky_model_qc_plot": False,
    }
    if inflationSetting is not _MISSING:
        skySettings[INFLATION_KEY] = inflationSetting
    subtractor.recipeSettings = {"sky-subtraction": skySettings}
    _stub_fit_stages(subtractor)
    return subtractor


def _line_in_second_half() -> np.ndarray:
    scatter = np.ones(FRAME_COLUMNS)
    scatter[FRAME_COLUMNS // 2 :] = LINE_SCATTER
    return scatter


def _subtracted_errors(log: Any, arm: str, inflationSetting: Any) -> np.ndarray:
    subtractor = _order_subtractor(log, arm=arm, inflationSetting=inflationSetting, scatter=_line_in_second_half())
    _, subtracted, _, _, _ = subtractor.subtract()
    return subtracted.uncertainty.array[0]


def test_nir_inflates_the_errors_where_the_sky_residuals_scatter_by_default(log: Any) -> None:
    # ACT
    errors = _subtracted_errors(log, "NIR", _MISSING)

    # ASSERT
    quietEnd = errors[: FRAME_COLUMNS // 4]
    lineEnd = errors[-FRAME_COLUMNS // 4 :]
    assert (quietEnd >= STATED_ERROR).all()
    assert quietEnd.max() < STATED_ERROR * NOISE_FACTOR_LIMIT
    assert lineEnd.mean() == pytest.approx(STATED_ERROR * LINE_SCATTER, rel=SCATTER_TOLERANCE)


@pytest.mark.parametrize(("arm", "inflationSetting"), [("NIR", False), ("VIS", _MISSING), ("UVB", _MISSING)])
def test_errors_are_untouched_when_the_inflation_is_off(log: Any, arm: str, inflationSetting: Any) -> None:
    # ACT
    errors = _subtracted_errors(log, arm, inflationSetting)

    # ASSERT
    np.testing.assert_array_equal(errors, np.full(FRAME_COLUMNS, STATED_ERROR))


def test_the_setting_turns_the_inflation_on_for_an_arm_that_defaults_to_off(log: Any) -> None:
    # ACT
    errors = _subtracted_errors(log, "VIS", True)

    # ASSERT
    assert errors[-FRAME_COLUMNS // 4 :].mean() == pytest.approx(STATED_ERROR * LINE_SCATTER, rel=SCATTER_TOLERANCE)


def _clipped_order_errors(log: Any, outlier: float) -> np.ndarray:
    """The errors of an order whose every tenth pixel is clipped, and whose clipped pixels are `outlier` high."""
    subtractor = _order_subtractor(log, arm="NIR", inflationSetting=True, scatter=_line_in_second_half())
    clipped = np.zeros(FRAME_COLUMNS, dtype=bool)
    clipped[::10] = True
    mapDF = subtractor.mapDF
    mapDF.loc[clipped, "flux"] += outlier
    mapDF.loc[clipped, "sky_subtracted_flux"] += outlier
    mapDF["flagged_all_clipped"] = clipped
    _, subtracted, _, _, _ = subtractor.subtract()
    return subtracted.uncertainty.array[0]


def test_clipped_pixels_do_not_set_the_scatter_but_still_get_the_factor(log: Any) -> None:
    # ACT
    withOutliers = _clipped_order_errors(log, outlier=1000.0)
    withoutOutliers = _clipped_order_errors(log, outlier=0.0)

    # ASSERT
    np.testing.assert_array_equal(withOutliers, withoutOutliers)
    assert withOutliers[-FRAME_COLUMNS // 4 :].mean() > STATED_ERROR * LINE_SCATTER / 2


def test_the_residuals_image_stays_in_units_of_the_stated_error(log: Any) -> None:
    # ARRANGE
    subtractor = _order_subtractor(log, arm="NIR", inflationSetting=True, scatter=_line_in_second_half())
    expected = (subtractor.mapDF["sky_subtracted_flux"] / subtractor.mapDF["error"]).to_numpy()

    # ACT
    _, _, residuals, _, _ = subtractor.subtract()

    # ASSERT
    np.testing.assert_allclose(residuals.data[0], expected)


@pytest.mark.parametrize("badValue", ["False", 1, None])
def test_a_setting_that_is_not_a_bool_raises_a_value_error_naming_the_key(log: Any, badValue: Any) -> None:
    # ARRANGE
    subtractor = _order_subtractor(log, arm="NIR", inflationSetting=badValue, scatter=np.ones(FRAME_COLUMNS))

    # ACT / ASSERT
    with pytest.raises(ValueError, match=INFLATION_KEY):
        subtractor.subtract()


def _errors_with_wing_light(log: Any, hasWingModel: bool) -> np.ndarray:
    """The errors of a skyline-free order carrying smooth object-wing light that the sky fit did not remove."""
    subtractor = _order_subtractor(log, arm="NIR", inflationSetting=True, scatter=np.ones(FRAME_COLUMNS))
    wingLight = WING_AMPLITUDE * np.sin(np.linspace(0.0, 40.0, FRAME_COLUMNS))
    mapDF = subtractor.mapDF
    mapDF["flux"] += wingLight
    mapDF["sky_subtracted_flux"] += wingLight
    if hasWingModel:
        mapDF["object_wing_model"] = wingLight
    _, subtracted, _, _, _ = subtractor.subtract()
    return subtracted.uncertainty.array[0]


def test_object_wing_light_removed_by_the_refit_does_not_inflate_the_errors(log: Any) -> None:
    # ACT
    errors = _errors_with_wing_light(log, hasWingModel=True)

    # ASSERT
    assert errors.max() < STATED_ERROR * NOISE_FACTOR_LIMIT


def test_object_wing_light_with_no_wing_model_inflates_the_errors(log: Any) -> None:
    # ACT
    errors = _errors_with_wing_light(log, hasWingModel=False)

    # ASSERT
    assert errors.mean() > 1.4 * STATED_ERROR
