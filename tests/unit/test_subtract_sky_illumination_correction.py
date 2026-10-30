"""Stare sky subtraction can multiply the sky model by a skyline-measured slit-illumination factor (DY-2054)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils.slit_illumination_correction as illumination
from soxspipe.commonutils.subtract_sky import ILLUMINATION_OBJECT_MARGIN_ARCSEC, subtract_sky

pytestmark = pytest.mark.unit

SETTING = "skyline_illumination_correction"
COLUMNS = 700
SLIT_ROWS = 41
SLIT_STEP = 0.25
READ_NOISE = 3.0
CONTINUUM = 20.0
LINE_AMPLITUDE = 6000.0
LINE_SPACING = 25
LINE_SIGMA_COLUMNS = 1.5
GRADIENT_PER_ARCSEC = 0.004
# THE SLOPE OF THE FRACTIONAL RESIDUAL OF THE BRIGHT PIXELS AFTER THE CORRECTION. THE NOISE ON A FIT TO 41 BIN MEANS IS
# ABOUT 0.0003 PER ARCSEC
RESIDUAL_SLOPE_LIMIT = 0.0010
# THE CURVE IS THE MEAN OF ABOUT 100 PIXELS PER BIN
BIN_TOLERANCE = 0.004


def _order_frame(isObjectContaminated: bool = False, lineAmplitude: float = LINE_AMPLITUDE) -> pd.DataFrame:
    """One order of a wavelength-only sky of narrow lines whose flux has a known gradient along the slit."""
    rng = np.random.default_rng(4)
    columns = np.arange(COLUMNS, dtype=float)
    rows = (np.arange(SLIT_ROWS) - SLIT_ROWS // 2) * SLIT_STEP
    sky = np.full(COLUMNS, CONTINUUM)
    for centre in range(LINE_SPACING, COLUMNS - LINE_SPACING, LINE_SPACING):
        sky += lineAmplitude * np.exp(-0.5 * ((columns - centre) / LINE_SIGMA_COLUMNS) ** 2)
    slitPosition = np.repeat(rows, COLUMNS)
    skyModel = np.tile(sky, SLIT_ROWS)
    error = np.sqrt(skyModel + READ_NOISE**2)
    flux = skyModel * (1.0 + GRADIENT_PER_ARCSEC * slitPosition) + rng.normal(0.0, 1.0, skyModel.size) * error
    if isObjectContaminated:
        flux = np.where((slitPosition > 0.0) & (slitPosition < 1.0), flux + 0.09 * skyModel, flux)
    return pd.DataFrame(
        {
            "order": 11,
            "x": np.tile(np.arange(COLUMNS), SLIT_ROWS),
            "y": np.repeat(np.arange(SLIT_ROWS), COLUMNS),
            "wavelength": np.tile(1500.0 + 0.05 * columns, SLIT_ROWS),
            "slit_position": slitPosition,
            "flux": flux,
            "error": error,
            "sky_model": skyModel,
            "sky_subtracted_flux": flux - skyModel,
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
            "flagged_edge_clipped": False,
            "flagged_bad_pixel_clipped": False,
        }
    )


def _stub_fit_stages(subtractor: subtract_sky) -> None:
    """Make the clip and fit stages hand the map pixels straight through, with their sky model already fitted."""
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.copy()
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = lambda frame: (
        frame.copy(),
        (None, None, 3),
        np.array([501.0]),
        np.array([1.0, 2.0]),
        0.5,
    )


def _subtractor(log: Any, mapDF: pd.DataFrame, settingValue: Any = None) -> subtract_sky:
    """A subtractor whose fit returns the given order, with the illumination setting set when a value is given."""
    subtractor = object.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "NIR"
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.debug = False
    subtractor.detectorParams = {"dispersion-axis": "x"}
    subtractor.dateObs = "2024-01-02T03:04:05"
    shape = (SLIT_ROWS, COLUMNS)
    subtractor.objectFrame = CCDData(
        np.full(shape, 10.0),
        unit=u.electron,
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(shape, 2.0)),
    )
    subtractor.mapDF = mapDF
    subtractor.qc = pd.DataFrame()
    subtractor.products = pd.DataFrame()
    skySettings = {
        "bspline_order": 3,
        "clip-slit-edge-fraction": 0.1,
        "aggressive_object_masking": False,
        "object_wing_correction": False,
        "sky_error_inflation": False,
        "sky_model_qc_plot": False,
    }
    if settingValue is not None:
        skySettings[SETTING] = settingValue
    subtractor.recipeSettings = {"sky-subtraction": skySettings}
    _stub_fit_stages(subtractor)
    return subtractor


def _reduce(log: Any, mapDF: pd.DataFrame, settingValue: Any = None, objectRanges: Any = None):
    subtractor = _subtractor(log, mapDF, settingValue)
    if objectRanges is not None:
        subtractor.objectSlitRanges = objectRanges
    return subtractor.subtract()


def _residual_slope(mapDF: pd.DataFrame, skySubtracted: np.ndarray) -> float:
    """The slope along the slit of the sky-subtracted flux as a fraction of the sky, over the bright lines."""
    isBright = (mapDF["sky_model"] > 1000.0).to_numpy()
    fractional = skySubtracted[mapDF["y"].to_numpy(), mapDF["x"].to_numpy()] / mapDF["sky_model"].to_numpy()
    return float(np.polyfit(mapDF["slit_position"].to_numpy()[isBright], fractional[isBright], 1)[0])


def test_the_sky_subtracted_frame_is_unchanged_when_the_setting_is_absent(log: Any, monkeypatch) -> None:
    # ARRANGE
    mapDF = _order_frame()
    monkeypatch.setattr(illumination, "measure_slit_illumination", _must_not_be_called)

    # ACT
    skyModel, skySubtracted, _, _, _ = _reduce(log, mapDF)

    # ASSERT
    np.testing.assert_array_equal(skySubtracted.data[mapDF["y"], mapDF["x"]], mapDF["sky_subtracted_flux"])
    np.testing.assert_array_equal(skyModel.data[mapDF["y"], mapDF["x"]], mapDF["sky_model"])


def _must_not_be_called(*args: Any, **kwargs: Any) -> None:
    raise AssertionError("the illumination correction ran with its setting off")


def test_the_sky_subtracted_frame_is_unchanged_when_the_setting_is_false(log: Any, monkeypatch) -> None:
    # ARRANGE
    mapDF = _order_frame()
    monkeypatch.setattr(illumination, "measure_slit_illumination", _must_not_be_called)

    # ACT
    _, skySubtracted, _, _, _ = _reduce(log, mapDF, settingValue=False)

    # ASSERT
    np.testing.assert_array_equal(skySubtracted.data[mapDF["y"], mapDF["x"]], mapDF["sky_subtracted_flux"])


def test_the_gradient_along_the_slit_is_removed_from_the_sky_subtracted_frame_when_on(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame()

    # ACT
    _, withoutCorrection, _, _, _ = _reduce(log, mapDF, settingValue=False)
    _, withCorrection, _, _, _ = _reduce(log, mapDF, settingValue=True)

    # ASSERT
    assert _residual_slope(mapDF, withoutCorrection.data) == pytest.approx(GRADIENT_PER_ARCSEC, rel=0.1)
    assert abs(_residual_slope(mapDF, withCorrection.data)) < RESIDUAL_SLOPE_LIMIT


def test_the_sky_model_image_is_the_fitted_model_times_the_measured_factor(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame()
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()

    # ACT
    skyModel, skySubtracted, _, _, _ = _reduce(log, mapDF, settingValue=True)

    # ASSERT
    factor = skyModel.data[rows, columns] / mapDF["sky_model"].to_numpy()
    expected = 1.0 + GRADIENT_PER_ARCSEC * mapDF["slit_position"].to_numpy()
    inside = np.abs(mapDF["slit_position"].to_numpy()) < 4.5
    np.testing.assert_allclose(factor[inside], expected[inside], atol=BIN_TOLERANCE)
    np.testing.assert_allclose(
        skySubtracted.data[rows, columns], mapDF["flux"].to_numpy() - skyModel.data[rows, columns]
    )


def test_an_order_with_too_faint_a_sky_is_left_as_it_is_when_the_setting_is_true(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame(lineAmplitude=40.0)

    # ACT
    skyModel, skySubtracted, _, _, _ = _reduce(log, mapDF, settingValue=True)

    # ASSERT
    np.testing.assert_array_equal(skyModel.data[mapDF["y"], mapDF["x"]], mapDF["sky_model"])
    np.testing.assert_array_equal(skySubtracted.data[mapDF["y"], mapDF["x"]], mapDF["sky_subtracted_flux"])


def test_pixels_inside_an_object_slit_range_do_not_set_the_curve(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame(isObjectContaminated=True)
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()
    slitPosition = mapDF["slit_position"].to_numpy()
    outside = (slitPosition < -0.5) | (slitPosition > 1.5)
    expected = 1.0 + GRADIENT_PER_ARCSEC * slitPosition

    # ACT
    withRange, _, _, _, _ = _reduce(log, mapDF, settingValue=True, objectRanges=[[-0.05, 1.05]])
    withoutRange, _, _, _, _ = _reduce(log, mapDF, settingValue=True, objectRanges=[])

    # ASSERT
    factorWithRange = withRange.data[rows, columns] / mapDF["sky_model"].to_numpy()
    factorWithoutRange = withoutRange.data[rows, columns] / mapDF["sky_model"].to_numpy()
    np.testing.assert_allclose(factorWithRange[outside], expected[outside], atol=BIN_TOLERANCE)
    assert factorWithoutRange[(slitPosition > 0.2) & (slitPosition < 0.8)].mean() > 1.03


def test_the_sky_model_across_an_object_slit_range_follows_the_curve_fitted_to_the_bins(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame()
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()
    slitPosition = mapDF["slit_position"].to_numpy()
    # A RANGE AWAY FROM THE SLIT CENTRE, WHERE THE CURVE IS ABOUT +1%. ITS ENDS ARE ON SAMPLED SLIT POSITIONS
    isInRange = (slitPosition >= 2.0) & (slitPosition <= 3.0)
    assert np.isclose(slitPosition, 2.0).any() and np.isclose(slitPosition, 3.0).any()
    expected = 1.0 + GRADIENT_PER_ARCSEC * slitPosition

    # ACT
    skyModel, skySubtracted, _, _, _ = _reduce(log, mapDF, settingValue=True, objectRanges=[[2.0, 3.0]])

    # ASSERT
    factor = skyModel.data[rows, columns] / mapDF["sky_model"].to_numpy()
    isInside = np.abs(slitPosition) < 4.5
    # THE CURVE IS ZERO ON AVERAGE OVER THE BINS USED, WHICH LEAVE OUT THE RANGE AND ITS MARGINS, SO ITS ZERO POINT IS
    # SHIFTED FROM THE SLIT CENTRE; ONLY ITS SHAPE IS COMPARED
    offset = np.mean(factor[isInside] - expected[isInside])
    np.testing.assert_allclose(factor[isInside], expected[isInside] + offset, atol=BIN_TOLERANCE)
    assert factor[isInRange].mean() - offset > 1.006
    np.testing.assert_allclose(
        skySubtracted.data[rows, columns], mapDF["flux"].to_numpy() - skyModel.data[rows, columns], rtol=1e-12
    )


def test_pixels_within_the_margin_around_an_object_slit_range_do_not_set_the_curve(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame()
    slitPosition = mapDF["slit_position"].to_numpy()
    # OBJECT LIGHT THAT REACHES 0.9 ARCSEC PAST THE RANGE END, OR A WING MODEL THAT TAKES TOO MUCH THERE, LIFTS THE BINS
    # BY 9%
    reach = 0.9
    assert reach < ILLUMINATION_OBJECT_MARGIN_ARCSEC
    isNearRange = (slitPosition > 1.05) & (slitPosition < 1.05 + reach)
    lifted = mapDF.assign(flux=np.where(isNearRange, mapDF["flux"] + 0.09 * mapDF["sky_model"], mapDF["flux"]))
    lifted["sky_subtracted_flux"] = lifted["flux"] - lifted["sky_model"]
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()
    expected = 1.0 + GRADIENT_PER_ARCSEC * slitPosition

    # ACT
    skyModel, _, _, _, _ = _reduce(log, lifted, settingValue=True, objectRanges=[[-0.05, 1.05]])

    # ASSERT
    factor = skyModel.data[rows, columns] / mapDF["sky_model"].to_numpy()
    isInside = np.abs(slitPosition) < 4.5
    offset = np.mean(factor[isInside] - expected[isInside])
    np.testing.assert_allclose(factor[isInside], expected[isInside] + offset, atol=BIN_TOLERANCE)


def test_the_object_wing_model_is_taken_from_the_flux_before_the_curve_is_measured(log: Any) -> None:
    # ARRANGE
    mapDF = _order_frame()
    wings = np.where(mapDF["slit_position"] > 1.5, 0.05 * mapDF["sky_model"], 0.0)
    withWings = mapDF.assign(flux=mapDF["flux"] + wings, object_wing_model=wings)
    withWings["sky_subtracted_flux"] = withWings["flux"] - withWings["sky_model"]
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()

    # ACT
    skyModel, _, _, _, _ = _reduce(log, withWings, settingValue=True)

    # ASSERT
    factor = skyModel.data[rows, columns] / mapDF["sky_model"].to_numpy()
    expected = 1.0 + GRADIENT_PER_ARCSEC * mapDF["slit_position"].to_numpy()
    inside = np.abs(mapDF["slit_position"].to_numpy()) < 4.5
    np.testing.assert_allclose(factor[inside], expected[inside], atol=BIN_TOLERANCE)


@pytest.mark.parametrize(
    "flagColumn",
    ["flagged_edge_clipped", "flagged_bad_pixel_clipped", "flagged_object_clipped", "flagged_all_clipped"],
)
def test_flagged_pixels_do_not_set_the_curve(log: Any, flagColumn: str) -> None:
    # ARRANGE
    mapDF = _order_frame()
    isFlagged = (mapDF["slit_position"] > 0.0) & (mapDF["slit_position"] < 1.0)
    flagged = mapDF.assign(flux=np.where(isFlagged, mapDF["flux"] + 0.09 * mapDF["sky_model"], mapDF["flux"]))
    flagged[flagColumn] = isFlagged
    clean = mapDF.assign(flux=flagged["flux"])
    clean[flagColumn] = False

    # ACT
    skyModelFlagged, _, _, _, _ = _reduce(log, flagged, settingValue=True)
    skyModelClean, _, _, _, _ = _reduce(log, clean, settingValue=True)

    # ASSERT
    rows, columns = mapDF["y"].to_numpy(), mapDF["x"].to_numpy()
    factorFlagged = skyModelFlagged.data[rows, columns] / mapDF["sky_model"].to_numpy()
    factorClean = skyModelClean.data[rows, columns] / mapDF["sky_model"].to_numpy()
    expected = 1.0 + GRADIENT_PER_ARCSEC * mapDF["slit_position"].to_numpy()
    isFlaggedCore = ((mapDF["slit_position"] > 0.2) & (mapDF["slit_position"] < 0.8)).to_numpy()
    np.testing.assert_allclose(factorFlagged[isFlaggedCore], expected[isFlaggedCore], atol=2 * BIN_TOLERANCE)
    assert factorClean[isFlaggedCore].mean() > 1.03


@pytest.mark.parametrize("badValue", ["True", 1, "false"])
def test_a_setting_that_is_not_a_bool_raises_a_value_error_naming_the_key(log: Any, badValue: Any) -> None:
    # ARRANGE
    subtractor = _subtractor(log, _order_frame(), settingValue=badValue)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=SETTING):
        subtractor.subtract()


def _flagged_orders() -> list[pd.DataFrame]:
    """Two orders with 70% of the pixels inside a 1.0 to 1.6 arcsec object flagged, plus 3% scattered flags."""
    count = 4000
    orders = []
    for offset, seed in enumerate((11, 12)):
        rng = np.random.default_rng(seed)
        slitPosition = rng.uniform(-5.5, 5.5, count)
        isObject = (slitPosition > 1.0) & (slitPosition < 1.6) & (rng.uniform(0, 1, count) < 0.7)
        orders.append(
            pd.DataFrame(
                {
                    "order": 10 + offset,
                    "slit_position": slitPosition,
                    "flagged_object_clipped": isObject | (rng.uniform(0, 1, count) < 0.03),
                    "flagged_all_clipped": False,
                    "flux_minus_smoothed_residual": 3.0 * np.sin(np.arange(count) / 5.0),
                    "flux_percentile_smoothed": 100.0 + np.cos(np.arange(count) / 50.0),
                }
            )
        )
    return orders


def test_the_clip_records_the_object_slit_ranges_for_the_correction(log: Any) -> None:
    # ARRANGE
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.ron = 3.3
    subtractor.recipeSettings = {
        "sky-subtraction": {"percentile_rolling_window_size": 15, "noise_rolling_window_size": 45}
    }

    # ACT
    subtractor.clip_object_slit_positions(_flagged_orders(), aggressive=True)
    aggressiveRanges = subtractor.objectSlitRanges
    subtractor.clip_object_slit_positions(_flagged_orders(), aggressive=False)

    # ASSERT
    assert len(aggressiveRanges) == 1
    lower, upper = aggressiveRanges[0]
    assert lower <= 1.0 and upper >= 1.6
    assert subtractor.objectSlitRanges == []


def test_the_object_slit_ranges_default_to_an_empty_list_of_one_type(log: Any) -> None:
    # ASSERT
    assert subtract_sky.objectSlitRanges == []
    assert isinstance(subtract_sky.objectSlitRanges, list)
