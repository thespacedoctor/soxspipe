"""The Kelson-style skyline slit-illumination correction measures a multiplicative factor along the slit (DY-2054)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from soxspipe.commonutils.slit_illumination_correction import (
    CONTINUUM_SKY_SNR,
    MAX_BIN_ERROR,
    MIN_SKY_SNR,
    measure_slit_illumination,
)

pytestmark = pytest.mark.unit

COLUMNS = 700
SLIT_ROWS = 41
SLIT_STEP = 0.25
READ_NOISE = 3.0
CONTINUUM = 20.0
LINE_AMPLITUDE = 6000.0
LINE_SIGMA_COLUMNS = 1.5
LINE_SPACING = 25
# THE FRACTIONAL ILLUMINATION GRADIENT INJECTED INTO THE DATA: 0.4% PER ARCSEC, +-2% ACROSS THE SLIT
GRADIENT_PER_ARCSEC = 0.004
# THE BIN MEAN OF ABOUT 100 PIXELS AT A SIGNAL-TO-NOISE OF 30 TO 80 IS GOOD TO ABOUT 0.5%
BIN_TOLERANCE = 0.004
SLOPE_TOLERANCE = 0.15


@dataclass
class SyntheticOrder:
    slitPosition: np.ndarray
    wavelength: np.ndarray
    skyModel: np.ndarray
    data: np.ndarray
    error: np.ndarray
    isUsable: np.ndarray


def _synthetic_order(
    gradientPerArcsec: float = GRADIENT_PER_ARCSEC,
    lineAmplitude: float = LINE_AMPLITUDE,
    seed: int = 4,
    lineSpacing: int = LINE_SPACING,
):
    """One order: a wavelength-only sky of narrow lines on a continuum, scaled along the slit by a known gradient."""
    rng = np.random.default_rng(seed)
    columns = np.arange(COLUMNS, dtype=float)
    rows = (np.arange(SLIT_ROWS) - SLIT_ROWS // 2) * SLIT_STEP
    sky = np.full(COLUMNS, CONTINUUM)
    for centre in range(lineSpacing, COLUMNS - lineSpacing, lineSpacing):
        sky += lineAmplitude * np.exp(-0.5 * ((columns - centre) / LINE_SIGMA_COLUMNS) ** 2)
    slitPosition = np.repeat(rows, COLUMNS)
    wavelength = np.tile(1500.0 + 0.05 * columns, SLIT_ROWS)
    skyModel = np.tile(sky, SLIT_ROWS)
    error = np.sqrt(skyModel + READ_NOISE**2)
    data = skyModel * (1.0 + gradientPerArcsec * slitPosition) + rng.normal(0.0, 1.0, skyModel.size) * error
    return SyntheticOrder(slitPosition, wavelength, skyModel, data, error, np.ones(skyModel.size, dtype=bool))


def _measure(order: SyntheticOrder, **overrides):
    return measure_slit_illumination(
        order.slitPosition,
        order.wavelength,
        order.skyModel,
        overrides.pop("data", order.data),
        order.error,
        overrides.pop("isUsable", order.isUsable),
        **overrides,
    )


def test_recovers_the_gradient_injected_into_the_skyline_pixels() -> None:
    # ARRANGE
    order = _synthetic_order()

    # ACT
    correction = _measure(order)

    # ASSERT
    assert correction.isApplied
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    interior = np.abs(order.slitPosition) < 4.5
    np.testing.assert_allclose(correction.factor[interior], expected[interior], atol=BIN_TOLERANCE)
    slope = np.polyfit(correction.binCentres, correction.binValues, 1)[0]
    assert slope == pytest.approx(GRADIENT_PER_ARCSEC, rel=SLOPE_TOLERANCE)


def test_factor_is_one_when_too_few_pixels_are_bright_enough() -> None:
    # ARRANGE
    order = _synthetic_order(lineAmplitude=40.0)

    # ACT
    correction = _measure(order)

    # ASSERT
    assert not correction.isApplied
    np.testing.assert_array_equal(correction.factor, np.ones(order.slitPosition.size))


def test_factor_is_one_when_every_slit_bin_mean_is_too_noisy() -> None:
    # ARRANGE
    order = _synthetic_order(lineSpacing=100, lineAmplitude=300.0)

    # ACT
    correction = _measure(order)

    # ASSERT
    assert correction.pixelCount > 0
    assert not correction.isApplied
    np.testing.assert_array_equal(correction.factor, np.ones(order.slitPosition.size))


def test_factor_is_one_when_no_pixel_is_usable() -> None:
    # ARRANGE
    order = _synthetic_order()

    # ACT
    correction = _measure(order, isUsable=np.zeros(order.slitPosition.size, dtype=bool))

    # ASSERT
    assert not correction.isApplied
    np.testing.assert_array_equal(correction.factor, np.ones(order.slitPosition.size))


def test_pixels_that_are_not_usable_do_not_set_the_curve() -> None:
    # ARRANGE
    order = _synthetic_order()
    isObject = (order.slitPosition > 0.0) & (order.slitPosition < 1.0)
    contaminated = np.where(isObject, order.data + 0.05 * order.skyModel, order.data)

    # ACT
    unmasked = _measure(order, data=contaminated)
    masked = _measure(order, data=contaminated, isUsable=~isObject)

    # ASSERT
    inObject = isObject & (order.wavelength > 0)
    assert unmasked.factor[inObject].mean() > 1.03
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    np.testing.assert_allclose(masked.factor[~isObject], expected[~isObject], atol=BIN_TOLERANCE)


def test_the_object_model_is_taken_from_the_data_before_the_curve_is_measured() -> None:
    # ARRANGE
    order = _synthetic_order()
    objectModel = np.where(order.slitPosition > 1.0, 0.03 * order.skyModel, 0.0)
    contaminated = order.data + objectModel

    # ACT
    withModel = _measure(order, data=contaminated - objectModel)
    withoutModel = _measure(order, data=contaminated)

    # ASSERT
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    inside = np.abs(order.slitPosition) < 4.5
    np.testing.assert_allclose(withModel.factor[inside], expected[inside], atol=BIN_TOLERANCE)
    isFarSide = order.slitPosition > 2.0
    assert withoutModel.factor[isFarSide].mean() > withModel.factor[isFarSide].mean() + 0.01


def test_pixels_with_a_large_fractional_residual_are_left_out() -> None:
    # ARRANGE
    order = _synthetic_order()
    isHot = np.random.default_rng(9).random(order.data.size) < 0.05
    hot = np.where(isHot, order.data + 0.5 * order.skyModel, order.data)

    # ACT
    correction = _measure(order, data=hot)

    # ASSERT
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    inside = np.abs(order.slitPosition) < 4.5
    np.testing.assert_allclose(correction.factor[inside], expected[inside], atol=BIN_TOLERANCE)


def test_the_ends_of_the_wavelength_range_are_left_out() -> None:
    # ARRANGE
    order = _synthetic_order()
    span = np.ptp(order.wavelength)
    isEnd = (order.wavelength < order.wavelength.min() + 0.05 * span) | (
        order.wavelength > order.wavelength.max() - 0.05 * span
    )
    biased = np.where(isEnd, order.data + 0.09 * order.skyModel, order.data)

    # ACT
    clean = _measure(order)
    withBiasedEnds = _measure(order, data=biased)

    # ASSERT
    np.testing.assert_array_equal(withBiasedEnds.factor, clean.factor)


def test_the_factor_follows_the_input_pixel_order_and_the_inputs_are_not_modified() -> None:
    # ARRANGE
    order = _synthetic_order()
    shuffle = np.random.default_rng(2).permutation(order.slitPosition.size)
    dataBefore = order.data.copy()

    # ACT
    inOrder = _measure(order)
    shuffled = measure_slit_illumination(
        order.slitPosition[shuffle],
        order.wavelength[shuffle],
        order.skyModel[shuffle],
        order.data[shuffle],
        order.error[shuffle],
        order.isUsable[shuffle],
    )

    # ASSERT
    np.testing.assert_allclose(shuffled.factor, inOrder.factor[shuffle])
    np.testing.assert_array_equal(order.data, dataBefore)


def test_non_finite_pixels_are_ignored_and_keep_a_factor_of_one() -> None:
    # ARRANGE
    order = _synthetic_order()
    noisy = order.data.copy()
    noisy[::13] = np.nan
    slitPosition = order.slitPosition.copy()
    slitPosition[5::17] = np.nan

    # ACT
    correction = measure_slit_illumination(
        slitPosition, order.wavelength, order.skyModel, noisy, order.error, order.isUsable
    )

    # ASSERT
    assert correction.isApplied
    assert np.isfinite(correction.factor).all()
    assert (correction.factor[np.isnan(slitPosition)] == 1.0).all()


def test_a_slit_bin_with_no_skyline_pixels_is_interpolated_across() -> None:
    # ARRANGE
    order = _synthetic_order()
    isGap = (order.slitPosition > 1.0) & (order.slitPosition < 2.0)

    # ACT
    correction = _measure(order, isUsable=~isGap)

    # ASSERT
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    np.testing.assert_allclose(correction.factor[isGap], expected[isGap], atol=2 * BIN_TOLERANCE)


def test_a_uniform_offset_of_the_skyline_pixels_is_left_to_the_wavelength_fit() -> None:
    # ARRANGE
    order = _synthetic_order(gradientPerArcsec=0.0)
    offset = order.data + 0.03 * order.skyModel

    # ACT
    correction = _measure(order, data=offset)

    # ASSERT
    np.testing.assert_allclose(correction.factor, 1.0, atol=BIN_TOLERANCE)


# THE SIGNAL-TO-NOISE OF THE STACKED SKYLINE PIXELS OF _dense_bins, AND THE PIXEL COUNT WHOSE BIN MEAN HAS AN ERROR OF
# EXACTLY MAX_BIN_ERROR: 1 / (DENSE_SNR * SQRT(400)) = 0.005
DENSE_SNR = 10.0
PIXELS_AT_MAX_BIN_ERROR = 400


def _dense_bins(positions: list[float], pixelsPerPosition: int):
    """Skyline pixels of signal-to-noise DENSE_SNR stacked at the given slit positions, the same number of faint usable
    pixels (so the order is not continuum-dominated), and two usable pixels that mark the wavelength range."""
    bright = np.concatenate([np.full(pixelsPerPosition, position) for position in positions])
    slit = np.concatenate([bright, bright, [positions[0]] * 2])
    wavelength = np.full(slit.size, 1550.0)
    wavelength[-2:] = [1500.0, 1600.0]
    skyModel = np.full(slit.size, 1000.0)
    skyModel[bright.size : 2 * bright.size] = 1.0
    skyModel[-2:] = 0.0
    error = np.full(slit.size, 1000.0 / DENSE_SNR)
    data = skyModel * (1.0 + 0.01 * slit)
    return slit, wavelength, skyModel, data, error, np.ones(slit.size, dtype=bool)


def test_factor_is_one_when_no_wavelength_is_finite() -> None:
    # ARRANGE
    order = _synthetic_order()
    wavelength = np.full(order.wavelength.size, np.nan)

    # ACT
    correction = measure_slit_illumination(
        order.slitPosition, wavelength, order.skyModel, order.data, order.error, order.isUsable
    )

    # ASSERT
    assert not correction.isApplied
    np.testing.assert_array_equal(correction.factor, np.ones(order.slitPosition.size))


def test_the_wavelength_range_is_that_of_the_usable_pixels() -> None:
    # ARRANGE
    order = _synthetic_order()
    isOutlier = order.slitPosition > 3.0
    wavelength = np.where(isOutlier, 9999.0, order.wavelength)
    isUsable = order.isUsable & ~isOutlier

    # ACT
    withOutliers = measure_slit_illumination(
        order.slitPosition, wavelength, order.skyModel, order.data, order.error, isUsable
    )
    without = measure_slit_illumination(
        order.slitPosition, order.wavelength, order.skyModel, order.data, order.error, isUsable
    )

    # ASSERT
    np.testing.assert_array_equal(withOutliers.factor, without.factor)


def test_a_wide_masked_gap_leaves_the_dense_bins_on_either_side_in_use() -> None:
    # ARRANGE
    order = _synthetic_order()
    isGap = np.abs(order.slitPosition) < 3.5

    # ACT
    correction = _measure(order, isUsable=~isGap)

    # ASSERT
    assert correction.isApplied
    expected = 1.0 + GRADIENT_PER_ARCSEC * order.slitPosition
    np.testing.assert_allclose(correction.factor[~isGap], expected[~isGap], atol=2 * BIN_TOLERANCE)


def test_a_bin_whose_mean_has_exactly_the_maximum_error_is_used() -> None:
    # ARRANGE
    positions = [0.0, 0.3, 0.55, 0.8]
    assert 1.0 / (DENSE_SNR * np.sqrt(PIXELS_AT_MAX_BIN_ERROR)) == MAX_BIN_ERROR

    # ACT
    atMaximum = measure_slit_illumination(*_dense_bins(positions, PIXELS_AT_MAX_BIN_ERROR))
    aboveMaximum = measure_slit_illumination(*_dense_bins(positions, PIXELS_AT_MAX_BIN_ERROR - 1))

    # ASSERT
    assert atMaximum.isApplied
    np.testing.assert_allclose(atMaximum.binErrors, MAX_BIN_ERROR)
    assert not aboveMaximum.isApplied


def test_a_pixel_at_the_far_edge_of_the_slit_falls_in_the_last_bin() -> None:
    # ARRANGE
    slit, *rest = _dense_bins([0.0, 0.3, 0.55, 0.8, 1.0], PIXELS_AT_MAX_BIN_ERROR)

    # ACT
    correction = measure_slit_illumination(slit, *rest)

    # ASSERT
    assert correction.binCentres.size == 4
    assert correction.binCounts[-1] == 2 * PIXELS_AT_MAX_BIN_ERROR
    assert correction.binCentres[-1] == pytest.approx(0.875)


def test_bin_means_are_weighted_by_the_inverse_variance_of_each_pixel() -> None:
    # ARRANGE
    slit, wavelength, skyModel, data, error, isUsable = _dense_bins([0.0, 1.0], PIXELS_AT_MAX_BIN_ERROR)
    # EVERY PIXEL AT TWICE DENSE_SNR, BUT IN THE FIRST BIN HALF THE SKYLINE PIXELS ARE TWICE AS NOISY (STILL AT THE
    # FLOOR) AND HIGH BY 4%: THEIR WEIGHT IS A QUARTER
    isNoisy = np.zeros(slit.size, dtype=bool)
    isNoisy[: PIXELS_AT_MAX_BIN_ERROR // 2] = True
    error = np.where(isNoisy, error, error / 2.0)
    data = np.where(isNoisy, data + 0.04 * skyModel, data)

    # ACT
    correction = measure_slit_illumination(slit, wavelength, skyModel, data, error, isUsable)

    # ASSERT
    # WEIGHTED: (0.25 * 0.04) / 1.25 = 0.008 ABOVE THE SECOND BIN'S 0.01 STEP. UNWEIGHTED IT WOULD BE 0.02
    step = correction.binValues[1] - correction.binValues[0]
    assert step == pytest.approx(0.01 - 0.008, abs=1e-9)


def test_a_continuum_dominated_order_measures_the_curve_from_its_skylines_only() -> None:
    # ARRANGE
    # A CONTINUUM AT ABOUT 20 SIGMA PER PIXEL: EVERY PIXEL PASSES THE MIN_SKY_SNR FLOOR, AS IN THE THERMAL NIR ORDER 10.
    # THE CONTINUUM HAS A STEEP OPPOSITE GRADIENT; THE LINES HAVE THE ILLUMINATION GRADIENT AND ARE BRIGHT ENOUGH THAT
    # THE CONTINUUM UNDER THEIR 30 SIGMA PIXELS CARRIES LITTLE OF THE WEIGHT
    rng = np.random.default_rng(5)
    order = _synthetic_order(gradientPerArcsec=0.0, lineAmplitude=60000.0)
    # CONTINUUM IS THE SYNTHETIC ORDER'S OWN BASELINE; IT IS REPLACED BY A BRIGHTER THERMAL CONTINUUM
    thermalContinuum = 400.0
    lines = order.skyModel - CONTINUUM
    skyModel = thermalContinuum + lines
    error = np.sqrt(skyModel + READ_NOISE**2)
    data = (
        thermalContinuum * (1.0 - 0.02 * order.slitPosition)
        + lines * (1.0 + GRADIENT_PER_ARCSEC * order.slitPosition)
        + rng.normal(0.0, 1.0, skyModel.size) * error
    )
    assert np.mean(skyModel >= MIN_SKY_SNR * error) > 0.9
    assert np.median(skyModel / error) < CONTINUUM_SKY_SNR

    # ACT
    correction = measure_slit_illumination(order.slitPosition, order.wavelength, skyModel, data, error, order.isUsable)

    # ASSERT
    assert correction.isApplied
    assert correction.skySnrFloor == CONTINUUM_SKY_SNR
    slope = np.polyfit(correction.binCentres, correction.binValues, 1)[0]
    assert slope == pytest.approx(GRADIENT_PER_ARCSEC, abs=0.0015)


def test_a_line_dominated_order_uses_the_minimum_signal_to_noise_floor() -> None:
    # ARRANGE
    order = _synthetic_order()

    # ACT
    correction = _measure(order)

    # ASSERT
    assert correction.skySnrFloor == MIN_SKY_SNR


def test_the_fill_zone_takes_a_quadratic_fitted_to_the_bins_across_it() -> None:
    # ARRANGE
    # A CONCAVE-UP ILLUMINATION, LOWEST AT THE SLIT CENTRE, AS ON THE NIR SLIT1.0 TARGET 3FGL
    curvature = 0.004
    order = _synthetic_order(gradientPerArcsec=0.0)
    data = order.data * (1.0 + curvature * order.slitPosition**2)
    isZone = np.abs(order.slitPosition) <= 1.5

    # ACT
    filled = _measure(order, data=data, isUsable=~isZone, isFillZone=isZone)
    interpolated = _measure(order, data=data, isUsable=~isZone)

    # ASSERT
    shape = curvature * order.slitPosition**2
    outside = ~isZone & (np.abs(order.slitPosition) < 4.5)
    offset = np.mean(filled.factor[outside] - shape[outside])
    np.testing.assert_allclose(filled.factor[isZone], offset + shape[isZone], atol=BIN_TOLERANCE)
    np.testing.assert_array_equal(filled.factor[~isZone], interpolated.factor[~isZone])
    # LINEAR INTERPOLATION ACROSS THE ZONE STAYS AT THE LEVEL OF ITS EDGES, ABOUT 0.009 ABOVE THE CURVE AT THE CENTRE
    isCentre = np.abs(order.slitPosition) < 0.3
    assert interpolated.factor[isCentre].mean() > filled.factor[isCentre].mean() + 0.006


def test_the_fill_zone_is_interpolated_when_too_few_bins_lie_near_the_slit_centre() -> None:
    # ARRANGE
    slit, wavelength, skyModel, data, error, isUsable = _dense_bins([0.0, 0.3, 0.8], PIXELS_AT_MAX_BIN_ERROR)
    # THE FAINT PIXELS OF _dense_bins MOVE INTO THE ZONE: NOT USABLE, SO THEY TAKE THEIR FACTOR FROM THE CURVE. THE
    # USABLE PIXELS ARE THEN MOSTLY BRIGHT, SO THE CONTINUUM FLOOR APPLIES: THE SKYLINES ARE MADE 40 SIGMA TO PASS IT
    isZone = skyModel == 1.0
    slit = np.where(isZone, 0.55, slit)
    error = error / 4.0

    # ACT
    filled = measure_slit_illumination(slit, wavelength, skyModel, data, error, isUsable & ~isZone, isFillZone=isZone)

    # ASSERT
    assert filled.binCentres.size == 3
    expected = 1.0 + np.interp(0.55, filled.binCentres, filled.binValues)
    np.testing.assert_allclose(filled.factor[isZone], expected)


def test_the_fill_zone_keeps_the_interpolated_curve_beyond_the_outermost_fitted_bins() -> None:
    # ARRANGE
    order = _synthetic_order(gradientPerArcsec=0.0)
    data = order.data * (1.0 + 0.004 * order.slitPosition**2)
    # AN OBJECT AT THE SLIT END: THE ZONE REACHES PAST THE LAST USABLE PIXEL, WHERE A QUADRATIC WOULD BE EXTRAPOLATED
    isUsable = order.slitPosition <= 3.5
    isZone = order.slitPosition > 3.5

    # ACT
    filled = _measure(order, data=data, isUsable=isUsable, isFillZone=isZone)
    interpolated = _measure(order, data=data, isUsable=isUsable)

    # ASSERT
    assert filled.binCentres.max() < 3.5
    np.testing.assert_array_equal(filled.factor, interpolated.factor)


def test_a_skyline_pixel_with_an_overflowing_weight_is_left_out() -> None:
    # ARRANGE
    order = _synthetic_order()
    error = order.error.copy()
    isTiny = order.skyModel > 1000.0
    tinyIndex = np.flatnonzero(isTiny)[::50]
    error[tinyIndex] = 1e-200

    # ACT
    clean = _measure(order)
    withTiny = measure_slit_illumination(
        order.slitPosition, order.wavelength, order.skyModel, order.data, error, order.isUsable
    )

    # ASSERT
    assert np.isfinite(withTiny.factor).all()
    assert withTiny.pixelCount < clean.pixelCount
    np.testing.assert_allclose(withTiny.factor, clean.factor, atol=0.001)


def test_inputs_of_different_lengths_raise_a_value_error() -> None:
    # ARRANGE
    order = _synthetic_order()

    # ACT / ASSERT
    with pytest.raises(ValueError, match="same length"):
        measure_slit_illumination(
            order.slitPosition, order.wavelength[:-1], order.skyModel, order.data, order.error, order.isUsable
        )
