"""The Kelson-style skyline slit-illumination correction measures a multiplicative factor along the slit (DY-2054)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pytest

from soxspipe.commonutils.slit_illumination_correction import (
    BIN_WIDTH_ARCSEC,
    MIN_PIXELS_PER_SLIT_BIN,
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


def test_factor_is_one_when_the_skylines_give_too_few_pixels_for_each_slit_bin() -> None:
    # ARRANGE
    order = _synthetic_order(lineSpacing=100)

    # ACT
    correction = _measure(order)

    # ASSERT
    pixelsPerBin = correction.pixelCount / (np.ptp(order.slitPosition) / BIN_WIDTH_ARCSEC)
    assert 20 < pixelsPerBin < MIN_PIXELS_PER_SLIT_BIN
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


def _dense_bins(positions: list[float], pixelsPerPosition: int):
    """Skyline pixels stacked at the given slit positions, plus two usable pixels that mark the wavelength range."""
    slit = np.concatenate([np.full(pixelsPerPosition, position) for position in positions] + [[positions[0]] * 2])
    wavelength = np.full(slit.size, 1550.0)
    wavelength[-2:] = [1500.0, 1600.0]
    skyModel = np.full(slit.size, 1000.0)
    skyModel[-2:] = 0.0
    error = np.ones(slit.size)
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


def test_a_bin_with_exactly_the_minimum_pixels_is_not_used() -> None:
    # ARRANGE
    positions = [0.0, 0.3, 0.55, 0.8]

    # ACT
    atMinimum = measure_slit_illumination(*_dense_bins(positions, MIN_PIXELS_PER_SLIT_BIN))
    aboveMinimum = measure_slit_illumination(*_dense_bins(positions, MIN_PIXELS_PER_SLIT_BIN + 1))

    # ASSERT
    assert not atMinimum.isApplied
    assert aboveMinimum.isApplied


def test_a_pixel_at_the_far_edge_of_the_slit_falls_in_the_last_bin() -> None:
    # ARRANGE
    slit, *rest = _dense_bins([0.0, 0.3, 0.55, 0.8, 1.0], MIN_PIXELS_PER_SLIT_BIN + 1)

    # ACT
    correction = measure_slit_illumination(slit, *rest)

    # ASSERT
    assert correction.binCentres.size == 4
    assert correction.binCounts[-1] == 2 * (MIN_PIXELS_PER_SLIT_BIN + 1)
    assert correction.binCentres[-1] == pytest.approx(0.875)


def test_inputs_of_different_lengths_raise_a_value_error() -> None:
    # ARRANGE
    order = _synthetic_order()

    # ACT / ASSERT
    with pytest.raises(ValueError, match="same length"):
        measure_slit_illumination(
            order.slitPosition, order.wavelength[:-1], order.skyModel, order.data, order.error, order.isUsable
        )
