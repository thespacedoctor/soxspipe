"""The Kelson-style two-pass offset-bin pre-filter of the sky pixels (DY-2055)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.sky_prefilter import (
    OFFSET_BIN_SIGMA,
    offset_bin_pass_masks,
    order_dispersion,
    robust_sigma,
    two_pass_offset_bin_outliers,
)

pytestmark = pytest.mark.unit

COLUMNS = 40
SLIT_PIXELS = 30
# ONE COLUMN IS ONE WAVELENGTH UNIT, SO A BIN OF FIVE COLUMNS IS FIVE UNITS WIDE
BIN_WIDTH = 5.0
SKY_LEVEL = 100.0
# BOUNDED NOISE: A UNIFORM +-1 HAS A SIGMA OF 0.58, SO NO PIXEL EVER SITS 3 SIGMA FROM ITS BIN MEDIAN
NOISE_HALF_RANGE = 1.0
LINE_AMPLITUDE = 200.0
SPIKE_AMPLITUDE = 500.0


def _sky_grid(seed: int = 3) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A COLUMNS x SLIT_PIXELS block of sky: the wavelength, the flux and the column of every pixel."""
    rng = np.random.default_rng(seed)
    column = np.repeat(np.arange(COLUMNS), SLIT_PIXELS)
    wavelength = column.astype(float)
    flux = SKY_LEVEL + rng.uniform(-NOISE_HALF_RANGE, NOISE_HALF_RANGE, column.size)
    return wavelength, flux, column


def test_a_narrow_line_straddling_a_bin_edge_is_kept() -> None:
    # ARRANGE: A THREE-COLUMN LINE ACROSS THE EDGE AT 15; EACH FIRST-PASS BIN HOLDS A MINORITY OF IT, THE SECOND-PASS
    # BIN CENTRED ON IT HOLDS A MAJORITY
    wavelength, flux, column = _sky_grid()
    isLine = np.isin(column, (14, 15, 16))
    flux[isLine] += LINE_AMPLITUDE

    # ACT
    firstPass, secondPass = offset_bin_pass_masks(wavelength, flux, BIN_WIDTH)
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT: THE FIRST PASS ALONE WOULD HAVE CUT THE LINE
    assert firstPass[isLine].all()
    assert not rejected[isLine].any()


def test_a_cosmic_ray_spike_is_rejected_and_its_neighbours_are_not() -> None:
    # ARRANGE
    wavelength, flux, column = _sky_grid()
    spike = 17 * SLIT_PIXELS + 11
    flux[spike] += SPIKE_AMPLITUDE

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert np.flatnonzero(rejected).tolist() == [spike]


def test_a_negative_spike_is_rejected_because_the_filter_is_two_sided() -> None:
    # ARRANGE
    wavelength, flux, _ = _sky_grid()
    cold = 5 * SLIT_PIXELS + 2
    flux[cold] -= SPIKE_AMPLITUDE

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert np.flatnonzero(rejected).tolist() == [cold]


def test_only_pixels_flagged_in_both_passes_are_rejected() -> None:
    # ARRANGE: THE STRADDLING LINE IS FLAGGED IN THE FIRST PASS ONLY, AND A SPIKE IN BOTH
    wavelength, flux, column = _sky_grid()
    isLine = np.isin(column, (14, 15, 16))
    flux[isLine] += LINE_AMPLITUDE
    spike = 30 * SLIT_PIXELS + 4
    flux[spike] += SPIKE_AMPLITUDE
    expectedRejected = np.zeros(wavelength.size, dtype=bool)
    expectedRejected[spike] = True

    # ACT
    firstPass, secondPass = offset_bin_pass_masks(wavelength, flux, BIN_WIDTH)
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert firstPass[isLine].all()
    assert not secondPass[isLine].any()
    assert firstPass[spike] and secondPass[spike]
    assert np.array_equal(rejected, expectedRejected)


def test_the_masks_follow_the_row_order_of_the_input() -> None:
    # ARRANGE
    wavelength, flux, _ = _sky_grid()
    spike = 17 * SLIT_PIXELS + 11
    flux[spike] += SPIKE_AMPLITUDE
    shuffle = np.random.default_rng(9).permutation(wavelength.size)

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength[shuffle], flux[shuffle], BIN_WIDTH)

    # ASSERT
    assert np.flatnonzero(rejected).tolist() == [int(np.flatnonzero(shuffle == spike)[0])]


def test_non_finite_pixels_are_never_rejected_and_do_not_set_the_bin_statistics() -> None:
    # ARRANGE
    wavelength, flux, _ = _sky_grid()
    flux[::7] = np.nan
    wavelength[3::11] = np.nan
    spike = 17 * SLIT_PIXELS + 11
    flux[spike] = SKY_LEVEL + SPIKE_AMPLITUDE

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert rejected[spike]
    assert not rejected[~np.isfinite(flux) | ~np.isfinite(wavelength)].any()


def test_a_bin_with_too_few_pixels_is_not_tested() -> None:
    # ARRANGE: THREE PIXELS IN ONE BIN, ONE OF THEM FAR OFF
    wavelength = np.array([0.0, 0.1, 0.2])
    flux = np.array([100.0, 101.0, 900.0])

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert not rejected.any()


def test_a_bin_whose_scale_is_zero_is_skipped_even_with_one_outlier() -> None:
    # ARRANGE: 60 IDENTICAL PIXELS AND ONE FAR OFF. THE MAD AND THE BIWEIGHT SCALE ARE ZERO, SO THE BIN IS NOT TESTED;
    # THE
    # IDL WOULD FLAG EVERY PIXEL OF IT ("GE" A ZERO LIMIT). THE SPIKE IN THE NEXT, NOISY BIN IS STILL REJECTED
    wavelength = np.repeat(np.arange(20, dtype=float), 12)
    flux = np.full(wavelength.size, SKY_LEVEL)
    outlier = 12 * 1 + 3
    flux[outlier] += SPIKE_AMPLITUDE
    noisy = wavelength >= 10
    flux[noisy] += np.tile(np.linspace(-NOISE_HALF_RANGE, NOISE_HALF_RANGE, 12), 10)
    spike = 15 * 12 + 5
    flux[spike] += SPIKE_AMPLITUDE
    expected = np.zeros(wavelength.size, dtype=bool)
    expected[spike] = True

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert np.array_equal(rejected, expected)


def test_flat_data_rejects_nothing() -> None:
    # ARRANGE: A SIGMA OF ZERO MUST NOT MAKE EVERY PIXEL "AT LEAST 3 SIGMA" FROM THE MEDIAN
    wavelength = np.repeat(np.arange(20, dtype=float), 12)
    flux = np.full(wavelength.size, SKY_LEVEL)

    # ACT
    rejected = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH)

    # ASSERT
    assert not rejected.any()


def test_an_empty_order_returns_empty_masks() -> None:
    # ACT
    rejected = two_pass_offset_bin_outliers(np.array([]), np.array([]), BIN_WIDTH)

    # ASSERT
    assert rejected.shape == (0,)
    assert rejected.dtype == bool


def test_a_wider_sigma_limit_rejects_fewer_pixels() -> None:
    # ARRANGE
    rng = np.random.default_rng(4)
    wavelength = np.repeat(np.arange(COLUMNS, dtype=float), SLIT_PIXELS)
    flux = SKY_LEVEL + rng.normal(0.0, 1.0, wavelength.size)

    # ACT
    tight = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH, sigma=2.0)
    wide = two_pass_offset_bin_outliers(wavelength, flux, BIN_WIDTH, sigma=OFFSET_BIN_SIGMA)

    # ASSERT
    assert 0 < wide.sum() < tight.sum()


def test_robust_sigma_recovers_the_scatter_of_gaussian_noise() -> None:
    # ARRANGE
    values = np.random.default_rng(1).normal(50.0, 4.0, 5000)

    # ACT
    sigma = robust_sigma(values)

    # ASSERT
    assert sigma == pytest.approx(4.0, rel=0.05)


def test_robust_sigma_ignores_a_minority_of_far_outliers() -> None:
    # ARRANGE
    values = np.random.default_rng(1).normal(50.0, 4.0, 5000)
    values[:800] += 400.0

    # ACT
    sigma = robust_sigma(values)

    # ASSERT: THE IDL FORMULA SCALES BY THE COUNT OF ALL VALUES, SO 16% OUTLIERS LEAVE IT ABOUT 12% HIGH; THE STANDARD
    # DEVIATION OF THESE VALUES IS ABOUT 150
    assert sigma == pytest.approx(4.0, rel=0.15)


def test_robust_sigma_of_constant_values_is_zero() -> None:
    assert robust_sigma(np.full(20, 7.0)) == 0.0


def test_order_dispersion_is_the_wavelength_step_of_one_dispersion_pixel() -> None:
    # ARRANGE: 100 COLUMNS RUN 1000 TO 1099 NM ALONG X; THE SLIT FILLS 12 ROWS
    x = np.repeat(np.arange(100), 12)
    y = np.tile(np.arange(12), 100)
    wavelength = 1000.0 + x.astype(float)

    # ACT
    step = order_dispersion(wavelength, x, y)

    # ASSERT
    assert step == pytest.approx(1.0)


def test_order_dispersion_finds_a_dispersion_axis_along_y() -> None:
    # ARRANGE: THE SAME ORDER ROTATED SO THAT WAVELENGTH RUNS ALONG Y AT 0.5 NM PER PIXEL
    y = np.repeat(np.arange(100), 12)
    x = np.tile(np.arange(12), 100)
    wavelength = 500.0 + 0.5 * y

    # ACT
    step = order_dispersion(wavelength, x, y)

    # ASSERT
    assert step == pytest.approx(0.5)


def test_a_slit_tilt_over_estimates_the_dispersion_by_the_tilt_span_over_the_order_length() -> None:
    # ARRANGE: THE WAVELENGTH ALSO RISES ALONG THE SLIT, BY 5.5 NM OVER 12 ROWS, SO THE SPAN GROWS BY 5.5 OVER 99
    # INTERVALS
    x = np.repeat(np.arange(100), 12)
    y = np.tile(np.arange(12), 100)
    wavelength = 1000.0 + x + 0.5 * y

    # ACT
    step = order_dispersion(wavelength, x, y)

    # ASSERT: THE ACCEPTED ERROR OF THE RULE, ABOUT +5.6%
    assert step == pytest.approx(1.0 + 5.5 / 99, rel=1e-9)


def test_a_column_gap_over_estimates_the_dispersion_by_the_gap_over_the_pixels_present() -> None:
    # ARRANGE: COLUMNS 50 TO 59 ARE ABSENT (MASKED OUT), SO 100 COLUMNS SPAN 109 NM
    columns = np.r_[np.arange(0, 50), np.arange(60, 110)]
    x = np.repeat(columns, 12)
    y = np.tile(np.arange(12), columns.size)
    wavelength = 1000.0 + x.astype(float)

    # ACT
    step = order_dispersion(wavelength, x, y)

    # ASSERT: THE ACCEPTED ERROR OF THE RULE, ABOUT +10%
    assert step == pytest.approx(109 / 99, rel=1e-9)


def test_binned_coordinates_give_the_wavelength_step_of_one_binned_pixel() -> None:
    # ARRANGE: 2X BINNED, SO COLUMN INDICES RUN 0 TO 99 AND EACH IS 2 NM
    x = np.repeat(np.arange(100), 6)
    y = np.tile(np.arange(6), 100)
    wavelength = 1000.0 + 2.0 * x

    # ACT
    step = order_dispersion(wavelength, x, y)

    # ASSERT
    assert step == pytest.approx(2.0)


def test_order_dispersion_of_a_single_pixel_order_is_zero() -> None:
    assert order_dispersion(np.array([500.0]), np.array([3]), np.array([4])) == 0.0
