"""Unit tests for the fixed log-wavelength grid and the flux-conserving resampler (DY-2348).

The grid points are lambda_k = 1000 nm * exp(k * dv / c), so every bin has the same width in km/s. The resampler
moves a density (flat-fielded electrons per pixel) onto the grid with normalized overlap weights, and the variance
of each bin is the sum of the squared weights times the input variances.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.log_wavelength_grid import (
    LOG_GRID_ORIGIN_NM,
    SPEED_OF_LIGHT_KMS,
    bin_edges,
    covered_bin_range,
    grid_index,
    grid_step_kms,
    grid_wavelength,
    overlap_weights,
    snap_to_edge_index,
)

pytestmark = pytest.mark.unit

STEP_KMS = 20.0


def test_grid_wavelength_is_the_origin_at_index_zero_and_grows_by_the_velocity_step() -> None:
    # ACT
    wavelength = grid_wavelength(np.array([0, 1, -1]), STEP_KMS)

    # ASSERT
    assert LOG_GRID_ORIGIN_NM == 1000.0
    np.testing.assert_allclose(wavelength[0], 1000.0, rtol=1e-15)
    np.testing.assert_allclose(wavelength[1] / 1000.0, np.exp(STEP_KMS / SPEED_OF_LIGHT_KMS), rtol=1e-15)
    np.testing.assert_allclose(wavelength[2] * wavelength[1], 1000.0**2, rtol=1e-15)


def test_grid_wavelength_is_the_same_for_every_caller() -> None:
    # ARRANGE: TWO FRAMES THAT SPAN DIFFERENT RANGES MUST SHARE EVERY GRID POINT THEY BOTH HOLD
    first = grid_wavelength(np.arange(100, 200), STEP_KMS)
    second = grid_wavelength(np.arange(150, 260), STEP_KMS)

    # ASSERT
    np.testing.assert_array_equal(first[50:], second[:50])


def test_bin_edges_hold_every_bin_at_the_same_width_in_velocity() -> None:
    # ACT
    edges = bin_edges(10, 14, STEP_KMS)

    # ASSERT: 5 BINS NEED 6 EDGES, AND EACH EDGE IS THE PREVIOUS ONE TIMES THE SAME FACTOR
    assert edges.shape == (6,)
    np.testing.assert_allclose(edges[1:] / edges[:-1], np.exp(STEP_KMS / SPEED_OF_LIGHT_KMS), rtol=1e-14)


def test_grid_points_are_the_geometric_centres_of_their_bins() -> None:
    # ACT
    edges = bin_edges(10, 14, STEP_KMS)
    centres = grid_wavelength(np.arange(10, 15), STEP_KMS)

    # ASSERT
    np.testing.assert_allclose(centres, np.sqrt(edges[:-1] * edges[1:]), rtol=1e-14)


def test_grid_step_multiplies_the_per_unbinned_pixel_step_by_the_binning() -> None:
    # ARRANGE
    detectorParams = {"log-wavelength-step-kms": 17.83}

    # ACT AND ASSERT
    assert grid_step_kms(detectorParams, 1) == pytest.approx(17.83)
    assert grid_step_kms(detectorParams, 2) == pytest.approx(35.66)


@pytest.mark.parametrize("badValue", [0, -3.0, float("nan"), float("inf"), "17.83", None, True])
def test_grid_step_refuses_a_value_that_is_not_a_positive_number(badValue: object) -> None:
    with pytest.raises(ValueError, match="log-wavelength-step-kms"):
        grid_step_kms({"log-wavelength-step-kms": badValue}, 1)


def test_grid_step_refuses_detector_parameters_without_the_key() -> None:
    with pytest.raises(ValueError, match="log-wavelength-step-kms"):
        grid_step_kms({"ron": 3.0}, 1)


@pytest.mark.parametrize("stepKms", [13.16, 17.83, 26.32, 35.66])
def test_grid_index_recovers_the_index_of_every_grid_point_exactly(stepKms: float) -> None:
    # ARRANGE
    indices = np.arange(-300, 2500)

    # ACT
    recovered = grid_index(grid_wavelength(indices, stepKms), stepKms)

    # ASSERT
    np.testing.assert_array_equal(recovered, indices)
    assert recovered.dtype.kind == "i"


def test_grid_index_gives_the_bin_that_holds_a_wavelength_off_the_grid_point() -> None:
    # ARRANGE: JUST INSIDE EACH EDGE OF BIN 40, AND JUST OUTSIDE IT
    lowerEdge, upperEdge = bin_edges(40, 40, STEP_KMS)
    inside = np.array([lowerEdge * (1 + 1e-9), upperEdge * (1 - 1e-9)])
    outside = np.array([lowerEdge * (1 - 1e-9), upperEdge * (1 + 1e-9)])

    # ACT / ASSERT
    np.testing.assert_array_equal(grid_index(inside, STEP_KMS), [40, 40])
    np.testing.assert_array_equal(grid_index(outside, STEP_KMS), [39, 41])


def test_grid_index_ignores_a_shift_of_a_grid_point_by_the_old_four_decimal_rounding() -> None:
    # ARRANGE: THE LARGEST SHIFT THE OLD ROUNDING GAVE A WAVELENGTH IS 5E-5 NM
    wavelength = grid_wavelength(np.arange(0, 800), STEP_KMS)

    # ACT / ASSERT
    np.testing.assert_array_equal(grid_index(wavelength + 5e-5, STEP_KMS), np.arange(0, 800))
    np.testing.assert_array_equal(grid_index(wavelength - 5e-5, STEP_KMS), np.arange(0, 800))


@pytest.mark.parametrize("badValue", [np.nan, np.inf, 0.0, -500.0])
def test_grid_index_refuses_a_wavelength_that_is_not_a_positive_finite_number(badValue: float) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match="positive, finite"):
        grid_index(np.array([500.0, badValue]), STEP_KMS)


def test_snap_to_edge_index_picks_the_nearest_edge_in_log_wavelength() -> None:
    # ARRANGE: EDGE m SITS AT 1000 * exp((m - 0.5) * dv / c)
    edges = bin_edges(3, 6, STEP_KMS)

    # ACT
    nearJustAbove = snap_to_edge_index(edges[2] * (1 + 0.2 * STEP_KMS / SPEED_OF_LIGHT_KMS), STEP_KMS)
    nearJustBelow = snap_to_edge_index(edges[2] * (1 - 0.2 * STEP_KMS / SPEED_OF_LIGHT_KMS), STEP_KMS)
    nearNext = snap_to_edge_index(edges[3] * (1 - 0.4 * STEP_KMS / SPEED_OF_LIGHT_KMS), STEP_KMS)

    # ASSERT: EDGES 3..7 ARE edges[0]..edges[4]
    assert nearJustAbove == 5
    assert nearJustBelow == 5
    assert nearNext == 6


def test_covered_bin_range_leaves_out_the_partly_covered_end_bins() -> None:
    # ARRANGE: THE DATA RUN FROM A QUARTER INTO BIN 20 TO THREE QUARTERS OF THE WAY THROUGH BIN 30
    edges = bin_edges(20, 30, STEP_KMS)
    lower = edges[0] * (1 + 0.25 * STEP_KMS / SPEED_OF_LIGHT_KMS)
    upper = edges[-2] * (1 + 0.75 * STEP_KMS / SPEED_OF_LIGHT_KMS)

    # ACT
    first, last = covered_bin_range(lower, upper, STEP_KMS)

    # ASSERT: BIN 20 IS ONLY PARTLY COVERED, BIN 30 STARTS AT edges[-2] AND ENDS BEYOND THE DATA
    assert (first, last) == (21, 29)


def test_covered_bin_range_keeps_a_bin_whose_edges_match_the_data_ends() -> None:
    # ARRANGE
    edges = bin_edges(20, 30, STEP_KMS)

    # ACT
    first, last = covered_bin_range(edges[0], edges[-1], STEP_KMS)

    # ASSERT
    assert (first, last) == (20, 30)


def test_covered_bin_range_is_empty_when_no_bin_fits_inside_the_data() -> None:
    # ARRANGE
    edges = bin_edges(20, 21, STEP_KMS)

    # ACT
    first, last = covered_bin_range(edges[0] * 1.00001, edges[1] * 0.99999, STEP_KMS)

    # ASSERT
    assert last < first


def _pixels(width: float, count: int, start: float) -> tuple[np.ndarray, np.ndarray]:
    """Edge-to-edge pixels of equal width: their lower and upper edges."""
    lower = start + width * np.arange(count)
    return lower, lower + width


def test_resampling_a_constant_density_keeps_the_constant() -> None:
    # ARRANGE: PIXELS 0.037 NM WIDE RESAMPLED ONTO A GRID OF A DIFFERENT, NON-COMMENSURATE WIDTH
    lower, upper = _pixels(0.037, 400, 1000.0)
    first, last = covered_bin_range(lower[0], upper[-1], 14.0)
    edges = bin_edges(first, last, 14.0)
    weights = overlap_weights(lower, upper, edges)

    # ACT
    resampled = weights.mean(np.full(400, 123.5))

    # ASSERT
    assert weights.hasData.all()
    np.testing.assert_allclose(resampled, 123.5, rtol=1e-12)


def test_resampling_keeps_the_integrated_flux() -> None:
    # ARRANGE: A DENSITY THAT VARIES FROM PIXEL TO PIXEL, ON PIXELS OF UNEQUAL WIDTH
    rng = np.random.default_rng(7)
    widths = rng.uniform(0.03, 0.05, 300)
    upper = 1000.0 + np.cumsum(widths)
    lower = upper - widths
    density = rng.uniform(50.0, 150.0, 300)
    first, last = covered_bin_range(lower[0], upper[-1], 15.0)
    edges = bin_edges(first, last, 15.0)
    weights = overlap_weights(lower, upper, edges)

    # ACT
    resampled = weights.mean(density)

    # ASSERT: THE FLUX IS THE DENSITY INTEGRATED OVER WAVELENGTH, AND THE PIXELS THE BINS DO NOT COVER ARE LEFT OUT
    binFlux = (resampled * np.diff(edges)).sum()
    pixelOverlap = np.clip(np.minimum(upper, edges[-1]) - np.maximum(lower, edges[0]), 0, None)
    assert binFlux == pytest.approx((density * pixelOverlap).sum(), rel=1e-12)


def test_resampling_weights_each_pixel_by_its_overlap_with_the_bin() -> None:
    # ARRANGE: ONE BIN SPANNING 3/4 OF PIXEL 0 AND 1/4 OF PIXEL 1
    binEdge = bin_edges(0, 0, 10.0)
    width = binEdge[1] - binEdge[0]
    lower = np.array([binEdge[0] - 0.25 * width, binEdge[0] + 0.75 * width])
    upper = np.array([binEdge[0] + 0.75 * width, binEdge[0] + 1.75 * width])
    weights = overlap_weights(lower, upper, binEdge)

    # ACT
    value = weights.mean(np.array([100.0, 200.0]))

    # ASSERT
    np.testing.assert_allclose(value, [125.0], rtol=1e-12)


def test_resampled_variance_is_the_sum_of_the_squared_weights_times_the_variances() -> None:
    # ARRANGE: THE SAME BIN, WEIGHTS 3/4 AND 1/4
    binEdge = bin_edges(0, 0, 10.0)
    width = binEdge[1] - binEdge[0]
    lower = np.array([binEdge[0] - 0.25 * width, binEdge[0] + 0.75 * width])
    upper = np.array([binEdge[0] + 0.75 * width, binEdge[0] + 1.75 * width])
    weights = overlap_weights(lower, upper, binEdge)

    # ACT
    variance = weights.propagate_variance(np.array([16.0, 64.0]))

    # ASSERT: (3/4)^2 * 16 + (1/4)^2 * 64 = 9 + 4
    np.testing.assert_allclose(variance, [13.0], rtol=1e-12)


def test_a_bin_inside_one_pixel_takes_that_pixel_and_its_variance() -> None:
    # ARRANGE: THE PIXEL IS 4 BINS WIDE
    edges = bin_edges(0, 3, 10.0)
    weights = overlap_weights(np.array([edges[0]]), np.array([edges[-1]]), edges)

    # ACT AND ASSERT
    np.testing.assert_allclose(weights.mean(np.array([7.0])), [7.0] * 4)
    np.testing.assert_allclose(weights.propagate_variance(np.array([9.0])), [9.0] * 4)


def test_a_bin_with_no_pixel_has_no_data_and_no_interpolated_value() -> None:
    # ARRANGE: TWO RUNS OF PIXELS WITH A GAP OF SEVERAL BINS BETWEEN THEM
    edges = bin_edges(0, 19, 10.0)
    lower = np.array([edges[0], edges[14]])
    upper = np.array([edges[4], edges[20]])
    weights = overlap_weights(lower, upper, edges)

    # ACT
    value = weights.mean(np.array([100.0, 300.0]))

    # ASSERT: BINS 4 TO 13 LIE IN THE GAP
    assert weights.hasData.tolist() == [True] * 4 + [False] * 10 + [True] * 6
    np.testing.assert_allclose(value[:4], 100.0)
    np.testing.assert_allclose(value[14:], 300.0)
    assert (value[4:14] == 0.0).all()
    assert (weights.propagate_variance(np.array([1.0, 1.0]))[4:14] == 0.0).all()


def test_a_dropped_pixel_leaves_a_hole_that_the_other_pixels_do_not_fill() -> None:
    # ARRANGE: THREE PIXELS, EACH ONE BIN WIDE, AND THE MIDDLE ONE IS DROPPED
    edges = bin_edges(0, 2, 10.0)
    lower, upper = edges[:-1], edges[1:]
    weights = overlap_weights(lower, upper, edges, keep=np.array([True, False, True]))

    # ACT
    value = weights.mean(np.array([10.0, 99999.0, 30.0]))

    # ASSERT
    assert weights.hasData.tolist() == [True, False, True]
    np.testing.assert_allclose(value, [10.0, 0.0, 30.0])


def test_a_bin_that_loses_only_some_pixels_keeps_the_others_as_a_weighted_mean() -> None:
    # ARRANGE: ONE BIN MADE OF TWO HALF-BIN PIXELS, THE SECOND DROPPED
    binEdge = bin_edges(0, 0, 10.0)
    middle = np.sqrt(binEdge[0] * binEdge[1])
    weights = overlap_weights(
        np.array([binEdge[0], middle]), np.array([middle, binEdge[1]]), binEdge, keep=np.array([True, False])
    )

    # ACT
    value = weights.mean(np.array([40.0, 1.0e9]))
    variance = weights.propagate_variance(np.array([4.0, 1.0e9]))

    # ASSERT
    assert weights.hasData.tolist() == [True]
    np.testing.assert_allclose(value, [40.0])
    np.testing.assert_allclose(variance, [4.0])


def test_nan_in_a_kept_pixel_reaches_only_the_bins_that_pixel_overlaps() -> None:
    # ARRANGE: SEE test_a_dropped_pixel_leaves_a_hole... BUT THE MIDDLE VALUE IS NAN AND THE PIXEL IS KEPT
    edges = bin_edges(0, 2, 10.0)
    weights = overlap_weights(edges[:-1], edges[1:], edges)

    # ACT
    value = weights.mean(np.array([10.0, np.nan, 30.0]))

    # ASSERT
    np.testing.assert_allclose(value[[0, 2]], [10.0, 30.0])
    assert np.isnan(value[1])


@pytest.mark.parametrize("badEdge", [np.nan, np.inf])
def test_a_pixel_with_a_non_finite_edge_is_treated_as_dropped(badEdge: float) -> None:
    # ARRANGE
    edges = bin_edges(0, 2, 10.0)
    lower, upper = edges[:-1].copy(), edges[1:].copy()
    upper[1] = badEdge

    # ACT
    weights = overlap_weights(lower, upper, edges)

    # ASSERT
    assert weights.hasData.tolist() == [True, False, True]


def test_pixels_beyond_the_bins_are_ignored() -> None:
    # ARRANGE: PIXELS RUN WELL PAST BOTH ENDS OF THE BINS
    edges = bin_edges(5, 7, 10.0)
    lower = np.array([edges[0] * 0.99, edges[0], edges[-1]])
    upper = np.array([edges[0], edges[-1], edges[-1] * 1.01])
    weights = overlap_weights(lower, upper, edges)

    # ACT
    value = weights.mean(np.array([1.0e6, 5.0, 1.0e6]))

    # ASSERT
    np.testing.assert_allclose(value, [5.0, 5.0, 5.0])
