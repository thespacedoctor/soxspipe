"""Unit tests for the order merge on the log-wavelength grid (DY-2348).

Each order is cut at its join with the next order, the join is snapped to a bin edge, and a bin that an order covers
only in part is dropped. There is no interpolation: a bin that no pixel overlaps has no data and flux 0.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.log_grid_merge import OrderPixels, merge_orders_on_log_grid
from soxspipe.commonutils.log_wavelength_grid import SPEED_OF_LIGHT_KMS, bin_edges, grid_wavelength

pytestmark = pytest.mark.unit

STEP_KMS = 13.16
PIXEL_NM = 0.05


def _order(
    start: float,
    count: int,
    value: float,
    variance: float = 4.0,
    keep: np.ndarray | None = None,
    sky: float | None = None,
) -> OrderPixels:
    wavelength = start + PIXEL_NM * np.arange(count)
    values = {"flux": np.full(count, value)}
    if sky is not None:
        values["sky"] = np.full(count, sky)
    return OrderPixels(
        wavelength=wavelength,
        pixelScale=np.full(count, PIXEL_NM),
        keep=np.ones(count, dtype=bool) if keep is None else keep,
        values=values,
        variance=np.full(count, variance),
    )


def test_one_order_gives_every_bin_it_covers_wholly_and_no_other() -> None:
    # ARRANGE
    order = _order(1000.0, 200, value=50.0)

    # ACT
    merged = merge_orders_on_log_grid([order], [], STEP_KMS)

    # ASSERT: THE DATA SPAN 999.975 TO 1009.975 NM, AND EVERY BIN INSIDE IT HAS THE ORDER'S VALUE
    assert merged.hasData.all()
    np.testing.assert_allclose(merged.values["flux"], 50.0, rtol=1e-12)
    assert merged.wavelength[0] > 999.975
    assert merged.wavelength[-1] < 1009.975
    edges = bin_edges(merged.firstIndex, merged.firstIndex + len(merged.wavelength) - 1, STEP_KMS)
    assert edges[0] >= 999.975 - 1e-9
    assert edges[-1] <= 1009.975 + 1e-9
    assert edges[0] < 999.975 + 0.05
    assert edges[-1] > 1009.975 - 0.05


def test_merged_wavelengths_are_grid_points() -> None:
    # ACT
    merged = merge_orders_on_log_grid([_order(1000.0, 200, 50.0)], [], STEP_KMS)

    # ASSERT
    expected = grid_wavelength(merged.firstIndex + np.arange(len(merged.wavelength)), STEP_KMS)
    np.testing.assert_array_equal(merged.wavelength, expected)


def test_each_bin_of_two_overlapping_orders_comes_from_one_order_only() -> None:
    # ARRANGE: RED ORDER 1005 TO 1015 NM AT 100 E-, BLUE ORDER 995 TO 1010 NM AT 200 E-, OVERLAPPING FROM 1005 TO 1010
    red = _order(1005.0, 200, value=100.0)
    blue = _order(995.0, 300, value=200.0)
    join = 1007.5

    # ACT
    merged = merge_orders_on_log_grid([red, blue], [join], STEP_KMS)

    # ASSERT: NOTHING BETWEEN THE TWO VALUES, AND ONE HARD STEP
    flux = merged.values["flux"].round(9)
    assert set(np.unique(flux)) == {100.0, 200.0}
    assert merged.hasData.all()
    stepAt = np.flatnonzero(np.diff(flux) != 0)
    assert stepAt.size == 1
    np.testing.assert_allclose(flux[: stepAt[0] + 1], 200.0)
    np.testing.assert_allclose(flux[stepAt[0] + 1 :], 100.0)


def test_the_join_is_snapped_to_the_nearest_bin_edge() -> None:
    # ARRANGE
    red = _order(1005.0, 200, value=100.0)
    blue = _order(995.0, 300, value=200.0)
    join = 1007.5

    # ACT
    merged = merge_orders_on_log_grid([red, blue], [join], STEP_KMS)

    # ASSERT: THE FIRST RED BIN STARTS AT THE SNAPPED JOIN, WITHIN HALF A BIN OF THE REQUESTED ONE
    snapped = merged.joinWavelengths[0]
    binWidth = join * STEP_KMS / SPEED_OF_LIGHT_KMS
    assert abs(snapped - join) <= binWidth / 2 * (1 + 1e-9)
    firstRed = merged.wavelength[merged.values["flux"].round(9) == 100.0][0]
    edgeBelowFirstRed = firstRed * np.exp(-0.5 * STEP_KMS / SPEED_OF_LIGHT_KMS)
    assert edgeBelowFirstRed == pytest.approx(snapped, rel=1e-12)


def test_a_gap_between_orders_gets_no_flux_and_no_variance() -> None:
    # ARRANGE: THE BLUE ORDER ENDS AT 1010 NM AND THE RED ORDER STARTS AT 1020 NM
    red = _order(1020.0, 200, value=100.0)
    blue = _order(1000.0, 200, value=200.0)
    join = 1015.0

    # ACT
    merged = merge_orders_on_log_grid([red, blue], [join], STEP_KMS)

    # ASSERT: ONLY THE BINS OF THE GAP HAVE NO DATA, AND THEY HOLD 0, NOT A LINE BETWEEN THE ORDERS
    gap = ~merged.hasData
    assert gap.any()
    assert (merged.values["flux"][gap] == 0.0).all()
    assert (merged.variance[gap] == 0.0).all()
    assert merged.wavelength[gap].min() > 1009.9
    assert merged.wavelength[gap].max() < 1020.1
    inGap = (merged.wavelength > 1010.1) & (merged.wavelength < 1019.9)
    assert (~merged.hasData[inGap]).all()
    np.testing.assert_allclose(merged.values["flux"][merged.hasData & (merged.wavelength < 1010.0)], 200.0)
    np.testing.assert_allclose(merged.values["flux"][merged.hasData & (merged.wavelength > 1020.0)], 100.0)


def test_the_bins_of_a_gap_form_one_unbroken_run_on_the_grid() -> None:
    # ARRANGE
    red = _order(1020.0, 200, value=100.0)
    blue = _order(1000.0, 200, value=200.0)

    # ACT
    merged = merge_orders_on_log_grid([red, blue], [1015.0], STEP_KMS)

    # ASSERT: THE GRID IS CONTIGUOUS, SO THE GAP BINS ARE PRESENT AND FLAGGED
    steps = merged.wavelength[1:] / merged.wavelength[:-1]
    np.testing.assert_allclose(steps, np.exp(STEP_KMS / SPEED_OF_LIGHT_KMS), rtol=1e-12)


def test_a_bin_that_an_order_covers_only_in_part_at_its_end_is_dropped() -> None:
    # ARRANGE: THREE PIXELS, EACH ONE BIN WIDE, THAT RUN FROM HALF WAY INTO BIN 10 TO HALF WAY INTO BIN 13
    edges = grid_wavelength(np.arange(10, 14), STEP_KMS)
    order = OrderPixels(
        wavelength=(edges[:-1] + edges[1:]) / 2,
        pixelScale=np.diff(edges),
        keep=np.ones(3, dtype=bool),
        values={"flux": np.array([1.0, 1.0, 1.0])},
        variance=np.ones(3),
    )

    # ACT
    merged = merge_orders_on_log_grid([order], [], STEP_KMS)

    # ASSERT: BINS 10 AND 13 ARE ONLY HALF COVERED, SO ONLY BINS 11 AND 12 REMAIN
    assert merged.firstIndex == 11
    assert len(merged.wavelength) == 2
    assert merged.hasData.all()


def test_a_pixel_dropped_before_the_merge_leaves_a_hole_and_no_fill() -> None:
    # ARRANGE: A RUN OF PIXELS WITH A STRETCH OF 40 PIXELS DROPPED, WIDER THAN SEVERAL BINS
    keep = np.ones(200, dtype=bool)
    keep[80:120] = False
    order = _order(1000.0, 200, value=50.0, keep=keep)

    # ACT
    merged = merge_orders_on_log_grid([order], [], STEP_KMS)

    # ASSERT
    hole = ~merged.hasData
    assert hole.sum() > 30
    assert (merged.values["flux"][hole] == 0.0).all()
    np.testing.assert_allclose(merged.values["flux"][merged.hasData], 50.0, rtol=1e-12)
    assert merged.wavelength[hole].min() > 1004.0
    assert merged.wavelength[hole].max() < 1006.1


def test_the_variance_of_a_bin_is_the_sum_of_the_squared_weights() -> None:
    # ARRANGE: PIXELS ABOUT AS WIDE AS THE BINS, SO A BIN TAKES SEVERAL PIXELS
    order = _order(1000.0, 400, value=50.0, variance=9.0)

    # ACT
    merged = merge_orders_on_log_grid([order], [], STEP_KMS)

    # ASSERT: FOR n EQUAL WEIGHTS 1/n THE VARIANCE IS 9/n, NEVER ABOVE THE PIXEL VARIANCE AND ABOVE 9/3
    assert (merged.variance <= 9.0 + 1e-9).all()
    assert (merged.variance > 3.0).all()


def test_the_extra_columns_are_resampled_with_the_same_weights_as_the_flux() -> None:
    # ARRANGE
    order = _order(1000.0, 200, value=50.0, sky=7.0)

    # ACT
    merged = merge_orders_on_log_grid([order], [], STEP_KMS)

    # ASSERT
    np.testing.assert_allclose(merged.values["sky"], 7.0, rtol=1e-12)


def test_the_merge_refuses_orders_that_cover_no_whole_bin() -> None:
    # ARRANGE: ONE PIXEL, NARROWER THAN A BIN
    order = OrderPixels(
        wavelength=np.array([1000.0]),
        pixelScale=np.array([0.01]),
        keep=np.ones(1, dtype=bool),
        values={"flux": np.array([50.0])},
        variance=np.ones(1),
    )

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="No bin of the log-wavelength grid"):
        merge_orders_on_log_grid([order], [], STEP_KMS)


def test_the_merge_needs_one_join_for_each_pair_of_neighbouring_orders() -> None:
    # ARRANGE
    orders = [_order(1005.0, 200, 100.0), _order(995.0, 300, 200.0)]

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="one join"):
        merge_orders_on_log_grid(orders, [], STEP_KMS)


def test_three_orders_are_cut_at_both_joins() -> None:
    # ARRANGE: RED, MIDDLE AND BLUE ORDERS THAT OVERLAP BY 5 NM EACH
    red = _order(1015.0, 400, value=1.0)
    middle = _order(1000.0, 400, value=2.0)
    blue = _order(985.0, 400, value=3.0)

    # ACT
    merged = merge_orders_on_log_grid([red, middle, blue], [1017.5, 1002.5], STEP_KMS)

    # ASSERT: THREE RUNS OF CONSTANT VALUE, BLUE TO RED, WITH NO OTHER VALUE BETWEEN THEM
    flux = merged.values["flux"].round(9)
    runs = flux[np.r_[True, np.diff(flux) != 0]]
    np.testing.assert_allclose(runs, [3.0, 2.0, 1.0])
    assert len(merged.joinWavelengths) == 2


def test_an_order_with_no_usable_pixel_is_skipped_and_the_others_still_merge() -> None:
    # ARRANGE: THE BLUE ORDER HAS NO FINITE PIXEL WIDTH, SO IT GIVES NO BIN
    good = _order(1000.0, 200, value=50.0)
    unusable = OrderPixels(
        wavelength=np.array([990.0, 990.05]),
        pixelScale=np.array([np.nan, np.nan]),
        keep=np.ones(2, dtype=bool),
        values={"flux": np.array([1.0, 1.0])},
        variance=np.ones(2),
    )

    # ACT
    merged = merge_orders_on_log_grid([good, unusable], [995.0], STEP_KMS)

    # ASSERT
    np.testing.assert_allclose(merged.values["flux"], 50.0, rtol=1e-12)
    assert merged.hasData.all()
