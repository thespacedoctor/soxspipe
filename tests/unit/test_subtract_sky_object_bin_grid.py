"""The object slit-range bins sit on a fixed grid, so the flagged slit extent cannot move them (DY-1982)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import OBJECT_BIN_WIDTH, subtract_sky

pytestmark = pytest.mark.unit

PIXELS_PER_ORDER = 4000
OBJECT_LOW = 0.9
OBJECT_HIGH = 1.7


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.ron = 3.3
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "percentile_rolling_window_size": 15,
            "noise_rolling_window_size": 45,
        }
    }
    return subtractor


def _order(seed: int, extraSlitPositions: tuple[float, ...] = (), isObject: bool = True) -> pd.DataFrame:
    """One order with scattered flags over -4.9 to 4.9 arcsec, an object at 0.9 to 1.7, and extra flagged pixels."""
    rng = np.random.default_rng(seed)
    slitPosition = rng.uniform(-4.9, 4.9, PIXELS_PER_ORDER)
    scattered = rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.03
    inObject = (slitPosition > OBJECT_LOW) & (slitPosition < OBJECT_HIGH) & (rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.7)
    flagged = scattered | (inObject if isObject else False)
    frame = pd.DataFrame(
        {
            "slit_position": slitPosition,
            "flagged_object_clipped": flagged,
            "flagged_all_clipped": False,
        }
    )
    extra = pd.DataFrame(
        {
            "slit_position": list(extraSlitPositions),
            "flagged_object_clipped": True,
            "flagged_all_clipped": False,
        }
    )
    frame = pd.concat([frame, extra], ignore_index=True)
    index = np.arange(len(frame))
    frame["flux_minus_smoothed_residual"] = 3.0 * np.sin(index / 5.0)
    frame["flux_percentile_smoothed"] = 100.0 + np.cos(index / 50.0)
    return frame


def _detected(log: Any, orders: list[pd.DataFrame]) -> tuple[list[list[float]], np.ndarray]:
    """Run the aggressive clip and return the object ranges it found and the bin edges it used."""
    subtractor = _subtractor(log)
    captured: dict[str, Any] = {}
    realRanges = subtractor._object_slit_ranges

    def record(counts: np.ndarray, binEdges: np.ndarray, edgeMargin: int) -> list[list[float]]:
        captured["edges"] = np.asarray(binEdges)
        captured["ranges"] = realRanges(counts, binEdges, edgeMargin)
        return captured["ranges"]

    subtractor._object_slit_ranges = record
    subtractor.clip_object_slit_positions(orders, aggressive=True)
    return captured["ranges"], captured["edges"]


@pytest.mark.parametrize("extraPosition", [-5.4, -5.04, -4.96])
def test_a_flagged_pixel_at_the_far_slit_edge_of_one_order_does_not_move_the_object_range(
    log: Any, extraPosition: float
) -> None:
    before, _ = _detected(log, [_order(11), _order(12)])

    after, _ = _detected(log, [_order(11, (extraPosition,)), _order(12)])

    assert len(before) == 1
    assert after == before


def test_a_flagged_pixel_at_the_red_slit_edge_of_one_order_does_not_move_the_object_range(log: Any) -> None:
    before, _ = _detected(log, [_order(11), _order(12)])
    assert len(before) == 1

    after, _ = _detected(log, [_order(11), _order(12, (5.37,))])

    assert after == before


@pytest.mark.parametrize("lowestFlag", [-5.04, -4.96, -5.23, -5.5])
def test_bin_edges_are_whole_multiples_of_the_bin_width_whatever_the_lowest_flag(log: Any, lowestFlag: float) -> None:
    _, edges = _detected(log, [_order(11, (lowestFlag,)), _order(12)])

    multiples = edges / OBJECT_BIN_WIDTH

    assert np.allclose(multiples, np.round(multiples), atol=1e-9, rtol=0)
    assert np.allclose(np.diff(edges), OBJECT_BIN_WIDTH)
    assert edges[0] <= lowestFlag < edges[0] + OBJECT_BIN_WIDTH + 1e-9


@pytest.mark.parametrize("lowest, highest", [(57.9, 58.4), (1.7, 2.2), (-1.2, 4.3), (-0.9, 8.1)])
def test_the_bin_edges_enclose_the_extent_whatever_the_float_rounding(log: Any, lowest: float, highest: float) -> None:
    edges = _subtractor(log)._object_bin_edges(lowest, highest)

    assert edges[0] <= lowest
    assert edges[-1] >= highest


def test_an_empty_flagged_set_gives_no_object_range(log: Any) -> None:
    orders = [_order(11, isObject=False), _order(12, isObject=False)]
    for order in orders:
        order["flagged_object_clipped"] = False

    subtractor = _subtractor(log)
    result = subtractor.clip_object_slit_positions(orders, aggressive=True)

    assert result is orders
    assert all(not order["flagged_all_clipped"].any() for order in result)


def test_a_flagged_set_with_no_finite_slit_position_gives_no_object_range(log: Any) -> None:
    orders = [_order(11, isObject=False)]
    orders[0]["slit_position"] = np.nan

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert not result[0]["flagged_all_clipped"].any()


def test_an_object_longer_than_the_minimum_run_is_detected_with_its_edges_on_the_grid(log: Any) -> None:
    ranges, _ = _detected(log, [_order(11), _order(12)])

    assert len(ranges) == 1
    lower, upper = ranges[0]
    assert OBJECT_LOW - OBJECT_BIN_WIDTH <= lower <= OBJECT_LOW
    assert OBJECT_HIGH <= upper <= OBJECT_HIGH + OBJECT_BIN_WIDTH
    assert lower / OBJECT_BIN_WIDTH == pytest.approx(round(lower / OBJECT_BIN_WIDTH), abs=1e-9)


def test_a_flat_distribution_of_flagged_pixels_gives_no_object_range(log: Any) -> None:
    ranges, _ = _detected(log, [_order(11, isObject=False), _order(12, isObject=False)])

    assert ranges == []


def test_a_uniform_distribution_has_one_normalised_count_per_bin(log: Any) -> None:
    positions = (np.arange(1000) + 0.5) * OBJECT_BIN_WIDTH - 50.0

    edges, counts = _subtractor(log)._normalised_object_counts(pd.Series(positions))

    assert len(counts) == len(edges) - 1 == 1000
    assert counts.to_numpy() == pytest.approx(1.0)


@pytest.mark.parametrize("positions", [[57.9], [1.7, 1.7], [-0.9, 4.3], [57.9, 58.4, 1.7], [0.05, 8.1]])
def test_every_flagged_position_is_counted_in_the_normalised_counts(log: Any, positions: list[float]) -> None:
    edges, counts = _subtractor(log)._normalised_object_counts(pd.Series(positions))

    assert counts.sum() == pytest.approx(len(edges) - 1)
    assert not counts.isna().any()


def test_slit_positions_that_are_not_finite_are_ignored_in_the_normalised_counts(log: Any) -> None:
    subtractor = _subtractor(log)

    edgesA, countsA = subtractor._normalised_object_counts(pd.Series([0.05, 0.15, np.nan, np.inf]))
    edgesB, countsB = subtractor._normalised_object_counts(pd.Series([0.05, 0.15]))

    assert np.array_equal(edgesA, edgesB)
    assert countsA.to_numpy() == pytest.approx(countsB.to_numpy())
