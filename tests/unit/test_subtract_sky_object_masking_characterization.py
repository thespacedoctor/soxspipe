"""Characterization of aggressive object slit-range masking in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

PIXELS_PER_ORDER = 4000
BIN_WIDTH = 11.0 / 99
NOISE_COLUMNS = ("residual_windowed_std", "residual_windowed_long_median", "flux_windowed_long_median")


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.recipeSettings = {
        "sky-subtraction": {
            "percentile_rolling_window_size": 15,
            "noise_rolling_window_size": 45,
        }
    }
    return subtractor


def _orders(objectLow: float, objectHigh: float) -> list[pd.DataFrame]:
    """Two orders with 70% of pixels inside the object slit range flagged, plus 3% scattered flags."""
    orders = []
    for offset, seed in enumerate((11, 12)):
        rng = np.random.default_rng(seed)
        slitPosition = rng.uniform(-5.5, 5.5, PIXELS_PER_ORDER)
        insideObject = (slitPosition > objectLow) & (slitPosition < objectHigh)
        objectPixels = insideObject & (rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.7)
        scatteredPixels = rng.uniform(0, 1, PIXELS_PER_ORDER) < 0.03
        index = np.arange(PIXELS_PER_ORDER)
        orders.append(
            pd.DataFrame(
                {
                    "order": 10 + offset,
                    "slit_position": slitPosition,
                    "flagged_object_clipped": objectPixels | scatteredPixels,
                    "flagged_all_clipped": False,
                    "flux_minus_smoothed_residual": 3.0 * np.sin(index / 5.0),
                    "flux_percentile_smoothed": 100.0 + np.cos(index / 50.0),
                }
            )
        )
    return orders


@pytest.fixture
def masked_ranges(monkeypatch: pytest.MonkeyPatch) -> list[tuple[Any, Any]]:
    # THE OBJECT RANGES ARE LOCAL TO THE METHOD; RECORD THEM WHERE THEY ARE APPLIED
    ranges: list[tuple[Any, Any]] = []
    realBetween = pd.Series.between

    def record_between(self: pd.Series, left: Any, right: Any, inclusive: str = "both") -> pd.Series:
        ranges.append((left, right))
        return realBetween(self, left, right, inclusive)

    monkeypatch.setattr(pd.Series, "between", record_between)
    return ranges


def _noise_summary(order: pd.DataFrame) -> list[tuple[Any, int]]:
    return [
        (pytest.approx(float(order[column].sum()), rel=1e-12, abs=0), int(order[column].isna().sum()))
        for column in NOISE_COLUMNS
    ]


def _expected_noise(values: list[tuple[float, int]]) -> list[tuple[Any, int]]:
    return [(pytest.approx(total, rel=1e-12, abs=0), nanCount) for total, nanCount in values]


def test_a_central_object_masks_its_slit_range_in_every_order(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """One range is found from both orders combined and applied twice per order.

    The range runs from the left edge of the first positive bin to the right edge
    of the last one, so it covers the whole 1.0 to 1.6 arcsec object.
    """
    orders = _orders(1.0, 1.6)

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert result is orders
    assert len(masked_ranges) == 4
    assert len(set(masked_ranges)) == 1
    lower, upper = masked_ranges[0]
    assert upper >= 1.6
    assert 1.0 - BIN_WIDTH <= lower <= 1.0
    assert lower == pytest.approx(0.9228259187436887, rel=1e-12, abs=0)
    assert upper == pytest.approx(1.6964064054220023, rel=1e-12, abs=0)
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [392, 338]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [392, 338]
    for order in result:
        insideRange = order["slit_position"].between(lower, upper)
        assert order.loc[insideRange, "flagged_object_clipped"].all()
    assert _noise_summary(result[0]) == _expected_noise(
        [(6385.0131472153225, 392), (137.46062579534546, 405), (359471.17853026406, 405)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6428.201453266465, 338), (-67.38329929430459, 351), (364825.4955009001, 351)]
    )


def test_an_object_at_the_blue_slit_edge_is_masked_from_the_first_examined_bin(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """A positive run starting at the first examined bin is recorded from that bin's left edge."""
    orders = _orders(-5.5, -4.0)
    flaggedBefore = [int(order["flagged_object_clipped"].sum()) for order in orders]

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert flaggedBefore == [501, 482]
    assert len(masked_ranges) == 4
    assert len(set(masked_ranges)) == 1
    lower, upper = masked_ranges[0]
    assert lower is not False
    assert lower == pytest.approx(-4.946292206383818, rel=1e-12, abs=0)
    assert upper == pytest.approx(-3.9505389518741967, rel=1e-12, abs=0)
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [631, 600]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [631, 600]
    assert _noise_summary(result[0]) == _expected_noise(
        [(6242.460744138139, 631), (19.521000344039592, 644), (335582.05770925974, 644)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6226.8607221658385, 600), (231.64411541533178, 613), (338658.9684423819, 613)]
    )


def test_without_an_object_range_aggressive_masking_leaves_object_flags_unpromoted(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """An object in the red edge margin yields no range, and object pixels stay out of `flagged_all_clipped`.

    The non-aggressive path promotes every object flag; the aggressive path only
    does so inside the per-range loop.
    """
    orders = _orders(4.5, 5.5)

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert masked_ranges == []
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [367, 349]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [0, 0]
    assert _noise_summary(result[0]) == _expected_noise(
        [(6411.643045030396, 367), (103.67070347766978, 380), (361976.18861746334, 380)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6405.264848625138, 349), (75.97294583409979, 362), (363737.83604636916, 362)]
    )


def test_the_same_red_edge_object_is_promoted_without_aggressive_masking(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """For contrast with the aggressive path, every object flag becomes an exclusion."""
    orders = _orders(4.5, 5.5)

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=False)

    assert masked_ranges == []
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [367, 349]


def test_an_object_running_into_the_red_margin_is_closed_after_the_loop(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """A positive run still open at the last examined bin is recorded once the loop ends.

    The range stops at 4.94, the right edge of the last examined bin, though the
    object continues to the slit end inside the edge margin.
    """
    orders = _orders(3.8, 5.5)
    flaggedBefore = [int(order["flagged_object_clipped"].sum()) for order in orders]

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    expectedRange = (
        pytest.approx(3.722386682854049, rel=1e-12, abs=0),
        pytest.approx(4.942886733348475, rel=1e-12, abs=0),
    )
    assert masked_ranges == [expectedRange] * 4
    assert flaggedBefore == [527, 532]
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [654, 672]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [654, 672]
    assert _noise_summary(result[0]) == _expected_noise(
        [(6212.313468960861, 654), (127.97333372969166, 667), (333259.15392645646, 667)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6143.968166835047, 672), (43.69186079400444, 685), (331434.72029450943, 685)]
    )


def _ranges(log: Any, counts: list[float], edgeMargin: int = 0) -> list[list[float]]:
    """Run detection on unit-width bins, so bin `i` spans `i` to `i + 1`."""
    binEdges = np.arange(len(counts) + 1, dtype=float)
    return _subtractor(log)._object_slit_ranges(np.array(counts), binEdges, edgeMargin)


def test_a_run_ends_at_the_right_edge_of_its_last_positive_bin(log: Any) -> None:
    ranges = _ranges(log, [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0])

    assert ranges == [[1.0, 6.0]]


def test_a_run_starting_at_the_first_examined_bin_starts_at_that_bins_left_edge(log: Any) -> None:
    ranges = _ranges(log, [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.0])

    assert ranges == [[0.0, 5.0]]


def test_a_run_of_four_positive_bins_is_not_recorded_but_five_are(log: Any) -> None:
    fourBins = _ranges(log, [0.0, 0.2, 0.2, 0.2, 0.2, 0.0])
    fiveBins = _ranges(log, [0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0])

    assert fourBins == []
    assert fiveBins == [[1.0, 6.0]]


def test_a_run_whose_peak_does_not_exceed_the_threshold_is_not_recorded_inside_the_loop(log: Any) -> None:
    ranges = _ranges(log, [0.0, 0.05, 0.01, 0.01, 0.01, 0.01, 0.0, 0.0])

    assert ranges == []


def test_a_final_run_whose_peak_does_not_exceed_the_threshold_is_not_recorded(log: Any) -> None:
    ranges = _ranges(log, [0.0, 0.0, 0.05, 0.01, 0.01, 0.01, 0.01])

    assert ranges == []


def test_a_final_run_with_a_peak_is_closed_at_the_right_edge_of_the_last_examined_bin(log: Any) -> None:
    ranges = _ranges(log, [0.0, 0.0, 0.01, 0.01, 0.2, 0.01, 0.01])

    assert ranges == [[2.0, 7.0]]


def test_the_peak_of_one_run_does_not_carry_over_to_the_next_run(log: Any) -> None:
    ranges = _ranges(log, [0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.01, 0.01, 0.01, 0.01, 0.01, 0.0])

    assert ranges == [[0.0, 5.0]]


def test_bins_inside_the_edge_margin_are_not_examined(log: Any) -> None:
    ranges = _ranges(log, [0.2, 0.2, 0.0, 0.0, 0.2, 0.2, 0.2, 0.2, 0.2, 0.0, 0.2, 0.2], edgeMargin=2)

    assert ranges == [[4.0, 9.0]]
