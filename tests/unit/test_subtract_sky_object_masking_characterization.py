"""Characterization of aggressive object slit-range masking in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

PIXELS_PER_ORDER = 4000
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

    The range ends at the left edge of the last positive bin, 1.586, short of
    the object's 1.6 arcsec edge.
    """
    # SUSPICIOUS (WRONG SCIENCE): RANGE TRIMS THE LAST OBJECT BIN, FILED AS DY-596
    orders = _orders(1.0, 1.6)

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert result is orders
    expectedRange = (
        pytest.approx(0.8123144206467865, rel=1e-12, abs=0),
        pytest.approx(1.5858949073251, rel=1e-12, abs=0),
    )
    assert masked_ranges == [expectedRange] * 4
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [406, 336]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [406, 336]
    for order in result:
        insideRange = order["slit_position"].between(0.8123144206467865, 1.5858949073251)
        assert order.loc[insideRange, "flagged_object_clipped"].all()
    assert _noise_summary(result[0]) == _expected_noise(
        [(6371.956674832487, 406), (88.6669295593187, 419), (358063.37340630684, 419)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6422.569136110474, 336), (-115.92489730639494, 349), (365026.1278999605, 349)]
    )


def test_an_object_at_the_blue_slit_edge_is_recorded_from_false_and_masks_nothing(
    log: Any, masked_ranges: list[tuple[Any, Any]]
) -> None:
    """A positive run starting at the first examined bin keeps `lower = False`, an empty range."""
    # SUSPICIOUS (WRONG SCIENCE): EDGE OBJECT RANGE STARTS AT False AND SELECTS NOTHING, FILED AS DY-596
    orders = _orders(-5.5, -4.0)
    flaggedBefore = [int(order["flagged_object_clipped"].sum()) for order in orders]

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=True)

    assert masked_ranges == [(False, pytest.approx(-4.061178202375266, rel=1e-12, abs=0))] * 4
    assert masked_ranges[0][0] is False
    assert flaggedBefore == [501, 482]
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [501, 482]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [501, 482]
    assert _noise_summary(result[0]) == _expected_noise(
        [(6322.006763896867, 501), (17.233931234827608, 514), (348574.90063246945, 514)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6316.503674671294, 482), (259.5343562330326, 495), (350454.6068999977, 495)]
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
