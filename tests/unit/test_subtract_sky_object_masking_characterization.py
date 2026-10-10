"""Characterization of object-flag promotion in subtract_sky without aggressive masking (DY-40).

The aggressive path builds each order's mask from its spatial profile; it is tested in
``test_subtract_sky_object_profile_mask.py`` (DY-1279).
"""

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


def _noise_summary(order: pd.DataFrame) -> list[tuple[Any, int]]:
    return [
        (pytest.approx(float(order[column].sum()), rel=1e-12, abs=0), int(order[column].isna().sum()))
        for column in NOISE_COLUMNS
    ]


def _expected_noise(values: list[tuple[float, int]]) -> list[tuple[Any, int]]:
    return [(pytest.approx(total, rel=1e-12, abs=0), nanCount) for total, nanCount in values]


def test_every_object_flag_is_promoted_without_aggressive_masking(log: Any) -> None:
    """Without aggressive masking, every object flag becomes an exclusion and no slit range is added."""
    orders = _orders(4.5, 5.5)
    flaggedBefore = [int(order["flagged_object_clipped"].sum()) for order in orders]

    result = _subtractor(log).clip_object_slit_positions(orders, aggressive=False)

    assert flaggedBefore == [367, 349]
    assert [int(order["flagged_object_clipped"].sum()) for order in result] == [367, 349]
    assert [int(order["flagged_all_clipped"].sum()) for order in result] == [367, 349]
    assert _noise_summary(result[0]) == _expected_noise(
        [(6411.643045030396, 367), (103.67070347766978, 380), (361976.18861746334, 380)]
    )
    assert _noise_summary(result[1]) == _expected_noise(
        [(6405.264848625138, 349), (75.97294583409979, 362), (363737.83604636916, 362)]
    )
