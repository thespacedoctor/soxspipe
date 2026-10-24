"""Faint-order detection, anchoring and the binned trace fit (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_trace import (
    VIS_BAND_GROUPS,
    anchor_neighbours,
    faint_trace_settings,
    find_faint_orders,
)
from soxspipe.commonutils.toolkit import SOXS_VIS_ORDERS_BLUE_TO_RED

pytestmark = pytest.mark.unit

SETTINGS = faint_trace_settings.from_block({"faint-order-trace": True})


def _order_rows(
    order: float,
    count: int = 100,
    detectedFraction: float = 1.0,
    scatterPx: float = 0.0,
    curvaturePx: float = 0.0,
    seed: int = 1941,
) -> pd.DataFrame:
    """Return the sampled rows of one order: `detectedFraction` of them found, with Gaussian centre-shift noise."""
    rng = np.random.default_rng(seed + int(order))
    wavelength = np.linspace(400.0, 500.0, count)
    unit = (wavelength - 450.0) / 50.0
    shift = (
        2.0 + curvaturePx * unit**2 + rng.normal(0.0, scatterPx, count) if scatterPx else 2.0 + curvaturePx * unit**2
    )
    isDetected = np.arange(count) < round(detectedFraction * count)
    fitX = np.full(count, 100.0 * order)
    return pd.DataFrame(
        {
            "order": order,
            "wavelength": wavelength,
            "fit_x": fitX,
            "fit_y": np.linspace(10.0, 10.0 + count - 1, count),
            "cont_x": np.where(isDetected, fitX - shift, np.nan),
            "cont_y": np.where(isDetected, np.linspace(10.0, 10.0 + count - 1, count), np.nan),
            "centre_shift": np.where(isDetected, shift, np.nan),
            "pre-clipped": ~isDetected,
        }
    )


def _table(*orderRows: pd.DataFrame) -> pd.DataFrame:
    return pd.concat(orderRows, ignore_index=True)


def test_find_faint_orders_flags_an_order_below_the_detection_fraction() -> None:
    # ARRANGE: ORDER 2 IS DETECTED IN 40% OF ITS SAMPLES, ORDER 3 IN ALL
    table = _table(_order_rows(2.0, detectedFraction=0.4), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert list(faint) == [2.0]
    assert faint[2.0].detectionFraction == pytest.approx(0.4)


def test_find_faint_orders_flags_an_order_whose_centre_shift_scatters_more_than_the_limit() -> None:
    # ARRANGE: ORDER 2 IS DETECTED EVERYWHERE, BUT ITS CENTRES SCATTER BY 3 PX
    table = _table(_order_rows(2.0, scatterPx=3.0), _order_rows(3.0, scatterPx=0.2))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert list(faint) == [2.0]
    assert faint[2.0].scatterPx == pytest.approx(3.0, rel=0.3)


def test_find_faint_orders_leaves_a_curved_but_clean_trace_alone() -> None:
    # ARRANGE: A 4 PX PARABOLA IS THE FIT'S OWN SHAPE, NOT SCATTER
    table = _table(_order_rows(3.0, curvaturePx=4.0, scatterPx=0.2))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert faint == {}


def test_find_faint_orders_flags_an_order_with_no_detection() -> None:
    # ARRANGE
    table = _table(_order_rows(2.0, detectedFraction=0.0), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert faint[2.0].detectionFraction == 0.0
    assert np.isnan(faint[2.0].scatterPx)


def test_find_faint_orders_flags_an_order_with_too_few_accepted_samples_for_a_fit() -> None:
    # ARRANGE: THREE ACCEPTED SAMPLES IN 6 ARE DETECTED, AND A DEGREE-2 FIT NEEDS 4
    table = _table(_order_rows(2.0, count=6, detectedFraction=0.5), _order_rows(3.0))

    # ACT
    faint = find_faint_orders(table, "x", SETTINGS)

    # ASSERT
    assert 2.0 in faint


def test_find_faint_orders_ignores_pre_clipped_samples_when_measuring_scatter() -> None:
    # ARRANGE: A CLEAN ORDER WHOSE FIVE OUTLIERS WERE PRE-CLIPPED BY THE SAMPLER
    rows = _order_rows(3.0, scatterPx=0.2)
    rows.loc[:4, "centre_shift"] = [30.0, -30.0, 25.0, -25.0, 20.0]
    rows.loc[:4, "pre-clipped"] = True

    # ACT
    faint = find_faint_orders(rows, "x", SETTINGS)

    # ASSERT
    assert faint == {}


def test_find_faint_orders_uses_the_limits_in_the_settings() -> None:
    # ARRANGE: 70% DETECTED IS FAINT UNDER A 0.8 LIMIT
    table = _table(_order_rows(2.0, detectedFraction=0.7), _order_rows(3.0))
    strict = faint_trace_settings.from_block({"faint-order-trace": True, "faint-order-min-detection-fraction": 0.8})

    # ACT
    faintDefault = find_faint_orders(table, "x", SETTINGS)
    faintStrict = find_faint_orders(table, "x", strict)

    # ASSERT
    assert faintDefault == {}
    assert list(faintStrict) == [2.0]


def _orders_table(*orders: float) -> pd.DataFrame:
    return pd.DataFrame({"order": [float(order) for order in orders for _ in range(3)]})


def test_band_groups_are_the_adjacent_pairs_of_the_blue_to_red_order_list() -> None:
    # ARRANGE
    blueToRed = SOXS_VIS_ORDERS_BLUE_TO_RED

    # ACT
    pairs = [tuple(blueToRed[start : start + 2]) for start in (0, 2)]

    # ASSERT
    assert tuple(pairs) == VIS_BAND_GROUPS


@pytest.mark.parametrize(
    ("faintOrder", "expected"),
    [(2.0, 3.0), (4.0, 1.0), (3.0, 2.0), (1.0, 4.0)],
    ids=["u-from-g", "i-from-r", "g-from-u", "r-from-i"],
)
def test_anchor_neighbours_takes_the_other_order_of_the_band_group(faintOrder: float, expected: float) -> None:
    # ARRANGE
    table = _orders_table(1, 2, 3, 4)

    # ACT
    anchors = anchor_neighbours(table, {faintOrder: None})

    # ASSERT
    assert anchors == {faintOrder: expected}


def test_anchor_neighbours_gives_no_anchor_when_both_orders_of_a_group_are_faint() -> None:
    # ARRANGE: U AND G ARE BOTH FAINT, AND R IS BRIGHT BUT ACROSS THE BAND-GROUP BOUNDARY FROM G
    table = _orders_table(1, 2, 3, 4)

    # ACT
    anchors = anchor_neighbours(table, {2.0: None, 3.0: None, 4.0: None})

    # ASSERT: ONLY I GETS AN ANCHOR
    assert anchors == {4.0: 1.0}


def test_anchor_neighbours_skips_an_order_whose_partner_is_not_in_the_table() -> None:
    # ARRANGE
    table = _orders_table(2, 1, 4)

    # ACT
    anchors = anchor_neighbours(table, {2.0: None})

    # ASSERT
    assert anchors == {}


def test_anchor_neighbours_ignores_orders_outside_the_vis_band_groups() -> None:
    # ARRANGE
    table = _orders_table(11, 12)

    # ACT
    anchors = anchor_neighbours(table, {11.0: None})

    # ASSERT
    assert anchors == {}
