"""Contracts for the both-lit overlap mask built by ``slit_illumination_mask.overlap_exclude_mask`` (DY-1280)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_illumination import (
    FLAT_EDGE_COLUMNS,
    OFFSET_COLUMNS,
    STATUS_BORROWED,
    STATUS_FAILED,
    STATUS_FITTED,
)
from soxspipe.commonutils.slit_illumination_mask import overlap_exclude_mask

pytestmark = pytest.mark.unit

# EACH ORDER IS ONE ROW OF THE TEST IMAGES, WITH SLIT POSITIONS FROM -6 TO +6 ARCSEC IN 0.1 ARCSEC STEPS
SLIT_STEPS = np.round(np.arange(-60, 61) * 0.1, 1)
ORDERS = (11, 12)
FLAT_LOWER = -5.0
FLAT_UPPER = 5.0


def _images(orders=ORDERS):
    """One image row per order, plus a last row of inter-order pixels with a NaN order and slit."""
    slit = np.vstack([SLIT_STEPS for _ in orders] + [np.full(SLIT_STEPS.size, np.nan)])
    orderMap = np.vstack([np.full(SLIT_STEPS.size, float(order)) for order in orders])
    orderMap = np.vstack([orderMap, np.full(SLIT_STEPS.size, np.nan)])
    return slit, orderMap


def _flat_edges(edges=None):
    edges = edges or dict.fromkeys(ORDERS, (FLAT_LOWER, FLAT_UPPER))
    rows = [
        {"order": order, "lower": lower, "upper": upper, "width": 0.15, "reason": ""}
        for order, (lower, upper) in edges.items()
    ]
    return pd.DataFrame(rows, columns=FLAT_EDGE_COLUMNS)


def _offsets(values, statuses=None):
    statuses = statuses or {}
    rows = [
        {
            "order": order,
            "offset": offset,
            "offsetErr": 0.02,
            "scienceWidth": 0.15,
            "nBins": 100,
            "status": statuses.get(order, STATUS_FITTED),
            "reason": "" if statuses.get(order, STATUS_FITTED) != STATUS_FAILED else "fit did not converge",
            "borrowedFrom": np.nan,
        }
        for order, offset in values.items()
    ]
    return pd.DataFrame(rows, columns=OFFSET_COLUMNS)


def _excluded_slit_positions(mask, row):
    return set(SLIT_STEPS[mask[row]].tolist())


def _expected_exclusions(keepLower, keepUpper):
    return {s for s in SLIT_STEPS.tolist() if s < keepLower - 1e-9 or s > keepUpper + 1e-9}


def test_a_positive_offset_trims_the_lower_end_by_the_offset_and_both_ends_by_the_margin(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 1.0, 12: 1.0}), 0.3)

    # BOTH-LIT RANGE IS [max(-5, -4), min(5, 6)] = [-4, 5], TRIMMED BY 0.3 TO [-3.7, 4.7]
    assert _excluded_slit_positions(mask, 0) == _expected_exclusions(-3.7, 4.7)
    assert _excluded_slit_positions(mask, 1) == _expected_exclusions(-3.7, 4.7)


def test_a_negative_offset_trims_the_upper_end_by_the_offset(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: -1.0, 12: -1.0}), 0.3)

    # BOTH-LIT RANGE IS [max(-5, -6), min(5, 4)] = [-5, 4], TRIMMED BY 0.3 TO [-4.7, 3.7]
    assert _excluded_slit_positions(mask, 0) == _expected_exclusions(-4.7, 3.7)


def test_each_order_is_masked_with_its_own_offset(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: -1.0, 12: 1.0}), 0.3)

    assert _excluded_slit_positions(mask, 0) == _expected_exclusions(-4.7, 3.7)
    assert _excluded_slit_positions(mask, 1) == _expected_exclusions(-3.7, 4.7)


def test_a_zero_margin_keeps_the_whole_both_lit_range(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 0.5, 12: 0.5}), 0.0)

    assert _excluded_slit_positions(mask, 0) == _expected_exclusions(-4.5, 5.0)


def test_each_order_uses_its_own_flat_edges(log):
    slit, orderMap = _images()
    flatEdges = _flat_edges({11: (-5.0, 5.0), 12: (-4.6, 5.2)})

    mask = overlap_exclude_mask(log, slit, orderMap, flatEdges, _offsets({11: 0.0, 12: 0.0}), 0.3)

    assert _excluded_slit_positions(mask, 0) == _expected_exclusions(-4.7, 4.7)
    assert _excluded_slit_positions(mask, 1) == _expected_exclusions(-4.3, 4.9)


def test_inter_order_pixels_are_never_excluded(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 1.0, 12: 1.0}), 0.3)

    assert not mask[-1].any()


def test_a_pixel_on_a_fractional_order_boundary_is_never_excluded(log):
    slit, orderMap = _images()
    # A BINNED MAP AVERAGES ORDER NUMBERS ACROSS A BOUNDARY, SO SUCH A PIXEL BELONGS TO NO ORDER
    orderMap[0, 0] = 11.5

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 1.0, 12: 1.0}), 0.3)

    assert not mask[0, 0]
    assert mask[0, 1]


def test_a_failed_order_is_not_masked_but_a_borrowed_one_is(log):
    slit, orderMap = _images()
    offsets = _offsets({11: np.nan, 12: 1.0}, statuses={11: STATUS_FAILED, 12: STATUS_BORROWED})

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), offsets, 0.3)

    assert not mask[0].any()
    assert _excluded_slit_positions(mask, 1) == _expected_exclusions(-3.7, 4.7)


def test_an_order_whose_status_is_not_fitted_or_borrowed_is_not_masked(log):
    slit, orderMap = _images()
    offsets = _offsets({11: 1.0, 12: 1.0}, statuses={11: "unknown"})

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), offsets, 0.3)

    assert not mask[0].any()
    assert mask[1].any()


def test_an_order_with_no_offset_row_or_no_flat_edges_is_not_masked(log):
    slit, orderMap = _images()
    flatEdges = _flat_edges({11: (np.nan, np.nan), 12: (FLAT_LOWER, FLAT_UPPER)})

    mask = overlap_exclude_mask(log, slit, orderMap, flatEdges, _offsets({11: 1.0}), 0.3)

    assert not mask.any()


def test_an_order_whose_both_lit_range_is_empty_is_not_masked_and_is_logged(log):
    slit, orderMap = _images()
    flatEdges = _flat_edges({11: (-0.2, 0.2), 12: (FLAT_LOWER, FLAT_UPPER)})

    mask = overlap_exclude_mask(log, slit, orderMap, flatEdges, _offsets({11: 0.0, 12: 0.0}), 0.3)

    assert not mask[0].any()
    assert mask[1].any()
    assert any(level == "warning" and "order 11" in message for level, message in log.messages)


def test_the_mask_is_boolean_and_the_shape_of_the_slit_image(log):
    slit, orderMap = _images()

    mask = overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 1.0, 12: 1.0}), 0.3)

    assert mask.dtype == bool
    assert mask.shape == slit.shape


def test_images_of_different_shapes_raise(log):
    slit, orderMap = _images()

    with pytest.raises(ValueError, match="shape"):
        overlap_exclude_mask(log, slit, orderMap[:-1], _flat_edges(), _offsets({11: 1.0}), 0.3)


@pytest.mark.parametrize("margin", [-0.1, np.nan, np.inf])
def test_a_margin_that_is_negative_or_not_finite_raises(log, margin):
    slit, orderMap = _images()

    with pytest.raises(ValueError, match="margin"):
        overlap_exclude_mask(log, slit, orderMap, _flat_edges(), _offsets({11: 1.0}), margin)


def test_does_not_modify_its_inputs(log):
    slit, orderMap = _images()
    flatEdges = _flat_edges()
    offsets = _offsets({11: 1.0, 12: -1.0})
    copies = (slit.copy(), orderMap.copy(), flatEdges.copy(), offsets.copy())

    overlap_exclude_mask(log, slit, orderMap, flatEdges, offsets, 0.3)

    np.testing.assert_array_equal(slit, copies[0])
    np.testing.assert_array_equal(orderMap, copies[1])
    pd.testing.assert_frame_equal(flatEdges, copies[2])
    pd.testing.assert_frame_equal(offsets, copies[3])
