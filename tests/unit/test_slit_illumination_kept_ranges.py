"""Contracts for ``slit_illumination_mask.flat_fit_ranges_from_result``, the slit range to fit the flat over.

The range is the both-lit range of the order (DY-2147). An order without a usable offset takes the median offset of
the frame's usable orders, and a frame with no usable offset takes an offset of 0 and an extra inset (DY-1977).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_illumination import (
    FLAT_EDGE_COLUMNS,
    OFFSET_COLUMNS,
    PROFILE_COLUMNS,
    STATUS_FAILED,
    STATUS_FITTED,
    slit_illumination_result,
)
from soxspipe.commonutils.slit_illumination_mask import (
    FLAT_FIT_CLEARANCE_ARCSEC,
    flat_fit_offsets,
    flat_fit_ranges_from_result,
    kept_slit_range,
)

pytestmark = pytest.mark.unit

MARGIN_ARCSEC = 0.3
MAX_OFFSET_ERR_ARCSEC = 0.25
# THE INSET THAT THE OVERLAP MARGIN LEAVES SHORT OF THE CLEARANCE; THE SHIPPED MARGIN LEAVES NONE
CLEARANCE_INSET_ARCSEC = FLAT_FIT_CLEARANCE_ARCSEC - MARGIN_ARCSEC


def _result(offsets, edges, isApplied=True, margin=MARGIN_ARCSEC, maxOffsetErr=MAX_OFFSET_ERR_ARCSEC):
    """A result whose ``offsets`` is ``{order: (offset, status)}`` and ``edges`` is ``{order: (lower, upper)}``."""
    offsetRows = [
        {
            "order": order,
            "offset": offset,
            "offsetErr": 0.02,
            "scienceWidth": 0.15,
            "nBins": 100,
            "status": status,
            "reason": "" if status == STATUS_FITTED else "fit did not converge",
            "borrowedFrom": np.nan,
        }
        for order, (offset, status) in offsets.items()
    ]
    edgeRows = [
        {"order": order, "lower": lower, "upper": upper, "width": 0.15, "reason": ""}
        for order, (lower, upper) in edges.items()
    ]
    return slit_illumination_result(
        offsets=pd.DataFrame(offsetRows, columns=OFFSET_COLUMNS),
        flatEdges=pd.DataFrame(edgeRows, columns=FLAT_EDGE_COLUMNS),
        profiles=pd.DataFrame(columns=PROFILE_COLUMNS),
        isApplied=isApplied,
        reason="" if isApplied else "no order has a usable offset",
        overlapMarginArcsec=margin,
        maxOffsetErrArcsec=maxOffsetErr,
    )


def _fit_range(lower, upper, offset, extraInset=0.0):
    """The flat fit range of one order by the rule: the both-lit range, the clearance inset, then the extra inset."""
    keepLower, keepUpper = kept_slit_range(lower, upper, offset, MARGIN_ARCSEC)
    keepLower, keepUpper = keepLower + CLEARANCE_INSET_ARCSEC, keepUpper - CLEARANCE_INSET_ARCSEC
    return (keepLower + extraInset, keepUpper - extraInset) if extraInset else (keepLower, keepUpper)


def test_an_order_with_a_usable_offset_gets_the_fit_range_of_its_own_offset():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (-0.5, STATUS_FITTED)},
        {11: (-5.0, 5.0), 12: (-5.5, 5.5)},
    )

    ranges = flat_fit_ranges_from_result(result)

    assert ranges == {11: _fit_range(-5.0, 5.0, 1.0), 12: _fit_range(-5.5, 5.5, -0.5)}


def test_an_order_with_no_usable_offset_takes_the_median_offset_of_the_usable_orders():
    # ARRANGE: OFFSETS 0.1, 0.4 AND 0.9 HAVE A MEDIAN OF 0.4 AND A MEAN OF 0.467
    result = _result(
        {11: (0.1, STATUS_FITTED), 12: (0.4, STATUS_FITTED), 13: (0.9, STATUS_FITTED), 14: (np.nan, STATUS_FAILED)},
        {11: (-5.0, 5.0), 12: (-5.0, 5.0), 13: (-5.0, 5.0), 14: (-5.5, 5.5), 15: (-5.0, 5.0)},
    )

    ranges = flat_fit_ranges_from_result(result)

    assert ranges[14] == _fit_range(-5.5, 5.5, 0.4)
    assert ranges[15] == _fit_range(-5.0, 5.0, 0.4)


def test_a_usable_order_keeps_its_range_when_another_order_is_unusable():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (np.nan, STATUS_FAILED)},
        {11: (-5.0, 5.0), 12: (-5.5, 5.5)},
    )

    assert flat_fit_ranges_from_result(result)[11] == _fit_range(-5.0, 5.0, 1.0)


def test_a_frame_with_no_usable_offset_gets_an_offset_of_zero_and_an_extra_inset_of_the_order_scatter():
    result = _result(
        {11: (np.nan, STATUS_FAILED), 12: (np.nan, STATUS_FAILED)},
        {11: (-5.0, 5.0), 12: (-5.5, 5.5)},
        isApplied=False,
    )

    ranges = flat_fit_ranges_from_result(result)

    assert ranges == {
        11: _fit_range(-5.0, 5.0, 0.0, MAX_OFFSET_ERR_ARCSEC),
        12: _fit_range(-5.5, 5.5, 0.0, MAX_OFFSET_ERR_ARCSEC),
    }


def test_a_frame_whose_fitted_offsets_are_all_not_finite_takes_the_offset_zero_fallback():
    result = _result({11: (np.nan, STATUS_FITTED)}, {11: (-5.0, 5.0)})

    assert flat_fit_ranges_from_result(result) == {11: _fit_range(-5.0, 5.0, 0.0, MAX_OFFSET_ERR_ARCSEC)}


def test_a_result_with_no_flat_edges_gives_no_ranges():
    result = _result({11: (np.nan, STATUS_FAILED)}, {}, isApplied=False)

    assert flat_fit_ranges_from_result(result) == {}


def test_leaves_out_an_order_whose_range_is_empty_or_whose_flat_edges_are_not_finite():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (0.0, STATUS_FITTED)},
        {11: (-0.1, 0.1), 12: (np.nan, 5.0)},
    )

    assert flat_fit_ranges_from_result(result) == {}


def test_the_extra_inset_leaves_out_an_order_it_empties():
    # ARRANGE: THE RANGE IS ABOUT 0.3 ARCSEC WIDE BEFORE THE EXTRA INSET, AND 2 * 0.25 ARCSEC LEAVES NOTHING
    edges = {11: (-0.6, 0.6)}
    unusable = _result({11: (np.nan, STATUS_FAILED)}, edges, isApplied=False)
    usable = _result({11: (0.0, STATUS_FITTED)}, edges)

    assert flat_fit_ranges_from_result(unusable) == {}
    assert flat_fit_ranges_from_result(usable) == {11: _fit_range(-0.6, 0.6, 0.0)}


def test_the_flat_fit_range_is_inset_by_what_the_overlap_margin_leaves_short_of_the_clearance():
    # ARRANGE: AN OVERLAP MARGIN 0.15 ARCSEC SHORT OF THE CLEARANCE, AS A WORKSPACE THAT KEPT THE OLD 0.3 HAS
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)})

    # ACT
    ranges = flat_fit_ranges_from_result(result)

    # ASSERT
    keepLower, keepUpper = kept_slit_range(-5.0, 5.0, 1.0, MARGIN_ARCSEC)
    assert CLEARANCE_INSET_ARCSEC > 0
    assert ranges == {
        11: (pytest.approx(keepLower + CLEARANCE_INSET_ARCSEC), pytest.approx(keepUpper - CLEARANCE_INSET_ARCSEC))
    }


@pytest.mark.parametrize("margin", [FLAT_FIT_CLEARANCE_ARCSEC, FLAT_FIT_CLEARANCE_ARCSEC + 0.2])
def test_the_flat_fit_range_is_the_kept_range_when_the_overlap_margin_reaches_the_clearance(margin):
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)}, margin=margin)

    assert flat_fit_ranges_from_result(result) == {11: kept_slit_range(-5.0, 5.0, 1.0, margin)}


def test_the_flat_fit_range_leaves_out_an_order_the_clearance_inset_empties():
    result = _result({11: (0.0, STATUS_FITTED)}, {11: (-0.3, 0.3)}, margin=0.0)

    assert flat_fit_ranges_from_result(result) == {}


def test_flat_fit_offsets_report_the_offset_extra_inset_and_source_of_each_order():
    result = _result(
        {11: (0.1, STATUS_FITTED), 12: (0.9, STATUS_FITTED), 13: (np.nan, STATUS_FAILED)},
        {11: (-5.0, 5.0), 12: (-5.0, 5.0), 13: (-5.0, 5.0), 14: (np.nan, 5.0)},
    )

    assert flat_fit_offsets(result) == {
        11: (0.1, 0.0, "measured"),
        12: (0.9, 0.0, "measured"),
        13: (0.5, 0.0, "median"),
    }


def test_flat_fit_offsets_of_a_frame_with_no_usable_offset_are_zero_with_the_order_scatter_inset():
    result = _result({11: (np.nan, STATUS_FAILED)}, {11: (-5.0, 5.0), 12: (-5.0, 5.0)}, isApplied=False)

    assert flat_fit_offsets(result) == {
        11: (0.0, MAX_OFFSET_ERR_ARCSEC, "zero"),
        12: (0.0, MAX_OFFSET_ERR_ARCSEC, "zero"),
    }
