"""Contracts for ``slit_illumination_mask.kept_ranges_from_result``, the both-lit range of each order (DY-2147)."""

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
    flat_fit_ranges_from_result,
    kept_ranges_from_result,
    kept_slit_range,
)

pytestmark = pytest.mark.unit

MARGIN_ARCSEC = 0.3


def _result(offsets, edges, isApplied=True, margin=MARGIN_ARCSEC):
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
    )


def test_returns_no_ranges_when_the_measurement_was_not_applied():
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)}, isApplied=False)

    assert kept_ranges_from_result(result) == {}


def test_returns_the_kept_slit_range_of_each_order_with_its_own_offset():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (-0.5, STATUS_FITTED)},
        {11: (-5.0, 5.0), 12: (-5.5, 5.5)},
    )

    ranges = kept_ranges_from_result(result)

    assert ranges == {
        11: kept_slit_range(-5.0, 5.0, 1.0, MARGIN_ARCSEC),
        12: kept_slit_range(-5.5, 5.5, -0.5, MARGIN_ARCSEC),
    }


def test_uses_a_zero_offset_for_an_order_whose_offset_is_not_usable():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (np.nan, STATUS_FAILED)},
        {11: (-5.0, 5.0), 12: (-5.5, 5.5), 13: (-5.0, 5.0)},
    )

    ranges = kept_ranges_from_result(result)

    assert ranges[12] == kept_slit_range(-5.5, 5.5, 0.0, MARGIN_ARCSEC)
    assert ranges[13] == kept_slit_range(-5.0, 5.0, 0.0, MARGIN_ARCSEC)


def test_leaves_out_an_order_whose_kept_range_is_empty_or_whose_flat_edges_are_not_finite():
    result = _result(
        {11: (1.0, STATUS_FITTED), 12: (0.0, STATUS_FITTED)},
        {11: (-0.1, 0.1), 12: (np.nan, 5.0)},
    )

    assert kept_ranges_from_result(result) == {}


def test_the_flat_fit_range_is_inset_by_what_the_overlap_margin_leaves_short_of_the_clearance():
    # ARRANGE: AN OVERLAP MARGIN 0.15 ARCSEC SHORT OF THE CLEARANCE, AS A WORKSPACE THAT KEPT THE OLD 0.3 HAS
    inset = FLAT_FIT_CLEARANCE_ARCSEC - MARGIN_ARCSEC
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)})

    # ACT
    ranges = flat_fit_ranges_from_result(result)

    # ASSERT
    keepLower, keepUpper = kept_slit_range(-5.0, 5.0, 1.0, MARGIN_ARCSEC)
    assert inset > 0
    assert ranges == {11: (pytest.approx(keepLower + inset), pytest.approx(keepUpper - inset))}


@pytest.mark.parametrize("margin", [FLAT_FIT_CLEARANCE_ARCSEC, FLAT_FIT_CLEARANCE_ARCSEC + 0.2])
def test_the_flat_fit_range_is_the_kept_range_when_the_overlap_margin_reaches_the_clearance(margin):
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)}, margin=margin)

    assert flat_fit_ranges_from_result(result) == kept_ranges_from_result(result)


def test_the_flat_fit_range_leaves_out_an_order_the_inset_empties():
    result = _result({11: (0.0, STATUS_FITTED)}, {11: (-0.3, 0.3)}, margin=0.0)

    assert kept_ranges_from_result(result) != {}
    assert flat_fit_ranges_from_result(result) == {}


def test_the_flat_fit_range_is_empty_when_the_measurement_was_not_applied():
    result = _result({11: (1.0, STATUS_FITTED)}, {11: (-5.0, 5.0)}, isApplied=False)

    assert flat_fit_ranges_from_result(result) == {}
