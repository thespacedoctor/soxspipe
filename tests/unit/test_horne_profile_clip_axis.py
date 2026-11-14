"""The Horne cosmic-ray clip runs along the slit, per wavelength column, weighted by distance from the centre (DY-1976).

In a rectified order image the rows are slit positions and the columns are wavelengths. A hot pixel far out on the slit
must be clipped, and a bright pixel on the trace at the slit centre must not be.
"""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.horne_extraction import _sigma_clip_and_mask

pytestmark = pytest.mark.unit

# THE 7 SIGMA UPPER LIMIT CAN ONLY BE REACHED IN A COLUMN OF ABOUT 50 ROWS OR MORE. THE SHIPPED SLIT LENGTHS ARE
# 10 TO 20 ROWS; THE STARE BLOCKS USE
SLIT_ROWS_SHIPPED = 15
SLIT_ROWS = 75
WAVELENGTH_COLUMNS = 200
CENTRE_ROW = SLIT_ROWS // 2
CENTRE_COLUMN = WAVELENGTH_COLUMNS // 2
FAR_ROW = 2
TRACE_PEAK = 1000.0
TRACE_SIGMA_ROWS = 5.0
SKY_LEVEL = 10.0
HOT_PIXEL = 5000.0


def _rectified_order() -> np.ndarray:
    """A smooth rectified order: a Gaussian trace at the slit centre on a flat sky, the same in every column."""
    rows = np.arange(SLIT_ROWS)
    profile = TRACE_PEAK * np.exp(-0.5 * ((rows - CENTRE_ROW) / TRACE_SIGMA_ROWS) ** 2) + SKY_LEVEL
    return np.tile(profile[:, None], (1, WAVELENGTH_COLUMNS))


def _clip(flux: np.ndarray) -> np.ndarray:
    return _sigma_clip_and_mask(flux, np.zeros(flux.shape, dtype=np.uint8))


def test_clean_rectified_order_has_no_clipped_pixels() -> None:
    mask = _clip(_rectified_order())

    assert not mask.any()


def test_hot_pixel_far_out_on_the_slit_is_clipped() -> None:
    flux = _rectified_order()
    flux[FAR_ROW, CENTRE_COLUMN] = HOT_PIXEL

    mask = _clip(flux)

    assert mask[FAR_ROW, CENTRE_COLUMN]
    assert mask.sum() == 1


def test_bright_pixel_at_the_slit_centre_on_the_trace_is_kept() -> None:
    flux = _rectified_order()
    flux[CENTRE_ROW, CENTRE_COLUMN] = 3 * TRACE_PEAK

    mask = _clip(flux)

    assert not mask[CENTRE_ROW, CENTRE_COLUMN]


def test_hot_pixel_is_not_clipped_on_a_rectified_slit_of_15_rows() -> None:
    """Pin the production geometry: the 7 sigma upper limit cannot be reached by one pixel in a 15-row column.

    The rectified order images of the shipped 15 px slit have 15 rows. One outlier among n values is at most
    sqrt(n - 1) standard deviations from their mean, 3.7 for n = 15, so the clip along the slit masks nothing.
    DY-1976 follow-up: lower ``sigma_upper`` (or lengthen the column) if the pre-fit clip is to catch cosmic rays.
    """
    flux = _rectified_order()[:SLIT_ROWS_SHIPPED]
    flux[1, CENTRE_COLUMN] = HOT_PIXEL

    mask = _clip(flux)

    assert not mask.any()
