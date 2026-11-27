"""``detrend`` sets the pixels where the master flat is unusable to NaN and masks them (DY-2279)."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from tests.factories.unusable_flat import (
    FLAT_SHAPE,
    USABLE_FLAT_VALUE,
    flat_with_unusable_pixels,
    unusable_flat_pixel_mask,
)
from tests.unit.test_base_recipe_detrend import _recipe

pytestmark = pytest.mark.unit

FRAME_VALUE = 10.0
FRAME_ERROR = 1.0

# A REAL VIS INTER-ORDER MASTER-FLAT PIXEL (20251012T101822 MFLAT, ROW 3667, COLUMN 27) AND A FRAME VALUE AND ERRORS OF
# THE SIZE A NOD DIFFERENCE AND THE FLAT CARRY THERE (DY-2440)
TINY_FLAT_VALUE = 1.4551915e-12
TINY_FLAT_ERROR = 2.2e-4
BRIGHT_FRAME_VALUE = 1000.0
BRIGHT_FRAME_ERROR = 30.0


def _frame_to_flat_correct() -> CCDData:
    return CCDData(
        np.full(FLAT_SHAPE, FRAME_VALUE, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full(FLAT_SHAPE, FRAME_ERROR, dtype=np.float32)),
        mask=np.zeros(FLAT_SHAPE, dtype=bool),
    )


def _detrend_without_warnings(log: Any) -> CCDData:
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        return _recipe(log).detrend(_frame_to_flat_correct(), master_flat=flat_with_unusable_pixels())


def test_detrend_sets_the_data_and_error_to_nan_where_the_flat_is_unusable(log: Any) -> None:
    # ARRANGE
    unusable = unusable_flat_pixel_mask()

    # ACT
    result = _detrend_without_warnings(log)

    # ASSERT
    assert np.isnan(result.data[unusable]).all()
    assert np.isnan(result.uncertainty.array[unusable]).all()


def test_detrend_masks_the_pixels_where_the_flat_is_unusable(log: Any) -> None:
    # ARRANGE
    unusable = unusable_flat_pixel_mask()

    # ACT
    result = _detrend_without_warnings(log)

    # ASSERT
    assert result.mask[unusable].all()


def test_detrend_still_divides_the_pixels_where_the_flat_is_usable(log: Any) -> None:
    # ARRANGE
    usable = ~unusable_flat_pixel_mask()

    # ACT
    result = _detrend_without_warnings(log)

    # ASSERT
    np.testing.assert_allclose(result.data[usable], FRAME_VALUE / USABLE_FLAT_VALUE)
    assert not result.mask[usable & ~flat_with_unusable_pixels().mask].any()


def _detrend_with_a_tiny_positive_flat_pixel(log: Any) -> CCDData:
    flat = CCDData(
        np.array([[USABLE_FLAT_VALUE, TINY_FLAT_VALUE]], dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full((1, 2), TINY_FLAT_ERROR, dtype=np.float32)),
        mask=np.zeros((1, 2), dtype=bool),
    )
    frame = CCDData(
        np.full((1, 2), BRIGHT_FRAME_VALUE, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full((1, 2), BRIGHT_FRAME_ERROR, dtype=np.float32)),
        mask=np.zeros((1, 2), dtype=bool),
    )
    return _recipe(log).detrend(frame, master_flat=flat)


def test_detrend_masks_a_pixel_where_the_flat_is_a_tiny_positive_residue(log: Any) -> None:
    # ARRANGE
    tinyFlatPixel = (0, 1)

    # ACT
    result = _detrend_with_a_tiny_positive_flat_pixel(log)

    # ASSERT
    assert result.mask[tinyFlatPixel]
    assert np.isnan(result.uncertainty.array[tinyFlatPixel])
    assert not result.mask[0, 0]


def test_the_float32_variance_of_a_frame_detrended_by_a_tiny_flat_pixel_does_not_overflow(log: Any) -> None:
    """The nod slit alignment and shift-and-add square the float32 uncertainty of a detrended frame (DY-2440)."""
    # ARRANGE
    result = _detrend_with_a_tiny_positive_flat_pixel(log)

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        variance = np.asarray(result.uncertainty.array, dtype=np.float32) ** 2

    # ASSERT
    assert not np.isinf(variance).any()
