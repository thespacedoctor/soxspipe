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
