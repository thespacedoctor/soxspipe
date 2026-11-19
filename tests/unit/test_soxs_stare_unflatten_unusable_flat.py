"""The stare recipe's unflattened frame is NaN and masked where the master flat is unusable (DY-2279)."""

from __future__ import annotations

import warnings

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
from tests.unit.test_soxs_stare_unflatten import _recipe

pytestmark = pytest.mark.unit

SKY_SUBTRACTED_VALUE = 10.0
SKY_SUBTRACTED_ERROR = 1.0


def _unflatten_without_warnings(log) -> CCDData:
    skySubtracted = CCDData(
        np.full(FLAT_SHAPE, SKY_SUBTRACTED_VALUE, dtype=np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full(FLAT_SHAPE, SKY_SUBTRACTED_ERROR, dtype=np.float32)),
        mask=np.zeros(FLAT_SHAPE, dtype=bool),
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        return _recipe(log)._unflatten_sky_subtracted_frame(skySubtracted, flat_with_unusable_pixels(), None)


def test_unflattening_sets_the_data_and_error_to_nan_where_the_flat_is_unusable(log) -> None:
    # ARRANGE
    unusable = unusable_flat_pixel_mask()

    # ACT
    result = _unflatten_without_warnings(log)

    # ASSERT
    assert np.isnan(result.data[unusable]).all()
    assert np.isnan(result.uncertainty.array[unusable]).all()


def test_unflattening_masks_the_pixels_where_the_flat_is_unusable(log) -> None:
    # ARRANGE
    unusable = unusable_flat_pixel_mask()

    # ACT
    result = _unflatten_without_warnings(log)

    # ASSERT
    assert result.mask[unusable].all()


def test_unflattening_still_multiplies_the_pixels_where_the_flat_is_usable(log) -> None:
    # ARRANGE
    usable = ~unusable_flat_pixel_mask()

    # ACT
    result = _unflatten_without_warnings(log)

    # ASSERT
    np.testing.assert_allclose(result.data[usable], SKY_SUBTRACTED_VALUE * USABLE_FLAT_VALUE)
    np.testing.assert_allclose(result.uncertainty.array[usable], SKY_SUBTRACTED_ERROR * USABLE_FLAT_VALUE)
