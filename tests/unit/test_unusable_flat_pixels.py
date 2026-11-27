"""The shared rule for unusable master-flat pixels (DY-2279)."""

from __future__ import annotations

import warnings

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.unusable_flat_pixels import (
    flag_unusable_pixels,
    neutralise_unusable_flat,
    unusable_flat_pixel_mask,
)
from tests.factories.unusable_flat import flat_with_unusable_pixels
from tests.factories.unusable_flat import unusable_flat_pixel_mask as expected_unusable

pytestmark = pytest.mark.unit


def test_a_flat_pixel_is_unusable_when_it_is_zero_negative_or_not_finite() -> None:
    # ARRANGE
    flatData = np.array([1.0, 1e-6, 0.0, -0.0, -0.5, np.nan, np.inf, -np.inf])

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        unusable = unusable_flat_pixel_mask(flatData)

    # ASSERT
    np.testing.assert_array_equal(unusable, [False, False, True, True, True, True, True, True])


def test_a_flat_pixel_below_float32_epsilon_is_unusable() -> None:
    """An inter-order master-flat pixel holds background-subtracted noise, so some hold a tiny positive residue.

    Dividing a frame by one gives an uncertainty above 1.8e19, whose float32 square overflows (DY-2440).
    """
    # ARRANGE
    epsilon = np.finfo(np.float32).eps
    flatData = np.array([1.46e-12, 7.45e-9, epsilon / 2, epsilon, 2 * epsilon, 1e-6], dtype=np.float32)

    # ACT
    unusable = unusable_flat_pixel_mask(flatData)

    # ASSERT
    np.testing.assert_array_equal(unusable, [True, True, True, False, False, False])


def test_neutralising_the_flat_sets_the_unusable_pixels_to_one_with_no_error() -> None:
    # ARRANGE
    flat = flat_with_unusable_pixels()
    unusable = expected_unusable()

    # ACT
    neutral = neutralise_unusable_flat(flat, unusable)

    # ASSERT
    np.testing.assert_array_equal(neutral.data[unusable], 1.0)
    np.testing.assert_array_equal(neutral.uncertainty.array[unusable], 0.0)
    np.testing.assert_array_equal(neutral.data[~unusable], flat.data[~unusable])
    np.testing.assert_array_equal(neutral.uncertainty.array[~unusable], flat.uncertainty.array[~unusable])


def test_neutralising_the_flat_leaves_the_original_flat_unchanged() -> None:
    # ARRANGE
    flat = flat_with_unusable_pixels()
    dataBefore = flat.data.copy()
    uncertaintyBefore = flat.uncertainty.array.copy()
    maskBefore = flat.mask.copy()

    # ACT
    neutralise_unusable_flat(flat, expected_unusable())

    # ASSERT
    np.testing.assert_array_equal(flat.data, dataBefore)
    np.testing.assert_array_equal(flat.uncertainty.array, uncertaintyBefore)
    np.testing.assert_array_equal(flat.mask, maskBefore)


def test_neutralising_a_flat_with_no_uncertainty_sets_only_the_data() -> None:
    # ARRANGE
    flat = CCDData(np.array([[1.0, 0.0]]), unit=u.electron)

    # ACT
    neutral = neutralise_unusable_flat(flat, np.array([[False, True]]))

    # ASSERT
    np.testing.assert_array_equal(neutral.data, [[1.0, 1.0]])
    assert neutral.uncertainty is None


def test_flagging_sets_data_and_uncertainty_to_nan_and_masks_the_unusable_pixels() -> None:
    # ARRANGE
    frame = CCDData(
        np.full((1, 3), 5.0),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full((1, 3), 0.5)),
        mask=np.array([[True, False, False]]),
    )
    unusable = np.array([[False, True, False]])

    # ACT
    flagged = flag_unusable_pixels(frame, unusable)

    # ASSERT
    assert np.isnan(flagged.data[0, 1])
    assert np.isnan(flagged.uncertainty.array[0, 1])
    np.testing.assert_array_equal(flagged.mask, [[True, True, False]])
    np.testing.assert_array_equal(flagged.data[0, [0, 2]], 5.0)
    np.testing.assert_array_equal(flagged.uncertainty.array[0, [0, 2]], 0.5)


def test_flagging_a_frame_with_no_mask_or_uncertainty_creates_the_mask() -> None:
    # ARRANGE
    frame = CCDData(np.full((1, 2), 5.0), unit=u.electron)

    # ACT
    flagged = flag_unusable_pixels(frame, np.array([[False, True]]))

    # ASSERT
    np.testing.assert_array_equal(flagged.mask, [[False, True]])
    assert np.isnan(flagged.data[0, 1])
    assert flagged.uncertainty is None


def test_flagging_leaves_the_original_frame_unchanged() -> None:
    # ARRANGE
    frame = CCDData(
        np.full((1, 3), 5.0),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full((1, 3), 0.5)),
        mask=np.array([[True, False, False]]),
    )

    # ACT
    flagged = flag_unusable_pixels(frame, np.array([[False, True, False]]))

    # ASSERT
    assert flagged is not frame
    np.testing.assert_array_equal(frame.data, [[5.0, 5.0, 5.0]])
    np.testing.assert_array_equal(frame.uncertainty.array, [[0.5, 0.5, 0.5]])
    np.testing.assert_array_equal(frame.mask, [[True, False, False]])


def test_flagging_an_integer_frame_promotes_it_to_float_and_leaves_the_original_integer() -> None:
    # ARRANGE
    frame = CCDData(
        np.array([[5, 6, 7]], dtype=np.int32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.array([[1, 2, 3]], dtype=np.int32)),
    )

    # ACT
    flagged = flag_unusable_pixels(frame, np.array([[False, True, False]]))

    # ASSERT
    assert flagged.data.dtype.kind == "f"
    assert flagged.uncertainty.array.dtype.kind == "f"
    assert np.isnan(flagged.data[0, 1]) and np.isnan(flagged.uncertainty.array[0, 1])
    np.testing.assert_array_equal(flagged.data[0, [0, 2]], [5.0, 7.0])
    np.testing.assert_array_equal(flagged.uncertainty.array[0, [0, 2]], [1.0, 3.0])
    assert frame.data.dtype == np.int32
    np.testing.assert_array_equal(frame.data, [[5, 6, 7]])
    np.testing.assert_array_equal(frame.uncertainty.array, [[1, 2, 3]])
