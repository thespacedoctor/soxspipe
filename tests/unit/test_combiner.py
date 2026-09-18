"""Characterization contracts for the error-map-safe ``Combiner`` subclass."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.nddata import CCDData

from soxspipe.commonutils.combiner import Combiner

pytestmark = pytest.mark.unit

NAN = np.nan
STACK_RTOL = 1e-12


def _stack(dtype: type | None = None) -> Combiner:
    """Build a 3-frame 2x2 stack with one partly and one fully blank pixel."""
    frameA = np.array([[1.0, 2.0], [NAN, 4.0]])
    frameB = np.array([[3.0, NAN], [NAN, 8.0]])
    frameC = np.array([[5.0, 6.0], [NAN, NAN]])
    frames = [CCDData(f, unit="adu") for f in (frameA, frameB, frameC)]
    if dtype is None:
        return Combiner(frames)
    return Combiner(frames, dtype=dtype)


def test_average_combine_returns_nan_ignoring_mean_per_pixel() -> None:
    """Pixels average over the frames that carry a finite value."""
    combined = _stack().average_combine()

    expected = np.array([[3.0, 4.0], [NAN, 6.0]])
    np.testing.assert_allclose(
        combined.data, expected, rtol=STACK_RTOL, atol=0, equal_nan=True
    )


def test_average_combine_masks_only_pixels_blank_in_every_frame() -> None:
    """The mask flags the pixel with no finite value in any frame."""
    combined = _stack().average_combine()

    assert combined.mask.tolist() == [[False, False], [True, False]]


def test_average_combine_records_frame_count_and_unit() -> None:
    """NCOMBINE holds the stack depth and the unit is carried through."""
    combined = _stack().average_combine()

    assert combined.meta["NCOMBINE"] == 3
    assert combined.unit == "adu"


def test_average_combine_adds_no_synthetic_uncertainty() -> None:
    """The override leaves uncertainty unset, unlike upstream ccdproc."""
    assert _stack().average_combine().uncertainty is None


def test_average_combine_honours_combiner_dtype() -> None:
    """The combined array takes the dtype the combiner was built with."""
    assert _stack().average_combine().data.dtype == np.float64
    assert _stack(np.float32).average_combine().data.dtype == np.float32
