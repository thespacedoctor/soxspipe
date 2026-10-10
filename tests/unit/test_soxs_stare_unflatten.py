"""The stare recipe's un-flattened sky-subtracted frame carries the error of un-flattened data (DY-1391)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.soxs_stare import soxs_stare

pytestmark = pytest.mark.unit

# A NON-UNIT FLAT THAT VARIES ACROSS THE FRAME, SMALL AT ONE CORNER LIKE AN ORDER EDGE
FLAT = np.array([[0.5, 2.0], [4.0, 0.25]])
FLAT_ERROR = np.array([[0.05, 0.1], [0.2, 0.02]])
SKY_SUBTRACTED = np.array([[10.0, 20.0], [30.0, 40.0]])
SKY_SUBTRACTED_ERROR = np.array([[1.0, 2.0], [3.0, 4.0]])


def _recipe(log: Any, *, subtractSky: bool = True) -> soxs_stare:
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.subtractSky = subtractSky
    return recipe


def _sky_subtracted_and_varying_flat() -> tuple[CCDData, CCDData]:
    skySubtracted = CCDData(
        SKY_SUBTRACTED.copy(),
        unit=u.electron,
        uncertainty=StdDevUncertainty(SKY_SUBTRACTED_ERROR.copy()),
        mask=np.array([[False, True], [False, False]]),
    )
    skySubtracted.header["MARKER"] = "skysub"
    masterFlat = CCDData(
        FLAT.copy(),
        unit=u.electron,
        uncertainty=StdDevUncertainty(FLAT_ERROR.copy()),
        mask=np.array([[False, False], [True, False]]),
    )
    return skySubtracted, masterFlat


def test_the_unflattened_error_is_the_flat_times_the_flat_corrected_error(log: Any) -> None:
    # ARRANGE
    skySubtracted, masterFlat = _sky_subtracted_and_varying_flat()

    # ACT
    result = _recipe(log)._unflatten_sky_subtracted_frame(skySubtracted, masterFlat, None)

    # ASSERT
    assert isinstance(result.uncertainty, StdDevUncertainty)
    np.testing.assert_allclose(result.uncertainty.array, [[0.5, 4.0], [12.0, 1.0]])


def test_the_unflattened_frame_carries_the_multiplied_data_the_combined_mask_and_the_sky_subtracted_header(
    log: Any,
) -> None:
    # ARRANGE
    skySubtracted, masterFlat = _sky_subtracted_and_varying_flat()

    # ACT
    result = _recipe(log)._unflatten_sky_subtracted_frame(skySubtracted, masterFlat, None)

    # ASSERT
    np.testing.assert_allclose(result.data, [[5.0, 40.0], [120.0, 10.0]])
    np.testing.assert_array_equal(result.mask, [[False, True], [True, False]])
    assert result.header["MARKER"] == "skysub"


def test_unflattening_leaves_the_sky_subtracted_frame_and_the_flat_unchanged(log: Any) -> None:
    # ARRANGE
    skySubtracted, masterFlat = _sky_subtracted_and_varying_flat()

    # ACT
    _recipe(log)._unflatten_sky_subtracted_frame(skySubtracted, masterFlat, None)

    # ASSERT
    np.testing.assert_allclose(skySubtracted.data, SKY_SUBTRACTED)
    np.testing.assert_allclose(skySubtracted.uncertainty.array, SKY_SUBTRACTED_ERROR)
    np.testing.assert_allclose(masterFlat.data, FLAT)
    np.testing.assert_allclose(masterFlat.uncertainty.array, FLAT_ERROR)


def test_without_a_flat_the_sky_subtracted_frame_is_returned_as_it_is(log: Any) -> None:
    # ARRANGE
    skySubtracted, _ = _sky_subtracted_and_varying_flat()

    # ACT
    result = _recipe(log)._unflatten_sky_subtracted_frame(skySubtracted, False, None)

    # ASSERT
    assert result is skySubtracted


def test_without_sky_subtraction_the_stacked_frame_is_returned(log: Any) -> None:
    # ARRANGE
    skySubtracted, masterFlat = _sky_subtracted_and_varying_flat()
    stacked = CCDData(np.zeros(FLAT.shape), unit=u.electron)

    # ACT
    result = _recipe(log, subtractSky=False)._unflatten_sky_subtracted_frame(skySubtracted, masterFlat, stacked)

    # ASSERT
    assert result is stacked
