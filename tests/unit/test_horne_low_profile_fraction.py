"""DY-2195: the optimal extraction drops a slice whose unmasked pixels carry too little of the object profile."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import (
    compute_extractions,
    read_minimum_profile_fraction,
)

pytestmark = pytest.mark.unit

# A FIVE-ROW SLIT PROFILE, NORMALISED TO 1, WITH THE OBJECT IN ROWS 1 TO 4 AND A FAINT WING IN ROW 0
PROFILE = np.array([0.002, 0.2, 0.3, 0.3, 0.198])
AMPLITUDE = 1000.0
# A THRESHOLD ABOVE THE 0.05 DEFAULT, SO THE DEFAULT-THRESHOLD TEST CHECKS SOMETHING THE OTHERS DO NOT
THRESHOLD = 0.1


def _three_slices() -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Three slices at 500, 501 and 502 nm with the same profile and flux, and no pixel masked."""
    profile = np.tile(PROFILE[:, np.newaxis], (1, 3))
    flux = profile * AMPLITUDE
    # A NEGATIVE NOISE SPIKE IN THE FAINT WING: ONLY ROW 0 SURVIVES IN THE LOW-FRACTION SLICE
    flux[0, 1] = -50.0
    images = {
        "fluxRaw": flux,
        "objectProfile": profile,
        "variance": np.full_like(profile, 4.0),
        "wavelength": np.tile(np.array([500.0, 501.0, 502.0]), (PROFILE.size, 1)),
        "mask": np.zeros_like(profile, dtype=bool),
    }
    slices = pd.DataFrame({"order": np.full(3, 24), "pixelScaleNm": np.full(3, 0.1)})
    return slices, images


def test_compute_extractions_drops_a_slice_whose_profile_pixels_are_masked() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["mask"][1:, 1] = True

    # ACT
    result = compute_extractions(slices, images, order=24, minimumProfileFraction=THRESHOLD)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 502.0])


def test_compute_extractions_keeps_a_normal_slice_and_a_slice_with_one_masked_profile_pixel() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["fluxRaw"][0, 1] = PROFILE[0] * AMPLITUDE
    images["mask"][2, 1] = True

    # ACT
    result = compute_extractions(slices, images, order=24, minimumProfileFraction=THRESHOLD)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 501.0, 502.0])
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE] * 3, rtol=1e-12)
    np.testing.assert_allclose(result["extractedFluxBoxcarRobust"], [AMPLITUDE, 700.0, AMPLITUDE], rtol=1e-12)


def test_compute_extractions_keeps_every_slice_when_the_threshold_is_zero() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["mask"][1:, 1] = True

    # ACT
    result = compute_extractions(slices, images, order=24, minimumProfileFraction=0.0)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 501.0, 502.0])
    # THE SURVIVING WING PIXEL ALONE SETS THE OPTIMAL FLUX: -50 / 0.002
    assert result["extractedFluxOptimal"].iloc[1] == pytest.approx(-25000.0)


def test_compute_extractions_drops_a_low_profile_slice_at_the_default_threshold() -> None:
    # ARRANGE: ONLY THE 0.002 WING OF THE MIDDLE SLICE IS UNMASKED
    slices, images = _three_slices()
    images["mask"][1:, 1] = True

    # ACT
    result = compute_extractions(slices, images, order=24)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 502.0])


def test_read_minimum_profile_fraction_uses_the_default_when_the_key_is_absent() -> None:
    assert read_minimum_profile_fraction({}) == 0.05


def test_read_minimum_profile_fraction_returns_the_setting() -> None:
    assert read_minimum_profile_fraction({"horne-extraction-minimum-profile-fraction": 0.2}) == 0.2


@pytest.mark.parametrize("value", [-0.1, 1.0, 2, True, "0.1", None])
def test_read_minimum_profile_fraction_rejects_a_value_outside_zero_to_one(value: object) -> None:
    with pytest.raises(ValueError, match="horne-extraction-minimum-profile-fraction"):
        read_minimum_profile_fraction({"horne-extraction-minimum-profile-fraction": value})
