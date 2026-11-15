"""DY-2221: the optimal extraction rejects a pixel whose residual from the profile model is too large."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import (
    compute_extractions,
    read_residual_clipping_sigma,
)

pytestmark = pytest.mark.unit

# A FIVE-ROW SLIT PROFILE, NORMALISED TO 1, WITH THE OBJECT IN ROWS 1 TO 4 AND A FAINT WING IN ROW 0
PROFILE = np.array([0.002, 0.2, 0.3, 0.3, 0.198])
AMPLITUDE = 1000.0
# EVERY PIXEL HAS A VARIANCE OF 4, SO ONE STANDARD DEVIATION IS 2 COUNTS
VARIANCE = 4.0
COSMIC_RAY = 500.0
SIGMA = 5.0


def _three_slices() -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """Three noiseless slices at 500, 501 and 502 nm with the same profile and flux, and no pixel masked."""
    profile = np.tile(PROFILE[:, np.newaxis], (1, 3))
    images = {
        "fluxRaw": profile * AMPLITUDE,
        "objectProfile": profile,
        "variance": np.full_like(profile, VARIANCE),
        "wavelength": np.tile(np.array([500.0, 501.0, 502.0]), (PROFILE.size, 1)),
        "mask": np.zeros_like(profile, dtype=bool),
    }
    slices = pd.DataFrame({"order": np.full(3, 24), "pixelScaleNm": np.full(3, 0.1)})
    return slices, images


def test_compute_extractions_rejects_a_cosmic_ray_on_the_trace() -> None:
    # ARRANGE: WITHOUT THE CLIP THE MIDDLE SLICE WOULD GIVE 1000 + 0.3 * 500 / 0.259208 = 1578.7
    slices, images = _three_slices()
    images["fluxRaw"][2, 1] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE] * 3, rtol=1e-12)


def test_compute_extractions_keeps_a_pixel_below_the_threshold() -> None:
    # ARRANGE: A +8 COUNT (4 SIGMA) EXCESS IN ROW 2 OF THE MIDDLE SLICE STAYS IN AT 5 SIGMA
    slices, images = _three_slices()
    images["fluxRaw"][2, 1] += 8.0

    # ACT
    result = compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT: HORNE FLUX WITH THE EXCESS KEPT, 1000 + 0.3 * 8 / 0.259208
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE, 1009.2589735, AMPLITUDE], rtol=1e-9)


def test_compute_extractions_drops_a_slice_that_the_clip_leaves_with_too_little_profile() -> None:
    # ARRANGE: THE MIDDLE SLICE KEEPS ROWS 0 TO 2 (PROFILE FRACTION 0.502) AND ROW 2 HOLDS A COSMIC RAY. ONE CLIPPED
    # PIXEL LEAVES AT MOST 0.302 OF THE PROFILE, BELOW THE 0.35 THRESHOLD
    slices, images = _three_slices()
    images["mask"][3:, 1] = True
    images["fluxRaw"][2, 1] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24, minimumProfileFraction=0.35, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 502.0])


def test_compute_extractions_keeps_the_cosmic_ray_in_the_boxcar_columns() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["fluxRaw"][2, 1] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_allclose(result["extractedFluxBoxcarRobust"], [AMPLITUDE, 1500.0, AMPLITUDE], rtol=1e-12)


def test_compute_extractions_clips_at_the_default_threshold() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["fluxRaw"][2, 1] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24)

    # ASSERT
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE] * 3, rtol=1e-12)


def test_read_residual_clipping_sigma_uses_the_default_when_the_key_is_absent() -> None:
    assert read_residual_clipping_sigma({}) == 5.0


def test_read_residual_clipping_sigma_returns_the_setting() -> None:
    assert read_residual_clipping_sigma({"horne-extraction-residual-clipping-sigma": 4}) == 4.0


@pytest.mark.parametrize("value", [0, -1.0, float("nan"), float("inf"), True, "5", None])
def test_read_residual_clipping_sigma_rejects_a_value_that_is_not_a_positive_number(value: object) -> None:
    with pytest.raises(ValueError, match="horne-extraction-residual-clipping-sigma"):
        read_residual_clipping_sigma({"horne-extraction-residual-clipping-sigma": value})


def test_compute_extractions_rejects_two_cosmic_rays_in_one_column() -> None:
    # ARRANGE: ONE PIXEL IS REJECTED PER PASS, SO THE SECOND COSMIC RAY NEEDS A SECOND PASS. THE HITS ARE IN ROWS 1 AND
    # 4: TWO HITS ON THE TWO BRIGHTEST ROWS (60% OF THE PROFILE) PULL THE FIRST FLUX SO HIGH THAT A CLEAN ROW HAS THE
    # LARGEST RESIDUAL, AND THE CLIP CANNOT REPAIR THE COLUMN
    slices, images = _three_slices()
    images["fluxRaw"][1, 1] += COSMIC_RAY
    images["fluxRaw"][4, 1] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE] * 3, rtol=1e-12)


def test_compute_extractions_leaves_a_column_with_a_non_finite_variance_alone() -> None:
    # ARRANGE: THE MIDDLE SLICE HAS A NaN VARIANCE, SO IT HAS NO OPTIMAL FLUX; THE COSMIC RAY BESIDE IT IS STILL CLIPPED
    slices, images = _three_slices()
    images["variance"][0, 1] = np.nan
    images["fluxRaw"][2, 2] += COSMIC_RAY

    # ACT
    result = compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_array_equal(result["wavelengthMean"], [500.0, 502.0])
    np.testing.assert_allclose(result["extractedFluxOptimal"], [AMPLITUDE] * 2, rtol=1e-12)


def test_compute_extractions_does_not_change_the_input_mask() -> None:
    # ARRANGE
    slices, images = _three_slices()
    images["fluxRaw"][2, 1] += COSMIC_RAY
    maskBefore = images["mask"].copy()

    # ACT
    compute_extractions(slices, images, order=24, residualClippingSigma=SIGMA)

    # ASSERT
    np.testing.assert_array_equal(images["mask"], maskBefore)
