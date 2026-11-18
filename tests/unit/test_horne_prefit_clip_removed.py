"""DY-2223: the optimal extraction has no pre-fit cosmic-ray clip, and the Horne residual clip rejects the cosmic ray.

In a rectified order image the rows are slit positions and the columns are wavelengths. ``generate_masks`` masks the
flagged bad pixels and the NaN pixels only. The residual clip in ``compute_extractions`` (DY-2221) rejects cosmic rays
on any slit length. The shipped slit lengths are 10 to 20 rows; the old pre-fit clip needed about 50 rows to fire.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import compute_extractions, generate_masks

pytestmark = pytest.mark.unit

SLIT_ROWS_SHIPPED = 15
SLIT_ROWS_LONG = 75
WAVELENGTH_COLUMNS = 40
CENTRE_COLUMN = WAVELENGTH_COLUMNS // 2
TRACE_AMPLITUDE = 1000.0
TRACE_SIGMA_ROWS = 1.5
VARIANCE = 4.0
COSMIC_RAY = 5000.0
OFF_TRACE_ROW = 2
RESIDUAL_SIGMA = 5.0


def _trace_profile(rows: int) -> np.ndarray:
    """A normalised Gaussian slit profile centred on the middle row."""
    positions = np.arange(rows)
    profile = np.exp(-0.5 * ((positions - rows // 2) / TRACE_SIGMA_ROWS) ** 2)
    return profile / profile.sum()


def _order_images(rows: int, badPixels: np.ndarray | None = None) -> tuple[pd.DataFrame, dict[str, np.ndarray]]:
    """A noiseless rectified order: the same trace in every column, no sky, no bad pixel unless given."""
    profile = np.tile(_trace_profile(rows)[:, np.newaxis], (1, WAVELENGTH_COLUMNS))
    wavelengths = np.linspace(500.0, 510.0, WAVELENGTH_COLUMNS)
    images = {
        "bpMask": np.zeros(profile.shape, dtype=np.uint8) if badPixels is None else badPixels,
        "fluxRaw": profile * TRACE_AMPLITUDE,
        "objectProfile": profile,
        "variance": np.full_like(profile, VARIANCE),
        "wavelength": np.tile(wavelengths, (rows, 1)),
    }
    slices = pd.DataFrame({"order": np.full(WAVELENGTH_COLUMNS, 12)}, index=range(WAVELENGTH_COLUMNS))
    return slices, images


def test_generate_masks_leaves_a_cosmic_ray_on_a_long_slit_unmasked() -> None:
    # ARRANGE: 75 ROWS IS LONG ENOUGH FOR THE OLD 7 SIGMA PRE-FIT CLIP TO MASK A HOT PIXEL FAR FROM THE TRACE
    slices, images = _order_images(SLIT_ROWS_LONG)
    images["fluxRaw"][OFF_TRACE_ROW, CENTRE_COLUMN] += COSMIC_RAY

    # ACT
    _, result = generate_masks(slices, images)

    # ASSERT
    assert not result["mask"].any()


def test_generate_masks_masks_nan_pixels_and_every_nonzero_bad_pixel_value() -> None:
    # ARRANGE
    badPixels = np.zeros((SLIT_ROWS_SHIPPED, WAVELENGTH_COLUMNS), dtype=np.uint8)
    badPixels[3, 4] = 1
    badPixels[5, 6] = 2
    slices, images = _order_images(SLIT_ROWS_SHIPPED, badPixels)
    images["fluxRaw"][0, 0] = np.nan

    # ACT
    _, result = generate_masks(slices, images)

    # ASSERT
    assert {tuple(index) for index in np.argwhere(result["mask"])} == {(0, 0), (3, 4), (5, 6)}


def test_generate_masks_leaves_the_callers_bad_pixel_mask_unchanged() -> None:
    # ARRANGE
    badPixels = np.zeros((SLIT_ROWS_SHIPPED, WAVELENGTH_COLUMNS), dtype=np.uint8)
    badPixels[5, 6] = 2
    slices, images = _order_images(SLIT_ROWS_SHIPPED, badPixels)
    before = badPixels.copy()

    # ACT
    generate_masks(slices, images)

    # ASSERT
    np.testing.assert_array_equal(badPixels, before)


def test_residual_clip_rejects_an_off_trace_cosmic_ray_in_a_15_row_order() -> None:
    # ARRANGE: THE PRE-FIT CLIP COULD NOT FIRE ON 15 ROWS, SO THE RESIDUAL CLIP IS THE ONLY DEFENCE
    slices, images = _order_images(SLIT_ROWS_SHIPPED)
    images["fluxRaw"][OFF_TRACE_ROW, CENTRE_COLUMN] += COSMIC_RAY

    # ACT
    slices, images = generate_masks(slices, images)
    result = compute_extractions(slices, images, order=12, residualClippingSigma=RESIDUAL_SIGMA)

    # ASSERT
    assert not images["mask"][OFF_TRACE_ROW, CENTRE_COLUMN]
    np.testing.assert_allclose(result["extractedFluxOptimal"], TRACE_AMPLITUDE, rtol=1e-9)
