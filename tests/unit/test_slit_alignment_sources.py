"""Where a slit-shifted pixel takes its value from, on a tilted and curved 2D map (DY-1692)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.nddata import CCDData

from soxspipe.commonutils.slit_alignment import (
    _interpolate,
    _newton_step,
    _order_sources,
    build_two_d_map_images,
    shift_frame_along_slit,
)
from tests.factories.slit_alignment import PIXEL_ARCSEC

pytestmark = pytest.mark.unit

SHAPE = (60, 120)
CENTRE_COLUMN = 59.5
# THE SLIT PLATE SCALE GROWS ALONG THE SLIT, AS AT THE RED END OF THE SOXS NIR ORDER 10
STRETCH = 0.0017
# ARCSEC PER ROW OF THE SLIT ISO-LINES' TILT
SLIT_PER_ROW = 0.01
# THE WAVELENGTH CHANGES BY THIS MANY NM PER ROW, PER COLUMN, AND PER COLUMN SQUARED, AND BY THE CROSS TERM
NM_PER_ROW = 0.02
NM_PER_COLUMN = 0.004
NM_PER_COLUMN_SQUARED = 4e-5
NM_PER_ROW_COLUMN = 2e-5
# A SHIFT OF SEVERAL ARCSEC, ABOUT 15 COLUMNS, LIKE A NOD THROW
SHIFT_ARCSEC = 3.0
# THE WAVELENGTH TOLERANCE IS ABOUT 1/40 OF A ROW AND ABOUT THE FLOAT32 PRECISION OF THE WAVELENGTH IMAGE
WAVELENGTH_TOLERANCE_NM = 5e-4
SLIT_TOLERANCE_ARCSEC = 2e-3


def _tilted_curved_map():
    """A one-order map whose wavelength depends on the column and whose slit coordinate is curved and tilted."""
    rows, columns = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]]
    offset = columns - CENTRE_COLUMN
    slit = PIXEL_ARCSEC * (offset + STRETCH * offset**2) + SLIT_PER_ROW * (rows - SHAPE[0] / 2)
    wavelength = (
        1000.0
        + NM_PER_ROW * rows
        + NM_PER_COLUMN * offset
        + NM_PER_COLUMN_SQUARED * offset**2
        + NM_PER_ROW_COLUMN * rows * offset
    )
    return build_two_d_map_images(wavelength, slit, np.ones(SHAPE))


def _interior_pixels(mapImages, targets, sources):
    """The mask of target pixels whose source lies well inside the detector."""
    rows, columns = sources
    return (rows > 2) & (rows < SHAPE[0] - 3) & (columns > 2) & (columns < SHAPE[1] - 3)


def test_a_shifted_source_keeps_the_target_wavelength_on_a_tilted_curved_map() -> None:
    # ARRANGE: THE WAVELENGTH GRADIENT HAS A COLUMN COMPONENT, SO A SLIT-ONLY STEP DRIFTS IN WAVELENGTH
    mapImages = _tilted_curved_map()

    # ACT
    targets, sources = _order_sources(mapImages, 1, SHIFT_ARCSEC)

    # ASSERT
    isInterior = _interior_pixels(mapImages, targets, sources)
    wavelengthError = _interpolate(mapImages.wavelength, sources) - mapImages.wavelength[targets]
    assert isInterior.sum() > 1000
    assert np.abs(wavelengthError[isInterior]).max() < WAVELENGTH_TOLERANCE_NM


def test_a_shifted_source_lies_the_shift_above_the_target_in_slit_coordinate_on_a_tilted_curved_map() -> None:
    # ARRANGE
    mapImages = _tilted_curved_map()

    # ACT
    targets, sources = _order_sources(mapImages, 1, SHIFT_ARCSEC)

    # ASSERT
    isInterior = _interior_pixels(mapImages, targets, sources)
    slitError = _interpolate(mapImages.slit, sources) - (mapImages.slit[targets] + SHIFT_ARCSEC)
    assert np.abs(slitError[isInterior]).max() < SLIT_TOLERANCE_ARCSEC


def test_a_pixel_with_parallel_wavelength_and_slit_gradients_is_flagged_bad_not_given_a_finite_wrong_value() -> None:
    # ARRANGE: THE WAVELENGTH AND SLIT IMAGES VARY ONLY ALONG THE COLUMNS, SO THE 2X2 JACOBIAN IS SINGULAR
    rows, columns = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]]
    mapImages = build_two_d_map_images(
        1000.0 + 0.02 * columns, PIXEL_ARCSEC * (columns - CENTRE_COLUMN), np.ones(SHAPE)
    )
    frame = CCDData(np.ones(SHAPE), unit="electron", mask=np.zeros(SHAPE, dtype=bool))

    # ACT
    with np.errstate(all="raise"):
        aligned = shift_frame_along_slit(frame, {1: SHIFT_ARCSEC}, mapImages)

    # ASSERT
    assert aligned.mask.all()


def test_a_well_conditioned_map_with_a_slit_aligned_to_the_columns_still_shifts_by_the_slit_step() -> None:
    # ARRANGE: THE SLIT RUNS ALONG THE COLUMNS AND THE WAVELENGTH ALONG THE ROWS, SO THE JACOBIAN IS DIAGONAL
    rows, columns = np.mgrid[0 : SHAPE[0], 0 : SHAPE[1]]
    mapImages = build_two_d_map_images(1000.0 + 0.02 * rows, PIXEL_ARCSEC * (columns - CENTRE_COLUMN), np.ones(SHAPE))

    # ACT
    targets, (sourceRows, sourceColumns) = _order_sources(mapImages, 1, SHIFT_ARCSEC)

    # ASSERT
    expectedColumns = targets[1] + SHIFT_ARCSEC / PIXEL_ARCSEC
    isOnDetector = expectedColumns < SHAPE[1] - 3
    np.testing.assert_allclose(sourceRows[isOnDetector], targets[0][isOnDetector], atol=1e-4)
    np.testing.assert_allclose(sourceColumns[isOnDetector], expectedColumns[isOnDetector], rtol=1e-5)


def test_a_newton_step_solves_both_residuals_with_the_inverse_of_the_gradient_matrix() -> None:
    # ARRANGE: A TILTED JACOBIAN WITH A WAVELENGTH COLUMN GRADIENT, AND BOTH RESIDUALS NON-ZERO
    gradients = [np.array([0.02]), np.array([0.004]), np.array([0.01]), np.array([0.2])]
    jacobian = np.array([[0.02, 0.004], [0.01, 0.2]])
    residuals = np.array([0.03, 1.5])

    # ACT
    newRows, newColumns = _newton_step(np.array([10.0]), np.array([20.0]), gradients, residuals[0], residuals[1])

    # ASSERT
    expected = np.linalg.solve(jacobian, residuals)
    np.testing.assert_allclose([newRows[0] - 10.0, newColumns[0] - 20.0], expected, rtol=1e-12)


def test_a_newton_step_ignores_the_wavelength_residual_where_the_gradient_matrix_is_ill_conditioned() -> None:
    # ARRANGE: THE SLIT GRADIENT IS PARALLEL TO THE WAVELENGTH GRADIENT TO WITHIN 1E-4, BELOW THE CONDITIONING LIMIT
    gradients = [np.array([0.02]), np.array([0.004]), np.array([0.02]), np.array([0.004 + 2e-6])]
    expectedRows, expectedColumns = _newton_step(np.array([10.0]), np.array([20.0]), gradients, 0.0, 1.5)

    # ACT
    newRows, newColumns = _newton_step(np.array([10.0]), np.array([20.0]), gradients, 0.03, 1.5)

    # ASSERT
    np.testing.assert_array_equal(newRows, expectedRows)
    np.testing.assert_array_equal(newColumns, expectedColumns)


def test_a_newton_step_gives_a_non_finite_coordinate_not_a_finite_wrong_one_for_a_zero_determinant() -> None:
    # ARRANGE
    gradients = [np.array([0.02]), np.array([0.004]), np.array([0.02]), np.array([0.004])]

    # ACT
    with np.errstate(all="raise"):
        newRows, newColumns = _newton_step(np.array([10.0]), np.array([20.0]), gradients, 0.03, 1.5)

    # ASSERT
    assert not np.isfinite(newRows).any()
    assert not np.isfinite(newColumns).any()
