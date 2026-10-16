"""Resampling a frame along the slit through the 2D map (DY-1612)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy.nddata import CCDData, StdDevUncertainty, VarianceUncertainty

from soxspipe.commonutils.slit_alignment import build_two_d_map_images, shift_frame_along_slit
from tests.factories.slit_alignment import PIXEL_ARCSEC, synthetic_alignment_map, synthetic_nod_frame

pytestmark = pytest.mark.unit

# ROWS AWAY FROM THE ENDS OF THE ORDERS, WHERE THE OBJECT PROFILE IS COMPLETE AFTER A SHIFT
MIDDLE_ROWS = slice(200, 400)


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def two_d_map(scene):
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def _centroid(scene, frame, order: int) -> float:
    """Return the flux-weighted mean slit coordinate of an object-only frame in the middle rows of an order."""
    isOrder = (scene.orderMap == order)[MIDDLE_ROWS]
    flux = np.where(isOrder, frame.data[MIDDLE_ROWS], 0.0)
    return float(np.sum(np.nan_to_num(scene.slit[MIDDLE_ROWS]) * flux) / np.sum(flux))


def test_a_positive_shift_moves_the_object_toward_lower_slit_coordinates_by_that_shift(scene, two_d_map) -> None:
    # ARRANGE: THE OBJECT SITS AT 1.0 ARCSEC IN ORDER 1, 0.6 ARCSEC ABOVE THE REFERENCE FRAME
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0, hasSky=False)

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 0.6, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT
    assert _centroid(scene, aligned, 1) == pytest.approx(0.4, abs=0.005)
    assert _centroid(scene, aligned, 2) == pytest.approx(1.0, abs=0.005)


def test_each_order_moves_by_its_own_shift(scene, two_d_map) -> None:
    # ARRANGE: THE SAME OBJECT POSITION IN EVERY ORDER, SHIFTED BY A DIFFERENT AMOUNT IN EACH
    frame = synthetic_nod_frame(scene, objectSlit=-2.0, mjd=60000.0, hasSky=False)
    shifts = {1: 0.6, 2: -0.4, 3: 0.0, 4: 1.1}

    # ACT
    aligned = shift_frame_along_slit(frame, shifts, two_d_map)

    # ASSERT
    for order, shift in shifts.items():
        assert _centroid(scene, aligned, order) == pytest.approx(-2.0 - shift, abs=0.005)


def test_the_sky_does_not_move_when_a_frame_is_shifted(scene, two_d_map) -> None:
    # ARRANGE: A SKY-ONLY FRAME. ITS LINES FOLLOW LINES OF CONSTANT WAVELENGTH, WHICH ARE TILTED ON THE DETECTOR
    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0, objectPeak=0.0)

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 0.6, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT: THE SHIFT FOLLOWS CONSTANT WAVELENGTH, SO ONLY THE INTERPOLATION BETWEEN SAMPLES ALONG A LINE REMAINS
    interior = (scene.orderMap == 1) & (np.abs(scene.slit) < 3.0)
    interior[:50] = interior[-50:] = False
    peakToPeak = np.ptp(frame.data[interior])
    assert np.max(np.abs(aligned.data[interior] - frame.data[interior])) < 0.05 * peakToPeak


def test_a_bad_pixel_moves_with_the_data(scene, two_d_map) -> None:
    # ARRANGE: ONE BAD PIXEL AT ROW 300, COLUMN 40 OF ORDER 1, THAT LIES 3 PIXELS TO THE HIGH SIDE OF THE SHIFT
    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0, badPixels=((300, 40),))
    shift = 3 * PIXEL_ARCSEC

    # ACT
    aligned = shift_frame_along_slit(frame, {1: shift, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT: THE PIXEL THAT SAW THE BAD PIXEL NOW SITS 3 COLUMNS LOWER, AND NO OTHER PIXEL INSIDE THE SLIT IS BAD
    # THE FIRST AND LAST ROWS ARE LEFT OUT: THE TILT OF THE MAP PUTS THEIR SOURCES JUST BEYOND THE DETECTOR EDGE
    flagged = np.argwhere(aligned.mask[5:-5, :80] & (np.abs(scene.slit[5:-5, :80]) < 4.5)) + [5, 0]
    assert len(flagged) >= 1
    assert np.all(np.abs(flagged[:, 0] - 300) <= 2)
    assert np.all(np.abs(flagged[:, 1] - 37) <= 2)


def test_the_uncertainty_is_interpolated_as_a_variance(scene, two_d_map) -> None:
    # ARRANGE: UNCERTAINTIES OF 1 AND 3 ON ALTERNATE COLUMNS, SHIFTED BY HALF A PIXEL ALONG THE COLUMNS
    from astropy.nddata import StdDevUncertainty

    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0)
    alternating = np.where(np.arange(scene.shape[1]) % 2 == 0, 1.0, 3.0)[None, :] * np.ones(scene.shape)
    frame.uncertainty = StdDevUncertainty(alternating, unit=frame.unit)

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 0.5 * PIXEL_ARCSEC, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT: HALF-WAY BETWEEN VARIANCES 1 AND 9 IS 5, SO THE UNCERTAINTY IS SQRT(5), NOT 2
    middle = (scene.orderMap == 1) & (np.abs(scene.slit) < 3.0)
    middle[:50] = middle[-50:] = False
    assert aligned.uncertainty.array[middle] == pytest.approx(np.sqrt(5.0), rel=0.03)
    assert aligned.uncertainty.uncertainty_type == "std"


def test_pixels_whose_source_lies_beyond_the_slit_edge_are_flagged_bad(scene, two_d_map) -> None:
    # ARRANGE
    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0)

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 1.0, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT: AT A SLIT COORDINATE OF 5.0 THE SOURCE (SLIT 6.0) IS OFF THE LIT SLIT; AT 3.0 IT IS NOT
    rows = slice(200, 400)
    edge = (scene.orderMap[rows] == 1) & (scene.slit[rows] > 4.8) & (scene.slit[rows] < 5.2)
    centre = (scene.orderMap[rows] == 1) & (np.abs(scene.slit[rows] - 3.0) < 0.3)
    assert aligned.mask[rows][edge].all()
    assert not aligned.mask[rows][centre].any()


def test_pixels_outside_every_order_and_the_input_frame_are_left_alone(scene, two_d_map) -> None:
    # ARRANGE
    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0)
    gap = ~np.isfinite(scene.orderMap)
    frame.data[gap] = 7.0
    dataBefore, maskBefore, errsBefore = frame.data.copy(), frame.mask.copy(), frame.uncertainty.array.copy()

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 0.6, 2: -0.4, 3: 0.2, 4: 0.9}, two_d_map)

    # ASSERT
    assert (aligned.data[gap] == 7.0).all()
    assert aligned.data.dtype == np.float32
    assert aligned is not frame
    np.testing.assert_array_equal(frame.data, dataBefore)
    np.testing.assert_array_equal(frame.mask, maskBefore)
    np.testing.assert_array_equal(frame.uncertainty.array, errsBefore)
    assert aligned.header["MJD-OBS"] == frame.header["MJD-OBS"]
    assert aligned.header is not frame.header


def test_a_zero_shift_returns_the_same_values(scene, two_d_map) -> None:
    # ARRANGE
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0, badPixels=((300, 40),))

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ASSERT
    np.testing.assert_array_equal(aligned.data, frame.data.astype(np.float32))
    np.testing.assert_array_equal(aligned.mask, frame.mask)
    np.testing.assert_array_equal(aligned.uncertainty.array, frame.uncertainty.array.astype(np.float32))


def test_a_frame_that_does_not_match_the_map_shape_is_refused(scene, two_d_map) -> None:
    # ARRANGE
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0)
    frame.data = frame.data[:, :-1]

    # ACT / ASSERT
    with pytest.raises(ValueError, match="differs from the 2D map shape"):
        shift_frame_along_slit(frame, {1: 0.6}, two_d_map)


def test_an_uncertainty_that_is_not_a_standard_deviation_is_refused(scene, two_d_map) -> None:
    # ARRANGE
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0)
    frame.uncertainty = VarianceUncertainty(frame.uncertainty.array**2, unit=frame.unit**2)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="not std"):
        shift_frame_along_slit(frame, {1: 0.6}, two_d_map)


def test_pixels_whose_source_lies_beyond_the_detector_edge_are_flagged_bad() -> None:
    # ARRANGE: ONE ORDER FILLS THE WHOLE DETECTOR, SO THE EDGE PIXEL IS IN THE ORDER. THE SLIT RUNS ALONG THE COLUMNS
    rows, columns = np.mgrid[0:40, 0:30]
    edgeMap = build_two_d_map_images(500.0 + 0.02 * rows, PIXEL_ARCSEC * (columns - 14.5), np.ones(rows.shape))
    frame = CCDData(np.ones(rows.shape), unit="electron", mask=np.zeros(rows.shape, dtype=bool))

    # ACT: A SHIFT OF 5 PIXELS, SO THE LAST 5 COLUMNS TAKE THEIR VALUE FROM BEYOND THE EDGE
    aligned = shift_frame_along_slit(frame, {1: 5 * PIXEL_ARCSEC}, edgeMap)

    # ASSERT
    assert aligned.mask[:, -5:].all()
    assert not aligned.mask[:, :-6].any()


def test_a_nan_pixel_the_input_mask_does_not_flag_is_flagged_in_the_output(scene, two_d_map) -> None:
    # ARRANGE: A NAN IN THE DATA AND ANOTHER IN THE UNCERTAINTY, NEITHER FLAGGED IN THE MASK
    frame = synthetic_nod_frame(scene, objectSlit=0.0, mjd=60000.0)
    errs = frame.uncertainty.array.copy()
    frame.data[300, 40] = np.nan
    errs[200, 40] = np.nan
    frame.uncertainty = StdDevUncertainty(errs, unit=frame.unit)

    # ACT
    aligned = shift_frame_along_slit(frame, {1: 3 * PIXEL_ARCSEC}, two_d_map)

    # ASSERT: THE PIXELS THAT READ THE NAN (3 COLUMNS LOWER) ARE BAD, AND NO PIXEL KEEPS A NON-FINITE VALUE UNFLAGGED
    assert aligned.mask[300, 36:38].any()
    assert aligned.mask[200, 36:38].any()
    isNonFinite = ~np.isfinite(aligned.data) | ~np.isfinite(aligned.uncertainty.array)
    assert aligned.mask[isNonFinite].all()
