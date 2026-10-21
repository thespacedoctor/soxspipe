"""Combining the A-B difference with the shifted B-A difference of a nod pair (DY-1674)."""

from __future__ import annotations

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.nod_shift_and_add import combine_nod_differences
from soxspipe.commonutils.slit_alignment import build_two_d_map_images
from tests.factories.slit_alignment import PIXEL_ARCSEC, synthetic_alignment_map, synthetic_nod_frame

pytestmark = pytest.mark.unit

# ROWS AWAY FROM THE ENDS OF THE ORDERS, WHERE EVERY RESAMPLED PROFILE IS COMPLETE
MIDDLE_ROWS = slice(200, 400)
# THE SLIT POSITION OF THE A AND B NODS IN THE SYNTHETIC SCENE, IN ARCSEC
NOD_A = 1.5
NOD_B = -1.5


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def two_d_map(scene):
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def _differences(scene, positionA: float, positionB: float):
    """Return the A-B and B-A differences of two object-only frames, with the A header on A-B."""
    frameA = synthetic_nod_frame(scene, objectSlit=positionA, mjd=60000.0, offset=3.0, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=positionB, mjd=60000.1, offset=-3.0, hasSky=False)
    abDifference = frameA.subtract(frameB)
    abDifference.header = frameA.header
    baDifference = frameB.subtract(frameA)
    baDifference.header = frameB.header
    return abDifference, baDifference


def _profile(scene, frame, order: int) -> tuple[np.ndarray, np.ndarray]:
    """Return the slit coordinates and the median counts across the middle rows of an order, along the slit."""
    isOrder = (scene.orderMap == order)[MIDDLE_ROWS]
    slit = np.where(isOrder, scene.slit[MIDDLE_ROWS], np.nan)
    data = np.where(isOrder, frame.data[MIDDLE_ROWS], np.nan)
    columns = np.flatnonzero(isOrder.all(axis=0))
    return np.median(slit[:, columns], axis=0), np.median(data[:, columns], axis=0)


def _extreme_near(scene, frame, order: int, slitPosition: float, *, isPositive: bool = True) -> float:
    """Return the highest (or lowest) median count within 0.4 arcsec of a slit position of an order."""
    slit, counts = _profile(scene, frame, order)
    window = counts[np.abs(slit - slitPosition) <= 0.4]
    return float(window.max() if isPositive else window.min())


def test_the_positive_trace_keeps_the_amplitude_of_one_difference_at_the_a_position(scene, two_d_map) -> None:
    # ARRANGE: A-B HAS ITS POSITIVE TRACE AT NOD A, AND B-A AT NOD B. B-A MOVES BY THE B MINUS A POSITION
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    orderShifts = dict.fromkeys(scene.orders, NOD_B - NOD_A)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT: THE PEAK IS THAT OF ONE DIFFERENCE, LESS THE SMALL LOSS OF A BILINEAR RESAMPLING
    for order in scene.orders:
        oneDifference = _extreme_near(scene, abDifference, order, NOD_A)
        assert _extreme_near(scene, combined, order, NOD_A) == pytest.approx(oneDifference, rel=0.06)


def test_the_negative_traces_are_halved_where_they_were_and_where_the_shift_puts_them(scene, two_d_map) -> None:
    # ARRANGE
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    orderShifts = dict.fromkeys(scene.orders, NOD_B - NOD_A)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT: THE NEGATIVE OF A-B STAYS AT NOD B, THE NEGATIVE OF B-A LANDS TWO NOD STEPS FROM NOD A
    oneDifference = _extreme_near(scene, abDifference, 1, NOD_B, isPositive=False)
    assert _extreme_near(scene, combined, 1, NOD_B, isPositive=False) == pytest.approx(oneDifference / 2, rel=0.02)
    assert _extreme_near(scene, combined, 1, 2 * NOD_A - NOD_B, isPositive=False) == pytest.approx(
        oneDifference / 2, rel=0.1
    )


def test_a_shift_of_the_wrong_sign_leaves_the_positive_trace_off_the_a_position(scene, two_d_map) -> None:
    # ARRANGE: THE SHIFT IS THE A MINUS B POSITION, THE REVERSE OF WHAT THE B-A DIFFERENCE NEEDS
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    orderShifts = dict.fromkeys(scene.orders, NOD_A - NOD_B)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT: THE POSITIVE TRACE OF B-A HAS MOVED TO 2B - A, SO NOD A HOLDS ONLY HALF OF ONE DIFFERENCE
    oneDifference = _extreme_near(scene, abDifference, 1, NOD_A)
    assert _extreme_near(scene, combined, 1, NOD_A) == pytest.approx(oneDifference / 2, rel=0.02)


def test_swapping_the_nod_positions_puts_the_positive_trace_at_the_swapped_a_position(scene, two_d_map) -> None:
    # ARRANGE: THE OBJECT IS AT -1.5 ARCSEC AT A AND AT +1.5 ARCSEC AT B
    abDifference, baDifference = _differences(scene, NOD_B, NOD_A)
    orderShifts = dict.fromkeys(scene.orders, NOD_A - NOD_B)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT
    oneDifference = _extreme_near(scene, abDifference, 1, NOD_B)
    assert _extreme_near(scene, combined, 1, NOD_B) == pytest.approx(oneDifference, rel=0.06)


def test_each_order_moves_by_its_own_shift(scene, two_d_map) -> None:
    # ARRANGE: THE B POSITION DIFFERS FROM ORDER TO ORDER
    positionsB = {1: -1.5, 2: -2.0, 3: -1.0, 4: -2.5}
    frameA = synthetic_nod_frame(scene, objectSlit=NOD_A, mjd=60000.0, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=positionsB, mjd=60000.1, offset=-3.0, hasSky=False)
    abDifference = frameA.subtract(frameB)
    abDifference.header = frameA.header
    baDifference = frameB.subtract(frameA)
    baDifference.header = frameB.header
    orderShifts = {order: position - NOD_A for order, position in positionsB.items()}

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT
    for order in scene.orders:
        oneDifference = _extreme_near(scene, abDifference, order, NOD_A)
        assert _extreme_near(scene, combined, order, NOD_A) == pytest.approx(oneDifference, rel=0.06)


def test_the_variance_is_a_quarter_of_the_summed_variances_and_the_a_header_is_kept(scene, two_d_map) -> None:
    # ARRANGE: UNIFORM ERRORS OF 3 AND 4 ELECTRONS, SO THE INTERPOLATION CANNOT CHANGE THEM
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    abDifference.uncertainty = StdDevUncertainty(np.full(abDifference.shape, 3.0), unit=u.electron)
    baDifference.uncertainty = StdDevUncertainty(np.full(baDifference.shape, 4.0), unit=u.electron)
    orderShifts = dict.fromkeys(scene.orders, NOD_B - NOD_A)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT
    isOrder = scene.orderMap == 1
    assert combined.uncertainty.uncertainty_type == "std"
    assert np.median(combined.uncertainty.array[isOrder]) == pytest.approx(np.sqrt((9.0 + 16.0) / 4.0), rel=1e-4)
    assert combined.header["MJD-OBS"] == abDifference.header["MJD-OBS"]
    assert combined.header["ESO SEQ CUMOFF Y"] == 3.0


def test_a_pixel_flagged_in_either_difference_is_flagged_in_the_combination(scene, two_d_map) -> None:
    # ARRANGE: A-B FLAGS THE PIXEL AT (300, 20), AND B-A THE PIXEL AT (300, 40) OF THE SAME ORDER
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    abMask = np.zeros(abDifference.shape, dtype=bool)
    baMask = np.zeros(baDifference.shape, dtype=bool)
    abMask[300, 20] = True
    baMask[300, 40] = True
    abDifference.mask, baDifference.mask = abMask, baMask
    shift = NOD_B - NOD_A
    orderShifts = dict.fromkeys(scene.orders, shift)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)

    # ASSERT: THE B-A FLAG MOVES BY THE SHIFT, 15 COLUMNS TOWARDS HIGHER COLUMNS, WITHIN ONE ROW OF TILT
    assert combined.mask[300, 20]
    assert combined.mask[299:302, 40 - round(shift / PIXEL_ARCSEC) - 1 : 40 - round(shift / PIXEL_ARCSEC) + 2].any()
    assert not combined.mask[300, 30]


def test_a_difference_without_a_mask_gives_a_combination_with_a_mask_of_the_shifted_pixels(scene, two_d_map) -> None:
    # ARRANGE
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    abDifference.mask = None
    baDifference.mask = None

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, dict.fromkeys(scene.orders, -3.0), two_d_map)

    # ASSERT
    assert combined.mask is not None
    assert not combined.mask[300, 30]


def _snapshot(frame) -> dict:
    """Return copies of everything a frame holds, to compare after a call."""
    return {
        "data": frame.data.copy(),
        "mask": None if frame.mask is None else frame.mask.copy(),
        "uncertainty": frame.uncertainty.array.copy(),
        "uncertaintyUnit": frame.uncertainty.unit,
        "unit": frame.unit,
        "header": dict(frame.header),
    }


def _assert_unchanged(frame, before: dict) -> None:
    after = _snapshot(frame)
    np.testing.assert_array_equal(after["data"], before["data"])
    np.testing.assert_array_equal(after["mask"], before["mask"])
    np.testing.assert_array_equal(after["uncertainty"], before["uncertainty"])
    assert after["uncertaintyUnit"] == before["uncertaintyUnit"]
    assert after["unit"] == before["unit"]
    assert after["header"] == before["header"]


def test_the_inputs_are_not_changed(scene, two_d_map) -> None:
    # ARRANGE: BOTH DIFFERENCES CARRY A MASK, AN UNCERTAINTY, A UNIT AND A HEADER
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    abDifference.mask = np.zeros(abDifference.shape, dtype=bool)
    abDifference.mask[300, 20] = True
    baDifference.mask = np.zeros(baDifference.shape, dtype=bool)
    baDifference.mask[300, 40] = True
    abBefore, baBefore = _snapshot(abDifference), _snapshot(baDifference)

    # ACT
    combined = combine_nod_differences(abDifference, baDifference, dict.fromkeys(scene.orders, -3.0), two_d_map)
    combined.header["ADDED"] = 1
    combined.mask[0, 0] = not combined.mask[0, 0]

    # ASSERT
    _assert_unchanged(abDifference, abBefore)
    _assert_unchanged(baDifference, baBefore)


def test_a_difference_with_no_uncertainty_is_refused(scene, two_d_map) -> None:
    # ARRANGE
    abDifference, baDifference = _differences(scene, NOD_A, NOD_B)
    baDifference = CCDData(baDifference.data, unit=u.electron, meta=baDifference.meta)

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="uncertainty"):
        combine_nod_differences(abDifference, baDifference, {1: -3.0}, two_d_map)


SKY_LEVEL = 50.0
# THE OBJECT PROFILE HAS A SIGMA OF 0.3 ARCSEC, SO THIS HALF-WIDTH HOLDS ALL OF ITS FLUX
APERTURE_HALF_WIDTH = 0.8


def _sky_and_object_frames(scene):
    """Return an A frame with its object at NOD_A and a B frame with its object at NOD_B, on one flat sky."""
    frameA = synthetic_nod_frame(scene, objectSlit=NOD_A, mjd=60000.0, offset=3.0, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=NOD_B, mjd=60000.1, offset=-3.0, hasSky=False)
    frameA.data[:] = np.where(np.isfinite(scene.orderMap), frameA.data + SKY_LEVEL, 0.0)
    frameB.data[:] = np.where(np.isfinite(scene.orderMap), frameB.data + SKY_LEVEL, 0.0)
    return frameA, frameB


def _aperture_mean(scene, frame, order: int, slitPosition: float) -> float:
    """Return the mean count of the middle rows inside an aperture centred on a slit position of an order."""
    inAperture = (scene.orderMap == order) & (np.abs(scene.slit - slitPosition) <= APERTURE_HALF_WIDTH)
    return float(frame.data[MIDDLE_ROWS][inAperture[MIDDLE_ROWS]].mean())


def test_the_sky_frame_holds_no_object_at_the_a_position_where_the_science_object_lands(scene, two_d_map) -> None:
    # ARRANGE: THE SAME SHIFT MOVES THE B-A TRACE FROM NOD B TO NOD A AND THE A FRAME'S SKY AT NOD B ONTO NOD A
    frameA, frameB = _sky_and_object_frames(scene)
    orderShifts = dict.fromkeys(scene.orders, NOD_B - NOD_A)
    abDifference = frameA.subtract(frameB)
    abDifference.header = frameA.header
    baDifference = frameB.subtract(frameA)
    baDifference.header = frameB.header

    # ACT
    science = combine_nod_differences(abDifference, baDifference, orderShifts, two_d_map)
    sky = combine_nod_differences(frameB, frameA, orderShifts, two_d_map)
    swappedSky = combine_nod_differences(frameA, frameB, orderShifts, two_d_map)

    # ASSERT: THE SCIENCE OBJECT IS AT NOD A, THE SKY APERTURE AT NOD A IS SKY ONLY, AND SWAPPING THE FRAMES WOULD PUT
    # THE OBJECT OF THE A FRAME INTO IT
    for order in scene.orders:
        oneDifference = _extreme_near(scene, abDifference, order, NOD_A)
        assert _extreme_near(scene, science, order, NOD_A) == pytest.approx(oneDifference, rel=0.06)
        assert _aperture_mean(scene, sky, order, NOD_A) == pytest.approx(SKY_LEVEL, abs=0.5)
        assert _aperture_mean(scene, swappedSky, order, NOD_A) > SKY_LEVEL + 10
