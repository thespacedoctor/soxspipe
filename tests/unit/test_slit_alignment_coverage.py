"""The slit range that a shifted frame covers, for leaving the uncovered pixels out of a stack (DY-2290)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.slit_alignment import build_two_d_map_images, shift_coverage_mask
from tests.factories.slit_alignment import FOOTPRINT_HALF, PIXEL_ARCSEC, synthetic_alignment_map

pytestmark = pytest.mark.unit

# THE SLIT LENGTH OF THE SYNTHETIC ORDERS IS 2 * FOOTPRINT_HALF, SO A SHIFT BEYOND IT LEAVES NO COMMON RANGE
BEYOND_SLIT_LENGTH = 2 * FOOTPRINT_HALF + 1.0
# THE ROWS AWAY FROM THE ORDER ENDS, WHERE THE TILTED SLIT COORDINATE IS CLOSE TO ITS NOMINAL VALUE
MIDDLE_ROWS = slice(200, 400)
MARGIN_ARCSEC = 0.5


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def two_d_map(scene):
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def test_an_unshifted_frame_leaves_nothing_uncovered(scene, two_d_map) -> None:
    # ARRANGE
    shifts = dict.fromkeys(scene.orders, 0.0)

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT
    assert uncovered.dtype == bool
    assert uncovered.shape == scene.shape
    assert not uncovered.any()


def test_a_positive_shift_leaves_the_high_end_of_the_slit_uncovered(scene, two_d_map) -> None:
    # ARRANGE: A SHIFT OF 4.6 ARCSEC IN ORDER 1 ONLY. PIXEL AT SLIT COORDINATE X TAKES ITS VALUE FROM X + 4.6
    shifts = {1: 4.6, 2: 0.0, 3: 0.0, 4: 0.0}
    isOrder = scene.orderMap == 1
    edge = FOOTPRINT_HALF - 4.6

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT: ABOVE THE EDGE THE PIXELS ARE UNCOVERED, BELOW IT COVERED, AND NO OTHER ORDER IS TOUCHED
    assert uncovered[MIDDLE_ROWS][(isOrder & (scene.slit > edge + MARGIN_ARCSEC))[MIDDLE_ROWS]].all()
    assert not uncovered[MIDDLE_ROWS][(isOrder & (scene.slit < edge - MARGIN_ARCSEC))[MIDDLE_ROWS]].any()
    assert not uncovered[scene.orderMap != 1].any()


def test_a_negative_shift_leaves_the_low_end_of_the_slit_uncovered(scene, two_d_map) -> None:
    # ARRANGE
    shifts = {1: -3.0, 2: 0.0, 3: 0.0, 4: 0.0}
    isOrder = scene.orderMap == 1
    edge = -FOOTPRINT_HALF + 3.0

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT
    assert uncovered[MIDDLE_ROWS][(isOrder & (scene.slit < edge - MARGIN_ARCSEC))[MIDDLE_ROWS]].all()
    assert not uncovered[MIDDLE_ROWS][(isOrder & (scene.slit > edge + MARGIN_ARCSEC))[MIDDLE_ROWS]].any()


def test_each_order_is_uncovered_by_its_own_shift(scene, two_d_map) -> None:
    # ARRANGE
    shifts = {1: 1.0, 2: 4.0, 3: 0.0, 4: -2.0}

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT: THE UNCOVERED SLIT RANGE GROWS WITH THE SIZE OF THE SHIFT, AND AN UNSHIFTED ORDER HAS NONE
    counts = {order: np.count_nonzero(uncovered & (scene.orderMap == order)) for order in scene.orders}
    assert counts[3] == 0
    assert 0 < counts[1] < counts[4] < counts[2]


def test_a_shift_longer_than_the_slit_leaves_the_whole_order_uncovered(scene, two_d_map) -> None:
    # ARRANGE
    shifts = {1: BEYOND_SLIT_LENGTH, 2: 0.0, 3: 0.0, 4: 0.0}

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT: ORDER 1 HAS NO COMMON RANGE AND THE OTHER ORDERS ARE UNTOUCHED
    assert uncovered[scene.orderMap == 1].all()
    assert not uncovered[scene.orderMap != 1].any()


def test_pixels_outside_every_order_are_never_uncovered(scene, two_d_map) -> None:
    # ARRANGE
    shifts = dict.fromkeys(scene.orders, BEYOND_SLIT_LENGTH)

    # ACT
    uncovered = shift_coverage_mask(shifts, two_d_map)

    # ASSERT
    assert not uncovered[np.isnan(scene.orderMap)].any()


def test_a_source_beyond_the_detector_edge_is_uncovered() -> None:
    # ARRANGE: ONE ORDER FILLS THE WHOLE DETECTOR, SO ONLY THE DETECTOR EDGE LIMITS THE COVERAGE. THE SLIT RUNS ALONG
    # THE COLUMNS
    rows, columns = np.mgrid[0:40, 0:30]
    edgeMap = build_two_d_map_images(500.0 + 0.02 * rows, PIXEL_ARCSEC * (columns - 14.5), np.ones(rows.shape))

    # ACT: A SHIFT OF 5 PIXELS, SO THE LAST 5 COLUMNS TAKE THEIR VALUE FROM BEYOND THE EDGE
    uncovered = shift_coverage_mask({1: 5 * PIXEL_ARCSEC}, edgeMap)

    # ASSERT
    assert uncovered[:, -5:].all()
    assert not uncovered[:, :-6].any()
