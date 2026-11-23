"""The slit range that every frame of a nod stack covers, for leaving the uncovered pixels out of it (DY-2290)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from soxspipe.commonutils.slit_alignment import build_two_d_map_images, two_d_map_images
from soxspipe.commonutils.slit_coverage import alignment_coverage_mask, orders_without_coverage
from tests.factories.slit_alignment import FOOTPRINT_HALF, PIXEL_ARCSEC, alignment_map, synthetic_alignment_map

pytestmark = pytest.mark.unit

# THE SLIT LENGTH OF THE SYNTHETIC ORDERS IS 2 * FOOTPRINT_HALF, SO A SHIFT BEYOND IT LEAVES NO COMMON RANGE
BEYOND_SLIT_LENGTH = 3 * FOOTPRINT_HALF
# THE NOD THROW OF A REAL STANDARD, AND A SMALLER ONE
NOD_THROW = 4.6
SHORT_THROW = 3.0
# THE JITTER OF A FRAME WITHIN ITS NOD POSITION
JITTER_ARCSEC = 2.0
# THE ROWS AWAY FROM THE ORDER ENDS, WHERE THE TILTED SLIT COORDINATE IS CLOSE TO ITS NOMINAL VALUE
MIDDLE_ROWS = slice(200, 400)
# A PIXEL THIS FAR FROM A COVERAGE EDGE IS SAFELY ON ONE SIDE OF IT
MARGIN_ARCSEC = 0.5


@pytest.fixture
def scene(tmp_path: Path) -> alignment_map:
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def two_d_map(scene: alignment_map) -> two_d_map_images:
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def _middle(isPixel: np.ndarray) -> np.ndarray:
    """Return ``isPixel`` limited to the middle rows."""
    limited = np.zeros(isPixel.shape, dtype=bool)
    limited[MIDDLE_ROWS] = isPixel[MIDDLE_ROWS]
    return limited


def test_frames_with_no_alignment_shift_leave_nothing_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ACT
    uncovered = alignment_coverage_mask([{}, {}], two_d_map)

    # ASSERT
    assert uncovered.dtype == bool
    assert uncovered.shape == scene.shape
    assert not uncovered.any()


def test_a_position_with_no_frame_shifts_leaves_nothing_uncovered(two_d_map: two_d_map_images) -> None:
    # ACT / ASSERT
    assert not alignment_coverage_mask([], two_d_map).any()


def test_a_frame_aligned_upward_leaves_the_high_end_of_the_slit_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE: THE REFERENCE FRAME IS UNSHIFTED AND THE SECOND IS ALIGNED BY +SHORT_THROW IN ORDER 1
    isOrder = scene.orderMap == 1
    edge = FOOTPRINT_HALF - SHORT_THROW

    # ACT
    uncovered = alignment_coverage_mask([{}, {1: SHORT_THROW}], two_d_map)

    # ASSERT: ABOVE THE EDGE THE PIXELS ARE UNCOVERED, BELOW IT COVERED, AND NO OTHER ORDER IS TOUCHED
    assert uncovered[_middle(isOrder & (scene.slit > edge + MARGIN_ARCSEC))].all()
    assert not uncovered[_middle(isOrder & (scene.slit < edge - MARGIN_ARCSEC))].any()
    assert not uncovered[scene.orderMap != 1].any()


def test_a_frame_aligned_downward_leaves_the_low_end_of_the_slit_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE
    isOrder = scene.orderMap == 1
    edge = -FOOTPRINT_HALF + SHORT_THROW

    # ACT
    uncovered = alignment_coverage_mask([{}, {1: -SHORT_THROW}], two_d_map)

    # ASSERT
    assert uncovered[_middle(isOrder & (scene.slit < edge - MARGIN_ARCSEC))].all()
    assert not uncovered[_middle(isOrder & (scene.slit > edge + MARGIN_ARCSEC))].any()


def test_the_widest_shift_of_any_frame_sets_each_edge(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ARRANGE: THREE FRAMES OF ONE NOD POSITION, ALIGNED BY -2.0, 0 AND +2.0 ARCSEC
    isOrder = scene.orderMap == 1
    highEdge = FOOTPRINT_HALF - JITTER_ARCSEC
    lowEdge = -FOOTPRINT_HALF + JITTER_ARCSEC

    # ACT
    uncovered = alignment_coverage_mask([{1: -JITTER_ARCSEC}, {}, {1: JITTER_ARCSEC}], two_d_map)

    # ASSERT
    assert uncovered[_middle(isOrder & (scene.slit > highEdge + MARGIN_ARCSEC))].all()
    assert uncovered[_middle(isOrder & (scene.slit < lowEdge - MARGIN_ARCSEC))].all()
    assert not uncovered[
        _middle(isOrder & (scene.slit < highEdge - MARGIN_ARCSEC) & (scene.slit > lowEdge + MARGIN_ARCSEC))
    ].any()


def test_a_frame_between_the_widest_shifts_adds_nothing(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ARRANGE
    widest = alignment_coverage_mask([{1: -JITTER_ARCSEC}, {}, {1: JITTER_ARCSEC}], two_d_map)

    # ACT
    withMiddleFrame = alignment_coverage_mask(
        [{1: -JITTER_ARCSEC}, {1: 0.5}, {}, {1: -0.5}, {1: JITTER_ARCSEC}], two_d_map
    )

    # ASSERT
    np.testing.assert_array_equal(withMiddleFrame, widest)


def test_each_order_is_uncovered_by_its_own_alignment_shift(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ACT
    uncovered = alignment_coverage_mask([{}, {1: 1.0, 2: 4.0, 3: 0.0, 4: -2.0}], two_d_map)

    # ASSERT: THE UNCOVERED SLIT RANGE GROWS WITH THE SIZE OF THE SHIFT, AND AN UNSHIFTED ORDER HAS NONE
    counts = {order: np.count_nonzero(uncovered & (scene.orderMap == order)) for order in scene.orders}
    assert counts[3] == 0
    assert 0 < counts[1] < counts[4] < counts[2]


def test_an_alignment_shift_longer_than_the_slit_leaves_the_whole_order_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ACT
    uncovered = alignment_coverage_mask([{}, {1: BEYOND_SLIT_LENGTH}], two_d_map)

    # ASSERT
    assert uncovered[scene.orderMap == 1].all()
    assert not uncovered[scene.orderMap != 1].any()
    assert orders_without_coverage(uncovered, two_d_map) == [1]


def test_pixels_outside_every_order_are_never_uncovered_by_an_alignment_shift(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ACT
    uncovered = alignment_coverage_mask([{}, dict.fromkeys(scene.orders, BEYOND_SLIT_LENGTH)], two_d_map)

    # ASSERT
    assert not uncovered[np.isnan(scene.orderMap)].any()


def test_an_alignment_source_beyond_the_detector_edge_is_uncovered() -> None:
    # ARRANGE: ONE ORDER FILLS THE WHOLE DETECTOR, SO ONLY THE DETECTOR EDGE LIMITS THE COVERAGE
    rows, columns = np.mgrid[0:40, 0:30]
    edgeMap = build_two_d_map_images(500.0 + 0.02 * rows, PIXEL_ARCSEC * (columns - 14.5), np.ones(rows.shape))

    # ACT: A SHIFT OF 5 PIXELS, SO THE LAST 5 COLUMNS TAKE THEIR VALUE FROM BEYOND THE EDGE
    uncovered = alignment_coverage_mask([{}, {1: 5 * PIXEL_ARCSEC}], edgeMap)

    # ASSERT
    assert uncovered[:, -5:].all()
    assert not uncovered[:, :-6].any()
