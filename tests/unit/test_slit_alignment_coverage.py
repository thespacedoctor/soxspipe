"""The slit range that every frame of a nod stack covers, for leaving the uncovered pixels out of it (DY-2290)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from soxspipe.commonutils.nod_shift_and_add import shift_and_add_plan
from soxspipe.commonutils.slit_alignment import (
    build_two_d_map_images,
    orders_without_coverage,
    shift_coverage_mask,
    two_d_map_images,
)
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


def _plan(shifts: dict[int, float], mapImages: two_d_map_images) -> shift_and_add_plan:
    return shift_and_add_plan(shifts, mapImages)


def _middle(isPixel: np.ndarray) -> np.ndarray:
    """Return ``isPixel`` limited to the middle rows."""
    limited = np.zeros(isPixel.shape, dtype=bool)
    limited[MIDDLE_ROWS] = isPixel[MIDDLE_ROWS]
    return limited


def test_an_unshifted_frame_leaves_nothing_uncovered(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ARRANGE
    plan = _plan(dict.fromkeys(scene.orders, 0.0), two_d_map)

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT
    assert uncovered.dtype == bool
    assert uncovered.shape == scene.shape
    assert not uncovered.any()


def test_a_positive_shift_leaves_the_high_end_of_the_slit_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE: A SHIFT IN ORDER 1 ONLY. A PIXEL AT SLIT COORDINATE X TAKES ITS VALUE FROM X + THE SHIFT
    plan = _plan({1: NOD_THROW, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)
    isOrder = scene.orderMap == 1
    edge = FOOTPRINT_HALF - NOD_THROW

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT: ABOVE THE EDGE THE PIXELS ARE UNCOVERED, BELOW IT COVERED, AND NO OTHER ORDER IS TOUCHED
    assert uncovered[_middle(isOrder & (scene.slit > edge + MARGIN_ARCSEC))].all()
    assert not uncovered[_middle(isOrder & (scene.slit < edge - MARGIN_ARCSEC))].any()
    assert not uncovered[scene.orderMap != 1].any()


def test_a_negative_shift_leaves_the_low_end_of_the_slit_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE
    plan = _plan({1: -SHORT_THROW, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)
    isOrder = scene.orderMap == 1
    edge = -FOOTPRINT_HALF + SHORT_THROW

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT
    assert uncovered[_middle(isOrder & (scene.slit < edge - MARGIN_ARCSEC))].all()
    assert not uncovered[_middle(isOrder & (scene.slit > edge + MARGIN_ARCSEC))].any()


def test_each_order_is_uncovered_by_its_own_shift(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ARRANGE
    plan = _plan({1: 1.0, 2: 4.0, 3: 0.0, 4: -2.0}, two_d_map)

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT: THE UNCOVERED SLIT RANGE GROWS WITH THE SIZE OF THE SHIFT, AND AN UNSHIFTED ORDER HAS NONE
    counts = {order: np.count_nonzero(uncovered & (scene.orderMap == order)) for order in scene.orders}
    assert counts[3] == 0
    assert 0 < counts[1] < counts[4] < counts[2]


def test_a_shift_longer_than_the_slit_leaves_the_whole_order_uncovered(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE
    plan = _plan({1: BEYOND_SLIT_LENGTH, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT: ORDER 1 HAS NO COMMON RANGE AND THE OTHER ORDERS ARE UNTOUCHED
    assert uncovered[scene.orderMap == 1].all()
    assert not uncovered[scene.orderMap != 1].any()
    assert orders_without_coverage(uncovered, two_d_map) == [1]


def test_orders_without_coverage_is_empty_when_every_order_keeps_some_slit(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE
    uncovered = shift_coverage_mask(_plan(dict.fromkeys(scene.orders, NOD_THROW), two_d_map))

    # ACT / ASSERT
    assert orders_without_coverage(uncovered, two_d_map) == []


def test_pixels_outside_every_order_are_never_uncovered(scene: alignment_map, two_d_map: two_d_map_images) -> None:
    # ARRANGE
    plan = _plan(dict.fromkeys(scene.orders, BEYOND_SLIT_LENGTH), two_d_map)

    # ACT
    uncovered = shift_coverage_mask(plan)

    # ASSERT
    assert not uncovered[np.isnan(scene.orderMap)].any()


def test_a_source_beyond_the_detector_edge_is_uncovered() -> None:
    # ARRANGE: ONE ORDER FILLS THE WHOLE DETECTOR, SO ONLY THE DETECTOR EDGE LIMITS THE COVERAGE. THE SLIT RUNS ALONG
    # THE COLUMNS
    rows, columns = np.mgrid[0:40, 0:30]
    edgeMap = build_two_d_map_images(500.0 + 0.02 * rows, PIXEL_ARCSEC * (columns - 14.5), np.ones(rows.shape))

    # ACT: A SHIFT OF 5 PIXELS, SO THE LAST 5 COLUMNS TAKE THEIR VALUE FROM BEYOND THE EDGE
    uncovered = shift_coverage_mask(_plan({1: 5 * PIXEL_ARCSEC}, edgeMap))

    # ASSERT
    assert uncovered[:, -5:].all()
    assert not uncovered[:, :-6].any()


def test_the_jitter_of_the_frames_within_a_nod_position_widens_the_uncovered_band(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE: THE B FRAMES ARE SHIFTED BY THE NOD THROW AND THE SECOND ONE IS ALSO ALIGNED BY +2.0 ARCSEC TO THE
    # FIRST, SO ITS TOTAL SHIFT IS 6.6. A FRAMES ARE ALIGNED BY 0 AND -2.0 ARCSEC
    plan = _plan({1: NOD_THROW, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)
    isOrder = scene.orderMap == 1
    withoutJitter = shift_coverage_mask(plan)
    plainEdge = FOOTPRINT_HALF - NOD_THROW
    highEdge = plainEdge - JITTER_ARCSEC
    lowEdge = -FOOTPRINT_HALF + JITTER_ARCSEC
    isInJitterBand = _middle(
        isOrder & (scene.slit > highEdge + MARGIN_ARCSEC) & (scene.slit < plainEdge - MARGIN_ARCSEC)
    )
    isInLowBand = _middle(isOrder & (scene.slit < lowEdge - MARGIN_ARCSEC))
    # WITHOUT THE JITTER THESE PIXELS ARE COVERED BY EVERY FRAME
    assert not withoutJitter[isInJitterBand].any()
    assert not withoutJitter[isInLowBand].any()

    # ACT
    uncovered = shift_coverage_mask(plan, [{}, {1: -JITTER_ARCSEC}], [{}, {1: JITTER_ARCSEC}])

    # ASSERT: THE WIDEST SHIFT OF ANY FRAME SETS EACH EDGE
    assert uncovered[isInJitterBand].all()
    assert uncovered[isInLowBand].all()
    assert not uncovered[
        _middle(isOrder & (scene.slit < highEdge - MARGIN_ARCSEC) & (scene.slit > lowEdge + MARGIN_ARCSEC))
    ].any()


def test_a_frame_that_jitters_back_toward_the_reference_adds_nothing(
    scene: alignment_map, two_d_map: two_d_map_images
) -> None:
    # ARRANGE: THE REFERENCE B FRAME, ALIGNED BY 0, ALREADY SETS THE EDGE AT THE NOD THROW
    plan = _plan({1: NOD_THROW, 2: 0.0, 3: 0.0, 4: 0.0}, two_d_map)

    # ACT
    uncovered = shift_coverage_mask(plan, [{}], [{}, {1: -JITTER_ARCSEC}])

    # ASSERT
    np.testing.assert_array_equal(uncovered, shift_coverage_mask(plan))
