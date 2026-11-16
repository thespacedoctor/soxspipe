"""The A/B split, the balance check, the per-order rule and the nominal throw of the nod frames (DY-2194)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy.nddata import CCDData

from soxspipe.commonutils.nod_frames import (
    aligns_each_order,
    check_ab_balance,
    nominal_throw_shift_plan,
    split_ab_frames,
)
from soxspipe.commonutils.slit_alignment import build_two_d_map_images, shift_frame_along_slit
from tests.factories.slit_alignment import FOOTPRINT_HALF, PIXEL_ARCSEC, synthetic_alignment_map, synthetic_nod_frame

pytestmark = pytest.mark.unit

OFFSET_KEYWORD = "ESO SEQ CUMOFF Y"


class FakeFrame:
    """A frame that carries only a header."""

    def __init__(self, offset: float) -> None:
        self.header = {OFFSET_KEYWORD: offset}


def _keyword(name: str) -> str:
    return {"NOD_CUMULATIVE_OFFSETY": OFFSET_KEYWORD}[name]


def test_split_ab_frames_puts_positive_offsets_in_a_and_the_rest_in_b() -> None:
    # ARRANGE
    frames = [FakeFrame(2.1), FakeFrame(-2.9), FakeFrame(2.1), FakeFrame(-2.9), FakeFrame(0.0)]
    names = ["a1", "b1", "a2", "b2", "zero"]

    # ACT
    allA, allB, offsetsA, offsetsB, namesA, namesB = split_ab_frames(frames, names, _keyword)

    # ASSERT
    assert allA == [frames[0], frames[2]]
    assert allB == [frames[1], frames[3], frames[4]]
    assert offsetsA == [2.1, 2.1]
    assert offsetsB == [-2.9, -2.9, 0.0]
    assert namesA == ["a1", "a2"]
    assert namesB == ["b1", "b2", "zero"]


def test_check_ab_balance_returns_the_distinct_a_offsets(log: Any) -> None:
    # ACT
    uniqueOffsets = check_ab_balance(log, [2.1, 2.1, 1.0], [-2.9, -2.9, -4.0])

    # ASSERT
    assert sorted(uniqueOffsets) == [1.0, 2.1]
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# PROCESSING 3 AB NODDING CYCLES WITH 2 UNIQUE PAIRS OF OFFSET LOCATIONS"]


def test_check_ab_balance_announces_a_single_pair_without_a_plural(log: Any) -> None:
    # ACT
    check_ab_balance(log, [2.1], [-2.9])

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# PROCESSING 1 AB NODDING CYCLES WITH 1 UNIQUE PAIR OF OFFSET LOCATIONS"]


def test_check_ab_balance_refuses_unequal_numbers_of_a_and_b_frames(log: Any) -> None:
    # ACT / ASSERT
    with pytest.raises(Exception, match="Found 2 A frames and 1 B frames"):
        check_ab_balance(log, [2.1, 2.1], [-2.9])

    assert any(level == "error" for level, _ in log.messages)


def test_check_ab_balance_refuses_a_set_with_no_a_frames(log: Any) -> None:
    # ACT / ASSERT
    with pytest.raises(Exception, match="Did not find frames with a positive offset"):
        check_ab_balance(log, [], [])


@pytest.mark.parametrize(
    ("instrument", "arm", "expected"),
    [
        ("SOXS", "VIS", True),
        ("soxs", "vis", True),
        ("SOXS", "NIR", False),
        ("SOXS", "UVB", False),
        ("XSHOOTER", "VIS", False),
    ],
)
def test_only_soxs_vis_aligns_each_order(instrument: str, arm: str, expected: bool) -> None:
    assert aligns_each_order(instrument, arm) is expected


# ============================================================================
# THE NOMINAL THROW (DY-2194): SHIFT B ONTO A BY THE COMMANDED THROW WHEN NO OBJECT POSITION WAS MEASURED
# ============================================================================

MJD_KEYWORD = "MJD-OBS"
# THE SLIT POSITION OF THE OBJECT IS MINUS THE CUMULATIVE OFFSET, PLUS A CONSTANT. MEASURED ON THE REAL VIS NOD SOF
# 20251212T045200 (A AT +2.1, B AT -2.9: POSITIVE TRACE OF A-B AT SLIT -2.4, NEGATIVE AT +1.9) AND ON NIR NOD STANDARD
# 20260113T012905 (A AT -1.38, B AT +3.38)
OFFSET_A = 2.1
OFFSET_B = -2.9
SLIT_A = -OFFSET_A
SLIT_B = -OFFSET_B


def _nominal_keyword(name: str) -> str:
    return {"NOD_CUMULATIVE_OFFSETY": OFFSET_KEYWORD, "MJDOBS": MJD_KEYWORD}[name]


class DatedFrame:
    """A frame with a header only: its observation date and cumulative offset."""

    def __init__(self, mjd: float, offset: float) -> None:
        self.header = {MJD_KEYWORD: mjd, OFFSET_KEYWORD: offset}


@pytest.fixture
def scene(tmp_path):
    return synthetic_alignment_map(tmp_path)


@pytest.fixture
def two_d_map(scene):
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


def test_the_nominal_throw_shifts_every_order_by_the_a_offset_minus_the_b_offset(log: Any, two_d_map: Any) -> None:
    # ARRANGE
    framesA = [DatedFrame(60000.0, OFFSET_A)]
    framesB = [DatedFrame(60000.1, OFFSET_B)]

    # ACT
    plan = nominal_throw_shift_plan(log, framesA, framesB, two_d_map, _nominal_keyword)

    # ASSERT
    assert set(plan.orderShifts) == set(two_d_map.orders)
    assert all(shift == pytest.approx(OFFSET_A - OFFSET_B) for shift in plan.orderShifts.values())
    assert plan.mapImages is two_d_map
    assert plan.fallbackCount == len(two_d_map.orders)


def test_the_nominal_throw_uses_the_first_frame_in_time_of_each_group(log: Any, two_d_map: Any) -> None:
    # ARRANGE: THE LATER FRAME OF EACH GROUP IS LISTED FIRST AND CARRIES A JITTERED OFFSET
    framesA = [DatedFrame(60000.4, 9.0), DatedFrame(60000.0, OFFSET_A)]
    framesB = [DatedFrame(60000.5, -9.0), DatedFrame(60000.1, OFFSET_B)]

    # ACT
    plan = nominal_throw_shift_plan(log, framesA, framesB, two_d_map, _nominal_keyword)

    # ASSERT
    assert list(plan.orderShifts.values()) == [pytest.approx(OFFSET_A - OFFSET_B)] * len(two_d_map.orders)


def test_the_nominal_throw_warns_that_it_was_used(log: Any, two_d_map: Any) -> None:
    # ACT
    nominal_throw_shift_plan(
        log, [DatedFrame(60000.0, OFFSET_A)], [DatedFrame(60000.1, OFFSET_B)], two_d_map, _nominal_keyword
    )

    # ASSERT
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "nominal" in warnings[0]
    assert "5.000" in warnings[0]


def _peak_slit(scene: Any, frame: CCDData, order: int) -> float:
    """Return the slit coordinate of the brightest median column profile of an order, over the middle rows."""
    middle = slice(200, 400)
    isOrder = (scene.orderMap == order)[middle]
    columns = np.flatnonzero(isOrder.all(axis=0))
    slit = np.median(np.where(isOrder, scene.slit[middle], np.nan)[:, columns], axis=0)
    counts = np.median(np.where(isOrder, frame.data[middle], np.nan)[:, columns], axis=0)
    return float(slit[np.nanargmax(counts)])


def test_the_nominal_throw_puts_the_b_object_on_the_a_object(log: Any, scene: Any, two_d_map: Any) -> None:
    """Pins the sign: the B frames move by A offset minus B offset, in the direction of increasing slit coordinate."""
    # ARRANGE
    frameA = synthetic_nod_frame(scene, objectSlit=SLIT_A, mjd=60000.0, offset=OFFSET_A, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=SLIT_B, mjd=60000.1, offset=OFFSET_B, hasSky=False)
    plan = nominal_throw_shift_plan(log, [frameA], [frameB], two_d_map, _nominal_keyword)

    # ACT
    shifted = shift_frame_along_slit(frameB, plan.orderShifts, plan.mapImages)

    # ASSERT
    for order in scene.orders:
        assert _peak_slit(scene, frameB, order) == pytest.approx(SLIT_B, abs=PIXEL_ARCSEC)
        assert _peak_slit(scene, shifted, order) == pytest.approx(SLIT_A, abs=PIXEL_ARCSEC)


def test_the_shift_masks_the_band_of_the_slit_whose_source_fell_off_it(log: Any, scene: Any, two_d_map: Any) -> None:
    # ARRANGE
    frameB = synthetic_nod_frame(scene, objectSlit=SLIT_B, mjd=60000.1, offset=OFFSET_B, hasSky=False)
    plan = nominal_throw_shift_plan(
        log, [DatedFrame(60000.0, OFFSET_A)], [DatedFrame(60000.1, OFFSET_B)], two_d_map, _nominal_keyword
    )

    # ACT
    shifted = shift_frame_along_slit(frameB, plan.orderShifts, plan.mapImages)

    # ASSERT: A TARGET PIXEL TAKES ITS VALUE FROM THE SLIT COORDINATE THAT IS THE THROW HIGHER
    order = scene.orders[0]
    isOrder = scene.orderMap == order
    throw = OFFSET_A - OFFSET_B
    middleRows = np.zeros(scene.shape, dtype=bool)
    middleRows[200:400] = True
    offSlit = isOrder & middleRows & (scene.slit + throw > FOOTPRINT_HALF + 0.4)
    onSlit = isOrder & middleRows & (scene.slit < 0.0) & (scene.slit > -3.0)
    assert offSlit.any() and onSlit.any()
    assert shifted.mask[offSlit].all()
    assert not shifted.mask[onSlit].any()
