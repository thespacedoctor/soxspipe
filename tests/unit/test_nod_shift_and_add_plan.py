"""The slit shift that puts the B-A difference onto the A-B difference, from the two reference positions (DY-1674)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.nod_shift_and_add import build_shift_and_add_plan
from soxspipe.commonutils.slit_alignment import SHIFT_COLUMNS, WHOLE_FRAME, slit_alignment_result

pytestmark = pytest.mark.unit

MAP_IMAGES = SimpleNamespace(orders=(1, 2, 3, 4))
NAN = float("nan")


def _result(rows: list[tuple[str, float, int, float, bool]], mapImages=MAP_IMAGES) -> slit_alignment_result:
    """Return a result whose shift table holds ``(frame, mjd, order, position, isReference)`` rows."""
    shifts = pd.DataFrame(
        [
            {
                "frame": frame,
                "mjd": mjd,
                "order": order,
                "position_arcsec": position,
                "shift_arcsec": 0.0,
                "is_fallback": False,
                "is_reference": isReference,
            }
            for frame, mjd, order, position, isReference in rows
        ]
    )[SHIFT_COLUMNS]
    return slit_alignment_result(frames=(), shifts=shifts, fallbackCount=0, mapImages=mapImages)


def _per_order(frame: str, positions: list[float], *, isReference: bool = True, mjd: float = 60000.0):
    return [(frame, mjd, order, position, isReference) for order, position in enumerate(positions, start=1)]


def test_each_order_shift_is_the_b_reference_position_minus_the_a_reference_position(log) -> None:
    # ARRANGE
    resultA = _result(_per_order("a0", [1.5, 1.4, 1.6, 1.5]))
    resultB = _result(_per_order("b0", [-1.5, -1.7, -1.2, -1.5]))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan.orderShifts == pytest.approx({1: -3.0, 2: -3.1, 3: -2.8, 4: -3.0})
    assert plan.mapImages is MAP_IMAGES


def test_only_the_reference_frame_of_each_group_sets_the_shift(log) -> None:
    # ARRANGE: THE SECOND FRAME OF EACH GROUP SITS ELSEWHERE, AND IS NOT THE FRAME THE STACK IS ALIGNED TO
    resultA = _result(_per_order("a0", [1.5] * 4) + _per_order("a1", [1.9] * 4, isReference=False, mjd=60000.1))
    resultB = _result(_per_order("b0", [-1.5] * 4) + _per_order("b1", [-1.1] * 4, isReference=False, mjd=60000.1))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan.orderShifts == pytest.approx({1: -3.0, 2: -3.0, 3: -3.0, 4: -3.0})


def test_a_whole_frame_position_gives_every_order_the_same_shift(log) -> None:
    # ARRANGE
    resultA = _result([("a0", 60000.0, WHOLE_FRAME, 1.5, True)])
    resultB = _result([("b0", 60000.1, WHOLE_FRAME, -1.4, True)])

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan.orderShifts == pytest.approx({1: -2.9, 2: -2.9, 3: -2.9, 4: -2.9})


def test_an_order_with_no_position_in_either_group_takes_the_median_shift_and_warns(log) -> None:
    # ARRANGE: ORDER 3 HAS NO A POSITION AND ORDER 4 HAS NO B POSITION. THE MEASURED SHIFTS ARE -3.0 AND -3.4
    resultA = _result(_per_order("a0", [1.5, 1.5, NAN, 1.5]))
    resultB = _result(_per_order("b0", [-1.5, -1.9, -1.5, NAN]))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan.orderShifts == pytest.approx({1: -3.0, 2: -3.4, 3: -3.2, 4: -3.2})
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("orders 3, 4" in message for message in warnings)


def test_no_position_in_any_order_gives_no_plan_and_a_warning(log) -> None:
    # ARRANGE
    resultA = _result(_per_order("a0", [NAN] * 4))
    resultB = _result(_per_order("b0", [-1.5] * 4))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan is None
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("extracted separately" in message for message in warnings)


def test_a_group_that_was_not_measured_gives_no_plan_and_a_warning(log) -> None:
    # ARRANGE: A GROUP OF ONE FRAME THAT WAS NOT MEASURED HAS NO MAP AND NO POSITION
    resultA = _result([("a0", 60000.0, WHOLE_FRAME, NAN, True)], mapImages=None)
    resultB = _result([("b0", 60000.1, WHOLE_FRAME, NAN, True)], mapImages=None)

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan is None
    assert any(level == "warning" for level, _ in log.messages)


def test_the_shifts_are_plain_floats(log) -> None:
    # ARRANGE
    resultA = _result(_per_order("a0", [1.5] * 4))
    resultB = _result(_per_order("b0", [-1.5] * 4))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert all(type(shift) is float for shift in plan.orderShifts.values())
    assert not any(np.isnan(shift) for shift in plan.orderShifts.values())
