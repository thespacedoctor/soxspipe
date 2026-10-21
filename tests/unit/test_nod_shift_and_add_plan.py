"""The slit shift that puts the B-A difference onto the A-B difference, from the two reference positions (DY-1674)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.nod_shift_and_add import (
    build_shift_and_add_plan,
    horne_aperture_pixels,
    plan_clears_aperture,
    shift_and_add_plan,
)
from soxspipe.commonutils.slit_alignment import (
    SHIFT_COLUMNS,
    WHOLE_FRAME,
    build_two_d_map_images,
    slit_alignment_result,
)
from tests.factories.slit_alignment import PIXEL_ARCSEC, synthetic_alignment_map

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


def test_the_plan_counts_the_orders_that_took_the_median_shift(log) -> None:
    # ARRANGE: ORDERS 3 AND 4 LACK A POSITION IN ONE GROUP
    resultA = _result(_per_order("a0", [1.5, 1.5, NAN, 1.5]))
    resultB = _result(_per_order("b0", [-1.5, -1.9, -1.5, NAN]))

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)

    # ASSERT
    assert plan.fallbackCount == 2


def test_a_plan_with_every_shift_measured_counts_no_fallback(log) -> None:
    # ARRANGE AND ACT
    plan = build_shift_and_add_plan(log, _result(_per_order("a0", [1.5] * 4)), _result(_per_order("b0", [-1.5] * 4)))

    # ASSERT
    assert plan.fallbackCount == 0


def test_the_shifts_of_a_plan_cannot_be_changed_through_the_plan_or_the_dict_it_was_built_from() -> None:
    # ARRANGE
    shifts = {1: -3.0, 2: -3.0}
    plan = shift_and_add_plan(shifts, MAP_IMAGES)

    # ACT
    shifts[1] = 99.0

    # ASSERT
    assert plan.orderShifts[1] == -3.0
    with pytest.raises(TypeError):
        plan.orderShifts[1] = 0.0


@pytest.mark.parametrize(
    ("binx", "biny", "dispersionAxis", "expected"),
    [(1, 1, "y", 15.0), (2, 4, "x", 7.5), (2, 4, "y", 3.75)],
)
def test_the_aperture_is_the_horne_slit_length_divided_by_the_binning_along_the_dispersion_axis(
    binx: int, biny: int, dispersionAxis: str, expected: float
) -> None:
    # ARRANGE AND ACT
    aperture = horne_aperture_pixels(15, dispersionAxis, binx, biny)

    # ASSERT
    assert aperture == pytest.approx(expected)


@pytest.fixture
def real_map(tmp_path):
    scene = synthetic_alignment_map(tmp_path)
    return build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)


# 15 PIXELS OF 0.2 ARCSEC ARE A 3 ARCSEC APERTURE
APERTURE_PIXELS = 15.0


@pytest.mark.parametrize(("throw", "expected"), [(-3.3, True), (3.3, True), (-2.7, False), (2.7, False)])
def test_a_throw_clears_the_aperture_only_when_it_is_at_least_the_aperture_length(
    log, real_map, throw: float, expected: bool
) -> None:
    # ARRANGE
    apertureArcsec = APERTURE_PIXELS * PIXEL_ARCSEC
    assert apertureArcsec == pytest.approx(3.0)
    plan = shift_and_add_plan(dict.fromkeys(real_map.orders, throw), real_map)

    # ACT
    clears = plan_clears_aperture(log, plan, APERTURE_PIXELS)

    # ASSERT
    assert clears is expected
    warnings = [message for level, message in log.messages if level == "warning"]
    assert bool(warnings) is (not expected)


def test_one_order_with_a_short_throw_stops_the_combination_and_the_warning_names_it(log, real_map) -> None:
    # ARRANGE
    shifts = dict.fromkeys(real_map.orders, -3.3)
    shifts[3] = -1.0
    plan = shift_and_add_plan(shifts, real_map)

    # ACT
    clears = plan_clears_aperture(log, plan, APERTURE_PIXELS)

    # ASSERT
    assert clears is False
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("orders 3," in message and "extracted separately" in message for message in warnings)
