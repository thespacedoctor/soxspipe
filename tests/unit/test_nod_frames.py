"""The A/B split, the balance check and the per-order rule of the nod frames (DY-2194)."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.commonutils.nod_frames import aligns_each_order, check_ab_balance, split_ab_frames

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
