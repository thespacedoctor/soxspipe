"""The target position of a merged spectrum, decoded from the telescope target cards (DY-2357).

``ESO TEL TARG ALPHA`` is a packed HHMMSS.TTT and ``ESO TEL TARG DELTA`` a packed +-DDMMSS.TTT. Neither is in degrees.
The ``RA`` and ``DEC`` cards of a raw frame are the pointing, which is not the target.
"""

from __future__ import annotations

import math
from typing import Any

import pytest

from soxspipe.commonutils.on_target_frames import on_target_frame
from soxspipe.commonutils.target_position import (
    decode_packed_dec,
    decode_packed_ra,
    target_position_cards,
)

pytestmark = pytest.mark.unit

POSITION_COMMENT = "[deg] Spectroscopic target position (J2000)"


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _frame(name: str, ra: float | None, dec: float | None) -> on_target_frame:
    return on_target_frame(
        sourceName=name,
        name=name,
        arcfile=name,
        mjdObs=None,
        exposureSeconds=None,
        obsId=None,
        progId=None,
        tplStart=None,
        targetRa=ra,
        targetDec=dec,
    )


# DECODING THE RIGHT ASCENSION


def test_the_real_packed_right_ascension_decodes_to_degrees() -> None:
    # ACT / ASSERT: 00:30:19.8 IS 7.5825 DEG
    assert decode_packed_ra(3019.8) == pytest.approx(7.5825, abs=1e-9)


def test_a_right_ascension_of_many_hours_decodes_to_degrees() -> None:
    # ACT / ASSERT: 12:00:00.0 IS 180 DEG
    assert decode_packed_ra(120000.0) == pytest.approx(180.0, abs=1e-9)


def test_a_right_ascension_given_as_a_string_decodes() -> None:
    # ACT / ASSERT
    assert decode_packed_ra("3019.8") == pytest.approx(7.5825, abs=1e-9)


@pytest.mark.parametrize("value", [240000.0, 6000.0, 3060.0, -3019.8, None, "text", math.nan, math.inf, True])
def test_an_invalid_right_ascension_decodes_to_none(value: Any) -> None:
    # ACT / ASSERT
    assert decode_packed_ra(value) is None


# DECODING THE DECLINATION


def test_the_real_packed_declination_decodes_to_degrees() -> None:
    # ACT / ASSERT: -33:15:01.08
    assert decode_packed_dec(-331501.08) == pytest.approx(-33.2503, abs=1e-9)


def test_a_positive_declination_decodes_to_degrees() -> None:
    # ACT / ASSERT: +12:30:00.0
    assert decode_packed_dec(123000.0) == pytest.approx(12.5, abs=1e-9)


@pytest.mark.parametrize("value", [-530.0, "-000530.0"])
def test_a_declination_between_minus_one_and_zero_degrees_keeps_its_sign(value: Any) -> None:
    # ACT / ASSERT: -00:05:30.0
    assert decode_packed_dec(value) == pytest.approx(-0.0916666667, abs=1e-9)


def test_a_declination_of_zero_decodes_to_zero() -> None:
    # ACT / ASSERT
    assert decode_packed_dec(-0.0) == 0.0


def test_a_declination_at_the_pole_decodes() -> None:
    # ACT / ASSERT
    assert decode_packed_dec(900000.0) == 90.0


@pytest.mark.parametrize("value", [910000.0, -910000.0, 906000.0, 6000.0, 3060.0, None, "text", math.nan, -math.inf])
def test_an_invalid_declination_decodes_to_none(value: Any) -> None:
    # ACT / ASSERT
    assert decode_packed_dec(value) is None


# THE CARDS OF A MERGED SPECTRUM


def test_the_earliest_frame_with_a_target_position_gives_the_cards(log: Any) -> None:
    # ARRANGE: THE FIRST FRAME HAS NO POSITION
    entries = (_frame("a", None, None), _frame("b", 7.5825, -33.2503), _frame("c", 99.0, 10.0))

    # ACT
    cards = target_position_cards(entries, log)

    # ASSERT
    assert cards == [
        ("RA", 7.5825, POSITION_COMMENT),
        ("DEC", -33.2503, POSITION_COMMENT),
        ("RADESYS", "FK5", "Reference frame of RA and DEC"),
        ("EQUINOX", 2000.0, "[yr] Equinox of RA and DEC"),
    ]
    assert not _warnings(log)


def test_a_frame_with_only_one_coordinate_is_not_used(log: Any) -> None:
    # ARRANGE
    entries = (_frame("a", 7.5825, None), _frame("b", 8.0, -30.0))

    # ACT
    cards = target_position_cards(entries, log)

    # ASSERT
    assert {keyword: value for keyword, value, _ in cards}["RA"] == 8.0


def test_no_frame_with_a_target_position_gives_no_cards_and_one_warning(log: Any) -> None:
    # ARRANGE
    entries = (_frame("a", None, None), _frame("b", None, None))

    # ACT
    cards = target_position_cards(entries, log)

    # ASSERT
    assert cards == []
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "ESO TEL TARG ALPHA" in warnings[0] and "DELTA" in warnings[0] and "arm fails" in warnings[0]
