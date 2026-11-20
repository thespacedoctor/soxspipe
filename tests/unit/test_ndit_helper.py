"""The NDIT of a frame is read from its header through the keyword lookup (DY-2317).

A NIR raw file is the mean of its `ESO DET NDIT` integrations. The reader returns that count as an integer of at
least 1: a missing card is 1 with no message (VIS and UVB frames carry none), and an invalid card is 1 with a warning.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.ndit import read_ndit

pytestmark = pytest.mark.unit

NDIT_CARD = "HIERARCH ESO DET NDIT"


def _kw(log: Any, instrument: str) -> Any:
    return keyword_lookup(log=log, settings={"instrument": instrument}).get


def _header(value: object | None = None) -> fits.Header:
    header = fits.Header()
    if value is not None:
        header[NDIT_CARD] = value
    return header


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
@pytest.mark.parametrize("ndit", [1, 2, 4, 600])
def test_returns_the_ndit_card_value_for_both_instruments(log: Any, instrument: str, ndit: int) -> None:
    # ARRANGE
    kw = _kw(log, instrument)

    # ACT
    result = read_ndit(_header(ndit), kw, log)

    # ASSERT
    assert result == ndit
    assert isinstance(result, int)
    assert not [message for level, message in log.messages if level == "warning"]


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
def test_returns_one_without_a_message_when_the_card_is_missing(log: Any, instrument: str) -> None:
    # ARRANGE
    kw = _kw(log, instrument)

    # ACT
    result = read_ndit(_header(), kw, log)

    # ASSERT
    assert result == 1
    assert not [message for level, message in log.messages if level == "warning"]


@pytest.mark.parametrize("value", [0, -2, 2.5, float("nan"), float("inf"), "two", "", True])
def test_returns_one_with_a_warning_when_the_card_is_invalid(log: Any, value: object) -> None:
    # ARRANGE
    kw = _kw(log, "soxs")

    # ACT
    result = read_ndit({kw("DET_NDIT"): value}, kw, log)

    # ASSERT
    assert result == 1
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "ESO DET NDIT" in warnings[0]


def test_accepts_a_whole_number_float_and_a_numpy_integer(log: Any) -> None:
    # ARRANGE
    kw = _kw(log, "soxs")

    # ACT
    fromFloat = read_ndit(_header(2.0), kw, log)
    fromNumpy = read_ndit({kw("DET_NDIT"): np.int64(4)}, kw, log)

    # ASSERT
    assert fromFloat == 2
    assert fromNumpy == 4
    assert isinstance(fromFloat, int)
    assert isinstance(fromNumpy, int)
