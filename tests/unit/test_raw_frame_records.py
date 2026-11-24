"""`PRO REC1 RAWi` lists the raw frames a product was made from, and only those (DY-2356)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from tests.factories import instrument_header
from tests.factories.keyword_recipes import keyword_recipe

pytestmark = pytest.mark.unit

RAW_NAMES = ["SOXS.A.fits", "SOXS.B.fits", "SOXS.C.fits"]


def _table() -> pd.DataFrame:
    """Return the summary of three raw frames. A raw row is one with no PRO TYPE, which reaches pandas as NaN."""
    return pd.DataFrame(
        {
            "filename": RAW_NAMES,
            "tag": ["OBJECT_VIS"] * 3,
            "ESO PRO TYPE": [float("nan")] * 3,
        }
    )


def _recorded_names(frame: CCDData) -> list[str]:
    return [
        frame.header[f"ESO PRO REC1 RAW{index} NAME"]
        for index in range(1, 4)
        if f"ESO PRO REC1 RAW{index} NAME" in frame.header
    ]


def _frame() -> CCDData:
    return CCDData(np.zeros((2, 2)), unit=u.electron, meta=instrument_header())


def test_a_second_stamp_replaces_the_raw_records_of_the_first(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE FIRST STAMP LISTS EVERY RAW FRAME, THE SECOND ONLY ONE
    recipe = keyword_recipe(tmp_path, log, rawFrames=None)
    frame = _frame()
    recipe._stamp_raw_frame_records(frame, _table())

    # ACT
    recipe._stamp_raw_frame_records(frame, _table(), rawFrames=["SOXS.B.fits"])

    # ASSERT
    assert _recorded_names(frame) == ["SOXS.B.fits"]
    assert "ESO PRO REC1 RAW2 CATG" not in frame.header


def test_a_requested_name_that_no_raw_frame_has_is_reported_with_a_warning(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    recipe = keyword_recipe(tmp_path, log, rawFrames=None)
    frame = _frame()

    # ACT
    recipe._stamp_raw_frame_records(frame, _table(), rawFrames=["SOXS.A.fits", "ARC-9.fits"])

    # ASSERT
    assert _recorded_names(frame) == ["SOXS.A.fits"]
    warnings = [message for level, message in log.messages if level == "warning"]
    assert len(warnings) == 1
    assert "ARC-9.fits" in warnings[0] and "SOXS.A.fits" not in warnings[0]


def test_every_raw_frame_is_listed_when_no_names_are_requested(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    recipe = keyword_recipe(tmp_path, log, rawFrames=None)
    frame = _frame()

    # ACT
    recipe._stamp_raw_frame_records(frame, _table())

    # ASSERT
    assert _recorded_names(frame) == RAW_NAMES
    assert not [message for level, message in log.messages if level == "warning"]
