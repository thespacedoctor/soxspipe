"""Contracts for the not-flattened Horne extraction reusing the flattened pass's object trace (DY-1650).

Both passes of a standard trace the same unflattened frame, so the not-flattened (response) pass must not run the
continuum detection a second time: it reads the trace the flattened pass wrote.
"""

from __future__ import annotations

import pickle

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import toolkit
from soxspipe.commonutils.horne_extraction import FlatPassHandoff, horne_extraction

pytestmark = pytest.mark.unit

FLAT_TRACE_PATH = "/products/STD_OBJTRACE.fits"
DETECTED_TRACE_PATH = "/products/STD_OBJTRACE_DETECTED.fits"


def _trace_table(path: str) -> pd.DataFrame:
    """A pixel table that names the trace file it was read from, so a test can tell which trace was used."""
    return pd.DataFrame({"order": [11, 12], "ycoord_centre": [100.0, 120.0], "source": [path, path]})


class _DetectorSpy:
    """Stands in for ``detect_continuum`` and records every construction."""

    constructions: list[dict] = []
    productPath: str | None = DETECTED_TRACE_PATH
    reads: list[str] = []

    def __init__(self, **kwargs: object) -> None:
        type(self).constructions.append(kwargs)
        self.qc = kwargs["qcTable"]
        self.products = kwargs["productsTable"]

    def get(self) -> tuple:
        path = type(self).productPath
        table = None if path is None else _trace_table(path)
        return path, self.qc, self.products, None, table, None


@pytest.fixture
def detector(monkeypatch: pytest.MonkeyPatch) -> type[_DetectorSpy]:
    monkeypatch.setattr(_DetectorSpy, "constructions", [])
    monkeypatch.setattr(_DetectorSpy, "productPath", DETECTED_TRACE_PATH)
    monkeypatch.setattr(commonutils, "detect_continuum", _DetectorSpy)
    reads: list[str] = []

    def unpack_spy(log: object, orderTablePath: str) -> tuple:
        reads.append(orderTablePath)
        return None, _trace_table(orderTablePath), None

    monkeypatch.setattr(toolkit, "unpack_order_table", unpack_spy)
    monkeypatch.setattr(_DetectorSpy, "reads", reads, raising=False)
    return _DetectorSpy


def _extractor(log: object, *, notFlattened: bool, flatPass: FlatPassHandoff | None = None) -> horne_extraction:
    """An extractor holding the state the trace step reads, without the frames the real ``__init__`` opens."""
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.settings = {}
    extractor.recipeSettings = {}
    extractor.unflattenedFrame = object()
    extractor.dispersionMap = "/products/DISP_MAP.fits"
    extractor.sofName = "STD"
    extractor.recipeName = "soxs-stare"
    extractor.qc = pd.DataFrame({"qc_name": ["EARLIER"]})
    extractor.products = pd.DataFrame()
    extractor.startNightDate = ""
    extractor.debug = False
    extractor.notFlattened = "_NOTFLAT" if notFlattened else ""
    extractor.flatPass = flatPass
    extractor.arm = "NIR"
    extractor.inst = "SOXS"
    extractor.orderNums = [11, 12]
    extractor.slitHalfLength = 30
    return extractor


def test_the_flattened_pass_detects_the_trace_once_and_hands_it_on(log: object, detector: type) -> None:
    # ARRANGE
    extractor = _extractor(log, notFlattened=False)

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is True
    assert len(detector.constructions) == 1
    assert detector.constructions[0]["traceFrame"] is extractor.unflattenedFrame
    assert detector.reads == [DETECTED_TRACE_PATH]
    assert extractor.flatPassHandoff.objectTracePath == DETECTED_TRACE_PATH


def test_the_not_flattened_pass_reuses_the_flat_trace_without_detecting(log: object, detector: type) -> None:
    # ARRANGE
    handoff = FlatPassHandoff(skyShiftsNm={11: 0.1}, objectTracePath=FLAT_TRACE_PATH)
    extractor = _extractor(log, notFlattened=True, flatPass=handoff)
    qcBefore = extractor.qc

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is True
    assert detector.constructions == []
    assert detector.reads == [FLAT_TRACE_PATH]
    assert extractor.orderPixelTable["source"].unique().tolist() == [FLAT_TRACE_PATH]
    # THE REUSED TRACE ADDS NO CONTINUUM QCS OF ITS OWN
    assert extractor.qc is qcBefore


def test_the_not_flattened_pass_finds_no_trace_when_the_flat_pass_found_none(log: object, detector: type) -> None:
    # ARRANGE
    handoff = FlatPassHandoff(skyShiftsNm=None, objectTracePath=None)
    extractor = _extractor(log, notFlattened=True, flatPass=handoff)

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is False
    assert detector.constructions == []
    assert detector.reads == []
    assert extractor.orderPixelTable is None


def test_the_flattened_pass_hands_on_that_it_found_no_trace(log: object, detector: type) -> None:
    # ARRANGE
    detector.productPath = None
    extractor = _extractor(log, notFlattened=False)

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is False
    assert len(detector.constructions) == 1
    assert extractor.flatPassHandoff == FlatPassHandoff(skyShiftsNm=None, objectTracePath=None)


def test_a_not_flattened_pass_without_a_hand_off_still_detects(log: object, detector: type) -> None:
    # ARRANGE
    extractor = _extractor(log, notFlattened=True, flatPass=None)

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is True
    assert len(detector.constructions) == 1
    assert detector.constructions[0]["notFlattened"] is True
    assert detector.reads == [DETECTED_TRACE_PATH]


def test_the_not_flattened_pass_falls_back_to_the_slit_centre_like_the_flat_pass(log: object, detector: type) -> None:
    # ARRANGE
    extractor = _extractor(log, notFlattened=True, flatPass=FlatPassHandoff(skyShiftsNm=None, objectTracePath=None))
    extractor.noTraceFallback = True

    # ACT
    hasTrace = extractor._find_object_trace(locationSetIndex=False)

    # ASSERT
    assert hasTrace is True
    assert detector.constructions == []
    assert extractor.traceFallback is True
    assert extractor.orderPixelTable["order"].tolist() == [11, 12]


def test_the_flattened_extraction_keeps_its_trace_in_the_tuned_hand_off(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, notFlattened=False)
    extractor.objectTracePath = FLAT_TRACE_PATH

    def tune(orders: pd.DataFrame, arm: str) -> pd.DataFrame:
        extractor.skyShiftsNm = {11: 0.2}
        return orders

    extractor._tune_or_skip = tune

    # ACT
    extractor._correct_wavelengths(pd.DataFrame({"order": [11]}), "NIR")

    # ASSERT
    assert extractor.flatPassHandoff == FlatPassHandoff(skyShiftsNm={11: 0.2}, objectTracePath=FLAT_TRACE_PATH)


def test_the_trace_survives_pickling_of_the_hand_off() -> None:
    # ARRANGE
    handoff = FlatPassHandoff(skyShiftsNm={16: 0.21}, objectTracePath=FLAT_TRACE_PATH)

    # ACT
    restored = pickle.loads(pickle.dumps(handoff))  # noqa: S301 - ROUND-TRIPS AN OBJECT THIS TEST BUILT

    # ASSERT
    assert restored == handoff
    assert restored.objectTracePath == FLAT_TRACE_PATH
