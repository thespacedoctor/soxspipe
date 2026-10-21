"""Contracts for the not-flattened Horne extraction reusing the flattened pass's sky shifts (DY-1680)."""

from __future__ import annotations

import copy
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import FlatPassHandoff
from tests.unit._horne_vis_scene import (
    I_LINE_OFFSET_NM,
    R_LINE_OFFSET_NM,
    applied_shift,
    build_extractor,
    build_orders,
    night_header,
)

pytestmark = pytest.mark.unit

FLAT_SHIFTS_NM = {11: 0.0, 16: 0.21, 21: 0.18, 24: 0.07}


def _spectrum(orders=(11, 16, 21, 24)) -> pd.DataFrame:
    frames = [
        pd.DataFrame({"order": order, "wavelengthMean": [1000.0 + order, 1000.05 + order], "skyFlux": [1.0, 1.0]})
        for order in orders
    ]
    return pd.concat(frames, ignore_index=True)


def _extractor(log: object, tmp_path: Path, notFlattened: bool, flatPass: FlatPassHandoff | None = None) -> object:
    extractor = build_extractor(log, tmp_path, night_header())
    extractor.notFlattened = "_NOTFLAT" if notFlattened else ""
    extractor.flatPass = flatPass
    extractor.calls = []

    def tune_spy(extractedOrdersDF: pd.DataFrame, arm: str, byOrder: bool = True) -> pd.DataFrame:
        extractor.calls.append(("tune", arm))
        extractor.skyShiftsNm = dict(FLAT_SHIFTS_NM)
        return extractedOrdersDF

    def measure_spy(*args: object, **kwargs: object) -> tuple[float, list]:
        extractor.calls.append(("measure", args[1] if len(args) > 1 else None))
        return 0.0, []

    def gate_spy() -> None:
        extractor.calls.append(("gate", None))

    extractor.tune_wavelength_calibration_to_skylines = tune_spy
    extractor._measure_order_shift = measure_spy
    extractor._arc_attachment = gate_spy
    return extractor


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_the_not_flattened_pass_applies_the_flat_shifts_without_peak_matching(
    log: object, tmp_path: Path, arm: str
) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, notFlattened=True, flatPass=FlatPassHandoff(skyShiftsNm=FLAT_SHIFTS_NM))
    before = _spectrum()

    # ACT
    after = extractor._correct_wavelengths(before.copy(), arm)

    # ASSERT
    assert extractor.calls == []
    for order, shift in FLAT_SHIFTS_NM.items():
        assert applied_shift(before, after, order) == pytest.approx(shift, abs=1e-12)


def test_the_not_flattened_pass_records_no_sky_qcs(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, notFlattened=True, flatPass=FlatPassHandoff(skyShiftsNm=FLAT_SHIFTS_NM))

    # ACT
    extractor._correct_wavelengths(_spectrum(), "NIR")

    # ASSERT
    assert extractor.qc.empty or not [name for name in extractor.qc["qc_name"] if name.startswith("SKY")]


def test_the_flattened_pass_tunes_and_hands_its_shifts_on(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, notFlattened=False)

    # ACT
    extractor._correct_wavelengths(_spectrum(), "NIR")

    # ASSERT
    assert ("tune", "NIR") in extractor.calls
    assert extractor.flatPassHandoff == FlatPassHandoff(skyShiftsNm=FLAT_SHIFTS_NM)


def test_an_order_the_flat_pass_did_not_shift_is_left_unshifted(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, notFlattened=True, flatPass=FlatPassHandoff(skyShiftsNm={16: 0.21}))
    before = _spectrum(orders=(16, 17))

    # ACT
    after = extractor._correct_wavelengths(before.copy(), "NIR")

    # ASSERT
    assert applied_shift(before, after, 16) == pytest.approx(0.21, abs=1e-12)
    assert applied_shift(before, after, 17) == 0.0


def test_a_skipped_flat_tuning_leaves_the_not_flattened_pass_untuned(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE FLAT PASS SKIPPED TUNING FOR AN ATTACHED ARC, SO IT HANDS ON NO SHIFTS
    extractor = _extractor(log, tmp_path, notFlattened=True, flatPass=FlatPassHandoff(skyShiftsNm=None))
    before = _spectrum()

    # ACT
    after = extractor._correct_wavelengths(before.copy(), "NIR")

    # ASSERT
    assert extractor.calls == []
    pd.testing.assert_frame_equal(after, before)


def test_a_not_flattened_pass_without_a_hand_off_runs_no_peak_matching(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, notFlattened=True, flatPass=None)
    before = _spectrum()

    # ACT
    after = extractor._correct_wavelengths(before.copy(), "NIR")

    # ASSERT
    assert extractor.calls == []
    pd.testing.assert_frame_equal(after, before)
    assert any(level == "warning" for level, _ in log.messages)


def _forbid_peak_matching(extractor: object) -> None:
    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("the not-flattened pass must not run skyline peak matching")

    extractor.tune_wavelength_calibration_to_skylines = forbidden
    extractor._measure_order_shift = forbidden
    extractor._match_peaks_to_skylines = forbidden
    extractor._get_local_skylines_for_order = forbidden


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_real_flat_tuning_hands_the_not_flattened_pass_the_same_wavelengths(
    log: object, tmp_path: Path, arm: str
) -> None:
    # ARRANGE: THE FLAT PASS REALLY MATCHES SKY LINES; THE NOT-FLAT PASS SEES A DIFFERENT SKY BUT MAY NOT MATCH IT
    for name in ("flat", "notflat"):
        (tmp_path / name).mkdir()
    flat = build_extractor(log, tmp_path / "flat", night_header())
    flat.notFlattened = ""
    notFlat = build_extractor(log, tmp_path / "notflat", night_header())
    notFlat.notFlattened = "_NOTFLAT"
    _forbid_peak_matching(notFlat)
    flatScene = build_orders()
    notFlatScene = build_orders().assign(skyFlux=1.0)

    # ACT
    flatAfter = flat._correct_wavelengths(flatScene.copy(), arm)
    notFlat.flatPass = flat.flatPassHandoff
    notFlatAfter = notFlat._correct_wavelengths(notFlatScene.copy(), arm)

    # ASSERT: THE FLAT PASS MEASURED BLUEWARD SHIFTS (UNDOING THE RED LINE OFFSETS); THE NOT-FLAT PASS REUSED THEM
    assert applied_shift(flatScene, flatAfter, 4) == pytest.approx(-I_LINE_OFFSET_NM, abs=0.01)
    if arm == "NIR":
        assert applied_shift(flatScene, flatAfter, 1) == pytest.approx(-R_LINE_OFFSET_NM, abs=0.01)
    assert set(flat.flatPassHandoff.skyShiftsNm) == {1, 2, 3, 4}
    np.testing.assert_array_equal(notFlatAfter["wavelengthMean"].to_numpy(), flatAfter["wavelengthMean"].to_numpy())
    assert not [name for name in notFlat.qc.get("qc_name", []) if name.startswith("SKY")]
    assert [name for name in flat.qc["qc_name"] if name.startswith("SKY SHIFT O4")]


def test_the_hand_off_cannot_be_changed_after_the_flat_pass(log: object, tmp_path: Path) -> None:
    # ARRANGE
    shifts = {16: 0.21}
    handoff = FlatPassHandoff(skyShiftsNm=shifts)

    # ACT
    shifts[16] = 9.0

    # ASSERT
    assert handoff.skyShiftsNm[16] == 0.21
    with pytest.raises(TypeError):
        handoff.skyShiftsNm[16] = 1.0


def test_the_hand_off_survives_pickling_and_deep_copy() -> None:
    # ARRANGE
    handoff = FlatPassHandoff(skyShiftsNm={16: 0.21, "all": 0.0})

    # ACT
    restored = pickle.loads(pickle.dumps(handoff))  # noqa: S301 - ROUND-TRIPS AN OBJECT THIS TEST BUILT
    copied = copy.deepcopy(FlatPassHandoff(skyShiftsNm=None))

    # ASSERT
    assert restored == handoff
    assert copied == FlatPassHandoff(skyShiftsNm=None)


def test_a_repeated_flat_extraction_does_not_hand_on_the_earlier_shifts(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE FIRST TUNING APPLIES SHIFTS; THE SECOND (A SPY THAT RECORDS NOTHING) APPLIES NONE
    extractor = _extractor(log, tmp_path, notFlattened=False)
    extractor._correct_wavelengths(_spectrum(), "NIR")
    extractor.tune_wavelength_calibration_to_skylines = lambda orders, arm, byOrder=True: orders

    # ACT
    extractor._correct_wavelengths(_spectrum(), "NIR")

    # ASSERT
    assert extractor.flatPassHandoff == FlatPassHandoff(skyShiftsNm=None)
