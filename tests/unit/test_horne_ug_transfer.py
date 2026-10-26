"""Contracts for the VIS u/g wavelength transfer inside the Horne-extraction skyline tuning (DY-1608)."""

from __future__ import annotations

from pathlib import Path

import pytest

import soxspipe.commonutils.ug_wavelength_transfer as ug_transfer
from tests.unit._horne_vis_scene import (
    CATALOGUE,
    I_LINE_OFFSET_NM,
    KEYWORDS,
    R_LINE_OFFSET_NM,
    RI_DISPERSION,
    SCIENCE_ABSROT,
    UG_DISPERSION,
    applied_shift,
    build_extractor,
    build_orders,
    night_header,
    qc_by_name,
)

pytestmark = pytest.mark.unit


class _PredictorSpy:
    def __init__(self, shifts: dict[int, float]) -> None:
        self.shifts = shifts
        self.calls: list[tuple] = []

    def __call__(self, *args: object) -> dict[int, float]:
        self.calls.append(args)
        return self.shifts


def test_ug_receive_the_model_shift_when_r_and_i_pass_the_guards(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    spy = _PredictorSpy({2: 2.0, 3: -1.0})
    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", spy)
    extractor = build_extractor(log, tmp_path, night_header())
    before = build_orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "VIS")

    # ASSERT
    assert len(spy.calls) == 1
    deltaR, deltaI, scienceAngle, calibrationAngle, coeffs = spy.calls[0]
    assert deltaR == pytest.approx(R_LINE_OFFSET_NM / RI_DISPERSION, abs=0.1)
    assert deltaI == pytest.approx(I_LINE_OFFSET_NM / RI_DISPERSION, abs=0.1)
    assert scienceAngle == pytest.approx(SCIENCE_ABSROT)
    # THE DISPERSION MAP CARRIES NO ABSROT, AS A DAYTIME ARC, SO THE PARK ANGLE IS USED
    assert calibrationAngle == 0.0
    assert coeffs[2] == {"m": -2.714, "cos": -4.624, "sin": -5.17}
    # A +2 px MOTION IN A FALLING ORDER READS 0.05 nm TOO BLUE; A -1 px MOTION READS 0.025 nm TOO RED
    assert applied_shift(before, after, 2) == pytest.approx(0.05)
    assert applied_shift(before, after, 3) == pytest.approx(-0.025)
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_value"] == "model"


def test_ug_stay_unshifted_with_a_logged_reason_when_the_science_angle_is_a_placeholder(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    spy = _PredictorSpy({2: 2.0, 3: -1.0})
    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", spy)
    placeholder = {KEYWORDS["ABSROT"]: 0.0, KEYWORDS["PAR_ANG_START"]: 0.0}
    extractor = build_extractor(log, tmp_path, placeholder)
    before = build_orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "VIS")

    # ASSERT
    assert spy.calls == []
    assert applied_shift(before, after, 2) == 0.0
    assert applied_shift(before, after, 3) == 0.0
    source = qc_by_name(extractor, "UG SHIFT SOURCE")
    assert source["qc_value"] == "none"
    assert source["qc_comment"] == "no usable science ABSROT"
    assert ("info", "u/g orders left unshifted: no usable science ABSROT") in log.messages


def test_ug_stay_unshifted_when_too_few_i_lines_survive(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    spy = _PredictorSpy({2: 2.0, 3: -1.0})
    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", spy)
    extractor = build_extractor(log, tmp_path, night_header())
    extractor._get_local_skylines_for_order = lambda wmin, wmax, order, col: (
        (CATALOGUE[:5], CATALOGUE[:5]) if order == 4 else (CATALOGUE, CATALOGUE)
    )
    before = build_orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "VIS")

    # ASSERT
    assert spy.calls == []
    assert applied_shift(before, after, 2) == 0.0
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_comment"] == "order 4: 5 sky lines (min 10)"


def test_sky_shift_qcs_are_written_in_nm_and_in_signed_detector_px(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    # ACT
    extractor.tune_wavelength_calibration_to_skylines(build_orders(), "VIS")

    # ASSERT
    shiftNm = qc_by_name(extractor, "SKY SHIFT O1")
    assert shiftNm["qc_unit"] == "nm"
    assert shiftNm["qc_value"] == pytest.approx(-R_LINE_OFFSET_NM, abs=0.005)
    shiftPx = qc_by_name(extractor, "SKY SHIFT O1 PX")
    assert shiftPx["qc_unit"] == "px"
    assert shiftPx["qc_value"] == pytest.approx(R_LINE_OFFSET_NM / RI_DISPERSION, abs=0.1)
    assert shiftPx["qc_order"] == 1
    assert shiftPx["to_header"]
    # THE MODELLED u SHIFT IS RECORDED IN BOTH UNITS, AND THE TWO AGREE THROUGH THE FALLING u DISPERSION
    uNm = qc_by_name(extractor, "SKY SHIFT O2")["qc_value"]
    uPx = qc_by_name(extractor, "SKY SHIFT O2 PX")["qc_value"]
    assert uNm != 0.0
    assert uNm == pytest.approx(-UG_DISPERSION * uPx, abs=0.001)


def test_arms_without_ug_coefficients_get_px_qcs_but_no_ug_source(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams={"dispersion-axis": "y"})
    before = build_orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "NIR")

    # ASSERT
    assert "UG SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert "SKY SHIFT O1 PX" in set(extractor.qc["qc_name"])
    assert applied_shift(before, after, 1) == pytest.approx(-R_LINE_OFFSET_NM, abs=0.005)
