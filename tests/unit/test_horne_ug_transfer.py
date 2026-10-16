"""Contracts for the VIS u/g wavelength transfer inside the Horne-extraction skyline tuning (DY-1608)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy.io import fits

import soxspipe
import soxspipe.commonutils.ug_wavelength_transfer as ug_transfer
from soxspipe.commonutils.horne_extraction import horne_extraction
from tests.factories import dispersion_map_fits

pytestmark = pytest.mark.unit

PIXEL_NM = 0.05
# SYNTHETIC DISPERSION MAP: y = 30000 + (30 o^2 - 150 o + 140) * wavelength, SO dy/dwavelength IS +20 px/nm IN
# ORDERS 1 AND 4 (r, i: +0.05 nm/px) AND -40 px/nm IN ORDERS 2 AND 3 (u, g: -0.025 nm/px)
RI_DISPERSION = 0.05
UG_DISPERSION = -0.025
# SKY LINES APPEAR THIS FAR RED OF THEIR CATALOGUE WAVELENGTHS: +0.6 px IN r AND +0.8 px IN i
R_LINE_OFFSET_NM = 0.03
I_LINE_OFFSET_NM = 0.04
SCIENCE_ABSROT = -19.76611
KEYWORDS = {"ABSROT": "ESO ADA ABSROT END", "PAR_ANG_START": "ESO TEL PARANG START"}
CATALOGUE = 501.0 + 1.5 * np.arange(12)


def _dispersion_map(tmp_path: Path) -> str:
    coefficients = pd.DataFrame(
        {
            "axis": ["x", "y"],
            "order_deg": [2, 2],
            "wavelength_deg": [1, 1],
            "slit_deg": [0, 0],
            "c000": [100.0, 30000.0],
            "c010": [0.0, 140.0],
            "c100": [0.0, 0.0],
            "c110": [0.0, -150.0],
            "c200": [0.0, 0.0],
            "c210": [0.0, 30.0],
        }
    )
    return str(dispersion_map_fits(tmp_path / "dispersion_map.fits", coefficients=coefficients))


def _sky(wavelength: np.ndarray, offsetNm: float) -> np.ndarray:
    lines = sum(np.exp(-0.5 * ((wavelength - (c + offsetNm)) / 0.08) ** 2) for c in CATALOGUE)
    return 1.0 + 20.0 * lines


def _orders() -> pd.DataFrame:
    wavelength = np.arange(500.0, 520.0, PIXEL_NM)
    frames = []
    for order, offset in [(1, R_LINE_OFFSET_NM), (2, None), (3, None), (4, I_LINE_OFFSET_NM)]:
        sky = np.ones_like(wavelength) if offset is None else _sky(wavelength, offset)
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": sky})
        )
    return pd.concat(frames, ignore_index=True)


def _vis_parameters() -> dict:
    path = Path(soxspipe.__file__).parent / "resources" / "soxs_detector_parameters.yaml"
    return yaml.safe_load(path.read_text())["VIS"]


def _extractor(
    log: object, tmp_path: Path, scienceHeader: dict, detectorParams: dict | None = None
) -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.debug = False
    extractor.axisB = "y"
    extractor.dispersionMap = _dispersion_map(tmp_path)
    extractor.detectorParams = _vis_parameters() if detectorParams is None else detectorParams
    extractor.skySubtractedFrame = SimpleNamespace(header=fits.Header(scienceHeader))
    extractor.kw = KEYWORDS.__getitem__
    extractor.qc = pd.DataFrame()
    extractor.recipeName = "soxs-stare"
    extractor.dateObs = "2026-07-05T05:26:29"
    extractor._get_local_skylines_for_order = lambda *args: (CATALOGUE, CATALOGUE)
    return extractor


def _night_header() -> dict:
    return {KEYWORDS["ABSROT"]: SCIENCE_ABSROT, KEYWORDS["PAR_ANG_START"]: -85.867}


def _qc(extractor: horne_extraction, name: str) -> pd.Series:
    rows = extractor.qc.loc[extractor.qc["qc_name"] == name]
    assert len(rows) == 1, name
    return rows.iloc[0]


def _applied_shift(before: pd.DataFrame, after: pd.DataFrame, order: int) -> float:
    mask = before["order"] == order
    return float((after.loc[mask, "wavelengthMean"] - before.loc[mask, "wavelengthMean"]).mean())


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
    extractor = _extractor(log, tmp_path, _night_header())
    before = _orders()

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
    assert _applied_shift(before, after, 2) == pytest.approx(0.05)
    assert _applied_shift(before, after, 3) == pytest.approx(-0.025)
    assert _qc(extractor, "UG SHIFT SOURCE")["qc_value"] == "model"


def test_ug_stay_unshifted_with_a_logged_reason_when_the_science_angle_is_a_placeholder(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    spy = _PredictorSpy({2: 2.0, 3: -1.0})
    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", spy)
    placeholder = {KEYWORDS["ABSROT"]: 0.0, KEYWORDS["PAR_ANG_START"]: 0.0}
    extractor = _extractor(log, tmp_path, placeholder)
    before = _orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "VIS")

    # ASSERT
    assert spy.calls == []
    assert _applied_shift(before, after, 2) == 0.0
    assert _applied_shift(before, after, 3) == 0.0
    source = _qc(extractor, "UG SHIFT SOURCE")
    assert source["qc_value"] == "none"
    assert source["qc_comment"] == "no usable science ABSROT"
    assert ("info", "u/g orders left unshifted: no usable science ABSROT") in log.messages


def test_ug_stay_unshifted_when_too_few_i_lines_survive(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    spy = _PredictorSpy({2: 2.0, 3: -1.0})
    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", spy)
    extractor = _extractor(log, tmp_path, _night_header())
    extractor._get_local_skylines_for_order = lambda wmin, wmax, order, col: (
        (CATALOGUE[:5], CATALOGUE[:5]) if order == 4 else (CATALOGUE, CATALOGUE)
    )
    before = _orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "VIS")

    # ASSERT
    assert spy.calls == []
    assert _applied_shift(before, after, 2) == 0.0
    assert _qc(extractor, "UG SHIFT SOURCE")["qc_comment"] == "order 4: 5 sky lines (min 10)"


def test_sky_shift_qcs_are_written_in_nm_and_in_signed_detector_px(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, _night_header())

    # ACT
    extractor.tune_wavelength_calibration_to_skylines(_orders(), "VIS")

    # ASSERT
    shiftNm = _qc(extractor, "SKY SHIFT O1")
    assert shiftNm["qc_unit"] == "nm"
    assert shiftNm["qc_value"] == pytest.approx(-R_LINE_OFFSET_NM, abs=0.005)
    shiftPx = _qc(extractor, "SKY SHIFT O1 PX")
    assert shiftPx["qc_unit"] == "px"
    assert shiftPx["qc_value"] == pytest.approx(R_LINE_OFFSET_NM / RI_DISPERSION, abs=0.1)
    assert shiftPx["qc_order"] == 1
    assert shiftPx["to_header"]
    # THE MODELLED u SHIFT IS RECORDED IN BOTH UNITS, AND THE TWO AGREE THROUGH THE FALLING u DISPERSION
    uNm = _qc(extractor, "SKY SHIFT O2")["qc_value"]
    uPx = _qc(extractor, "SKY SHIFT O2 PX")["qc_value"]
    assert uNm != 0.0
    assert uNm == pytest.approx(-UG_DISPERSION * uPx, abs=0.001)


def test_arms_without_ug_coefficients_get_px_qcs_but_no_ug_source(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, _night_header(), detectorParams={"dispersion-axis": "y"})
    before = _orders()

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(before.copy(), "NIR")

    # ASSERT
    assert "UG SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert "SKY SHIFT O1 PX" in set(extractor.qc["qc_name"])
    assert _applied_shift(before, after, 1) == pytest.approx(-R_LINE_OFFSET_NM, abs=0.005)
