"""Contracts for the SOXS VIS sky-line signal-to-noise gate in the Horne-extraction skyline tuning (DY-1663)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.ug_wavelength_transfer import MIN_R_SECOND_LINE_SNR, MIN_SKY_LINE_SNR
from tests.unit._horne_vis_scene import (
    CATALOGUE,
    I_LINE_OFFSET_NM,
    PIXEL_NM,
    R_LINE_OFFSET_NM,
    applied_shift,
    build_extractor,
    night_header,
    qc_by_name,
)

pytestmark = pytest.mark.unit

# THE SCENE SKY IS 1 e- BETWEEN LINES; WITH THE 3.8 e- VIS READ NOISE A LINE OF HEIGHT h HAS SNR h / SQRT(1 + 3.8^2)
PHOTON_NOISE = np.sqrt(1.0 + 3.8**2)
BRIGHT_HEIGHT = 20.0
# 12 LINES OF HEIGHT 2 COMBINE TO SNR 2 * SQRT(12) / 3.93 = 1.8, FAR BELOW THE GATE
FAINT_HEIGHT = 2.0
# THE i LINES SIT 0.04 nm FROM THEIR CATALOGUE WAVELENGTHS, SO THE NEAREST 0.05 nm SAMPLE IS 0.01 nm OFF EACH PEAK
SAMPLED_PEAK_FRACTION = np.exp(-0.5 * (0.01 / 0.08) ** 2)


def _expected_i_snr(height: float) -> float:
    return height * SAMPLED_PEAK_FRACTION * np.sqrt(len(CATALOGUE)) / PHOTON_NOISE


def _sky(wavelength: np.ndarray, centres: np.ndarray, height: float | np.ndarray) -> np.ndarray:
    heights = np.broadcast_to(height, np.shape(centres))
    return 1.0 + sum(h * np.exp(-0.5 * ((wavelength - c) / 0.08) ** 2) for c, h in zip(centres, heights, strict=True))


def _scene(
    rHeight: float | np.ndarray = BRIGHT_HEIGHT,
    iHeight: float = BRIGHT_HEIGHT,
    rCentres: np.ndarray = CATALOGUE + R_LINE_OFFSET_NM,
) -> pd.DataFrame:
    wavelength = np.arange(500.0, 520.0, PIXEL_NM)
    frames = []
    for order, skyFlux in [
        (1, _sky(wavelength, rCentres, rHeight)),
        (2, np.ones_like(wavelength)),
        (3, np.ones_like(wavelength)),
        (4, _sky(wavelength, CATALOGUE + I_LINE_OFFSET_NM, iHeight)),
    ]:
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": skyFlux})
        )
    return pd.concat(frames, ignore_index=True)


def _tune(extractor: horne_extraction, scene: pd.DataFrame) -> pd.DataFrame:
    return extractor.tune_wavelength_calibration_to_skylines(scene.copy(), "VIS")


def test_a_bright_i_sky_is_shifted_and_its_snr_recorded(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene()

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, after, 4) == pytest.approx(-I_LINE_OFFSET_NM, abs=0.005)
    snr = qc_by_name(extractor, "SKY LINE SNR O4")
    assert snr["qc_value"] == pytest.approx(_expected_i_snr(BRIGHT_HEIGHT), abs=0.05)
    assert snr["qc_order"] == 4


def test_a_faint_i_sky_with_the_same_lines_gets_no_shift(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(iHeight=FAINT_HEIGHT)
    expectedSnr = _expected_i_snr(FAINT_HEIGHT)

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, after, 4) == 0.0
    # THE QC IS ROUNDED TO ONE DECIMAL
    assert qc_by_name(extractor, "SKY LINE SNR O4")["qc_value"] == pytest.approx(expectedSnr, abs=0.05)
    assert (
        "info",
        f"order 4: sky-line SNR {expectedSnr:.1f} (min {MIN_SKY_LINE_SNR}); no shift applied",
    ) in log.messages


def test_a_faint_i_sky_sends_r_to_its_own_measurement_and_leaves_u_and_g_unshifted(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(iHeight=FAINT_HEIGHT)

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, after, 1) == pytest.approx(-R_LINE_OFFSET_NM, abs=0.005)
    rSource = qc_by_name(extractor, "R SHIFT SOURCE")
    assert rSource["qc_value"] == "measured"
    assert any("r shift not predicted from i: order 4: sky-line SNR" in message for _, message in log.messages)
    ugSource = qc_by_name(extractor, "UG SHIFT SOURCE")
    assert ugSource["qc_value"] == "none"
    assert ugSource["qc_comment"].startswith("order 4: sky-line SNR")
    assert applied_shift(scene, after, 2) == 0.0
    assert applied_shift(scene, after, 3) == 0.0


def test_a_faint_r_sky_measured_on_its_own_gets_no_shift(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(rHeight=FAINT_HEIGHT, iHeight=FAINT_HEIGHT)

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    assert applied_shift(scene, after, 1) == 0.0
    rSource = qc_by_name(extractor, "R SHIFT SOURCE")
    assert rSource["qc_value"] == "none"
    assert rSource["qc_comment"].startswith("order 1: sky-line SNR")


def test_r_measured_on_its_own_from_one_bright_and_one_noise_level_line_gets_no_shift(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    # A BRIGHT LINE (SNR 7.4) CARRIES THE COMBINED SNR PAST THE GATE, BUT THE SECOND LINE (SNR 5 / 3.93 = 1.27, AND
    # 1.23 AT THE NEAREST SAMPLE, 0.02 nm OFF ITS PEAK) IS AT NOISE LEVEL, AS ON 20260626T044348, WHERE IT PULLED THE
    # TWO-LINE MEDIAN 1.2 px OFF
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(rHeight=np.array([30.0, 5.0]), iHeight=FAINT_HEIGHT, rCentres=CATALOGUE[[2, 7]] + R_LINE_OFFSET_NM)

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    assert qc_by_name(extractor, "SKY LINE SNR O1")["qc_value"] > MIN_SKY_LINE_SNR
    assert applied_shift(scene, after, 1) == 0.0
    rSource = qc_by_name(extractor, "R SHIFT SOURCE")
    assert rSource["qc_value"] == "none"
    assert rSource["qc_comment"].startswith("order 1: 2nd sky line SNR 1.2")
    assert f"(min {MIN_R_SECOND_LINE_SNR})" in rSource["qc_comment"]


def test_a_faint_r_sky_falls_back_to_the_i_prediction(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(rHeight=FAINT_HEIGHT)

    # ACT
    _tune(extractor, scene)

    # ASSERT
    rSource = qc_by_name(extractor, "R SHIFT SOURCE")
    assert rSource["qc_value"] == "predicted"
    assert "sky-line SNR" in rSource["qc_comment"]


def test_an_i_order_whose_wavelengths_are_all_nan_gets_no_shift(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene()
    scene.loc[scene["order"] == 4, "wavelengthMean"] = np.nan

    # ACT
    after = _tune(extractor, scene)

    # ASSERT
    iRows = scene["order"] == 4
    assert after.loc[iRows, "wavelengthMean"].isna().all()
    assert ("info", "order 4: 0 sky lines kept (min 2); no shift applied") in log.messages


def test_the_gate_leaves_the_nir_arm_alone(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams={"dispersion-axis": "y"})
    scene = _scene(rHeight=FAINT_HEIGHT, iHeight=FAINT_HEIGHT)

    # ACT
    after = extractor.tune_wavelength_calibration_to_skylines(scene.copy(), "NIR")

    # ASSERT
    assert applied_shift(scene, after, 4) == pytest.approx(-I_LINE_OFFSET_NM, abs=0.005)
    assert not any(name.startswith("SKY LINE SNR") for name in extractor.qc["qc_name"])
