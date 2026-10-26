"""Contracts for predicting the VIS r sky shift from the i shift in the Horne-extraction skyline tuning (DY-1645)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import soxspipe.commonutils.ug_wavelength_transfer as ug_transfer
from soxspipe.commonutils.horne_extraction import R_SEEDED_TOLERANCE_NM, horne_extraction
from soxspipe.commonutils.ug_wavelength_transfer import MIN_MEASURED_SIGMA_PX, nm_shift_to_detector_px
from tests.unit._horne_vis_scene import (
    CATALOGUE,
    I_LINE_OFFSET_NM,
    KEYWORDS,
    PIXEL_NM,
    RI_DISPERSION,
    SCIENCE_ABSROT,
    applied_shift,
    build_extractor,
    night_header,
    qc_by_name,
    sky_from_lines,
    sky_spectrum,
    vis_parameters,
)

pytestmark = pytest.mark.unit

# THE i SKY LINES SIT 0.8 px RED OF THEIR CATALOGUE WAVELENGTHS, WHICH THE r-FROM-i MODEL TURNS INTO ABOUT 0.89 px FOR r
I_SHIFT_PX = I_LINE_OFFSET_NM / RI_DISPERSION
R_COEFFICIENTS = ug_transfer.load_r_transfer_coefficients(vis_parameters())
PREDICTED_R_PX = ug_transfer.predict_r_detector_shift(I_SHIFT_PX, SCIENCE_ABSROT, 0.0, R_COEFFICIENTS)
PREDICTED_R_NM = -RI_DISPERSION * PREDICTED_R_PX
FIRST_PASS_TOLERANCE_NM = 0.5
JUNK_PEAK_OFFSET_NM = 8 * RI_DISPERSION


def _scene(rLineCentres: np.ndarray, iLineOffsetNm: float = I_LINE_OFFSET_NM) -> pd.DataFrame:
    wavelength = np.arange(500.0, 520.0, PIXEL_NM)
    frames = []
    for order, skyFlux in [
        (1, sky_from_lines(wavelength, rLineCentres)),
        (2, np.ones_like(wavelength)),
        (3, np.ones_like(wavelength)),
        (4, sky_spectrum(wavelength, iLineOffsetNm)),
    ]:
        frames.append(
            pd.DataFrame({"order": order, "wavelengthMean": wavelength, "pixelScaleNm": PIXEL_NM, "skyFlux": skyFlux})
        )
    return pd.concat(frames, ignore_index=True)


class _LookupSpy:
    """Records the orders whose catalogue lines are requested and serves per-order catalogues."""

    def __init__(self, catalogues: dict[int, np.ndarray]) -> None:
        self.catalogues = catalogues
        self.orders: list[int] = []

    def __call__(self, wmin: float, wmax: float, order: int | str, col: str) -> tuple[np.ndarray, np.ndarray]:
        self.orders.append(order)
        return self.catalogues[order], self.catalogues[order]

    @property
    def measured_orders(self) -> list[int]:
        return list(dict.fromkeys(self.orders))


def _tune(
    extractor: horne_extraction,
    scene: pd.DataFrame,
    catalogues: dict[int | str, np.ndarray],
    arm: str = "VIS",
    byOrder: bool = True,
) -> tuple[pd.DataFrame, _LookupSpy, list[float]]:
    lookup = _LookupSpy({1: CATALOGUE, 2: CATALOGUE, 3: CATALOGUE, 4: CATALOGUE, "all": CATALOGUE} | catalogues)
    extractor._get_local_skylines_for_order = lookup
    tolerances: list[float] = []
    realMatch = extractor._match_peaks_to_skylines

    def record_tolerance(peakPositions: object, calibrationSkylines: object, tolerance: float = 15) -> tuple:
        tolerances.append(tolerance)
        return realMatch(peakPositions, calibrationSkylines, tolerance=tolerance)

    extractor._match_peaks_to_skylines = record_tolerance
    after = extractor.tune_wavelength_calibration_to_skylines(scene.copy(), arm, byOrder=byOrder)
    return after, lookup, tolerances


def test_i_is_measured_before_r(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + 0.03)

    # ACT
    _, lookup, _ = _tune(extractor, scene, {})

    # ASSERT
    assert lookup.measured_orders == [4, 1]


def test_r_shift_is_predicted_when_fewer_than_two_seeded_r_lines_survive_and_a_junk_peak_is_ignored(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    # ONE REAL r LINE AT THE PREDICTED POSITION AND A JUNK PEAK 8 px FROM ANOTHER CATALOGUE LINE, WHICH THE 0.5 nm FIRST
    # PASS WOULD MATCH
    extractor = build_extractor(log, tmp_path, night_header())
    rLines = np.array([CATALOGUE[3] - PREDICTED_R_NM, CATALOGUE[7] + JUNK_PEAK_OFFSET_NM])
    scene = _scene(rLines)
    before = scene.copy()

    # ACT
    after, _, tolerances = _tune(extractor, scene, {})

    # ASSERT
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] + [R_SEEDED_TOLERANCE_NM] * 3
    assert applied_shift(before, after, 1) == pytest.approx(PREDICTED_R_NM, abs=0.002)
    source = qc_by_name(extractor, "R SHIFT SOURCE")
    assert source["qc_value"] == "predicted"
    assert "1 sky line" in source["qc_comment"]


def _expected_combination(measuredPx: float, sigmaPx: float) -> float:
    """The inverse-variance mean of a measured and the predicted r shift, worked out from the definition."""
    weightMeasured = 1 / sigmaPx**2
    weightPredicted = 1 / R_COEFFICIENTS["rms_px"] ** 2
    return (weightMeasured * measuredPx + weightPredicted * PREDICTED_R_PX) / (weightMeasured + weightPredicted)


def test_identical_r_line_shifts_are_combined_with_the_floored_sigma(log: object, tmp_path: Path) -> None:
    # ARRANGE
    # FOUR LINES AGREE TO WITHIN A FRACTION OF A PIXEL, SO THE MEASURED SIGMA FALLS BELOW THE 0.05 px FLOOR. THE
    # MEASUREMENT IS 0.6 px BEYOND THE PREDICTION, INSIDE THE 1 px SEEDED WINDOW
    measuredPx = PREDICTED_R_PX + 0.6
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE[[2, 5, 8, 10]] + measuredPx * RI_DISPERSION)
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {})

    # ASSERT
    expectedPx = _expected_combination(measuredPx, MIN_MEASURED_SIGMA_PX)
    assert expectedPx - PREDICTED_R_PX == pytest.approx(0.6 * 400 / (400 + 1 / 0.364**2), abs=1e-9)
    assert nm_shift_to_detector_px(applied_shift(before, after, 1), RI_DISPERSION) == pytest.approx(
        expectedPx, abs=0.003
    )
    assert applied_shift(before, after, 1) == pytest.approx(-RI_DISPERSION * expectedPx, abs=0.0002)
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "combined"


def test_scattered_r_line_shifts_pull_the_combined_value_towards_the_prediction(log: object, tmp_path: Path) -> None:
    # ARRANGE
    # LINES ALTERNATE 0.5 px EITHER SIDE OF A MEAN 0.2 px BEYOND THE PREDICTION: SCATTER 0.5 px, SIGMA 0.5 / SQRT(4)
    meanPx = PREDICTED_R_PX + 0.2
    linePx = meanPx + np.array([-0.5, 0.5, -0.5, 0.5])
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE[[2, 5, 8, 10]] + linePx * RI_DISPERSION)
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {})

    # ASSERT
    expectedPx = _expected_combination(meanPx, 0.25)
    assert expectedPx - PREDICTED_R_PX == pytest.approx(0.2 * (1 / 0.25**2) / (1 / 0.25**2 + 1 / 0.364**2), abs=1e-9)
    assert nm_shift_to_detector_px(applied_shift(before, after, 1), RI_DISPERSION) == pytest.approx(
        expectedPx, abs=0.01
    )
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "combined"


def test_u_and_g_consume_the_final_predicted_r_shift(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    received: list[tuple] = []

    def predictor_spy(*args: object) -> dict[int, float]:
        received.append(args)
        return {2: 0.0, 3: 0.0}

    monkeypatch.setattr(ug_transfer, "predict_ug_detector_shifts", predictor_spy)
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(np.array([CATALOGUE[3] - PREDICTED_R_NM]))

    # ACT
    _tune(extractor, scene, {})

    # ASSERT
    assert len(received) == 1
    deltaR, deltaI = received[0][:2]
    assert deltaR == pytest.approx(PREDICTED_R_PX, abs=0.1)
    assert deltaI == pytest.approx(I_SHIFT_PX, abs=0.1)
    assert qc_by_name(extractor, "UG SHIFT SOURCE")["qc_value"] == "model"


def test_r_is_measured_independently_when_i_fails_its_guards(log: object, tmp_path: Path) -> None:
    # ARRANGE
    rOffsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + rOffsetNm)
    before = scene.copy()

    # ACT
    after, _, tolerances = _tune(extractor, scene, {4: CATALOGUE[:5]})

    # ASSERT
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] * 2
    assert applied_shift(before, after, 1) == pytest.approx(-rOffsetNm, abs=0.005)
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"


def test_r_is_measured_independently_when_the_science_angle_is_unusable(log: object, tmp_path: Path) -> None:
    # ARRANGE
    rOffsetNm = 0.03
    placeholder = {KEYWORDS["ABSROT"]: 0.0, KEYWORDS["PAR_ANG_START"]: 0.0}
    extractor = build_extractor(log, tmp_path, placeholder)
    scene = _scene(CATALOGUE + rOffsetNm)
    before = scene.copy()

    # ACT
    after, _, tolerances = _tune(extractor, scene, {})

    # ASSERT
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] * 2
    assert applied_shift(before, after, 1) == pytest.approx(-rOffsetNm, abs=0.005)
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"


def test_r_gets_no_shift_when_i_fails_and_fewer_than_two_r_lines_are_kept(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + 0.03)
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {4: CATALOGUE[:5], 1: CATALOGUE[:1]})

    # ASSERT
    assert applied_shift(before, after, 1) == 0.0
    source = qc_by_name(extractor, "R SHIFT SOURCE")
    assert source["qc_value"] == "none"
    assert "1 sky lines kept (min 2)" in source["qc_comment"]
    assert ("info", "order 1: 1 sky lines kept (min 2); no shift applied") in log.messages


def test_any_vis_order_with_fewer_than_two_kept_lines_gets_no_shift(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + 0.03)
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {4: CATALOGUE[:1]})

    # ASSERT
    assert applied_shift(before, after, 4) == 0.0
    assert ("info", "order 4: 1 sky lines kept (min 2); no shift applied") in log.messages


def test_the_nir_arm_still_applies_a_shift_from_a_single_line_and_records_no_r_source(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    offsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams={"dispersion-axis": "y"})
    scene = _scene(CATALOGUE + offsetNm)
    before = scene.copy()

    # ACT
    after, lookup, tolerances = _tune(extractor, scene, {1: CATALOGUE[:1]}, arm="NIR")

    # ASSERT
    assert lookup.measured_orders == [1, 2, 3, 4]
    assert tolerances[:3] == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1]
    assert applied_shift(before, after, 1) == pytest.approx(-offsetNm, abs=0.005)
    assert "R SHIFT SOURCE" not in set(extractor.qc["qc_name"])


def test_r_is_measured_independently_when_the_dispersion_is_unknown(log: object, tmp_path: Path) -> None:
    # ARRANGE
    rOffsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header())
    extractor._signed_order_dispersions = lambda *args: {}
    scene = _scene(CATALOGUE + rOffsetNm)
    before = scene.copy()

    # ACT
    after, _, tolerances = _tune(extractor, scene, {})

    # ASSERT
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] * 2
    assert applied_shift(before, after, 1) == pytest.approx(-rOffsetNm, abs=0.005)
    assert ("info", "r shift not predicted from i: order 4: no dispersion") in log.messages
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"


def test_an_order_without_usable_sky_data_gets_no_shift(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + 0.03)
    scene.loc[scene["order"] == 4, "skyFlux"] = np.nan
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {})

    # ASSERT
    assert applied_shift(before, after, 4) == 0.0
    assert ("info", "order 4: 0 sky lines kept (min 2); no shift applied") in log.messages


def test_global_tuning_shifts_every_order_together_and_looks_up_all_orders(log: object, tmp_path: Path) -> None:
    # ARRANGE
    offsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header(), detectorParams={"dispersion-axis": "y"})
    scene = _scene(CATALOGUE + offsetNm)
    scene = scene.loc[scene["order"] == 1].reset_index(drop=True)
    before = scene.copy()

    # ACT
    after, lookup, _ = _tune(extractor, scene, {}, arm="NIR", byOrder=False)

    # ASSERT
    assert set(lookup.orders) == {"all"}
    assert applied_shift(before, after, 1) == pytest.approx(-offsetNm, abs=0.005)


def test_xshooter_vis_keeps_the_old_per_order_behaviour(log: object, tmp_path: Path) -> None:
    # ARRANGE
    # THE SOXS VIS PARAMETERS CARRY THE r-FROM-i BLOCK, SO ONLY THE INSTRUMENT GATE KEEPS THE NEW BEHAVIOUR OFF
    offsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header(), inst="XSHOOTER")
    scene = _scene(CATALOGUE + offsetNm)
    before = scene.copy()

    # ACT
    after, lookup, tolerances = _tune(extractor, scene, {1: CATALOGUE[:1], 4: CATALOGUE[:1]})

    # ASSERT
    assert lookup.measured_orders == [1, 4]
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] * 2
    assert applied_shift(before, after, 1) == pytest.approx(-offsetNm, abs=0.005)
    assert applied_shift(before, after, 4) == pytest.approx(-I_LINE_OFFSET_NM, abs=0.005)
    assert "R SHIFT SOURCE" not in set(extractor.qc["qc_name"])
    assert not any("no shift applied" in message for _, message in log.messages)


def test_global_tuning_in_soxs_vis_applies_a_single_line_shift_with_the_old_tolerances(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    offsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header())
    scene = _scene(CATALOGUE + offsetNm)
    scene = scene.loc[scene["order"] == 1].reset_index(drop=True)
    before = scene.copy()

    # ACT
    after, lookup, tolerances = _tune(extractor, scene, {"all": CATALOGUE[:1]}, byOrder=False)

    # ASSERT
    assert set(lookup.orders) == {"all"}
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1]
    assert applied_shift(before, after, 1) == pytest.approx(-offsetNm, abs=0.005)
    assert extractor.qc.empty


def test_an_unreadable_dispersion_map_leaves_r_to_be_measured_alone(log: object, tmp_path: Path) -> None:
    # ARRANGE
    rOffsetNm = 0.03
    extractor = build_extractor(log, tmp_path, night_header())
    extractor.dispersionMap = str(tmp_path / "missing_dispersion_map.fits")
    scene = _scene(CATALOGUE + rOffsetNm)
    before = scene.copy()

    # ACT
    after, _, tolerances = _tune(extractor, scene, {})

    # ASSERT
    assert tolerances == [FIRST_PASS_TOLERANCE_NM, 0.1, 0.1] * 2
    assert applied_shift(before, after, 1) == pytest.approx(-rOffsetNm, abs=0.005)
    assert any(level == "info" and "dispersion" in message for level, message in log.messages)
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"


@pytest.mark.parametrize("mapContent", [b"not a fits file", None], ids=["corrupt-file", "primary-hdu-only"])
def test_a_corrupt_or_tableless_dispersion_map_leaves_r_to_be_measured_alone(
    log: object, tmp_path: Path, mapContent: bytes | None
) -> None:
    # ARRANGE
    from astropy.io import fits

    mapPath = tmp_path / "bad_dispersion_map.fits"
    if mapContent is None:
        fits.PrimaryHDU().writeto(mapPath)
    else:
        mapPath.write_bytes(mapContent)
    extractor = build_extractor(log, tmp_path, night_header())
    extractor.dispersionMap = str(mapPath)
    scene = _scene(CATALOGUE + 0.03)
    before = scene.copy()

    # ACT
    after, _, _ = _tune(extractor, scene, {})

    # ASSERT
    assert applied_shift(before, after, 1) == pytest.approx(-0.03, abs=0.005)
    assert qc_by_name(extractor, "R SHIFT SOURCE")["qc_value"] == "measured"


@pytest.mark.parametrize("error", [KeyError("axis"), TypeError("bad operand")], ids=["key-error", "type-error"])
def test_an_unexpected_error_in_the_dispersion_computation_propagates(
    log: object, tmp_path: Path, error: Exception
) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path, night_header())

    def fail(*args: object) -> dict:
        raise error

    extractor._signed_order_dispersions = fail

    # ACT / ASSERT
    with pytest.raises(type(error)):
        _tune(extractor, _scene(CATALOGUE + 0.03), {})
