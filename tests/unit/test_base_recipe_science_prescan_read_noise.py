"""The error plane of a UVB/VIS science frame uses the read noise of its own prescan (DY-2437).

The master bias read noise (`QC RAW RON`) is measured on the morning biases, and the night frames can have more.
A science recipe measures each raw frame's prescan before the frame is trimmed, builds that frame's uncertainty with
it, and writes the value to the product header as `QC SCI RON`. Calibration recipes and the NIR arm are unchanged.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.commonutils.read_noise import SOURCE_SCIENCE_PRESCAN, ReadNoise
from tests.factories.prepared_frames import (
    KEYWORDS,
    TEST_CRYOSTAT_TEMP,
    TEST_GAIN,
    TEST_NIR_DETECTOR_TEMP,
    TEST_RON,
    TEST_VIS_DETECTOR_TEMP,
)
from tests.factories.prepared_frames import prepare_frame_recipe as _recipe

pytestmark = pytest.mark.unit

ROWS = 200
PRESCAN_COLUMNS = 30
SCIENCE_COLUMNS = 40
SKY_ADU = 100.0
PRESCAN_SIGMA_ADU = 5.0
# THE PRESCAN IS READ AS COLUMNS 4 TO 28 OF EVERY ROW
PRESCAN_WINDOW = {"rows": {"start": 0, "end": ROWS}, "columns": {"start": 4, "end": 28}}
MASTER_BIAS_RON = TEST_RON
PRESCAN_RON = PRESCAN_SIGMA_ADU * TEST_GAIN
TOLERANCE = 0.04


def _raw_science_frame(tmpPath: Path, *, prescanSigma: float = PRESCAN_SIGMA_ADU, name: str = "raw.fits") -> str:
    """A raw VIS science frame in ADU: a flat science region beside a noisy prescan."""
    rng = np.random.default_rng(7)
    data = np.full((ROWS, PRESCAN_COLUMNS + SCIENCE_COLUMNS), SKY_ADU, dtype=np.float32)
    data[:, 4:28] = rng.normal(590.0, prescanSigma, (ROWS, 24)) if prescanSigma else 590.0
    frame = CCDData(data, unit=u.adu)
    frame.header[KEYWORDS["DPR_TYPE"]] = "OBJECT"
    frame.header[KEYWORDS["OBS_NAME"]] = "synthetic"
    frame.header[KEYWORDS["CP_TEMP_C"]] = TEST_CRYOSTAT_TEMP
    frame.header[KEYWORDS["NIR_TEMP_K"]] = TEST_NIR_DETECTOR_TEMP
    frame.header[KEYWORDS["VIS_TEMP_C"]] = TEST_VIS_DETECTOR_TEMP
    path = tmpPath / name
    frame.write(str(path), overwrite=True)
    return str(path)


def _prescan_recipe(log: Any, tmpPath: Path, *, recipeName: str = "soxs-stare", arm: str = "VIS", window: Any = None):
    """A recipe for the 200 x 40 science pixels beside the 30-column prescan strip."""
    calibrations = tmpPath / "calibrations"
    calibrations.mkdir(parents=True, exist_ok=True)
    CCDData(np.zeros((ROWS, SCIENCE_COLUMNS), dtype=np.int16), unit=u.dimensionless_unscaled).write(
        str(calibrations / "bpm.fits"), overwrite=True
    )
    recipe = _recipe(log, tmpPath, recipeName=recipeName, arm=arm, recipeSettings={"use_flat": False})
    recipe.dateObs = None
    recipe.detectorParams["science-pixels"] = {
        "rows": {"start": 0, "end": ROWS},
        "columns": {"start": PRESCAN_COLUMNS, "end": PRESCAN_COLUMNS + SCIENCE_COLUMNS},
    }
    recipe.detectorParams["prescan-pixels"] = PRESCAN_WINDOW if window is None else window
    return recipe


def _prepared_errors(recipe: Any, rawPath: str) -> np.ndarray:
    preparedPath = recipe._prepare_single_frame(rawPath)
    prepared = CCDData.read(preparedPath, hdu=0, unit=u.electron, hdu_uncertainty="ERRS")
    return prepared.uncertainty.array


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_a_science_frame_error_plane_uses_the_read_noise_of_its_prescan(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    rawPath = _raw_science_frame(tmp_path)

    # ACT
    errors = _prepared_errors(recipe, rawPath)

    # ASSERT: sqrt(PHOTON ELECTRONS + THE PRESCAN READ NOISE SQUARED)
    expected = np.sqrt(SKY_ADU * TEST_GAIN + PRESCAN_RON**2)
    assert errors.shape == (ROWS, SCIENCE_COLUMNS)
    assert errors[10, 10] == pytest.approx(expected, rel=TOLERANCE)
    assert errors[10, 10] > np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2)


def test_a_standard_star_recipe_also_uses_the_prescan(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, recipeName="soxs-nod-std")

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + PRESCAN_RON**2), rel=TOLERANCE)


def test_a_calibration_recipe_keeps_the_master_bias_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, recipeName="soxs-mflat")

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2), rel=1e-6)
    assert recipe.frameReadNoises == ()


def test_a_nir_science_recipe_keeps_its_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path, arm="NIR")

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2), rel=1e-6)
    assert recipe.frameReadNoises == ()
    assert _warnings(log) == []


def test_a_detector_with_no_prescan_window_keeps_the_master_bias_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    del recipe.detectorParams["prescan-pixels"]

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2), rel=1e-6)
    assert _warnings(log) == []


def test_a_frame_with_no_prescan_falls_back_to_the_master_bias_with_a_warning(log: Any, tmp_path: Path) -> None:
    # ARRANGE: THE PRESCAN WINDOW LIES OUTSIDE THE FRAME
    window = {"rows": {"start": 0, "end": ROWS}, "columns": {"start": 200, "end": 240}}
    recipe = _prescan_recipe(log, tmp_path, window=window)

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path, name="raw_a.fits"))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2), rel=1e-6)
    assert recipe.frameReadNoises == ((MASTER_BIAS_RON, False),)
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "raw_a.fits" in warnings[0]
    assert "outside the frame" in warnings[0]


def test_a_prescan_with_no_scatter_falls_back_to_the_master_bias_with_a_warning(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)

    # ACT
    errors = _prepared_errors(recipe, _raw_science_frame(tmp_path, prescanSigma=0.0))

    # ASSERT
    assert errors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + MASTER_BIAS_RON**2), rel=1e-6)
    assert "not a positive number" in _warnings(log)[0]


def test_each_frame_error_plane_keeps_its_own_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    quiet = _raw_science_frame(tmp_path, prescanSigma=2.0, name="quiet.fits")
    noisy = _raw_science_frame(tmp_path, prescanSigma=8.0, name="noisy.fits")

    # ACT
    quietErrors = _prepared_errors(recipe, quiet)
    noisyErrors = _prepared_errors(recipe, noisy)

    # ASSERT
    assert quietErrors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + (2.0 * TEST_GAIN) ** 2), rel=TOLERANCE)
    assert noisyErrors[10, 10] == pytest.approx(np.sqrt(SKY_ADU * TEST_GAIN + (8.0 * TEST_GAIN) ** 2), rel=TOLERANCE)
    assert [isMeasured for _, isMeasured in recipe.frameReadNoises] == [True, True]


def test_the_reduction_takes_the_root_mean_square_of_the_frame_read_noises(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.readNoise = ReadNoise(MASTER_BIAS_RON, 5, MASTER_BIAS_RON * np.sqrt(1.2), "master bias")
    recipe.frameReadNoises = ((4.0, True), (6.0, True))

    # ACT
    recipe._apply_science_read_noise()

    # ASSERT
    expected = float(np.sqrt((16.0 + 36.0) / 2.0))
    assert recipe.detectorParams["ron"].value == pytest.approx(expected)
    assert recipe.detectorParams["ron"].unit == u.electron
    assert recipe.readNoise.source == SOURCE_SCIENCE_PRESCAN
    assert recipe.readNoise.skyFit == pytest.approx(expected * np.sqrt(1.2))


def test_the_science_read_noise_is_recorded_as_a_qc_in_the_product_header(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.frameReadNoises = ((5.1, True),)

    # ACT
    recipe._apply_science_read_noise()

    # ASSERT
    row = recipe.qc[recipe.qc["qc_name"] == "SCI RON"].iloc[0]
    assert row["qc_value"] == pytest.approx(5.1)
    assert row["qc_unit"] == "electrons"
    assert bool(row["to_header"]) is True


def test_the_qc_can_be_recorded_again_after_the_tables_are_emptied(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    import pandas as pd

    recipe = _prescan_recipe(log, tmp_path)
    recipe.frameReadNoises = ((5.1, True),)
    recipe._apply_science_read_noise()
    recipe.qc, recipe.products = recipe._empty_qc_and_product_tables(pd)

    # ACT
    recipe._add_science_read_noise_qc()

    # ASSERT
    assert list(recipe.qc["qc_name"]) == ["SCI RON"]


def test_frames_that_all_fell_back_leave_the_read_noise_and_the_qc_unchanged(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.readNoise = ReadNoise(MASTER_BIAS_RON, 5, MASTER_BIAS_RON * np.sqrt(1.2), "master bias")
    recipe.frameReadNoises = ((MASTER_BIAS_RON, False),)

    # ACT
    recipe._apply_science_read_noise()
    recipe._add_science_read_noise_qc()

    # ASSERT
    assert recipe.detectorParams["ron"].value == MASTER_BIAS_RON
    assert recipe.readNoise.source == "master bias"
    assert len(recipe.qc) == 0


def test_a_frame_that_fell_back_still_counts_in_the_mean_read_noise(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _prescan_recipe(log, tmp_path)
    recipe.frameReadNoises = ((6.0, True), (MASTER_BIAS_RON, False))

    # ACT
    recipe._apply_science_read_noise()

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(float(np.sqrt((36.0 + MASTER_BIAS_RON**2) / 2.0)))
    assert recipe.readNoise.nCombined is None


def _stare_with_science_read_noise(log: Any, monkeypatch: pytest.MonkeyPatch) -> tuple[Any, list[list[str]]]:
    """A stare recipe whose read noise came from the prescans, and a record of the QCs each pass starts with."""
    import pandas as pd

    from soxspipe.recipes.soxs_stare import soxs_stare

    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.sofName = "synthetic"
    recipe.filenameTemplate = "synthetic.fits"
    recipe.productNameMarker = ""
    recipe.onTargetFrames = None
    recipe.recipeSettings = {}
    recipe.subtractSky = False
    recipe.recipeName = "soxs-stare"
    recipe.dateObs = None
    recipe.qc, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.readNoise = ReadNoise(5.1, None, 5.1, SOURCE_SCIENCE_PRESCAN)
    startingQcs: list[list[str]] = []

    def reduce_one(frames: list[Any], rawFrames: Any = False) -> tuple[str, bool]:
        startingQcs.append(list(recipe.qc["qc_name"]))
        return "product.fits", False

    monkeypatch.setattr(recipe, "_reduce_object_frames", reduce_one)
    monkeypatch.setattr(recipe, "report_output", lambda: recipe.qc.copy())
    return recipe, startingQcs


def test_each_frame_extracted_on_its_own_starts_with_the_science_read_noise_qc(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, startingQcs = _stare_with_science_read_noise(log, monkeypatch)
    recipe._add_science_read_noise_qc()

    # ACT
    recipe._extract_frames_on_their_own(["frame_1", "frame_2", "frame_3"], ["a.fits", "b.fits", "c.fits"])

    # ASSERT
    assert startingQcs == [["SCI RON"], ["SCI RON"], ["SCI RON"]]


def test_the_stacked_pass_before_the_frame_passes_starts_with_the_qc_and_so_does_each_frame(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, startingQcs = _stare_with_science_read_noise(log, monkeypatch)
    recipe._add_science_read_noise_qc()
    monkeypatch.setattr(recipe, "clean_up", lambda **_: None)

    # ACT
    recipe._reduce_stacked_then_each_object_frame(["frame_1", "frame_2"], ["a.fits", "b.fits"])

    # ASSERT
    assert startingQcs == [["SCI RON"], ["SCI RON"], ["SCI RON"]]


def test_the_cycle_inspection_pass_of_a_nodding_recipe_starts_with_the_qc(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    import pandas as pd

    from soxspipe.recipes.soxs_nod import soxs_nod

    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.sofName = "synthetic"
    recipe.filenameTemplate = "synthetic.fits"
    recipe.productNameMarker = ""
    recipe.extractionNotes = ()
    recipe.onTargetFrames = None
    recipe.onTargetCandidates = None
    recipe.recipeName = "soxs-nod"
    recipe.dateObs = None
    recipe.qc, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.readNoise = ReadNoise(5.1, None, 5.1, SOURCE_SCIENCE_PRESCAN)
    startingQcs: list[list[str]] = []

    def inspect(*_: object) -> Any:
        startingQcs.append(list(recipe.qc["qc_name"]))
        return recipe.qc.copy()

    monkeypatch.setattr(recipe, "_inspect_each_cycle", inspect)

    # ACT
    recipe._extract_each_cycle_as_well(([], [], [], []), recipe.qc.copy(), "order.fits", False)

    # ASSERT
    assert startingQcs == [["SCI RON"]]


class _PreparationStoppedError(Exception):
    """Raised by a stub to stop `prepare_frames` once the read noise has been applied."""


def test_preparing_the_frames_applies_the_prescan_read_noise_of_this_call_only(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    from types import SimpleNamespace

    recipe = _prescan_recipe(log, tmp_path)
    recipe.inputFrames = SimpleNamespace(files_filtered=lambda include_path: ["raw.fits"])
    recipe.verbose = False
    recipe.frameReadNoises = ((20.0, True),)

    def prepare_one(frame: str, save: bool) -> str:
        recipe.frameReadNoises = (*recipe.frameReadNoises, (5.1, True))
        return "raw_pre.fits"

    def stop() -> None:
        raise _PreparationStoppedError

    monkeypatch.setattr(recipe, "_prepare_single_frame", prepare_one)
    monkeypatch.setattr(recipe, "_write_cosmic_ray_mask_product", stop)

    # ACT
    with pytest.raises(_PreparationStoppedError):
        recipe.prepare_frames()

    # ASSERT: THE VALUE LEFT OVER FROM AN EARLIER CALL DOES NOT COUNT
    assert recipe.detectorParams["ron"].value == pytest.approx(5.1)
    assert list(recipe.qc["qc_name"]) == ["SCI RON"]
