"""Characterization of `base_recipe._prepare_single_frame`.

These tests pin what the method does today, before it is split. The method is
225 lines long, runs once per input frame, and most of its uncovered lines are
branches the offline suite never reaches: the two early returns for a frame
that cannot be read as an image, the dummy bad-pixel map it writes when the
real one is missing, the cosmic-ray cleaning branch, and the SOXS temperature
keywords it copies onto the recipe.

The numerical assertions pin the gain conversion and both uncertainty maps,
which is what a split of this method could silently change.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

# THE GAIN AND READNOISE THE FAKE DETECTOR PARAMETERS CARRY. BOTH ARE READ
# BACK IN THE ASSERTIONS, SO THEY ARE NAMED RATHER THAN REPEATED.
TEST_GAIN = 2.0
TEST_RON = 3.0

# THE TEMPERATURES EVERY SYNTHETIC FRAME CARRIES. A SOXS FRAME READS THE
# CRYOSTAT TEMPERATURE FOR EVERY ARM AND THE DETECTOR TEMPERATURE FOR TWO OF
# THEM, SO ALL THREE KEYWORDS ARE PRESENT ON EVERY FRAME.
TEST_CRYOSTAT_TEMP = -40.5
TEST_NIR_DETECTOR_TEMP = 80.5
TEST_VIS_DETECTOR_TEMP = -110.25

# THE FITS KEYWORDS THE FAKE LOOKUP SERVES. THEY ARE SHORT NAMES BECAUSE
# ASTROPY REQUIRES `HIERARCH` FOR ANYTHING LONGER THAN EIGHT CHARACTERS, AND
# THE METHOD UNDER TEST ONLY EVER SEES WHAT `self.kw` RETURNS.
KEYWORDS = {
    "DPR_TYPE": "DPRTYPE",
    "OBS_NAME": "OBSNAME",
    "WIN_BINX": "WINBINX",
    "WIN_BINY": "WINBINY",
    "CP_TEMP_C": "CPTEMP",
    "NIR_TEMP_K": "NIRTEMP",
    "VIS_TEMP_C": "VISTEMP",
}

QC_COLUMNS = [
    "soxspipe_recipe",
    "qc_name",
    "qc_value",
    "qc_unit",
    "qc_order",
    "qc_comment",
    "obs_date_utc",
    "reduction_date_utc",
    "to_header",
]


def _detector_params(badPixelMapName: str) -> dict[str, Any]:
    """Return detector parameters covering an unrotated, unbinned 8x8 detector."""
    return {
        "clockwise-rotation": 0.0,
        "binning": [1, 1],
        "science-pixels": {
            "rows": {"start": 0, "end": 8},
            "columns": {"start": 0, "end": 8},
        },
        "gain": TEST_GAIN * u.electron / u.adu,
        "ron": TEST_RON * u.electron,
        "bad-pixel map": {"1x1": badPixelMapName},
    }


def _raw_frame(tmpPath: Path, *, dprType: str = "BIAS", **headerExtras: Any) -> str:
    """Write a single-extension raw frame in ADU and return its path."""
    data = np.arange(64, dtype=np.float32).reshape(8, 8)
    frame = CCDData(data, unit=u.adu)
    frame.header[KEYWORDS["DPR_TYPE"]] = dprType
    frame.header[KEYWORDS["OBS_NAME"]] = "synthetic"
    frame.header[KEYWORDS["CP_TEMP_C"]] = TEST_CRYOSTAT_TEMP
    frame.header[KEYWORDS["NIR_TEMP_K"]] = TEST_NIR_DETECTOR_TEMP
    frame.header[KEYWORDS["VIS_TEMP_C"]] = TEST_VIS_DETECTOR_TEMP
    for keyword, value in headerExtras.items():
        frame.header[keyword] = value
    framePath = tmpPath / f"raw_{dprType.replace(',', '_')}.fits"
    frame.write(str(framePath), overwrite=True)
    return str(framePath)


def _bad_pixel_map(calibrationPath: Path, name: str = "bpm.fits") -> str:
    """Write an all-good bad-pixel map and return its file name."""
    calibrationPath.mkdir(parents=True, exist_ok=True)
    mask = CCDData(np.zeros((8, 8), dtype=np.int16), unit=u.dimensionless_unscaled)
    mask.write(str(calibrationPath / name), overwrite=True)
    return name


def _recipe(
    log: Any,
    tmpPath: Path,
    *,
    badPixelMapName: str = "bpm.fits",
    arm: str = "VIS",
    inst: str = "SOXS",
    recipeName: str = "soxs-mbias",
    recipeSettings: dict[str, Any] | None = None,
) -> base_recipe:
    """Build a recipe carrying only the attributes frame preparation reads."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = lambda name: KEYWORDS.get(name, name)
    recipe.settings = {"instrument": "soxs", "workspace-root-dir": str(tmpPath)}
    recipe.recipeName = recipeName
    recipe.recipeSettings = {} if recipeSettings is None else recipeSettings
    recipe.arm = arm
    recipe.inst = inst
    recipe.debug = False
    recipe.detectorParams = _detector_params(badPixelMapName)
    recipe.calibrationRootPath = str(tmpPath / "calibrations")
    recipe.workspaceRootPath = str(tmpPath / "workspace")
    recipe.outDir = str(tmpPath / "tmp")
    recipe.sofName = "synthetic"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = pd.DataFrame({column: [] for column in QC_COLUMNS})
    recipe.generateReponseCurve = False
    return recipe


def test_a_bias_frame_is_gain_corrected_and_given_a_flat_readnoise_uncertainty(
    log: Any,
    tmp_path: Path,
) -> None:
    """A BIAS frame's uncertainty is the readnoise alone, at every pixel."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)

    # ACT
    preparedPath = recipe._prepare_single_frame(rawPath)

    # ASSERT
    prepared = CCDData.read(
        preparedPath,
        hdu=0,
        unit=u.electron,
        hdu_uncertainty="ERRS",
        hdu_mask="QUAL",
    )
    assert prepared.data[1, 1] == pytest.approx(9.0 * TEST_GAIN, rel=1e-12)
    assert np.allclose(prepared.uncertainty.array, TEST_RON, rtol=1e-12)
    assert prepared.mask.sum() == 0


def test_a_non_bias_frame_takes_a_poisson_plus_readnoise_uncertainty(
    log: Any,
    tmp_path: Path,
) -> None:
    """A non-BIAS frame's uncertainty is `sqrt(counts + readnoise squared)`."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path, dprType="LAMP,FLAT")

    # ACT
    preparedPath = recipe._prepare_single_frame(rawPath)

    # ASSERT
    prepared = CCDData.read(preparedPath, hdu=0, unit=u.electron, hdu_uncertainty="ERRS")
    expected = np.sqrt(9.0 * TEST_GAIN + TEST_RON**2)
    assert prepared.uncertainty.array[1, 1] == pytest.approx(expected, rel=1e-6)


def test_an_already_prepared_frame_is_returned_untouched(
    log: Any,
    tmp_path: Path,
) -> None:
    """A frame carrying `SXSPRE` is handed back as its own path, byte for byte unchanged."""
    # ARRANGE
    import hashlib

    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path, SXSPRE="2024-01-02T03:04:05.678901")
    digestBefore = hashlib.sha256(Path(rawPath).read_bytes()).hexdigest()

    # ACT
    result = recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert result == rawPath
    # THE FRAME IS NOT REWRITTEN, AND NOTHING IS WRITTEN BESIDE IT EITHER, SO A
    # REFACTOR THAT PREPARED THE FRAME ANYWAY AND THEN RETURNED THE OLD PATH
    # WOULD FAIL HERE.
    assert hashlib.sha256(Path(rawPath).read_bytes()).hexdigest() == digestBefore
    assert not (tmp_path / "tmp").exists()


def test_a_corrupted_frame_is_dropped_from_the_reduction(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A `buffer is too small` read failure returns `None` and warns."""
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)

    def failing_read(*_: object, **__: object) -> CCDData:
        raise TypeError("buffer is too small for requested array")

    monkeypatch.setattr(CCDData, "read", failing_read)

    # ACT
    result = recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert result is None
    assert any("likely corrupted" in str(message) for message in log.messages)


def test_a_binary_table_is_passed_through_by_path(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Any other read `TypeError` returns the path, unprepared."""
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)

    def failing_read(*_: object, **__: object) -> CCDData:
        raise TypeError("no data in this HDU")

    monkeypatch.setattr(CCDData, "read", failing_read)

    # ACT
    result = recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert result == rawPath
    assert any("FITS Binary Table" in str(message) for message in log.messages)


def test_a_missing_bad_pixel_map_is_written_as_zeros_and_then_fails_the_frame(
    log: Any,
    tmp_path: Path,
) -> None:
    """The method writes a dummy map to the missing path, then raises `OSError`.

    Pinned as it stands: writing the file it has just declared missing and
    then failing anyway is the current behavior, and a later run finds the
    zeroed map in place.
    """
    # ARRANGE
    (tmp_path / "calibrations").mkdir(parents=True, exist_ok=True)
    recipe = _recipe(log, tmp_path, badPixelMapName="absent.fits")
    rawPath = _raw_frame(tmp_path)

    # ACT / ASSERT
    with pytest.raises(OSError, match="does not exist on this machine"):
        recipe._prepare_single_frame(rawPath)

    writtenMap = tmp_path / "calibrations" / "absent.fits"
    assert writtenMap.exists()
    assert CCDData.read(str(writtenMap), unit=u.adu).data.sum() == 0


def test_cosmic_ray_cleaning_runs_when_the_recipe_settings_ask_for_it(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`use_lacosmic` with an arm sigma sends the frame through `cosmicray_lacosmic`."""
    # ARRANGE
    import ccdproc

    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(
        log,
        tmp_path,
        recipeSettings={"use_lacosmic": True, "vis": {"lacosmic-sigma": 4.5}},
    )
    rawPath = _raw_frame(tmp_path)
    calls: list[dict[str, Any]] = []

    def spying_lacosmic(frame: CCDData, **kwargs: Any) -> CCDData:
        calls.append(kwargs)
        return frame

    monkeypatch.setattr(ccdproc, "cosmicray_lacosmic", spying_lacosmic)

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert len(calls) == 1
    assert calls[0]["sigclip"] == pytest.approx(4.5)
    assert calls[0]["gain_apply"] is False
    assert calls[0]["niter"] == 2
    assert calls[0]["cleantype"] == "meanmask"


def test_a_standard_star_frame_switches_on_the_response_curve(
    log: Any,
    tmp_path: Path,
) -> None:
    """A `STD,FLUX` frame in a standard recipe sets `generateReponseCurve`."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(
        log,
        tmp_path,
        recipeName="soxs-stare-std",
        recipeSettings={"use_flat": True},
    )
    rawPath = _raw_frame(tmp_path, dprType="STD,FLUX")

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert recipe.generateReponseCurve is True


@pytest.mark.parametrize(
    ("arm", "temperatureValue"),
    [
        ("NIR", TEST_NIR_DETECTOR_TEMP),
        ("VIS", TEST_VIS_DETECTOR_TEMP),
    ],
)
def test_soxs_frames_record_the_arm_detector_temperature(
    log: Any,
    tmp_path: Path,
    arm: str,
    temperatureValue: float,
) -> None:
    """A SOXS frame copies the cryostat and arm detector temperatures onto the recipe."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path, arm=arm)
    rawPath = _raw_frame(tmp_path)

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert recipe.cptemp == pytest.approx(TEST_CRYOSTAT_TEMP)
    assert recipe.detectorTemp == pytest.approx(temperatureValue)


def test_a_uvb_frame_records_no_detector_temperature(
    log: Any,
    tmp_path: Path,
) -> None:
    """An arm with no temperature keyword records `None`."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path, arm="UVB")
    rawPath = _raw_frame(tmp_path)

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert recipe.detectorTemp is None


def test_an_xshooter_frame_records_no_temperatures_at_all(
    log: Any,
    tmp_path: Path,
) -> None:
    """Outside SOXS the temperature keywords are never read."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path, inst="XSHOOTER")
    rawPath = _raw_frame(tmp_path)

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert not hasattr(recipe, "cptemp")
    assert not hasattr(recipe, "detectorTemp")


def test_saving_redirects_the_prepared_frame_to_the_workspace_root(
    log: Any,
    tmp_path: Path,
) -> None:
    """`save=True` writes beside the workspace rather than into the scratch directory."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)

    # ACT
    preparedPath = recipe._prepare_single_frame(rawPath, save=True)

    # ASSERT
    assert Path(preparedPath).parent == tmp_path / "workspace"
    assert Path(preparedPath).name.endswith("_pre.fits")


def test_the_prepared_frame_carries_a_microsecond_resolution_pre_timestamp(
    log: Any,
    tmp_path: Path,
) -> None:
    """`SXSPRE` records when the frame was prepared, to microseconds."""
    # ARRANGE
    import re

    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)

    # ACT
    preparedPath = recipe._prepare_single_frame(rawPath)

    # ASSERT
    prepared = CCDData.read(preparedPath, hdu=0, unit=u.electron)
    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{6}$", prepared.header["SXSPRE"])
