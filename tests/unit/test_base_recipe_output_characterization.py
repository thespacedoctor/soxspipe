"""Characterization of read-out noise, header stamping and frame writing."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import MaskedColumn, Table

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories import instrument_header, pipeline_settings, qc_row, qc_table

pytestmark = pytest.mark.unit


class StubInputFrames:
    """Input-frame collection boundary returning fixed CCDData objects."""

    def __init__(self, ccdList: list[CCDData]) -> None:
        self.files = [f"synthetic_{index}.fits" for index in range(len(ccdList))]
        self._ccdList = ccdList
        self.ccdKwargs: list[dict[str, Any]] = []

    def ccds(self, ccd_kwargs: dict[str, Any] | None = None) -> list[CCDData]:
        """Return the fixed frames and record the requested extension names."""
        self.ccdKwargs.append(dict(ccd_kwargs or {}))
        return list(self._ccdList)


def _frame(data: np.ndarray, mask: np.ndarray | None = None) -> CCDData:
    shape = data.shape
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros(shape, dtype=bool) if mask is None else mask.copy(),
        uncertainty=StdDevUncertainty(np.full(shape, 1.0), unit=u.electron),
    )


def _recipe(log: Any) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = qc_table().iloc[0:0].copy()
    return recipe


def _noise_frames(*, maskFirstPixel: bool = False, extraFrame: bool = False) -> list[CCDData]:
    """Return frames whose difference has a known standard deviation."""
    rng = np.random.default_rng(11)
    first = rng.normal(loc=100.0, scale=3.0, size=(16, 16))
    second = rng.normal(loc=100.0, scale=3.0, size=(16, 16))
    mask = np.zeros((16, 16), dtype=bool)
    if maskFirstPixel:
        mask[0, 0] = True
    frames = [_frame(first, mask), _frame(second, mask)]
    if extraFrame:
        # A THIRD FRAME WIDE ENOUGH TO MOVE THE RESULT IF IT WERE MEASURED
        frames.append(_frame(rng.normal(loc=100.0, scale=90.0, size=(16, 16))))
    return frames


def test_qc_ron_measures_raw_noise_from_the_first_two_input_frames(log: Any) -> None:
    """Raw read-out noise is half the clipped spread of the first frame pair."""
    # ARRANGE
    recipe = _recipe(log)
    # THE THIRD FRAME IS NINETY TIMES NOISIER, SO A MEASUREMENT THAT READ IT
    # COULD NOT RETURN THE SAME VALUE AS THE TWO-FRAME CASE BELOW
    recipe.inputFrames = StubInputFrames(_noise_frames(extraFrame=True))

    # ACT
    rawRon, masterRon = recipe.qc_ron(frameType="MBIAS", frameName="master bias")

    # ASSERT
    assert rawRon == pytest.approx(2.979132461386683, rel=1e-12)
    assert masterRon is None
    assert recipe.inputFrames.ccdKwargs == [
        {
            "hdu_uncertainty": "ERRS",
            "hdu_mask": "QUAL",
            "hdu_flags": "FLAGS",
            "key_uncertainty_type": "UTYPE",
        }
    ]
    rawRow = recipe.qc.set_index("qc_name").loc["RAW RON"]
    assert rawRow["qc_value"] == pytest.approx(2.979132461386683, rel=1e-12)
    assert rawRow["qc_comment"] == "[e-] RON in single BIAS"
    assert rawRow["qc_unit"] == "electrons"


def test_qc_ron_measures_master_noise_against_the_raw_pair_mask(log: Any) -> None:
    """Master read-out noise is the spread of the master frame under the raw mask."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inputFrames = StubInputFrames(_noise_frames(maskFirstPixel=True))
    masterFrame = _frame(np.full((16, 16), 5.0))
    # THIS PIXEL IS MASKED IN BOTH RAW FRAMES, SO THE MASK THE RAW MEASUREMENT
    # RETURNS MUST EXCLUDE IT. AN UNMASKED MEASUREMENT WOULD RETURN ABOUT 62.
    masterFrame.data[0, 0] = 1000.0
    masterFrame.data[1, 1] = 9.0

    # ACT
    rawRon, masterRon = recipe.qc_ron(
        frameType="MBIAS",
        frameName="master bias",
        masterFrame=masterFrame,
    )

    # ASSERT
    assert rawRon == pytest.approx(2.9832380152474314, rel=1e-12)
    assert masterRon == pytest.approx(0.24999807765504675, rel=1e-12)
    # THE MASTER BRANCH RECORDS NO QC ROW OF ITS OWN, ONLY THE RAW ROW
    assert recipe.qc["qc_name"].tolist() == ["RAW RON"]


def test_qc_ron_rejects_input_frames_with_no_measurable_spread(log: Any) -> None:
    """Identical input frames raise rather than report zero read-out noise."""
    # ARRANGE
    recipe = _recipe(log)
    flat = np.full((8, 8), 12.0)
    recipe.inputFrames = StubInputFrames([_frame(flat), _frame(flat)])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="appear to be corrupted"):
        recipe.qc_ron(frameType="MBIAS", frameName="master bias")


def test_qc_ron_reports_nothing_when_a_single_input_frame_is_supplied(
    log: Any,
) -> None:
    """A single input frame yields no raw measurement and no QC rows."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.inputFrames = StubInputFrames([_frame(np.full((4, 4), 3.0))])

    # ACT
    rawRon, masterRon = recipe.qc_ron(frameType="MBIAS", frameName="master bias")

    # ASSERT
    assert rawRon is False
    assert masterRon is None
    assert recipe.qc.empty


def _keyword_recipe(tmp_path: Path, log: Any, rawFrames: Table) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.imageType = "BIAS"
    recipe.detectorParams = {}
    recipe.recipeName = "soxs-mbias"
    recipe.rawFrames = rawFrames
    return recipe


def _raw_frame_table(tmp_path: Path, calibrationPath: Path) -> Table:
    """Return a raw-frame summary carrying one raw and one calibration row."""
    # A RAW ROW IS ONE WHOSE PRO TYPE IS MISSING, WHICH REACHES PANDAS AS NaN
    proType = MaskedColumn(["", "REDUCED"], name="ESO PRO TYPE", mask=[True, False])
    return Table(
        {
            "filename": ["raw_one_pre.fits", "master_bias.fits"],
            "file": [str(tmp_path / "raw_one.fits"), str(calibrationPath)],
            "TYPE": ["BIAS", "MASTER_BIAS"],
            "ARM": ["VIS", "VIS"],
            "ESO PRO TYPE": proType,
            "ESO PRO CATG": ["", "MASTER_BIAS_VIS"],
        },
        masked=True,
    )


def test_update_fits_keywords_stamps_raw_calibration_and_parameter_records(
    tmp_path: Path,
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Raw, calibration and recipe-parameter records reach the product header."""
    # ARRANGE
    calibrationPath = tmp_path / "master_bias.fits"
    calibrationPath.write_bytes(b"synthetic calibration payload")
    recipe = _keyword_recipe(tmp_path, log, _raw_frame_table(tmp_path, calibrationPath))
    monkeypatch.setattr(
        base_recipe,
        "get_recipe_settings",
        lambda self: {
            "stacked-clipping-sigma": 3.0,
            "uvb": {"ignored": True},
            "clipping": {"iterations": 2},
            "lamps": ["QTH", "D2"],
        },
        raising=False,
    )
    frame = CCDData(np.zeros((2, 2)), unit=u.electron, meta=instrument_header())

    # ACT
    result = recipe.update_fits_keywords(frame)

    # ASSERT
    assert result is None
    header = frame.header
    assert header["ESO SEQ ARM"] == "VIS"
    assert header["ESO PRO TYPE"] == "REDUCED"
    assert header["ESO PRO CATG"] == "MASTER_BIAS_VIS"
    assert header["ESO PRO TECH"] == "IMAGE"
    assert header["ESO PRO REC1 RAW1 NAME"] == "raw_one.fits"
    assert header["ESO PRO REC1 RAW1 CATG"] == "BIAS_VIS"
    assert "ESO PRO REC1 RAW2 NAME" not in header
    assert header["ESO PRO REC1 CAL1 NAME"] == "master_bias.fits"
    assert header["ESO PRO REC1 CAL1 CATG"] == "MASTER_BIAS_VIS"
    assert header["ESO PRO REC1 CAL1 DATAMD5"] == "2b6ee51a5e887abfe41e7630fb17f214"
    assert header["ESO PRO REC1 PARAM1 NAME"] == "stacked-clipping-sigma"
    assert header["ESO PRO REC1 PARAM1 VALUE"] == 3.0
    assert header["ESO PRO REC1 PARAM2 NAME"] == "iterations"
    assert header["ESO PRO REC1 PARAM2 VALUE"] == 2
    assert header["ESO PRO REC1 PARAM3 NAME"] == "lamps"
    assert header["ESO PRO REC1 PARAM3 VALUE"] == "QTH, D2"
    assert header["ESO PRO REC1 ID"] == "soxs-mbias"
    assert header["ESO PRO REC1 PIPE ID"].startswith("soxspipe/v")


def test_update_fits_keywords_limits_raw_records_to_the_requested_frames(
    tmp_path: Path,
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A raw-frame filter removes unlisted raw records but keeps calibrations."""
    # ARRANGE
    calibrationPath = tmp_path / "master_bias.fits"
    calibrationPath.write_bytes(b"synthetic calibration payload")
    recipe = _keyword_recipe(tmp_path, log, _raw_frame_table(tmp_path, calibrationPath))
    monkeypatch.setattr(base_recipe, "get_recipe_settings", lambda self: {}, raising=False)
    frame = CCDData(np.zeros((2, 2)), unit=u.electron, meta=instrument_header())

    # ACT
    recipe.update_fits_keywords(frame, rawFrames=["other_frame.fits"])

    # ASSERT
    assert "ESO PRO REC1 RAW1 NAME" not in frame.header
    assert frame.header["ESO PRO REC1 CAL1 NAME"] == "master_bias.fits"


def _write_recipe(tmp_path: Path, log: Any) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.sofName = "synthetic_product"
    recipe.recipeName = "soxs-mbias"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = qc_table().iloc[0:0].copy()
    return recipe


def test_write_stamps_flagged_qc_values_and_builds_the_product_path(
    tmp_path: Path,
    log: Any,
) -> None:
    """A product write lands under the recipe directory carrying header QC rows."""
    # ARRANGE
    recipe = _write_recipe(tmp_path, log)
    recipe.qc = pd.concat(
        [
            qc_row(
                recipeName="soxs-mbias",
                name="RON",
                value=1.25,
                unit="electron",
                comment="[e-] read noise",
            ),
            qc_row(
                recipeName="soxs-mbias",
                name="BIAS STRUCTURE",
                value=2.5,
                unit="electron",
                comment="[e-] structure",
                toHeader=False,
            ),
        ],
        ignore_index=True,
    )
    frame = _frame(np.full((4, 4), 7.0))

    # ACT
    filepath = recipe._write(frame, str(tmp_path))

    # ASSERT
    expectedDirectory = tmp_path / "reduced" / "2024-01-02" / "soxs-mbias"
    assert Path(filepath) == expectedDirectory / "synthetic_product.fits"
    assert Path(filepath).exists()
    written = CCDData.read(filepath, hdu_uncertainty="ERRS", hdu_mask="QUAL")
    assert written.header["ESO QC RON"] == pytest.approx(1.25)
    assert written.header.comments["ESO QC RON"] == "[e-] read noise"
    assert "ESO QC BIAS STRUCTURE" not in written.header
    assert written.data[0, 0] == pytest.approx(7.0)


def test_write_can_save_a_named_non_product_frame_with_masked_pixels_set(
    tmp_path: Path,
    log: Any,
) -> None:
    """A non-product write keeps the given directory and rewrites masked pixels."""
    # ARRANGE
    recipe = _write_recipe(tmp_path, log)
    mask = np.zeros((4, 4), dtype=bool)
    mask[1, 1] = True
    frame = _frame(np.full((4, 4), 7.0), mask)

    # ACT
    filepath = recipe._write(
        frame,
        str(tmp_path),
        filename="intermediate.fits",
        product=False,
        maskToZero=True,
    )

    # ASSERT
    assert Path(filepath) == tmp_path / "intermediate.fits"
    written = CCDData.read(filepath, hdu_uncertainty="ERRS", hdu_mask="QUAL")
    # THE METHOD WRITES 1, NOT 0, INTO MASKED PIXELS DESPITE ITS ARGUMENT NAME
    assert written.data[1, 1] == pytest.approx(1.0)
    assert written.data[0, 0] == pytest.approx(7.0)
    # THE MASK ITSELF SURVIVES THE ROUND TRIP, IN ITS OWN QUAL EXTENSION
    assert bool(written.mask[1, 1]) is True
    assert bool(written.mask[0, 0]) is False
