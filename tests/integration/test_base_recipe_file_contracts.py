"""Synthetic FITS contracts owned by the base recipe."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories import pipeline_settings, prepared_fits, raw_fits, synthetic_ccd

pytestmark = pytest.mark.integration


def _recipe(tmp_path: Path, log: object) -> base_recipe:
    recipe = object.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.qc = __import__("pandas").DataFrame(
        columns=["qc_name", "qc_value", "qc_comment", "to_header"]
    )
    recipe.sofName = False
    recipe.startNightDate = "2024-01-02"
    recipe.recipeName = "soxs-mbias"
    recipe.outDir = str(tmp_path / "prepared")
    recipe.workspaceRootPath = str(tmp_path)
    recipe.inst = "XSHOOTER"
    return recipe


def test_update_fits_keywords_records_recipe_provenance(
    log: object,
    tmp_path: Path,
) -> None:
    """Recipe outputs retain stable product, provenance, and parameter cards."""
    calibrationPath = tmp_path / "MASTER_DARK_VIS.fits"
    calibrationPath.write_bytes(b"synthetic calibration")
    inventory = pd.DataFrame(
        [
            {
                "filename": "raw-science_pre.fits",
                "TYPE": "BIAS",
                "ARM": "VIS",
                "ESO PRO TYPE": float("nan"),
                "ESO PRO CATG": None,
                "file": str(tmp_path / "raw-science_pre.fits"),
            },
            {
                "filename": calibrationPath.name,
                "TYPE": "DARK",
                "ARM": "VIS",
                "ESO PRO TYPE": "REDUCED",
                "ESO PRO CATG": "MASTER_DARK_VIS",
                "file": str(calibrationPath),
            },
        ]
    )
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    recipe.detectorParams = {}
    recipe.imageType = "BIAS"
    recipe.recipeName = "soxs-mbias"
    recipe.settings = pipeline_settings(
        tmp_path,
        overrides={
            "soxs-mbias": {
                "frame-clipping-sigma": 3.0,
                "stacking": {"method": "median"},
            }
        },
    )
    recipe.rawFrames = type(
        "SyntheticRawFrames",
        (),
        {"to_pandas": lambda self: inventory.copy(deep=True)},
    )()
    outputFrame = synthetic_ccd(seed=101, prepared=True)

    recipe.update_fits_keywords(outputFrame, rawFrames=["raw-science.fits"])

    header = outputFrame.header
    assert header["ESO SEQ ARM"] == "VIS"
    assert header["ESO PRO TYPE"] == "REDUCED"
    assert header["ESO PRO CATG"] == "MASTER_BIAS_VIS"
    assert header["ESO PRO TECH"] == "IMAGE"
    assert header["ESO PRO REC1 ID"] == "soxs-mbias"
    assert str(header["ESO PRO REC1 PIPE ID"]).startswith("soxspipe/v")
    assert header["ESO PRO REC1 RAW1 NAME"] == "raw-science.fits"
    assert header["ESO PRO REC1 RAW1 CATG"] == "BIAS_VIS"
    assert header["ESO PRO REC1 CAL1 NAME"] == calibrationPath.name
    assert header["ESO PRO REC1 CAL1 CATG"] == "MASTER_DARK_VIS"
    assert header["ESO PRO REC1 CAL1 DATAMD5"] == "70237715eb696d0cf9dc67431f353b44"
    assert header["ESO PRO REC1 PARAM1 NAME"] == "frame-clipping-sigma"
    assert header["ESO PRO REC1 PARAM1 VALUE"] == 3.0
    assert header["ESO PRO REC1 PARAM2 NAME"] == "method"
    assert header["ESO PRO REC1 PARAM2 VALUE"] == "median"


def test_prepare_frames_preserves_processing_order_then_sorts_by_observation_time(
    tmp_path: Path,
    log: object,
) -> None:
    recipe = _recipe(tmp_path, log)
    recipe.settings = pipeline_settings(
        tmp_path,
        overrides={
            "summary-keys": {
                "default": [
                    "MJDOBS",
                    "SEQ_ARM",
                    "DPR_TYPE",
                    "LAMP1",
                    "LAMP2",
                    "LAMP3",
                    "LAMP4",
                    "LAMP5",
                    "LAMP6",
                    "LAMP7",
                ],
                "verbose": [],
                "nodding_extras": [],
            },
            "soxs-mbias": {"use_lacosmic": False},
        },
    )
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    laterPath = prepared_fits(tmp_path / "later.fits", seed=2)
    earlierPath = prepared_fits(tmp_path / "earlier.fits", seed=1)
    fits.setval(laterPath, recipe.kw("MJDOBS"), value=60_001.0)
    fits.setval(earlierPath, recipe.kw("MJDOBS"), value=60_000.0)
    sourcePaths = [str(laterPath), str(earlierPath)]
    recipe.inputFrames = SimpleNamespace(
        files_filtered=lambda include_path: sourcePaths.copy()
    )
    recipe.recipeSettings = {"use_lacosmic": False}
    recipe.verbose = False
    recipe.arm = "VIS"
    processedPaths: list[str] = []

    def prepare(frame: str, save: bool) -> str:
        processedPaths.append(frame)
        return frame

    recipe._prepare_single_frame = prepare

    collection = recipe.prepare_frames()

    assert processedPaths == sourcePaths
    assert list(collection.summary[recipe.kw("MJDOBS")]) == [60_000.0, 60_001.0]
    assert list(collection.summary["filename"]) == ["earlier.fits", "later.fits"]


def test_write_round_trip_preserves_prepared_extensions(
    tmp_path: Path,
    log: object,
) -> None:
    recipe = _recipe(tmp_path, log)
    frame = synthetic_ccd(prepared=True)

    outputPath = recipe._write(
        frame=frame,
        filedir=str(tmp_path),
        filename="prepared.fits",
        product=False,
    )

    with fits.open(outputPath) as hdus:
        assert [hdu.name for hdu in hdus] == ["FLUX", "QUAL", "ERRS"]
    restored = CCDData.read(
        outputPath,
        hdu="FLUX",
        unit=u.electron,
        hdu_uncertainty="ERRS",
        hdu_mask="QUAL",
        key_uncertainty_type="UTYPE",
    )
    assert restored.shape == (32, 32)
    assert restored.unit == u.electron
    assert restored.header["ESO SEQ ARM"] == "VIS"
    assert restored.mask[0, 0]
    np.testing.assert_allclose(restored.data, frame.data)
    assert np.allclose(restored.uncertainty.array, 1.5)


def test_prepare_single_frame_uses_shape_matched_bad_pixel_map(
    tmp_path: Path,
    log: object,
) -> None:
    recipe = _recipe(tmp_path, log)
    calibrationPath = tmp_path / "calibrations"
    calibrationPath.mkdir()
    badPixelData = np.zeros((32, 32), dtype=np.uint8)
    badPixelData[1, 2] = 1
    fits.PrimaryHDU(badPixelData).writeto(calibrationPath / "bitmap.fits")
    recipe.calibrationRootPath = str(calibrationPath)
    recipe.detectorParams = {
        "gain": np.float32(2.0),
        "ron": np.float32(3.0),
        "bad-pixel map": {"1x1": "bitmap.fits"},
    }
    recipe.arm = "VIS"
    recipe.recipeSettings = {"use_lacosmic": False}
    recipe.debug = False
    recipe.xsh2soxs = lambda frame: frame
    recipe._trim_frame = lambda frame: frame
    inputPath = raw_fits(tmp_path / "raw.fits")

    preparedPath = recipe._prepare_single_frame(str(inputPath))
    prepared = CCDData.read(
        preparedPath,
        hdu="FLUX",
        unit=u.electron,
        hdu_uncertainty="ERRS",
        hdu_mask="QUAL",
        key_uncertainty_type="UTYPE",
    )

    assert preparedPath == str(tmp_path / "prepared" / "raw_pre.fits")
    assert prepared.shape == (32, 32)
    assert prepared.unit == u.electron
    assert prepared.mask[1, 2]
    assert "SXSPRE" in prepared.header
    assert np.issubdtype(prepared.uncertainty.array.dtype, np.floating)
    assert prepared.uncertainty.array.dtype.itemsize == np.dtype(np.float32).itemsize
    np.testing.assert_allclose(prepared.uncertainty.array, 3.0)
    np.testing.assert_allclose(prepared.data, fits.getdata(inputPath) * 2.0)


def test_prepare_single_frame_leaves_existing_prepared_layout_unchanged(
    tmp_path: Path,
    log: object,
) -> None:
    """Prepared FITS input remains a stable path instead of being prepared twice."""
    recipe = _recipe(tmp_path, log)
    recipe.detectorParams = {}
    inputPath = prepared_fits(tmp_path / "already-prepared.fits", seed=7)

    result = recipe._prepare_single_frame(str(inputPath))

    assert result == str(inputPath)
    with fits.open(inputPath) as hdus:
        assert [hdu.name for hdu in hdus] == ["FLUX", "QUAL", "ERRS"]
