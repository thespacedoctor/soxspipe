"""Synthetic FITS contracts owned by the base recipe."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
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
