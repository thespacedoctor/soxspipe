"""Contracts for the cosmic-ray mask QC product (DY-1721).

Frame preparation runs L.A.Cosmic on each raw frame. The pixels it newly flags are collected on the recipe and
written once, as one multi-extension FITS file with one image extension per cleaned raw frame, together with a
QC row holding the total number of flagged pixels.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.recipes import base_recipe
from tests.unit.test_base_recipe_prepare_frame_characterization import (
    KEYWORDS,
    QC_COLUMNS,
)
from tests.unit.test_base_recipe_prepare_frame_characterization import (
    _recipe as _stub_recipe,
)

pytestmark = pytest.mark.unit

FRAME_SIZE = 64
SKY_LEVEL = 100.0
NOISE_LEVEL = 3.0
SPIKE_LEVEL = 60000.0
BAD_PIXEL = (10, 10)
FIRST_RAY = (30, 40)
SECOND_RAY = (50, 20)
MASK_PRODUCT_LABEL = "LACOSMIC_MASK"
MASK_QC_NAME = "COSMIC RAY PIX NUM"
PRODUCT_COLUMNS = [
    "soxspipe_recipe",
    "product_label",
    "file_name",
    "file_type",
    "obs_date_utc",
    "reduction_date_utc",
    "file_path",
    "label",
]
LACOSMIC_SETTINGS = {"use_lacosmic": True, "vis": {"lacosmic-sigma": 25}}


def _raw_object_frame(
    tmpPath: Path,
    name: str,
    rays: tuple[tuple[int, int], ...],
    mjdObs: float,
    seed: int,
) -> str:
    """Write a noisy 64x64 raw frame with a spike on the bad pixel and on each ray position."""
    data = np.random.default_rng(seed).normal(SKY_LEVEL, NOISE_LEVEL, (FRAME_SIZE, FRAME_SIZE)).astype(np.float32)
    data[BAD_PIXEL] = SPIKE_LEVEL
    for ray in rays:
        data[ray] = SPIKE_LEVEL
    frame = CCDData(data, unit=u.adu)
    frame.header[KEYWORDS["DPR_TYPE"]] = "OBJECT"
    frame.header[KEYWORDS["OBS_NAME"]] = "synthetic"
    frame.header[KEYWORDS["CP_TEMP_C"]] = -40.0
    frame.header[KEYWORDS["NIR_TEMP_K"]] = 80.0
    frame.header[KEYWORDS["VIS_TEMP_C"]] = -110.0
    frame.header["MJD-OBS"] = mjdObs
    frame.header["DATE-OBS"] = "2024-01-02T03:04:05.000"
    path = tmpPath / name
    frame.write(str(path), overwrite=True)
    return str(path)


def _bad_pixel_map(calibrationPath: Path) -> None:
    """Write a 64x64 bad-pixel map with one bad pixel."""
    calibrationPath.mkdir(parents=True, exist_ok=True)
    mask = np.zeros((FRAME_SIZE, FRAME_SIZE), dtype=np.int16)
    mask[BAD_PIXEL] = 1
    CCDData(mask, unit=u.dimensionless_unscaled).write(str(calibrationPath / "bpm.fits"), overwrite=True)


def _recipe(log: Any, tmpPath: Path, recipeSettings: dict[str, Any] | None = None) -> base_recipe:
    """Build a stub recipe that can prepare frames, collect cosmic-ray masks and register QC products."""
    _bad_pixel_map(tmpPath / "calibrations")
    recipe = _stub_recipe(
        log,
        tmpPath,
        recipeName="soxs-stare",
        recipeSettings=LACOSMIC_SETTINGS if recipeSettings is None else recipeSettings,
    )
    recipe.detectorParams["science-pixels"] = {
        "rows": {"start": 0, "end": FRAME_SIZE},
        "columns": {"start": 0, "end": FRAME_SIZE},
    }
    keywords = {**KEYWORDS, "MJDOBS": "MJD-OBS", "DATE_OBS": "DATE-OBS"}
    recipe.kw = lambda name: keywords.get(name, name)
    recipe.sofName = "OBJECT_VIS"
    recipe.qcDir = str(tmpPath / "qc")
    Path(recipe.qcDir).mkdir()
    recipe.dateObs = None
    recipe.cosmicRayMasks = []
    recipe.products = pd.DataFrame({column: [] for column in PRODUCT_COLUMNS})
    recipe.qc = pd.DataFrame({column: [] for column in QC_COLUMNS})
    return recipe


def _mask_path(recipe: base_recipe) -> Path:
    return Path(recipe.qcDir) / f"{recipe.sofName}_LACOSMIC_MASK.fits"


def test_a_cosmic_ray_is_written_to_the_qc_mask_and_the_bad_pixel_is_not(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_object_frame(tmp_path, "raw_1.fits", (FIRST_RAY,), mjdObs=60000.5, seed=1)
    recipe._prepare_single_frame(rawPath)

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    with fits.open(_mask_path(recipe)) as hdus:
        assert len(hdus) == 2
        mask = hdus[1].data
        assert mask.dtype == np.uint8
        assert mask[FIRST_RAY] == 1
        assert mask[BAD_PIXEL] == 0
        assert mask.sum() == 1
        header = hdus[1].header
        assert header["ORIGFILE"] == "raw_1.fits"
        assert header["MJD-OBS"] == pytest.approx(60000.5)
        assert header["DATE-OBS"] == "2024-01-02T03:04:05.000"
        assert header["NCRPIX"] == 1


def test_the_mask_product_and_its_qc_value_are_registered(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_1.fits", (FIRST_RAY,), 60000.5, seed=1))

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    assert len(recipe.products) == 1
    product = recipe.products.iloc[0]
    assert product["product_label"] == MASK_PRODUCT_LABEL
    assert product["label"] == "QC"
    assert product["file_name"] == "OBJECT_VIS_LACOSMIC_MASK.fits"
    assert product["file_path"] == str(_mask_path(recipe))
    assert product["obs_date_utc"] == "2024-01-02T03:04:05.000"
    assert recipe.qc["qc_name"].tolist() == [MASK_QC_NAME]
    assert recipe.qc["qc_value"].tolist() == [1]


def test_a_multi_frame_sof_gives_one_file_with_one_extension_per_frame_in_order(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_1.fits", (FIRST_RAY,), 60000.5, seed=1))
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_2.fits", (SECOND_RAY,), 60000.6, seed=2))

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    assert [path.name for path in Path(recipe.qcDir).iterdir()] == ["OBJECT_VIS_LACOSMIC_MASK.fits"]
    with fits.open(_mask_path(recipe)) as hdus:
        assert len(hdus) == 3
        assert [hdu.header["ORIGFILE"] for hdu in hdus[1:]] == ["raw_1.fits", "raw_2.fits"]
        assert hdus[1].data[FIRST_RAY] == 1
        assert hdus[1].data[SECOND_RAY] == 0
        assert hdus[2].data[SECOND_RAY] == 1
        assert hdus[2].data[FIRST_RAY] == 0
    assert recipe.qc["qc_value"].tolist() == [2]


def test_no_file_and_no_qc_value_are_written_when_lacosmic_is_off(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path, recipeSettings={"use_lacosmic": False, "vis": {"lacosmic-sigma": 25}})
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_1.fits", (FIRST_RAY,), 60000.5, seed=1))

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    assert list(Path(recipe.qcDir).iterdir()) == []
    assert recipe.products.empty
    assert recipe.qc.empty


def test_a_frame_with_no_cosmic_rays_still_gets_an_empty_extension_and_a_zero_count(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_1.fits", (), 60000.5, seed=1))

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    with fits.open(_mask_path(recipe)) as hdus:
        assert hdus[1].data.sum() == 0
        assert hdus[1].header["NCRPIX"] == 0
    assert recipe.qc["qc_value"].tolist() == [0]


def test_extensions_follow_observation_time_whatever_the_order_the_frames_were_prepared_in(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    recipe = _recipe(log, tmp_path)
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_late.fits", (SECOND_RAY,), 60000.9, seed=2))
    recipe._prepare_single_frame(_raw_object_frame(tmp_path, "raw_early.fits", (FIRST_RAY,), 60000.1, seed=1))

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    with fits.open(_mask_path(recipe)) as hdus:
        assert [hdu.header["ORIGFILE"] for hdu in hdus[1:]] == ["raw_early.fits", "raw_late.fits"]


def test_an_extension_omits_the_observation_time_cards_when_the_frame_has_none() -> None:
    # ARRANGE
    from soxspipe.commonutils.cosmic_ray_mask import cosmic_ray_mask_hdulist

    mask = np.zeros((4, 4), dtype=np.uint8)
    mask[1, 2] = 1

    # ACT
    hdus = cosmic_ray_mask_hdulist([{"fileName": "raw.fits", "mjdObs": None, "dateObs": None, "mask": mask}])

    # ASSERT
    header = hdus[1].header
    assert hdus[1].name == "CRMASK"
    assert "MJD-OBS" not in header
    assert "DATE-OBS" not in header
    assert header["NCRPIX"] == 1


def test_only_the_pixels_cleaning_newly_flagged_are_in_the_new_mask() -> None:
    # ARRANGE
    from soxspipe.commonutils.cosmic_ray_mask import new_cosmic_ray_mask

    before = np.array([[True, False], [False, False]])
    after = np.array([[True, True], [False, False]])

    # ACT
    newMask = new_cosmic_ray_mask(before, after)

    # ASSERT
    assert newMask.dtype == np.uint8
    assert newMask.tolist() == [[0, 1], [0, 0]]
