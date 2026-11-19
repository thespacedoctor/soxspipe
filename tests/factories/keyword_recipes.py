"""Builders for a bare `base_recipe` that stamps product-header keywords."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from astropy.table import MaskedColumn, Table

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories.settings import pipeline_settings


def keyword_recipe(tmp_path: Path, log: Any, rawFrames: Table) -> base_recipe:
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


def raw_frame_table_with_calibration(tmp_path: Path, calibrationPath: Path) -> Table:
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
