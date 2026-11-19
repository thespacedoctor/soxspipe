"""Builders for a bare `base_recipe` that can prepare a synthetic 8x8 raw frame."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes import base_recipe

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


def detector_params(badPixelMapName: str) -> dict[str, Any]:
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


def raw_frame_file(tmpPath: Path, *, dprType: str = "BIAS", **headerExtras: Any) -> str:
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


def bad_pixel_map_file(calibrationPath: Path, name: str = "bpm.fits") -> str:
    """Write an all-good bad-pixel map and return its file name."""
    calibrationPath.mkdir(parents=True, exist_ok=True)
    mask = CCDData(np.zeros((8, 8), dtype=np.int16), unit=u.dimensionless_unscaled)
    mask.write(str(calibrationPath / name), overwrite=True)
    return name


def prepare_frame_recipe(
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
    recipe.detectorParams = detector_params(badPixelMapName)
    recipe.calibrationRootPath = str(tmpPath / "calibrations")
    recipe.workspaceRootPath = str(tmpPath / "workspace")
    recipe.outDir = str(tmpPath / "tmp")
    recipe.sofName = "synthetic"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = pd.DataFrame({column: [] for column in QC_COLUMNS})
    recipe.generateReponseCurve = False
    recipe.cosmicRayMasks = []
    return recipe
