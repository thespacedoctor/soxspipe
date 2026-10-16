"""A dark-subtracted NIR frame's ERRS is calibrated when it uses the measured read noise (DY-1581).

The raw frame carries read noise of 15 e- while the detector yaml says 7.3 e-. The master dark is the mean of four
frames, so its ERRS is 15 / sqrt(4) = 7.5 e-. After preparation and detrending, the residual over the truth, compared
with ERRS, has unit scatter only when the preparation used the measured 15 e-, with no 1/N term of its own: the master
dark's ERRS supplies it.

The scatter is the pooled ratio sqrt(mean(residual^2) / mean(ERRS^2)); see test_measured_read_noise_calibration.py.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.read_noise import resolve_read_noise
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

SHAPE = (400, 400)
DARK_LEVEL = 40.0
SKY = 30.0
YAML_RON = 7.3
MEASURED_RON = 15.0
N_COMBINED = 4
Z_SCATTER_TOLERANCE = 0.02
# THE SCATTER AN UNDERSTATED READ NOISE LEAVES: sqrt((70 + 225 + 56.25) / (70 + 7.3^2 + 56.25)) = 1.40
UNDERSTATED_SCATTER_FLOOR = 1.3
HEADER = {"DPR_TYPE": "OBJECT", "EXPTIME": 60.0}


def _recipe(log: Any, ron: float) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.arm = "NIR"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: keyword
    recipe.detectorParams = {"ron": ron * u.electron}
    recipe.recipeSettings = {"subtract_background": False}
    recipe.darkDetrendWarningIssued1 = False
    recipe.darkDetrendWarningIssued2 = False
    return recipe


def _raw_frame(rng: np.random.Generator) -> CCDData:
    """A gain-corrected raw frame: Poisson dark level and sky, and read noise of the measured size."""
    data = rng.poisson(DARK_LEVEL + SKY, SHAPE) + rng.normal(0.0, MEASURED_RON, SHAPE)
    return CCDData(data.astype(np.float32), unit=u.electron, meta=fits.Header(HEADER))


def _master_dark(rng: np.random.Generator) -> CCDData:
    """The mean of four dark frames: the dark level plus noise of RON/sqrt(N), with that ERRS."""
    error = MEASURED_RON / np.sqrt(N_COMBINED)
    data = DARK_LEVEL + rng.normal(0.0, error, SHAPE)
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        meta=fits.Header(HEADER),
        uncertainty=StdDevUncertainty(np.full(SHAPE, error, dtype=np.float32), unit=u.electron),
        mask=np.zeros(SHAPE, dtype=bool),
    )


def _z_scatter(recipe: base_recipe, seed: int) -> float:
    rng = np.random.default_rng(seed)
    prepared = recipe._add_uncertainty_map(_raw_frame(rng))
    prepared.mask = np.zeros(SHAPE, dtype=bool)

    detrended = recipe.detrend(prepared, dark=_master_dark(rng))

    residual = detrended.data.astype(np.float64) - SKY
    variance = detrended.uncertainty.array.astype(np.float64) ** 2
    return float(np.sqrt(np.mean(residual**2) / np.mean(variance)))


def test_the_measured_read_noise_gives_unit_scatter_in_the_dark_subtracted_residuals(log: Any) -> None:
    # ARRANGE
    header = fits.Header({"HIERARCH ESO QC RAW RON": MEASURED_RON, "HIERARCH ESO PRO DATANCOM": N_COMBINED})
    kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    readNoise = resolve_read_noise(header, None, YAML_RON, kw, log, masterLabel="master dark")
    recipe = _recipe(log, readNoise.singleFrame)

    # ACT
    scatter = _z_scatter(recipe, seed=5)

    # ASSERT
    assert scatter == pytest.approx(1.0, abs=Z_SCATTER_TOLERANCE)


def test_the_yaml_read_noise_understates_the_dark_subtracted_error(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log, YAML_RON)

    # ACT
    scatter = _z_scatter(recipe, seed=5)

    # ASSERT
    assert scatter > UNDERSTATED_SCATTER_FLOOR
