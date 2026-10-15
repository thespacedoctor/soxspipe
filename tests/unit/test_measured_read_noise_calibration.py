"""A detrended VIS frame's ERRS is calibrated when it uses the measured read noise (DY-1458).

The raw frame carries read noise of 4.0 e- while its header nominally says 3.3 e-. The master bias is the mean of
four frames, so its ERRS is 4.0 / sqrt(4) = 2.0 e-. After preparation and detrending, the residual over the truth,
compared with ERRS, has unit scatter only when the preparation used the measured 4.0 e-.

The scatter is the pooled ratio sqrt(mean(residual^2) / mean(ERRS^2)). ERRS is built from each pixel's own counts, so
a per-pixel (residual / ERRS) scatter reads about 3% high at this sky level whatever the read noise (Jensen's
inequality on the noisy variance estimate); the pooled ratio is free of that and is exactly 1 for a calibrated ERRS.
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
PEDESTAL = 650.0
SKY = 20.0
NOMINAL_RON = 3.3
MEASURED_RON = 4.0
N_COMBINED = 4
Z_SCATTER_TOLERANCE = 0.02
# THE SCATTER AN UNDERSTATED READ NOISE LEAVES: sqrt((20 + 16 + 4) / (20 + 3.3^2 + 3.3^2 / 4)) = 1.07
UNDERSTATED_SCATTER_FLOOR = 1.05


def _recipe(log: Any, ron: float) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: keyword
    recipe.detectorParams = {"ron": ron * u.electron}
    recipe.recipeSettings = {"subtract_background": False}
    return recipe


def _raw_frame(rng: np.random.Generator) -> CCDData:
    """A gain-corrected raw frame: pedestal, Poisson sky and read noise of the measured size."""
    data = PEDESTAL + rng.poisson(SKY, SHAPE) + rng.normal(0.0, MEASURED_RON, SHAPE)
    return CCDData(data.astype(np.float32), unit=u.electron, meta=fits.Header({"DPR_TYPE": "OBJECT"}))


def _master_bias(rng: np.random.Generator) -> CCDData:
    """The mean of four bias frames: pedestal plus noise of RON/sqrt(N), with that ERRS."""
    error = MEASURED_RON / np.sqrt(N_COMBINED)
    data = PEDESTAL + rng.normal(0.0, error, SHAPE)
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full(SHAPE, error, dtype=np.float32), unit=u.electron),
        mask=np.zeros(SHAPE, dtype=bool),
    )


def _z_scatter(recipe: base_recipe, seed: int) -> float:
    rng = np.random.default_rng(seed)
    prepared = recipe._add_uncertainty_map(_raw_frame(rng))
    prepared.mask = np.zeros(SHAPE, dtype=bool)

    detrended = recipe.detrend(prepared, master_bias=_master_bias(rng))

    residual = detrended.data.astype(np.float64) - SKY
    variance = detrended.uncertainty.array.astype(np.float64) ** 2
    return float(np.sqrt(np.mean(residual**2) / np.mean(variance)))


def test_the_measured_read_noise_gives_unit_scatter_in_the_detrended_residuals(log: Any) -> None:
    # ARRANGE
    header = fits.Header({"HIERARCH ESO QC RAW RON": MEASURED_RON, "HIERARCH ESO PRO DATANCOM": N_COMBINED})
    kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    readNoise = resolve_read_noise(header, NOMINAL_RON, 3.8, kw, log)
    recipe = _recipe(log, readNoise.singleFrame)

    # ACT
    scatter = _z_scatter(recipe, seed=5)

    # ASSERT
    assert scatter == pytest.approx(1.0, abs=Z_SCATTER_TOLERANCE)


def test_the_nominal_read_noise_understates_the_detrended_error(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log, NOMINAL_RON)

    # ACT
    scatter = _z_scatter(recipe, seed=5)

    # ASSERT
    assert scatter > UNDERSTATED_SCATTER_FLOOR
