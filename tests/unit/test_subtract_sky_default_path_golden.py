"""On the default pre-fit clip, ``subtract_sky.subtract`` gives the same products as before the DY-2055 prototype.

The golden file ``data/subtract_sky_default_golden_a151d6c5.npz`` holds the sky model, the sky-subtracted frame, its
``ERRS`` and its mask for each arm, as develop a151d6c5 (the commit before the prototype) wrote them for the synthetic
frame built here. The frame has a star, sky lines and cosmic rays, so the rolling clip, the object masking, the DY-1972
skyline restore, the object-wing correction and (NIR) the DY-1973 error inflation all act on it.

The file pins the pre-prototype behaviour. To re-record it, run ``record_golden(path, log, workDir)`` with ``soxspipe``
imported from a checkout of the commit that is the new reference and this module imported from its file path. Re-record
it only on purpose: a failure here means the default sky subtraction changed.
"""

from __future__ import annotations

import copy
import importlib
from contextlib import ExitStack
from pathlib import Path
from typing import Any
from unittest import mock

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

GOLDEN_FILE = Path(__file__).parent / "data" / "subtract_sky_default_golden_a151d6c5.npz"
ARMS = ("NIR", "VIS")
PRODUCTS = ("model", "subtracted", "errs", "mask")
STAR_ORDER = 12
SKY_ORDER = 13
COLUMNS = 400
SLIT_PIXELS = 24
SLIT_HALF_LENGTH = 5.5
SKY_LEVEL = 20.0
READ_NOISE = 3.3
STAR_CENTRE = 0.4
STAR_FWHM = 1.0
STAR_FLUX_PER_COLUMN = 3000.0
MOFFAT_BETA = 4.765
# (CENTRE COLUMN, PEAK IN E-) OF THE GAUSSIAN SKY LINES, AND THEIR SIGMA IN COLUMNS
SKY_LINES = ((60.0, 400.0), (150.5, 1500.0), (240.0, 120.0), (330.2, 800.0))
SKY_LINE_SIGMA_COLUMNS = 1.5
COSMIC_COUNT = 25
COSMIC_FLUX = 600.0
FRAME_SHAPE = (2 * SLIT_PIXELS, COLUMNS)
# THE GOLDEN FILE STORES FLOAT32 TO STAY SMALL. THE TOLERANCE IS FAR BELOW THE NOISE (ABOUT 5 E-) AND WIDE ENOUGH FOR
# THAT ROUNDING AND FOR FLOATING-POINT DIFFERENCES OF THE LINEAR ALGEBRA BETWEEN MACHINES
ABSOLUTE_TOLERANCE = 1e-4
RELATIVE_TOLERANCE = 1e-6

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"
subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")


def _sky_settings(arm: str) -> dict[str, Any]:
    """The shipped SOXS stare sky-subtraction settings of the arm, with the QC plot turned off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare"][arm.lower()]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    return skySettings


def _order(order: int, starFluxPerColumn: float, firstRow: int, seed: int) -> pd.DataFrame:
    """One synthetic order: a flat sky, Gaussian sky lines, a Moffat star along the slit and noise."""
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitGrid = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)
    slitPosition = slitGrid[slitPixel]
    # THE SLIT IS TILTED A LITTLE AGAINST THE WAVELENGTH AXIS, SO NEIGHBOURING COLUMNS INTERLEAVE IN WAVELENGTH
    tiltedColumn = column + 0.2 * slitPosition
    sky = np.full(column.size, SKY_LEVEL)
    for centre, peak in SKY_LINES:
        sky += peak * np.exp(-0.5 * ((tiltedColumn - centre) / SKY_LINE_SIGMA_COLUMNS) ** 2)
    alpha = STAR_FWHM / (2.0 * np.sqrt(2.0 ** (1.0 / MOFFAT_BETA) - 1.0))
    profile = (1.0 + ((slitGrid - STAR_CENTRE) / alpha) ** 2) ** -MOFFAT_BETA
    star = starFluxPerColumn * (profile / profile.sum())[slitPixel]
    noiseless = sky + star
    error = np.sqrt(noiseless + READ_NOISE**2)
    flux = noiseless + rng.standard_normal(noiseless.size) * error
    cosmics = rng.choice(flux.size, size=COSMIC_COUNT, replace=False)
    flux[cosmics] += COSMIC_FLUX
    return pd.DataFrame(
        {
            # THE 2D MAP STORES THE ORDER NUMBER AS A FLOAT
            "order": float(order),
            "x": column,
            "y": slitPixel + firstRow,
            "wavelength": 1000.0 + 0.02 * tiltedColumn + 20.0 * (order - STAR_ORDER),
            "slit_position": slitPosition,
            "flux": flux,
            "error": error,
            "mask": False,
        }
    )


def _map_table() -> pd.DataFrame:
    star = _order(STAR_ORDER, STAR_FLUX_PER_COLUMN, firstRow=0, seed=2055)
    sky = _order(SKY_ORDER, 0.0, firstRow=SLIT_PIXELS, seed=2056)
    return pd.concat([star, sky], ignore_index=True)


def _object_frame(mapTable: pd.DataFrame, arm: str) -> CCDData:
    data = np.zeros(FRAME_SHAPE)
    data[mapTable["y"], mapTable["x"]] = mapTable["flux"]
    error = np.full(FRAME_SHAPE, 5.0)
    error[mapTable["y"], mapTable["x"]] = mapTable["error"]
    slitKeyword = "ESO INS NISE NAME" if arm == "NIR" else "ESO INS VISE NAME"
    return CCDData(
        data,
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(error),
        meta=instrument_header(arm=arm, overrides={slitKeyword: "SLIT1.0"}),
    )


def run_default_subtract(log: Any, arm: str, workDir: Path) -> dict[str, np.ndarray]:
    """Run ``subtract`` on the synthetic frame with the shipped settings of the arm; return its products."""
    mapTable = _map_table()
    with ExitStack() as stack:
        stack.enter_context(
            mock.patch.object(
                subtractSkyModule,
                "twoD_disp_map_image_to_dataframe",
                lambda **_: (mapTable.copy(), np.zeros(FRAME_SHAPE, dtype=bool)),
            )
        )
        stack.enter_context(mock.patch.object(subtractSkyModule, "quicklook_image", lambda **_: None))
        stack.enter_context(mock.patch.object(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits"))
        stack.enter_context(
            mock.patch.object(toolkit, "utility_setup", lambda **_: (str(workDir / "qc"), str(workDir / "products")))
        )
        subtractor = subtract_sky(
            log=log,
            settings={"instrument": "soxs"},
            recipeSettings={"sky-subtraction": _sky_settings(arm)},
            objectFrame=_object_frame(mapTable, arm),
            twoDMap="two-d-map.fits",
            qcTable=pd.DataFrame(),
            productsTable=pd.DataFrame(),
        )
        model, subtracted, _, _, _ = subtractor.subtract()
    return {
        "model": np.asarray(model.data, dtype=float),
        "subtracted": np.asarray(subtracted.data, dtype=float),
        "errs": np.asarray(subtracted.uncertainty.array, dtype=float),
        "mask": np.asarray(subtracted.mask, dtype=bool),
    }


def record_golden(path: Path, log: Any, workDir: Path) -> None:
    """Write the golden products of both arms to ``path``; run it only from the reference commit."""
    arrays = {}
    for arm in ARMS:
        for name, values in run_default_subtract(log, arm, workDir / arm).items():
            arrays[f"{arm}_{name}"] = values if values.dtype == bool else values.astype(np.float32)
    np.savez_compressed(path, **arrays)


@pytest.fixture(scope="module")
def golden() -> dict[str, np.ndarray]:
    with np.load(GOLDEN_FILE) as stored:
        return {key: stored[key] for key in stored.files}


@pytest.mark.parametrize("arm", ARMS)
def test_the_default_subtract_matches_the_pre_prototype_products(
    log: Any, tmp_path: Path, golden: dict[str, np.ndarray], arm: str
) -> None:
    # ARRANGE
    expected = {name: golden[f"{arm}_{name}"] for name in PRODUCTS}

    # ACT
    products = run_default_subtract(log, arm, tmp_path)

    # ASSERT
    np.testing.assert_array_equal(products["mask"], expected["mask"])
    for name in ("model", "subtracted", "errs"):
        np.testing.assert_allclose(
            products[name],
            expected[name],
            rtol=RELATIVE_TOLERANCE,
            atol=ABSOLUTE_TOLERANCE,
            equal_nan=True,
            err_msg=name,
        )


def test_the_golden_frame_exercises_the_sky_fit(golden: dict[str, np.ndarray]) -> None:
    # A GOLDEN FILE OF EMPTY PRODUCTS WOULD PASS ANY CODE: THE MODEL MUST HOLD THE SKY LINES AND THE ROWS MUST BE FITTED
    for arm in ARMS:
        model = golden[f"{arm}_model"]
        assert np.nanmax(model) > 10 * SKY_LEVEL
        assert np.nanmedian(model) == pytest.approx(SKY_LEVEL, rel=0.1)
        assert not golden[f"{arm}_mask"].all()


def test_the_nir_errs_are_inflated_at_the_sky_lines_and_the_vis_errs_are_not(golden: dict[str, np.ndarray]) -> None:
    # THE FRAME REACHES THE DY-1973 INFLATION: NIR ERRS RISE ABOVE THE INPUT ERROR SOMEWHERE, VIS ERRS DO NOT
    mapTable = _map_table()
    inputError = np.full(FRAME_SHAPE, 5.0)
    inputError[mapTable["y"], mapTable["x"]] = mapTable["error"]
    assert (golden["NIR_errs"] > inputError * 1.01).any()
    assert not (golden["VIS_errs"] > inputError * 1.01).any()
