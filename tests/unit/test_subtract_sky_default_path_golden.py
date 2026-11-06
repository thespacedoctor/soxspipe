"""With the shipped stare settings, ``subtract_sky.subtract`` gives the products recorded for DY-2079 on f33c0e7b.

The golden file ``data/subtract_sky_default_golden_f33c0e7b_dy2079.npz`` holds the sky model, the sky-subtracted frame,
its ``ERRS`` and its mask for each arm, for the synthetic frame built here. It was recorded from develop cff556ae and
re-recorded for DY-2081 (the VIS arrays: the pre-fit clip no longer retries once it clips about half the pixels
Gaussian sky noise would; the fix does not change the NIR arrays, the VIS ``ERRS`` or the VIS mask) and for DY-2079
(below). The frame has two orders. Order 12 holds a star, sky lines and cosmic rays, and order 13 holds sky lines and
cosmic rays, so the rolling clip, the object masking, the DY-1972 skyline restore, the object-wing correction and
(NIR) the DY-1973 error inflation all act on it.

The frame has 96 slit pixels per column and 100 cosmic rays per order. It had 24 slit pixels and 25 cosmic rays when the
file was recorded from a151d6c5. At 24 pixels the pre-fit clip can flag every pixel of a sky-line core in the VIS arm,
and the sky model then loses the line. Which lines it loses changes with the clip settings, so the old file pinned a
lossy state (DY-2118). At 96 pixels per column every line is held.

The geometry changes every array, the NIR arrays and ``ERRS`` included: the frame is 192 rows (was 48), and the noise
draw and the cosmic rays are new. The sizes stay close. The NIR ``ERRS`` are at most 1.47 times the input error (was
1.39) and 48% of the pixels are above 1.01 times it (was 41%). The VIS ``ERRS`` equal the input error, as before. The
median sky model is 20.1 e- in NIR and 20.1 e- in VIS (was 19.9 and 20.0 e-), and the highest is 1520 e- in NIR and
1525 e- in VIS (was 1515 and 1513 e-). The golden file is 1.3 MB (was 0.33 MB).

The DY-2081 fix changes the VIS sky model and sky-subtracted frame only: by 0.58 e- at the median and 106 e- at most
(1690 of 76800 pixels move by more than 10 e-). In order 13 the mean of (model - true sky) over the input error goes
from -0.210 to -0.021, because the clip no longer removes the high tail of the noise. Order 12 moves from +0.004 to
+0.015.

No VIS order retries its pre-fit clip on this frame under the fix, so the unit tests pin the retry, not this file. The
DY-1972 skyline restore is NIR-only by default, so the VIS arm has no restore.

DY-2079 re-recorded the file on f33c0e7b (develop with the DY-2081 clip fix) because the object-wing model now measures
its profile in up to 6 wavelength chunks. Only order 12, the star order, changes, in both arms. Against the f33c0e7b
file (pixels counted where the change exceeds 1e-3, both medians signed as new minus old):

- NIR sky model and sky-subtracted frame: 37909 and 37911 of 76800 pixels change, by at most 0.28 e- (+0.23 up, -0.28
  down). The median change over the changed pixels is +0.021 e- in the model and -0.021 e- in the subtracted frame.
- NIR ``ERRS``: 23626 pixels change, by at most 2.46 e- (+2.03 up, -2.46 down). The median signed change is -0.0023 e-
  and the median size is 0.098 e-, because the DY-1973 inflation follows the residuals of the new sky model.
- VIS sky model and sky-subtracted frame: 38055 and 38054 pixels change, by at most 0.25 e- (+0.25 up, -0.20 down). The
  median signed change is -0.003 e- in the model and +0.003 e- in the subtracted frame.
- The VIS ``ERRS`` and both masks do not change. The model still holds every sky line in both orders.

The file pins the default sky subtraction. To re-record it, run ``record_golden(path, log, workDir)`` with ``soxspipe``
imported from a checkout of the commit that is the new reference and this module imported from its file path. ``log`` is
the ``log`` fixture of ``tests/conftest.py`` or any object with ``print``, ``debug`` and ``warning``. Recording from the
DY-2081 fix commit reproduces every array of the file. Re-record it only on purpose: a failure here means the default
sky subtraction changed.
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

GOLDEN_FILE = Path(__file__).parent / "data" / "subtract_sky_default_golden_f33c0e7b_dy2079.npz"
ARMS = ("NIR", "VIS")
PRODUCTS = ("model", "subtracted", "errs", "mask")
STAR_ORDER = 12
SKY_ORDER = 13
COLUMNS = 400
SLIT_PIXELS = 96
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
COSMIC_COUNT = 100
COSMIC_FLUX = 600.0
FRAME_SHAPE = (2 * SLIT_PIXELS, COLUMNS)
# ERROR OF THE FRAME PIXELS THAT NO ORDER COVERS
OFF_ORDER_ERROR = 5.0
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


def _true_sky(tiltedColumn: np.ndarray) -> np.ndarray:
    """The noiseless sky (e-) at each tilted column: the flat level plus the Gaussian sky lines."""
    sky = np.full(tiltedColumn.size, SKY_LEVEL)
    for centre, peak in SKY_LINES:
        sky += peak * np.exp(-0.5 * ((tiltedColumn - centre) / SKY_LINE_SIGMA_COLUMNS) ** 2)
    return sky


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
    sky = _true_sky(tiltedColumn)
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


def _input_error(mapTable: pd.DataFrame) -> np.ndarray:
    error = np.full(FRAME_SHAPE, OFF_ORDER_ERROR)
    error[mapTable["y"], mapTable["x"]] = mapTable["error"]
    return error


def _object_frame(mapTable: pd.DataFrame, arm: str) -> CCDData:
    data = np.zeros(FRAME_SHAPE)
    data[mapTable["y"], mapTable["x"]] = mapTable["flux"]
    error = _input_error(mapTable)
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
def test_the_default_subtract_matches_the_recorded_products(
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


def test_the_golden_model_holds_every_sky_line_in_both_orders(golden: dict[str, np.ndarray]) -> None:
    # THE PRE-FIT CLIP CAN FLAG A WHOLE LINE CORE AND THE FIT THEN LOSES THE LINE (DY-2118): THE OLD 24-PIXEL FRAME DID
    mapTable = _map_table()
    tiltedColumn = mapTable["x"].to_numpy(dtype=float) + 0.2 * mapTable["slit_position"].to_numpy(dtype=float)
    trueSky = _true_sky(tiltedColumn)
    for arm in ARMS:
        model = golden[f"{arm}_model"][mapTable["y"].to_numpy(), mapTable["x"].to_numpy()]
        for order in (STAR_ORDER, SKY_ORDER):
            for centre, _ in SKY_LINES:
                isCore = (mapTable["order"].to_numpy() == order) & (np.abs(tiltedColumn - centre) < 1.0)
                held = model[isCore].sum() / trueSky[isCore].sum()
                assert 0.95 < held < 1.05, (
                    f"{arm} order {order} line at column {centre}: model over true sky {held:.3f}"
                )


def test_the_nir_errs_are_inflated_at_the_sky_lines_and_the_vis_errs_are_not(golden: dict[str, np.ndarray]) -> None:
    # THE FRAME REACHES THE DY-1973 INFLATION: NIR ERRS RISE ABOVE THE INPUT ERROR SOMEWHERE, VIS ERRS DO NOT
    mapTable = _map_table()
    inputError = _input_error(mapTable)
    assert (golden["NIR_errs"] > inputError * 1.01).any()
    assert not (golden["VIS_errs"] > inputError * 1.01).any()
