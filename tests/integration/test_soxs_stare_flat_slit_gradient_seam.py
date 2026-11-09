"""The shipped NIR stare settings reach the real sky fitter through the real ``soxs_stare`` constructor (DY-2147).

The recipe's own constructor runs for real, with the shipped SOXS settings: it merges the NIR block, so the
``flat-slit-gradient`` key that the sky fitter reads is the one in ``soxs_default_settings.yaml``. The recipe's
``_subtract_sky`` then builds the real ``subtract_sky`` with the kept slit ranges the recipe holds, and fits a
flat-fielded synthetic order whose master flat tilts along the slit. Only the inherited boundaries (the session lookup,
the set of files, the frame preparation) and the sky fitter's own file reads (the 2D map, the plots) are replaced.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe
import soxspipe.commonutils.set_of_files as set_of_files_module
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_stare import soxs_stare
from tests.factories import instrument_header
from tests.unit.test_soxs_stare_characterization import FrameInventory, _isolate, _stub_basics
from tests.unit.test_subtract_sky_flat_slit_gradient import KEPT_RANGE, ORDER, _tilted_order

pytestmark = pytest.mark.integration

stareModule = importlib.import_module("soxspipe.recipes.soxs_stare")
subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")
DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"
SETTING_KEY = "flat-slit-gradient"
# THE FLAT TILTS BY 10 % ALONG THE SLIT, SO THE FLAT-FIELDED SKY FALLS BY ABOUT 9 % FROM THE BOTTOM TO THE TOP OF THE
# KEPT RANGE. A MODEL THAT IGNORES IT IS LEVEL ALONG THE SLIT
MIN_TILT_FOLLOWED = 0.05


def _shipped_settings(tmpPath: Path) -> dict[str, Any]:
    """The shipped SOXS settings with the workspace keys the recipe constructor reads."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    return {
        **settings,
        "workspace-root-dir": str(tmpPath),
        "save-intermediate-products": False,
    }


def _real_nir_stare_recipe(log: Any, monkeypatch: pytest.MonkeyPatch, tmpPath: Path) -> soxs_stare:
    """Run the real ``soxs_stare.__init__`` for the NIR arm against a stubbed set of files."""
    _isolate(monkeypatch, tmpPath)
    _stub_basics(monkeypatch, arm="NIR")

    class StubSetOfFiles:
        def __init__(self, **_: object) -> None:
            return None

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            return FrameInventory(), {}

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)

    def prepare_frames(self: soxs_stare, save: bool = False) -> object:
        # THE REAL PREPARE_FRAMES ENDS BY MERGING THE ARM'S SETTINGS INTO THE RECIPE'S, NOW THE ARM IS KNOWN
        self.recipeSettings = self.get_recipe_settings()
        return object()

    monkeypatch.setattr(base_recipe, "prepare_frames", prepare_frames)
    return soxs_stare(
        log=log,
        settings=_shipped_settings(tmpPath),
        inputFrames=str(tmpPath / "2024-01-02_stare.sof"),
        turnOffMP=True,
    )


def _frame_and_map(order: pd.DataFrame) -> tuple[CCDData, CCDData, pd.DataFrame]:
    """A flat-fielded object frame, its master flat and the 2D map table, with the order's columns along x."""
    mapTable = order.assign(x=order["y"], y=order["x"])
    shape = (int(mapTable["y"].max()) + 1, int(mapTable["x"].max()) + 1)
    data, error, flat = np.zeros(shape), np.ones(shape), np.ones(shape)
    rows, columns = mapTable["y"].to_numpy(), mapTable["x"].to_numpy()
    data[rows, columns], error[rows, columns], flat[rows, columns] = mapTable[["flux", "error", "flat"]].to_numpy().T
    frame = CCDData(
        data,
        unit=u.electron,
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(error),
        meta=instrument_header(arm="NIR", overrides={"ESO INS NISE NAME": "SLIT1.0"}),
    )
    return frame, CCDData(flat, unit=u.electron), mapTable


def _model_tilt(skyModel: CCDData, mapTable: pd.DataFrame) -> float:
    """The sky model of the top of the kept range over that of the bottom, as the mean over the columns."""
    model = np.asarray(skyModel.data, dtype=float)
    slit = mapTable["slit_position"].to_numpy()
    values = model[mapTable["y"].to_numpy(), mapTable["x"].to_numpy()]
    bottom = values[(slit > KEPT_RANGE[0]) & (slit < KEPT_RANGE[0] + 1)].mean()
    top = values[(slit < KEPT_RANGE[1]) & (slit > KEPT_RANGE[1] - 1)].mean()
    return float(top / bottom)


def _subtract(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmpPath: Path, isOn: bool | None
) -> tuple[soxs_stare, CCDData, pd.DataFrame]:
    """Build the real recipe, set the key (None leaves the shipped value), and run its ``_subtract_sky``."""
    recipe = _real_nir_stare_recipe(log, monkeypatch, tmpPath)
    if isOn is not None:
        recipe.recipeSettings["sky-subtraction"][SETTING_KEY] = isOn
    frame, flat, mapTable = _frame_and_map(_tilted_order())
    monkeypatch.setattr(
        subtractSkyModule,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (mapTable.copy(), np.zeros(frame.shape, bool)),
    )
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits")
    monkeypatch.setattr(stareModule, "generic_quality_checks", lambda **kwargs: kwargs["qcTable"])
    monkeypatch.setattr(stareModule, "spectroscopic_image_quality_checks", lambda **kwargs: kwargs["qcTable"])
    recipe.subtractSky = True
    recipe.qc = pd.DataFrame()
    recipe.products = pd.DataFrame()
    recipe.readNoise = None
    recipe.debug = False
    recipe._write_sky_products = lambda *_: "product.fits"
    recipe.slitIlluminationFitRanges = {ORDER: KEPT_RANGE}
    skyModel, _, _ = recipe._subtract_sky(frame, "map.fits", "disp.fits", "orders.fits", None, flat)
    return recipe, skyModel, mapTable


def test_the_shipped_nir_setting_reaches_the_real_sky_fit_and_the_model_follows_the_flat_tilt(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ACT
    recipe, skyModel, mapTable = _subtract(log, monkeypatch, tmp_path, None)

    # ASSERT
    assert recipe.recipeSettings["sky-subtraction"][SETTING_KEY] is True
    assert _model_tilt(skyModel, mapTable) < 1 - MIN_TILT_FOLLOWED


def test_the_real_sky_fit_leaves_the_model_level_along_the_slit_with_the_setting_off(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ACT
    _, skyModel, mapTable = _subtract(log, monkeypatch, tmp_path, False)

    # ASSERT
    assert _model_tilt(skyModel, mapTable) == pytest.approx(1.0, abs=0.02)
