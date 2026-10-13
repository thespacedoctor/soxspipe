"""Object light must stay out of the sky model: the object-wing correction refits the sky without it (DY-1279).

The synthetic orders have a flat sky and a bright star with Moffat wings along the slit. The first tests drive an order
through the per-order steps of ``subtract_sky.subtract``, because the share of object light in the model is read from
the clipping flags of the fitted pixels. The last tests run ``subtract`` itself through the real ``__init__``.
"""

from __future__ import annotations

import copy
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
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")

STAR_ORDER = 12
SKY_ORDER = 13
COLUMNS = 800
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
SKY_LEVEL = 20.0
READ_NOISE = 3.3
OBJECT_CENTRE = 0.4
SEEING_FWHM = 1.0
BRIGHT_STAR_FLUX_PER_COLUMN = 10000.0
# A STAR AS BRIGHT AS THE SPECTROPHOTOMETRIC STANDARDS THE CORRECTION IS MEANT FOR
STANDARD_STAR_FLUX_PER_COLUMN = 100000.0
# THE ATMOSPHERIC-TURBULENCE MOFFAT INDEX (TRUJILLO ET AL. 2001)
TURBULENCE_BETA = 4.765
# THE FRACTION OF A STANDARD STAR'S FLUX THAT MAY GO INTO THE SKY MODEL (DY-1279 ACCEPTANCE)
MAX_OBJECT_FLUX_IN_MODEL = 0.001
# THE LIMIT FOR THE FAINTER BRIGHT STAR. THE PROFILE COMES FROM THE PIXELS THE UPPER-CLIPPED FIT USED, WHICH UNDERSTATES
# THE WINGS IN THIS NOISE-LIMITED SCENE: THE CORRECTION TAKES OUT 70 TO 85 PERCENT OF THE OBJECT LIGHT, NOT ALL OF IT.
# ON REAL STANDARDS A PROFILE FROM EVERY PIXEL TAKES OUT 20 TO 50 PERCENT TOO MUCH, SO THIS IS THE ACCEPTED TRADE
MAX_BRIGHT_STAR_FLUX_IN_MODEL = 0.002
# THE SHARE OF THE OBJECT FLUX THAT THE CORRECTION MAY TAKE OUT OF THE SKY MODEL BEYOND THE OBJECT LIGHT IN IT. A
# PROFILE BUILT FROM PIXELS THE UPPER-CLIPPED FIT DID NOT USE OVERSTATES THE WINGS AND BREAKS THIS
MAX_OVER_CORRECTION = 0.0002
# THE SYNTHETIC STAR'S WINGS PUT ABOUT 3% OF THE SKY LEVEL INTO THE FIRST MODEL; A GATE OF 0 MAKES THESE TESTS MEASURE
# THE CORRECTION ITSELF, AND THE GATE IS TESTED WITH THE WING MODEL (tests/unit/test_object_wing_model.py)
OPEN_GATE = 0.0

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"


def _sky_settings(**overrides: Any) -> dict[str, Any]:
    """The shipped SOXS VIS stare-std sky-subtraction settings, with the QC plot turned off."""
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare-std"]["vis"]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    skySettings.update(overrides)
    return skySettings


def _bare_subtractor(log: Any, **overrides: Any) -> subtract_sky:
    """A subtract_sky with only the attributes the per-order steps read."""
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = "VIS"
    subtractor.debug = False
    subtractor.binx = 1
    subtractor.biny = 1
    subtractor.axisA = "x"
    subtractor.axisB = "y"
    subtractor.ron = READ_NOISE
    subtractor.recipeSettings = {
        "sky-subtraction": _sky_settings(**{"object_wing_gate_fraction": OPEN_GATE, **overrides})
    }
    subtractor.bspline_order = subtractor.recipeSettings["sky-subtraction"]["bspline_order"]
    subtractor.qcPlotOrder = -1
    return subtractor


def _moffat_slit_profile(slitPosition: np.ndarray) -> np.ndarray:
    """A Moffat profile along the slit, normalised to sum to 1 over one wavelength column."""
    alpha = SEEING_FWHM / (2.0 * np.sqrt(2.0 ** (1.0 / TURBULENCE_BETA) - 1.0))
    profile = (1.0 + ((slitPosition - OBJECT_CENTRE) / alpha) ** 2) ** -TURBULENCE_BETA
    return profile / profile.sum()


def _order(order: int, objectFluxPerColumn: float, firstRow: int = 0, seed: int = 1279) -> pd.DataFrame:
    """One synthetic order on a flat sky; ``object_flux`` holds each pixel's noiseless object flux.

    The slit is tilted slightly against the wavelength axis, so neighbouring columns interleave in wavelength.
    """
    rng = np.random.default_rng(seed)
    column, slitPixel = np.meshgrid(np.arange(COLUMNS), np.arange(SLIT_PIXELS), indexing="ij")
    column = column.ravel()
    slitPixel = slitPixel.ravel()
    slitGrid = np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)
    slitPosition = slitGrid[slitPixel]
    objectFlux = objectFluxPerColumn * _moffat_slit_profile(slitGrid)[slitPixel]
    noiseless = SKY_LEVEL + objectFlux
    error = np.sqrt(noiseless + READ_NOISE**2)
    frame = pd.DataFrame(
        {
            # THE 2D MAP STORES THE ORDER NUMBER AS A FLOAT
            "order": float(order),
            "x": column,
            "y": slitPixel + firstRow,
            "wavelength": 400.0 + 0.01 * column + 0.002 * slitPosition + 100.0 * (order - STAR_ORDER),
            "slit_position": slitPosition,
            "flux": noiseless + rng.standard_normal(noiseless.size) * error,
            "error": error,
            "mask": False,
            "object_flux": objectFlux,
        }
    )
    return frame.sort_values("wavelength", kind="stable").reset_index(drop=True)


def _modelled(subtractor: subtract_sky, order: pd.DataFrame) -> pd.DataFrame:
    """Run the per-order steps of ``subtract`` on one order and return the modelled order."""
    settings = subtractor.recipeSettings["sky-subtraction"]
    clipped = subtractor.get_over_sampled_sky_from_order(
        order.copy(), clipBPs=True, clipSlitEdge=settings["clip-slit-edge-fraction"]
    )
    clipped = subtractor.clip_object_slit_positions([clipped], aggressive=settings["aggressive_object_masking"])
    fittedOrders, _, _ = subtractor._fit_sky_to_orders(clipped)
    return fittedOrders[0]


def _object_flux_fraction_in_model(modelled: pd.DataFrame) -> float:
    """The share of the object flux in the sky model.

    The object spectrum is flat, so the object light in the model is the mean object light left in the fitted
    pixels, at every pixel of the order. The correction removes its wing model from the fitted pixels first.
    """
    isFitted = ~modelled["flagged_all_clipped"]
    wingModel = modelled.get("object_wing_model", 0.0)
    objectLeft = (modelled["object_flux"] - wingModel)[isFitted]
    return float(objectLeft.mean() * len(modelled) / modelled["object_flux"].sum())


@pytest.mark.parametrize("seed", [1279, 1, 2, 3])
def test_with_the_correction_a_standard_star_puts_under_a_tenth_of_a_percent_of_its_flux_into_the_sky_model(
    log: Any, seed: int
) -> None:
    modelled = _modelled(_bare_subtractor(log), _order(STAR_ORDER, STANDARD_STAR_FLUX_PER_COLUMN, seed=seed))

    assert abs(_object_flux_fraction_in_model(modelled)) < MAX_OBJECT_FLUX_IN_MODEL


@pytest.mark.parametrize("seed", [1279, 1, 2, 3])
def test_with_the_correction_a_bright_star_puts_under_a_fifth_of_a_percent_of_its_flux_into_the_sky_model(
    log: Any, seed: int
) -> None:
    modelled = _modelled(_bare_subtractor(log), _order(STAR_ORDER, BRIGHT_STAR_FLUX_PER_COLUMN, seed=seed))

    assert abs(_object_flux_fraction_in_model(modelled)) < MAX_BRIGHT_STAR_FLUX_IN_MODEL


def test_without_the_correction_the_same_bright_star_puts_more_than_that_into_the_sky_model(log: Any) -> None:
    modelled = _modelled(
        _bare_subtractor(log, object_wing_correction=False), _order(STAR_ORDER, BRIGHT_STAR_FLUX_PER_COLUMN)
    )

    assert _object_flux_fraction_in_model(modelled) > MAX_BRIGHT_STAR_FLUX_IN_MODEL


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_the_correction_does_not_take_more_than_the_object_light_out_of_the_sky_model(log: Any, seed: int) -> None:
    modelled = _modelled(_bare_subtractor(log), _order(STAR_ORDER, BRIGHT_STAR_FLUX_PER_COLUMN, seed=seed))

    assert _object_flux_fraction_in_model(modelled) > -MAX_OVER_CORRECTION


def test_the_corrected_order_keeps_its_real_flux_and_subtracts_the_refitted_sky(log: Any) -> None:
    order = _order(STAR_ORDER, BRIGHT_STAR_FLUX_PER_COLUMN)

    modelled = _modelled(_bare_subtractor(log), order)

    realFlux = order.set_index(["x", "y"])["flux"]
    np.testing.assert_array_equal(modelled.set_index(["x", "y"])["flux"].reindex(realFlux.index), realFlux)
    np.testing.assert_allclose(modelled["sky_subtracted_flux"], modelled["flux"] - modelled["sky_model"])
    assert (modelled["object_wing_model"] != 0).any()


def test_a_sky_only_order_is_left_exactly_as_without_the_correction(log: Any) -> None:
    corrected = _modelled(_bare_subtractor(log), _order(SKY_ORDER, 0.0))
    uncorrected = _modelled(_bare_subtractor(log, object_wing_correction=False), _order(SKY_ORDER, 0.0))

    assert "object_wing_model" not in corrected
    pd.testing.assert_frame_equal(corrected, uncorrected)


def test_settings_without_the_object_wing_keys_use_the_shipped_defaults(log: Any) -> None:
    subtractor = _bare_subtractor(log)
    del subtractor.recipeSettings["sky-subtraction"]["object_wing_correction"]
    del subtractor.recipeSettings["sky-subtraction"]["object_wing_gate_fraction"]

    assert subtractor._object_wing_settings() == (True, 0.1)


@pytest.mark.parametrize("gateFraction", [-0.1, float("nan"), float("inf")])
def test_an_invalid_gate_fraction_raises(log: Any, gateFraction: float) -> None:
    subtractor = _bare_subtractor(log, object_wing_gate_fraction=gateFraction)

    with pytest.raises(ValueError, match="object_wing_gate_fraction"):
        subtractor._object_wing_settings()


# ----- subtract() THROUGH THE REAL __init__ -----

FRAME_SHAPE = (2 * SLIT_PIXELS, COLUMNS)


def _map_table() -> pd.DataFrame:
    star = _order(STAR_ORDER, BRIGHT_STAR_FLUX_PER_COLUMN)
    sky = _order(SKY_ORDER, 0.0, firstRow=SLIT_PIXELS, seed=1280)
    return pd.concat([star, sky], ignore_index=True)


@pytest.fixture
def map_table() -> pd.DataFrame:
    return _map_table()


@pytest.fixture
def collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, map_table: pd.DataFrame) -> None:
    """Replace the file-reading and plotting collaborators; the keyword and detector lookups stay real."""
    monkeypatch.setattr(
        subtractSkyModule,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (map_table.drop(columns="object_flux"), np.zeros(FRAME_SHAPE, dtype=bool)),
    )
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))


def _object_frame(map_table: pd.DataFrame) -> CCDData:
    data = np.zeros(FRAME_SHAPE)
    data[map_table["y"], map_table["x"]] = map_table["flux"]
    return CCDData(
        data,
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(FRAME_SHAPE, 5.0)),
        meta=instrument_header(arm="VIS", overrides={"ESO INS VISE NAME": "SLIT1.0"}),
    )


def _subtract(log: Any, map_table: pd.DataFrame, **overrides: Any) -> tuple[CCDData, CCDData]:
    subtractor = subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={"sky-subtraction": _sky_settings(**{"object_wing_gate_fraction": OPEN_GATE, **overrides})},
        objectFrame=_object_frame(map_table),
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
    )
    model, subtracted, _, _, _ = subtractor.subtract()
    return model, subtracted


@pytest.mark.usefixtures("collaborators")
def test_subtract_corrects_only_the_star_order_and_logs_it_once(log: Any, map_table: pd.DataFrame) -> None:
    modelOn, subtractedOn = _subtract(log, map_table)
    messages = [message for _, message in log.messages if "OBJECT-WING CORRECTION" in message]
    modelOff, subtractedOff = _subtract(log, map_table, object_wing_correction=False)

    starRows, skyRows = slice(0, SLIT_PIXELS), slice(SLIT_PIXELS, None)
    np.testing.assert_array_equal(modelOn.data[skyRows], modelOff.data[skyRows])
    np.testing.assert_array_equal(subtractedOn.data[skyRows], subtractedOff.data[skyRows])
    assert not np.array_equal(modelOn.data[starRows], modelOff.data[starRows])
    assert len(messages) == 1
    assert f"ORDER {STAR_ORDER}:" in messages[0]


@pytest.mark.usefixtures("collaborators")
def test_subtract_returns_the_object_frame_minus_the_corrected_sky_model(log: Any, map_table: pd.DataFrame) -> None:
    model, subtracted = _subtract(log, map_table)

    rows, columns = map_table["y"].to_numpy(), map_table["x"].to_numpy()
    np.testing.assert_allclose(
        subtracted.data[rows, columns], map_table["flux"].to_numpy() - model.data[rows, columns], rtol=1e-6, atol=1e-3
    )
