"""The spline sky fit on a nod A-B difference frame removes a residual skyline and leaves both traces alone (DY-2057).

The synthetic order is what is left of an A-B difference: no sky continuum, a residual skyline at three wavelengths
(the same at every slit position, one of them negative), a positive trace at the A position and a negative trace at the
B position. ``subtract`` runs through the real ``__init__``.
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
from soxspipe.commonutils.difference_sky_fit import difference_fit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")

DEFAULT_SETTINGS = Path(soxspipe.__file__).parent / "soxs_default_settings.yaml"
SLIT_KEYWORD = {"VIS": "ESO INS VISE NAME", "NIR": "ESO INS NISE NAME"}
ORDER = 12
COLUMNS = 1000
SLIT_PIXELS = 40
SLIT_HALF_LENGTH = 5.5
SKY_LEVEL = 150.0
READ_NOISE = 10.0
# THE NOISE OF ONE A-B DIFFERENCE: TWICE THE SKY AND READ-NOISE VARIANCE OF ONE FRAME
NOISE = float(np.sqrt(2.0 * (SKY_LEVEL + READ_NOISE**2)))
# (CENTRE COLUMN, AMPLITUDE IN NOISE UNITS); THE LAST ONE IS A SKY THAT WAS FAINTER IN A THAN IN B
SKYLINES = [(250, 8.0), (550, 5.0), (800, -6.0)]
LINE_WIDTH_COLUMNS = 2.5
POSITIVE_CENTRE, POSITIVE_PEAK = 1.5, 6.0
NEGATIVE_CENTRE, NEGATIVE_PEAK = -2.0, 4.0
TRACE_SIGMA = 0.4
FRAME_SHAPE = (SLIT_PIXELS, COLUMNS)
# A TRACE PIXEL IS WITHIN THIS MANY TRACE SIGMA OF ITS CENTRE
TRACE_CORE_SIGMA = 1.5
# THE FIT MAY LEAVE THIS FRACTION OF A RESIDUAL SKYLINE AT ITS CORE
MAX_LINE_LEFTOVER = 0.15
# A COLUMN IS ON A SKYLINE CORE WHEN THE RESIDUAL IS ABOVE THIS FRACTION OF THE LINE'S PEAK, AND LINE-FREE BELOW
# LINE_FREE_RESIDUAL NOISE UNITS
LINE_CORE_FRACTION = 0.8
LINE_FREE_RESIDUAL = 0.01
# THE MEAN OF A TRACE'S SKY-SUBTRACTED PIXELS MAY DIFFER FROM THE TRUE TRACE BY THIS MANY NOISE UNITS
MAX_TRACE_BIAS_SIGMA = 0.12


def _slit_grid() -> np.ndarray:
    return np.linspace(-SLIT_HALF_LENGTH, SLIT_HALF_LENGTH, SLIT_PIXELS)


def _residual_sky(column: np.ndarray) -> np.ndarray:
    sky = np.zeros(column.shape)
    for centre, amplitude in SKYLINES:
        sky = sky + amplitude * NOISE * np.exp(-0.5 * ((column - centre) / LINE_WIDTH_COLUMNS) ** 2)
    return sky


def _trace_flux(slitPosition: np.ndarray, isBothSigns: bool = True) -> np.ndarray:
    trace = POSITIVE_PEAK * NOISE * np.exp(-0.5 * ((slitPosition - POSITIVE_CENTRE) / TRACE_SIGMA) ** 2)
    if isBothSigns:
        trace = trace - NEGATIVE_PEAK * NOISE * np.exp(-0.5 * ((slitPosition - NEGATIVE_CENTRE) / TRACE_SIGMA) ** 2)
    return trace


def _difference_order(seed: int = 2057, hasSkyline: bool = True) -> pd.DataFrame:
    """One order of an A-B difference. ``residual`` and ``trace`` hold the noiseless parts of ``flux``."""
    rng = np.random.default_rng(seed)
    row, column = np.meshgrid(np.arange(SLIT_PIXELS), np.arange(COLUMNS), indexing="ij")
    row, column = row.ravel(), column.ravel()
    slitPosition = _slit_grid()[row]
    residual = _residual_sky(column) if hasSkyline else np.zeros(column.shape)
    trace = _trace_flux(slitPosition)
    table = pd.DataFrame(
        {
            "order": float(ORDER),
            "x": column,
            "y": row,
            "wavelength": 1500.0 + 0.01 * column + 0.002 * slitPosition,
            "slit_position": slitPosition,
            "flux": residual + trace + rng.normal(0.0, NOISE, residual.size),
            "error": NOISE,
            "mask": False,
            "residual": residual,
            "trace": trace,
        }
    )
    return table.sort_values("wavelength", kind="stable").reset_index(drop=True)


def _sky_settings(arm: str) -> dict[str, Any]:
    with DEFAULT_SETTINGS.open() as stream:
        settings = yaml.safe_load(stream)
    skySettings = copy.deepcopy(settings["soxs-stare"][arm.lower()]["sky-subtraction"])
    skySettings["sky_model_qc_plot"] = False
    return skySettings


@pytest.fixture
def map_table() -> pd.DataFrame:
    return _difference_order()


def _patch_collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, table: pd.DataFrame) -> None:
    """Replace the file-reading and plotting collaborators; the keyword and detector lookups stay real."""
    monkeypatch.setattr(
        subtractSkyModule,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (table.drop(columns=["residual", "trace"]), np.zeros(FRAME_SHAPE, dtype=bool)),
    )
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))


@pytest.fixture
def collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, map_table: pd.DataFrame) -> None:
    _patch_collaborators(monkeypatch, tmp_path, map_table)


def _object_frame(map_table: pd.DataFrame, arm: str) -> CCDData:
    data = np.zeros(FRAME_SHAPE)
    data[map_table["y"], map_table["x"]] = map_table["flux"]
    return CCDData(
        data,
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(FRAME_SHAPE, NOISE)),
        meta=instrument_header(arm=arm, overrides={SLIT_KEYWORD[arm]: "SLIT1.0"}),
    )


def _sky_level_frame() -> CCDData:
    return CCDData(np.full(FRAME_SHAPE, SKY_LEVEL), unit=u.electron)


def _subtract(
    log: Any, map_table: pd.DataFrame, arm: str = "NIR", *, differenceFit: difference_fit | None = None
) -> tuple[CCDData, CCDData]:
    subtractor = subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={"sky-subtraction": _sky_settings(arm)},
        objectFrame=_object_frame(map_table, arm),
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        ron=READ_NOISE,
        differenceFit=differenceFit,
    )
    model, subtracted, _, _, _ = subtractor.subtract()
    return model, subtracted


def _difference_fit() -> difference_fit:
    return difference_fit(noiseSkyFrame=_sky_level_frame())


def _is_trace_core(table: pd.DataFrame, centre: float) -> np.ndarray:
    return ((table["slit_position"] - centre).abs() < TRACE_CORE_SIGMA * TRACE_SIGMA).to_numpy()


def _frame_values(frame: CCDData, table: pd.DataFrame) -> np.ndarray:
    return frame.data[table["y"].to_numpy(), table["x"].to_numpy()]


def _line_core(table: pd.DataFrame, centre: int, amplitude: float) -> np.ndarray:
    peak = abs(amplitude) * NOISE
    return (table["residual"].abs() > LINE_CORE_FRACTION * peak).to_numpy() & (
        (table["x"] - centre).abs() <= LINE_WIDTH_COLUMNS
    ).to_numpy()


@pytest.mark.usefixtures("collaborators")
@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_the_fit_removes_each_residual_skyline_from_the_pixels_off_the_traces(
    log: Any, map_table: pd.DataFrame, arm: str
) -> None:
    # ARRANGE
    isOffTrace = ~(_is_trace_core(map_table, POSITIVE_CENTRE) | _is_trace_core(map_table, NEGATIVE_CENTRE))

    # ACT
    model, _ = _subtract(log, map_table, arm, differenceFit=_difference_fit())

    # ASSERT
    fitted = _frame_values(model, map_table)
    for centre, amplitude in SKYLINES:
        isCore = _line_core(map_table, centre, amplitude) & isOffTrace
        trueMean = map_table.loc[isCore, "residual"].mean()
        assert isCore.sum() > 50
        assert abs(fitted[isCore].mean() - trueMean) < MAX_LINE_LEFTOVER * abs(trueMean)


@pytest.mark.usefixtures("collaborators")
@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_the_traces_keep_their_flux_in_the_columns_with_no_skyline(log: Any, map_table: pd.DataFrame, arm: str) -> None:
    # ARRANGE
    isLineFree = (map_table["residual"].abs() < LINE_FREE_RESIDUAL * NOISE).to_numpy()

    # ACT
    _, subtracted = _subtract(log, map_table, arm, differenceFit=_difference_fit())

    # ASSERT
    measured = _frame_values(subtracted, map_table)
    for centre in (POSITIVE_CENTRE, NEGATIVE_CENTRE):
        isTrace = _is_trace_core(map_table, centre) & isLineFree
        change = np.mean(measured[isTrace] - map_table["flux"].to_numpy()[isTrace])
        assert isTrace.sum() > 1000
        assert abs(change) < MAX_TRACE_BIAS_SIGMA * NOISE


@pytest.mark.usefixtures("collaborators")
@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_the_traces_lose_only_the_skyline_in_the_columns_with_a_skyline(
    log: Any, map_table: pd.DataFrame, arm: str
) -> None:
    # ARRANGE
    isTrace = _is_trace_core(map_table, POSITIVE_CENTRE) | _is_trace_core(map_table, NEGATIVE_CENTRE)

    # ACT
    _, subtracted = _subtract(log, map_table, arm, differenceFit=_difference_fit())

    # ASSERT
    measured = _frame_values(subtracted, map_table)
    for centre, amplitude in SKYLINES:
        isCore = _line_core(map_table, centre, amplitude) & isTrace
        leftover = np.mean(measured[isCore] - map_table["trace"].to_numpy()[isCore])
        assert isCore.sum() > 5
        assert abs(leftover) < MAX_LINE_LEFTOVER * abs(amplitude) * NOISE + 3.0 * NOISE / np.sqrt(isCore.sum())


@pytest.mark.parametrize("arm", ["NIR", "VIS"])
def test_a_difference_order_with_no_sky_residual_gets_no_sky_model(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, arm: str
) -> None:
    # ARRANGE
    quietTable = _difference_order(seed=7, hasSkyline=False)
    _patch_collaborators(monkeypatch, tmp_path, quietTable)

    # ACT
    model, subtracted = _subtract(log, quietTable, arm, differenceFit=_difference_fit())

    # ASSERT
    assert not model.data.any()
    np.testing.assert_array_equal(_frame_values(subtracted, quietTable), quietTable["flux"].to_numpy())


# ----- THE NOISE MODEL (DY-2057) -----


def _weights_subtractor(log: Any, differenceFit: difference_fit | None) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.ron = READ_NOISE
    subtractor.differenceFit = differenceFit
    return subtractor


def _noise_table(flat: float | None = None) -> pd.DataFrame:
    table = pd.DataFrame({"noise_sky": [0.0, 50.0, 400.0, -20.0]})
    if flat is not None:
        table["flat"] = flat
    return table


def test_an_ordinary_frame_is_weighted_by_its_sky_model_as_before(log: Any) -> None:
    # ARRANGE
    subtractor = _weights_subtractor(log, None)
    skyModel = np.array([0.0, 50.0, 400.0, -20.0])

    # ACT
    weights = subtractor._pixel_weights(_noise_table(), skyModel)

    # ASSERT
    np.testing.assert_array_equal(weights, subtract_sky._inverse_noise_weights(skyModel, READ_NOISE))


def test_a_difference_frame_has_twice_the_sky_and_read_noise_variance_and_no_bias_term(log: Any) -> None:
    # ARRANGE
    subtractor = _weights_subtractor(log, difference_fit(varianceFactor=2.0))
    table = _noise_table()

    # ACT
    weights = subtractor._pixel_weights(table, np.array([9e9, 9e9, 9e9, 9e9]))

    # ASSERT
    expectedVariance = 2.0 * (np.clip(table["noise_sky"].to_numpy(), 0.0, None) + READ_NOISE**2)
    np.testing.assert_allclose(weights, 1.0 / np.sqrt(expectedVariance))


def test_the_variance_factor_scales_both_noise_terms_alike(log: Any) -> None:
    # ARRANGE
    table = _noise_table()
    atTwo = _weights_subtractor(log, difference_fit(varianceFactor=2.0))
    atOne = _weights_subtractor(log, difference_fit(varianceFactor=1.0))

    # ACT
    ratio = atOne._pixel_weights(table, None) / atTwo._pixel_weights(table, None)

    # ASSERT
    np.testing.assert_allclose(ratio, np.sqrt(2.0))


def test_a_difference_frame_on_a_flat_gives_a_low_flat_pixel_a_low_weight(log: Any) -> None:
    # ARRANGE
    subtractor = _weights_subtractor(log, difference_fit())
    flat = 0.5
    table = _noise_table(flat=flat)

    # ACT
    weights = subtractor._pixel_weights(table, None)

    # ASSERT
    sky = np.clip(table["noise_sky"].to_numpy(), 0.0, None)
    expected = flat / np.sqrt(2.0 * sky * flat + 2.0 * READ_NOISE**2)
    np.testing.assert_allclose(weights, expected)


def test_a_difference_fit_with_no_noise_sky_frame_is_refused_when_the_subtractor_is_built(
    log: Any, map_table: pd.DataFrame, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    _patch_collaborators(monkeypatch, tmp_path, map_table)

    # ACT, ASSERT
    with pytest.raises(ValueError, match="noiseSkyFrame"):
        _subtract(log, map_table, "NIR", differenceFit=difference_fit())


def test_weights_for_a_difference_frame_with_no_noise_sky_column_raise_a_clear_error(log: Any) -> None:
    # ARRANGE
    subtractor = _weights_subtractor(log, difference_fit())

    # ACT, ASSERT
    with pytest.raises(ValueError, match="noise_sky"):
        subtractor._pixel_weights(pd.DataFrame({"flux": [1.0]}), None)
