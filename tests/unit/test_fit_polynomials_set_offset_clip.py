"""The dispersion-map fit drops multi-pinhole sets that sit off the model in every pinhole (DY-1526)."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from soxspipe.commonutils.pinhole_set_offsets import OFFSET_SET_CLIP_SIGMA

pytestmark = pytest.mark.unit

PINHOLES = 9
NOISE_OFFSETS = np.linspace(-0.3, 0.3, 11)
OFFSET_WAVELENGTH = 511.0
# THE SOXS DEFAULTS: A 7-SIGMA PER-PINHOLE CLIP THAT A CONSISTENT +2 PX SET SURVIVES
SOXS_CLIPPING_SIGMA = 7.0
WELL_FITTING_LIMIT = 0.5
LINE_LIST_COLUMNS = [
    "wavelength",
    "order",
    "slit_index",
    "slit_position",
    "detector_x",
    "detector_y",
    "observed_x",
    "observed_y",
    "x_diff",
    "y_diff",
    "fit_x",
    "fit_y",
    "residuals_x",
    "residuals_y",
    "residuals_xy",
    "sigma_clipped",
    "sharpness",
    "roundness1",
    "roundness2",
    "npix",
    "sky",
    "peak",
    "flux",
    "fwhm_pin_px",
    "R_pin",
    "pixelScaleNm",
    "detector_x_shifted",
    "detector_y_shifted",
    "R_slit",
    "fwhm_slit_px",
]


def _mapper(log: object, **settings: object) -> create_dispersion_map:
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.recipeName = "soxs-spat-solution"
    mapper.settings = {}
    mapper.detectorParams = {}
    mapper.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": SOXS_CLIPPING_SIGMA,
        "poly-clipping-iteration-limit": 3,
        "poly-clipping-pinhole-sets": False,
        **settings,
    }
    return mapper


def _multi_pinhole_lines(offsetsY) -> pd.DataFrame:
    rows = []
    for setIndex, offsetY in enumerate(offsetsY):
        for pinhole in range(PINHOLES):
            rows.append(
                {
                    "order": 10,
                    "wavelength": 500.0 + setIndex,
                    "slit_index": pinhole,
                    "slit_position": pinhole - PINHOLES // 2,
                    "observed_x": 100.0 + setIndex,
                    "observed_y": 20.0 + pinhole,
                    "R_pin": 4000.0,
                    "dropped": False,
                    "offset_y": float(offsetY),
                }
            )
    return pd.DataFrame(rows)


def _stub_fit(
    mapper: create_dispersion_map,
    monkeypatch: pytest.MonkeyPatch,
    tightenAfterFirstFit: bool = False,
) -> None:
    """Fix the fit at zero and give each set its own residual offset, with a small within-set pattern.

    With ``tightenAfterFirstFit`` the sets that fit well (offset under 0.5 px) halve their offsets from the second
    residual calculation on, as a refit without the outliers does.
    """
    polynomialsModule = importlib.import_module("soxspipe.commonutils.polynomials")
    monkeypatch.setattr(
        polynomialsModule,
        "solve_linear_poly_coefficients",
        lambda **kwargs: np.zeros(kwargs["nCoeff"]),
    )
    pattern = np.linspace(-0.1, 0.1, PINHOLES)
    calls = []

    def residuals(**kwargs: object) -> tuple[float, float, float, pd.DataFrame]:
        calls.append(1)
        table = kwargs["orderPixelTable"].copy()
        offsets = table["offset_y"]
        if tightenAfterFirstFit and len(calls) > 1:
            offsets = offsets.where(offsets.abs() >= WELL_FITTING_LIMIT, offsets / 2)
        table["residuals_x"] = pattern[table["slit_index"].to_numpy()]
        table["residuals_y"] = offsets + pattern[::-1][table["slit_index"].to_numpy()]
        table["residuals_xy"] = np.hypot(table["residuals_x"], table["residuals_y"])
        return 0.0, 0.0, 0.0, table

    mapper.calculate_residuals = residuals


@pytest.mark.parametrize("clipOnSets", [False, True])
def test_a_set_offset_in_every_pinhole_is_dropped_whichever_clip_branch_runs(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    clipOnSets: bool,
) -> None:
    # ARRANGE
    mapper = _mapper(log, **{"poly-clipping-pinhole-sets": clipOnSets})
    _stub_fit(mapper, monkeypatch)
    lines = _multi_pinhole_lines([*NOISE_OFFSETS, 2.0])

    # ACT
    _, _, retained, clipped = mapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    assert OFFSET_WAVELENGTH not in set(retained["wavelength"])
    offsetRows = clipped.loc[clipped["wavelength"] == OFFSET_WAVELENGTH]
    assert len(offsetRows) == PINHOLES
    assert offsetRows["set_offset_clipped"].all()
    assert not retained["set_offset_clipped"].any()


def test_the_existing_seven_sigma_clip_alone_keeps_the_offset_set(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: A CUT FAR ABOVE THE OFFSET SWITCHES THE NEW TEST OFF
    mapper = _mapper(log, **{"poly-clipping-set-offset-sigma": 100})
    _stub_fit(mapper, monkeypatch)
    lines = _multi_pinhole_lines([*NOISE_OFFSETS, 2.0])

    # ACT
    _, _, retained, clipped = mapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    assert (retained["wavelength"] == OFFSET_WAVELENGTH).sum() == PINHOLES
    assert len(clipped) == 0


def test_the_setting_moves_the_cut(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: A +1 PX SET IS ABOUT 3.7 MAD-SIGMA OUT, INSIDE THE DEFAULT CUT OF 5
    lines = _multi_pinhole_lines([*NOISE_OFFSETS, 1.0])
    defaultMapper = _mapper(log)
    tightMapper = _mapper(log, **{"poly-clipping-set-offset-sigma": 3})
    _stub_fit(defaultMapper, monkeypatch)
    _stub_fit(tightMapper, monkeypatch)

    # ACT
    _, _, retainedDefault, _ = defaultMapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)
    _, _, retainedTight, _ = tightMapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    assert OFFSET_WAVELENGTH in set(retainedDefault["wavelength"])
    assert OFFSET_WAVELENGTH not in set(retainedTight["wavelength"])


def test_a_set_that_only_crosses_the_cut_after_a_refit_is_kept(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: THE +2 PX SET GOES ON THE FIRST FIT, WHICH STARTS A SECOND ITERATION. THE REFIT TIGHTENS THE SCATTER OF
    # THE REST, SO THE +1 PX SET (ABOUT 3.7 MAD-SIGMA ON THE FIRST FIT) WOULD LATER READ AS AN OUTLIER
    mapper = _mapper(log)
    _stub_fit(mapper, monkeypatch, tightenAfterFirstFit=True)
    lines = _multi_pinhole_lines([*NOISE_OFFSETS, 1.0, 2.0])

    # ACT
    _, _, retained, clipped = mapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    assert set(clipped["wavelength"]) == {OFFSET_WAVELENGTH + 1}
    assert OFFSET_WAVELENGTH in set(retained["wavelength"])


def test_the_dropped_set_count_is_logged(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    mapper = _mapper(log)
    _stub_fit(mapper, monkeypatch)
    lines = _multi_pinhole_lines([*NOISE_OFFSETS, 2.0])

    # ACT
    mapper.fit_polynomials(lines, wavelengthDeg=0, orderDeg=0, slitDeg=0)

    # ASSERT
    infoMessages = [message for level, message in log.messages if level == "info"]
    assert any("1 multi-pinhole set" in message and "order 10: 1" in message for message in infoMessages)


def test_fitted_lines_output_carries_the_set_offset_flag(log: object) -> None:
    # ARRANGE: A LINE DROPPED BEFORE THE FIT HAS NO FLAG COLUMN AT ALL
    mapper = _mapper(log)
    good = pd.DataFrame([dict.fromkeys(LINE_LIST_COLUMNS, 0.0)]).assign(wavelength=500.0, set_offset_clipped=False)
    clipped = pd.DataFrame([dict.fromkeys(LINE_LIST_COLUMNS, 0.0)] * 2).assign(
        wavelength=[501.0, 502.0], set_offset_clipped=[True, np.nan]
    )

    # ACT
    combined, goodOut = mapper._prepare_line_list_columns(good, clipped)

    # ASSERT
    flags = dict(zip(combined["wavelength"], combined["set_offset_clipped"], strict=True))
    assert flags == {500.0: False, 501.0: True, 502.0: False}
    assert combined["set_offset_clipped"].dtype == bool
    assert goodOut["set_offset_clipped"].tolist() == [False]


def test_lines_dropped_before_the_fit_are_written_as_not_set_offset_clipped(log: object) -> None:
    # ARRANGE: NEITHER TABLE HAS THE FLAG COLUMN, AS WHEN EVERY LINE WAS DROPPED BEFORE THE FIT
    mapper = _mapper(log)
    good = pd.DataFrame([dict.fromkeys(LINE_LIST_COLUMNS, 0.0)]).assign(sigma_clipped=False)
    clipped = pd.DataFrame([dict.fromkeys(LINE_LIST_COLUMNS, 0.0)]).assign(sigma_clipped=True)

    # ACT
    combined, goodOut = mapper._prepare_line_list_columns(good, clipped)

    # ASSERT
    assert combined["set_offset_clipped"].tolist() == [False, False]
    assert goodOut["set_offset_clipped"].tolist() == [False]


@pytest.mark.parametrize(
    ("settingsFile", "arms"),
    [("soxs_default_settings.yaml", ["vis", "nir"]), ("xsh_default_settings.yaml", ["uvb", "vis", "nir"])],
)
def test_default_settings_set_the_cut_to_five_mad_sigma(settingsFile: str, arms: list[str]) -> None:
    # ARRANGE
    settings = yaml.safe_load((Path(soxspipe.__file__).parent / settingsFile).read_text())

    # ACT
    cuts = {arm: settings["soxs-spat-solution"][arm]["poly-clipping-set-offset-sigma"] for arm in arms}

    # ASSERT
    assert cuts == dict.fromkeys(arms, OFFSET_SET_CLIP_SIGMA)
