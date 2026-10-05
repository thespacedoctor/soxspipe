"""Characterization of the subtract_sky constructor branches (DY-40)."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")

MAP_TABLE = pd.DataFrame(
    {
        "order": [11, 12, 13, 14],
        "wavelength": [1500.0, 1600.0, 1700.0, 1800.0],
        "slit_position": [0.0, 0.1, 0.2, 0.3],
    }
)


@pytest.fixture
def collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> dict[str, list]:
    """Record the file-reading and plotting collaborators; the keyword and detector lookups stay real."""
    calls: dict[str, list] = {"map": [], "quicklook": [], "filenamer": [], "setup": []}

    def unpack_map(**kwargs: Any) -> tuple[pd.DataFrame, np.ndarray]:
        calls["map"].append(kwargs)
        return MAP_TABLE.copy(), np.zeros((3, 3), dtype=bool)

    def name_frame(**kwargs: Any) -> str:
        calls["filenamer"].append(kwargs)
        return "SYNTHETIC_FROM_HEADER.fits"

    def setup_directories(**kwargs: Any) -> tuple[str, str]:
        calls["setup"].append(kwargs)
        return str(tmp_path / "init_qc"), str(tmp_path / "init_products")

    monkeypatch.setattr(subtractSkyModule, "twoD_disp_map_image_to_dataframe", unpack_map)
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **kwargs: calls["quicklook"].append(kwargs))
    monkeypatch.setattr(subtractSkyModule, "filenamer", name_frame)
    monkeypatch.setattr(toolkit, "utility_setup", setup_directories)
    return calls


def _construct(log: Any, frame: CCDData, **kwargs: Any) -> subtract_sky:
    return subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={"sky-subtraction": {"bspline_order": 3}},
        objectFrame=frame,
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        dispMap="dispersion-map.fits",
        **kwargs,
    )


def test_a_nir_jh_frame_keeps_orders_above_12_swaps_axes_and_names_from_the_header(
    log: Any, collaborators: dict[str, list]
) -> None:
    """The JH blocking filter drops orders 12 and below; NIR disperses along y and is never binned."""
    header = instrument_header(
        arm="NIR",
        overrides={"ESO INS NISE NAME": "JH0.9", "ESO DET BINX": 2, "ESO DET BINY": 2},
    )
    frame = CCDData(np.ones((3, 3)), unit=u.electron, meta=header)

    subtractor = _construct(log, frame, recipeName="soxs-nod", startNightDate="2024-01-02")

    assert subtractor.mapDF["order"].tolist() == [13, 14]
    assert subtractor.slit == "JH0.9"
    assert (subtractor.axisA, subtractor.axisB) == ("y", "x")
    assert subtractor.filenameTemplate == "SYNTHETIC_FROM_HEADER.fits"
    assert collaborators["filenamer"] == [{"log": log, "frame": frame, "settings": {"instrument": "soxs"}}]
    # NIR IGNORES THE BINNING KEYWORDS
    assert (subtractor.binx, subtractor.biny) == (1, 1)
    assert subtractor.dateObs == "2024-01-02T03:04:05.678"
    assert subtractor.stopSubtraction is False
    [mapCall] = collaborators["map"]
    assert mapCall["slit_length"] == 12
    assert mapCall["dispAxis"] == "y"
    assert mapCall["twoDMapPath"] == "two-d-map.fits"
    assert mapCall["associatedFrame"] is frame
    [quicklookCall] = collaborators["quicklook"]
    assert {key: value for key, value in quicklookCall.items() if key not in ("log", "CCDObject", "settings")} == {
        "show": False,
        "ext": False,
        "stdWindow": 0.1,
        "title": "science frame awaiting sky-subtraction",
        "surfacePlot": False,
        "dispMap": "dispersion-map.fits",
        "dispMapImage": "two-d-map.fits",
        "skylines": True,
    }
    assert collaborators["setup"] == [
        {
            "log": log,
            "settings": {"instrument": "soxs"},
            "recipeName": "soxs-nod",
            "startNightDate": "2024-01-02",
        }
    ]


def test_a_vis_frame_without_binning_keywords_defaults_to_unbinned(
    log: Any, collaborators: dict[str, list]
) -> None:
    """Without `ESO DET BINX` the binning falls back to 1x1 and every order is kept."""
    header = instrument_header(arm="VIS", overrides={"ESO INS VISE NAME": "SLIT1.0"})
    frame = CCDData(np.ones((3, 3)), unit=u.electron, meta=header)

    subtractor = _construct(log, frame, sofName="SYNTHETIC_SOF")

    assert subtractor.mapDF["order"].tolist() == [11, 12, 13, 14]
    assert (subtractor.axisA, subtractor.axisB) == ("x", "y")
    assert (subtractor.binx, subtractor.biny) == (1, 1)
    assert subtractor.filenameTemplate == "SYNTHETIC_SOF.fits"
    assert collaborators["filenamer"] == []
    assert collaborators["map"][0]["dispAxis"] == "x"
    assert collaborators["setup"][0]["recipeName"] == "soxs-stare"
