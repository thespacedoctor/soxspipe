"""Initialization contracts for the shared utility base class."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe import commonutils
from soxspipe.commonutils import toolkit
from soxspipe.commonutils.base_util import base_util

pytestmark = pytest.mark.unit


def _patch_lookups(monkeypatch: pytest.MonkeyPatch, dispersionAxis: str) -> None:
    monkeypatch.setattr(
        commonutils,
        "keyword_lookup",
        lambda **kwargs: SimpleNamespace(get=lambda name: name),
    )
    monkeypatch.setattr(
        commonutils,
        "detector_lookup",
        lambda **kwargs: SimpleNamespace(
            get=lambda arm: {
                "dispersion-axis": dispersionAxis,
                "slit_length": 4,
            }
        ),
    )
    monkeypatch.setattr(
        toolkit,
        "get_skylines_dataframe",
        lambda *args, **kwargs: pd.DataFrame({"wavelength": [500.0]}),
    )


@pytest.mark.parametrize(
    ("dispersionAxis", "expectedAxes"), [("x", ("x", "y")), ("y", ("y", "x"))]
)
def test_base_util_initializes_orientation_and_default_binning(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    dispersionAxis: str,
    expectedAxes: tuple[str, str],
) -> None:
    _patch_lookups(monkeypatch, dispersionAxis)
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02"},
    )

    result = base_util(log, {}, associatedFrame=frame)

    assert (result.axisA, result.axisB) == expectedAxes
    assert (result.binx, result.biny) == (1, 1)
    assert result.skylinesDF["wavelength"].tolist() == [500.0]


def test_base_util_reads_format_and_rebins_two_d_map(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    monkeypatch.setattr(
        toolkit,
        "read_spectral_format",
        lambda **kwargs: ([10], [400.0], [500.0], [1], [2]),
    )
    mapDataframe = pd.DataFrame({"order": [10]})
    interOrderMask = np.array([[False, True], [False, False]])
    monkeypatch.setattr(
        toolkit,
        "twoD_disp_map_image_to_dataframe",
        lambda **kwargs: (mapDataframe, interOrderMask),
    )
    mapPath = tmp_path / "map.fits"
    primary = fits.PrimaryHDU(np.ones((4, 4)))
    primary.header["WIN_BINX"] = 1
    primary.header["WIN_BINY"] = 1
    hdus = [primary]
    for name, value in (("WAVELENGTH", 500.0), ("SLIT", 0.5), ("ORDER", 10.0)):
        hdu = fits.ImageHDU(np.full((4, 4), value), name=name)
        hdus.append(hdu)
    fits.HDUList(hdus).writeto(mapPath)
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2024-01-02",
            "WIN_BINX": 2,
            "WIN_BINY": 2,
        },
    )

    result = base_util(
        log,
        {},
        associatedFrame=frame,
        dispersionMap="solution.fits",
        twoDMapPath=str(mapPath),
    )

    assert result.orderNums == [10]
    assert result.mapDF is mapDataframe
    assert np.isnan(frame.data[0, 1])
    assert result.imageMap.shape == (4, 6)
    assert result.imageMap["order"].unique().tolist() == [10.0]
    result.twoDMap.close()
