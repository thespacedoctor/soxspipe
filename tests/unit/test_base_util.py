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


def _write_two_d_map(
    mapPath: Path, shape: tuple[int, int], binning: tuple[int, int] | None
) -> None:
    """Write a 2D dispersion map file, optionally without its binning keywords."""
    primary = fits.PrimaryHDU(np.ones(shape))
    if binning is not None:
        primary.header["WIN_BINX"] = binning[0]
        primary.header["WIN_BINY"] = binning[1]
    hdus = [primary]
    for name, value in (("WAVELENGTH", 500.0), ("SLIT", 0.5), ("ORDER", 10.0)):
        hdus.append(fits.ImageHDU(np.full(shape, value), name=name))
    fits.HDUList(hdus).writeto(mapPath)


def _patch_two_d_map_helpers(
    monkeypatch: pytest.MonkeyPatch, interOrderMask: np.ndarray
) -> None:
    """Stub the toolkit helpers the 2D map branch of the constructor calls."""
    monkeypatch.setattr(
        toolkit,
        "read_spectral_format",
        lambda **kwargs: ([10], [400.0], [500.0], [1], [2]),
    )
    monkeypatch.setattr(
        toolkit,
        "twoD_disp_map_image_to_dataframe",
        lambda **kwargs: (pd.DataFrame({"order": [10]}), interOrderMask),
    )


def test_base_util_reads_binning_from_the_associated_frame_header(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2024-01-02",
            "WIN_BINX": "2",
            "WIN_BINY": "3",
        },
    )

    result = base_util(log, {}, associatedFrame=frame)

    assert (result.binx, result.biny) == (2, 3)
    assert result.arm == "VIS"
    assert result.dateObs == "2024-01-02"


def test_base_util_defaults_binning_when_the_map_has_no_binning_keywords(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_without_binning.fits"
    _write_two_d_map(mapPath, (2, 2), binning=None)
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02"},
    )

    result = base_util(log, {}, associatedFrame=frame, twoDMapPath=str(mapPath))

    assert result.imageMap.shape == (4, 6)
    assert result.twoDMap["WAVELENGTH"].data.shape == (2, 2)
    assert any(
        "dpBinx = self.twoDMap[0].header" in message
        for level, message in log.messages
        if level == "debug"
    )
    result.twoDMap.close()


def test_base_util_builds_the_image_map_without_rebinning_when_the_ratios_are_one(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_unbinned.fits"
    _write_two_d_map(mapPath, (2, 2), binning=(1, 1))
    frame = CCDData(
        np.arange(4, dtype=float).reshape((2, 2)),
        unit="adu",
        meta={
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2024-01-02",
            "WIN_BINX": 1,
            "WIN_BINY": 1,
        },
    )

    result = base_util(log, {}, associatedFrame=frame, twoDMapPath=str(mapPath))

    assert list(result.imageMap.columns) == [
        "x",
        "y",
        "wavelength",
        "slit_position",
        "order",
        "flux",
    ]
    assert result.imageMap["x"].tolist() == [0, 1, 0, 1]
    assert result.imageMap["y"].tolist() == [0, 0, 1, 1]
    assert result.imageMap["wavelength"].tolist() == [500.0] * 4
    assert result.imageMap["slit_position"].tolist() == [0.5] * 4
    assert result.imageMap["flux"].tolist() == [0.0, 1.0, 2.0, 3.0]
    assert [str(dtype) for dtype in result.imageMap.dtypes.tolist()] == [
        "int64",
        "int64",
        "float32",
        "float32",
        "float32",
        "float32",
    ]
    assert result.twoDMap["WAVELENGTH"].data.shape == (2, 2)
    result.twoDMap.close()


def test_base_util_drops_image_map_rows_masked_out_of_every_map_plane(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_with_gap.fits"
    primary = fits.PrimaryHDU(np.ones((2, 2)))
    primary.header["WIN_BINX"] = 1
    primary.header["WIN_BINY"] = 1
    hdus = [primary]
    for name, value in (("WAVELENGTH", 500.0), ("SLIT", 0.5), ("ORDER", 10.0)):
        plane = np.full((2, 2), value)
        plane[1, 1] = np.nan
        hdus.append(fits.ImageHDU(plane, name=name))
    fits.HDUList(hdus).writeto(mapPath)
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02"},
    )

    result = base_util(log, {}, associatedFrame=frame, twoDMapPath=str(mapPath))

    assert result.imageMap.shape == (3, 6)
    assert result.imageMap.index.tolist() == [0, 1, 2]
    result.twoDMap.close()


def test_base_util_raises_when_no_associated_frame_supplies_the_arm(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")

    with pytest.raises(AttributeError):
        base_util(log, {})


def test_base_util_emits_only_its_instantiation_debug_message(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02", "WIN_BINX": 1, "WIN_BINY": 1},
    )

    base_util(log, {}, associatedFrame=frame)

    assert log.messages == [("debug", "instansiating a new 'base_util' object")]
