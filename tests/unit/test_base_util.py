"""Initialization contracts for the shared utility base class."""

from __future__ import annotations

import gc
import warnings
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


@pytest.mark.parametrize(("dispersionAxis", "expectedAxes"), [("x", ("x", "y")), ("y", ("y", "x"))])
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
    assert result.twoDMap["WAVELENGTH"].data.shape == (2, 2)
    assert result.twoDMap["WAVELENGTH"].data.tolist() == [[500.0, 500.0], [500.0, 500.0]]
    assert result.imageMap["wavelength"].tolist() == [500.0] * 4
    assert result.imageMap["slit_position"].tolist() == [0.5] * 4


def test_base_util_imports_the_toolkit_helpers_before_it_branches(monkeypatch: pytest.MonkeyPatch, log: object) -> None:
    _patch_lookups(monkeypatch, "x")
    monkeypatch.delattr(toolkit, "read_spectral_format")
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02"},
    )

    with pytest.raises(ImportError):
        base_util(log, {}, associatedFrame=frame)


def _write_two_d_map(
    mapPath: Path,
    shape: tuple[int, int],
    binning: tuple[int, int] | None,
    planes: tuple[tuple[str, float], ...] = (("WAVELENGTH", 500.0), ("SLIT", 0.5), ("ORDER", 10.0)),
) -> None:
    """Write a 2D dispersion map file, optionally without its binning keywords or some of its planes."""
    primary = fits.PrimaryHDU(np.ones(shape))
    if binning is not None:
        primary.header["WIN_BINX"] = binning[0]
        primary.header["WIN_BINY"] = binning[1]
    hdus = [primary]
    for name, value in planes:
        hdus.append(fits.ImageHDU(np.full(shape, value), name=name))
    fits.HDUList(hdus).writeto(mapPath)


def _patch_two_d_map_helpers(monkeypatch: pytest.MonkeyPatch, interOrderMask: np.ndarray) -> None:
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


def test_base_util_reads_binning_from_the_associated_frame_header(monkeypatch: pytest.MonkeyPatch, log: object) -> None:
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
    assert any("dpBinx = self.twoDMap[0].header" in message for level, message in log.messages if level == "debug")


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


def test_base_util_raises_type_error_when_associated_frame_is_omitted(
    monkeypatch: pytest.MonkeyPatch, log: object
) -> None:
    _patch_lookups(monkeypatch, "x")

    with pytest.raises(TypeError):
        base_util(log, {})


def test_base_util_emits_only_its_instantiation_debug_message(monkeypatch: pytest.MonkeyPatch, log: object) -> None:
    _patch_lookups(monkeypatch, "x")
    frame = CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02", "WIN_BINX": 1, "WIN_BINY": 1},
    )

    base_util(log, {}, associatedFrame=frame)

    assert log.messages == [("debug", "instansiating a new 'base_util' object")]


def _resource_warnings_naming(path: Path, caught: list[warnings.WarningMessage]) -> list[str]:
    """Return the messages of the caught ResourceWarnings that name ``path``."""
    return [str(item.message) for item in caught if item.category is ResourceWarning and str(path) in str(item.message)]


def _binned_frame(binning: int) -> CCDData:
    """Build a 2x2 VIS frame with the same binning on both axes."""
    return CCDData(
        np.ones((2, 2)),
        unit="adu",
        meta={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02", "WIN_BINX": binning, "WIN_BINY": binning},
    )


def test_base_util_closes_the_two_d_map_file_before_the_constructor_returns(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    # ARRANGE
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_closed.fits"
    _write_two_d_map(mapPath, (2, 2), binning=(1, 1))

    # ACT
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        result = base_util(log, {}, associatedFrame=_binned_frame(1), twoDMapPath=str(mapPath))
        wavelength = result.twoDMap["WAVELENGTH"].data.tolist()
        wavelengthHeader = result.twoDMap["WAVELENGTH"].header
        del result
        gc.collect()

    # ASSERT: THE MAP PLANES AND HEADERS STAY READABLE AFTER CONSTRUCTION
    assert wavelength == [[500.0, 500.0], [500.0, 500.0]]
    assert wavelengthHeader["EXTNAME"] == "WAVELENGTH"

    # ASSERT: NO FILE HANDLE WAS LEFT FOR THE GARBAGE COLLECTOR TO CLOSE
    assert _resource_warnings_naming(mapPath, caught) == []


def test_base_util_closes_the_two_d_map_file_when_it_block_reduces_the_map(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    # ARRANGE
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_block_reduced.fits"
    _write_two_d_map(mapPath, (4, 4), binning=(1, 1))

    # ACT
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        result = base_util(log, {}, associatedFrame=_binned_frame(2), twoDMapPath=str(mapPath))
        planeShapes = [result.twoDMap[name].data.shape for name in ("WAVELENGTH", "SLIT", "ORDER")]
        slitValues = result.twoDMap["SLIT"].data.tolist()
        del result
        gc.collect()

    # ASSERT: THE PLANES KEEP THEIR BLOCK-REDUCED SHAPE AND VALUES
    assert planeShapes == [(2, 2), (2, 2), (2, 2)]
    assert slitValues == [[0.5, 0.5], [0.5, 0.5]]

    # ASSERT: NO FILE HANDLE WAS LEFT FOR THE GARBAGE COLLECTOR TO CLOSE
    assert _resource_warnings_naming(mapPath, caught) == []


def test_base_util_closes_the_two_d_map_file_when_construction_raises(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: object
) -> None:
    # ARRANGE: A MAP WITH NO ORDER EXTENSION, SO THE CONSTRUCTOR RAISES AFTER IT OPENS THE FILE
    _patch_lookups(monkeypatch, "x")
    _patch_two_d_map_helpers(monkeypatch, np.zeros((2, 2), dtype=bool))
    mapPath = tmp_path / "map_without_order.fits"
    _write_two_d_map(mapPath, (2, 2), binning=(1, 1), planes=(("WAVELENGTH", 500.0), ("SLIT", 0.5)))

    # ACT
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ResourceWarning)
        with pytest.raises(KeyError, match="ORDER"):
            base_util(log, {}, associatedFrame=_binned_frame(1), twoDMapPath=str(mapPath))
        gc.collect()

    # ASSERT
    assert _resource_warnings_naming(mapPath, caught) == []
