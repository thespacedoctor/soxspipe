"""Characterization contracts pinning the toolkit branches DY-82 splits."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import Table
from numpy.testing import assert_allclose, assert_array_equal

from soxspipe.commonutils import toolkit
from tests.factories import order_table_fits

pytestmark = pytest.mark.integration

KEYWORD_MAP = {"SEQ_ARM": "ARM", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}


def _x_dispersion_order_table(tmp_path: Path) -> Path:
    """Write an order table whose polynomials are keyed on the x-axis."""
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 1,
                "degx_cent": 1,
                "cent_00": 1.0,
                "cent_01": 2.0,
                "cent_10": 3.0,
                "cent_11": 4.0,
                "degorder_edgeup": 1,
                "degx_edgeup": 1,
                "edgeup_00": 5.0,
                "edgeup_01": 0.5,
                "edgeup_10": 0.25,
                "edgeup_11": 0.125,
                "degorder_edgelow": 1,
                "degx_edgelow": 1,
                "edgelow_00": -5.0,
                "edgelow_01": -0.5,
                "edgelow_10": -0.25,
                "edgelow_11": -0.125,
            }
        ]
    )
    metadata = pd.DataFrame(
        {
            "order": [10, 11],
            "xmin": [0.0, 2.0],
            "xmax": [4.0, 6.0],
        }
    )
    return order_table_fits(
        tmp_path / "x_orders.fits",
        polynomials=polynomials,
        metadata=metadata,
    )


def test_unpack_order_table_evaluates_edge_polynomials_on_an_x_axis_table(
    tmp_path: Path,
    log: object,
) -> None:
    orderPath = _x_dispersion_order_table(tmp_path)

    unpacked = toolkit.unpack_order_table(
        log=log,
        orderTablePath=str(orderPath),
        pixelDelta=1,
    )

    assert len(unpacked) == 3
    _, pixelTable, metadataTable = unpacked
    assert list(pixelTable["xcoord"]) == [0, 1, 2, 3, 2, 3, 4, 5]
    assert "std" not in pixelTable.columns

    xcoord = pixelTable["xcoord"]
    order = pixelTable["order"]
    assert_allclose(
        pixelTable["ycoord_centre"],
        1.0 + 2.0 * xcoord + 3.0 * order + 4.0 * order * xcoord,
        rtol=1e-12,
        atol=0,
    )
    assert_allclose(
        pixelTable["ycoord_edgeup"],
        5.0 + 0.5 * xcoord + 0.25 * order + 0.125 * order * xcoord,
        rtol=1e-12,
        atol=0,
    )
    assert_allclose(
        pixelTable["ycoord_edgelow"],
        -5.0 - 0.5 * xcoord - 0.25 * order - 0.125 * order * xcoord,
        rtol=1e-12,
        atol=0,
    )
    assert list(metadataTable["xmin"]) == [0.0, 2.0]


def test_unpack_order_table_scales_edge_coordinates_by_the_axis_a_binning(
    tmp_path: Path,
    log: object,
) -> None:
    orderPath = _x_dispersion_order_table(tmp_path)

    _, unbinned, _ = toolkit.unpack_order_table(log=log, orderTablePath=str(orderPath))
    _, binned, _ = toolkit.unpack_order_table(
        log=log,
        orderTablePath=str(orderPath),
        biny=2,
    )

    for column in ["ycoord_centre", "ycoord_edgeup", "ycoord_edgelow"]:
        assert_allclose(binned[column], unbinned[column] / 2, rtol=1e-12, atol=0)
    assert list(binned["xcoord"]) == list(unbinned["xcoord"])


def _binned_2d_map(path: Path) -> Path:
    """Write a 2D dispersion map whose extensions are twice the frame size."""
    wavelength = np.arange(16, dtype=np.float32).reshape(4, 4) + 500.0
    slit = np.zeros((4, 4), dtype=np.float32)
    order = np.full((4, 4), 10.0, dtype=np.float32)
    fits.HDUList(
        [
            fits.PrimaryHDU(np.zeros((2, 2), dtype=np.float32)),
            fits.ImageHDU(wavelength, name="WAVELENGTH"),
            fits.ImageHDU(slit, name="SLIT"),
            fits.ImageHDU(order, name="ORDER"),
        ]
    ).writeto(path)
    return path


def test_two_d_disp_map_dataframe_block_reduces_a_binned_frame(
    tmp_path: Path,
    log: object,
) -> None:
    mapPath = _binned_2d_map(tmp_path / "binned_map.fits")
    frame = CCDData(
        np.arange(4, dtype=float).reshape(2, 2),
        unit=u.electron,
        mask=np.zeros((2, 2), dtype=bool),
        uncertainty=StdDevUncertainty(np.full((2, 2), 0.5)),
        meta={"ARM": "VIS", "BINX": 2, "BINY": 2},
    )

    result, interOrderMask = toolkit.twoD_disp_map_image_to_dataframe(
        log,
        slit_length=4,
        twoDMapPath=str(mapPath),
        kw=lambda key: KEYWORD_MAP[key],
        associatedFrame=frame,
        dispAxis="y",
    )

    assert list(result["x"]) == [0, 1, 0, 1]
    assert list(result["y"]) == [0, 0, 1, 1]
    assert_allclose(
        result.sort_values(["y", "x"], kind="stable")["wavelength"],
        [502.5, 504.5, 510.5, 512.5],
        rtol=1e-12,
        atol=0,
    )
    assert_allclose(
        result.sort_values(["y", "x"], kind="stable")["min"],
        [500.0, 502.0, 508.0, 510.0],
        rtol=1e-12,
        atol=0,
    )
    assert list(result["flux"]) == [0.0, 1.0, 2.0, 3.0]
    assert_array_equal(interOrderMask, np.zeros((2, 2)))


def _unbinned_2d_map(path: Path) -> Path:
    """Write a 2D dispersion map whose extensions match the primary HDU size."""
    wavelength = np.arange(4, dtype=np.float32).reshape(2, 2) + 500.0
    fits.HDUList(
        [
            fits.PrimaryHDU(np.zeros((2, 2), dtype=np.float32)),
            fits.ImageHDU(wavelength, name="WAVELENGTH"),
            fits.ImageHDU(np.zeros((2, 2), dtype=np.float32), name="SLIT"),
            fits.ImageHDU(np.full((2, 2), 10.0, dtype=np.float32), name="ORDER"),
        ]
    ).writeto(path)
    return path


def test_two_d_disp_map_dataframe_expands_a_home_relative_path(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
) -> None:
    _unbinned_2d_map(tmp_path / "home_map.fits")
    monkeypatch.setenv("HOME", str(tmp_path))

    result, _ = toolkit.twoD_disp_map_image_to_dataframe(
        log,
        slit_length=4,
        twoDMapPath="~/home_map.fits",
        dispAxis="y",
    )

    assert_allclose(result["wavelength"], [500.0, 501.0, 502.0, 503.0], rtol=1e-12, atol=0)
    assert list(result["min"]) == [1.0, 1.0, 1.0, 1.0]
    assert "flux" not in result.columns


def _spectral_format_detector(tmp_path: Path, dispersionAxis: str) -> dict:
    formatPath = tmp_path / "format.fits"
    if not formatPath.exists():
        Table(
            {
                "ORDER": [10, 11],
                "WLMINFUL": [490.0, 590.0],
                "WLMAXFUL": [560.0, 660.0],
            }
        ).write(formatPath)
    return {
        "science-pixels": {"rows": {"end": 20}, "columns": {"end": 20}},
        "spectral format table": formatPath.name,
        "dispersion-axis": dispersionAxis,
    }


def _patch_spectral_format_lookups(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    detector: dict,
) -> None:
    monkeypatch.setattr(
        toolkit,
        "detector_lookup",
        lambda **kwargs: SimpleNamespace(get=lambda arm: detector),
    )
    monkeypatch.setattr(toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **kwargs: str(tmp_path))
    monkeypatch.setattr(
        toolkit,
        "dispersion_map_to_pixel_arrays",
        lambda log, dispersionMapPath, orderPixelTable, **kwargs: orderPixelTable.assign(
            fit_x=[-4.0, 8.0, 12.0, 40.0],
            fit_y=[-6.0, 9.0, 14.0, 44.0],
        ),
    )


@pytest.mark.parametrize(
    ("dispersionAxis", "binx", "biny", "expectedMins", "expectedMaxs"),
    [
        ("y", 2, 1, [0.0, 6.0], [4.0, 10.0]),
        ("x", 1, 4, [0.0, 3.5], [2.25, 5.0]),
    ],
)
def test_read_spectral_format_clamps_and_bins_the_dispersion_pixel_limits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    dispersionAxis: str,
    binx: int,
    biny: int,
    expectedMins: list[float],
    expectedMaxs: list[float],
) -> None:
    detector = _spectral_format_detector(tmp_path, dispersionAxis)
    _patch_spectral_format_lookups(monkeypatch, tmp_path, detector)

    result = toolkit.read_spectral_format(
        log,
        {},
        "VIS",
        dispersionMap=str(tmp_path / "map.fits"),
        binx=binx,
        biny=biny,
    )

    assert len(result) == 5
    orderNums, waveLengthMin, waveLengthMax, amins, amaxs = result
    assert_array_equal(orderNums, [10, 11])
    assert_array_equal(waveLengthMin, [490.0, 590.0])
    assert_array_equal(waveLengthMax, [560.0, 660.0])
    assert_allclose(amins, expectedMins, rtol=1e-12, atol=0)
    assert_allclose(amaxs, expectedMaxs, rtol=1e-12, atol=0)
