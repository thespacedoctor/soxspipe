"""Toolkit integration contracts using real FITS and Astropy tables."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from astropy.table import Table
from numpy.testing import assert_allclose, assert_array_equal

from soxspipe.commonutils import toolkit

pytestmark = pytest.mark.integration


def _write_2d_map(path: Path) -> Path:
    wavelength = np.array([[500.0, 501.0], [0.0, 503.0]], dtype=np.float32)
    slit = np.array([[-1.0, 0.0], [0.0, 1.0]], dtype=np.float32)
    order = np.array([[10.0, 10.0], [np.nan, 11.0]], dtype=np.float32)
    fits.HDUList(
        [
            fits.PrimaryHDU(np.zeros((2, 2), dtype=np.float32)),
            fits.ImageHDU(wavelength, name="WAVELENGTH"),
            fits.ImageHDU(slit, name="SLIT"),
            fits.ImageHDU(order, name="ORDER"),
        ]
    ).writeto(path)
    return path


def test_two_d_map_dataframe_preserves_associated_frame_arrays(tmp_path: Path, log: object) -> None:
    mapPath = _write_2d_map(tmp_path / "map.fits")
    frame = CCDData(
        np.array([[1.0, 2.0], [3.0, 4.0]]),
        unit=u.electron,
        mask=np.array([[False, True], [False, False]]),
        uncertainty=StdDevUncertainty(np.full((2, 2), 0.5)),
        meta={"ARM": "NIR"},
    )

    result, interOrderMask = toolkit.twoD_disp_map_image_to_dataframe(
        log,
        slit_length=4,
        twoDMapPath=str(mapPath),
        kw=lambda key: {"SEQ_ARM": "ARM", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key],
        associatedFrame=frame,
        removeMaskedPixels=True,
        dispAxis="x",
    )

    assert result["wavelength"].tolist() == [500.0, 503.0]
    assert result["flux"].tolist() == [1.0, 4.0]
    assert result["error"].tolist() == [0.5, 0.5]
    assert result["mask"].tolist() == [False, False]
    assert_array_equal(interOrderMask, [[0.0, 0.0], [1.0, 0.0]])


def test_extinction_correction_reads_fits_table_and_uses_next_sample(
    tmp_path: Path,
) -> None:
    extinctionPath = tmp_path / "extinction.fits"
    Table({"WAVE": [4000.0, 5000.0, 6000.0], "MAG_AIRMASS": [0.1, 0.2, 0.3]}).write(extinctionPath)

    factors = toolkit.extinction_correction_factor(np.array([400.0, 450.0, 600.0]), str(extinctionPath), airmass=2.0)

    assert_allclose(factors, 10 ** (0.8 * np.array([0.1, 0.2, 0.3])))


def test_read_spectral_format_selects_full_or_reduced_wavelengths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
) -> None:
    formatPath = tmp_path / "format.fits"
    Table(
        {
            "ORDER": [10, 11],
            "WLMIN": [500.0, 600.0],
            "WLMAX": [550.0, 650.0],
            "WLMINFUL": [490.0, 590.0],
            "WLMAXFUL": [560.0, 660.0],
        }
    ).write(formatPath)
    detector = {
        "science-pixels": {"rows": {"end": 31}, "columns": {"end": 31}},
        "spectral format table": formatPath.name,
        "dispersion-axis": "x",
    }
    monkeypatch.setattr(
        toolkit,
        "detector_lookup",
        lambda **kwargs: SimpleNamespace(get=lambda arm: detector),
    )
    monkeypatch.setattr(toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **kwargs: str(tmp_path))

    full = toolkit.read_spectral_format(log, {}, "VIS")
    reduced = toolkit.read_spectral_format(log, {}, "VIS", extended=False)

    assert_array_equal(full[0], [10, 11])
    assert_array_equal(full[1], [490.0, 590.0])
    assert_array_equal(full[2], [560.0, 660.0])
    assert_array_equal(reduced[1], [500.0, 600.0])
    assert_array_equal(reduced[2], [550.0, 650.0])


@pytest.mark.parametrize(
    ("arm", "threshold", "expected"),
    [("VIS", 5.0, [501.0, 502.0]), ("NIR", 100.0, [501.0, 502.0])],
)
def test_get_skylines_dataframe_filters_real_fits_table(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    arm: str,
    threshold: float,
    expected: list[float],
) -> None:
    skyPath = tmp_path / "skylines.fits"
    Table(
        {
            "WAVELENGTH": [500.0, 501.0, 502.0],
            "FLUX": [threshold, threshold + 1.0, 150.0],
            "ISOLATED": [0, 1, 1],
        }
    ).write(skyPath)
    monkeypatch.setattr(
        toolkit,
        "detector_lookup",
        lambda **kwargs: SimpleNamespace(get=lambda value: {"skylines": skyPath.name}),
    )
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **kwargs: str(tmp_path))

    result = toolkit.get_skylines_dataframe(
        log,
        {},
        arm,
        minBrightnessVIS=5,
        minBrightnessNIR=100,
    )

    assert result["WAVELENGTH"].tolist() == expected
    assert result["WAVELENGTH"].dtype.kind == "f"
    assert result["FLUX"].dtype.kind == "f"
    assert result["ISOLATED"].dtype.kind == "b"


def _record_opened_hdulists(monkeypatch: pytest.MonkeyPatch) -> list:
    """Record every HDUList that ``fits.open`` returns, so a test can inspect its file handle afterwards."""
    opened: list = []
    realOpen = fits.open

    def recording_open(*args, **kwargs):
        hdul = realOpen(*args, **kwargs)
        opened.append(hdul)
        return hdul

    monkeypatch.setattr(fits, "open", recording_open)
    return opened


def _write_map(path: Path, *, hasOrder: bool = True, hasOrderData: bool = True) -> Path:
    """Write a small 2D map, optionally without its ORDER extension or without data in it."""
    image = np.ones((6, 4), dtype="float32")
    hdus = [fits.PrimaryHDU(image, header=fits.Header({"EXTNAME": "WAVELENGTH"})), fits.ImageHDU(image, name="SLIT")]
    if hasOrder:
        hdus.append(fits.ImageHDU(image if hasOrderData else None, name="ORDER"))
    fits.HDUList(hdus).writeto(path)
    return path


@pytest.mark.parametrize(
    ("mapOptions", "error"),
    [({"hasOrder": False}, KeyError), ({"hasOrderData": False}, AttributeError)],
)
def test_open_disp_map_hdus_closes_the_file_when_the_map_is_malformed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, mapOptions: dict, error: type
) -> None:
    path = _write_map(tmp_path / "map.fits", **mapOptions)
    opened = _record_opened_hdulists(monkeypatch)

    with pytest.raises(error):
        toolkit._open_disp_map_hdus(str(path), 1, 1)

    assert len(opened) == 1
    assert opened[0]._file.closed


def test_open_disp_map_hdus_leaves_the_file_open_for_the_caller_when_the_map_is_valid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_map(tmp_path / "map.fits")
    opened = _record_opened_hdulists(monkeypatch)

    hdul, minimumBinnedPixelValue, binned = toolkit._open_disp_map_hdus(str(path), 1, 1)

    assert hdul is opened[0]
    assert not hdul._file.closed
    assert minimumBinnedPixelValue is None
    assert not binned
    hdul.close()


def test_open_disp_map_hdus_closes_the_file_when_the_binning_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import astropy.nddata

    def refuse(*args, **kwargs):
        raise ValueError("cannot block-reduce")

    path = _write_map(tmp_path / "map.fits")
    opened = _record_opened_hdulists(monkeypatch)
    monkeypatch.setattr(astropy.nddata, "block_reduce", refuse)

    with pytest.raises(ValueError):
        toolkit._open_disp_map_hdus(str(path), 2, 2)

    assert opened[0]._file.closed
