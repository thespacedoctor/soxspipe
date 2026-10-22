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

LASILLA_EXTINCTION_PATH = (
    Path(toolkit.__file__).resolve().parent.parent / "resources/static_calibrations/soxs/extinction_lasilla.fits"
)
# A STEP OF THIS FRACTION BETWEEN NEIGHBOURING DENSE-GRID SAMPLES IS A JUMP, NOT THE CURVE'S SLOPE
MAX_DENSE_GRID_STEP = 1e-3


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


def test_extinction_correction_reads_fits_table_and_interpolates_linearly(
    tmp_path: Path,
) -> None:
    # ARRANGE
    extinctionPath = tmp_path / "extinction.fits"
    Table({"WAVE": [4000.0, 5000.0, 6000.0], "MAG_AIRMASS": [0.1, 0.2, 0.3]}).write(extinctionPath)

    # ACT
    factors = toolkit.extinction_correction_factor(np.array([400.0, 450.0, 600.0]), str(extinctionPath), airmass=2.0)

    # ASSERT: 450 NM IS HALFWAY BETWEEN THE 0.1 AND 0.2 MAG/AIRMASS NODES
    assert_allclose(factors, 10 ** (0.8 * np.array([0.1, 0.15, 0.3])))


def test_extinction_correction_has_no_jump_at_table_wavelengths_and_hits_table_values() -> None:
    # ARRANGE
    table = Table.read(LASILLA_EXTINCTION_PATH)
    nodesNm = np.asarray(table["WAVE"], float) / 10
    magPerAirmass = np.asarray(table["MAG_AIRMASS"], float)
    airmass = 1.37
    epsilonNm = 1e-6
    interiorNodes = nodesNm[1:-1]

    # ACT
    atNodes = toolkit.extinction_correction_factor(nodesNm, str(LASILLA_EXTINCTION_PATH), airmass)
    justBelow = toolkit.extinction_correction_factor(interiorNodes - epsilonNm, str(LASILLA_EXTINCTION_PATH), airmass)
    justAbove = toolkit.extinction_correction_factor(interiorNodes + epsilonNm, str(LASILLA_EXTINCTION_PATH), airmass)
    denseWave = np.linspace(nodesNm[0], nodesNm[-1], 200001)
    denseFactors = toolkit.extinction_correction_factor(denseWave, str(LASILLA_EXTINCTION_PATH), airmass)

    # ASSERT: EQUAL TO THE TABLE AT EVERY NODE
    assert_allclose(atNodes, 10 ** (0.4 * magPerAirmass * airmass), rtol=1e-6)
    # ASSERT: CONTINUOUS ACROSS EVERY INTERIOR NODE (THE OLD STAIRCASE JUMPED BY OVER 50% AT 320 NM)
    assert_allclose(justBelow, justAbove, rtol=1e-5)
    # ASSERT: NO STEP ANYWHERE ON A DENSE GRID (0.003 NM SPACING; STEEPEST SLOPE IS ~0.07 %/0.003 NM AT 310-320 NM)
    assert np.max(np.abs(np.diff(denseFactors) / denseFactors[:-1])) < MAX_DENSE_GRID_STEP


def test_extinction_correction_holds_end_values_outside_the_table(tmp_path: Path) -> None:
    # ARRANGE
    extinctionPath = tmp_path / "extinction.fits"
    Table({"WAVE": [4000.0, 5000.0, 6000.0], "MAG_AIRMASS": [0.1, 0.2, 0.3]}).write(extinctionPath)
    # A NAN WAVELENGTH GIVES NAN, SO A MASKED PIXEL STAYS MASKED
    wave = np.array([300.0, 399.0, 601.0, 900.0, 1020.0, np.nan])

    # ACT
    factors = toolkit.extinction_correction_factor(wave, str(extinctionPath), airmass=2.0)

    # ASSERT: BELOW THE TABLE HOLDS THE FIRST VALUE; ABOVE THE TABLE HOLDS THE LAST VALUE (DY-1688)
    assert_allclose(factors[:2], 10 ** (0.8 * 0.1))
    assert_allclose(factors[2:5], 10 ** (0.8 * 0.3))
    assert np.isnan(factors[5])


def test_extinction_correction_is_finite_above_the_lasilla_table_end() -> None:
    # ARRANGE
    table = Table.read(LASILLA_EXTINCTION_PATH)
    lastMagPerAirmass = float(np.asarray(table["MAG_AIRMASS"], float)[-1])
    airmass = 1.3
    xshooterVisWaveNm = np.linspace(900.5, 1020.0, 50)

    # ACT
    factors = toolkit.extinction_correction_factor(xshooterVisWaveNm, str(LASILLA_EXTINCTION_PATH), airmass)

    # ASSERT: THE X-SHOOTER VIS RANGE ABOVE 900 NM IS NOT BLANKED (DY-1688)
    assert_allclose(factors, 10 ** (0.4 * lastMagPerAirmass * airmass))


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
        toolkit._open_disp_map_hdus(str(path), 1, 1, kw=None)

    assert len(opened) == 1
    assert opened[0]._file.closed


def test_open_disp_map_hdus_leaves_the_file_open_for_the_caller_when_the_map_is_valid(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_map(tmp_path / "map.fits")
    opened = _record_opened_hdulists(monkeypatch)

    hdul, minimumBinnedPixelValue, binned = toolkit._open_disp_map_hdus(str(path), 1, 1, kw=None)

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
        toolkit._open_disp_map_hdus(str(path), 2, 2, kw=None)

    assert opened[0]._file.closed


def test_open_disp_map_hdus_closes_the_file_when_the_map_binning_does_not_divide_the_frame_binning(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = _write_map(tmp_path / "map.fits")
    fits.setval(path, "BINX", value=2, ext=0)
    fits.setval(path, "BINY", value=2, ext=0)
    opened = _record_opened_hdulists(monkeypatch)

    with pytest.raises(ValueError, match="cannot be block-reduced"):
        toolkit._open_disp_map_hdus(str(path), 3, 2, kw={"SEQ_ARM": "ARM", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}.get)

    assert opened[0]._file.closed
