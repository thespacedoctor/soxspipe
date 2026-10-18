"""Public response-function behavior with synthetic spectra."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.response_function import response_function
from tests.factories import product_table, qc_table

pytestmark = pytest.mark.unit


def test_response_constructor_reads_and_normalizes_synthetic_standard(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Initialize response state from temporary extracted and reference spectra."""
    extractedPath = tmp_path / "standard.fits"
    unflattenedPath = tmp_path / "standard-unflattened.fits"
    standardsPath = tmp_path / "standards.fits"
    header = fits.Header(
        {
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2024-01-02T03:04:05",
            "EXPTIME": 30.0,
            "OBJECT": "CD-30 17706",
            "HIERARCH ESO TEL AIRM START": 1.1,
            "HIERARCH ESO TEL AIRM END": 1.3,
        }
    )
    extractedTable = Table({"WAVE": [600.0, 500.0], "FLUX_COUNTS": [20.0, 10.0]})
    unflattenedTable = Table({"WAVE": [500.0, 600.0], "FLUX_DENSITY_COUNTS": [10.0, 20.0]})
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(extractedTable)]).writeto(extractedPath)
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(unflattenedTable)]).writeto(unflattenedPath)
    Table({"WAVE": [500.0, 600.0], "LTT7987": [1.0, 2.0]}).write(standardsPath)

    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit

    monkeypatch.setattr(
        commonutils,
        "keyword_lookup",
        lambda **_: SimpleNamespace(get=lambda key: key),
    )
    monkeypatch.setattr(
        commonutils,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"flux-standards": "standards.fits"}),
    )
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: str(tmp_path))
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")),
    )

    response = response_function(
        log=log,
        stdExtractionPath=str(extractedPath),
        recipeName="soxs-response",
        settings={"instrument": "soxs", "soxs-response": {"vis": {}}},
        qcTable=qc_table().iloc[0:0],
        productsTable=product_table().iloc[0:0],
        sofName="SYNTHETIC",
        startNightDate="2024-01-02",
        stdNotFlatExtractionPath=str(unflattenedPath),
        orderJoins=["12"],
    )

    assert response.arm == "VIS"
    assert response.std_objName == "LTT7987"
    assert response.stdExtractionDF["WAVE"].to_list() == [500.0, 600.0]
    assert response.stdExtractionNotFlatDF["WAVE"].to_list() == [500.0, 600.0]
    assert response.airmass == pytest.approx(1.2)
    assert response.qcDir == str(tmp_path / "qc")
    assert response.productDir == str(tmp_path / "products")


def _response_writer(tmp_path: Path, log: object) -> response_function:
    """Return a VIS response_function with only the state the RESP writer reads."""
    response = object.__new__(response_function)
    response.log = log
    response.arm = "VIS"
    response.kw = lambda key: key
    response.sofName = "SYNTHETIC"
    response.header = fits.Header({"OBJECT": "Synthetic standard"})
    response.productDir = str(tmp_path)
    response.qc = qc_table()
    response.products = product_table().iloc[0:0]
    response.recipeName = "soxs-response"
    response.dateObs = "2024-01-02T03:04:05"
    return response


def test_write_response_function_to_file_preserves_fits_and_product_contract(
    tmp_path: Path,
    log: object,
) -> None:
    """Write response coefficients, header QC values, and product metadata."""
    response = _response_writer(tmp_path, log)

    filePath = response.write_response_function_to_file(
        responseFuncCoeffs=np.array([0.5, 2.0]),
        polyOrder=1,
        wavelength=np.array([500.0, 600.0]),
        response=np.array([250.5, 302.0]),
    )

    assert filePath == str(tmp_path / "SYNTHETIC_RESP.fits")
    with fits.open(filePath, checksum=True) as hdul:
        assert hdul[0].header["OBJECT"] == "Synthetic standard"
        assert hdul[0].header["SEQ_ARM"] == "VIS"
        assert hdul[0].header["PRO_TYPE"] == "REDUCED"
        assert hdul[0].header["PRO_CATG"] == "RESP_TAB_VIS"
        assert hdul[0].header["ESO QC RON"] == 3.2
        assert hdul[0].header.comments["ESO QC RON"] == "Synthetic read noise"
        assert hdul[1].data["polyOrder"][0] == 1
        np.testing.assert_allclose(hdul[1].data["c0"], [0.5])
        np.testing.assert_allclose(hdul[1].data["c1"], [2.0])

    assert response.products.to_dict("records") == [
        {
            "soxspipe_recipe": "soxs-response",
            "product_label": "RESPONSE_FUNC",
            "file_name": "SYNTHETIC_RESP.fits",
            "file_type": "FITS",
            "obs_date_utc": "2024-01-02T03:04:05",
            "reduction_date_utc": response.products.loc[0, "reduction_date_utc"],
            "product_desc": "Response function coeffs.",
            "file_path": filePath,
            "label": "PROD",
        }
    ]


def _plot_fit_and_settings(wavelengths: np.ndarray):
    from soxspipe.commonutils.response_fit import fit_response_curve, resolve_response_settings

    fitSettings = resolve_response_settings({}, "VIS", "soxs")
    ratio = 2.0 * (1.0 + 0.001 * np.sin(wavelengths))
    ratio[100] = 5.0
    return fit_response_curve(wavelengths, ratio, fitSettings), fitSettings


@pytest.mark.parametrize(
    ("includeEfficiency", "qcWarnings"),
    [(False, []), (True, []), (False, ["response data gap of 99 nm exceeds 40 nm"])],
    ids=["plain", "efficiency", "qc-warning"],
)
def test_plot_response_curve_writes_qc_pdf_and_product_row(
    tmp_path: Path,
    log: object,
    includeEfficiency: bool,
    qcWarnings: list[str],
) -> None:
    """Render each response QC layout into the isolated QC directory."""
    wavelengths = np.linspace(400.0, 800.0, 1000)
    fit, fitSettings = _plot_fit_and_settings(wavelengths)
    response = object.__new__(response_function)
    response.log = log
    response.std_wavelength_to_abs_flux = lambda wave: np.full_like(wave, 2.0e17, dtype=float)
    response.std_objName = "SYNTHETIC_STAR"
    response.sofName = "SYNTHETIC"
    response.qcDir = str(tmp_path)
    response.products = product_table().iloc[0:0]
    response.recipeName = "soxs-response"
    response.dateObs = "2024-01-02T03:04:05"
    efficiency = np.linspace(0.2, 0.4, 1000) if includeEfficiency else None

    plotPath = response.plot_response_curve(
        stdExtWave=wavelengths,
        stdExtWaveNotFlat=wavelengths,
        stdExtFlux=np.full(1000, 2.0e17),
        fit=fit,
        fitSettings=fitSettings,
        stdEfficiencyEstimate=efficiency,
        qcWarnings=qcWarnings,
    )

    assert plotPath == str(tmp_path / "SYNTHETIC_RESPONSE.pdf")
    assert Path(plotPath).read_bytes().startswith(b"%PDF-")
    assert response.products.to_dict("records") == [
        {
            "soxspipe_recipe": "soxs-response",
            "product_label": "RESPONSE_QC_PLOT",
            "file_name": "SYNTHETIC_RESPONSE.pdf",
            "file_type": "PDF",
            "obs_date_utc": "2024-01-02T03:04:05",
            "reduction_date_utc": response.products.loc[0, "reduction_date_utc"],
            "product_desc": "Response curve QC plot.",
            "file_path": plotPath,
            "label": "QC",
        }
    ]


def test_write_response_function_to_file_adds_sampled_response_curve_extension(
    tmp_path: Path,
    log: object,
) -> None:
    """The RESP product carries the fitted curve as a WAVE/RESPONSE table after the coefficient table (DY-1378)."""
    # ARRANGE
    response = _response_writer(tmp_path, log)

    # ACT
    filePath = response.write_response_function_to_file(
        responseFuncCoeffs=np.array([0.5, 2.0]),
        polyOrder=1,
        wavelength=np.array([600.0, np.nan, 500.0, 550.0]),
        response=np.array([302.0, 999.0, 252.0, 277.0]),
    )

    # ASSERT
    with fits.open(filePath, checksum=True) as hdul:
        assert [hdu.name for hdu in hdul] == ["PRIMARY", "", "RESPONSE"]
        assert hdul[1].data["polyOrder"][0] == 1
        np.testing.assert_allclose(hdul[1].data["c0"], [0.5])
        np.testing.assert_allclose(hdul[1].data["c1"], [2.0])
    curve = Table.read(filePath, hdu="RESPONSE")
    assert curve.colnames == ["WAVE", "RESPONSE"]
    assert str(curve["WAVE"].unit) == "nm"
    assert curve["RESPONSE"].unit == u.Unit("1e-17 erg / (cm2 s Angstrom)") / (u.ct / u.s)
    np.testing.assert_array_equal(curve["WAVE"], [500.0, 550.0, 600.0])
    np.testing.assert_allclose(curve["RESPONSE"], [252.0, 277.0, 302.0], rtol=1e-14, atol=0)


def test_write_response_function_to_file_skips_nan_qc_values(
    tmp_path: Path,
    log: object,
) -> None:
    """A NaN QC value from an upstream recipe is left out of the header instead of crashing the write (DY-1259)."""
    # ARRANGE
    response = _response_writer(tmp_path, log)
    response.qc = pd.concat(
        [qc_table(), qc_table().assign(qc_name="INNER ORDER PIX MEAN", qc_value=float("nan"))],
        ignore_index=True,
    )

    # ACT
    filePath = response.write_response_function_to_file(
        responseFuncCoeffs=np.array([0.5, 2.0]),
        polyOrder=1,
        wavelength=np.array([500.0, 600.0]),
        response=np.array([250.5, 302.0]),
    )

    # ASSERT
    with fits.open(filePath, checksum=True) as hdul:
        assert hdul[0].header["ESO QC RON"] == 3.2
        assert "ESO QC INNER ORDER PIX MEAN" not in hdul[0].header


def test_get_writes_response_when_optional_efficiency_is_unavailable(
    log: object,
) -> None:
    """Produce a NIR response without an efficiency product for a flat synthetic standard."""
    wavelengths = np.linspace(1000.0, 2400.0, 3000)
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame({"WAVE": wavelengths, "FLUX_COUNTS": np.full(3000, 10.0)})
    # A THREE-SAMPLE UNFLATTENED SPECTRUM IS TOO SHORT FOR THE EFFICIENCY SMOOTHING, SO NO EFFICIENCY IS PRODUCED
    response.stdExtractionNotFlatDF = pd.DataFrame(
        {"WAVE": wavelengths[[0, 1500, 2999]], "FLUX_DENSITY_COUNTS": np.full(3, 10.0)}
    )
    response.stdAbsFluxDF = pd.DataFrame({"WAVE": wavelengths, "SYNTHETIC_STAR": np.full(3000, 2.0)})
    response.std_objName = "SYNTHETIC_STAR"
    response.arm = "NIR"
    response.instrument = "soxs"
    response.texp = 10.0
    response.recipeSettings = {"nir": {"max_iteration": 5, "poly_order": 2}}
    response.calibrationRootPath = "/calibration"
    response.detectorParams = {"extinction": "extinction.fits"}
    response.airmass = 1.0
    response.qc = pd.DataFrame()
    response.products = pd.DataFrame()
    response.recipeName = "soxs-response"
    response.dateObs = "2024-01-02T03:04:05"
    response.orderJoins = []
    captured: dict[str, object] = {}
    response.write_response_function_to_file = lambda **kwargs: captured.update({"response_write": kwargs})
    response.plot_response_curve = lambda **kwargs: captured.update({"plot": kwargs})

    qc, products, recipeError = response.get()

    assert recipeError is False
    assert "RESP BIN CLIP FRAC" in qc["qc_name"].tolist()
    assert products.empty
    # THE OLD-STYLE poly_order IN THE WORKSPACE BLOCK IS HONOURED AS THE LEGACY POLYNOMIAL ORDER
    assert captured["response_write"]["polyOrder"] == 2
    np.testing.assert_allclose(captured["response_write"]["response"], 2.0e17, rtol=1e-6)
    np.testing.assert_array_equal(captured["response_write"]["wavelength"], wavelengths)
    assert captured["plot"]["stdEfficiencyEstimate"] is None


def _build_real_nir_response(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, **tables: object):
    """Construct a NIR response_function through its real constructor from synthetic files."""
    import soxspipe.commonutils as commonutils
    import soxspipe.commonutils.toolkit as toolkit

    wavelengths = np.linspace(1000.0, 2400.0, 3000)
    header = fits.Header(
        {
            "SEQ_ARM": "NIR",
            "DATE_OBS": "2024-01-02T03:04:05",
            "EXPTIME": 10.0,
            "OBJECT": "LTT7987",
            "HIERARCH ESO TEL AIRM START": 1.0,
            "HIERARCH ESO TEL AIRM END": 1.0,
        }
    )
    extractedPath = tmp_path / "standard.fits"
    unflattenedPath = tmp_path / "standard-unflattened.fits"
    extracted = Table({"WAVE": wavelengths, "FLUX_COUNTS": np.full(3000, 10.0)})
    unflattened = Table({"WAVE": wavelengths, "FLUX_DENSITY_COUNTS": np.full(3000, 10.0)})
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(extracted)]).writeto(extractedPath)
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(unflattened)]).writeto(unflattenedPath)
    Table({"WAVE": wavelengths, "LTT7987": np.full(3000, 2.0)}).write(tmp_path / "standards.fits")
    (tmp_path / "qc").mkdir()
    (tmp_path / "products").mkdir()

    monkeypatch.setattr(commonutils, "keyword_lookup", lambda **_: SimpleNamespace(get=lambda key: key))
    monkeypatch.setattr(
        commonutils,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"flux-standards": "standards.fits"}),
    )
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: str(tmp_path))
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))
    return response_function(
        log=log,
        stdExtractionPath=str(extractedPath),
        recipeName="soxs-response",
        settings={"instrument": "soxs", "soxs-response": {"nir": {}}},
        sofName="SYNTHETIC",
        startNightDate="2024-01-02",
        stdNotFlatExtractionPath=str(unflattenedPath),
        orderJoins={"1-2": 1500.0},
        **tables,
    )


def test_get_builds_empty_qc_and_product_tables_when_none_are_passed(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Constructing without qcTable and productsTable still yields QC rows and product rows from get()."""
    # ARRANGE
    import soxspipe.commonutils.phase3 as phase3

    monkeypatch.setattr(phase3, "write_fits_table_to_disk", lambda **_: None)
    response = _build_real_nir_response(tmp_path, log, monkeypatch)

    # ACT
    qc, products, recipeError = response.get()

    # ASSERT
    assert recipeError is False
    assert "RESP BIN CLIP FRAC" in qc["qc_name"].tolist()
    assert {"RESPONSE_FUNC", "RESPONSE_QC_PLOT", "EFFICIENCY"} <= set(products["product_label"])
    assert list(qc.columns) == list(qc_table().columns)
    assert list(products.columns) == list(product_table().columns)


def test_plot_response_curve_closes_only_its_own_figure(tmp_path: Path, log: object) -> None:
    """A figure the caller has open survives the QC plot, and the QC figure itself is closed."""
    import matplotlib.pyplot as plt

    # ARRANGE
    wavelengths = np.linspace(400.0, 800.0, 1000)
    fit, fitSettings = _plot_fit_and_settings(wavelengths)
    response = object.__new__(response_function)
    response.log = log
    response.std_wavelength_to_abs_flux = lambda wave: np.full_like(wave, 2.0e17, dtype=float)
    response.std_objName = "SYNTHETIC_STAR"
    response.sofName = "SYNTHETIC"
    response.qcDir = str(tmp_path)
    response.products = product_table().iloc[0:0]
    response.recipeName = "soxs-response"
    response.dateObs = "2024-01-02T03:04:05"
    plt.close("all")
    callerFigure = plt.figure()

    # ACT
    response.plot_response_curve(
        stdExtWave=wavelengths,
        stdExtWaveNotFlat=wavelengths,
        stdExtFlux=np.full(1000, 2.0e17),
        fit=fit,
        fitSettings=fitSettings,
        stdEfficiencyEstimate=None,
        qcWarnings=[],
    )

    # ASSERT
    assert plt.get_fignums() == [callerFigure.number]
    plt.close(callerFigure)
