"""Public response-function behavior with synthetic spectra."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

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
    extractedTable = Table(
        {"WAVE": [600.0, 500.0], "FLUX_COUNTS": [20.0, 10.0]}
    )
    unflattenedTable = Table(
        {"WAVE": [500.0, 600.0], "FLUX_DENSITY_COUNTS": [10.0, 20.0]}
    )
    fits.HDUList(
        [fits.PrimaryHDU(header=header), fits.BinTableHDU(extractedTable)]
    ).writeto(extractedPath)
    fits.HDUList(
        [fits.PrimaryHDU(header=header), fits.BinTableHDU(unflattenedTable)]
    ).writeto(unflattenedPath)
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


def test_write_response_function_to_file_preserves_fits_and_product_contract(
    tmp_path: Path,
    log: object,
) -> None:
    """Write response coefficients, header QC values, and product metadata."""
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

    filePath = response.write_response_function_to_file(
        responseFuncCoeffs=np.array([0.5, 2.0]),
        polyOrder=1,
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


@pytest.mark.parametrize("includeEfficiency", [False, True])
def test_plot_response_curve_writes_qc_pdf_and_product_row(
    tmp_path: Path,
    log: object,
    includeEfficiency: bool,
) -> None:
    """Render each response QC layout into the isolated QC directory."""
    wavelengths = np.array([500.0, 550.0, 600.0])
    response = object.__new__(response_function)
    response.log = log
    response.std_wavelength_to_abs_flux = lambda wave: np.full_like(
        wave, 2.0e17, dtype=float
    )
    response.std_objName = "SYNTHETIC_STAR"
    response.sofName = "SYNTHETIC"
    response.qcDir = str(tmp_path)
    response.products = product_table().iloc[0:0]
    response.recipeName = "soxs-response"
    response.dateObs = "2024-01-02T03:04:05"
    efficiency = np.array([0.2, 0.3, 0.4]) if includeEfficiency else None

    plotPath = response.plot_response_curve(
        stdExtWave=wavelengths,
        stdExtWaveNotFlat=wavelengths,
        stdExtFlux=np.full(3, 2.0e17),
        binCentreWave=wavelengths,
        binCentreWaveOriginal=wavelengths,
        binIntegratedFlux=np.array([1.0, 1.1, 1.2]),
        absToExtFluxRatio=None,
        responseFuncCoeffs=np.array([0.0, 1.0]),
        stdEfficiencyEstimate=efficiency,
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


def test_get_writes_response_when_optional_efficiency_is_unavailable(
    log: object,
) -> None:
    """Produce a NIR response without an efficiency product for a short spectrum."""
    wavelengths = np.array([900.0, 1000.0, 1100.0])
    response = object.__new__(response_function)
    response.log = log
    response.stdExtractionDF = pd.DataFrame(
        {"WAVE": wavelengths, "FLUX_COUNTS": np.full(3, 10.0)}
    )
    response.stdExtractionNotFlatDF = pd.DataFrame(
        {"WAVE": wavelengths, "FLUX_DENSITY_COUNTS": np.full(3, 10.0)}
    )
    response.stdAbsFluxDF = pd.DataFrame(
        {"WAVE": wavelengths, "SYNTHETIC_STAR": np.full(3, 2.0)}
    )
    response.std_objName = "SYNTHETIC_STAR"
    response.arm = "NIR"
    response.instrument = "soxs"
    response.texp = 10.0
    response.recipeSettings = {"nir": {"poly_order": 0, "max_iteration": 2}}
    response.calibrationRootPath = "/calibration"
    response.detectorParams = {"extinction": "extinction.fits"}
    response.airmass = 1.0
    response.qc = pd.DataFrame()
    response.products = pd.DataFrame()
    response.orderJoins = []
    captured: dict[str, object] = {}
    response.write_response_function_to_file = lambda **kwargs: captured.update(
        {"response_write": kwargs}
    )
    response.plot_response_curve = lambda **kwargs: captured.update({"plot": kwargs})

    qc, products, recipeError = response.get()

    assert recipeError is False
    assert qc.empty
    assert products.empty
    np.testing.assert_allclose(
        captured["response_write"]["responseFuncCoeffs"], [2.0e17]
    )
    assert captured["response_write"]["polyOrder"] == 0
    assert captured["plot"]["stdEfficiencyEstimate"] is None
