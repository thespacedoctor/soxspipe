"""Integration contracts for calibrated extracted-spectrum products."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.flux_calibration import flux_calibration
from tests.factories import pipeline_settings

pytestmark = pytest.mark.integration


def test_constructor_resolves_keywords_and_isolated_output_directories(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Initialize calibration state without creating a real workspace."""
    import soxspipe.commonutils.toolkit as toolkit

    settings = pipeline_settings(tmp_path)
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")),
    )

    calibration = flux_calibration(
        log=log,
        responseFunction="response.fits",
        extractedSpectrum=pd.DataFrame({"WAVE": [500.0], "FLUX_COUNTS": [20.0]}),
        settings=settings,
        airmass=1.2,
        exptime=30.0,
        extinctionPath="extinction.fits",
        arm="VIS",
        header=fits.Header({"DATE-OBS": "2024-01-02T03:04:05"}),
        recipeName="soxs-stare",
        startNightDate="2024-01-02",
        sofName="SYNTHETIC",
    )

    assert calibration.arm == "VIS"
    assert calibration.airmass == 1.2
    assert calibration.exptime == 30.0
    assert calibration.qcDir == str(tmp_path / "qc")
    assert calibration.productDir == str(tmp_path / "products")
    assert calibration.products.empty


def test_calibrate_writes_the_analytic_nir_flux_and_product_row(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    responsePath = tmp_path / "response.fits"
    Table(
        {
            "polyOrder": [1],
            "c0": [0.01],
            "c1": [2.0],
        }
    ).write(responsePath)
    originalHeader = fits.Header({"DATE-OBS": "2024-01-02T03:04:05.678"})
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame(
        {"WAVE": [500.0, 600.0], "FLUX_COUNTS": [20.0, 40.0]}
    )
    calibration.exptime = 10.0
    calibration.airmass = 1.0
    calibration.extinctionPath = ""
    calibration.arm = "NIR"
    calibration.debug = False
    calibration.header = originalHeader
    calibration.recipeName = "soxs-stare"
    calibration.sofName = "synthetic"
    calibration.productDir = str(tmp_path)
    calibration.products = pd.DataFrame()
    keywordNames = {
        "DPR_CATG": "ESO DPR CATG",
        "DPR_TYPE": "ESO DPR TYPE",
        "DET_READ_SPEED": "ESO DET READ CURID",
        "CONAD": "ESO DET OUT1 CONAD",
        "GAIN": "ESO DET OUT1 GAIN",
        "RON": "ESO DET OUT1 RON",
        "SEQ_ARM": "ESO SEQ ARM",
        "PRO_TYPE": "ESO PRO TYPE",
        "PRO_CATG": "ESO PRO CATG",
    }
    calibration.kw = keywordNames.__getitem__
    captured: dict[str, object] = {}

    def capture_write(**kwargs: object) -> None:
        captured.update(kwargs)

    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        capture_write,
    )

    outputPath, products = calibration.calibrate()

    assert outputPath == str(tmp_path / "synthetic_FLUXCAL.fits")
    outputTable = captured["tables"][0]
    expected = (
        calibration.extractedSpectrum["FLUX_COUNTS"]
        / calibration.exptime
        * np.polyval([0.01, 2.0], calibration.extractedSpectrum["WAVE"])
        * 1e-17
    )
    np.testing.assert_allclose(
        outputTable["FLUX_CALIBRATED"], expected, rtol=1e-14, atol=0
    )
    np.testing.assert_array_equal(outputTable["WAVE"], [500.0, 600.0])
    assert originalHeader == fits.Header({"DATE-OBS": "2024-01-02T03:04:05.678"})
    assert products.loc[0, "product_label"] == "EXTRACTED_FLUXCAL_SPECTRUM"
    assert products.loc[0, "file_path"] == outputPath


def test_calibrate_applies_vis_extinction_factor(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Apply the wavelength-dependent extinction correction for VIS spectra."""
    import soxspipe.commonutils.toolkit as toolkit

    responsePath = tmp_path / "response.fits"
    Table({"polyOrder": [0], "c0": [2.0]}).write(responsePath)
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame(
        {"WAVE": [500.0, 600.0], "FLUX_COUNTS": [20.0, 40.0]}
    )
    calibration.exptime = 10.0
    calibration.airmass = 1.4
    calibration.extinctionPath = "extinction.fits"
    calibration.arm = "VIS"
    calibration.debug = False
    calibration.header = fits.Header({"DATE-OBS": "2024-01-02T03:04:05"})
    calibration.recipeName = "soxs-stare"
    calibration.sofName = "synthetic"
    calibration.productDir = str(tmp_path)
    calibration.products = pd.DataFrame()
    calibration.kw = {
        "DPR_CATG": "ESO DPR CATG",
        "DPR_TYPE": "ESO DPR TYPE",
        "DET_READ_SPEED": "ESO DET READ CURID",
        "CONAD": "ESO DET OUT1 CONAD",
        "GAIN": "ESO DET OUT1 GAIN",
        "RON": "ESO DET OUT1 RON",
        "SEQ_ARM": "ESO SEQ ARM",
        "PRO_TYPE": "ESO PRO TYPE",
        "PRO_CATG": "ESO PRO CATG",
    }.__getitem__
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        toolkit,
        "extinction_correction_factor",
        lambda wavelengths, path, airmass: np.array([1.5, 2.0]),
    )
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.update(kwargs),
    )

    outputPath, _ = calibration.calibrate()

    outputTable = captured["tables"][0]
    np.testing.assert_allclose(
        outputTable["FLUX_CALIBRATED"], [6.0e-17, 16.0e-17], rtol=0, atol=0
    )
    assert outputPath == str(tmp_path / "synthetic_FLUXCAL.fits")
