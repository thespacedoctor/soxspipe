"""Integration contracts for calibrated extracted-spectrum products."""

from __future__ import annotations

import warnings
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
    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([7.0, 8.0]))
    originalHeader = fits.Header({"DATE-OBS": "2024-01-02T03:04:05.678"})
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame({"WAVE": [500.0, 600.0], "FLUX_COUNTS": [20.0, 40.0]})
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
        "PRO_AIRM_MEAN": "ESO PRO AIRM MEAN",
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
    expected = calibration.extractedSpectrum["FLUX_COUNTS"] / calibration.exptime * np.array([7.0, 8.0]) * 1e-17
    np.testing.assert_allclose(outputTable["FLUX_CALIBRATED"], expected, rtol=1e-14, atol=0)
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

    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([2.0, 2.0]))
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame({"WAVE": [500.0, 600.0], "FLUX_COUNTS": [20.0, 40.0]})
    calibration.exptime = 10.0
    calibration.airmass = 1.4
    calibration.extinctionPath = "extinction.fits"
    calibration.arm = "VIS"
    calibration.debug = False
    calibration.header = fits.Header({"DATE-OBS": "2024-01-02T03:04:05", "ESO DET BINY": 1})
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
        "PRO_AIRM_MEAN": "ESO PRO AIRM MEAN",
        "WIN_BINY": "ESO DET BINY",
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
    np.testing.assert_allclose(outputTable["FLUX_CALIBRATED"], [6.0e-17, 16.0e-17], rtol=1e-14, atol=0)
    assert outputPath == str(tmp_path / "synthetic_FLUXCAL.fits")


def test_calibrate_records_the_airmass_it_used_in_the_product_header(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The flux-calibrated product header records the mean airmass of the extinction correction (DY-1695)."""
    # ARRANGE
    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([2.0, 2.0]))
    calibration = _nir_calibration(tmp_path, log, responsePath)
    calibration.airmass = 1.37
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.update(kwargs),
    )

    # ACT
    calibration.calibrate()

    # ASSERT
    assert captured["header"]["HIERARCH ESO PRO AIRM MEAN"] == 1.37
    assert "ESO PRO AIRM MEAN" not in calibration.header


def test_calibrate_writes_no_airmass_card_when_no_airmass_is_known(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When no contributing frame has a readable airmass, the product header gets no airmass card (DY-1695)."""
    # ARRANGE
    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([2.0, 2.0]))
    calibration = _nir_calibration(tmp_path, log, responsePath)
    calibration.airmass = None
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.update(kwargs),
    )

    # ACT
    calibration.calibrate()

    # ASSERT
    assert "ESO PRO AIRM MEAN" not in captured["header"]


def test_calibrate_refuses_a_vis_spectrum_with_no_airmass(
    tmp_path: Path,
    log: object,
) -> None:
    """UVB and VIS need an airmass for the extinction correction, so none is a clear error (DY-1695)."""
    # ARRANGE
    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([2.0, 2.0]))
    calibration = _nir_calibration(tmp_path, log, responsePath)
    calibration.arm = "VIS"
    calibration.airmass = None

    # ACT / ASSERT
    with pytest.raises(ValueError, match="airmass"):
        calibration.calibrate()


def _nir_calibration(tmp_path: Path, log: object, responsePath: Path) -> flux_calibration:
    """Return a NIR flux_calibration that needs no extinction file."""
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame({"WAVE": [500.0, 600.0], "FLUX_COUNTS": [20.0, 40.0]})
    calibration.exptime = 10.0
    calibration.airmass = 1.0
    calibration.extinctionPath = ""
    calibration.arm = "NIR"
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
        "SEQ_ARM": "ESO SEQ ARM",
        "PRO_TYPE": "ESO PRO TYPE",
        "PRO_CATG": "ESO PRO CATG",
        "PRO_AIRM_MEAN": "ESO PRO AIRM MEAN",
    }.__getitem__
    return calibration


def _write_response_product(tmp_path: Path, log: object, wave: np.ndarray, response: np.ndarray) -> Path:
    """Write a RESP product with its RESPONSE curve extension."""
    from soxspipe.commonutils.response_function import response_function
    from tests.factories import product_table, qc_table

    writer = object.__new__(response_function)
    writer.log = log
    writer.arm = "NIR"
    writer.kw = lambda key: key
    writer.sofName = "SYNTHETIC"
    writer.header = fits.Header({"OBJECT": "Synthetic standard"})
    writer.productDir = str(tmp_path)
    writer.qc = qc_table().iloc[0:0]
    writer.products = product_table().iloc[0:0]
    writer.recipeName = "soxs-stare"
    writer.dateObs = "2024-01-02T03:04:05"
    return Path(writer.write_response_function_to_file(wavelength=wave, response=response))


def test_calibrate_uses_the_response_table_and_logs_no_fallback_warning(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A RESP product with the RESPONSE extension is calibrated from the table, not the polynomial (DY-1614)."""
    # ARRANGE
    responsePath = _write_response_product(tmp_path, log, np.array([500.0, 600.0]), np.array([252.0, 302.0]))
    captured: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.append(kwargs["tables"][0]),
    )

    # ACT
    with warnings.catch_warnings():
        warnings.filterwarnings("error", message=".*multiple tables are present.*")
        _nir_calibration(tmp_path, log, responsePath).calibrate()

    # ASSERT
    np.testing.assert_allclose(captured[0]["FLUX_CALIBRATED"].to_numpy(), [20 / 10 * 252e-17, 40 / 10 * 302e-17])
    assert not any(level == "warning" for level, _ in log.messages)


def test_calibrate_sets_samples_outside_the_response_table_to_nan_and_logs_the_count(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Science wavelengths beyond the response curve are NaN, never extrapolated (DY-1614)."""
    # ARRANGE
    responsePath = _write_response_product(tmp_path, log, np.array([520.0, 580.0]), np.array([100.0, 100.0]))
    captured: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.append(kwargs["tables"][0]),
    )

    # ACT
    _nir_calibration(tmp_path, log, responsePath).calibrate()

    # ASSERT: BOTH SCIENCE SAMPLES (500 AND 600 NM) LIE OUTSIDE 520-580 NM
    assert captured[0]["FLUX_CALIBRATED"].isna().all()
    assert any(level == "info" and "2 " in message for level, message in log.messages)


def test_calibrate_counts_only_science_samples_outside_the_valid_response_range(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Unusable response rows are dropped first, so only science wavelengths beyond the good rows are NaN (DY-1614)."""
    # ARRANGE: GOOD ROWS SPAN 520-580 NM; THE BAD ROWS SIT AT 480 NM, 540 NM AND 620 NM
    responsePath = _write_response_product(
        tmp_path,
        log,
        np.array([480.0, 520.0, 540.0, 560.0, 580.0, 620.0]),
        np.array([0.0, 100.0, -5.0, 100.0, 100.0, np.nan]),
    )
    captured: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.append(kwargs["tables"][0]),
    )

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        _nir_calibration(tmp_path, log, responsePath).calibrate()

    # ASSERT: THE SCIENCE SAMPLES AT 500 AND 600 NM ARE OUTSIDE 520-580 NM
    assert captured[0]["FLUX_CALIBRATED"].isna().all()
    assert any(level == "info" and message.startswith("2 science samples") for level, message in log.messages)
