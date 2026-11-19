"""The columns the flux-calibrated table carries from the merged spectrum (DY-1361)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.flux_calibration import (
    _calculate_flux_calibration,
    _response_factors_from_table,
    flux_calibration,
)
from soxspipe.commonutils.response_basis import RESPONSE_BASIS, RESPONSE_BASIS_KEYWORD
from tests.factories import merged_spectrum

pytestmark = pytest.mark.unit

WAVES = [500.0, 600.0]
FLUX = [20.0, 40.0]
VARIANCE = [16.0, 64.0]
SKY = [3.0, 5.0]
RESPONSE = [7.0, 8.0]
EXPOSURE_TIME = 10.0
KEYWORDS = {
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
}


def _write_response(tmp_path: Path, wave: list[float], response: list[float]) -> Path:
    """Write a RESP product holding a RESPONSE curve with the current basis card."""
    curve = fits.table_to_hdu(Table({"WAVE": wave, "RESPONSE": response}))
    curve.name = "RESPONSE"
    primary = fits.PrimaryHDU(header=fits.Header({RESPONSE_BASIS_KEYWORD: RESPONSE_BASIS}))
    path = tmp_path / "SYNTHETIC_RESP.fits"
    fits.HDUList([primary, curve]).writeto(path)
    return path


def _calibrate(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    spectrum: object,
    *,
    arm: str = "NIR",
    binning: int = 1,
    responseWave: list[float] | None = None,
    extinction: np.ndarray | None = None,
) -> Table:
    """Run ``calibrate`` on ``spectrum`` and return the table it hands to the writer."""
    import soxspipe.commonutils.toolkit as toolkit

    # A RESPONSE COVERING OTHER WAVELENGTHS IS FLAT, SO ONLY ITS RANGE MATTERS
    responsePath = (
        _write_response(tmp_path, WAVES, RESPONSE)
        if responseWave is None
        else _write_response(tmp_path, responseWave, [7.0, 7.0])
    )
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = spectrum
    calibration.exptime = EXPOSURE_TIME
    calibration.airmass = 1.4
    calibration.extinctionPath = "extinction.fits"
    calibration.arm = arm
    calibration.debug = False
    calibration.header = fits.Header({"DATE-OBS": "2024-01-02T03:04:05", "ESO DET BINY": binning})
    calibration.recipeName = "soxs-stare"
    calibration.sofName = "synthetic"
    calibration.productDir = str(tmp_path)
    calibration.products = pd.DataFrame()
    calibration.kw = KEYWORDS.__getitem__
    captured: list[object] = []
    monkeypatch.setattr(
        "soxspipe.commonutils.phase3.write_fits_table_to_disk",
        lambda **kwargs: captured.append(kwargs["tables"][0]),
    )
    if extinction is not None:
        monkeypatch.setattr(toolkit, "extinction_correction_factor", lambda wavelengths, path, airmass: extinction)
    calibration.calibrate()
    return captured[0]


def test_fluxcal_table_has_the_seven_columns_in_order(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    assert list(table.columns) == [
        "WAVE",
        "FLUX_CALIBRATED",
        "ERR_CALIBRATED",
        "FLUX_COUNTS",
        "VARIANCE",
        "SKY_COUNTS",
        "SNR",
    ]


def test_fluxcal_table_has_no_sky_column_when_the_spectrum_has_none(tmp_path, log, monkeypatch):
    # ARRANGE: A NOD OR OFFSET SPECTRUM CARRIES NO SKY MODEL
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, asTable=True)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    assert list(table.columns) == ["WAVE", "FLUX_CALIBRATED", "ERR_CALIBRATED", "FLUX_COUNTS", "VARIANCE", "SNR"]


def test_fluxcal_flux_is_bit_identical_to_the_previous_formula(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)
    responseFactors = _response_factors_from_table(WAVES, WAVES, RESPONSE)
    expected = _calculate_flux_calibration(np.asarray(FLUX, dtype=float), EXPOSURE_TIME, responseFactors, 1.0)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    np.testing.assert_array_equal(np.asarray(table["FLUX_CALIBRATED"]), expected)


def test_fluxcal_error_is_sqrt_variance_times_the_flux_factor(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: SQRT(16) = 4 AND SQRT(64) = 8 COUNTS, CALIBRATED LIKE THE FLUX
    np.testing.assert_allclose(table["ERR_CALIBRATED"], [4 / 10 * 7e-17, 8 / 10 * 8e-17], rtol=1e-14, atol=0)


def test_fluxcal_error_includes_the_extinction_factor_like_the_flux(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum, arm="VIS", extinction=np.array([1.5, 2.0]))

    # ASSERT
    np.testing.assert_allclose(table["ERR_CALIBRATED"], [4 / 10 * 1.5 * 7e-17, 8 / 10 * 2.0 * 8e-17], rtol=1e-14)
    np.testing.assert_allclose(
        np.asarray(table["ERR_CALIBRATED"]) / np.asarray(table["FLUX_CALIBRATED"]),
        np.sqrt(VARIANCE) / FLUX,
        rtol=1e-6,
    )


def test_fluxcal_error_is_scaled_by_dispersion_binning_like_the_flux(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum, arm="VIS", binning=2, extinction=np.array([1.0, 1.0]))

    # ASSERT: THE COUNTS AND THE ERROR BOTH PER UNBINNED PIXEL, SO HALF THE UNBINNED-BASIS VALUES
    np.testing.assert_allclose(table["ERR_CALIBRATED"], [4 / 2 / 10 * 7e-17, 8 / 2 / 10 * 8e-17], rtol=1e-14)
    np.testing.assert_allclose(table["FLUX_CALIBRATED"], [20 / 2 / 10 * 7e-17, 40 / 2 / 10 * 8e-17], rtol=1e-14)


def test_fluxcal_copies_counts_variance_sky_and_snr_unchanged(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: THE SKY IS NOT CALIBRATED; IT STAYS IN E-/PIX
    np.testing.assert_array_equal(np.asarray(table["FLUX_COUNTS"]), FLUX)
    np.testing.assert_array_equal(np.asarray(table["VARIANCE"]), VARIANCE)
    np.testing.assert_array_equal(np.asarray(table["SKY_COUNTS"]), SKY)
    np.testing.assert_array_equal(np.asarray(table["SNR"]), np.asarray(FLUX) / np.sqrt(VARIANCE))


def test_fluxcal_flux_and_error_are_nan_outside_the_response_but_snr_and_sky_are_not(tmp_path, log, monkeypatch):
    # ARRANGE: THE RESPONSE COVERS ONLY 550-560 NM, SO BOTH SAMPLES LIE OUTSIDE IT
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum, responseWave=[550.0, 560.0])

    # ASSERT
    assert np.isnan(np.asarray(table["FLUX_CALIBRATED"])).all()
    assert np.isnan(np.asarray(table["ERR_CALIBRATED"])).all()
    assert np.isfinite(np.asarray(table["SNR"])).all()
    assert np.isfinite(np.asarray(table["SKY_COUNTS"])).all()


def test_fluxcal_carries_nan_sky_counts_unchanged_when_the_extraction_had_no_sky_model(tmp_path, log, monkeypatch):
    # ARRANGE: A STARE FRAME WITH NO SKY MODEL HAS NAN SKY_COUNTS IN THE MERGED SPECTRUM (DY-1948)
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, [np.nan, np.nan])

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: NAN STAYS NAN, NOT 0, AND THE FLUX AND ERROR ARE UNAFFECTED
    assert np.isnan(np.asarray(table["SKY_COUNTS"])).all()
    assert np.isfinite(np.asarray(table["FLUX_CALIBRATED"])).all()
    assert np.isfinite(np.asarray(table["ERR_CALIBRATED"])).all()


@pytest.mark.parametrize("asTable", [False, True], ids=["quantity-dataframe", "astropy-table"])
def test_fluxcal_accepts_a_quantity_dataframe_and_an_astropy_table(tmp_path, log, monkeypatch, asTable):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY, asTable=asTable)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    np.testing.assert_allclose(table["ERR_CALIBRATED"], [4 / 10 * 7e-17, 8 / 10 * 8e-17], rtol=1e-14)
    np.testing.assert_array_equal(np.asarray(table["SKY_COUNTS"]), SKY)


def test_fluxcal_refuses_a_spectrum_without_variance(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY).drop(columns=["VARIANCE"])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="VARIANCE"):
        _calibrate(tmp_path, log, monkeypatch, spectrum)


def test_fluxcal_error_column_is_float64(tmp_path, log, monkeypatch):
    # ARRANGE: A FAINT SOURCE'S CALIBRATED VARIANCE IS NEAR 1E-38 AND WOULD UNDERFLOW FLOAT32
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = _calibrate(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    assert np.asarray(table["ERR_CALIBRATED"]).dtype == np.float64
