"""The columns the flux-calibrated table carries from the merged spectrum (DY-1361)."""

from __future__ import annotations


import numpy as np
import pytest
from astropy.io import fits
from astropy.table import MaskedColumn, Table

from soxspipe.commonutils.flux_calibration import (
    _calculate_flux_calibration,
    _response_factors_from_table,
)
from tests.factories import merged_spectrum
from tests.factories.flux_calibration import (
    EXPOSURE_TIME,
    RESPONSE,
    WAVES,
    calibrate_and_capture_write,
    calibrate_spectrum,
)

pytestmark = pytest.mark.unit

FLUX = [20.0, 40.0]
VARIANCE = [16.0, 64.0]
SKY = [3.0, 5.0]


def test_fluxcal_table_has_the_seven_soxs_columns_in_order(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    assert list(table.columns) == [
        "WAVE",
        "FLUX",
        "ERR",
        "FLUX_COUNTS",
        "ERR_COUNTS",
        "SKY_COUNTS",
        "SNR",
    ]


def test_fluxcal_table_has_no_sky_column_when_the_spectrum_has_none(tmp_path, log, monkeypatch):
    # ARRANGE: A NOD OR OFFSET SPECTRUM CARRIES NO SKY MODEL
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, asTable=True)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    assert list(table.columns) == ["WAVE", "FLUX", "ERR", "FLUX_COUNTS", "ERR_COUNTS", "SNR"]


def test_fluxcal_flux_is_bit_identical_to_the_previous_formula(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)
    responseFactors = _response_factors_from_table(WAVES, WAVES, RESPONSE)
    expected = _calculate_flux_calibration(np.asarray(FLUX, dtype=float), EXPOSURE_TIME, responseFactors, 1.0)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    np.testing.assert_array_equal(np.asarray(table["FLUX"]), expected)


def test_fluxcal_error_is_the_counts_error_times_the_flux_factor(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: THE ERRORS ARE 4 AND 8 COUNTS (THE ROOTS OF THE VARIANCES 16 AND 64), CALIBRATED LIKE THE FLUX
    np.testing.assert_allclose(table["ERR"], [4 / 10 * 7e-17, 8 / 10 * 8e-17], rtol=1e-14, atol=0)


def test_fluxcal_error_includes_the_extinction_factor_like_the_flux(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, arm="VIS", extinction=np.array([1.5, 2.0]))

    # ASSERT
    np.testing.assert_allclose(table["ERR"], [4 / 10 * 1.5 * 7e-17, 8 / 10 * 2.0 * 8e-17], rtol=1e-14)
    np.testing.assert_allclose(
        np.asarray(table["ERR"]) / np.asarray(table["FLUX"]),
        np.sqrt(VARIANCE) / FLUX,
        rtol=1e-6,
    )


def test_fluxcal_error_is_scaled_by_dispersion_binning_like_the_flux(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(
        tmp_path, log, monkeypatch, spectrum, arm="VIS", binning=2, extinction=np.array([1.0, 1.0])
    )

    # ASSERT: THE COUNTS AND THE ERROR BOTH PER UNBINNED PIXEL, SO HALF THE UNBINNED-BASIS VALUES
    np.testing.assert_allclose(table["ERR"], [4 / 2 / 10 * 7e-17, 8 / 2 / 10 * 8e-17], rtol=1e-14)
    np.testing.assert_allclose(table["FLUX"], [20 / 2 / 10 * 7e-17, 40 / 2 / 10 * 8e-17], rtol=1e-14)


@pytest.mark.parametrize(
    ("arm", "binning", "extinction"),
    [("NIR", 1, None), ("VIS", 1, np.array([1.5, 2.0])), ("VIS", 2, np.array([1.5, 2.0]))],
    ids=["nir", "vis-extinction", "vis-binned-extinction"],
)
def test_fluxcal_err_is_err_counts_scaled_by_exactly_the_factors_of_the_flux(
    tmp_path, log, monkeypatch, arm, binning, extinction
):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, arm=arm, binning=binning, extinction=extinction)

    # ASSERT: FLUX / FLUX_COUNTS AND ERR / ERR_COUNTS ARE THE SAME FACTOR AT EVERY WAVELENGTH
    fluxFactor = np.asarray(table["FLUX"]) / np.asarray(table["FLUX_COUNTS"], dtype=float)
    errFactor = np.asarray(table["ERR"]) / np.asarray(table["ERR_COUNTS"], dtype=float)
    np.testing.assert_allclose(errFactor, fluxFactor, rtol=1e-14, atol=0)


def test_fluxcal_keeps_the_variance_layout_of_an_xshooter_spectrum(tmp_path, log, monkeypatch):
    # ARRANGE: THE X-SHOOTER MERGE STILL WRITES VARIANCE, AND ITS PRODUCT KEEPS THE COLUMNS OF BEFORE DY-2350
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY, hasVarianceColumn=True)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, settings={"instrument": "xsh"})

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
    np.testing.assert_allclose(table["ERR_CALIBRATED"], [4 / 10 * 7e-17, 8 / 10 * 8e-17], rtol=1e-14, atol=0)
    np.testing.assert_array_equal(np.asarray(table["VARIANCE"]), VARIANCE)


def test_fluxcal_refuses_a_soxs_run_given_a_variance_spectrum(tmp_path, log, monkeypatch):
    # ARRANGE: A SOXS RUN HOLDING THE LAYOUT OF BEFORE DY-2350, SUCH AS A STALE FILE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY, hasVarianceColumn=True)

    # ACT / ASSERT: NO TABLE WITH X-SHOOTER NAMES IS WRITTEN FOR A SOXS RUN
    with pytest.raises(ValueError, match=r"holds VARIANCE, but SOXS spectra hold ERR_COUNTS"):
        calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, settings={"instrument": "soxs"})


def test_fluxcal_refuses_an_xshooter_run_given_an_err_counts_spectrum(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=r"holds ERR_COUNTS, but X-Shooter spectra hold VARIANCE"):
        calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, settings={"instrument": "xsh"})


def test_fluxcal_copies_counts_error_sky_and_snr_unchanged(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: THE SKY IS NOT CALIBRATED; IT STAYS IN E-/PIX
    np.testing.assert_array_equal(np.asarray(table["FLUX_COUNTS"]), FLUX)
    np.testing.assert_array_equal(np.asarray(table["ERR_COUNTS"]), np.sqrt(VARIANCE))
    np.testing.assert_array_equal(np.asarray(table["SKY_COUNTS"]), SKY)
    np.testing.assert_array_equal(np.asarray(table["SNR"]), np.asarray(FLUX) / np.sqrt(VARIANCE))


def test_fluxcal_flux_and_error_are_nan_outside_the_response_but_snr_and_sky_are_not(tmp_path, log, monkeypatch):
    # ARRANGE: THE RESPONSE COVERS ONLY 550-560 NM, SO BOTH SAMPLES LIE OUTSIDE IT
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum, responseWave=[550.0, 560.0])

    # ASSERT
    assert np.isnan(np.asarray(table["FLUX"])).all()
    assert np.isnan(np.asarray(table["ERR"])).all()
    assert np.isfinite(np.asarray(table["SNR"])).all()
    assert np.isfinite(np.asarray(table["SKY_COUNTS"])).all()


def test_fluxcal_carries_nan_sky_counts_unchanged_when_the_extraction_had_no_sky_model(tmp_path, log, monkeypatch):
    # ARRANGE: A STARE FRAME WITH NO SKY MODEL HAS NAN SKY_COUNTS IN THE MERGED SPECTRUM (DY-1948)
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, [np.nan, np.nan])

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT: NAN STAYS NAN, NOT 0, AND THE FLUX AND ERROR ARE UNAFFECTED
    assert np.isnan(np.asarray(table["SKY_COUNTS"])).all()
    assert np.isfinite(np.asarray(table["FLUX"])).all()
    assert np.isfinite(np.asarray(table["ERR"])).all()


def test_fluxcal_turns_masked_cells_into_nan_instead_of_exposing_their_fill_values(tmp_path, log, monkeypatch):
    # ARRANGE: A TABLE WITH A MASKED SKY CELL, WHOSE UNDERLYING DATA IS A FILL VALUE
    spectrum = Table(
        {
            "WAVE": WAVES,
            "FLUX_COUNTS": np.asarray(FLUX, dtype=np.float32),
            "ERR_COUNTS": np.sqrt(np.asarray(VARIANCE, dtype=np.float32)),
            "SKY_COUNTS": MaskedColumn(np.asarray([3.0, 1.0e20], dtype=np.float32), mask=[False, True]),
            "SNR": np.asarray(FLUX) / np.sqrt(VARIANCE),
        }
    )

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    sky = np.asarray(table["SKY_COUNTS"])
    assert sky[0] == 3.0
    assert np.isnan(sky[1])
    assert sky.dtype == np.float32


@pytest.mark.parametrize("asTable", [False, True], ids=["pandas-dataframe", "astropy-table"])
def test_fluxcal_accepts_a_pandas_dataframe_and_an_astropy_table(tmp_path, log, monkeypatch, asTable):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY, asTable=asTable)

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    np.testing.assert_allclose(table["ERR"], [4 / 10 * 7e-17, 8 / 10 * 8e-17], rtol=1e-14)
    np.testing.assert_array_equal(np.asarray(table["SKY_COUNTS"]), SKY)


@pytest.mark.parametrize("asTable", [False, True], ids=["pandas-dataframe", "astropy-table"])
def test_fluxcal_keeps_the_float32_dtype_of_the_counts_error_and_sky(tmp_path, log, monkeypatch, asTable):
    # ARRANGE: THE MERGE STORES COUNTS, ERROR AND SKY AS FLOAT32, AND THE COPIED COLUMNS MUST NOT BE WIDENED
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY, asTable=asTable)
    assert np.asarray(spectrum["SKY_COUNTS"]).dtype == np.float32

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    for name in ("FLUX_COUNTS", "ERR_COUNTS", "SKY_COUNTS"):
        assert np.asarray(table[name]).dtype == np.float32, name
    np.testing.assert_array_equal(np.asarray(table["SKY_COUNTS"]), SKY)


def test_fluxcal_refuses_a_spectrum_without_an_error_column(tmp_path, log, monkeypatch):
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY).drop(columns=["ERR_COUNTS"])

    # ACT / ASSERT
    with pytest.raises(ValueError, match="ERR_COUNTS"):
        calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)


def test_fluxcal_error_column_is_float64_and_does_not_underflow(tmp_path, log, monkeypatch):
    # ARRANGE: AN ERROR OF 1E-30 COUNTS GIVES A CALIBRATED ERROR OF ABOUT 7E-48, WHICH FLOAT32 WOULD ROUND TO ZERO
    spectrum = Table(
        {
            "WAVE": WAVES,
            "FLUX_COUNTS": FLUX,
            "ERR_COUNTS": [1.0e-30, 2.0e-30],
            "SNR": [1.0, 1.0],
        }
    )

    # ACT
    table = calibrate_spectrum(tmp_path, log, monkeypatch, spectrum)

    # ASSERT
    error = np.asarray(table["ERR"])
    assert error.dtype == np.float64
    np.testing.assert_allclose(error, [1.0e-30 / 10 * 7e-17, 2.0e-30 / 10 * 8e-17], rtol=1e-12)
    assert (error > 0).all()


def test_the_calibrated_spectrum_header_says_absolute_and_keeps_the_target_position(tmp_path, log, monkeypatch):
    # ARRANGE: THE HEADER IS THE ONE THE RECIPE STAMPED, WITH AN UNCALIBRATED STATUS FROM THE COUNTS FILE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, SKY)
    header = fits.Header(
        {"DATE-OBS": "2024-01-02T03:04:05", "ESO DET BINY": 1, "RA": 7.5825, "DEC": -33.2503, "FLUXCAL": "UNCALIBRATED"}
    )
    header["FLUXERR"] = -1.0

    # ACT
    written = calibrate_and_capture_write(tmp_path, log, monkeypatch, spectrum, header=header)

    # ASSERT
    stamped = written["header"]
    assert (stamped["FLUXCAL"], stamped["TOT_FLUX"], stamped["FLUXERR"]) == ("ABSOLUTE", False, -2.0)
    assert (stamped["RA"], stamped["DEC"]) == (7.5825, -33.2503)
    assert header["FLUXCAL"] == "UNCALIBRATED"
