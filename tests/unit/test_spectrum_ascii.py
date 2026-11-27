"""The ASCII spectrum is written once, from the ``_SPECTRUM`` FITS table of the pass (DY-1361, DY-2346)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from astropy.table import MaskedColumn, Table

from soxspipe.commonutils.spectrum_ascii import (
    ascii_spectrum_path,
    write_ascii_spectrum,
    write_ascii_spectrum_product,
)
from tests.factories import merged_spectrum, product_table

pytestmark = pytest.mark.unit

MERGED_HEADER = "WAVE FLUX_COUNTS ERR_COUNTS SKY_COUNTS SNR"
FLUXCAL_HEADER = "WAVE FLUX ERR FLUX_COUNTS ERR_COUNTS SKY_COUNTS SNR"


def _write_merged_fits(directory: Path, *, withSky: bool = True) -> Path:
    """Write a merged-spectrum FITS table and return its path."""
    sky = [3.0, 5.0] if withSky else None
    table = Table(merged_spectrum([500.0, 500.02], [20.0, 40.0], [16.0, 64.0], sky, asTable=True))
    path = directory / "STEM_SPECTRUM.fits"
    table.write(path, format="fits")
    return path


def _write_fluxcal_fits(directory: Path) -> Path:
    """Write a flux-calibrated FITS table with the seven columns and return its path."""
    table = Table(
        {
            "WAVE": [500.0, 500.02],
            "FLUX": [1.4e-16, 3.2e-16],
            "ERR": [2.8e-17, 6.4e-17],
            "FLUX_COUNTS": [20.0, 40.0],
            "ERR_COUNTS": [4.0, 8.0],
            "SKY_COUNTS": [3.0, 5.0],
            "SNR": [5.0, 5.0],
        }
    )
    path = directory / "STEM_SPECTRUM.fits"
    table.write(path, format="fits")
    return path


def test_ascii_path_is_the_spectrum_name_with_a_txt_extension(tmp_path):
    # ARRANGE
    spectrum = tmp_path / "STEM_SPECTRUM.fits"

    # ACT
    path = ascii_spectrum_path(str(spectrum))

    # ASSERT
    assert path == str(tmp_path / "STEM_SPECTRUM.txt")


def test_ascii_of_a_flux_calibrated_spectrum_holds_the_seven_columns(tmp_path):
    # ARRANGE
    spectrum = _write_fluxcal_fits(tmp_path)

    # ACT
    products = write_ascii_spectrum_product(str(spectrum), True, product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT: THE HEADER HOLDS THE SEVEN COLUMNS AND WAVE IS IN ANGSTROM
    lines = (tmp_path / "STEM_SPECTRUM.txt").read_text().splitlines()
    assert lines[0] == FLUXCAL_HEADER
    assert lines[1].split()[0] == "5000.00"
    assert products.loc[0, "file_name"] == "STEM_SPECTRUM.txt"
    assert products.loc[0, "product_desc"] == "ASCII version of the flux-calibrated extracted spectrum"


def test_ascii_of_an_uncalibrated_spectrum_is_in_counts(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)

    # ACT
    products = write_ascii_spectrum_product(str(merged), False, product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT
    lines = (tmp_path / "STEM_SPECTRUM.txt").read_text().splitlines()
    assert lines[0] == MERGED_HEADER
    assert lines[1].split()[0] == "5000.00"
    assert products.loc[0, "file_name"] == "STEM_SPECTRUM.txt"
    assert products.loc[0, "product_desc"] == "ASCII version of the extracted spectrum (not flux calibrated)"


def test_ascii_of_a_spectrum_with_no_sky_has_no_sky_column(tmp_path):
    # ARRANGE: A NOD OR OFFSET SPECTRUM
    merged = _write_merged_fits(tmp_path, withSky=False)

    # ACT
    write_ascii_spectrum_product(str(merged), False, product_table().iloc[0:0], "soxs-nod", "2024")

    # ASSERT
    header = (tmp_path / "STEM_SPECTRUM.txt").read_text().splitlines()[0]
    assert header == "WAVE FLUX_COUNTS ERR_COUNTS SNR"


def test_ascii_records_one_product_row_with_the_spectrum_ascii_label(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)
    existing = product_table()

    # ACT
    products = write_ascii_spectrum_product(str(merged), False, existing, "soxs-stare", "2024-01-02T03:04:05")

    # ASSERT
    newRow = products.iloc[-1]
    assert len(products) == len(existing) + 1
    assert len(existing) == len(product_table())
    assert newRow["product_label"] == "SPECTRUM_ASCII"
    assert newRow["file_type"] == "TXT"
    assert newRow["label"] == "PROD"
    assert newRow["obs_date_utc"] == "2024-01-02T03:04:05"
    assert newRow["file_path"] == str(tmp_path / "STEM_SPECTRUM.txt")


def test_ascii_leaves_the_source_fits_table_unchanged(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)
    before = Table.read(merged)

    # ACT
    write_ascii_spectrum_product(str(merged), False, pd.DataFrame(), "soxs-stare", "2024")

    # ASSERT: THE FITS KEEPS ITS WAVELENGTHS IN NANOMETRES
    assert list(Table.read(merged)["WAVE"]) == list(before["WAVE"])
    assert Table.read(merged)["WAVE"][0] == 500.0


def test_ascii_keeps_a_nan_sky_as_nan_not_zero(tmp_path):
    # ARRANGE: A FLUX-CALIBRATED TABLE WHOSE EXTRACTION HAD NO SKY MODEL
    fluxcal = Table(
        {
            "WAVE": [500.0],
            "FLUX": [1.4e-16],
            "ERR": [2.8e-17],
            "FLUX_COUNTS": [20.0],
            "ERR_COUNTS": [4.0],
            "SKY_COUNTS": [float("nan")],
            "SNR": [5.0],
        }
    )
    fluxcalPath = tmp_path / "STEM_SPECTRUM.fits"
    fluxcal.write(fluxcalPath, format="fits")

    # ACT
    write_ascii_spectrum_product(str(fluxcalPath), True, product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT
    row = (tmp_path / "STEM_SPECTRUM.txt").read_text().splitlines()[1].split()
    assert row[FLUXCAL_HEADER.split().index("SKY_COUNTS")].lower() == "nan"


def test_ascii_writes_a_masked_cell_as_nan_not_blank(tmp_path, monkeypatch):
    # ARRANGE: THE FITS READER HANDS BACK A MASKED CELL
    maskedTable = Table(
        {
            "WAVE": [500.0, 500.02],
            "FLUX_COUNTS": MaskedColumn([20.0, 1.0e20], mask=[False, True]),
            "SNR": [5.0, 5.0],
        }
    )
    monkeypatch.setattr(Table, "read", classmethod(lambda cls, *args, **kwargs: maskedTable))

    # ACT
    write_ascii_spectrum(str(tmp_path / "ignored.fits"), str(tmp_path / "STEM.txt"))

    # ASSERT
    lines = (tmp_path / "STEM.txt").read_text().splitlines()
    assert lines[2].split() == ["5000.20", "nan", "5.0"]
