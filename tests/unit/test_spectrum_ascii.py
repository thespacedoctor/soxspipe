"""The ASCII spectrum is written once, from the flux-calibrated FITS table when there is one (DY-1361)."""

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

MERGED_HEADER = "WAVE FLUX_COUNTS VARIANCE SKY_COUNTS SNR"
FLUXCAL_HEADER = "WAVE FLUX_CALIBRATED ERR_CALIBRATED FLUX_COUNTS VARIANCE SKY_COUNTS SNR"


def _write_merged_fits(directory: Path, *, withSky: bool = True) -> Path:
    """Write a merged-spectrum FITS table and return its path."""
    sky = [3.0, 5.0] if withSky else None
    table = Table(merged_spectrum([500.0, 500.02], [20.0, 40.0], [16.0, 64.0], sky, asTable=True))
    path = directory / "STEM_EXTRACTED_MERGED.fits"
    table.write(path, format="fits")
    return path


def _write_fluxcal_fits(directory: Path) -> Path:
    """Write a flux-calibrated FITS table with the seven columns and return its path."""
    table = Table(
        {
            "WAVE": [500.0, 500.02],
            "FLUX_CALIBRATED": [1.4e-16, 3.2e-16],
            "ERR_CALIBRATED": [2.8e-17, 6.4e-17],
            "FLUX_COUNTS": [20.0, 40.0],
            "VARIANCE": [16.0, 64.0],
            "SKY_COUNTS": [3.0, 5.0],
            "SNR": [5.0, 5.0],
        }
    )
    path = directory / "STEM_FLUXCAL.fits"
    table.write(path, format="fits")
    return path


def test_ascii_path_is_the_fluxcal_name_when_the_spectrum_was_calibrated(tmp_path):
    # ARRANGE
    merged = tmp_path / "STEM_EXTRACTED_MERGED.fits"
    fluxcal = tmp_path / "STEM_FLUXCAL.fits"

    # ACT
    path = ascii_spectrum_path(str(merged), str(fluxcal))

    # ASSERT
    assert path == str(tmp_path / "STEM_FLUXCAL.txt")


def test_ascii_path_keeps_the_merged_name_without_a_fluxcal_table(tmp_path):
    # ARRANGE
    merged = tmp_path / "STEM_EXTRACTED_MERGED.fits"

    # ACT
    path = ascii_spectrum_path(str(merged), None)

    # ASSERT
    assert path == str(tmp_path / "STEM_EXTRACTED_MERGED.txt")


def test_ascii_is_written_from_the_fluxcal_table_when_one_exists(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)
    fluxcal = _write_fluxcal_fits(tmp_path)

    # ACT
    products = write_ascii_spectrum_product(str(merged), str(fluxcal), product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT: THE HEADER HOLDS THE SEVEN COLUMNS AND WAVE IS IN ANGSTROM
    lines = (tmp_path / "STEM_FLUXCAL.txt").read_text().splitlines()
    assert lines[0] == FLUXCAL_HEADER
    assert lines[1].split()[0] == "5000.00"
    assert not (tmp_path / "STEM_EXTRACTED_MERGED.txt").exists()
    assert products.loc[0, "file_name"] == "STEM_FLUXCAL.txt"
    assert products.loc[0, "product_desc"] == "ASCII version of the flux-calibrated extracted spectrum"


def test_ascii_falls_back_to_the_merged_spectrum_without_fluxcal(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)

    # ACT
    products = write_ascii_spectrum_product(str(merged), None, product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT
    lines = (tmp_path / "STEM_EXTRACTED_MERGED.txt").read_text().splitlines()
    assert lines[0] == MERGED_HEADER
    assert lines[1].split()[0] == "5000.00"
    assert products.loc[0, "file_name"] == "STEM_EXTRACTED_MERGED.txt"
    assert products.loc[0, "product_desc"] == "ASCII version of the extracted spectrum (not flux calibrated)"


def test_ascii_of_a_spectrum_with_no_sky_has_no_sky_column(tmp_path):
    # ARRANGE: A NOD OR OFFSET SPECTRUM
    merged = _write_merged_fits(tmp_path, withSky=False)

    # ACT
    write_ascii_spectrum_product(str(merged), None, product_table().iloc[0:0], "soxs-nod", "2024")

    # ASSERT
    header = (tmp_path / "STEM_EXTRACTED_MERGED.txt").read_text().splitlines()[0]
    assert header == "WAVE FLUX_COUNTS VARIANCE SNR"


def test_ascii_records_one_product_row_with_the_merged_ascii_label(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)
    existing = product_table()

    # ACT
    products = write_ascii_spectrum_product(str(merged), None, existing, "soxs-stare", "2024-01-02T03:04:05")

    # ASSERT
    newRow = products.iloc[-1]
    assert len(products) == len(existing) + 1
    assert len(existing) == len(product_table())
    assert newRow["product_label"] == "EXTRACTED_MERGED_ASCII"
    assert newRow["file_type"] == "TXT"
    assert newRow["label"] == "PROD"
    assert newRow["obs_date_utc"] == "2024-01-02T03:04:05"
    assert newRow["file_path"] == str(tmp_path / "STEM_EXTRACTED_MERGED.txt")


def test_ascii_leaves_the_source_fits_table_unchanged(tmp_path):
    # ARRANGE
    merged = _write_merged_fits(tmp_path)
    before = Table.read(merged)

    # ACT
    write_ascii_spectrum_product(str(merged), None, pd.DataFrame(), "soxs-stare", "2024")

    # ASSERT: THE FITS KEEPS ITS WAVELENGTHS IN NANOMETRES
    assert list(Table.read(merged)["WAVE"]) == list(before["WAVE"])
    assert Table.read(merged)["WAVE"][0] == 500.0


def test_ascii_keeps_a_nan_sky_as_nan_not_zero(tmp_path):
    # ARRANGE: A FLUX-CALIBRATED TABLE WHOSE EXTRACTION HAD NO SKY MODEL
    fluxcal = Table(
        {
            "WAVE": [500.0],
            "FLUX_CALIBRATED": [1.4e-16],
            "ERR_CALIBRATED": [2.8e-17],
            "FLUX_COUNTS": [20.0],
            "VARIANCE": [16.0],
            "SKY_COUNTS": [float("nan")],
            "SNR": [5.0],
        }
    )
    fluxcalPath = tmp_path / "STEM_FLUXCAL.fits"
    fluxcal.write(fluxcalPath, format="fits")
    mergedPath = _write_merged_fits(tmp_path)

    # ACT
    write_ascii_spectrum_product(str(mergedPath), str(fluxcalPath), product_table().iloc[0:0], "soxs-stare", "2024")

    # ASSERT
    row = (tmp_path / "STEM_FLUXCAL.txt").read_text().splitlines()[1].split()
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
