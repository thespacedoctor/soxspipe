"""Unit tests for the file, HDU-structure, checksum and ancillary-link checks of the Phase 3 checker (DY-2345)."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.phase3_check import FAILURE, check_phase3
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    SPECTRUM_LENGTH,
    finding_messages,
    phase3_spectrum_file,
    rewrite_product,
)

pytestmark = pytest.mark.unit


def _variable_length_table(hdus: fits.HDUList) -> fits.HDUList:
    columns = fits.ColDefs([fits.Column(name="WAVE", format="PD()", array=[np.array([600.0, 600.5])])])
    return fits.HDUList([hdus[0], fits.BinTableHDU.from_columns(columns, name="SPECTRUM")])


def test_a_second_extension_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    rewrite_product(product, lambda hdus: hdus.append(fits.BinTableHDU(name="OTHER")))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("exactly one extension" in message for message in finding_messages(findings, FAILURE))


def test_data_in_the_primary_hdu_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def add_data(hdus: fits.HDUList) -> None:
        hdus[0].data = np.zeros((2, 2))

    rewrite_product(product, add_data)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("primary HDU" in message and "NAXIS" in message for message in finding_messages(findings, FAILURE))


def test_an_extension_that_is_not_a_binary_table_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def image_extension(hdus: fits.HDUList) -> fits.HDUList:
        return fits.HDUList([hdus[0], fits.ImageHDU(np.zeros(3), name="SPECTRUM")])

    rewrite_product(product, image_extension)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("BINTABLE" in message for message in finding_messages(findings, FAILURE))


def test_a_table_with_two_rows_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    twoRows = np.zeros((2, SPECTRUM_LENGTH))
    product = phase3_spectrum_file(
        tmp_path / PRODUCT_NAME,
        log,
        columns={name: twoRows + 1 + np.arange(SPECTRUM_LENGTH) for name in ("WAVE", "FLUX", "ERR")},
    )

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("NAXIS2" in message for message in finding_messages(findings, FAILURE))


def test_a_variable_length_array_column_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    rewrite_product(product, _variable_length_table)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("variable-length" in message for message in finding_messages(findings, FAILURE))


def test_a_corrupted_header_fails_the_checksum(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    product.write_bytes(product.read_bytes().replace(b"SYNTHETIC TARGET", b"SYNTHETIC TARGEX", 1))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("CHECKSUM" in message and "primary" in message for message in finding_messages(findings, FAILURE))


def test_a_corrupted_table_fails_the_datasum(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    firstWave, changedWave = (np.array([value], dtype=">f8").tobytes() for value in (600.0, 601.0))
    product.write_bytes(product.read_bytes().replace(firstWave, changedWave, 1))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("DATASUM" in message and "extension" in message for message in finding_messages(findings, FAILURE))


def test_a_file_name_that_does_not_end_in_dot_fits_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / "SOXS_20260102T030405_STARE_ph3.fit", log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(".fits" in message for message in finding_messages(findings, FAILURE))


def test_a_phase3_file_name_over_68_characters_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    longName = "SOXS_" + "X" * 60 + "_ph3.fits"
    product = phase3_spectrum_file(tmp_path / longName, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("68" in message and "file name" in message for message in finding_messages(findings, FAILURE))


def test_an_internal_spectrum_file_name_over_68_characters_is_allowed(tmp_path: Path, log: object) -> None:
    # ARRANGE
    longName = "soxs_stare_" + "x" * 80 + "_SPECTRUM.fits"
    product = phase3_spectrum_file(tmp_path / longName, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert not any("file name" in message for message in finding_messages(findings, FAILURE))


PREVIEW_NAME = "SOXS_20260102T030405_STARE_ph3_preview.png"
TWO_D_NAME = "SOXS_20260102T030405_STARE_ph3_2D.fits"
PREVIEW_BYTES = b"not really a png"


def _ancillary_files(folder: Path, *, fitsCards: dict[str, object] | None = None) -> dict[str, object]:
    """Write a preview and a 2D spectrum beside the product, and return the primary cards that link them."""
    (folder / PREVIEW_NAME).write_bytes(PREVIEW_BYTES)
    fits.PrimaryHDU(np.zeros((2, 2)), fits.Header({"PRODCATG": "ANCILLARY.2DSPECTRUM", **(fitsCards or {})})).writeto(
        folder / TWO_D_NAME
    )
    return {
        "ASSON1": PREVIEW_NAME,
        "ASSOC1": "ANCILLARY.PREVIEW",
        "ASSOM1": hashlib.md5(PREVIEW_BYTES, usedforsecurity=False).hexdigest(),
        "ASSON2": TWO_D_NAME,
    }


def test_a_product_with_a_preview_and_a_2d_spectrum_linked_is_not_a_finding(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path)
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == []


def test_an_associated_file_that_is_not_beside_the_product_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path)
    (tmp_path / TWO_D_NAME).unlink()
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "not found" in message for message in failures)


@pytest.mark.parametrize(
    ("change", "expected"),
    [
        ({"ASSOC1": None}, "ASSOC1"),
        ({"ASSOM1": None}, "ASSOM1"),
        ({"ASSOM1": "0" * 32}, "ASSOM1"),
        ({"ASSOC1": "PREVIEW"}, "ASSOC1"),
    ],
)
def test_a_non_fits_associated_file_needs_its_category_and_md5(
    tmp_path: Path, log: object, change: dict[str, object], expected: str
) -> None:
    # ARRANGE
    links = {**_ancillary_files(tmp_path), **change}
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(expected in message for message in failures)


def test_an_associated_fits_file_needs_an_ancillary_prodcatg(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path, fitsCards={"PRODCATG": "SCIENCE.IMAGE"})
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "ANCILLARY" in message for message in failures)


def test_an_associated_fits_file_may_not_link_files_in_turn(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path, fitsCards={"ASSON1": PREVIEW_NAME})
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "nested" in message for message in failures)


def test_an_assom_without_an_asson_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"ASSOM3": "0" * 32})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSOM3" in message and "ASSON3" in message for message in failures)


def _links_to(name: str, folder: Path) -> dict[str, object]:
    """Return the primary cards that link one non-FITS file by the given ``ASSON1`` value."""
    return {"ASSON1": name, "ASSOC1": "ANCILLARY.PREVIEW", "ASSOM1": "0" * 32}


@pytest.mark.parametrize("name", ["../outside.png", "sub/outside.png", ".", "..", "/etc/hosts"])
def test_an_asson_that_is_not_a_bare_file_name_is_a_failure_that_names_the_value(
    tmp_path: Path, log: object, name: str
) -> None:
    # ARRANGE
    folder = tmp_path / "batch"
    folder.mkdir()
    (tmp_path / "outside.png").write_bytes(PREVIEW_BYTES)
    product = phase3_spectrum_file(folder / PRODUCT_NAME, log, primary=_links_to(name, folder))

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON1" in message and repr(name) in message and "bare file name" in message for message in failures)


def test_an_asson_that_is_a_symlink_to_a_file_outside_the_folder_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    folder = tmp_path / "batch"
    folder.mkdir()
    outside = tmp_path / "secret.png"
    outside.write_bytes(PREVIEW_BYTES)
    (folder / PREVIEW_NAME).symlink_to(outside)
    links = {
        "ASSON1": PREVIEW_NAME,
        "ASSOC1": "ANCILLARY.PREVIEW",
        "ASSOM1": hashlib.md5(PREVIEW_BYTES, usedforsecurity=False).hexdigest(),
    }
    product = phase3_spectrum_file(folder / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON1" in message and repr(PREVIEW_NAME) in message and "outside" in message for message in failures)


def test_an_asson_that_is_a_symlink_to_a_file_in_the_same_folder_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path)
    real = tmp_path / "real_preview.png"
    (tmp_path / PREVIEW_NAME).rename(real)
    (tmp_path / PREVIEW_NAME).symlink_to(real)
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == []


def test_the_md5_of_a_linked_file_is_read_in_chunks_not_all_at_once(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path)
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    def refuse(self: Path) -> bytes:
        raise AssertionError("the whole file was read into memory")

    monkeypatch.setattr(Path, "read_bytes", refuse)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == []


def test_a_linked_fits_file_that_cannot_be_read_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path)
    (tmp_path / TWO_D_NAME).write_text("this is not FITS")
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "cannot be opened" in message for message in failures)


def test_a_failed_phase3_file_name_over_68_characters_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / ("SOXS_" + "X" * 55 + "_ph3_FAILED.fits"), log)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("68" in message and "file name" in message for message in failures)
