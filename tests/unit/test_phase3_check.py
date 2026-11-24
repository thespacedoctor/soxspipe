"""Unit tests for the Phase 3 compliance checker `soxspipe.commonutils.phase3_check` (DY-2345).

Every good product is written by the real Phase 3 table writer, then one fault is put in to test one rule.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import sys
from collections.abc import Callable
from importlib import import_module
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from docopt import docopt as parse_docopt

from soxspipe.commonutils import phase3_check, phase3_check_framework
from soxspipe.commonutils.phase3_check import FAILURE, WARNING, check_phase3
from tests.factories.phase3_products import SPECTRUM_LENGTH, phase3_spectrum_file

pytestmark = pytest.mark.unit

PRODUCT_NAME = "SOXS_20260102T030405_STARE_ph3.fits"


def _messages(findings: list, severity: str) -> list[str]:
    return [finding.message for finding in findings if finding.severity == severity]


def _checks(findings: list, severity: str) -> list[str]:
    return [finding.check for finding in findings if finding.severity == severity]


def test_a_good_product_has_no_failures_and_no_warnings(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == []


@pytest.mark.parametrize("keyword", ["PRODCATG", "OBJECT", "RA", "DEC", "EXPTIME", "MJD-OBS", "FLUXCAL", "SPEC_RES"])
def test_a_missing_mandatory_primary_keyword_is_a_failure(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    failures = _messages(findings, FAILURE)
    assert any(keyword in message and "primary" in message for message in failures)
    assert _messages(findings, WARNING) == []


@pytest.mark.parametrize("keyword", ["NELEM", "VOCLASS", "TITLE", "TELAPSE", "SPEC_BW"])
def test_a_missing_mandatory_extension_keyword_is_a_failure(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, extension={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message and "extension" in message for message in _messages(findings, FAILURE))


@pytest.mark.parametrize("keyword", ["SPEC_ERR", "SPEC_SYE", "RA_ERR"])
def test_a_missing_recommended_keyword_is_a_warning(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message for message in _messages(findings, WARNING))
    assert _messages(findings, FAILURE) == []


@pytest.mark.parametrize("keyword", ["LAMRMS", "LAMNLIN"])
def test_a_missing_optional_keyword_is_a_warning(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message for message in _messages(findings, WARNING))
    assert _messages(findings, FAILURE) == []


def test_a_mandatory_string_keyword_that_is_empty_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"OBJECT": ""})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("OBJECT" in message and "empty" in message for message in _messages(findings, FAILURE))


def test_an_empty_referenc_is_not_a_finding(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert not any("REFERENC" in finding.message for finding in findings)


@pytest.mark.parametrize(
    ("keyword", "value"),
    [("EXPTIME", "long"), ("NCOMBINE", 2.5), ("TOT_FLUX", "F"), ("OBJECT", 7), ("OBID1", True)],
)
def test_a_value_of_the_wrong_kind_is_a_failure(tmp_path: Path, log: object, keyword: str, value: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: value})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message and "type" in message for message in _messages(findings, FAILURE))


def test_a_file_that_is_not_fits_gives_one_failure(tmp_path: Path) -> None:
    # ARRANGE
    notFits = tmp_path / PRODUCT_NAME
    notFits.write_text("this is not a FITS file")

    # ACT
    findings = check_phase3(notFits)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]
    assert "cannot be opened" in findings[0].message


@pytest.mark.parametrize("name", ["", "."])
def test_a_path_that_is_not_a_file_gives_one_failure(
    tmp_path: Path, name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    monkeypatch.chdir(tmp_path)

    # ACT
    findings = check_phase3(name)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


def test_a_truncated_fits_file_gives_one_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    product.write_bytes(product.read_bytes()[:100])

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


def test_a_missing_file_gives_one_failure(tmp_path: Path) -> None:
    # ARRANGE
    missing = tmp_path / PRODUCT_NAME

    # ACT
    findings = check_phase3(missing)

    # ASSERT
    assert [finding.severity for finding in findings] == [FAILURE]


def _rewrite(path: Path, change: Callable[[fits.HDUList], fits.HDUList | None]) -> Path:
    """Read a written product, let ``change`` break it, and write it back with fresh checksums."""
    with fits.open(path) as hdus:
        broken = change(hdus) or hdus
        broken.writeto(path, checksum=True, overwrite=True)
    return path


def _variable_length_table(hdus: fits.HDUList) -> fits.HDUList:
    columns = fits.ColDefs([fits.Column(name="WAVE", format="PD()", array=[np.array([600.0, 600.5])])])
    return fits.HDUList([hdus[0], fits.BinTableHDU.from_columns(columns, name="SPECTRUM")])


def test_a_second_extension_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    _rewrite(product, lambda hdus: hdus.append(fits.BinTableHDU(name="OTHER")))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("exactly one extension" in message for message in _messages(findings, FAILURE))


def test_data_in_the_primary_hdu_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def add_data(hdus: fits.HDUList) -> None:
        hdus[0].data = np.zeros((2, 2))

    _rewrite(product, add_data)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("primary HDU" in message and "NAXIS" in message for message in _messages(findings, FAILURE))


def test_an_extension_that_is_not_a_binary_table_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def image_extension(hdus: fits.HDUList) -> fits.HDUList:
        return fits.HDUList([hdus[0], fits.ImageHDU(np.zeros(3), name="SPECTRUM")])

    _rewrite(product, image_extension)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("BINTABLE" in message for message in _messages(findings, FAILURE))


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
    assert any("NAXIS2" in message for message in _messages(findings, FAILURE))


def test_a_variable_length_array_column_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    _rewrite(product, _variable_length_table)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("variable-length" in message for message in _messages(findings, FAILURE))


def test_a_corrupted_header_fails_the_checksum(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    product.write_bytes(product.read_bytes().replace(b"SYNTHETIC TARGET", b"SYNTHETIC TARGEX", 1))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("CHECKSUM" in message and "primary" in message for message in _messages(findings, FAILURE))


def test_a_corrupted_table_fails_the_datasum(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    firstWave, changedWave = (np.array([value], dtype=">f8").tobytes() for value in (600.0, 601.0))
    product.write_bytes(product.read_bytes().replace(firstWave, changedWave, 1))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("DATASUM" in message and "extension" in message for message in _messages(findings, FAILURE))


def test_a_file_name_that_does_not_end_in_dot_fits_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / "SOXS_20260102T030405_STARE_ph3.fit", log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(".fits" in message for message in _messages(findings, FAILURE))


def test_a_phase3_file_name_over_68_characters_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    longName = "SOXS_" + "X" * 60 + "_ph3.fits"
    product = phase3_spectrum_file(tmp_path / longName, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("68" in message and "file name" in message for message in _messages(findings, FAILURE))


def test_an_internal_spectrum_file_name_over_68_characters_is_allowed(tmp_path: Path, log: object) -> None:
    # ARRANGE
    longName = "soxs_stare_" + "x" * 80 + "_SPECTRUM.fits"
    product = phase3_spectrum_file(tmp_path / longName, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert not any("file name" in message for message in _messages(findings, FAILURE))


def _failures(tmp_path: Path, log: object, **factoryArguments: object) -> list[str]:
    return _messages(check_phase3(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, **factoryArguments)), FAILURE)


def test_the_first_three_columns_must_be_wave_flux_and_err_in_order(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, dropColumns=("ERR",))

    # ASSERT
    assert any("WAVE, FLUX, ERR" in message for message in failures)


def test_a_column_of_another_length_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columns={"FLUX_COUNTS": np.ones(SPECTRUM_LENGTH - 1)})

    # ASSERT
    assert any("FLUX_COUNTS" in message and "elements" in message for message in failures)


def test_nelem_that_differs_from_the_array_length_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, extension={"NELEM": SPECTRUM_LENGTH + 1})

    # ASSERT
    assert any("NELEM" in message for message in failures)


@pytest.mark.parametrize("wave", [np.array([600.0, 600.5, 600.5, 601.0, 601.5, 602.0]), np.arange(6.0)[::-1] + 600])
def test_wave_that_does_not_strictly_increase_is_a_failure(tmp_path: Path, log: object, wave: np.ndarray) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columns={"WAVE": wave})

    # ASSERT
    assert any("WAVE" in message and "strictly increas" in message for message in failures)


@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_a_non_finite_value_in_a_column_is_a_failure(tmp_path: Path, log: object, value: float) -> None:
    # ARRANGE
    flux = np.full(SPECTRUM_LENGTH, 1.0e-16)
    flux[2] = value

    # ACT
    failures = _failures(tmp_path, log, columns={"FLUX": flux})

    # ASSERT
    assert any("FLUX" in message and "NaN" in message for message in failures)


def test_a_float_qual_column_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columns={"QUAL": np.zeros(SPECTRUM_LENGTH)})

    # ASSERT
    assert any("QUAL" in message and "integer" in message for message in failures)


def test_a_repeated_column_name_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def repeat_name(hdus: fits.HDUList) -> None:
        hdus[1].header["TTYPE6"] = "FLUX_COUNTS"

    _rewrite(product, repeat_name)

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("FLUX_COUNTS" in message and "more than once" in message for message in failures)


def test_a_wave_ucd_that_is_neither_air_nor_vacuum_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(
        tmp_path, log, columnCards={"WAVE": ("nm", "em.wl;obs.vacuum", "spec:Data.SpectralAxis.Value")}
    )

    # ASSERT
    assert any("TUCD1" in message and "em.wl" in message for message in failures)


def test_a_vacuum_wave_ucd_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columnCards={"WAVE": ("nm", "em.wl", "spec:Data.SpectralAxis.Value")})

    # ASSERT
    assert failures == []


@pytest.mark.parametrize(
    ("column", "cards", "expected"),
    [
        ("FLUX", ("erg / (s cm2 Angstrom)", "phot.flux.density;em.wl", "spec:Data.FluxAxis.Value"), "meta.main"),
        (
            "ERR",
            ("erg / (s cm2 Angstrom)", "stat.error;phot.flux.density;em.ql;meta.main", "spec:Data.FluxAxis.Value"),
            "TUTYP3",
        ),
        ("WAVE", ("nm", "em.wl;obs.atmos", "Spectrum.Data.SpectralAxis.Value"), "TUTYP1"),
    ],
)
def test_the_main_columns_need_the_v2_utypes_and_the_main_flux_ucds(
    tmp_path: Path, log: object, column: str, cards: tuple[str, str, str], expected: str
) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columnCards={column: cards})

    # ASSERT
    assert any(expected in message for message in failures)


@pytest.mark.parametrize(
    ("cards", "expected"),
    [
        ((None, "phot.flux.density;em.wl", "eso:Data.FluxAxis.Value"), "TUNIT5"),
        (("adu", None, "eso:Data.FluxAxis.Value"), "TUCD5"),
        (("adu", "phot.flux.density;em.wl", None), "TUTYP5"),
    ],
)
def test_a_column_without_unit_ucd_or_utype_is_a_failure(
    tmp_path: Path, log: object, cards: tuple[str | None, str | None, str | None], expected: str
) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, columnCards={"FLUX_COUNTS": cards})

    # ASSERT
    assert any(expected in message and "missing" in message for message in failures)


def test_err_in_another_unit_than_flux_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(
        tmp_path,
        log,
        columnCards={
            "ERR": ("adu", "stat.error;phot.flux.density;em.ql;meta.main", "spec:Data.FluxAxis.Accuracy.StatError")
        },
    )

    # ASSERT
    assert any("ERR" in message and "unit of FLUX" in message for message in failures)


def test_a_bgflux_column_must_match_the_flux_unit_and_ucd(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(
        tmp_path,
        log,
        columns={"BGFLUX": np.ones(SPECTRUM_LENGTH)},
        columnCards={"BGFLUX": ("adu", "phot.flux.density;em.wl", "eso:Data.FluxAxis.Value")},
    )

    # ASSERT
    assert any("BGFLUX" in message and "FLUX" in message for message in failures)


def test_a_missing_tdmin_or_tdmax_of_wave_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, extension={"TDMIN1": None})

    # ASSERT
    assert any("TDMIN1" in message for message in failures)


def test_wave_outside_tdmin_and_tdmax_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, extension={"TDMAX1": 601.0})

    # ASSERT
    assert any("TDMAX1" in message and "WAVE" in message for message in failures)


def test_an_exposure_column_needs_its_ucd_and_an_empty_utype(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(
        tmp_path,
        log,
        columns={"EXPOSURE": np.full(SPECTRUM_LENGTH, 600.0)},
        columnCards={"EXPOSURE": ("s", "time.duration", "eso:Data.Exposure")},
    )

    # ASSERT
    assert any("EXPOSURE" in message and "time.duration;obs.exposure" in message for message in failures)
    assert any("EXPOSURE" in message and "empty" in message for message in failures)


def test_a_good_exposure_column_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(
        tmp_path,
        log,
        columns={"EXPOSURE": np.full(SPECTRUM_LENGTH, 600.0)},
        columnCards={"EXPOSURE": ("s", "time.duration;obs.exposure", "")},
    )

    # ASSERT
    assert failures == []


def _add_cards(
    path: Path, primary: dict[str, object] | None = None, extension: dict[str, object] | None = None
) -> Path:
    def add(hdus: fits.HDUList) -> None:
        hdus[0].header.update(primary or {})
        hdus[1].header.update(extension or {})

    return _rewrite(path, add)


def test_origfile_that_is_not_the_file_name_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    written = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    renamed = written.rename(tmp_path / "SOXS_20260102T030405_STARE_ph3_FAILED.fits")

    # ACT
    failures = _messages(check_phase3(renamed), FAILURE)

    # ASSERT
    assert any("ORIGFILE" in message and PRODUCT_NAME in message for message in failures)


def test_origfile_in_the_extension_header_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), extension={"ORIGFILE": PRODUCT_NAME})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ORIGFILE" in message and "extension" in message for message in failures)


def test_an_extname_in_the_primary_header_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={"EXTNAME": "SPECTRUM"})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "primary" in message for message in failures)


def test_a_table_not_named_spectrum_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), extension={"EXTNAME": "TABLE"})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "SPECTRUM" in message for message in failures)


@pytest.mark.parametrize("keyword", ["BUNIT", "SXSPRE", "SUBBIAS", "FLATCOR", "CCDPROC"])
@pytest.mark.parametrize("header", ["primary", "extension"])
def test_a_card_that_the_writer_drops_from_a_table_is_a_failure(
    tmp_path: Path, log: object, keyword: str, header: str
) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), **{header: {keyword: "x"}})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(keyword in message and header in message for message in failures)


def test_a_cd_matrix_card_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={"CD1_1": 1.0})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("CD1_1" in message for message in failures)


@pytest.mark.parametrize(
    ("keyword", "value"),
    [("ORIGIN", "NOT-ESO"), ("PROCSOFT", "somepipeline 1.0")],
)
def test_the_identification_cards_must_have_the_values_the_writer_sets(
    tmp_path: Path, log: object, keyword: str, value: str
) -> None:
    # ARRANGE
    product = _add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={keyword: value})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(keyword in message for message in failures)


@pytest.mark.parametrize("length", [69, 150])
def test_a_phase3_header_value_over_68_characters_is_a_failure(tmp_path: Path, log: object, length: int) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"OBSTECH": "X" * length})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("OBSTECH" in message for message in failures)


def test_an_internal_spectrum_may_hold_a_value_over_68_characters(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / "soxs_stare_SPECTRUM.fits", log, primary={"OBSTECH": "X" * 150})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert failures == []


@pytest.mark.parametrize(
    ("primary", "expected"),
    [
        ({"PRODCATG": "SCIENCE.IMAGE"}, "PRODCATG"),
        ({"RADESYS": "GALACTIC"}, "RADESYS"),
        ({"RADESYS": "FK5", "EQUINOX": 1950.0}, "EQUINOX"),
        ({"RADESYS": "FK5", "EQUINOX": None}, "EQUINOX"),
        ({"RADESYS": "ICRS", "EQUINOX": 1950.0}, "EQUINOX"),
        ({"TIMESYS": "TT"}, "TIMESYS"),
        ({"FLUXCAL": "CALIBRATED"}, "FLUXCAL"),
        ({"FLUXCAL": "UNCALIBRATED", "FLUXERR": -2.0}, "FLUXERR"),
        ({"SPECSYS": "TOPO"}, "SPECSYS"),
        ({"FLUXCAL": "UNCALIBRATED", "FLUXERR": -1.0, "TOT_FLUX": True}, "TOT_FLUX"),
    ],
)
def test_a_value_the_standard_does_not_allow_is_a_failure(
    tmp_path: Path, log: object, primary: dict, expected: str
) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, primary=primary)

    # ASSERT
    assert any(expected in message for message in failures)


@pytest.mark.parametrize(
    "primary",
    [
        {"RADESYS": "ICRS", "EQUINOX": None},
        {"RADESYS": "ICRS", "EQUINOX": 2000.0},
        {"TIMESYS": "TAI"},
        {"TIMESYS": "UTC"},
        {"FLUXCAL": "UNCALIBRATED", "FLUXERR": -1.0},
        {"SPECSYS": "TOPOCENT"},
    ],
)
def test_a_value_the_standard_allows_is_not_a_finding(tmp_path: Path, log: object, primary: dict) -> None:
    # ARRANGE, ACT
    findings = check_phase3(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=primary))

    # ASSERT
    assert findings == []


@pytest.mark.parametrize("keyword", ["OBJECT", "RA", "DEC"])
def test_a_primary_value_that_differs_from_the_extension_is_a_failure(
    tmp_path: Path, log: object, keyword: str
) -> None:
    # ARRANGE
    other = "OTHER TARGET" if keyword == "OBJECT" else 10.0
    failures = _failures(tmp_path, log, extension={keyword: other})

    # ASSERT
    assert any(keyword in message and "differs" in message for message in failures)


def test_texptime_over_telapse_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, primary={"TEXPTIME": 700.0})

    # ASSERT
    assert any("TEXPTIME" in message and "TELAPSE" in message for message in failures)


def test_texptime_equal_to_telapse_is_not_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = _failures(tmp_path, log, primary={"TEXPTIME": 691.2})

    # ASSERT
    assert failures == []


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
    failures = _messages(check_phase3(product), FAILURE)

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
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(expected in message for message in failures)


def test_an_associated_fits_file_needs_an_ancillary_prodcatg(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path, fitsCards={"PRODCATG": "SCIENCE.IMAGE"})
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "ANCILLARY" in message for message in failures)


def test_an_associated_fits_file_may_not_link_files_in_turn(tmp_path: Path, log: object) -> None:
    # ARRANGE
    links = _ancillary_files(tmp_path, fitsCards={"ASSON1": PREVIEW_NAME})
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=links)

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSON2" in message and "nested" in message for message in failures)


def test_an_assom_without_an_asson_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"ASSOM3": "0" * 32})

    # ACT
    failures = _messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ASSOM3" in message and "ASSON3" in message for message in failures)


@pytest.fixture
def isolated_registry(monkeypatch: pytest.MonkeyPatch) -> dict:
    """Give a test its own copy of the check registry, so a check it registers is gone afterwards."""
    registry = dict(phase3_check_framework.CHECKS)
    monkeypatch.setattr(phase3_check_framework, "CHECKS", registry)
    return registry


def test_a_registered_check_runs_and_its_finding_carries_the_check_name(
    tmp_path: Path, log: object, isolated_registry: dict
) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_object_is_not_a_test_target(product: phase3_check.Phase3Product) -> list:
        return [phase3_check.fail(f"{product.primary['OBJECT']} is a test target")]

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == [
        phase3_check.Finding(FAILURE, "SYNTHETIC TARGET is a test target", "check_object_is_not_a_test_target")
    ]


def test_a_registered_check_may_return_a_warning(tmp_path: Path, log: object, isolated_registry: dict) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_for_a_doubt(product: phase3_check.Phase3Product) -> list:
        return [phase3_check.warn("doubtful")]

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [(finding.severity, finding.message) for finding in findings] == [(WARNING, "doubtful")]


def test_registering_two_checks_with_one_name_is_an_error(isolated_registry: dict) -> None:
    # ARRANGE
    def check_something(product: phase3_check.Phase3Product) -> list:
        return []

    phase3_check.register_check(check_something)

    # ACT, ASSERT
    with pytest.raises(ValueError, match="check_something"):
        phase3_check.register_check(check_something)


def test_a_check_that_raises_is_a_failure_and_the_other_checks_still_run(
    tmp_path: Path, log: object, isolated_registry: dict
) -> None:
    # ARRANGE
    @phase3_check.register_check
    def check_that_crashes(product: phase3_check.Phase3Product) -> list:
        raise KeyError("EXPTIME")

    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"PRODCATG": None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    crash = [finding for finding in findings if finding.check == "check_that_crashes"]
    assert len(crash) == 1
    assert crash[0].severity == FAILURE
    assert "KeyError" in crash[0].message
    assert any("PRODCATG" in finding.message for finding in findings)


def test_a_table_check_is_skipped_when_the_file_has_no_table(tmp_path: Path, isolated_registry: dict) -> None:
    # ARRANGE
    calls: list[str] = []

    @phase3_check.register_check(needsTable=True)
    def check_that_reads_the_table(product: phase3_check.Phase3Product) -> list:
        calls.append("ran")
        return []

    emptyFile = tmp_path / PRODUCT_NAME
    fits.HDUList([fits.PrimaryHDU()]).writeto(emptyFile)

    # ACT
    findings = check_phase3(emptyFile)

    # ASSERT
    assert calls == []
    assert any("exactly one extension" in finding.message for finding in findings)


@pytest.mark.parametrize("extension", [None, "image"])
def test_no_check_crashes_on_a_file_without_a_binary_table(tmp_path: Path, extension: str | None) -> None:
    # ARRANGE
    product = tmp_path / PRODUCT_NAME
    hdus = [fits.PrimaryHDU()] + ([fits.ImageHDU(np.zeros(3), name="SPECTRUM")] if extension else [])
    fits.HDUList(hdus).writeto(product)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert not [finding for finding in findings if "crashed" in finding.message]
    assert not [finding for finding in findings if finding.check in ("check_column_names", "check_column_descriptors")]


def test_the_module_docstring_examples_register_and_run(tmp_path: Path, log: object, isolated_registry: dict) -> None:
    # ARRANGE
    namespace: dict = {}
    exampleStart = phase3_check.__doc__.index("```python") + len("```python")
    exec(phase3_check.__doc__[exampleStart : phase3_check.__doc__.index("```", exampleStart)], namespace)  # noqa: S102
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, columns={"FLUX": np.zeros(SPECTRUM_LENGTH)})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert [(finding.check, finding.severity) for finding in findings] == [("check_flux_is_not_all_zero", WARNING)]


def test_the_command_prints_each_file_and_returns_0_when_there_is_no_failure(
    tmp_path: Path, log: object, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"LAMRMS": None})

    # ACT
    status = phase3_check.check_phase3_command([str(product)])

    # ASSERT
    output = capsys.readouterr().out
    assert status == 0
    assert f"{product}: PASS (0 failures, 1 warning)" in output
    assert "warning [check_mandatory_and_recommended_keywords] optional keyword LAMRMS is missing" in output


def test_the_command_returns_1_when_any_file_has_a_failure(
    tmp_path: Path, log: object, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})

    # ACT
    status = phase3_check.check_phase3_command([str(good), str(bad)])

    # ASSERT
    output = capsys.readouterr().out
    assert status == 1
    assert f"{good}: PASS (0 failures, 0 warnings)" in output
    assert f"{bad}: FAIL (1 failure, 0 warnings)" in output
    assert "failure [check_mandatory_and_recommended_keywords] mandatory keyword PRODCATG is missing" in output


def test_the_command_returns_1_for_a_file_that_does_not_exist(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE, ACT
    status = phase3_check.check_phase3_command([str(tmp_path / "missing.fits")])

    # ASSERT
    assert status == 1
    assert "FAIL" in capsys.readouterr().out


def _run_soxspipe(monkeypatch: pytest.MonkeyPatch, *arguments: str) -> int:
    clUtils = import_module("soxspipe.cl_utils")
    monkeypatch.setattr(sys, "argv", ["soxspipe", *arguments])
    with pytest.raises(SystemExit) as exit_info:
        clUtils.main()
    return exit_info.value.code


def test_soxspipe_check_phase3_exits_0_for_a_good_product_and_writes_nothing(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    products = tmp_path / "products"
    products.mkdir()
    product = phase3_spectrum_file(products / PRODUCT_NAME, log)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    monkeypatch.chdir(workspace)

    # ACT
    code = _run_soxspipe(monkeypatch, "check-phase3", str(product))

    # ASSERT
    assert code == 0
    assert "PASS" in capsys.readouterr().out
    assert list(workspace.iterdir()) == []


def test_soxspipe_check_phase3_exits_1_when_any_of_several_files_fails(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})
    monkeypatch.chdir(tmp_path)

    # ACT
    code = _run_soxspipe(monkeypatch, "check-phase3", str(good), str(bad))

    # ASSERT
    assert code == 1
    output = capsys.readouterr().out
    assert "PASS" in output
    assert "FAIL" in output


def test_the_usage_text_that_fundamentals_reads_has_only_names_it_can_unpack() -> None:
    # ARRANGE
    clUtils = import_module("soxspipe.cl_utils")

    # ACT
    arguments = parse_docopt(clUtils.TOOLS_DOC, argv=["prep"])

    # ASSERT
    unpackedNames = [
        argument.replace("-", "") + "Flag" if argument[0] == "-" else argument.replace("<", "").replace(">", "")
        for argument in arguments
    ]
    assert all(name.isidentifier() for name in unpackedNames)
    assert "check-phase3" not in arguments


def test_the_full_usage_text_still_has_the_check_phase3_command() -> None:
    # ARRANGE
    clUtils = import_module("soxspipe.cl_utils")

    # ACT
    arguments = parse_docopt(clUtils.__doc__, argv=["check-phase3", "a.fits", "b.fits"])

    # ASSERT
    assert arguments["check-phase3"] is True
    assert arguments["<file>"] == ["a.fits", "b.fits"]


def test_the_real_command_line_exits_1_for_a_bad_product_and_0_for_a_good_one(tmp_path: Path, log: object) -> None:
    # ARRANGE
    good = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    bad = phase3_spectrum_file(tmp_path / "SOXS_20260102T030406_STARE_ph3.fits", log, primary={"PRODCATG": None})
    driver = "import sys; sys.argv = ['soxspipe', *sys.argv[1:]]; from soxspipe.cl_utils import main; main()"
    packageRoot = str(Path(phase3_check.__file__).resolve().parents[2])

    def run(*files: Path) -> subprocess.CompletedProcess:
        return subprocess.run(  # noqa: S603
            [sys.executable, "-W", "ignore", "-c", driver, "check-phase3", *map(str, files)],
            cwd=tmp_path,
            env={**os.environ, "PYTHONPATH": packageRoot},
            capture_output=True,
            text=True,
            check=False,
        )

    # ACT
    goodRun, badRun = run(good), run(good, bad)

    # ASSERT
    assert (goodRun.returncode, badRun.returncode) == (0, 1)
    assert "PASS" in goodRun.stdout
    assert "mandatory keyword PRODCATG is missing" in badRun.stdout
    assert goodRun.stderr == ""
