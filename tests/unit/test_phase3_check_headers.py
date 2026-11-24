"""Unit tests for the keyword and header checks of the Phase 3 checker (DY-2345).

Every good product is written by the real Phase 3 table writer, then one fault is put in.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from soxspipe.commonutils.phase3_check import FAILURE, WARNING, check_phase3
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    add_cards,
    failure_messages,
    finding_messages,
    phase3_spectrum_file,
)

pytestmark = pytest.mark.unit


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
    failures = finding_messages(findings, FAILURE)
    assert any(keyword in message and "primary" in message for message in failures)
    assert finding_messages(findings, WARNING) == []


@pytest.mark.parametrize("keyword", ["NELEM", "VOCLASS", "TITLE", "TELAPSE", "SPEC_BW"])
def test_a_missing_mandatory_extension_keyword_is_a_failure(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, extension={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message and "extension" in message for message in finding_messages(findings, FAILURE))


@pytest.mark.parametrize("keyword", ["SPEC_ERR", "SPEC_SYE", "RA_ERR"])
def test_a_missing_recommended_keyword_is_a_warning(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message for message in finding_messages(findings, WARNING))
    assert finding_messages(findings, FAILURE) == []


@pytest.mark.parametrize("keyword", ["LAMRMS", "LAMNLIN"])
def test_a_missing_optional_keyword_is_a_warning(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={keyword: None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any(keyword in message for message in finding_messages(findings, WARNING))
    assert finding_messages(findings, FAILURE) == []


def test_a_mandatory_string_keyword_that_is_empty_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"OBJECT": ""})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert any("OBJECT" in message and "empty" in message for message in finding_messages(findings, FAILURE))


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
    assert any(keyword in message and "type" in message for message in finding_messages(findings, FAILURE))


def test_origfile_that_is_not_the_file_name_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    written = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    renamed = written.rename(tmp_path / "SOXS_20260102T030405_STARE_ph3_FAILED.fits")

    # ACT
    failures = finding_messages(check_phase3(renamed), FAILURE)

    # ASSERT
    assert any("ORIGFILE" in message and PRODUCT_NAME in message for message in failures)


def test_origfile_in_the_extension_header_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), extension={"ORIGFILE": PRODUCT_NAME})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ORIGFILE" in message and "extension" in message for message in failures)


def test_an_extname_in_the_primary_header_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={"EXTNAME": "SPECTRUM"})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "primary" in message for message in failures)


def test_a_table_not_named_spectrum_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), extension={"EXTNAME": "TABLE"})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "SPECTRUM" in message for message in failures)


@pytest.mark.parametrize("keyword", ["BUNIT", "SXSPRE", "SUBBIAS", "FLATCOR", "CCDPROC"])
@pytest.mark.parametrize("header", ["primary", "extension"])
def test_a_card_that_the_writer_drops_from_a_table_is_a_failure(
    tmp_path: Path, log: object, keyword: str, header: str
) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), **{header: {keyword: "x"}})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(keyword in message and header in message for message in failures)


def test_a_cd_matrix_card_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={"CD1_1": 1.0})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

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
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={keyword: value})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(keyword in message for message in failures)


@pytest.mark.parametrize("length", [69, 150])
def test_a_phase3_header_value_over_68_characters_is_a_failure(tmp_path: Path, log: object, length: int) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"OBSTECH": "X" * length})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("OBSTECH" in message for message in failures)


def test_an_internal_spectrum_may_hold_a_value_over_68_characters(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / "soxs_stare_SPECTRUM.fits", log, primary={"OBSTECH": "X" * 150})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

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
    failures = failure_messages(tmp_path, log, primary=primary)

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
    failures = failure_messages(tmp_path, log, extension={keyword: other})

    # ASSERT
    assert any(keyword in message and "differs" in message for message in failures)


def test_texptime_over_telapse_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, primary={"TEXPTIME": 700.0})

    # ASSERT
    assert any("TEXPTIME" in message and "TELAPSE" in message for message in failures)


def test_texptime_equal_to_telapse_is_not_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, primary={"TEXPTIME": 691.2})

    # ASSERT
    assert failures == []
