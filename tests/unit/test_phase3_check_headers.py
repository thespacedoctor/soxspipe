"""Unit tests for the keyword and header checks of the Phase 3 checker (DY-2345).

Every good product is written by the real Phase 3 table writer, then one fault is put in.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from astropy.io import fits

from soxspipe.commonutils.phase3 import RENAMED_ESO_CARDS, SCRUBBED_ESO_CARDS
from soxspipe.commonutils.phase3_check import FAILURE, WARNING, check_phase3
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    add_cards,
    failure_messages,
    finding_messages,
    phase3_spectrum_file,
    rewrite_product,
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


def test_a_missing_tcomm_is_a_warning_for_each_column(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    rewrite_product(product, lambda hdus: hdus[1].header.remove("TCOMM2") or hdus[1].header.remove("TCOMM4"))

    # ACT
    findings = check_phase3(product)

    # ASSERT
    warnings = finding_messages(findings, WARNING)
    assert sorted(message.split()[2] for message in warnings) == ["TCOMM2", "TCOMM4"]
    assert all("optional" in message and "extension" in message for message in warnings)
    assert finding_messages(findings, FAILURE) == []


def test_a_multi_prog_id_needs_progid1(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, primary={"PROG_ID": "MULTI"})

    # ASSERT
    assert any("PROGID1" in message and "MULTI" in message for message in failures)


def test_a_multi_prog_id_with_progid1_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, primary={"PROG_ID": "MULTI", "PROGID1": "115.ABCD.001"})

    # ASSERT
    assert failures == []


@pytest.mark.parametrize(
    ("hdu", "keyword", "value"),
    [
        ("primary", "EXPTIME", 0.0),
        ("primary", "EXPTIME", -5.0),
        ("primary", "TEXPTIME", 0.0),
        ("extension", "TELAPSE", 0.0),
        ("primary", "NCOMBINE", 0),
        ("primary", "SNR", -1.0),
        ("primary", "WAVELMIN", 0.0),
        ("primary", "WAVELMIN", 650.0),
        ("primary", "SPEC_RES", 0.0),
    ],
)
def test_a_value_outside_its_allowed_range_is_a_failure(
    tmp_path: Path, log: object, hdu: str, keyword: str, value: float
) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, **{hdu: {keyword: value}})

    # ASSERT
    assert any(keyword in message and "range" in message for message in failures)


@pytest.mark.parametrize(
    ("hdu", "keyword", "value", "source"),
    [
        ("primary", "EXPTIME", 0.0, "(standard section 5)"),
        ("primary", "TEXPTIME", 0.0, "(house rule)"),
        ("extension", "TELAPSE", 0.0, "(house rule)"),
        ("primary", "NCOMBINE", 0, "(house rule)"),
        ("primary", "SNR", -1.0, "(house rule)"),
        ("primary", "WAVELMIN", 0.0, "(house rule)"),
        ("primary", "WAVELMIN", 650.0, "(house rule)"),
        ("primary", "SPEC_RES", 0.0, "(house rule)"),
    ],
)
def test_the_range_failure_names_where_the_rule_comes_from(
    tmp_path: Path, log: object, hdu: str, keyword: str, value: float, source: str
) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, **{hdu: {keyword: value}})

    # ASSERT
    assert any(keyword in message and f"allowed range {source}" in message for message in failures)


@pytest.mark.parametrize(("keyword", "value"), [("SNR", 0.0), ("NCOMBINE", 1)])
def test_a_value_at_the_edge_of_its_allowed_range_is_accepted(
    tmp_path: Path, log: object, keyword: str, value: float
) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, primary={keyword: value})

    # ASSERT
    assert all("range" not in message for message in failures)


def test_a_missing_extname_is_labelled_as_a_house_rule(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)
    rewrite_product(product, lambda hdus: hdus[1].header.remove("EXTNAME"))

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "house rule" in message and "DY-2393" in message for message in failures)


def test_a_wrong_extname_is_labelled_as_a_house_rule(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), extension={"EXTNAME": "TABLE"})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("EXTNAME" in message and "house rule" in message for message in failures)


@pytest.mark.parametrize("keyword", ["HIERARCH ESO DPR TYPE", "HIERARCH ESO DPR TECH", "HIERARCH ESO DPR CATG"])
def test_a_dpr_card_in_a_phase3_product_is_a_failure(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={keyword: "OBJECT"})

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any(keyword.replace("HIERARCH ", "") in message and "writer removes" in message for message in failures)


def test_the_equinox_card_that_the_writer_renames_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(
        phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary={"HIERARCH ESO TEL TARG EQUINOX": 2000.0}
    )

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("ESO TEL TARG EQUINOX" in message and "writer" in message for message in failures)


def test_other_eso_cards_are_kept_in_a_per_arm_product(tmp_path: Path, log: object) -> None:
    # ARRANGE
    cards = {
        "HIERARCH ESO OBS ID": 42,
        "HIERARCH ESO QC RON": 3.2,
        "HIERARCH ESO DET NDIT": 2,
        "HIERARCH ESO TEL TARG ALPHA": 123456.789,
    }
    product = add_cards(phase3_spectrum_file(tmp_path / PRODUCT_NAME, log), primary=cards)

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert findings == []


def test_a_dpr_card_in_an_internal_spectrum_is_not_checked(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = add_cards(
        phase3_spectrum_file(tmp_path / "soxs_stare_SPECTRUM.fits", log), primary={"HIERARCH ESO DPR TYPE": "OBJECT"}
    )

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert all("DPR" not in message for message in failures)


def test_the_writer_removes_and_renames_the_cards_the_checker_looks_for(tmp_path: Path, log: object) -> None:
    # ARRANGE
    cards = {f"HIERARCH {keyword}": "x" for keyword in SCRUBBED_ESO_CARDS}
    cards.update({f"HIERARCH {keyword}": 2000.0 for keyword in RENAMED_ESO_CARDS})

    # ACT
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary=cards)

    # ASSERT
    header = fits.getheader(product, 0)
    assert not [keyword for keyword in (*SCRUBBED_ESO_CARDS, *RENAMED_ESO_CARDS) if keyword in header]
    assert header["EQUINOX"] == 2000.0
