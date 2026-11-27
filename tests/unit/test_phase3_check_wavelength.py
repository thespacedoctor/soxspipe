"""Unit tests for the wavelength-calibration and resolution checks of the Phase 3 checker (DY-2358).

Every good product is written by the real Phase 3 table writer, then one fault is put in. The presence and type of
`SPEC_RES`, `SPEC_ERR`, `SPEC_SYE`, `LAMRMS` and `LAMNLIN` are rules of the standard, which the header checks hold.
These tests cover that the severities follow the standard for these keywords, and the house rules of DY-2358 on their
values.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from soxspipe.commonutils.phase3_check import check_phase3
from tests.factories.phase3_products import PRODUCT_NAME, failure_messages, finding_messages, phase3_spectrum_file

pytestmark = pytest.mark.unit

GOOD_RMS = 0.004
GOOD_LINES = 120
GOOD_ERR = GOOD_RMS / math.sqrt(GOOD_LINES)


def _mentioning(messages: list[str], keyword: str) -> list[str]:
    return [message for message in messages if keyword in message]


def test_a_product_with_the_resolution_and_accuracy_cards_has_no_finding(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"LAMRMS": GOOD_RMS, "LAMNLIN": GOOD_LINES})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert finding_messages(findings, "failure") == []


def test_a_missing_spec_res_fails_the_product(tmp_path: Path, log: object) -> None:
    # ACT
    messages = failure_messages(tmp_path, log, primary={"SPEC_RES": None})

    # ASSERT
    assert any("SPEC_RES" in message and "missing" in message for message in messages)


@pytest.mark.parametrize("keyword", ["LAMRMS", "LAMNLIN", "SPEC_ERR"])
def test_a_missing_accuracy_card_is_a_warning_not_a_failure(tmp_path: Path, log: object, keyword: str) -> None:
    # ARRANGE: AN OLD DISP_TAB LEAVES ALL THREE OUT, SO THE PRODUCT HAS NONE OF THEM
    product = phase3_spectrum_file(
        tmp_path / PRODUCT_NAME, log, primary={"LAMRMS": None, "LAMNLIN": None, "SPEC_ERR": None}
    )

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert _mentioning(finding_messages(findings, "failure"), keyword) == []
    assert any(keyword in message and "missing" in message for message in finding_messages(findings, "warning"))


@pytest.mark.parametrize(
    ("primary", "keyword"),
    [
        ({"LAMRMS": 0.0}, "LAMRMS"),
        ({"LAMRMS": -0.004}, "LAMRMS"),
        ({"LAMNLIN": 0}, "LAMNLIN"),
        ({"SPEC_ERR": 0.0}, "SPEC_ERR"),
        ({"SPEC_ERR": -1e-4}, "SPEC_ERR"),
    ],
)
def test_a_value_outside_its_range_fails(tmp_path: Path, log: object, primary: dict, keyword: str) -> None:
    # ARRANGE: KEEP SPEC_ERR CONSISTENT WITH THE OTHER TWO WHERE THE FAULT IS NOT IN SPEC_ERR
    cards = {"LAMRMS": GOOD_RMS, "LAMNLIN": GOOD_LINES, "SPEC_ERR": GOOD_ERR, **primary}

    # ACT
    messages = _mentioning(failure_messages(tmp_path, log, primary=cards), keyword)

    # ASSERT
    assert any("outside its allowed range" in message for message in messages)


def test_a_spec_err_that_is_not_lamrms_over_the_root_of_lamnlin_fails(tmp_path: Path, log: object) -> None:
    # ARRANGE: ESO'S OWN X-SHOOTER PRODUCTS USE LAMRMS / LAMNLIN, WHICH THE STANDARD DOES NOT
    wrong = GOOD_RMS / GOOD_LINES

    # ACT
    messages = _mentioning(
        failure_messages(tmp_path, log, primary={"LAMRMS": GOOD_RMS, "LAMNLIN": GOOD_LINES, "SPEC_ERR": wrong}),
        "SPEC_ERR",
    )

    # ASSERT
    assert any("sqrt(LAMNLIN)" in message and "house rule, DY-2358" in message for message in messages)


def test_a_spec_err_that_follows_from_lamrms_and_lamnlin_passes(tmp_path: Path, log: object) -> None:
    # ACT
    messages = failure_messages(
        tmp_path, log, primary={"LAMRMS": GOOD_RMS, "LAMNLIN": GOOD_LINES, "SPEC_ERR": GOOD_ERR}
    )

    # ASSERT
    assert _mentioning(messages, "SPEC_ERR") == []


def test_the_spec_err_rule_is_not_applied_when_lamrms_is_missing(tmp_path: Path, log: object) -> None:
    # ACT
    messages = failure_messages(tmp_path, log, primary={"LAMRMS": None, "LAMNLIN": GOOD_LINES, "SPEC_ERR": 0.5})

    # ASSERT
    assert _mentioning(messages, "SPEC_ERR") == []


def test_a_non_integer_line_count_fails(tmp_path: Path, log: object) -> None:
    # ACT
    messages = failure_messages(tmp_path, log, primary={"LAMNLIN": 120.5, "SPEC_ERR": None})

    # ASSERT
    assert any("LAMNLIN" in message and "must be int" in message for message in messages)


def test_a_negative_spec_sye_fails(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _mentioning(failure_messages(tmp_path, log, primary={"SPEC_SYE": -0.001}), "SPEC_SYE")

    # ASSERT
    assert any("outside its allowed range" in message and "house rule, DY-2359" in message for message in messages)


@pytest.mark.parametrize("value", [0.0, 0.0285])
def test_a_spec_sye_of_zero_or_more_passes(tmp_path: Path, log: object, value: float) -> None:
    # ARRANGE: THE STANDARD ALLOWS 0 WHEN THE SYSTEMATIC ERROR IS KNOWN TO BE INSIGNIFICANT

    # ACT
    messages = failure_messages(tmp_path, log, primary={"SPEC_SYE": value})

    # ASSERT
    assert _mentioning(messages, "SPEC_SYE") == []


def test_a_missing_spec_sye_is_a_warning_and_the_range_rule_is_not_applied(tmp_path: Path, log: object) -> None:
    # ARRANGE: A NIR PRODUCT HAS NO SPEC_SYE UNTIL ITS ERROR FORMULAS EXIST
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log, primary={"SPEC_SYE": None})

    # ACT
    findings = check_phase3(product)

    # ASSERT
    assert _mentioning(finding_messages(findings, "failure"), "SPEC_SYE") == []
    assert any("SPEC_SYE" in message and "missing" in message for message in finding_messages(findings, "warning"))
