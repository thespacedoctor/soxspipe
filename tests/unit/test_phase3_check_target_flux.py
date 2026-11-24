"""Unit tests for the target position and flux-calibration checks of the Phase 3 checker (DY-2357).

Every good product is written by the real Phase 3 table writer, then one fault is put in. The rules that the standard
already gives these keywords (mandatory, type, ``RADESYS``, ``EQUINOX``, ``FLUXERR`` of an uncalibrated spectrum) are
tested with the other header checks. These tests cover the house rules of DY-2393 only.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.factories.phase3_products import failure_messages

pytestmark = pytest.mark.unit


def _failures_mentioning(tmp_path: Path, log: object, keyword: str, **factoryArguments: object) -> list[str]:
    return [message for message in failure_messages(tmp_path, log, **factoryArguments) if keyword in message]


@pytest.mark.parametrize(
    "primary",
    [
        {"FLUXCAL": "ABSOLUTE", "FLUXERR": -2.0},
        {"FLUXCAL": "UNCALIBRATED", "FLUXERR": -1.0},
        {"RA": 0.0, "DEC": -90.0},
        {"RA": 359.9999, "DEC": 90.0},
    ],
)
def test_a_good_target_and_flux_product_has_no_failure(tmp_path: Path, log: object, primary: dict) -> None:
    # ARRANGE: THE EXTENSION REPEATS THE POSITION OF THE PRIMARY HEADER
    extension = {keyword: value for keyword, value in primary.items() if keyword in ("RA", "DEC")}

    # ACT / ASSERT
    assert failure_messages(tmp_path, log, primary=primary, extension=extension) == []


def test_an_absolute_spectrum_must_have_a_flux_error_of_minus_two(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "FLUXERR", primary={"FLUXCAL": "ABSOLUTE", "FLUXERR": -1.0})

    # ASSERT
    assert any("-2" in message and "ABSOLUTE" in message and "house rule, DY-2393" in message for message in messages)


def test_an_absolute_spectrum_with_a_positive_flux_error_fails(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "FLUXERR", primary={"FLUXCAL": "ABSOLUTE", "FLUXERR": 0.5})

    # ASSERT
    assert any("-2" in message for message in messages)


def test_tot_flux_must_be_false(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "TOT_FLUX", primary={"TOT_FLUX": True})

    # ASSERT
    assert any("False" in message and "house rule, DY-2393" in message for message in messages)


@pytest.mark.parametrize("ra", [360.0, 400.0, -0.5])
def test_the_right_ascension_must_be_in_zero_to_less_than_360(tmp_path: Path, log: object, ra: float) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "RA is", primary={"RA": ra}, extension={"RA": ra})

    # ASSERT
    assert messages and all("house rule" in message for message in messages)


@pytest.mark.parametrize("dec", [90.5, -91.0])
def test_the_declination_must_be_in_minus_90_to_90(tmp_path: Path, log: object, dec: float) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "DEC is", primary={"DEC": dec}, extension={"DEC": dec})

    # ASSERT
    assert messages and all("house rule" in message for message in messages)


def test_a_position_that_is_not_a_number_is_left_to_the_type_check(tmp_path: Path, log: object) -> None:
    # ACT
    messages = failure_messages(tmp_path, log, primary={"RA": "07:35", "DEC": "text"})

    # ASSERT
    assert not [message for message in messages if "RA is" in message or "DEC is" in message]
    assert any("RA" in message for message in messages)
