"""Unit tests for the frame-accounting checks of the Phase 3 checker (DY-2356).

Every good product is written by the real Phase 3 table writer, then one fault is put in. A per-arm product carries the
``ESO SEQ ARM`` card, which the joined writer drops, so the rule that ``EXPTIME`` equals ``TEXPTIME`` reads that card
(Dave, 2026-11-24: the joined product has the median ``EXPTIME`` and the maximum ``TEXPTIME``).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.factories.phase3_products import failure_messages

pytestmark = pytest.mark.unit

PER_ARM_CARD = {"HIERARCH ESO SEQ ARM": "VIS"}


def _failures_mentioning(tmp_path: Path, log: object, keyword: str, **factoryArguments: object) -> list[str]:
    return [message for message in failure_messages(tmp_path, log, **factoryArguments) if keyword in message]


def test_a_good_product_has_no_frame_accounting_failure(tmp_path: Path, log: object) -> None:
    # ACT / ASSERT
    assert failure_messages(tmp_path, log) == []


def test_a_good_per_arm_product_has_no_failure(tmp_path: Path, log: object) -> None:
    # ACT / ASSERT
    assert failure_messages(tmp_path, log, primary=PER_ARM_CARD) == []


def test_ncombine_must_equal_the_number_of_provenance_cards(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "NCOMBINE", primary={"NCOMBINE": 3})

    # ASSERT
    assert any("PROV" in message and "3" in message and "2" in message for message in messages)


def test_the_provenance_cards_must_run_from_one_without_a_gap(tmp_path: Path, log: object) -> None:
    # ARRANGE: PROV2 IS MISSING AND PROV3 IS THERE, SO THERE ARE TWO CARDS FOR NCOMBINE 2
    primary = {"PROV2": None, "PROV3": "SOXS.2026-01-02T03:04:07.000"}

    # ACT
    messages = _failures_mentioning(tmp_path, log, "PROV", primary=primary)

    # ASSERT
    assert any("PROV2" in message for message in messages)


def test_the_observation_block_cards_must_run_from_one_without_a_gap(tmp_path: Path, log: object) -> None:
    # ARRANGE: OBID1 IS THERE, OBID2 IS MISSING AND OBID3 IS THERE
    primary = {"OBID3": 44}

    # ACT
    messages = _failures_mentioning(tmp_path, log, "OBID", primary=primary)

    # ASSERT
    assert any("OBID2" in message for message in messages)


def test_the_end_of_the_observation_must_be_after_its_start(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "MJD-END", primary={"MJD-END": 61000.1})

    # ASSERT
    assert any("MJD-OBS" in message for message in messages)


def test_a_provenance_card_must_not_end_in_fits(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "PROV1", primary={"PROV1": "SOXS.2026-01-02T03:04:05.678.fits"})

    # ASSERT
    assert any(".fits" in message for message in messages)


def test_a_per_arm_exptime_must_equal_its_texptime(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "TEXPTIME", primary={**PER_ARM_CARD, "TEXPTIME": 590.0})

    # ASSERT
    assert any("EXPTIME" in message and "590" in message for message in messages)


def test_a_joined_product_may_have_an_exptime_that_differs_from_its_texptime(tmp_path: Path, log: object) -> None:
    # ARRANGE: THE JOINED PRODUCT HAS THE MEDIAN EXPTIME AND THE MAXIMUM TEXPTIME, AND NO ESO SEQ CARD
    # ACT / ASSERT
    assert _failures_mentioning(tmp_path, log, "TEXPTIME", primary={"TEXPTIME": 650.0}) == []


def test_a_per_arm_exptime_that_is_not_a_number_is_left_to_the_type_check(tmp_path: Path, log: object) -> None:
    # ACT
    messages = _failures_mentioning(tmp_path, log, "is not TEXPTIME", primary={**PER_ARM_CARD, "EXPTIME": "long"})

    # ASSERT
    assert messages == []
