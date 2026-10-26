"""The real-data gate's chain-wide QC-verdict check, run against hand-built databases (DY-908)."""

from __future__ import annotations

import sqlite3
from collections.abc import Iterator

import pytest

from tests.real_data.qc_verdict import AllowedQcFailure, qc_verdict_problems

pytestmark = pytest.mark.unit


@pytest.fixture
def connection() -> Iterator[sqlite3.Connection]:
    """An in-memory database holding an empty `quality_control` table."""
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "CREATE TABLE quality_control "
        "(soxspipe_recipe, sof_name, qc_name, qc_value, qc_value_min, qc_value_max, qc_flag)"
    )
    yield connection
    connection.close()


def _insert(connection: sqlite3.Connection, *rows: tuple) -> None:
    connection.executemany(
        "INSERT INTO quality_control "
        "(soxspipe_recipe, sof_name, qc_name, qc_value, qc_value_min, qc_value_max, qc_flag) "
        "VALUES (?, ?, ?, ?, ?, ?, ?)",
        rows,
    )


_PASSING_ROW = ("soxs-mbias", "VIS_MBIAS.sof", "MBIAS MEDIAN", 644.5, None, None, "pass")
_FAILING_ROW = ("soxs-mflat", "VIS_MFLAT.sof", "COLDPIX FRAC", 0.31, 0.0, 0.1, "fail")


def test_reports_every_detail_of_a_failing_qc_row(connection: sqlite3.Connection) -> None:
    # ARRANGE
    _insert(connection, _PASSING_ROW, _FAILING_ROW)

    # ACT
    problems = qc_verdict_problems(connection, allowedFailures=())

    # ASSERT
    assert len(problems) == 1
    for detail in ("soxs-mflat", "VIS_MFLAT.sof", "COLDPIX FRAC", "0.31", "0.0", "0.1"):
        assert detail in problems[0]


def test_passes_when_the_failing_row_is_allowed(connection: sqlite3.Connection) -> None:
    # ARRANGE
    _insert(connection, _PASSING_ROW, _FAILING_ROW)
    allowed = (AllowedQcFailure("soxs-mflat", "VIS_MFLAT.sof", "COLDPIX FRAC", reason="known vignetting"),)

    # ACT
    problems = qc_verdict_problems(connection, allowedFailures=allowed)

    # ASSERT
    assert problems == []


def test_an_allowance_for_another_sof_does_not_hide_the_failure(connection: sqlite3.Connection) -> None:
    # ARRANGE
    _insert(connection, _FAILING_ROW)
    allowed = (AllowedQcFailure("soxs-mflat", "NIR_MFLAT.sof", "COLDPIX FRAC", reason="another arm"),)

    # ACT
    problems = qc_verdict_problems(connection, allowedFailures=allowed)

    # ASSERT
    assert len(problems) == 1


def test_an_empty_table_is_a_problem_not_a_pass(connection: sqlite3.Connection) -> None:
    # ACT
    problems = qc_verdict_problems(connection, allowedFailures=())

    # ASSERT
    assert problems == ["the quality_control table holds no rows, so no QC verdict was recorded"]


def test_a_missing_table_is_a_problem_not_a_pass() -> None:
    # ARRANGE
    connection = sqlite3.connect(":memory:")

    # ACT
    problems = qc_verdict_problems(connection, allowedFailures=())
    connection.close()

    # ASSERT
    assert problems == ["the workspace database has no quality_control table"]
