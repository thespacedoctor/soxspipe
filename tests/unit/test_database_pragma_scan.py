"""Source scan: production code assigns `synchronous` and `journal_mode` only literal durable values (DY-1706).

`PRAGMA synchronous = OFF` (or NORMAL) lets an operating-system crash or a power loss corrupt `soxspipe.db`, and WAL is
unsafe on the network volumes a workspace can sit on. The scan reads every module under `soxspipe/` and fails on a
`synchronous` assignment that is not a literal FULL, EXTRA, 2 or 3, which includes a value only known at run time such
as `f"PRAGMA synchronous = {level}"`. It fails on a `journal_mode` assignment that is not a literal DELETE. A query
(`PRAGMA synchronous;`) is not an assignment and passes.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_ROOT = Path(soxspipe.__file__).resolve().parent
# POSSESSIVE QUANTIFIERS (`*+`) STOP THE MATCH BACKTRACKING INTO THE SPACE BEFORE THE VALUE, WHICH WOULD HIDE IT FROM
# THE LOOKAHEAD
NON_DURABLE_SYNCHRONOUS = re.compile(r"synchronous\s*+(?:=|\()\s*+['\"]?+\s*+(?!(?:FULL|EXTRA|2|3)\b)", re.IGNORECASE)
NON_ROLLBACK_JOURNAL_MODE = re.compile(r"journal_mode\s*+(?:=|\()\s*+['\"]?+\s*+(?!DELETE\b)", re.IGNORECASE)


def _production_python_files() -> list[Path]:
    return sorted(path for path in PACKAGE_ROOT.rglob("*.py") if "tests" not in path.relative_to(PACKAGE_ROOT).parts)


def _matches(pattern: re.Pattern[str]) -> list[str]:
    found = []
    for path in _production_python_files():
        for lineNumber, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if pattern.search(line):
                found.append(f"{path.relative_to(PACKAGE_ROOT)}:{lineNumber}: {line.strip()}")
    return found


def test_production_code_scan_covers_the_database_modules() -> None:
    """The scan reads the modules that open the database, so an empty result cannot come from an empty scan."""
    scanned = {path.name for path in _production_python_files()}

    assert {"base_recipe.py", "reducer.py", "data_organiser.py", "workspace_database.py"} <= scanned


@pytest.mark.parametrize(
    "line",
    [
        'c.execute("PRAGMA synchronous = OFF")',
        'c.execute("PRAGMA synchronous=0")',
        'c.execute("PRAGMA synchronous = NORMAL")',
        'c.execute("PRAGMA synchronous(1)")',
        'c.execute(f"PRAGMA synchronous = {level}")',
        "conn.execute('pragma synchronous = off')",
    ],
)
def test_synchronous_pattern_flags_a_value_that_is_not_a_literal_durable_one(line: str) -> None:
    """OFF, NORMAL, 0, 1 and a value only known at run time are all flagged."""
    assert NON_DURABLE_SYNCHRONOUS.search(line)


@pytest.mark.parametrize(
    "line",
    [
        'c.execute("PRAGMA synchronous = FULL")',
        'c.execute("PRAGMA synchronous = EXTRA")',
        'c.execute("PRAGMA synchronous = 2")',
        'c.execute("PRAGMA synchronous = 3")',
        'c.execute("PRAGMA synchronous;")',
        'c.execute("PRAGMA synchronous").fetchone()',
    ],
)
def test_synchronous_pattern_passes_full_extra_and_a_query(line: str) -> None:
    """FULL, EXTRA, 2 and 3 are durable, and a bare query assigns nothing."""
    assert NON_DURABLE_SYNCHRONOUS.search(line) is None


@pytest.mark.parametrize(
    "line",
    [
        'c.execute("PRAGMA journal_mode = WAL")',
        'c.execute("PRAGMA journal_mode=WAL;")',
        'c.execute(f"PRAGMA journal_mode = {mode}")',
    ],
)
def test_journal_mode_pattern_flags_a_value_that_is_not_delete(line: str) -> None:
    """WAL and a value only known at run time are flagged."""
    assert NON_ROLLBACK_JOURNAL_MODE.search(line)


def test_journal_mode_pattern_passes_delete() -> None:
    """The rollback journal is the one allowed assignment."""
    assert NON_ROLLBACK_JOURNAL_MODE.search('target.execute("PRAGMA journal_mode=DELETE;")') is None


def test_no_production_code_assigns_a_non_durable_synchronous_value() -> None:
    """Nothing under `soxspipe/` assigns `synchronous` anything but a literal FULL, EXTRA, 2 or 3."""
    assert _matches(NON_DURABLE_SYNCHRONOUS) == []


def test_no_production_code_assigns_a_journal_mode_other_than_delete() -> None:
    """Nothing under `soxspipe/` switches the journal to WAL or any other mode but the rollback journal."""
    assert _matches(NON_ROLLBACK_JOURNAL_MODE) == []
