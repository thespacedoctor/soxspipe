"""No production module weakens SQLite's `synchronous` setting (DY-1706).

`PRAGMA synchronous = OFF` makes an operating-system crash or a power loss able to corrupt `soxspipe.db`. The
database keeps SQLite's default, `FULL`, so nothing under `soxspipe/` may set it to a weaker value.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_ROOT = Path(soxspipe.__file__).resolve().parent
WEAKENED_SYNCHRONOUS = re.compile(
    r"synchronous\s*(?:=|\()\s*['\"]?\s*(?!FULL\b|2\b)\w+",
    re.IGNORECASE,
)
JOURNAL_MODE_OTHER_THAN_DELETE = re.compile(
    r"journal_mode\s*(?:=|\()\s*['\"]?\s*(?!DELETE\b)\w+",
    re.IGNORECASE,
)


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

    assert {"base_recipe.py", "reducer.py", "data_organiser.py"} <= scanned


def test_no_production_code_weakens_the_synchronous_pragma() -> None:
    """Nothing under `soxspipe/` sets `synchronous` to a value other than `FULL`."""
    offenders = _matches(WEAKENED_SYNCHRONOUS)

    assert offenders == []


def test_no_production_code_changes_the_journal_mode_away_from_delete() -> None:
    """Nothing under `soxspipe/` switches the journal to WAL or any other mode but the rollback journal."""
    offenders = _matches(JOURNAL_MODE_OTHER_THAN_DELETE)

    assert offenders == []
