"""*tests for the changelog duplicate-bullet checker*

The checker is pure text processing, so these tests feed it changelog text directly and only
touch the filesystem for the CLI path.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "check_changelog.py"


def _load_check_changelog() -> ModuleType:
    """*import the checker from tools/, which is not an installed package*

    **Return:**

    - ``module`` -- the imported `check_changelog` module
    """
    specification = importlib.util.spec_from_file_location("check_changelog", MODULE_PATH)
    if specification is None or specification.loader is None:
        raise AssertionError(f"the changelog checker is not importable from {MODULE_PATH}")

    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    return module


check_changelog = _load_check_changelog()


def _write_changelog(tmp_path: Path, text: str) -> Path:
    """*write changelog text to a file under the test's temporary directory*

    **Key Arguments:**

    - ``tmp_path`` -- the test's temporary directory
    - ``text`` -- the changelog text

    **Return:**

    - ``path`` -- the path of the file written
    """
    path = tmp_path / "CHANGES.md"
    path.write_text(text, encoding="utf-8")

    return path


CLEAN_CHANGELOG = """# Release Notes

* **FIXED**: first fix
* **FIXED**: second fix

## v1.0.0 - January 1, 2026

* **FEATURE**: a feature
* **DOCS**: some docs
"""

UNRELEASED_DUPLICATE_CHANGELOG = """# Release Notes

* **FIXED**: first fix
* **TEST**: repeated entry
* **TEST**: repeated entry

## v1.0.0 - January 1, 2026

* **FEATURE**: a feature
"""

VERSIONED_DUPLICATE_CHANGELOG = """# Release Notes

* **FIXED**: first fix

## v1.0.0 - January 1, 2026

* **FEATURE**: a feature
* **DOCS**: some docs
* **FEATURE**: a feature
"""

CROSS_BLOCK_CHANGELOG = """# Release Notes

* **FIXED**: shared entry

## v1.0.0 - January 1, 2026

* **FIXED**: shared entry

## v0.9.0 - December 1, 2025

* **FIXED**: shared entry
"""


def test_returns_no_duplicates_for_a_clean_changelog():
    duplicates = check_changelog.find_duplicates(CLEAN_CHANGELOG)

    assert duplicates == []


def test_reports_a_duplicate_in_the_unreleased_block_with_both_line_numbers():
    duplicates = check_changelog.find_duplicates(UNRELEASED_DUPLICATE_CHANGELOG)

    assert len(duplicates) == 1
    assert duplicates[0].text == "* **TEST**: repeated entry"
    assert (duplicates[0].firstLine, duplicates[0].secondLine) == (4, 5)
    assert duplicates[0].heading == check_changelog.UNRELEASED_HEADING


def test_reports_a_duplicate_in_a_versioned_block_with_its_heading():
    duplicates = check_changelog.find_duplicates(VERSIONED_DUPLICATE_CHANGELOG)

    assert len(duplicates) == 1
    assert duplicates[0].heading == "## v1.0.0 - January 1, 2026"
    assert (duplicates[0].firstLine, duplicates[0].secondLine) == (7, 9)


def test_allows_the_same_bullet_in_different_release_blocks():
    duplicates = check_changelog.find_duplicates(CROSS_BLOCK_CHANGELOG)

    assert duplicates == []


def test_ignores_trailing_whitespace_when_comparing_bullets():
    text = "* **DOCS:** brought up to date  \n* **DOCS:** brought up to date\n"

    duplicates = check_changelog.find_duplicates(text)

    assert [(d.firstLine, d.secondLine) for d in duplicates] == [(1, 2)]


def test_ignores_repeated_lines_that_are_not_bullets():
    text = "## v1.0.0\n\nsame prose\nsame prose\n\n    * indented bullet\n    * indented bullet\n"

    duplicates = check_changelog.find_duplicates(text)

    assert duplicates == []


def test_reports_each_repeat_of_a_bullet_against_its_first_occurrence():
    text = "## v1.0.0\n\n* same\n* same\n* same\n"

    duplicates = check_changelog.find_duplicates(text)

    assert [(d.firstLine, d.secondLine) for d in duplicates] == [(3, 4), (3, 5)]


def test_formats_a_duplicate_with_heading_and_both_line_numbers():
    duplicates = check_changelog.find_duplicates(VERSIONED_DUPLICATE_CHANGELOG)

    message = str(duplicates[0])

    assert "## v1.0.0 - January 1, 2026" in message
    assert "7" in message and "9" in message
    assert "* **FEATURE**: a feature" in message


def test_main_returns_zero_for_a_clean_file(tmp_path, capsys):
    path = _write_changelog(tmp_path, CLEAN_CHANGELOG)

    status = check_changelog.main([str(path)])

    assert status == 0
    assert capsys.readouterr().err == ""


def test_main_returns_one_and_prints_each_duplicate(tmp_path, capsys):
    path = _write_changelog(tmp_path, UNRELEASED_DUPLICATE_CHANGELOG)

    status = check_changelog.main([str(path)])

    output = capsys.readouterr().out
    assert status == 1
    assert "repeated entry" in output
    assert "4" in output and "5" in output


def test_main_returns_two_when_the_file_is_unreadable(tmp_path, capsys):
    missing = tmp_path / "missing.md"

    status = check_changelog.main([str(missing)])

    assert status == 2
    assert "missing.md" in capsys.readouterr().err


def test_main_returns_two_when_the_file_is_not_utf8(tmp_path):
    path = tmp_path / "CHANGES.md"
    path.write_bytes(b"* \xff\xfe broken\n")

    status = check_changelog.main([str(path)])

    assert status == 2


def test_main_defaults_to_changes_md_in_the_working_directory(tmp_path, monkeypatch):
    _write_changelog(tmp_path, UNRELEASED_DUPLICATE_CHANGELOG)
    monkeypatch.chdir(tmp_path)

    status = check_changelog.main([])

    assert status == 1


def test_main_reads_a_release_heading_on_the_first_line_of_a_bom_prefixed_file(tmp_path, capsys):
    path = tmp_path / "CHANGES.md"
    path.write_bytes("﻿## v1.0.0\n* **FIXED**: a\n* **FIXED**: a\n".encode())

    status = check_changelog.main([str(path)])

    assert status == 1
    assert "## v1.0.0: lines 2 and 3" in capsys.readouterr().out
