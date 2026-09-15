"""*tests for the ruff lint ratchet's diff and finding intersection*

The intersection logic is pure, so these tests feed it diff text and ruff JSON directly
rather than shelling out to git or ruff.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import ModuleType
from typing import NamedTuple

import pytest

pytestmark = pytest.mark.unit

REPO_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO_ROOT / "tools" / "lint_ratchet.py"


def _load_lint_ratchet() -> ModuleType:
    """*import the ratchet from tools/, which is not an installed package*

    **Return:**

    - ``module`` -- the imported `lint_ratchet` module
    """
    specification = importlib.util.spec_from_file_location("lint_ratchet", MODULE_PATH)
    if specification is None or specification.loader is None:
        raise AssertionError(f"the lint ratchet is not importable from {MODULE_PATH}")

    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)

    return module


lint_ratchet = _load_lint_ratchet()


class _CompletedCommand(NamedTuple):
    """*the parts of `subprocess.CompletedProcess` the ratchet reads*"""

    returncode: int
    stdout: str
    stderr: str


def _ruff_report(*findings: tuple[str, int, str, str]) -> str:
    """*a ruff JSON report holding the given findings*

    **Key Arguments:**

    - ``findings`` -- (relative path, row, code, message) tuples

    **Return:**

    - ``ruffJson`` -- the report as ruff would print it
    """
    return json.dumps(
        [
            {
                "filename": str(REPO_ROOT / path),
                "location": {"row": row, "column": 1},
                "end_location": {"row": row, "column": 2},
                "code": code,
                "message": message,
            }
            for path, row, code, message in findings
        ]
    )


def test_changed_line_with_a_finding_fails_the_ratchet():
    # ARRANGE
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -10,0 +11,2 @@\n+x = 1\n+y = 2\n"
    changedLines = lint_ratchet.parse_changed_lines(diffText)
    findings = lint_ratchet.parse_ruff_findings(
        _ruff_report(("soxspipe/recipes/soxs_mbias.py", 12, "E501", "Line too long (256 > 120)")), REPO_ROOT
    )

    # ACT
    gated = lint_ratchet.findings_on_changed_lines(findings, changedLines)

    # ASSERT
    assert [(finding.path, finding.line, finding.code) for finding in gated] == [
        ("soxspipe/recipes/soxs_mbias.py", 12, "E501")
    ]


def test_unchanged_line_with_a_finding_passes_the_ratchet():
    # ARRANGE
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -10,0 +11,1 @@\n+x = 1\n"
    changedLines = lint_ratchet.parse_changed_lines(diffText)
    findings = lint_ratchet.parse_ruff_findings(
        _ruff_report(("soxspipe/recipes/soxs_mbias.py", 400, "E722", "Do not use bare `except`")), REPO_ROOT
    )

    # ACT
    gated = lint_ratchet.findings_on_changed_lines(findings, changedLines)

    # ASSERT
    assert gated == []


def test_changed_line_with_no_finding_passes_the_ratchet():
    # ARRANGE
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -10,0 +11,1 @@\n+x = 1\n"
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ACT
    gated = lint_ratchet.findings_on_changed_lines([], changedLines)

    # ASSERT
    assert gated == []


def test_a_renamed_file_is_gated_under_its_new_path():
    # ARRANGE
    diffText = (
        "diff --git a/soxspipe/commonutils/old_name.py b/soxspipe/commonutils/new_name.py\n"
        "similarity index 96%\n"
        "rename from soxspipe/commonutils/old_name.py\n"
        "rename to soxspipe/commonutils/new_name.py\n"
        "--- a/soxspipe/commonutils/old_name.py\n"
        "+++ b/soxspipe/commonutils/new_name.py\n"
        "@@ -5,0 +6,1 @@\n"
        "+import os, sys\n"
    )
    changedLines = lint_ratchet.parse_changed_lines(diffText)
    findings = lint_ratchet.parse_ruff_findings(
        _ruff_report(
            ("soxspipe/commonutils/new_name.py", 6, "E401", "Multiple imports on one line"),
            ("soxspipe/commonutils/old_name.py", 6, "E401", "Multiple imports on one line"),
        ),
        REPO_ROOT,
    )

    # ACT
    gated = lint_ratchet.findings_on_changed_lines(findings, changedLines)

    # ASSERT
    assert list(changedLines) == ["soxspipe/commonutils/new_name.py"]
    assert [finding.path for finding in gated] == ["soxspipe/commonutils/new_name.py"]


def test_a_pure_rename_with_no_edits_changes_no_lines():
    # ARRANGE
    diffText = (
        "diff --git a/soxspipe/commonutils/old_name.py b/soxspipe/commonutils/new_name.py\n"
        "similarity index 100%\n"
        "rename from soxspipe/commonutils/old_name.py\n"
        "rename to soxspipe/commonutils/new_name.py\n"
    )

    # ACT
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ASSERT
    assert changedLines == {}


def test_a_deleted_file_contributes_no_changed_lines():
    # ARRANGE
    diffText = "--- a/soxspipe/commonutils/gone.py\n+++ /dev/null\n@@ -1,5 +0,0 @@\n-x = 1\n"

    # ACT
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ASSERT
    assert changedLines == {}


def test_a_pure_deletion_hunk_marks_no_lines_changed():
    # ARRANGE
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -20,3 +19,0 @@\n-x = 1\n-y = 2\n-z = 3\n"

    # ACT
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ASSERT
    assert changedLines == {"soxspipe/recipes/soxs_mbias.py": set()}


def test_a_hunk_header_without_a_length_covers_one_line():
    # ARRANGE
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -10 +11 @@\n-x = 1\n+x = 2\n"

    # ACT
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ASSERT
    assert changedLines == {"soxspipe/recipes/soxs_mbias.py": {11}}


def test_several_hunks_in_one_file_accumulate():
    # ARRANGE
    diffText = (
        "+++ b/soxspipe/recipes/soxs_mbias.py\n"
        "@@ -10,0 +11,2 @@\n"
        "+x = 1\n"
        "+y = 2\n"
        "@@ -40,0 +43,1 @@\n"
        "+z = 3\n"
    )

    # ACT
    changedLines = lint_ratchet.parse_changed_lines(diffText)

    # ASSERT
    assert changedLines == {"soxspipe/recipes/soxs_mbias.py": {11, 12, 43}}


def test_findings_are_reported_with_enough_detail_to_fix_without_rerunning():
    # ARRANGE
    findings = lint_ratchet.parse_ruff_findings(
        _ruff_report(("soxspipe/recipes/soxs_mbias.py", 12, "E501", "Line too long (256 > 120)")), REPO_ROOT
    )

    # ACT
    rendered = str(findings[0])

    # ASSERT
    assert rendered == "soxspipe/recipes/soxs_mbias.py:12: E501 Line too long (256 > 120)"


def test_gated_findings_are_sorted_by_path_then_line():
    # ARRANGE
    diffText = (
        "+++ b/soxspipe/commonutils/toolkit.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
        "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -1,0 +2,2 @@\n+y = 1\n+z = 2\n"
    )
    changedLines = lint_ratchet.parse_changed_lines(diffText)
    findings = lint_ratchet.parse_ruff_findings(
        _ruff_report(
            ("soxspipe/recipes/soxs_mbias.py", 3, "F401", "Unused import"),
            ("soxspipe/recipes/soxs_mbias.py", 2, "E501", "Line too long"),
            ("soxspipe/commonutils/toolkit.py", 2, "E501", "Line too long"),
        ),
        REPO_ROOT,
    )

    # ACT
    gated = lint_ratchet.findings_on_changed_lines(findings, changedLines)

    # ASSERT
    assert [(finding.path, finding.line) for finding in gated] == [
        ("soxspipe/commonutils/toolkit.py", 2),
        ("soxspipe/recipes/soxs_mbias.py", 2),
        ("soxspipe/recipes/soxs_mbias.py", 3),
    ]


def test_an_empty_ruff_report_parses_to_no_findings():
    # ACT
    findings = lint_ratchet.parse_ruff_findings("", REPO_ROOT)

    # ASSERT
    assert findings == []


def test_ruff_is_not_run_when_the_change_touches_no_python_file():
    # ACT
    ruffJson = lint_ratchet.run_ruff([], REPO_ROOT)

    # ASSERT
    assert json.loads(ruffJson) == []


def test_the_diff_covers_every_extension_ruff_lints():
    # ASSERT
    assert lint_ratchet.LINTABLE_PATHSPECS == ("*.py", "*.pyi", "*.pyw", "*.ipynb")


def test_a_failing_command_raises_rather_than_returning_empty_output(monkeypatch):
    # ARRANGE
    def fail(command, **keywordArguments):
        return _CompletedCommand(returncode=128, stdout="", stderr="fatal: bad revision 'origin/nope'")

    monkeypatch.setattr(lint_ratchet.subprocess, "run", fail)

    # ACT, ASSERT
    with pytest.raises(lint_ratchet.RatchetError) as raised:
        lint_ratchet.run_git_diff("origin/nope", REPO_ROOT)

    assert "bad revision" in str(raised.value)


def test_a_missing_executable_raises_a_ratchet_error(monkeypatch):
    # ARRANGE
    def missing(command, **keywordArguments):
        raise OSError("No such file or directory: 'ruff'")

    monkeypatch.setattr(lint_ratchet.subprocess, "run", missing)

    # ACT, ASSERT
    with pytest.raises(lint_ratchet.RatchetError) as raised:
        lint_ratchet.run_ruff(["soxspipe/recipes/soxs_mbias.py"], REPO_ROOT)

    assert "could not run ruff" in str(raised.value)


def test_main_reports_a_tool_error_when_git_fails(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)

    def fail(compareBranch, repoRoot):
        raise lint_ratchet.RatchetError("git diff failed with status 128")

    monkeypatch.setattr(lint_ratchet, "run_git_diff", fail)

    # ACT
    status = lint_ratchet.main(["--compare-branch", "origin/develop"])

    # ASSERT
    assert status == lint_ratchet.EXIT_TOOL_ERROR
    assert "git diff failed" in capsys.readouterr().err


def test_main_reports_a_tool_error_when_the_ruff_report_is_malformed(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "+++ b/soxspipe/x.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
    )
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: "{not json")

    # ACT
    status = lint_ratchet.main(["--compare-branch", "origin/develop"])

    # ASSERT
    assert status == lint_ratchet.EXIT_TOOL_ERROR
    assert "could not read the ruff report" in capsys.readouterr().err


def test_main_exits_clean_when_no_finding_lands_on_a_changed_line(monkeypatch):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "+++ b/soxspipe/x.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
    )
    report = _ruff_report(("soxspipe/x.py", 99, "E501", "Line too long"))
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: report)

    # ACT
    status = lint_ratchet.main(["--staged"])

    # ASSERT
    assert status == lint_ratchet.EXIT_CLEAN


def test_main_fails_when_a_finding_lands_on_a_changed_line(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "+++ b/soxspipe/x.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
    )
    report = _ruff_report(("soxspipe/x.py", 2, "E501", "Line too long"))
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: report)

    # ACT
    status = lint_ratchet.main(["--staged"])

    # ASSERT
    assert status == lint_ratchet.EXIT_FINDINGS
    assert "soxspipe/x.py:2: E501 Line too long" in capsys.readouterr().out


def test_deletion_only_files_are_never_handed_to_ruff(monkeypatch):
    # ARRANGE
    inspectedPaths = []
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet,
        "run_git_diff",
        lambda compareBranch, repoRoot: (
            "+++ b/soxspipe/kept.py\n@@ -20,3 +19,0 @@\n-x = 1\n" "+++ b/soxspipe/edited.py\n@@ -1,0 +2,1 @@\n+y = 2\n"
        ),
    )

    def recording_ruff(paths, repoRoot):
        inspectedPaths.extend(paths)
        return "[]"

    monkeypatch.setattr(lint_ratchet, "run_ruff", recording_ruff)

    # ACT
    lint_ratchet.main(["--staged"])

    # ASSERT
    assert inspectedPaths == ["soxspipe/edited.py"]
