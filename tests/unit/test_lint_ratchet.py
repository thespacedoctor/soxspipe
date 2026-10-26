"""*tests for the ruff lint ratchet's diff and finding intersection*

The intersection logic is pure, so these tests feed it diff text and ruff JSON directly
rather than shelling out to git or ruff.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
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


@pytest.fixture(autouse=True)
def stub_the_hard_rule_scan(monkeypatch):
    """*keep the unit suite hermetic: `main` would otherwise run ruff over the whole tree*

    A test that needs the hard-rule path sets its own stub, which replaces this one.
    """
    monkeypatch.setattr(lint_ratchet, "run_ruff_selected", lambda codes, paths, repoRoot: "[]")


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
    diffText = "+++ b/soxspipe/recipes/soxs_mbias.py\n@@ -10,0 +11,2 @@\n+x = 1\n+y = 2\n@@ -40,0 +43,1 @@\n+z = 3\n"

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
            "+++ b/soxspipe/kept.py\n@@ -20,3 +19,0 @@\n-x = 1\n+++ b/soxspipe/edited.py\n@@ -1,0 +2,1 @@\n+y = 2\n"
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


def test_hard_rule_findings_keeps_only_the_absolute_rules():
    # ARRANGE
    findings = [
        lint_ratchet.Finding("soxspipe/a.py", 10, "E722", "Do not use bare `except`"),
        lint_ratchet.Finding("soxspipe/b.py", 20, "E501", "Line too long"),
        lint_ratchet.Finding("tools/c.py", 30, "S110", "`try`-`except`-`pass` detected"),
    ]

    # ACT
    hard = lint_ratchet.hard_rule_findings(findings)

    # ASSERT
    assert [finding.code for finding in hard] == ["E722", "S110"]


def test_hard_rule_findings_are_sorted_by_path_and_line():
    # ARRANGE
    findings = [
        lint_ratchet.Finding("tools/c.py", 30, "S112", "`try`-`except`-`continue` detected"),
        lint_ratchet.Finding("soxspipe/a.py", 40, "E722", "Do not use bare `except`"),
        lint_ratchet.Finding("soxspipe/a.py", 10, "E722", "Do not use bare `except`"),
    ]

    # ACT
    hard = lint_ratchet.hard_rule_findings(findings)

    # ASSERT
    assert [(finding.path, finding.line) for finding in hard] == [
        ("soxspipe/a.py", 10),
        ("soxspipe/a.py", 40),
        ("tools/c.py", 30),
    ]


def test_main_fails_on_a_bare_except_the_change_never_touched(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "+++ b/soxspipe/x.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
    )
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: "[]")
    untouched = _ruff_report(("soxspipe/elsewhere.py", 900, "E722", "Do not use bare `except`"))
    monkeypatch.setattr(lint_ratchet, "run_ruff_selected", lambda codes, paths, repoRoot: untouched)

    # ACT
    status = lint_ratchet.main(["--staged"])

    # ASSERT
    assert status == lint_ratchet.EXIT_FINDINGS
    assert "soxspipe/elsewhere.py:900: E722" in capsys.readouterr().out


def test_main_names_the_hard_rules_when_one_fails(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "")
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: "[]")
    untouched = _ruff_report(("tools/thing.py", 12, "S110", "`try`-`except`-`pass` detected"))
    monkeypatch.setattr(lint_ratchet, "run_ruff_selected", lambda codes, paths, repoRoot: untouched)

    # ACT
    lint_ratchet.main(["--staged"])

    # ASSERT
    output = capsys.readouterr().out
    assert "tools/thing.py:12: S110" in output
    assert "1 hard rule finding" in output


def test_main_checks_the_hard_rules_over_the_whole_repository(monkeypatch):
    # ARRANGE
    requested = {}
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "")
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: "[]")

    def recording_selected(codes, paths, repoRoot):
        requested["codes"] = codes
        requested["paths"] = paths
        return "[]"

    monkeypatch.setattr(lint_ratchet, "run_ruff_selected", recording_selected)

    # ACT
    status = lint_ratchet.main(["--staged"])

    # ASSERT
    assert status == lint_ratchet.EXIT_CLEAN
    assert set(requested["codes"]) == {"E722", "S110", "S112"}
    # THE REPOSITORY ROOT, SO A PYTHON FILE OUTSIDE soxspipe, tests AND tools IS GATED TOO
    assert requested["paths"] == ["."]


def test_a_hard_rule_finding_on_a_changed_line_is_not_reported_twice(monkeypatch, capsys):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)
    monkeypatch.setattr(
        lint_ratchet, "run_git_diff", lambda compareBranch, repoRoot: "+++ b/soxspipe/x.py\n@@ -1,0 +2,1 @@\n+x = 1\n"
    )
    bareExcept = _ruff_report(("soxspipe/x.py", 2, "E722", "Do not use bare `except`"))
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: bareExcept)
    monkeypatch.setattr(lint_ratchet, "run_ruff_selected", lambda codes, paths, repoRoot: bareExcept)

    # ACT
    status = lint_ratchet.main(["--staged"])

    # ASSERT
    output = capsys.readouterr().out
    assert status == lint_ratchet.EXIT_FINDINGS
    assert output.count("soxspipe/x.py:2: E722") == 1
    assert "No findings on changed lines." in output


def test_ruff_is_run_through_the_interpreter_when_it_is_installed_beside_it(monkeypatch):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet.importlib.util, "find_spec", lambda name: object())

    # ACT
    command = lint_ratchet.ruff_command()

    # ASSERT
    assert command == [lint_ratchet.sys.executable, "-m", "ruff"]


def test_ruff_falls_back_to_the_path_when_it_is_not_installed_beside_the_interpreter(monkeypatch):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet.importlib.util, "find_spec", lambda name: None)

    # ACT
    command = lint_ratchet.ruff_command()

    # ASSERT
    assert command == ["ruff"]


def test_the_ruff_runners_do_not_depend_on_the_path(monkeypatch):
    # ARRANGE
    # A FRESH MODULE, SO THE AUTOUSE STUB OF run_ruff_selected DOES NOT HIDE ITS COMMAND
    ratchet = _load_lint_ratchet()
    calls = []

    def _record(command, **keywordArguments):
        calls.append(command)

        return _CompletedCommand(returncode=0, stdout="[]", stderr="")

    monkeypatch.setattr(ratchet.subprocess, "run", _record)
    monkeypatch.setattr(ratchet, "ruff_command", lambda: ["/somewhere/python", "-m", "ruff"])

    # ACT
    ratchet.run_ruff(["soxspipe/recipes/soxs_mbias.py"], REPO_ROOT)
    ratchet.run_ruff_selected(ratchet.HARD_RULE_CODES, ["."], REPO_ROOT)

    # ASSERT
    assert [command[:4] for command in calls] == [
        ["/somewhere/python", "-m", "ruff", "check"],
        ["/somewhere/python", "-m", "ruff", "check"],
    ]


GIT_IDENTITY = {
    "GIT_AUTHOR_NAME": "Test",
    "GIT_AUTHOR_EMAIL": "test@example.com",
    "GIT_COMMITTER_NAME": "Test",
    "GIT_COMMITTER_EMAIL": "test@example.com",
}


def _git(repoPath: Path, *arguments: str) -> None:
    """*run git in a temporary repository with a fixed identity*

    **Key Arguments:**

    - ``repoPath`` -- the repository to run git in
    - ``arguments`` -- the git arguments
    """
    environment = {**os.environ, **GIT_IDENTITY}
    subprocess.run(["git", *arguments], cwd=repoPath, env=environment, check=True, capture_output=True)  # noqa: S603, S607


@pytest.fixture
def worktree_repo(tmp_path) -> Path:
    """*a git repository with a `main` branch holding a.py and a .gitignore, checked out on `feature`*

    **Return:**

    - ``repoPath`` -- the repository root
    """
    repoPath = tmp_path / "repo"
    repoPath.mkdir()
    _git(repoPath, "init", "-b", "main")
    (repoPath / "a.py").write_text("one = 1\ntwo = 2\nthree = 3\n", encoding="utf-8")
    (repoPath / ".gitignore").write_text("ignored.py\n", encoding="utf-8")
    _git(repoPath, "add", "-A")
    _git(repoPath, "commit", "-m", "base")
    _git(repoPath, "checkout", "-b", "feature")

    return repoPath


def test_worktree_lines_marks_an_uncommitted_edit_as_changed(worktree_repo):
    # ARRANGE
    (worktree_repo / "a.py").write_text("one = 1\ntwo = 22\nthree = 3\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"a.py": {2}}


def test_worktree_lines_marks_an_untracked_file_as_fully_added(worktree_repo):
    # ARRANGE
    (worktree_repo / "new.py").write_text("x = 1\ny = 2\nz = 3\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"new.py": {1, 2, 3}}


def test_worktree_lines_skips_an_ignored_file(worktree_repo):
    # ARRANGE
    (worktree_repo / "ignored.py").write_text("x = 1\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {}


def test_worktree_lines_includes_committed_work_since_the_merge_base(worktree_repo):
    # ARRANGE
    (worktree_repo / "a.py").write_text("one = 1\ntwo = 2\nthree = 33\n", encoding="utf-8")
    _git(worktree_repo, "commit", "-am", "edit three")
    (worktree_repo / "a.py").write_text("one = 11\ntwo = 2\nthree = 33\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"a.py": {1, 3}}


def test_worktree_lines_ignores_untracked_files_ruff_does_not_lint(worktree_repo):
    # ARRANGE
    (worktree_repo / "notes.txt").write_text("hello\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {}


def test_worktree_lines_gives_an_empty_untracked_file_no_changed_lines(worktree_repo):
    # ARRANGE
    (worktree_repo / "empty.py").write_text("", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {}


def test_worktree_lines_raises_when_the_base_branch_does_not_exist(worktree_repo):
    # ACT, ASSERT
    with pytest.raises(lint_ratchet.RatchetError):
        lint_ratchet.worktree_changed_lines("origin/nope", worktree_repo)


def test_main_gates_a_finding_in_an_untracked_file_in_worktree_mode(worktree_repo, monkeypatch, capsys):
    # ARRANGE
    (worktree_repo / "new.py").write_text("x = 1\n", encoding="utf-8")
    monkeypatch.chdir(worktree_repo)
    report = json.dumps(
        [
            {
                "filename": str(worktree_repo.resolve() / "new.py"),
                "location": {"row": 1, "column": 1},
                "end_location": {"row": 1, "column": 2},
                "code": "E501",
                "message": "Line too long",
            }
        ]
    )
    monkeypatch.setattr(lint_ratchet, "run_ruff", lambda paths, repoRoot: report)

    # ACT
    status = lint_ratchet.main(["--worktree", "main"])

    # ASSERT
    assert status == lint_ratchet.EXIT_FINDINGS
    assert "new.py:1: E501 Line too long" in capsys.readouterr().out


def test_worktree_mode_cannot_be_combined_with_the_other_modes():
    # ACT, ASSERT
    with pytest.raises(SystemExit):
        lint_ratchet.main(["--worktree", "origin/develop", "--staged"])
    with pytest.raises(SystemExit):
        lint_ratchet.main(["--worktree", "origin/develop", "--compare-branch", "origin/develop"])


def test_worktree_lines_skips_an_untracked_path_that_is_not_a_regular_file(worktree_repo):
    # ARRANGE
    (worktree_repo / "dangling.py").symlink_to(worktree_repo / "nowhere.py")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {}


def test_worktree_lines_skip_a_symlink_to_an_existing_file_outside_the_repository(worktree_repo, tmp_path):
    # ARRANGE
    outside = tmp_path / "outside.py"
    outside.write_text("x = 1\ny = 2\n", encoding="utf-8")
    (worktree_repo / "link.py").symlink_to(outside)

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {}


def test_worktree_lines_keep_an_untracked_filename_that_contains_a_newline(worktree_repo):
    # ARRANGE
    oddName = "we\nird.py"
    try:
        (worktree_repo / oddName).write_text("x = 1\ny = 2\n", encoding="utf-8")
    except OSError:
        pytest.skip("this filesystem does not allow a newline in a filename")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {oddName: {1, 2}}


def test_worktree_lines_keep_an_untracked_filename_with_non_ascii_characters(worktree_repo):
    # ARRANGE
    (worktree_repo / "naïve.py").write_text("x = 1\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"naïve.py": {1}}


def test_worktree_lines_merge_an_edit_and_an_untracked_file(worktree_repo):
    # ARRANGE
    (worktree_repo / "a.py").write_text("one = 1\ntwo = 22\nthree = 3\n", encoding="utf-8")
    (worktree_repo / "new.py").write_text("x = 1\n", encoding="utf-8")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"a.py": {2}, "new.py": {1}}


def test_worktree_lines_count_lines_by_the_endings_ruff_counts(worktree_repo):
    # ARRANGE
    # A FORM FEED, A VERTICAL TAB, FILE SEPARATORS AND NEXT-LINE ARE NOT LINE ENDINGS FOR RUFF
    (worktree_repo / "odd.py").write_bytes(b"x = 1\x0c\x0b\x1c\x1d\x1e\xc2\x85y = 2\n")
    (worktree_repo / "mixed.py").write_bytes(b"a = 1\r\nb = 2\rc = 3\nd = 4")

    # ACT
    changedLines = lint_ratchet.worktree_changed_lines("main", worktree_repo)

    # ASSERT
    assert changedLines == {"odd.py": {1}, "mixed.py": {1, 2, 3, 4}}


@pytest.mark.parametrize("mode", [["--staged"], ["--worktree", "main"]])
def test_main_reports_a_tool_error_when_the_git_output_is_not_utf8(monkeypatch, capsys, mode):
    # ARRANGE
    monkeypatch.setattr(lint_ratchet, "repository_root", lambda start: REPO_ROOT)

    def undecodable(*arguments):
        raise UnicodeDecodeError("utf-8", b"\xff", 0, 1, "invalid start byte")

    monkeypatch.setattr(lint_ratchet, "run_git_diff", undecodable)
    monkeypatch.setattr(lint_ratchet, "run_worktree_diff", undecodable)

    # ACT
    status = lint_ratchet.main(mode)

    # ASSERT
    errorLines = capsys.readouterr().err.strip().splitlines()
    assert status == lint_ratchet.EXIT_TOOL_ERROR
    assert len(errorLines) == 1
    assert "not valid UTF-8" in errorLines[0]
