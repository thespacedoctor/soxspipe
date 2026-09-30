"""*tests that the hard lint rules really are hard, by running ruff for real*

These belong in the integration suite rather than beside the ratchet's unit tests: each
one shells out to ruff, and one of them scans the whole repository.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

pytestmark = pytest.mark.integration

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


def _hard_rule_findings_for(paths: list[str], repoRoot: Path) -> list:
    """*run the real hard-rule scan over the given paths*

    **Key Arguments:**

    - ``paths`` -- the paths to scan
    - ``repoRoot`` -- the directory to run ruff in

    **Return:**

    - ``findings`` -- the hard-rule findings ruff reported
    """
    report = lint_ratchet.run_ruff_selected(lint_ratchet.HARD_RULE_CODES, paths, repoRoot)

    return lint_ratchet.hard_rule_findings(lint_ratchet.parse_ruff_findings(report, repoRoot))


def test_the_repository_itself_passes_the_hard_rules():
    # ARRANGE
    repoRoot = lint_ratchet.repository_root(REPO_ROOT)

    # ACT
    findings = _hard_rule_findings_for(list(lint_ratchet.HARD_RULE_PATHS), repoRoot)

    # ASSERT
    assert findings == []


def test_the_hard_rules_reach_every_tracked_python_file():
    # ARRANGE
    repoRoot = lint_ratchet.repository_root(REPO_ROOT)
    tracked = {(repoRoot / path).as_posix() for path in lint_ratchet._run(["git", "ls-files", "*.py"], repoRoot).split()}

    # ACT
    command = [*lint_ratchet.ruff_command(), "check", "--show-files", "--force-exclude"]
    command += ["--", *lint_ratchet.HARD_RULE_PATHS]
    scanned = {Path(line).as_posix() for line in lint_ratchet._run(command, repoRoot, allowedStatuses=(0, 1)).splitlines()}

    # ASSERT
    assert tracked - scanned == set()


def test_a_noqa_comment_cannot_silence_a_hard_rule(tmp_path):
    # ARRANGE
    offender = tmp_path / "suppressed.py"
    offender.write_text("def f():\n    try:\n        pass\n    except:  # noqa: E722, S110\n        pass\n")

    # ACT
    findings = _hard_rule_findings_for([offender.name], tmp_path)

    # ASSERT
    assert sorted({finding.code for finding in findings}) == ["E722", "S110"]


def test_a_per_file_ignore_cannot_silence_a_hard_rule(tmp_path):
    # ARRANGE
    offender = tmp_path / "ignored.py"
    offender.write_text("def f():\n    try:\n        pass\n    except:\n        pass\n")
    (tmp_path / "pyproject.toml").write_text(
        '[tool.ruff.lint.per-file-ignores]\n"ignored.py" = ["E722", "S110"]\n',
    )

    # ACT
    findings = _hard_rule_findings_for([offender.name], tmp_path)

    # ASSERT
    assert sorted({finding.code for finding in findings}) == ["E722", "S110"]
