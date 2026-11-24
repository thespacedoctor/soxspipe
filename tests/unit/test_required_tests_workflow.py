"""Contracts for the steps of the required-tests workflow."""

import re
from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.unit

WORKFLOW_PATH = Path(__file__).parents[2] / ".github" / "workflows" / "tests.yml"


CHANGELOG_CHECK = "python tools/check_changelog.py CHANGES.md"


def load_workflow() -> dict:
    """Return the required-tests workflow as nested str/list/dict values."""
    # BASELOADER BUILDS ONLY str/list/dict; safe_load WOULD TURN THE `on` KEY INTO True
    return yaml.load(WORKFLOW_PATH.read_text(), Loader=yaml.BaseLoader)  # noqa: S506


def load_step_commands(job: str = "tests") -> list[str]:
    """Return the `run` command of every step in one job, in order."""
    return [step.get("run", "").strip() for step in load_workflow()["jobs"][job]["steps"]]


def test_required_tests_gate_the_whole_repository_on_ruff_format() -> None:
    """Every file in the tree must stay formatted, not only the changed lines."""
    commands = load_step_commands()

    assert "ruff format --check ." in commands


def test_format_gate_runs_after_install_and_before_the_suite() -> None:
    """The pinned ruff comes from the test extras, and a format failure should end the job early."""
    commands = load_step_commands()
    formatIndex = commands.index("ruff format --check .")

    assert any("pip install" in command for command in commands[:formatIndex])
    assert all("pytest" not in command for command in commands[:formatIndex])


def test_required_suite_runs_under_xdist_on_every_runner_core() -> None:
    """The required suite runs in parallel (`-n auto`) to stay inside the job's 30-minute cap."""
    commands = load_step_commands()
    suiteCommand = next((command for command in commands if "-m pytest" in command), None)

    assert suiteCommand is not None, "no workflow step runs the required suite"
    assert re.search(r"(?:^|\s)-n auto(?:\s|$)", suiteCommand)


def test_changelog_check_runs_in_its_own_ungated_job() -> None:
    """A changelog-only change skips the suite, so the duplicate-entry check must not live in it."""
    workflow = load_workflow()

    assert CHANGELOG_CHECK in load_step_commands("changelog")
    assert "if" not in workflow["jobs"]["changelog"]
    assert CHANGELOG_CHECK not in load_step_commands("tests")


def test_required_suite_runs_only_when_the_changes_job_says_so() -> None:
    """The `tests` job keeps its name (a required check on main) and is skipped by `if`, not by a path filter."""
    workflow = load_workflow()

    assert workflow["jobs"]["tests"]["needs"] == "changes"
    assert workflow["jobs"]["tests"]["if"] == (
        "${{ !cancelled() && (needs.changes.result != 'success' || needs.changes.outputs.code == 'true') }}"
    )
    for trigger in workflow["on"].values():
        assert "paths" not in trigger
        assert "paths-ignore" not in trigger


def test_changes_job_checks_out_the_full_history() -> None:
    """The diff needs both pull request tips and their merge base; a shallow clone would always run the suite."""
    checkoutStep = load_workflow()["jobs"]["changes"]["steps"][0]

    assert checkoutStep["uses"].startswith("actions/checkout@")
    assert checkoutStep["with"]["fetch-depth"] == "0"
