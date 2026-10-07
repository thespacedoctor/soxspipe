"""Contracts for the steps of the required-tests workflow."""

from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.unit

WORKFLOW_PATH = Path(__file__).parents[2] / ".github" / "workflows" / "tests.yml"


def load_step_commands() -> list[str]:
    """Return the `run` command of every step in the required-tests job, in order."""
    # BASELOADER BUILDS ONLY str/list/dict; safe_load WOULD TURN THE `on` KEY INTO True
    workflow = yaml.load(WORKFLOW_PATH.read_text(), Loader=yaml.BaseLoader)  # noqa: S506
    return [step.get("run", "").strip() for step in workflow["jobs"]["tests"]["steps"]]


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
