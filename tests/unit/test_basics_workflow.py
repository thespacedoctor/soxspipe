"""Contracts for the basic Python checks workflow triggers."""

from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.unit

WORKFLOW_PATH = Path(__file__).parents[2] / ".github" / "workflows" / "basics.yml"


def test_basic_checks_skip_changes_that_touch_only_the_changelog() -> None:
    """A change to CHANGES.md alone has no Python to check, so the workflow does not start."""
    # BASELOADER BUILDS ONLY str/list/dict; safe_load WOULD TURN THE `on` KEY INTO True
    workflow = yaml.load(WORKFLOW_PATH.read_text(), Loader=yaml.BaseLoader)  # noqa: S506

    assert workflow["on"] == {
        "push": {"paths-ignore": ["CHANGES.md"]},
        "pull_request": {"paths-ignore": ["CHANGES.md"]},
    }
