"""Contracts for the legacy integration-test workflow triggers."""

from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.unit


def test_legacy_integration_runs_for_main_prs_and_weekly_schedule() -> None:
    """Keep the expensive legacy suite off routine feature pull requests."""
    workflowPath = (
        Path(__file__).parents[2] / ".github" / "workflows" / "integration-tests.yml"
    )

    workflow = yaml.load(workflowPath.read_text(), Loader=yaml.BaseLoader)

    assert workflow["on"] == {
        "pull_request": {"branches": ["main"]},
        "schedule": [{"cron": "17 3 * * 1"}],
    }
