"""Contracts for the real-data acceptance workflow triggers."""

from pathlib import Path

import pytest
import yaml

pytestmark = pytest.mark.unit

WORKFLOW_PATH = Path(__file__).parents[2] / ".github" / "workflows" / "real-data-tests.yml"


def load_workflow() -> dict:
    """Load the workflow with every scalar kept as a string, so `on` stays a key."""
    # BASELOADER BUILDS ONLY str/list/dict; safe_load WOULD TURN THE `on` KEY INTO True
    return yaml.load(WORKFLOW_PATH.read_text(), Loader=yaml.BaseLoader)  # noqa: S506


def test_real_data_runs_only_for_main_prs_and_manual_dispatch() -> None:
    """Keep the 120-minute real-data gate off develop pull requests and schedules."""
    triggers = load_workflow()["on"]

    assert set(triggers) == {"pull_request", "workflow_dispatch"}
    assert triggers["pull_request"]["branches"] == ["main"]
    # A PULL REQUEST CREATED READY NEVER FIRES ready_for_review, SO opened IS NEEDED TOO
    assert {"opened", "ready_for_review"} <= set(triggers["pull_request"]["types"])


def test_real_data_still_skips_draft_and_fork_pull_requests() -> None:
    """The job-level guard keeps drafts and forks off the runner but lets dispatch through."""
    guard = load_workflow()["jobs"]["real-data-tests"]["if"]

    assert "github.event_name != 'pull_request'" in guard
    assert "github.event.pull_request.draft == false" in guard
    assert "github.event.pull_request.head.repo.fork == false" in guard
