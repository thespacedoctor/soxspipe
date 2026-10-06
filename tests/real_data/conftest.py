"""Shared fixtures for the real-data acceptance checks."""

from __future__ import annotations

import os
from pathlib import Path

import pytest


@pytest.fixture
def reduced_workspace() -> Path:
    """Return the workflow-created real-data workspace, or skip when absent."""
    workspace_value = os.environ.get("SOXSPIPE_REAL_DATA_DIR")
    if workspace_value is None:
        pytest.skip("SOXSPIPE_REAL_DATA_DIR is required for real-data tests")
    workspace_path = Path(workspace_value).resolve()
    if not workspace_path.is_dir():
        pytest.skip("SOXSPIPE_REAL_DATA_DIR does not name a workspace directory")
    return workspace_path
