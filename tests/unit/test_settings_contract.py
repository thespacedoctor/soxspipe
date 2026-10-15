"""Characterization tests for existing required settings."""

from __future__ import annotations

import pytest

from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit


def test_the_packaged_advanced_settings_default_the_workspace_root_to_the_current_directory() -> None:
    """Session settings name no workspace root, so the packaged default must (DY-1601)."""
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)

    # ACT
    advancedSettings = recipe._advanced_settings()

    # ASSERT
    assert advancedSettings["workspace-root-dir"] == "./"
