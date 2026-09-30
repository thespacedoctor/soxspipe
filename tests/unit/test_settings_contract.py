"""Characterization tests for existing required settings."""

from __future__ import annotations

import pytest

from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit


def test_base_recipe_requires_workspace_root_before_side_effects(log: object) -> None:
    with pytest.raises(KeyError, match="workspace-root-dir"):
        base_recipe(log=log, settings={}, inputFrames=False)
