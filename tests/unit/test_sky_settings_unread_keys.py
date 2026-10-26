"""No ``sky-subtraction`` settings block carries a key that no code reads (DY-1287)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
SETTINGS_FILES = [
    "soxs_default_settings.yaml",
    "xsh_default_settings.yaml",
    "test_settings_soxs_sim.yaml",
    "test_settings_xsh.yaml",
]
UNREAD_KEYS = ["skyline_significance"]


def _sky_subtraction_blocks(node):
    """*every ``sky-subtraction`` block nested anywhere in a settings tree*"""
    if not isinstance(node, dict):
        return []
    blocks = [node["sky-subtraction"]] if "sky-subtraction" in node else []
    for value in node.values():
        blocks.extend(_sky_subtraction_blocks(value))
    return blocks


@pytest.mark.parametrize("fileName", SETTINGS_FILES)
def test_sky_subtraction_blocks_carry_no_unread_keys(fileName):
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ACT
    blocks = _sky_subtraction_blocks(settings)

    # ASSERT
    assert blocks
    for block in blocks:
        assert not set(UNREAD_KEYS) & set(block)
