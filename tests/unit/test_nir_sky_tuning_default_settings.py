"""The shipped SOXS NIR extraction recipes carry a valid sky-tuning block, and no other settings file does (DY-1968)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.nir_sky_refit import DEFAULT_SETTINGS, read_sky_tuning_settings

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
EXTRACTION_RECIPES = ["soxs-stare", "soxs-stare-std", "soxs-nod", "soxs-nod-std", "soxs-offset", "soxs-offset-std"]


def _settings(fileName: str) -> dict:
    return yaml.safe_load((PACKAGE_DIR / fileName).read_text())


def _blocks(node) -> list:
    """*every ``sky-tuning`` block nested anywhere in a settings tree*"""
    if not isinstance(node, dict):
        return []
    blocks = [node["sky-tuning"]] if "sky-tuning" in node else []
    for value in node.values():
        blocks.extend(_blocks(value))
    return blocks


@pytest.mark.parametrize("recipe", EXTRACTION_RECIPES)
def test_every_soxs_nir_extraction_recipe_ships_the_default_sky_tuning_block(recipe):
    # ARRANGE
    nirSettings = _settings("soxs_default_settings.yaml")[recipe]["nir"]

    # ACT
    block = nirSettings["sky-tuning"]

    # ASSERT
    assert block == DEFAULT_SETTINGS
    assert read_sky_tuning_settings(block).globalFit is True


@pytest.mark.parametrize("recipe", EXTRACTION_RECIPES)
def test_no_soxs_vis_arm_carries_a_sky_tuning_block(recipe):
    # ARRANGE
    settings = _settings("soxs_default_settings.yaml")[recipe]

    # ACT
    blocks = _blocks(settings["vis"])

    # ASSERT
    assert blocks == []


def test_the_soxs_default_settings_carry_exactly_one_block_per_extraction_recipe():
    # ACT
    blocks = _blocks(_settings("soxs_default_settings.yaml"))

    # ASSERT
    assert len(blocks) == len(EXTRACTION_RECIPES)


@pytest.mark.parametrize("fileName", ["xsh_default_settings.yaml", "test_settings_xsh.yaml"])
def test_xshooter_settings_carry_no_sky_tuning_block(fileName):
    # ACT
    blocks = _blocks(_settings(fileName))

    # ASSERT
    assert blocks == []
