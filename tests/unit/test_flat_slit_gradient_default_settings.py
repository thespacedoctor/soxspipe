"""The shipped ``flat-slit-gradient`` defaults: on for the SOXS NIR stare blocks, off everywhere else (DY-2147)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
KEY = "flat-slit-gradient"
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare", "nir", True),
    ("soxs_default_settings.yaml", "soxs-stare-std", "nir", True),
    ("soxs_default_settings.yaml", "soxs-stare", "vis", False),
    ("soxs_default_settings.yaml", "soxs-stare-std", "vis", False),
    ("xsh_default_settings.yaml", "soxs-stare", "uvb", False),
    ("xsh_default_settings.yaml", "soxs-stare", "vis", False),
    ("xsh_default_settings.yaml", "soxs-stare", "nir", False),
]
# NOD AND OFFSET REMOVE THE SKY AT THE SAME PIXEL, SO THE GRADIENT CANCELS AND THEY CARRY NO KEY
NOD_AND_OFFSET_RECIPES = ["soxs-nod", "soxs-nod-std", "soxs-offset", "soxs-offset-std"]


@pytest.mark.parametrize(("fileName", "recipe", "arm", "isOn"), CASES)
def test_each_stare_arm_sets_the_flat_slit_gradient_explicitly(fileName, recipe, arm, isOn):
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    value = settings[recipe][arm]["sky-subtraction"][KEY]

    assert value is isOn


def test_no_nod_or_offset_block_carries_the_key():
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())

    for recipe in NOD_AND_OFFSET_RECIPES:
        for armSettings in settings[recipe].values():
            if isinstance(armSettings, dict):
                assert KEY not in armSettings.get("sky-subtraction", {})
