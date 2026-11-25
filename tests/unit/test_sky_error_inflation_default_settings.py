"""The shipped ``sky_error_inflation`` defaults: on for SOXS VIS and NIR stare, off for X-Shooter VIS (DY-2436)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.subtract_sky import ERROR_INFLATION_ARMS

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
KEY = "sky_error_inflation"
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare", "vis", True),
    ("soxs_default_settings.yaml", "soxs-stare-std", "vis", True),
    ("soxs_default_settings.yaml", "soxs-stare", "nir", True),
    ("soxs_default_settings.yaml", "soxs-stare-std", "nir", True),
    ("xsh_default_settings.yaml", "soxs-stare", "vis", False),
]
# THE VIS SETTING IS A STOPGAP FOR THE LOW VIS ERRS PLANE. NOD AND OFFSET NEVER READ THE KEY, SO THEY CARRY NONE
NOD_AND_OFFSET_RECIPES = ["soxs-nod", "soxs-nod-std", "soxs-offset", "soxs-offset-std"]


@pytest.mark.parametrize(("fileName", "recipe", "arm", "isOn"), CASES)
def test_each_stare_arm_sets_the_sky_error_inflation_explicitly(fileName, recipe, arm, isOn):
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    value = settings[recipe][arm]["sky-subtraction"][KEY]

    assert value is isOn


def test_the_missing_key_fallback_matches_the_shipped_soxs_defaults():
    assert frozenset({"VIS", "NIR"}) == ERROR_INFLATION_ARMS


def test_no_nod_or_offset_block_carries_the_key():
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())

    for recipe in NOD_AND_OFFSET_RECIPES:
        for armSettings in settings[recipe].values():
            if isinstance(armSettings, dict):
                assert KEY not in armSettings.get("sky-subtraction", {})
