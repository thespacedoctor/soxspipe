"""Each stare arm has a ``slit-illumination`` default block: on for SOXS, off for XSH (DY-1280)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
BLOCK_NUMBERS = {
    "overlap-margin-arcsec": 0.3,
    "profile-bin-arcsec": 0.1,
    "max-offset-arcsec": 3.0,
    "max-order-scatter-arcsec": 0.3,
}
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare", "vis", True),
    ("soxs_default_settings.yaml", "soxs-stare", "nir", True),
    ("soxs_default_settings.yaml", "soxs-stare-std", "vis", True),
    ("soxs_default_settings.yaml", "soxs-stare-std", "nir", True),
    ("xsh_default_settings.yaml", "soxs-stare", "uvb", False),
    ("xsh_default_settings.yaml", "soxs-stare", "vis", False),
    ("xsh_default_settings.yaml", "soxs-stare", "nir", False),
]


@pytest.mark.parametrize(("fileName", "recipe", "arm", "isOn"), CASES)
def test_each_stare_arm_has_the_slit_illumination_block(fileName, recipe, arm, isOn):
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    block = settings[recipe][arm]["slit-illumination"]

    assert block["measure-offset"] is isOn
    for key, value in BLOCK_NUMBERS.items():
        assert block[key] == value
    assert set(block) == {"measure-offset", *BLOCK_NUMBERS}


def test_no_slit_illumination_qc_has_an_acceptable_range():
    for fileName in ("soxs_default_settings.yaml", "xsh_default_settings.yaml"):
        settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())
        for recipeSettings in settings.values():
            if not isinstance(recipeSettings, dict):
                continue
            for armSettings in recipeSettings.values():
                ranges = armSettings.get("qc-acceptable-ranges", {}) if isinstance(armSettings, dict) else {}
                assert not [name for name in ranges if "illum" in name.lower()]
