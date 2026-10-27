"""Each stare arm's ``sky-subtraction`` block turns the object-wing correction on and sets its gates."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
MIN_LEAK_SIGMA = 0.03
MIN_SKY_SIGMA = 1.0
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare", "vis"),
    ("soxs_default_settings.yaml", "soxs-stare", "nir"),
    ("soxs_default_settings.yaml", "soxs-stare-std", "vis"),
    ("soxs_default_settings.yaml", "soxs-stare-std", "nir"),
    ("xsh_default_settings.yaml", "soxs-stare", "uvb"),
    ("xsh_default_settings.yaml", "soxs-stare", "vis"),
    ("xsh_default_settings.yaml", "soxs-stare", "nir"),
]


@pytest.mark.parametrize(("fileName", "recipe", "arm"), CASES)
def test_each_stare_arm_turns_the_object_wing_correction_on(fileName, recipe, arm):
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    block = settings[recipe][arm]["sky-subtraction"]

    assert block["object_wing_correction"] is True
    assert block["object_wing_min_leak_sigma"] == MIN_LEAK_SIGMA
    assert block["object_wing_min_sky_sigma"] == MIN_SKY_SIGMA
    # THE SKY-FRACTION GATE SKIPPED REAL LEAKS OF 4 TO 7 PERCENT OF THE SKY AND IS GONE (DY-1964)
    assert "object_wing_gate_fraction" not in block


@pytest.mark.parametrize("fileName", ["test_settings_soxs_sim.yaml", "test_settings_xsh.yaml"])
def test_the_test_settings_carry_the_object_wing_keys(fileName):
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    blocks = []

    def collect(node):
        if isinstance(node, dict):
            if "sky-subtraction" in node:
                blocks.append(node["sky-subtraction"])
            for value in node.values():
                collect(value)

    collect(settings)

    assert blocks
    for block in blocks:
        assert block["object_wing_correction"] is True
        assert block["object_wing_min_leak_sigma"] == MIN_LEAK_SIGMA
        assert block["object_wing_min_sky_sigma"] == MIN_SKY_SIGMA
        # THE SKY-FRACTION GATE SKIPPED REAL LEAKS OF 4 TO 7 PERCENT OF THE SKY AND IS GONE (DY-1964)
        assert "object_wing_gate_fraction" not in block
