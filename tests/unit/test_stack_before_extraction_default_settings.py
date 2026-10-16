"""Every stare settings block stacks the object frames of a SOF before extraction by default (DY-1611)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare"),
    ("soxs_default_settings.yaml", "soxs-stare-std"),
    ("xsh_default_settings.yaml", "soxs-stare"),
    ("test_settings_soxs_sim.yaml", "soxs-stare"),
    ("test_settings_xsh.yaml", "soxs-stare"),
]


@pytest.mark.parametrize(("fileName", "recipe"), CASES)
def test_each_stare_block_stacks_before_extraction(fileName, recipe):
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ASSERT
    assert settings[recipe]["stack_before_extraction"] is True
