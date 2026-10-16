"""Every nod and offset recipe block ships ``stack_before_extraction: True`` (DY-1612)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
CASES = [
    ("soxs_default_settings.yaml", "soxs-nod"),
    ("soxs_default_settings.yaml", "soxs-nod-std"),
    ("soxs_default_settings.yaml", "soxs-offset"),
    ("soxs_default_settings.yaml", "soxs-offset-std"),
    ("xsh_default_settings.yaml", "soxs-nod"),
    ("test_settings_xsh.yaml", "soxs-nod"),
]


@pytest.mark.parametrize(("fileName", "recipe"), CASES)
def test_each_nod_and_offset_block_stacks_before_extraction_by_default(fileName, recipe):
    # ARRANGE
    settingsText = (PACKAGE_DIR / fileName).read_text()

    # ACT
    settings = yaml.safe_load(settingsText)

    # ASSERT
    assert settings[recipe]["stack_before_extraction"] is True
