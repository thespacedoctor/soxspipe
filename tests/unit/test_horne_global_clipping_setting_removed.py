"""No shipped settings file carries the removed ``horne-extraction-profile-global-clipping-sigma`` key (DY-1976).

Nothing read the key, so it was removed. A workspace ``soxspipe.yaml`` that still holds it is ignored: see
``test_constructor_ignores_the_removed_global_clipping_sigma_key_in_a_workspace_settings_file``.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
SETTINGS_FILES = sorted(path.name for path in PACKAGE_DIR.glob("*settings*.yaml"))
REMOVED_KEY = "horne-extraction-profile-global-clipping-sigma"


def _keys(node: object) -> set[str]:
    """*every mapping key nested anywhere in a settings tree*"""
    if not isinstance(node, dict):
        return set()
    keys = set(node)
    for value in node.values():
        keys |= _keys(value)
    return keys


def test_settings_files_are_found() -> None:
    assert "soxs_default_settings.yaml" in SETTINGS_FILES
    assert "xsh_default_settings.yaml" in SETTINGS_FILES


@pytest.mark.parametrize("fileName", SETTINGS_FILES)
def test_shipped_settings_carry_no_global_clipping_sigma_key(fileName: str) -> None:
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    assert REMOVED_KEY not in _keys(settings)
