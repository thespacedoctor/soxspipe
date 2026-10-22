"""Every ``detect-continuum`` block sets the length of sample-free trace that triggers a warning (DY-1716)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.detect_continuum import DEFAULT_MAX_UNSAMPLED_TRACE_FRACTION

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
KEY = "max-unsampled-trace-fraction"


def _detect_continuum_blocks(node: object) -> list[dict]:
    """Return every ``detect-continuum`` block in a nested settings tree."""
    blocks = []
    if isinstance(node, dict):
        if "detect-continuum" in node:
            blocks.append(node["detect-continuum"])
        for value in node.values():
            blocks.extend(_detect_continuum_blocks(value))
    return blocks


@pytest.mark.parametrize("fileName", ["soxs_default_settings.yaml", "xsh_default_settings.yaml"])
def test_each_detect_continuum_block_carries_the_extrapolation_threshold(fileName: str) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ACT
    blocks = _detect_continuum_blocks(settings)

    # ASSERT: THE FALLBACK IN THE CODE MATCHES THE SHIPPED VALUE
    assert blocks
    for block in blocks:
        assert block[KEY] == DEFAULT_MAX_UNSAMPLED_TRACE_FRACTION
