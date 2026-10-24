"""The faint-order trace settings: code fallbacks and the shipped blocks (DY-1941)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.faint_order_trace import FAINT_TRACE_DEFAULTS, faint_trace_settings

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
ON_RECIPES = ("soxs-stare", "soxs-nod", "soxs-offset")
OFF_RECIPES = ("soxs-stare-std", "soxs-nod-std", "soxs-offset-std", "soxs-order-centres")


def _shipped_block(recipe: str, arm: str) -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings[recipe][arm]["detect-continuum"]


def test_from_block_is_off_with_the_shipped_numbers_when_the_block_has_no_faint_keys() -> None:
    # ARRANGE: A SETTINGS FILE WRITTEN BEFORE THE FEATURE EXISTED
    block = {"slice-length": 70}

    # ACT
    settings = faint_trace_settings.from_block(block)

    # ASSERT
    assert settings.isEnabled is False
    assert settings.minDetectionFraction == 0.5
    assert settings.maxScatterPx == 1.5
    assert settings.maxAnchorScatterPx == 3.0
    assert settings.binCount == 15
    assert settings.minBinSnr == 3.0
    assert settings.maxJoinArcsec == 2.0


def test_from_block_reads_every_faint_key() -> None:
    # ARRANGE
    block = {
        "faint-order-trace": True,
        "faint-order-min-detection-fraction": 0.4,
        "faint-order-max-scatter-px": 2.5,
        "faint-order-max-anchor-scatter-px": 4.5,
        "faint-order-bin-count": 9,
        "faint-order-min-bin-snr": 4.0,
        "same-source-max-join-arcsec": 3.0,
    }

    # ACT
    settings = faint_trace_settings.from_block(block)

    # ASSERT
    assert settings == faint_trace_settings(
        isEnabled=True,
        minDetectionFraction=0.4,
        maxScatterPx=2.5,
        maxAnchorScatterPx=4.5,
        binCount=9,
        minBinSnr=4.0,
        maxJoinArcsec=3.0,
    )


def test_from_block_is_off_when_the_block_is_missing() -> None:
    # ACT
    settings = faint_trace_settings.from_block(None)

    # ASSERT
    assert settings.isEnabled is False


@pytest.mark.parametrize("recipe", ON_RECIPES)
def test_the_science_vis_blocks_switch_the_trace_on(recipe: str) -> None:
    assert _shipped_block(recipe, "vis")["faint-order-trace"] is True


@pytest.mark.parametrize("recipe", ON_RECIPES)
def test_the_science_vis_blocks_carry_the_anchor_scatter_limit(recipe: str) -> None:
    assert _shipped_block(recipe, "vis")["faint-order-max-anchor-scatter-px"] == 3.0


def test_the_code_fallbacks_include_the_anchor_scatter_limit() -> None:
    assert FAINT_TRACE_DEFAULTS["faint-order-max-anchor-scatter-px"] == 3.0


@pytest.mark.parametrize("recipe", ON_RECIPES + OFF_RECIPES)
@pytest.mark.parametrize("arm", ["vis", "nir"])
def test_every_block_carries_the_shipped_numbers_or_none_at_all(recipe: str, arm: str) -> None:
    # ARRANGE
    block = _shipped_block(recipe, arm)

    # ACT
    settings = faint_trace_settings.from_block(block)

    # ASSERT: A SHIPPED KEY HAS THE VALUE THE CODE FALLBACK USES
    for key, default in FAINT_TRACE_DEFAULTS.items():
        if key in block and key != "faint-order-trace":
            assert block[key] == default
    assert settings.isEnabled is (arm == "vis" and recipe in ON_RECIPES)
