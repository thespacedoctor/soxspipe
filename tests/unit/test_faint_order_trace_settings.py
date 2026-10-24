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


BOOLEAN_KEY = "faint-order-trace"
NUMBER_KEYS = (
    "faint-order-min-detection-fraction",
    "faint-order-max-scatter-px",
    "faint-order-max-anchor-scatter-px",
    "faint-order-min-bin-snr",
    "same-source-max-join-arcsec",
)
INTEGER_KEY = "faint-order-bin-count"


@pytest.mark.parametrize("value", ["False", "true", 0, 1, 1.0, [True]])
def test_from_block_rejects_a_switch_that_is_not_a_real_bool_and_names_the_key(value: object) -> None:
    # ACT / ASSERT: A QUOTED "False" IS A NON-EMPTY STRING AND WOULD TURN THE TRACE ON IF IT WERE CAST
    with pytest.raises(ValueError, match=BOOLEAN_KEY):
        faint_trace_settings.from_block({BOOLEAN_KEY: value})


@pytest.mark.parametrize("key", NUMBER_KEYS)
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), 0, 0.0, -1.5, "1.5", True, [1.0]])
def test_from_block_rejects_a_number_that_is_not_finite_and_positive_and_names_the_key(key: str, value: object) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match=key):
        faint_trace_settings.from_block({key: value})


@pytest.mark.parametrize("value", [0, -3, 9.5, 15.0, "15", True, float("nan")])
def test_from_block_rejects_a_bin_count_that_is_not_a_positive_integer_and_names_the_key(value: object) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match=INTEGER_KEY):
        faint_trace_settings.from_block({INTEGER_KEY: value})


@pytest.mark.parametrize("key", [BOOLEAN_KEY, *NUMBER_KEYS, INTEGER_KEY])
def test_from_block_uses_the_default_for_a_key_set_to_none(key: str) -> None:
    # ACT
    settings = faint_trace_settings.from_block({key: None})

    # ASSERT
    assert settings == faint_trace_settings.from_block({})


def test_from_block_accepts_whole_numbers_for_the_float_settings() -> None:
    # ACT
    settings = faint_trace_settings.from_block({"faint-order-max-scatter-px": 2, "faint-order-bin-count": 8})

    # ASSERT
    assert settings.maxScatterPx == 2.0
    assert isinstance(settings.maxScatterPx, float)
    assert settings.binCount == 8
