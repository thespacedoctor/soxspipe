"""Contracts for `settings_without_stacking`, the settings override behind the `--no-stack` option (DY-1723)."""

from __future__ import annotations

import copy

import pytest

from soxspipe.commonutils.stacking_override import settings_without_stacking

pytestmark = pytest.mark.unit

SOXS_BLOCKS = (
    "soxs-stare",
    "soxs-stare-std",
    "soxs-nod",
    "soxs-nod-std",
    "soxs-offset",
    "soxs-offset-std",
)


def _soxs_settings() -> dict[str, object]:
    settings: dict[str, object] = {"instrument": "soxs", "soxs-mbias": {"clip": 3}}
    for block in SOXS_BLOCKS:
        settings[block] = {"stack_before_extraction": True, "use_flat": True}
    return settings


def test_returns_a_new_dict_and_leaves_the_input_unchanged() -> None:
    # ARRANGE
    settings = _soxs_settings()
    snapshot = copy.deepcopy(settings)

    # ACT
    result = settings_without_stacking(settings)

    # ASSERT
    assert result is not settings
    assert settings == snapshot
    assert result["soxs-stare"] is not settings["soxs-stare"]


@pytest.mark.parametrize("block", SOXS_BLOCKS)
def test_turns_stacking_off_in_every_stare_nod_and_offset_block(block: str) -> None:
    # ARRANGE
    settings = _soxs_settings()

    # ACT
    result = settings_without_stacking(settings)

    # ASSERT
    assert result[block]["stack_before_extraction"] is False
    assert result[block]["use_flat"] is True


def test_marks_the_override_as_requested_on_the_command_line() -> None:
    # ACT
    result = settings_without_stacking(_soxs_settings())

    # ASSERT
    assert result["no-stack"] is True


def test_leaves_other_recipe_blocks_alone() -> None:
    # ACT
    result = settings_without_stacking(_soxs_settings())

    # ASSERT
    assert result["soxs-mbias"] == {"clip": 3}
    assert "stack_before_extraction" not in result["soxs-mbias"]


def test_handles_xshooter_settings_that_lack_the_std_and_offset_blocks() -> None:
    # ARRANGE
    settings = {
        "instrument": "xsh",
        "soxs-stare": {"stack_before_extraction": True},
        "soxs-nod": {"stack_before_extraction": True},
    }

    # ACT
    result = settings_without_stacking(settings)

    # ASSERT
    assert result["soxs-stare"]["stack_before_extraction"] is False
    assert result["soxs-nod"]["stack_before_extraction"] is False
    assert "soxs-offset" not in result
    assert "soxs-stare-std" not in result


def test_turns_stacking_off_in_arm_sub_blocks_that_carry_the_key() -> None:
    # ARRANGE
    settings = {
        "soxs-nod": {
            "stack_before_extraction": True,
            "nir": {"stack_before_extraction": True, "use_flat": False},
            "vis": {"use_flat": True},
        }
    }

    # ACT
    result = settings_without_stacking(settings)

    # ASSERT
    assert result["soxs-nod"]["nir"]["stack_before_extraction"] is False
    assert result["soxs-nod"]["nir"]["use_flat"] is False
    assert "stack_before_extraction" not in result["soxs-nod"]["vis"]


def test_adds_the_setting_to_a_block_that_does_not_carry_it() -> None:
    # ACT
    result = settings_without_stacking({"soxs-stare": {"use_flat": True}})

    # ASSERT
    assert result["soxs-stare"]["stack_before_extraction"] is False


@pytest.mark.parametrize(
    ("settings", "expected"),
    [({"no-stack": True}, True), ({"no-stack": False}, False), ({}, False)],
    ids=["set", "false", "absent"],
)
def test_a_recipe_sees_the_override_only_when_the_command_line_set_it(settings: dict, expected: bool) -> None:
    from soxspipe.recipes.base_recipe import base_recipe

    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.settings = settings

    # ACT
    requested = recipe._stacking_override_requested()

    # ASSERT
    assert requested is expected
