"""Tests for the hand-off of a nodding SOF to the stare recipe (DY-2194)."""

from __future__ import annotations

from importlib import import_module
from typing import Any

import pytest

pytestmark = pytest.mark.unit

KEY = "reduce-as-stare"
SCIENCE_SOF = "/workspace/20251212T045200_VIS_1X1_1_NOD_OBJ_SLIT1_SOXS.sof"
STANDARD_SOF = "/workspace/20251212T045200_VIS_1X1_1_NOD_STD_SLIT1_SOXS.sof"


def _stare_mode() -> Any:
    return import_module("soxspipe.commonutils.stare_mode")


@pytest.mark.parametrize(
    ("inputFrames", "expected"),
    [
        (STANDARD_SOF, True),
        (SCIENCE_SOF, False),
        (["/workspace/a_STD_b.fits"], False),
        (False, False),
    ],
)
def test_is_standard_input_follows_the_std_marker_in_a_sof_name(inputFrames: object, expected: bool) -> None:
    assert _stare_mode().is_standard_input(inputFrames) is expected


@pytest.mark.parametrize(
    ("settings", "isStandard", "isRequested", "expected"),
    [
        ({}, False, False, False),
        ({}, True, False, False),
        ({"soxs-nod": {KEY: True}}, False, False, True),
        ({"soxs-nod": {KEY: False}}, False, False, False),
        ({"soxs-nod-std": {KEY: True}}, True, False, True),
        ({"soxs-nod": {KEY: True}}, True, False, False),
        ({"soxs-nod-std": {KEY: True}}, False, False, False),
        ({}, False, True, True),
        ({"soxs-nod": {KEY: False}}, False, True, True),
        ({"soxs-nod-std": {KEY: False}}, True, True, True),
        ({"soxs-nod": None}, False, False, False),
        ({"soxs-nod": {}}, False, False, False),
    ],
)
def test_reduces_nod_as_stare_lets_the_request_win_over_the_block_key(
    settings: dict, isStandard: bool, isRequested: bool, expected: bool
) -> None:
    assert _stare_mode().reduces_nod_as_stare(settings, isStandard, isRequested) is expected


def test_nod_marker_names_the_stare_mode_products() -> None:
    assert _stare_mode().NOD_STARE_MARKER == "_STAREMODE"
    assert _stare_mode().REDUCE_AS_STARE_KEY == KEY


class RecordingRecipe:
    """Stand in for a recipe constructor and keep every call."""

    calls: list[dict[str, object]]

    def __init__(self, **kwargs: object) -> None:
        self.calls.append(kwargs)


def _recorders(monkeypatch: pytest.MonkeyPatch) -> tuple[list[dict], list[dict]]:
    recipes = import_module("soxspipe.recipes")
    nodCalls: list[dict] = []
    stareCalls: list[dict] = []
    monkeypatch.setattr(recipes, "soxs_nod", type("nod", (RecordingRecipe,), {"calls": nodCalls}))
    monkeypatch.setattr(recipes, "soxs_stare", type("stare", (RecordingRecipe,), {"calls": stareCalls}))
    return nodCalls, stareCalls


RECIPE_KWARGS = {"verbose": True, "overwrite": False, "command": "soxspipe nod x.sof", "debug": False}


def test_nod_recipe_builds_a_nod_recipe_with_todays_arguments_by_default(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    nodCalls, stareCalls = _recorders(monkeypatch)

    # ACT
    recipe = _stare_mode().nod_recipe(log="log", settings={"a": 1}, inputFrames=SCIENCE_SOF, **RECIPE_KWARGS)

    # ASSERT
    assert type(recipe).__name__ == "nod"
    assert stareCalls == []
    assert nodCalls == [{"log": "log", "settings": {"a": 1}, "inputFrames": SCIENCE_SOF, **RECIPE_KWARGS}]


def test_nod_recipe_builds_a_stare_recipe_in_nod_mode_when_requested(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    nodCalls, stareCalls = _recorders(monkeypatch)

    # ACT
    recipe = _stare_mode().nod_recipe(
        log="log", settings={"a": 1}, inputFrames=SCIENCE_SOF, isRequested=True, **RECIPE_KWARGS
    )

    # ASSERT
    assert type(recipe).__name__ == "stare"
    assert nodCalls == []
    assert stareCalls == [
        {"log": "log", "settings": {"a": 1}, "inputFrames": SCIENCE_SOF, "reduceNodAsStare": True, **RECIPE_KWARGS}
    ]


def test_nod_recipe_follows_the_settings_key_without_a_request(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    nodCalls, stareCalls = _recorders(monkeypatch)
    settings = {"soxs-nod": {KEY: True}}

    # ACT
    _stare_mode().nod_recipe(log="log", settings=settings, inputFrames=SCIENCE_SOF, **RECIPE_KWARGS)

    # ASSERT
    assert nodCalls == []
    assert len(stareCalls) == 1


def test_nod_recipe_reads_the_standard_block_for_a_standard_sof(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    nodCalls, stareCalls = _recorders(monkeypatch)
    settings = {"soxs-nod": {KEY: True}, "soxs-nod-std": {KEY: False}}

    # ACT
    _stare_mode().nod_recipe(log="log", settings=settings, inputFrames=STANDARD_SOF, **RECIPE_KWARGS)

    # ASSERT
    assert stareCalls == []
    assert len(nodCalls) == 1


@pytest.mark.parametrize(
    ("settings", "inputFrames", "isRequested", "recipeName", "marker"),
    [
        ({}, SCIENCE_SOF, False, "soxs-nod", ""),
        ({}, STANDARD_SOF, False, "soxs-nod-std", ""),
        ({"soxs-nod": {KEY: True}}, SCIENCE_SOF, False, "soxs-stare", "_STAREMODE"),
        ({"soxs-nod-std": {KEY: True}}, STANDARD_SOF, False, "soxs-stare-std", "_STAREMODE"),
        ({}, SCIENCE_SOF, True, "soxs-stare", "_STAREMODE"),
    ],
)
def test_predicted_nod_product_path_follows_the_mode_the_recipe_will_run_in(
    monkeypatch: pytest.MonkeyPatch,
    settings: dict,
    inputFrames: str,
    isRequested: bool,
    recipeName: str,
    marker: str,
) -> None:
    # ARRANGE
    calls: list[tuple] = []
    toolkit = import_module("soxspipe.commonutils.toolkit")
    monkeypatch.setattr(
        toolkit, "predict_product_path", lambda *args: calls.append(args) or ("predicted.fits", "night")
    )

    # ACT
    path = _stare_mode().predicted_nod_product_path(settings, inputFrames, isRequested)

    # ASSERT
    assert path == "predicted.fits"
    assert calls == [(inputFrames, recipeName, marker)]
