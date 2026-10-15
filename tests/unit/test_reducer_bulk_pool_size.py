"""Tests that bulk reduction hands the right pool size to the worker pool (DY-1377)."""

from __future__ import annotations

import sqlite3
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from tests.factories import logging_settings

pytestmark = pytest.mark.unit

MEMORY_HEAVY_RECIPES = ["mflat", "stare_obj", "nod_std", "offset_obj"]
LIGHT_RECIPES = ["mbias", "mdark", "disp_solution", "order_centres", "spat_solution"]


def _empty_product_frames_connection() -> sqlite3.Connection:
    """Return an in-memory database with the table `run_recipe_bulk` updates."""
    connection = sqlite3.connect(":memory:")
    connection.execute("create table product_frames (sof text, status_night1 text, error_message text)")
    return connection


def _pool_size_for(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    recipe: str,
    **poolArguments: Any,
) -> object:
    """Run `run_recipe_bulk` against a recording pool and return the `poolSize` it was handed."""
    seenPoolSizes: list[object] = []

    def recording_multiprocess(**kwargs: Any) -> list[dict[str, Any]]:
        seenPoolSizes.append(kwargs["poolSize"])
        return []

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", recording_multiprocess)
    import_module("soxspipe.commonutils.reducer").run_recipe_bulk(
        log=log,
        recipe=recipe,
        sofList=[str(tmp_path / "sof" / "first.sof")],
        commandList=[f"soxspipe {recipe} first.sof"],
        settings=logging_settings(tmp_path / "x.log"),
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        conn=_empty_product_frames_connection(),
        sessionId="night1",
        **poolArguments,
    )
    assert len(seenPoolSizes) == 1
    return seenPoolSizes[0]


@pytest.mark.parametrize("recipe", MEMORY_HEAVY_RECIPES)
def test_memory_heavy_recipes_keep_a_pool_of_three_without_a_user_pool_size(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    recipe: str,
) -> None:
    """Without a user value, the memory-heavy recipes keep the built-in cap of 3."""
    # ACT
    poolSize = _pool_size_for(monkeypatch, tmp_path, log, recipe)

    # ASSERT
    assert poolSize == 3


@pytest.mark.parametrize("recipe", LIGHT_RECIPES)
def test_light_recipes_keep_the_cpu_default_without_a_user_pool_size(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    recipe: str,
) -> None:
    """Without a user value, other recipes leave the pool size to `fmultiprocess` (one per CPU)."""
    # ACT
    poolSize = _pool_size_for(monkeypatch, tmp_path, log, recipe)

    # ASSERT
    assert poolSize is False


@pytest.mark.parametrize("recipe", MEMORY_HEAVY_RECIPES + LIGHT_RECIPES)
def test_a_user_pool_size_overrides_the_built_in_choice_for_every_recipe(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    recipe: str,
) -> None:
    """A user-set pool size replaces both the cap of 3 and the CPU default."""
    # ACT
    poolSize = _pool_size_for(monkeypatch, tmp_path, log, recipe, poolSize=2)

    # ASSERT
    assert poolSize == 2


def test_the_progress_message_names_the_user_pool_size(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The user can see the number of worker processes in use."""
    # ACT
    _pool_size_for(monkeypatch, tmp_path, log, "mbias", poolSize=5)

    # ASSERT
    assert "with a pool size of 5" in capsys.readouterr().out


def test_the_progress_message_names_the_cpu_default_pool_size(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Without a user value, a light recipe's message names the one-per-CPU pool size."""
    # ARRANGE
    monkeypatch.setattr(import_module("psutil"), "cpu_count", lambda: 7)

    # ACT
    _pool_size_for(monkeypatch, tmp_path, log, "mbias")

    # ASSERT
    assert "with a pool size of 7" in capsys.readouterr().out
