"""Unit tests for the scratch-directory clean-up of the single-recipe CLI commands (DY-2093).

``soxspipe mbias``, ``stare`` and the other single-recipe commands build a recipe and call
``produce_product`` themselves. A recipe that raises must not leave ``<workspace>/tmp/tmpXXXXXXXX``
behind, as it does not under ``reducer.run_recipe`` (DY-1700). Both paths share the one scope
``scratch.removing_new_scratch``, and these tests cover every command through the CLI.
"""

from __future__ import annotations

import shutil
import tempfile
import threading
from collections.abc import Iterator
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from soxspipe.commonutils import scratch as scratch_store
from tests.unit.test_cli_dispatch import _run_cli

pytestmark = pytest.mark.unit

# EVERY SINGLE-RECIPE COMMAND: ITS ARGV AND THE RECIPE CLASS THE CLI BUILDS FOR IT
RECIPE_COMMANDS = [
    (["soxspipe", "mbias", "frames"], "soxs_mbias"),
    (["soxspipe", "mdark", "frames"], "soxs_mdark"),
    (["soxspipe", "disp_solution", "frames"], "soxs_disp_solution"),
    (["soxspipe", "order_centres", "frames"], "soxs_order_centres"),
    (["soxspipe", "spat_solution", "frames"], "soxs_spatial_solution"),
    (["soxspipe", "mflat", "frames"], "soxs_mflat"),
    (["soxspipe", "stare", "frames"], "soxs_stare"),
    (["soxspipe", "stare_std", "frames"], "soxs_stare"),
    (["soxspipe", "nod", "frames"], "soxs_nod"),
    (["soxspipe", "nod_std", "frames"], "soxs_nod"),
    (["soxspipe", "offset", "frames"], "soxs_offset"),
]
COMMAND_IDS = [argv[1] for argv, _ in RECIPE_COMMANDS]


@pytest.fixture(autouse=True)
def isolated_registry(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Give each test its own scratch registry, so no test reads or leaves entries in the real one."""
    monkeypatch.setattr(scratch_store, "_REGISTRY", threading.local())
    yield


def _install_scratch_recipe(
    monkeypatch: pytest.MonkeyPatch, recipeAttribute: str, workspace: Path, failAt: str | None, error: Exception
) -> None:
    """Replace a recipe with one that makes and registers a scratch directory as ``base_recipe`` does.

    ``failAt`` is ``"constructor"`` or ``"produce_product"``, the stage at which the recipe raises ``error``.
    When ``failAt`` is *None* the recipe returns a product and removes its scratch, as ``clean_up`` does.
    """
    recipesModule = import_module("soxspipe.recipes")

    class ScratchRecipe:
        def __init__(self, **kwargs: object) -> None:
            self.outDir = tempfile.mkdtemp(dir=workspace / "tmp")
            scratch_store.register_scratch_directory(self.outDir)
            (Path(self.outDir) / "frame_pre.fits").write_bytes(b"\0" * 2880)
            if failAt == "constructor":
                raise error

        def produce_product(self) -> str:
            if failAt == "produce_product":
                raise error
            shutil.rmtree(self.outDir)
            scratch_store.unregister_scratch_directory(self.outDir)
            return "product.fits"

    monkeypatch.setattr(recipesModule, recipeAttribute, ScratchRecipe)


def _workspace_with_scratch_root(tmp_path: Path) -> Path:
    (tmp_path / "tmp").mkdir()
    return tmp_path


def _scratch_names(workspace: Path) -> set[str]:
    return {path.name for path in (workspace / "tmp").iterdir()}


@pytest.mark.parametrize("failAt", ["constructor", "produce_product"])
@pytest.mark.parametrize(("argv", "recipeAttribute"), RECIPE_COMMANDS, ids=COMMAND_IDS)
def test_a_failing_recipe_command_leaves_no_scratch_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, argv: list[str], recipeAttribute: str, failAt: str
) -> None:
    """A recipe that raises in its constructor or in ``produce_product`` leaves ``<workspace>/tmp/`` empty."""
    # ARRANGE
    workspace = _workspace_with_scratch_root(tmp_path)
    _install_scratch_recipe(monkeypatch, recipeAttribute, workspace, failAt, RuntimeError("recipe failed"))

    # ACT
    with pytest.raises(SystemExit):
        _run_cli(monkeypatch, tmp_path, argv)

    # ASSERT
    assert _scratch_names(workspace) == set()
    assert scratch_store.created_scratch_directories() == []


@pytest.mark.parametrize("failAt", ["constructor", "produce_product"])
@pytest.mark.parametrize(("argv", "recipeAttribute"), RECIPE_COMMANDS, ids=COMMAND_IDS)
def test_a_failing_recipe_command_raises_the_recipes_own_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, argv: list[str], recipeAttribute: str, failAt: str
) -> None:
    """The caller sees the recipe's exception, unchanged, as the cause of the CLI's exit."""
    # ARRANGE
    workspace = _workspace_with_scratch_root(tmp_path)
    error = ValueError("the recipe's own message")
    _install_scratch_recipe(monkeypatch, recipeAttribute, workspace, failAt, error)

    # ACT
    with pytest.raises(SystemExit) as exitInfo:
        _run_cli(monkeypatch, tmp_path, argv)

    # ASSERT
    assert exitInfo.value.code == 1
    assert exitInfo.value.__cause__ is error


@pytest.mark.parametrize(("argv", "recipeAttribute"), RECIPE_COMMANDS, ids=COMMAND_IDS)
def test_a_failing_recipe_command_keeps_scratch_directories_that_existed_before(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, argv: list[str], recipeAttribute: str
) -> None:
    """A scratch directory made before the command, and still registered, survives the failed command."""
    # ARRANGE
    workspace = _workspace_with_scratch_root(tmp_path)
    earlier = workspace / "tmp" / "tmpearlier1"
    earlier.mkdir()
    scratch_store.register_scratch_directory(str(earlier))
    _install_scratch_recipe(monkeypatch, recipeAttribute, workspace, "produce_product", RuntimeError("recipe failed"))

    # ACT
    with pytest.raises(SystemExit):
        _run_cli(monkeypatch, tmp_path, argv)

    # ASSERT
    assert _scratch_names(workspace) == {"tmpearlier1"}
    assert scratch_store.created_scratch_directories() == [str(earlier)]


@pytest.mark.parametrize(("argv", "recipeAttribute"), RECIPE_COMMANDS, ids=COMMAND_IDS)
def test_a_successful_recipe_command_runs_to_the_end(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, argv: list[str], recipeAttribute: str
) -> None:
    """A recipe that finishes is produced once and the command returns without raising."""
    # ARRANGE
    workspace = _workspace_with_scratch_root(tmp_path)
    _install_scratch_recipe(monkeypatch, recipeAttribute, workspace, None, RuntimeError("unused"))

    # ACT
    _run_cli(monkeypatch, tmp_path, argv)

    # ASSERT
    assert _scratch_names(workspace) == set()
    assert scratch_store.created_scratch_directories() == []


# ---------------------------------------------------------------------------
# THE SHARED SEAM: scratch.removing_new_scratch
# ---------------------------------------------------------------------------


def _new_scratch(workspace: Path, name: str) -> Path:
    path = workspace / "tmp" / name
    path.mkdir(parents=True)
    scratch_store.register_scratch_directory(str(path))
    return path


def test_removing_new_scratch_removes_directories_made_inside_the_scope_when_the_body_raises(
    log: Any, tmp_path: Path
) -> None:
    """The scope removes the directories registered inside it and re-raises the body's own exception."""
    # ARRANGE
    error = KeyError("body failed")
    settings = {"workspace-root-dir": str(tmp_path)}

    # ACT
    with pytest.raises(KeyError) as raised, scratch_store.removing_new_scratch(log, settings):
        inside = _new_scratch(tmp_path, "tmpinside1")
        raise error

    # ASSERT
    assert raised.value is error
    assert not inside.exists()
    assert scratch_store.created_scratch_directories() == []


def test_removing_new_scratch_keeps_directories_registered_before_the_scope(log: Any, tmp_path: Path) -> None:
    """A directory registered before the scope is neither removed nor forgotten."""
    # ARRANGE
    settings = {"workspace-root-dir": str(tmp_path)}
    before = _new_scratch(tmp_path, "tmpbefore1")

    # ACT
    with scratch_store.removing_new_scratch(log, settings):
        inside = _new_scratch(tmp_path, "tmpinside1")

    # ASSERT
    assert before.exists()
    assert not inside.exists()
    assert scratch_store.created_scratch_directories() == [str(before)]
