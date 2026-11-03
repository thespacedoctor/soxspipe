"""Unit tests for the scratch-directory clean-up of the serial reducer (DY-1700).

The end-to-end behaviour is pinned by ``test_a_recipe_that_raises_removes_its_scratch_directory`` and
``test_the_serial_reducer_removes_scratch_left_by_an_interrupted_recipe``. These tests cover what
those pins do not: the exception that reaches the caller, a recipe whose constructor raises, a
clean-up that itself fails, and the bounds of every deletion.
"""

from __future__ import annotations

import shutil
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
from tests.integration.test_recipe_from_frame_list import _bias_frames, _IsolatedOrganiser, _settings

pytestmark = pytest.mark.unit


def _reducer_module() -> Any:
    """Return the reducer module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.commonutils.reducer")


def _base_recipe_module() -> Any:
    """Return the base_recipe module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.recipes.base_recipe")


def _install_recipe(
    monkeypatch: pytest.MonkeyPatch, scratch: Path, error: Exception | None, *, raiseInConstructor: bool = False
) -> None:
    """Replace ``soxs_mbias`` with a recipe that registers ``scratch`` as base_recipe does."""
    recipesModule = import_module("soxspipe.recipes")
    baseRecipe = _base_recipe_module()

    class FakeRecipe:
        def __init__(self, **_: object) -> None:
            baseRecipe.CREATED_SCRATCH_DIRECTORIES.append(str(scratch))
            if raiseInConstructor and error is not None:
                raise error

        def produce_product(self) -> tuple[str, pd.DataFrame]:
            if error is not None:
                raise error
            return "/workspace/product.fits", pd.DataFrame()

    monkeypatch.setattr(recipesModule, "soxs_mbias", FakeRecipe)


def _run_mbias(log: Any, workspace: Path) -> tuple[str, pd.DataFrame]:
    return _reducer_module().run_recipe(
        log=log,
        recipe="mbias",
        sof="input.sof",
        settings={"workspace-root-dir": str(workspace)},
        overwrite=False,
        turnOffMP=True,
    )


def _scratch_with_frame(workspace: Path, name: str = "tmpabcdefgh") -> Path:
    scratch = workspace / "tmp" / name
    scratch.mkdir(parents=True)
    (scratch / "frame_pre.fits").write_bytes(b"\0" * 2880)
    return scratch


def _refuse(error: Exception) -> Any:
    def refuse(path: object, *args: object, **kwargs: object) -> None:
        raise error

    return refuse


# ---------------------------------------------------------------------------
# THE RAISE PATH: run_recipe
# ---------------------------------------------------------------------------


def test_run_recipe_reraises_the_same_exception_after_removing_the_scratch_directory(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The caller still sees the recipe's own exception object, and the scratch directory is gone."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = RuntimeError("recipe failed")
    _install_recipe(monkeypatch, scratch, failure)

    # ACT
    with pytest.raises(RuntimeError) as caught:
        _run_mbias(log, tmp_path)

    # ASSERT
    assert caught.value is failure
    assert not scratch.exists()


def test_run_recipe_removes_the_scratch_directory_when_the_recipe_constructor_raises(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A recipe that fails while it is being built never reaches ``produce_product``, and still cleans up."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = TypeError("Input frames are a mix of arms")
    _install_recipe(monkeypatch, scratch, failure, raiseInConstructor=True)

    # ACT
    with pytest.raises(TypeError) as caught:
        _run_mbias(log, tmp_path)

    # ASSERT
    assert caught.value is failure
    assert not scratch.exists()


def test_run_recipe_leaves_the_scratch_directories_of_earlier_recipes_alone(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Only the directories created during this call go; a sibling already in ``tmp/`` stays."""
    # ARRANGE
    sibling = _scratch_with_frame(tmp_path, "tmpsibling")
    baseRecipe = _base_recipe_module()
    baseRecipe.CREATED_SCRATCH_DIRECTORIES.append(str(sibling))
    scratch = _scratch_with_frame(tmp_path, "tmpmine")
    _install_recipe(monkeypatch, scratch, RuntimeError("recipe failed"))

    # ACT
    try:
        with pytest.raises(RuntimeError):
            _run_mbias(log, tmp_path)
    finally:
        baseRecipe.CREATED_SCRATCH_DIRECTORIES.remove(str(sibling))

    # ASSERT
    assert not scratch.exists()
    assert sibling.exists()


@pytest.mark.parametrize("error", [PermissionError("cannot remove"), RuntimeError("not an OSError")])
def test_run_recipe_raises_the_recipes_error_when_the_scratch_clean_up_fails(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, error: Exception
) -> None:
    """A clean-up error of any type is logged and never replaces the recipe's own exception."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = ValueError("recipe failed")
    _install_recipe(monkeypatch, scratch, failure)
    monkeypatch.setattr(shutil, "rmtree", _refuse(error))

    # ACT
    with pytest.raises(ValueError) as caught:
        _run_mbias(log, tmp_path)

    # ASSERT
    assert caught.value is failure
    assert any(level == "warning" and str(scratch) in message for level, message in log.messages)


def test_run_recipe_returns_the_product_when_the_recipe_already_removed_its_scratch_directory(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A recipe that finishes has run ``clean_up`` already, so the missing directory is not an error."""
    # ARRANGE
    (tmp_path / "tmp").mkdir()
    _install_recipe(monkeypatch, tmp_path / "tmp" / "already-removed", None)

    # ACT
    productPath, _ = _run_mbias(log, tmp_path)

    # ASSERT
    assert productPath == "/workspace/product.fits"
    assert not [item for item in log.messages if item[0] == "warning"]


def test_run_recipe_does_not_remove_a_scratch_directory_outside_the_workspace_tmp(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A registered directory whose parent is not ``<workspace>/tmp`` survives, and the refusal is logged."""
    # ARRANGE
    workspace = tmp_path / "workspace"
    (workspace / "tmp").mkdir(parents=True)
    outside = tmp_path / "elsewhere" / "tmpabcdefgh"
    outside.mkdir(parents=True)
    (outside / "keep.fits").write_bytes(b"keep")
    _install_recipe(monkeypatch, outside, RuntimeError("recipe failed"))

    # ACT
    with pytest.raises(RuntimeError):
        _run_mbias(log, workspace)

    # ASSERT
    assert (outside / "keep.fits").read_bytes() == b"keep"
    assert any(level == "warning" and str(outside) in message for level, message in log.messages)


def test_run_recipe_does_not_remove_the_target_of_a_scratch_link_that_points_outside_the_workspace(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A link inside ``tmp/`` is judged by where it leads, not where it sits."""
    # ARRANGE
    workspace = tmp_path / "workspace"
    (workspace / "tmp").mkdir(parents=True)
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep.txt").write_text("keep")
    link = workspace / "tmp" / "tmplinked"
    link.symlink_to(outside, target_is_directory=True)
    _install_recipe(monkeypatch, link, RuntimeError("recipe failed"))

    # ACT
    with pytest.raises(RuntimeError):
        _run_mbias(log, workspace)

    # ASSERT
    assert (outside / "keep.txt").read_text() == "keep"


def test_a_real_recipe_that_raises_in_its_constructor_removes_its_scratch_directory(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Real ``soxs_mbias`` writes its prepared frames in its constructor; a failure after that must not leak them."""
    # ARRANGE
    from soxspipe.recipes import base_recipe, soxs_mbias

    monkeypatch.setattr(commonutils, "data_organiser", _IsolatedOrganiser)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    framePaths = _bias_frames(tmp_path / "raw")
    seenAtFailure: dict[str, list[str]] = {}
    realPrepare = base_recipe.prepare_frames

    def prepare_then_fail(self: Any, *args: object, **kwargs: object) -> None:
        realPrepare(self, *args, **kwargs)
        seenAtFailure[self.outDir] = sorted(path.name for path in Path(self.outDir).iterdir())
        raise TypeError("synthetic failure after frame preparation")

    monkeypatch.setattr(soxs_mbias, "prepare_frames", prepare_then_fail)

    # ACT
    with pytest.raises(TypeError, match="synthetic failure after frame preparation"):
        _reducer_module().run_recipe(
            log, "mbias", framePaths, settings=_settings(workspace), overwrite=False, turnOffMP=True
        )

    # ASSERT: THE CONSTRUCTOR REALLY PREPARED FRAMES INTO ITS OWN SCRATCH DIRECTORY
    assert len(seenAtFailure) == 1
    (scratchDirectory, preparedFrames) = next(iter(seenAtFailure.items()))
    assert workspace in Path(scratchDirectory).parents
    assert any(name.endswith("_pre.fits") for name in preparedFrames)

    # ASSERT: AND REMOVED IT DESPITE THE FAILURE
    assert not Path(scratchDirectory).exists()


def test_run_recipe_logs_and_raises_the_recipes_error_when_the_scratch_root_cannot_be_found(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A settings value that breaks the scratch-root lookup leaves the directory and the recipe's error intact."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = ValueError("recipe failed")
    _install_recipe(monkeypatch, scratch, failure)

    # ACT
    with pytest.raises(ValueError) as caught:
        _reducer_module().run_recipe(
            log, "mbias", "input.sof", settings={"workspace-root-dir": ""}, overwrite=False, turnOffMP=True
        )

    # ASSERT
    assert caught.value is failure
    assert scratch.exists()
    assert any(level == "warning" and "scratch root" in message for level, message in log.messages)


def test_remove_scratch_tree_does_nothing_for_an_empty_path(log: Any) -> None:
    """A recipe with no scratch directory passes an empty value, and nothing happens."""
    # ACT
    _reducer_module().remove_scratch_tree(log, "")

    # ASSERT
    assert log.messages == []


def test_scratch_root_expands_a_home_directory_workspace_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Recipes and the reducer read ``~/workspace`` the same way, so they agree on the scratch root."""
    # ARRANGE
    monkeypatch.setenv("HOME", str(tmp_path))

    # ACT
    root = _base_recipe_module().scratch_root({"workspace-root-dir": "~/workspace"})

    # ASSERT
    assert root == f"{tmp_path}/workspace/tmp"


# ---------------------------------------------------------------------------
# THE KILL PATH: clear_workspace_scratch
# ---------------------------------------------------------------------------


def test_clear_workspace_scratch_removes_the_tmp_directory_and_keeps_the_rest_of_the_workspace(
    log: Any, tmp_path: Path
) -> None:
    """Only ``<workspace>/tmp/`` goes; its siblings stay."""
    # ARRANGE
    _scratch_with_frame(tmp_path)
    (tmp_path / "soxspipe.db").write_bytes(b"db")
    (tmp_path / "reduced").mkdir()
    (tmp_path / "reduced" / "product.fits").write_bytes(b"product")

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert not (tmp_path / "tmp").exists()
    assert (tmp_path / "soxspipe.db").read_bytes() == b"db"
    assert (tmp_path / "reduced" / "product.fits").read_bytes() == b"product"


def test_clear_workspace_scratch_does_nothing_when_there_is_no_tmp_directory(log: Any, tmp_path: Path) -> None:
    """A workspace that never ran a recipe has no ``tmp/``, and that is not an error."""
    # ARRANGE
    entriesBefore = sorted(tmp_path.iterdir())

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert sorted(tmp_path.iterdir()) == entriesBefore
    assert not [item for item in log.messages if item[0] == "warning"]


def test_clear_workspace_scratch_skips_an_empty_workspace_path(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An empty path would otherwise name ``tmp`` under the current directory; nothing is removed."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    monkeypatch.chdir(tmp_path)

    # ACT
    _reducer_module().clear_workspace_scratch(log, "")

    # ASSERT
    assert scratch.exists()
    assert any(level == "warning" for level, _ in log.messages)


def test_clear_workspace_scratch_skips_the_home_directory(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A workspace path that is the home directory is almost certainly a mistake, and its ``tmp/`` stays."""
    # ARRANGE
    monkeypatch.setenv("HOME", str(tmp_path))
    scratch = _scratch_with_frame(tmp_path)

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert scratch.exists()
    assert any(level == "warning" for level, _ in log.messages)


def test_clear_workspace_scratch_does_not_follow_a_tmp_symlink_out_of_the_workspace(log: Any, tmp_path: Path) -> None:
    """A ``tmp`` link that points outside the workspace is left alone, and the refusal is logged."""
    # ARRANGE
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep.txt").write_text("keep")
    (workspace / "tmp").symlink_to(outside, target_is_directory=True)

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(workspace))

    # ASSERT
    assert (workspace / "tmp").is_symlink()
    assert (outside / "keep.txt").read_text() == "keep"
    assert any(level == "warning" for level, _ in log.messages)


def test_clear_workspace_scratch_logs_and_continues_when_removal_fails(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A sweep that cannot delete must not stop the reduction that called it."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    monkeypatch.setattr(shutil, "rmtree", _refuse(PermissionError("cannot remove")))

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert scratch.exists()
    assert any(level == "warning" and "cannot remove" in message for level, message in log.messages)
