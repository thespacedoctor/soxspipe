"""Unit tests for the scratch-directory clean-up of the serial reducer (DY-1700).

The end-to-end behaviour is pinned by ``test_a_recipe_that_raises_removes_its_scratch_directory`` and
``test_the_serial_reducer_removes_scratch_left_by_an_interrupted_recipe``. These tests cover what
those pins do not: the exception that reaches the caller, a clean-up that itself fails, and the
bounds of the workspace sweep.
"""

from __future__ import annotations

import shutil
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytestmark = pytest.mark.unit


def _reducer_module() -> Any:
    """Return the reducer module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.commonutils.reducer")


def _install_recipe(monkeypatch: pytest.MonkeyPatch, outDir: Path, error: Exception | None) -> None:
    """Replace ``soxs_mbias`` with a recipe that owns ``outDir`` and raises ``error`` when it is given."""
    recipesModule = import_module("soxspipe.recipes")

    class FakeRecipe:
        def __init__(self, **_: object) -> None:
            self.outDir = str(outDir)

        def produce_product(self) -> tuple[str, pd.DataFrame]:
            if error is not None:
                raise error
            return "/workspace/product.fits", pd.DataFrame()

    monkeypatch.setattr(recipesModule, "soxs_mbias", FakeRecipe)


def _run_mbias(log: Any) -> tuple[str, pd.DataFrame]:
    return _reducer_module().run_recipe(
        log=log, recipe="mbias", sof="input.sof", settings={}, overwrite=False, turnOffMP=True
    )


def _scratch_with_frame(workspace: Path) -> Path:
    scratch = workspace / "tmp" / "tmpabcdefgh"
    scratch.mkdir(parents=True)
    (scratch / "frame_pre.fits").write_bytes(b"\0" * 2880)
    return scratch


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
        _run_mbias(log)

    # ASSERT
    assert caught.value is failure
    assert not scratch.exists()


def test_run_recipe_raises_the_recipes_error_when_the_scratch_clean_up_fails(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A clean-up error is logged and never replaces the recipe's own exception."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = ValueError("recipe failed")
    _install_recipe(monkeypatch, scratch, failure)

    def refuse(path: object, *args: object, **kwargs: object) -> None:
        raise PermissionError("cannot remove")

    monkeypatch.setattr(shutil, "rmtree", refuse)

    # ACT
    with pytest.raises(ValueError) as caught:
        _run_mbias(log)

    # ASSERT
    assert caught.value is failure
    assert any(level == "warning" and str(scratch) in message for level, message in log.messages)


def test_run_recipe_returns_the_product_when_the_recipe_already_removed_its_scratch_directory(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A recipe that finishes has run ``clean_up`` already, so the missing directory is not an error."""
    # ARRANGE
    _install_recipe(monkeypatch, tmp_path / "tmp" / "already-removed", None)

    # ACT
    productPath, _ = _run_mbias(log)

    # ASSERT
    assert productPath == "/workspace/product.fits"
    assert not [item for item in log.messages if item[0] == "warning"]


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
    entriesBefore = sorted(tmp_path.iterdir())

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert sorted(tmp_path.iterdir()) == entriesBefore
    assert not [item for item in log.messages if item[0] == "warning"]


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
    assert (outside / "keep.txt").read_text() == "keep"
    assert any(level == "warning" for level, _ in log.messages)


def test_clear_workspace_scratch_logs_and_continues_when_removal_fails(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A sweep that cannot delete must not stop the reduction that called it."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)

    def refuse(path: object, *args: object, **kwargs: object) -> None:
        raise PermissionError("cannot remove")

    monkeypatch.setattr(shutil, "rmtree", refuse)

    # ACT
    _reducer_module().clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert scratch.exists()
    assert any(level == "warning" and "cannot remove" in message for level, message in log.messages)
