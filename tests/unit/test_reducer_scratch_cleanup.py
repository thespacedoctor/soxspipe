"""Unit tests for the scratch-directory clean-up of the serial reducer (DY-1700).

The end-to-end behaviour is pinned by ``test_a_recipe_that_raises_removes_its_scratch_directory`` and
``test_the_serial_reducer_removes_scratch_left_by_an_interrupted_recipe``. These tests cover what
those pins do not: the exception that reaches the caller, a recipe whose constructor raises, a
clean-up that itself fails, and the bounds of every deletion.
"""

from __future__ import annotations

import shutil
import sqlite3
import tempfile
import threading
from collections.abc import Callable, Iterator
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
from soxspipe.commonutils import scratch as scratch_store
from tests.factories import logging_settings
from tests.integration.test_recipe_from_frame_list import _bias_frames, _IsolatedOrganiser, _settings

pytestmark = pytest.mark.unit


def _reducer_module() -> Any:
    """Return the reducer module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.commonutils.reducer")


@pytest.fixture(autouse=True)
def isolated_registry(monkeypatch: pytest.MonkeyPatch) -> Iterator[None]:
    """Give each test its own scratch registry, so no test reads or leaves entries in the real one."""
    monkeypatch.setattr(scratch_store, "_REGISTRY", threading.local())
    yield


def _install_recipes(monkeypatch: pytest.MonkeyPatch, behaviours: dict[str, Callable[[], Any]]) -> None:
    """Replace ``soxs_mbias`` and ``soxs_mdark`` with recipes that register a scratch directory as base_recipe does.

    ``behaviours`` maps a SOF name to ``(scratch, action)``: the recipe registers ``scratch`` in its
    constructor, and ``produce_product`` calls ``action`` and returns its result.
    """
    recipesModule = import_module("soxspipe.recipes")

    class FakeRecipe:
        def __init__(self, **kwargs: object) -> None:
            self.scratch, self.action, self.constructorError = behaviours[str(kwargs["inputFrames"])]
            scratch_store.register_scratch_directory(str(self.scratch))
            if self.constructorError is not None:
                raise self.constructorError

        def produce_product(self) -> tuple[str, pd.DataFrame]:
            return self.action() or ("/workspace/product.fits", pd.DataFrame())

    monkeypatch.setattr(recipesModule, "soxs_mbias", FakeRecipe)
    monkeypatch.setattr(recipesModule, "soxs_mdark", FakeRecipe)


def _install_recipe(
    monkeypatch: pytest.MonkeyPatch, scratchPath: Path, error: Exception | None, *, raiseInConstructor: bool = False
) -> None:
    """Replace ``soxs_mbias`` with a recipe that registers ``scratchPath`` and raises ``error`` when it is given."""

    def action() -> None:
        if error is not None:
            raise error

    constructorError = error if raiseInConstructor else None
    _install_recipes(monkeypatch, {"input.sof": (scratchPath, action, constructorError)})


def _run_mbias(log: Any, workspace: Path, recipe: str = "mbias", sof: str = "input.sof") -> tuple[str, pd.DataFrame]:
    return _reducer_module().run_recipe(
        log=log,
        recipe=recipe,
        sof=sof,
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
    scratch_store.register_scratch_directory(str(sibling))
    mine = _scratch_with_frame(tmp_path, "tmpmine")
    _install_recipe(monkeypatch, mine, RuntimeError("recipe failed"))

    # ACT
    with pytest.raises(RuntimeError):
        _run_mbias(log, tmp_path)

    # ASSERT
    assert not mine.exists()
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
    """Settings that break the scratch-root lookup leaves the directory and the recipe's error intact."""
    # ARRANGE
    scratch = _scratch_with_frame(tmp_path)
    failure = ValueError("recipe failed")
    _install_recipe(monkeypatch, scratch, failure)

    # ACT
    with pytest.raises(ValueError) as caught:
        _reducer_module().run_recipe(log, "mbias", "input.sof", settings=None, overwrite=False, turnOffMP=True)

    # ASSERT
    assert caught.value is failure
    assert scratch.exists()
    assert any(level == "warning" and "scratch root" in message for level, message in log.messages)


def test_remove_scratch_tree_does_nothing_for_an_empty_path(log: Any) -> None:
    """A recipe with no scratch directory passes an empty value, and nothing happens."""
    # ACT
    scratch_store.remove_scratch_tree(log, "")

    # ASSERT
    assert log.messages == []


def test_scratch_root_expands_a_home_directory_workspace_path(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Recipes and the reducer read ``~/workspace`` the same way, so they agree on the scratch root."""
    # ARRANGE
    monkeypatch.setenv("HOME", str(tmp_path))

    # ACT
    root = scratch_store.scratch_root({"workspace-root-dir": "~/workspace"})

    # ASSERT
    assert root == f"{tmp_path}/workspace/tmp"


# ---------------------------------------------------------------------------
# THE KILL PATH: clear_workspace_scratch
# ---------------------------------------------------------------------------


def test_clear_workspace_scratch_removes_the_mkdtemp_directories_and_keeps_everything_else(
    log: Any, tmp_path: Path
) -> None:
    """Only the ``tmp*`` directories in ``<workspace>/tmp/`` go; ``tmp/`` itself and other entries stay."""
    # ARRANGE
    leftover = _scratch_with_frame(tmp_path)
    (tmp_path / "tmp" / "notes.txt").write_text("keep")
    (tmp_path / "tmp" / "tmpfile.txt").write_text("keep")
    (tmp_path / "tmp" / "keepme").mkdir()
    (tmp_path / "tmp" / "keepme" / "data.fits").write_bytes(b"keep")
    (tmp_path / "soxspipe.db").write_bytes(b"db")
    (tmp_path / "reduced").mkdir()
    (tmp_path / "reduced" / "product.fits").write_bytes(b"product")

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert not leftover.exists()
    assert (tmp_path / "tmp" / "notes.txt").read_text() == "keep"
    assert (tmp_path / "tmp" / "tmpfile.txt").read_text() == "keep"
    assert (tmp_path / "tmp" / "keepme" / "data.fits").read_bytes() == b"keep"
    assert (tmp_path / "soxspipe.db").read_bytes() == b"db"
    assert (tmp_path / "reduced" / "product.fits").read_bytes() == b"product"


def test_clear_workspace_scratch_does_nothing_when_there_is_no_tmp_directory(log: Any, tmp_path: Path) -> None:
    """A workspace that never ran a recipe has no ``tmp/``, and that is not an error."""
    # ARRANGE
    entriesBefore = sorted(tmp_path.iterdir())

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

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
    scratch_store.clear_workspace_scratch(log, "")

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
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

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
    scratch_store.clear_workspace_scratch(log, str(workspace))

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
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert scratch.exists()
    assert any(level == "warning" and "cannot remove" in message for level, message in log.messages)


def test_clear_workspace_scratch_does_not_follow_a_tmp_prefixed_link_out_of_the_scratch_directory(
    log: Any, tmp_path: Path
) -> None:
    """A ``tmp*`` link inside ``tmp/`` that leads outside it is refused and its target kept."""
    # ARRANGE
    (tmp_path / "tmp").mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (outside / "keep.txt").write_text("keep")
    (tmp_path / "tmp" / "tmplnk00001").symlink_to(outside, target_is_directory=True)

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert (outside / "keep.txt").read_text() == "keep"
    assert any(level == "warning" for level, _ in log.messages)


def test_clear_workspace_scratch_keeps_a_tmp_prefixed_directory_that_mkdtemp_did_not_make(
    log: Any, tmp_path: Path
) -> None:
    """``tmp_backup`` starts with ``tmp`` but is not a ``mkdtemp`` name, so it stays."""
    # ARRANGE
    backup = tmp_path / "tmp" / "tmp_backup"
    backup.mkdir(parents=True)
    (backup / "data.fits").write_bytes(b"keep")

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert (backup / "data.fits").read_bytes() == b"keep"


def test_clear_workspace_scratch_removes_a_directory_made_by_mkdtemp(log: Any, tmp_path: Path) -> None:
    """A name that the standard library really generates is swept."""
    # ARRANGE
    root = Path(scratch_store.scratch_root({"workspace-root-dir": str(tmp_path)}))
    root.mkdir()
    made = Path(tempfile.mkdtemp(dir=root))

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert not made.exists()
    assert root.is_dir()


# ---------------------------------------------------------------------------
# THREADS AND NESTED CALLS: THE REGISTRY IS PER THREAD AND PER CALL
# ---------------------------------------------------------------------------


def test_a_failing_run_recipe_in_one_thread_leaves_the_live_scratch_of_another_thread_alone(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The failing call takes its snapshot, then another thread registers scratch, then the failing call raises.

    A registry shared by both threads would show the live thread's directory as new to the failing call, which
    would remove it. A per-thread registry cannot.
    """
    # ARRANGE
    liveScratch = _scratch_with_frame(tmp_path, "tmplive")
    failedScratch = _scratch_with_frame(tmp_path, "tmpfailed")
    liveHasRegistered = threading.Event()
    releaseLive = threading.Event()
    liveSawItsScratch: list[bool] = []

    def live_action() -> None:
        liveHasRegistered.set()
        assert releaseLive.wait(timeout=30)
        liveSawItsScratch.append(liveScratch.exists())

    liveThread = threading.Thread(target=_run_mbias, args=(log, tmp_path, "mbias", "live.sof"))

    def failing_action() -> None:
        # THE FAILING CALL IS ALREADY INSIDE run_recipe, SO ITS SNAPSHOT PRECEDES THE LIVE THREAD'S REGISTRATION
        liveThread.start()
        assert liveHasRegistered.wait(timeout=30)
        raise RuntimeError("other thread failed")

    _install_recipes(
        monkeypatch,
        {
            "live.sof": (liveScratch, live_action, None),
            "failed.sof": (failedScratch, failing_action, None),
        },
    )

    # ACT
    try:
        with pytest.raises(RuntimeError, match="other thread failed"):
            _run_mbias(log, tmp_path, "mdark", "failed.sof")
        liveWasIntactAfterTheFailure = liveScratch.exists()
    finally:
        releaseLive.set()
        liveThread.join(timeout=30)

    # ASSERT
    assert not liveThread.is_alive()
    assert not failedScratch.exists()
    assert liveWasIntactAfterTheFailure
    assert liveSawItsScratch == [True]


def test_run_recipe_removes_a_failed_recipes_scratch_after_an_outer_recipe_forgot_its_own(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """An entry that leaves the registry before this call's entries must not hide them from the clean-up."""
    # ARRANGE
    outerScratch = _scratch_with_frame(tmp_path, "tmpouter")
    failedScratch = _scratch_with_frame(tmp_path, "tmpfailed")
    scratch_store.register_scratch_directory(str(outerScratch))

    def failing_action() -> None:
        # THE OUTER RECIPE FINISHES, AND ITS clean_up FORGETS ITS DIRECTORY, WHILE THIS CALL IS RUNNING
        scratch_store.unregister_scratch_directory(str(outerScratch))
        raise RuntimeError("recipe failed")

    _install_recipes(monkeypatch, {"input.sof": (failedScratch, failing_action, None)})

    # ACT
    with pytest.raises(RuntimeError, match="recipe failed"):
        _run_mbias(log, tmp_path)

    # ASSERT
    assert not failedScratch.exists()
    assert scratch_store.created_scratch_directories() == []


def test_a_nested_run_recipe_removes_only_its_own_scratch_and_the_outer_call_removes_the_rest(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The inner call's mark starts after the outer call's directory, so the inner clean-up cannot reach it."""
    # ARRANGE
    outerScratch = _scratch_with_frame(tmp_path, "tmpouter")
    innerScratch = _scratch_with_frame(tmp_path, "tmpinner")
    seenAfterInner: dict[str, bool] = {}

    def inner_action() -> None:
        raise ValueError("inner failed")

    def outer_action() -> None:
        with pytest.raises(ValueError, match="inner failed"):
            _run_mbias(log, tmp_path, "mdark", "inner.sof")
        seenAfterInner["inner"] = innerScratch.exists()
        seenAfterInner["outer"] = outerScratch.exists()
        raise RuntimeError("outer failed")

    _install_recipes(
        monkeypatch,
        {"outer.sof": (outerScratch, outer_action, None), "inner.sof": (innerScratch, inner_action, None)},
    )

    # ACT
    with pytest.raises(RuntimeError, match="outer failed"):
        _run_mbias(log, tmp_path, "mbias", "outer.sof")

    # ASSERT
    assert seenAfterInner == {"inner": False, "outer": True}
    assert not outerScratch.exists()
    assert scratch_store.created_scratch_directories() == []


def test_run_recipe_trims_the_registry_after_a_run(log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """A long serial reduction does not grow the registry."""
    # ARRANGE
    _install_recipe(monkeypatch, tmp_path / "tmp" / "tmpdone", None)

    # ACT
    _run_mbias(log, tmp_path)

    # ASSERT
    assert scratch_store.created_scratch_directories() == []


# ---------------------------------------------------------------------------
# THE SCRATCH ROOT AND THE REGISTRY
# ---------------------------------------------------------------------------


def test_scratch_root_is_absolute_and_does_not_follow_a_later_change_of_directory(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A relative workspace root is fixed when the root is read, so a recipe that changes directory cannot move it."""
    # ARRANGE
    monkeypatch.chdir(tmp_path)
    (tmp_path / "elsewhere").mkdir()

    # ACT
    root = scratch_store.scratch_root({"workspace-root-dir": "./workspace"})
    monkeypatch.chdir(tmp_path / "elsewhere")

    # ASSERT
    assert root == str(tmp_path / "workspace" / "tmp")


def test_each_thread_has_its_own_scratch_registry() -> None:
    """A directory registered in one thread is invisible to another."""
    # ARRANGE
    scratch_store.register_scratch_directory("/workspace/tmp/tmpmine")
    seenByOtherThread: list[list[str]] = []

    # ACT
    thread = threading.Thread(
        target=lambda: seenByOtherThread.append(list(scratch_store.created_scratch_directories()))
    )
    thread.start()
    thread.join(timeout=30)

    # ASSERT
    assert seenByOtherThread == [[]]
    assert scratch_store.created_scratch_directories() == ["/workspace/tmp/tmpmine"]


def test_unregistering_a_directory_that_was_never_registered_is_ignored() -> None:
    """``clean_up`` can run on a recipe whose directory was already forgotten."""
    # ACT
    scratch_store.unregister_scratch_directory("/workspace/tmp/tmpnever")

    # ASSERT
    assert scratch_store.created_scratch_directories() == []


# ---------------------------------------------------------------------------
# base_recipe.clean_up USES THE SAME GUARDED REMOVAL
# ---------------------------------------------------------------------------


def _finished_recipe(log: Any, workspace: Path, outDir: Path) -> Any:
    """Return a recipe object that is ready for ``clean_up`` and has no database connection."""
    from soxspipe.recipes import base_recipe

    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["FLUX MEDIAN"]})
    recipe.status = None
    recipe.conn = None
    recipe.settings = {"workspace-root-dir": str(workspace)}
    recipe.workspaceRootPath = str(workspace)
    recipe.outDir = str(outDir)
    return recipe


def test_clean_up_removes_its_scratch_directory_and_forgets_it(log: Any, tmp_path: Path) -> None:
    """A finished recipe leaves nothing in ``tmp/`` and nothing in the registry."""
    # ARRANGE
    outDir = _scratch_with_frame(tmp_path)
    scratch_store.register_scratch_directory(str(outDir))
    recipe = _finished_recipe(log, tmp_path, outDir)

    # ACT
    recipe.clean_up()

    # ASSERT
    assert not outDir.exists()
    assert scratch_store.created_scratch_directories() == []


def test_clean_up_judges_its_out_dir_against_the_scratch_root_made_at_creation(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A relative ``workspaceRootPath`` and a later directory change must not make ``clean_up`` refuse its scratch."""
    # ARRANGE
    workspace = tmp_path / "workspace"
    outDir = _scratch_with_frame(workspace)
    (tmp_path / "elsewhere").mkdir()
    monkeypatch.chdir(tmp_path / "elsewhere")
    recipe = _finished_recipe(log, workspace, outDir)
    recipe.workspaceRootPath = "./"

    # ACT
    recipe.clean_up()

    # ASSERT
    assert not outDir.exists()


def test_clean_up_does_not_remove_an_out_dir_outside_the_workspace_tmp(log: Any, tmp_path: Path) -> None:
    """``clean_up`` follows the same deletion policy as ``run_recipe``."""
    # ARRANGE
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    outside = tmp_path / "elsewhere"
    outside.mkdir()
    (outside / "keep.fits").write_bytes(b"keep")
    recipe = _finished_recipe(log, workspace, outside)

    # ACT
    recipe.clean_up()

    # ASSERT
    assert (outside / "keep.fits").read_bytes() == b"keep"
    assert any(level == "warning" and str(outside) in message for level, message in log.messages)


# ---------------------------------------------------------------------------
# run_recipe_bulk SWEEPS EVEN WHEN ITS RESULT PROCESSING RAISES
# ---------------------------------------------------------------------------


def test_run_recipe_bulk_clears_the_scratch_when_the_pool_raises(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A failure after the workers started must not leave their scratch behind."""
    # ARRANGE
    leftover = _scratch_with_frame(tmp_path)
    reducerModule = _reducer_module()

    def failing_pool(**_: Any) -> list[dict[str, Any]]:
        raise RuntimeError("pool failed")

    monkeypatch.setattr(reducerModule, "_run_in_spawn_pool", failing_pool)
    connection = sqlite3.connect(":memory:")

    # ACT
    with pytest.raises(RuntimeError, match="pool failed"):
        reducerModule.run_recipe_bulk(
            log=log,
            recipe="mbias",
            sofList=[str(tmp_path / "first.sof")],
            commandList=["soxspipe mbias first.sof"],
            settings=logging_settings(tmp_path / "x.log"),
            overwrite=False,
            workspaceDirectory=str(tmp_path),
            conn=connection,
            sessionId="night1",
        )

    # ASSERT
    assert not leftover.exists()


def test_clear_workspace_scratch_logs_and_continues_when_tmp_cannot_be_listed(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A ``tmp/`` that cannot be read stops the sweep with a warning, not an error."""
    # ARRANGE
    leftover = _scratch_with_frame(tmp_path)

    def refuse_listing(path: object) -> list[str]:
        raise PermissionError("cannot list")

    monkeypatch.setattr("os.listdir", refuse_listing)

    # ACT
    scratch_store.clear_workspace_scratch(log, str(tmp_path))

    # ASSERT
    assert leftover.exists()
    assert any(level == "warning" and "cannot list" in message for level, message in log.messages)


def test_absolute_path_returns_an_empty_path_unchanged() -> None:
    """An empty path is not a home-directory path, and reading its first character must not raise."""
    # ACT
    result = scratch_store.absolute_path("")

    # ASSERT
    assert result == ""
