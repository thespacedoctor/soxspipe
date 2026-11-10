"""Unit tests for detaching a recipe's log handlers from the shared logger when the recipe ends (DY-1704).

``add_recipe_logger`` attaches the ``recipeLog`` and ``recipeErr`` file handlers to the logger every recipe shares. The
reducer detaches and closes them once the recipe has finished or failed and the run trailer (``Recipe Command`` and
``Recipe Run Time``) is written, so the trailer reaches the recipe's own ``.log`` and a later recipe's error never
lands in its ``_ERROR.log``.
"""

from __future__ import annotations

import logging
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

RECIPE_HANDLER_NAMES = {"recipeLog", "recipeErr"}


def _reducer_module() -> Any:
    """Return the reducer module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.commonutils.reducer")


@pytest.fixture(name="sharedLog")
def shared_log() -> Any:
    """Return an unregistered logger with a ``print`` method, and release its handlers afterwards."""
    sharedLogger = logging.Logger("dy1704", level=logging.DEBUG)
    sharedLogger.print = lambda message, *args, **kwargs: sharedLogger.log(logging.INFO + 5, message)
    sharedLogger.addHandler(logging.NullHandler())
    yield sharedLogger
    for handler in list(sharedLogger.handlers):
        handler.close()
        sharedLogger.removeHandler(handler)


def _recipe_handlers(logger: logging.Logger) -> list[logging.Handler]:
    return [handler for handler in logger.handlers if handler.get_name() in RECIPE_HANDLER_NAMES]


def _finish_like_a_recipe(log: Any, productPath: Path, forceFail: bool | str = False) -> tuple[str, pd.DataFrame]:
    """Attach the recipe handlers, then run the real ``clean_up`` as a recipe that finishes does."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = toolkit.add_recipe_logger(log, str(productPath))
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    recipe.status = None
    recipe.conn = None
    recipe.settings = {"workspace-root-dir": str(productPath.parent)}
    recipe.outDir = str(productPath.parent / "tmp" / "intermediate")
    Path(recipe.outDir).mkdir(parents=True)
    recipe.clean_up(forceFail=forceFail)
    return str(productPath), recipe.qc


def _finishing_run_recipe(productPath: Path) -> Any:
    """Return a ``run_recipe`` stand-in for a recipe that finishes."""

    def fake_run_recipe(log: Any, *args: object, **kwargs: object) -> tuple[str, pd.DataFrame]:
        return _finish_like_a_recipe(log, productPath)

    return fake_run_recipe


def _failing_run_recipe(productPath: Path, error: Exception) -> Any:
    """Return a ``run_recipe`` stand-in that attaches the recipe handlers as a constructor does, then raises."""

    def fake_run_recipe(log: Any, *args: object, **kwargs: object) -> None:
        toolkit.add_recipe_logger(log, str(productPath))
        raise error

    return fake_run_recipe


class _StubOrganiser:
    def __init__(self, **kwargs: object) -> None:
        pass

    def session_refresh(self, **kwargs: object) -> None:
        pass

    def get_incomplete_sets_report(self) -> tuple[pd.DataFrame, pd.DataFrame]:
        return pd.DataFrame(), pd.DataFrame()

    def close(self) -> None:
        pass


def _serial_reducer(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: Any, fakeRunRecipe: Any) -> Any:
    """Return a reducer with one mbias SOF to reduce serially, with ``run_recipe`` replaced by ``fakeRunRecipe``."""
    reducerModule = _reducer_module()
    selected = [
        pd.DataFrame([{"recipe": "mbias", "sof": "A.sof", "command": "soxspipe mbias A.sof"}]),
        pd.DataFrame(columns=["recipe", "sof", "command"]),
    ]
    runner = reducerModule.reducer.__new__(reducerModule.reducer)
    runner.log = log
    runner.settings = {}
    runner.workspaceDirectory = str(tmp_path)
    runner.reductionTarget = "all"
    runner.recipeList = ["mbias"]
    runner.sessionId = "night1"
    runner.sessionPath = str(tmp_path / "sessions" / "night1")
    runner.sessionDB = str(tmp_path / "soxspipe.db")
    runner.overwrite = False
    runner.quitOnFail = False
    runner.daemon = True
    runner.verbose = False
    runner.pathToSettings = False
    monkeypatch.setattr(runner, "select_sof_files_to_process", lambda **kwargs: selected.pop(0))
    monkeypatch.setattr(commonutils, "data_organiser", _StubOrganiser)
    monkeypatch.setattr(reducerModule, "run_recipe", fakeRunRecipe)
    return runner


def _bulk_worker_result(tmp_path: Path, log: Any) -> dict[str, Any]:
    return _reducer_module()._bulk_worker(
        {"sof": "A.sof", "command": "soxspipe mbias A.sof"},
        log=log,
        recipe="mbias",
        settings={},
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        logQueue=None,
    )


def test_remove_recipe_logger_detaches_and_closes_only_the_recipe_handlers(tmp_path: Path, sharedLog: Any) -> None:
    # ARRANGE
    toolkit.add_recipe_logger(sharedLog, str(tmp_path / "A" / "A.fits"))
    recipeHandlers = _recipe_handlers(sharedLog)
    otherHandlers = [handler for handler in sharedLog.handlers if handler not in recipeHandlers]
    assert len(recipeHandlers) == 2

    # ACT
    toolkit.remove_recipe_logger(sharedLog)

    # ASSERT
    assert _recipe_handlers(sharedLog) == []
    assert sharedLog.handlers == otherHandlers
    assert all(handler.stream is None for handler in recipeHandlers if isinstance(handler, logging.FileHandler))


def test_remove_recipe_logger_does_nothing_when_no_recipe_handler_is_attached(sharedLog: Any) -> None:
    # ARRANGE
    handlersBefore = list(sharedLog.handlers)

    # ACT
    toolkit.remove_recipe_logger(sharedLog)

    # ASSERT
    assert sharedLog.handlers == handlersBefore


def test_clean_up_keeps_the_recipe_handlers_so_the_caller_can_write_the_run_trailer(
    tmp_path: Path, sharedLog: Any
) -> None:
    # ACT
    _finish_like_a_recipe(sharedLog, tmp_path / "A" / "A.fits", forceFail="the recipe's own failure")

    # ASSERT
    assert "the recipe's own failure" in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert len(_recipe_handlers(sharedLog)) == 2


def test_a_serial_reducer_success_writes_the_trailer_to_the_recipe_log_then_detaches(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, _finishing_run_recipe(tmp_path / "A" / "A.fits"))

    # ACT
    runner.reduce(batch=10, multiprocess=False)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert "Recipe Command: soxspipe mbias A.sof" in (tmp_path / "A" / "A.log").read_text()
    assert not (tmp_path / "A" / "A_ERROR.log").exists()
    assert _recipe_handlers(sharedLog) == []


def test_a_serial_reducer_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    fakeRunRecipe = _failing_run_recipe(tmp_path / "A" / "A.fits", RuntimeError("boom"))
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, fakeRunRecipe)

    # ACT
    runner.reduce(batch=10, multiprocess=False)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert "boom" in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert "a later recipe failed" not in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert "Recipe Command: soxspipe mbias A.sof" in (tmp_path / "A" / "A.log").read_text()
    assert _recipe_handlers(sharedLog) == []


def test_a_serial_reducer_detaches_when_a_recipe_raises_file_exists_error(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    fakeRunRecipe = _failing_run_recipe(tmp_path / "A" / "A.fits", FileExistsError("the product exists"))
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, fakeRunRecipe)

    # ACT
    runner.reduce(batch=10, multiprocess=False)

    # ASSERT
    assert _recipe_handlers(sharedLog) == []


def test_a_bulk_worker_success_detaches_the_recipe_handlers(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    monkeypatch.setattr(_reducer_module(), "run_recipe", _finishing_run_recipe(tmp_path / "A" / "A.fits"))

    # ACT
    result = _bulk_worker_result(tmp_path, sharedLog)

    # ASSERT
    assert result["status"] == "pass"
    assert _recipe_handlers(sharedLog) == []


def test_a_bulk_worker_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    errorLog = tmp_path / "A" / "A_ERROR.log"
    fakeRunRecipe = _failing_run_recipe(tmp_path / "A" / "A.fits", RuntimeError("worker boom"))
    monkeypatch.setattr(_reducer_module(), "run_recipe", fakeRunRecipe)

    # ACT
    result = _bulk_worker_result(tmp_path, sharedLog)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert result["status"] == "fail"
    assert "worker boom" in errorLog.read_text()
    assert "a later recipe failed" not in errorLog.read_text()
    assert _recipe_handlers(sharedLog) == []
