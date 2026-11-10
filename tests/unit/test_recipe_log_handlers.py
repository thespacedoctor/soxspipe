"""Unit tests for detaching a recipe's log handlers from the shared logger when the recipe ends (DY-1704).

``add_recipe_logger`` attaches the ``recipeLog`` and ``recipeErr`` file handlers to the logger every recipe shares. A
recipe that finishes (``clean_up``) or fails (the reducer's ``except`` block) must detach and close them, after its own
error is logged, so a later recipe's error never lands in its ``_ERROR.log``.
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
    sharedLogger.print = lambda *args, **kwargs: None
    sharedLogger.addHandler(logging.NullHandler())
    yield sharedLogger
    for handler in list(sharedLogger.handlers):
        handler.close()
        sharedLogger.removeHandler(handler)


def _recipe_handlers(logger: logging.Logger) -> list[logging.Handler]:
    return [handler for handler in logger.handlers if handler.get_name() in RECIPE_HANDLER_NAMES]


def _fail_after_attaching_recipe_log(productPath: Path, error: Exception) -> Any:
    """Return a ``run_recipe`` stand-in that attaches the recipe handlers as a recipe constructor does, then raises."""

    def fake_run_recipe(log: Any, *args: object, **kwargs: object) -> None:
        toolkit.add_recipe_logger(log, str(productPath))
        raise error

    return fake_run_recipe


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


def test_clean_up_writes_a_failure_message_to_the_recipe_error_log_before_detaching(
    tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    errorLog = tmp_path / "A" / "A_ERROR.log"
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = toolkit.add_recipe_logger(sharedLog, str(tmp_path / "A" / "A.fits"))
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    recipe.status = None
    recipe.conn = None
    recipe.settings = {"workspace-root-dir": str(tmp_path)}
    recipe.outDir = str(tmp_path / "tmp" / "intermediate")
    Path(recipe.outDir).mkdir(parents=True)

    # ACT
    recipe.clean_up(forceFail="the recipe's own failure")

    # ASSERT
    assert "the recipe's own failure" in errorLog.read_text()
    assert _recipe_handlers(sharedLog) == []


def test_a_bulk_worker_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    reducerModule = _reducer_module()
    errorLog = tmp_path / "A" / "A_ERROR.log"
    fakeRunRecipe = _fail_after_attaching_recipe_log(tmp_path / "A" / "A.fits", RuntimeError("worker boom"))
    monkeypatch.setattr(reducerModule, "run_recipe", fakeRunRecipe)

    # ACT
    result = reducerModule._bulk_worker(
        {"sof": "A.sof", "command": "soxspipe mbias A.sof"},
        log=sharedLog,
        recipe="mbias",
        settings={},
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        logQueue=None,
    )
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert result["status"] == "fail"
    assert "worker boom" in errorLog.read_text()
    assert "a later recipe failed" not in errorLog.read_text()
    assert _recipe_handlers(sharedLog) == []


def test_a_serial_reducer_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, sharedLog: Any
) -> None:
    # ARRANGE
    reducerModule = _reducer_module()
    errorLog = tmp_path / "A" / "A_ERROR.log"
    selected = [
        pd.DataFrame([{"recipe": "mbias", "sof": "A.sof", "command": "soxspipe mbias A.sof"}]),
        pd.DataFrame(columns=["recipe", "sof", "command"]),
    ]

    class StubOrganiser:
        def __init__(self, **kwargs: object) -> None:
            pass

        def session_refresh(self, **kwargs: object) -> None:
            pass

        def get_incomplete_sets_report(self) -> tuple[pd.DataFrame, pd.DataFrame]:
            return pd.DataFrame(), pd.DataFrame()

        def close(self) -> None:
            pass

    runner = reducerModule.reducer.__new__(reducerModule.reducer)
    runner.log = sharedLog
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
    monkeypatch.setattr(commonutils, "data_organiser", StubOrganiser)
    monkeypatch.setattr(
        reducerModule, "run_recipe", _fail_after_attaching_recipe_log(tmp_path / "A" / "A.fits", RuntimeError("boom"))
    )

    # ACT
    runner.reduce(batch=10, multiprocess=False)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert "boom" in errorLog.read_text()
    assert "a later recipe failed" not in errorLog.read_text()
    assert _recipe_handlers(sharedLog) == []
