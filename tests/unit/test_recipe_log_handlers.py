"""Unit tests for detaching a recipe's log handlers from the shared logger when the recipe ends (DY-1704).

``add_recipe_logger`` attaches the ``recipeLog`` and ``recipeErr`` file handlers to the logger every recipe shares. The
reducer detaches and closes them once the recipe has finished or failed and the run trailer (``Recipe Command`` and
``Recipe Run Time``) is written, so the trailer reaches the recipe's own ``.log`` and a later recipe's error never
lands in its ``_ERROR.log``.
"""

import logging
from importlib import import_module
from pathlib import Path

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

RECIPE_HANDLER_NAMES = {"recipeLog", "recipeErr"}


def _reducer_module():
    """Return the reducer module; the package re-exports a class of the same name over it."""
    return import_module("soxspipe.commonutils.reducer")


@pytest.fixture(name="sharedLog")
def shared_log():
    """Return an unregistered logger with a ``print`` method, and release its handlers afterwards."""
    sharedLogger = logging.Logger("dy1704", level=logging.DEBUG)
    sharedLogger.print = lambda message, *args, **kwargs: sharedLogger.log(logging.INFO + 5, message)
    sharedLogger.addHandler(logging.NullHandler())
    yield sharedLogger
    for handler in list(sharedLogger.handlers):
        handler.close()
        sharedLogger.removeHandler(handler)


@pytest.fixture(name="attachedHandlers")
def attached_handlers(monkeypatch):
    """Return the list of recipe handlers that ``add_recipe_logger`` attached during the test."""
    attached = []
    addRecipeLogger = toolkit.add_recipe_logger

    def recording_add_recipe_logger(log, productPath):
        result = addRecipeLogger(log, productPath)
        attached.extend(_recipe_handlers(result))
        return result

    monkeypatch.setattr(toolkit, "add_recipe_logger", recording_add_recipe_logger)
    return attached


def _recipe_handlers(logger):
    return [handler for handler in logger.handlers if handler.get_name() in RECIPE_HANDLER_NAMES]


def _assert_closed(handlers):
    assert handlers
    assert all(handler.stream is None for handler in handlers if isinstance(handler, logging.FileHandler))


def _finish_like_a_recipe(log, productPath, forceFail=False):
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


def _finishing_run_recipe(productPath):
    """Return a ``run_recipe`` stand-in for a recipe that finishes."""

    def fake_run_recipe(log, *args, **kwargs):
        return _finish_like_a_recipe(log, productPath)

    return fake_run_recipe


def _failing_run_recipe(productPath, error):
    """Return a ``run_recipe`` stand-in that attaches the recipe handlers as a constructor does, then raises."""

    def fake_run_recipe(log, *args, **kwargs):
        toolkit.add_recipe_logger(log, str(productPath))
        raise error

    return fake_run_recipe


class _StubOrganiser:
    def __init__(self, **kwargs):
        pass

    def session_refresh(self, **kwargs):
        pass

    def get_incomplete_sets_report(self):
        return pd.DataFrame(), pd.DataFrame()

    def close(self):
        pass


def _serial_reducer(monkeypatch, tmp_path, log, fakeRunRecipe):
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


def _bulk_worker_result(tmp_path, log):
    return _reducer_module()._bulk_worker(
        {"sof": "A.sof", "command": "soxspipe mbias A.sof"},
        log=log,
        recipe="mbias",
        settings={},
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        logQueue=None,
    )


def test_remove_recipe_logger_detaches_and_closes_only_the_recipe_handlers(tmp_path, sharedLog):
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
    _assert_closed(recipeHandlers)


def test_remove_recipe_logger_does_nothing_when_no_recipe_handler_is_attached(sharedLog):
    # ARRANGE
    handlersBefore = list(sharedLog.handlers)

    # ACT
    toolkit.remove_recipe_logger(sharedLog)

    # ASSERT
    assert sharedLog.handlers == handlersBefore


def test_clean_up_keeps_the_recipe_handlers_so_the_caller_can_write_the_run_trailer(tmp_path, sharedLog):
    # ACT
    _finish_like_a_recipe(sharedLog, tmp_path / "A" / "A.fits", forceFail="the recipe's own failure")

    # ASSERT
    assert "the recipe's own failure" in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert len(_recipe_handlers(sharedLog)) == 2


def test_a_serial_reducer_success_writes_the_trailer_to_the_recipe_log_then_detaches(
    monkeypatch, tmp_path, sharedLog, attachedHandlers
):
    # ARRANGE
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, _finishing_run_recipe(tmp_path / "A" / "A.fits"))

    # ACT
    runner.reduce(batch=10, multiprocess=False)
    sharedLog.error("a later recipe failed")

    # ASSERT
    recipeLogText = (tmp_path / "A" / "A.log").read_text()
    assert "Recipe Command: soxspipe mbias A.sof" in recipeLogText
    assert "Recipe Run Time:" in recipeLogText
    assert not (tmp_path / "A" / "A_ERROR.log").exists()
    assert _recipe_handlers(sharedLog) == []
    _assert_closed(attachedHandlers)


def test_a_serial_reducer_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch, tmp_path, sharedLog
):
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


def test_a_serial_reducer_failure_detaches_when_quit_on_fail_exits(monkeypatch, tmp_path, sharedLog):
    # ARRANGE
    fakeRunRecipe = _failing_run_recipe(tmp_path / "A" / "A.fits", RuntimeError("boom"))
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, fakeRunRecipe)
    runner.quitOnFail = True

    # ACT
    with pytest.raises(SystemExit):
        runner.reduce(batch=10, multiprocess=False)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert "boom" in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert "a later recipe failed" not in (tmp_path / "A" / "A_ERROR.log").read_text()
    assert _recipe_handlers(sharedLog) == []


def test_a_serial_reducer_detaches_when_a_recipe_raises_file_exists_error(monkeypatch, tmp_path, sharedLog):
    # ARRANGE
    fakeRunRecipe = _failing_run_recipe(tmp_path / "A" / "A.fits", FileExistsError("the product exists"))
    runner = _serial_reducer(monkeypatch, tmp_path, sharedLog, fakeRunRecipe)

    # ACT
    runner.reduce(batch=10, multiprocess=False)

    # ASSERT
    assert _recipe_handlers(sharedLog) == []


def test_a_bulk_worker_success_detaches_the_recipe_handlers(monkeypatch, tmp_path, sharedLog):
    # ARRANGE
    monkeypatch.setattr(_reducer_module(), "run_recipe", _finishing_run_recipe(tmp_path / "A" / "A.fits"))

    # ACT
    result = _bulk_worker_result(tmp_path, sharedLog)
    sharedLog.error("a later recipe failed")

    # ASSERT
    assert result["status"] == "pass"
    assert not (tmp_path / "A" / "A_ERROR.log").exists()
    assert _recipe_handlers(sharedLog) == []


def test_a_bulk_worker_failure_is_logged_to_the_recipe_error_log_then_the_handlers_detach(
    monkeypatch, tmp_path, sharedLog
):
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
