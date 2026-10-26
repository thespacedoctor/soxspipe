"""Tests for the bulk reduction preconditions and the worker logger that feeds the parent process."""

from __future__ import annotations

import logging
import queue
import sqlite3
from collections.abc import Iterator
from importlib import import_module
from logging.handlers import QueueHandler
from pathlib import Path
from typing import Any

import multiprocess
import pytest

from tests.factories import logging_settings

pytestmark = pytest.mark.unit


def _empty_product_frames_connection() -> sqlite3.Connection:
    """Return an in-memory database with the table `run_recipe_bulk` updates."""
    connection = sqlite3.connect(":memory:")
    connection.execute("create table product_frames (sof text, status_night1 text, error_message text)")
    return connection


def _run_bulk(tmp_path: Path, log: Any, settings: dict[str, Any] | None = None) -> None:
    """Call `run_recipe_bulk` for one unsupported-recipe SOF."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    reducerModule.run_recipe_bulk(
        log=log,
        recipe="unsupported",
        sofList=[str(tmp_path / "sof" / "first.sof")],
        commandList=["soxspipe unsupported first.sof"],
        settings=logging_settings(tmp_path / "x.log") if settings is None else settings,
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        conn=_empty_product_frames_connection(),
        sessionId="night1",
    )


def test_run_recipe_bulk_raises_before_starting_the_pool_when_logging_settings_are_missing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """A spawned worker takes its log level from `logging settings`, so the parent checks for them first."""
    # ARRANGE
    poolWasStarted = False

    def stub_pool(**kwargs: Any) -> list[dict[str, Any]]:
        nonlocal poolWasStarted
        poolWasStarted = True
        return []

    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "_run_in_spawn_pool", stub_pool)

    # ACT
    with pytest.raises(ValueError, match="logging settings"):
        _run_bulk(tmp_path, log, settings={"workspace-root-dir": str(tmp_path)})

    # ASSERT
    assert poolWasStarted is False


def test_run_recipe_bulk_leaves_the_global_start_method_alone(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """The pool takes its own spawn context, so no process-wide start method is switched."""
    # ARRANGE
    multiprocess.set_start_method("fork", force=True)
    seenInsidePool: list[str | None] = []

    def stub_pool(**kwargs: Any) -> list[dict[str, Any]]:
        seenInsidePool.append(multiprocess.get_start_method(allow_none=True))
        return []

    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "_run_in_spawn_pool", stub_pool)

    # ACT
    _run_bulk(tmp_path, log)

    # ASSERT
    assert seenInsidePool == ["fork"]
    assert multiprocess.get_start_method(allow_none=True) == "fork"


@pytest.fixture
def added_log_levels(monkeypatch: pytest.MonkeyPatch) -> Iterator[list[tuple[str, int]]]:
    """Keep `_worker_log` from changing process-wide logging state, and restore the root logger afterwards.

    Yield the list of ``(name, number)`` arguments that `addLoggingLevel` was called with.
    """
    rootLog = logging.getLogger()
    handlersBefore, levelBefore = list(rootLog.handlers), rootLog.level
    addedLevels: list[tuple[str, int]] = []
    monkeypatch.setattr(import_module("fundamentals.logs"), "addLoggingLevel", lambda *args: addedLevels.append(args))
    monkeypatch.setattr(logging, "captureWarnings", lambda capture: None)
    yield addedLevels
    rootLog.handlers[:] = handlersBefore
    rootLog.setLevel(levelBefore)


def test_worker_log_sends_a_bare_loggers_records_to_the_queue_and_adds_no_file_handler(
    tmp_path: Path,
    added_log_levels: list[tuple[str, int]],
) -> None:
    """A bare logger gets the print level and a queue handler at the settings' root level, and no file handler."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    logQueue: queue.Queue[logging.LogRecord] = queue.Queue()
    handlersBefore = list(logging.getLogger().handlers)

    # ACT
    workerLog = reducerModule._worker_log(logging.Logger("dy1534-bare"), logging_settings(tmp_path / "x.log"), logQueue)
    logging.getLogger("dy1534.child").info("queued-marker")
    logging.getLogger("dy1534.child").debug("below-root-level")

    # ASSERT
    assert workerLog.name == "fundamentals.logs"
    assert added_log_levels == [("PRINT", logging.INFO + 5)]
    assert logQueue.get_nowait().getMessage() == "queued-marker"
    assert logQueue.empty()
    addedHandlers = [handler for handler in logging.getLogger().handlers if handler not in handlersBefore]
    assert [type(handler) for handler in addedHandlers] == [QueueHandler]


def test_worker_log_returns_a_working_logger_unchanged(
    tmp_path: Path,
    added_log_levels: list[tuple[str, int]],
) -> None:
    """A logger that has handlers and a `print` method needs no queue handler."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    workingLog = logging.Logger("dy1534-working")
    workingLog.addHandler(logging.NullHandler())
    workingLog.print = lambda message: None  # type: ignore[attr-defined]
    handlersBefore = list(logging.getLogger().handlers)

    # ACT
    workerLog = reducerModule._worker_log(workingLog, logging_settings(tmp_path / "x.log"), queue.Queue())

    # ASSERT
    assert workerLog is workingLog
    assert logging.getLogger().handlers == handlersBefore
    assert added_log_levels == []


def test_worker_log_returns_a_bare_logger_unchanged_when_settings_have_no_logging_block(
    added_log_levels: list[tuple[str, int]],
) -> None:
    """Without a `logging settings` block there is nothing to set the root level from."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    bareLog = logging.Logger("dy1534-bare")
    handlersBefore = list(logging.getLogger().handlers)

    # ACT
    workerLog = reducerModule._worker_log(bareLog, {}, queue.Queue())

    # ASSERT
    assert workerLog is bareLog
    assert logging.getLogger().handlers == handlersBefore
    assert added_log_levels == []
