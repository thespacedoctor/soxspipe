"""Bulk reduction workers log through the parent process, so one handler writes and rotates the log (DY-1534)."""

from __future__ import annotations

import io
import logging
import multiprocessing
import os
import queue
import signal
import sqlite3
from collections import Counter
from collections.abc import Iterator
from importlib import import_module
from logging.handlers import QueueListener
from pathlib import Path
from typing import Any

import pytest

from tests.factories import logging_settings

pytestmark = pytest.mark.unit

WORKER_COUNT = 4
SOF_COUNT = 8
LINES_PER_SOF = 150
# SMALL ENOUGH THAT THE WORKERS ROLL THE LOG OVER DOZENS OF TIMES, LARGE ENOUGH FOR ONE RECORD TO FIT
TINY_MAX_BYTES = 4096
# NO BACKUP IS EVER DISCARDED, SO EVERY RECORD CAN BE COUNTED
KEEP_ALL_BACKUPS = 500
CRASH_COMMAND = "crash"
TRACEBACK_COMMAND = "traceback"
RECORD_PREFIX = "RECORD"
POOL_TIMEOUT_SECONDS = 120


def log_lines(inputDict: dict[str, Any], log: Any, settings: dict[str, Any], **kwargs: Any) -> dict[str, Any]:
    """Log a numbered block of lines through the worker logger, as a bulk worker does."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    workerLog = reducerModule._worker_log(log, settings, kwargs["logQueue"])
    if inputDict["command"] == TRACEBACK_COMMAND:
        try:
            os.stat("/dy1534/no/such/path")
        except FileNotFoundError:
            workerLog.error(f"{RECORD_PREFIX} {inputDict['sof']} failed", exc_info=True)
        return {
            "status": "previous-pass",
            "sof": inputDict["sof"],
            "productPath": None,
            "qcTable": None,
            "error_message": None,
        }
    if inputDict["command"] == CRASH_COMMAND:
        workerLog.info(f"{RECORD_PREFIX} {inputDict['sof']} before-crash")
        os.kill(os.getpid(), signal.SIGSEGV)
    for lineNumber in range(LINES_PER_SOF):
        workerLog.info(f"{RECORD_PREFIX} {inputDict['sof']} line {lineNumber:04d} {'x' * 40}")
    # A PREVIOUS-PASS RESULT CARRIES NO QC TABLE, SO THE PARENT HAS NOTHING TO WRITE TO THE DATABASE
    return {
        "status": "previous-pass",
        "sof": inputDict["sof"],
        "productPath": None,
        "qcTable": None,
        "error_message": None,
    }


@pytest.fixture(autouse=True)
def bounded_pool(monkeypatch: pytest.MonkeyPatch) -> None:
    """Give the pool a deadline well under the ten hours of production, so a hung pool fails the test."""
    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "POOL_TIMEOUT_SECONDS", POOL_TIMEOUT_SECONDS)


@pytest.fixture
def parent_log_file(tmp_path: Path) -> Iterator[Path]:
    """Attach the shipped rotating file handler to the root logger of this process and return its log path."""
    from fundamentals.logs import GroupWriteRotatingFileHandler

    logPath = tmp_path / "workspace.log"
    handler = GroupWriteRotatingFileHandler(logPath, maxBytes=TINY_MAX_BYTES, backupCount=KEEP_ALL_BACKUPS, delay=True)
    handler.setLevel(logging.INFO)
    handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logging.getLogger().addHandler(handler)
    yield logPath
    logging.getLogger().removeHandler(handler)
    handler.close()


@pytest.fixture
def parent_log(monkeypatch: pytest.MonkeyPatch) -> logging.Logger:
    """Return the logger a real run hands to `run_recipe_bulk`, with the `print` method the shipped setup adds."""
    monkeypatch.setattr(logging.Logger, "print", lambda self, message: None, raising=False)
    return logging.getLogger("fundamentals.logs")


def _run_bulk(tmp_path: Path, log: Any, logPath: Path, commands: list[str], poolSize: int = WORKER_COUNT) -> None:
    """Run `run_recipe_bulk` for one SOF per command, with the shared log at ``logPath``."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    settings = logging_settings(logPath)
    settings["logging settings"]["handlers"]["file"].update(maxBytes=TINY_MAX_BYTES, backupCount=KEEP_ALL_BACKUPS)
    connection = sqlite3.connect(":memory:")
    connection.execute("create table product_frames (sof text, status_night1 text, error_message text)")
    reducerModule.run_recipe_bulk(
        log=log,
        recipe="unsupported",
        sofList=[str(tmp_path / "sof" / f"{index}.sof") for index in range(len(commands))],
        commandList=commands,
        settings=settings,
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        conn=connection,
        sessionId="night1",
        poolSize=poolSize,
    )


def _logged_records(logPath: Path) -> list[str]:
    """Return every record line from the log file and all of its rolled-over backups."""
    lines: list[str] = []
    for path in [logPath, *sorted(logPath.parent.glob(f"{logPath.name}.*"))]:
        if path.exists():
            lines.extend(line for line in path.read_text().splitlines() if RECORD_PREFIX in line)
    return lines


def _log_text(logPath: Path) -> str:
    """Return the text of the log file and all of its rolled-over backups."""
    return "".join(
        path.read_text() for path in [logPath, *sorted(logPath.parent.glob(f"{logPath.name}.*"))] if path.exists()
    )


def _sof_path(tmp_path: Path, index: int) -> str:
    """Return the SOF path `_run_bulk` gives the SOF at ``index``."""
    return str(tmp_path / "sof" / f"{index}.sof")


def _expected_records(tmp_path: Path, indexes: range | list[int]) -> list[str]:
    """Return the record lines `log_lines` writes for the SOFs at ``indexes``."""
    return [
        f"INFO: {RECORD_PREFIX} {_sof_path(tmp_path, index)} line {lineNumber:04d} {'x' * 40}"
        for index in indexes
        for lineNumber in range(LINES_PER_SOF)
    ]


def test_workers_that_force_many_rollovers_lose_no_records(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    parent_log: logging.Logger,
    parent_log_file: Path,
) -> None:
    # ARRANGE
    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "_bulk_worker", log_lines)

    # ACT
    _run_bulk(tmp_path, parent_log, parent_log_file, ["go"] * SOF_COUNT)

    # ASSERT
    assert list(parent_log_file.parent.glob(f"{parent_log_file.name}.*")), "the log never rolled over"
    assert sorted(_logged_records(parent_log_file)) == sorted(_expected_records(tmp_path, range(SOF_COUNT)))


def test_a_worker_killed_mid_batch_keeps_each_record_of_the_other_sofs_exactly_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    parent_log: logging.Logger,
    parent_log_file: Path,
) -> None:
    # ARRANGE: ONE WORKER RUNS THE SOFS IN TURN, SO THE CRASH FALLS ON THE SECOND AND THE THIRD NEVER STARTED
    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "_bulk_worker", log_lines)
    expected = Counter(
        [
            *_expected_records(tmp_path, [0, 2]),
            f"INFO: {RECORD_PREFIX} {_sof_path(tmp_path, 1)} before-crash",
        ]
    )

    # ACT
    _run_bulk(tmp_path, parent_log, parent_log_file, ["go", CRASH_COMMAND, "go"], poolSize=1)

    # ASSERT
    assert Counter(_logged_records(parent_log_file)) == expected


def test_a_worker_record_with_a_traceback_reaches_the_parent_file_with_the_traceback_text_once(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    parent_log: logging.Logger,
    parent_log_file: Path,
) -> None:
    # ARRANGE: THE QUEUE HANDLER FOLDS THE TRACEBACK INTO THE MESSAGE AND DROPS exc_info
    monkeypatch.setattr(import_module("soxspipe.commonutils.reducer"), "_bulk_worker", log_lines)

    # ACT
    _run_bulk(tmp_path, parent_log, parent_log_file, [TRACEBACK_COMMAND])

    # ASSERT
    logText = _log_text(parent_log_file)
    assert logText.startswith(
        f"ERROR: {RECORD_PREFIX} {_sof_path(tmp_path, 0)} failed\nTraceback (most recent call last):"
    )
    assert logText.count("Traceback (most recent call last):") == 1
    assert logText.rstrip().endswith("No such file or directory: '/dy1534/no/such/path'")


def test_a_listener_that_fails_to_stop_is_logged_and_the_bulk_run_still_finishes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    parent_log_file: Path,
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    monkeypatch.setattr(reducerModule, "_run_in_spawn_pool", lambda **kwargs: [])

    def fail_to_stop(listener: QueueListener) -> None:
        raise OSError("the manager is gone")

    monkeypatch.setattr(QueueListener, "stop", fail_to_stop)

    # ACT
    _run_bulk(tmp_path, log, parent_log_file, ["go"])

    # ASSERT
    errors = [message for level, message in log.messages if level == "error"]
    assert len(errors) == 1
    assert "stop the worker log listener" in errors[0]
    assert "the manager is gone" in errors[0]
    assert log.messages[-1] == ("debug", "completed the ``run_recipe_bulk`` method")


def test_a_failure_to_stop_the_listener_does_not_replace_the_pool_exception(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    parent_log_file: Path,
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")

    def failing_pool(**kwargs: Any) -> list[dict[str, Any]]:
        raise RuntimeError("the pool failed")

    def fail_to_stop(listener: QueueListener) -> None:
        raise OSError("the manager is gone")

    monkeypatch.setattr(reducerModule, "_run_in_spawn_pool", failing_pool)
    monkeypatch.setattr(QueueListener, "stop", fail_to_stop)

    # ACT
    with pytest.raises(RuntimeError, match="the pool failed"):
        _run_bulk(tmp_path, log, parent_log_file, ["go"])

    # ASSERT
    assert any("stop the worker log listener" in message for level, message in log.messages if level == "error")


def test_a_manager_that_fails_to_shut_down_is_logged_and_not_raised(
    monkeypatch: pytest.MonkeyPatch,
    log: Any,
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")

    class DeadManager:
        def Queue(self) -> queue.Queue:  # noqa: N802
            return queue.Queue()

        def shutdown(self) -> None:
            raise OSError("the manager process was killed")

    class FakeContext:
        def Manager(self) -> DeadManager:  # noqa: N802
            return DeadManager()

    monkeypatch.setattr(multiprocessing, "get_context", lambda method: FakeContext())

    # ACT
    with reducerModule._worker_log_queue(log):
        pass

    # ASSERT
    errors = [message for level, message in log.messages if level == "error"]
    assert len(errors) == 1
    assert "shut down the worker log manager" in errors[0]
    assert "the manager process was killed" in errors[0]


def test_the_listener_writes_queued_records_then_stops_when_the_pool_raises(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    parent_log: logging.Logger,
    parent_log_file: Path,
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    queueUsedByTheWorkers: list[Any] = []
    stopCalls: list[QueueListener] = []
    originalStop = QueueListener.stop

    def spy_on_stop(listener: QueueListener) -> None:
        stopCalls.append(listener)
        originalStop(listener)

    def failing_pool(**kwargs: Any) -> list[dict[str, Any]]:
        queueUsedByTheWorkers.append(kwargs["logQueue"])
        lastWords = {"msg": f"{RECORD_PREFIX} last words", "levelno": logging.INFO, "levelname": "INFO"}
        kwargs["logQueue"].put(logging.makeLogRecord(lastWords))
        raise RuntimeError("the pool failed")

    monkeypatch.setattr(reducerModule, "_run_in_spawn_pool", failing_pool)
    monkeypatch.setattr(QueueListener, "stop", spy_on_stop)

    # ACT
    with pytest.raises(RuntimeError, match="the pool failed"):
        _run_bulk(tmp_path, parent_log, parent_log_file, ["go"])

    # ASSERT
    assert len(stopCalls) == 1
    assert _logged_records(parent_log_file) == [f"INFO: {RECORD_PREFIX} last words"]
    with pytest.raises((OSError, EOFError)):
        queueUsedByTheWorkers[0].put("the manager is gone")


@pytest.fixture
def forwarding_logger(tmp_path: Path) -> Iterator[tuple[logging.Logger, Path, io.StringIO]]:
    """Return a parent logger with a file handler and a terminal handler, plus the file path and terminal buffer."""
    parentLog = logging.getLogger("dy1534.forward")
    parentLog.propagate = False
    filePath = tmp_path / "forward.log"
    fileHandler = logging.FileHandler(filePath)
    terminal = io.StringIO()
    terminalHandler = logging.StreamHandler(terminal)
    for handler in (fileHandler, terminalHandler):
        handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        parentLog.addHandler(handler)
    yield parentLog, filePath, terminal
    for handler in (fileHandler, terminalHandler):
        parentLog.removeHandler(handler)
        handler.close()
    parentLog.propagate = True


def test_forwarder_gives_file_handlers_every_record_and_terminal_handlers_only_warnings(
    forwarding_logger: tuple[logging.Logger, Path, io.StringIO],
) -> None:
    # ARRANGE: THE RECORDS COME FROM A CHILD LOGGER, SO THEY MUST REACH THE HANDLERS OF ITS PARENT
    _, filePath, terminal = forwarding_logger
    forwarder = import_module("soxspipe.commonutils.reducer")._parent_log_forwarder()
    records = [
        logging.makeLogRecord({"name": "dy1534.forward.child", "msg": text, "levelno": level, "levelname": name})
        for text, level, name in (("quiet", logging.INFO, "INFO"), ("loud", logging.WARNING, "WARNING"))
    ]

    # ACT
    for record in records:
        forwarder.handle(record)

    # ASSERT
    assert filePath.read_text().splitlines() == ["INFO: quiet", "WARNING: loud"]
    assert terminal.getvalue().splitlines() == ["WARNING: loud"]


def test_forwarder_respects_the_level_of_each_parent_handler(
    forwarding_logger: tuple[logging.Logger, Path, io.StringIO],
) -> None:
    # ARRANGE
    parentLog, filePath, _ = forwarding_logger
    fileHandler = next(handler for handler in parentLog.handlers if isinstance(handler, logging.FileHandler))
    fileHandler.setLevel(logging.ERROR)
    forwarder = import_module("soxspipe.commonutils.reducer")._parent_log_forwarder()
    record = logging.makeLogRecord(
        {"name": "dy1534.forward", "msg": "warned", "levelno": logging.WARNING, "levelname": "WARNING"}
    )

    # ACT
    forwarder.handle(record)

    # ASSERT
    assert not filePath.exists() or filePath.read_text() == ""
