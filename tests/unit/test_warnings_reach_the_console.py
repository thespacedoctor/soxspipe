"""Python warnings print on the terminal once, and are written to ``soxspipe.log`` once (DY-2284).

`logging.captureWarnings` sends each warning to the ``py.warnings`` logger, which propagates to the root logger and
no further. The terminal handler of the ``logging settings`` block sits on the ``fundamentals.logs`` logger, not on
root, so `route_warnings_to_console` attaches it to ``py.warnings`` as well.
"""

from __future__ import annotations

import io
import logging
import warnings
from collections.abc import Iterator
from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from soxspipe.commonutils import toolkit
from tests.unit.test_cli_dispatch import CliLogger, _run_cli

pytestmark = pytest.mark.unit

WARNING_TEXT = "invalid value encountered in sqrt"


class SetupLogs:
    """The terminal stream, the ``soxspipe.log`` path and the logger that mirror the shipped logging setup."""

    def __init__(self, log: logging.Logger, terminal: io.StringIO, filePath: Path) -> None:
        self.log = log
        self.terminal = terminal
        self.filePath = filePath

    def file_text(self) -> str:
        """Return what the root file handler has written."""
        return self.filePath.read_text()


@pytest.fixture(name="setupLogs")
def setup_logs(tmp_path: Path) -> Iterator[SetupLogs]:
    """Build the shipped layout: the file handler on root, the named terminal handler on a child logger."""
    rootLogger = logging.getLogger()
    warningsLogger = logging.getLogger(toolkit.WARNINGS_LOGGER_NAME)
    log = logging.getLogger("dy2284.fundamentals.logs")
    log.setLevel(logging.DEBUG)
    terminal = io.StringIO()
    terminalHandler = logging.StreamHandler(terminal)
    terminalHandler.set_name("consoleLog")
    log.addHandler(terminalHandler)
    filePath = tmp_path / "soxspipe.log"
    fileHandler = logging.FileHandler(filePath)
    rootLogger.addHandler(fileHandler)
    handlersBefore = list(warningsLogger.handlers)
    wasCapturing = warnings.showwarning is logging._showwarning
    logging.captureWarnings(False)
    logging.captureWarnings(True)
    yield SetupLogs(log, terminal, filePath)
    logging.captureWarnings(wasCapturing)
    rootLogger.removeHandler(fileHandler)
    fileHandler.close()
    log.removeHandler(terminalHandler)
    warningsLogger.handlers[:] = handlersBefore


def _raise_a_runtime_warning() -> None:
    """Raise numpy's `RuntimeWarning` with the warning filter switched to ``always``."""
    with warnings.catch_warnings():
        warnings.simplefilter("always")
        np.sqrt(np.array([-1.0]))


def test_a_warning_prints_once_on_the_terminal_and_is_written_once_to_the_log_file(setupLogs: SetupLogs) -> None:
    # ARRANGE
    toolkit.route_warnings_to_console(setupLogs.log)

    # ACT
    _raise_a_runtime_warning()

    # ASSERT
    assert setupLogs.terminal.getvalue().count(WARNING_TEXT) == 1
    assert setupLogs.file_text().count(WARNING_TEXT) == 1


def test_routing_the_console_twice_still_prints_a_warning_once(setupLogs: SetupLogs) -> None:
    # ARRANGE
    toolkit.route_warnings_to_console(setupLogs.log)
    toolkit.route_warnings_to_console(setupLogs.log)

    # ACT
    _raise_a_runtime_warning()

    # ASSERT
    assert setupLogs.terminal.getvalue().count(WARNING_TEXT) == 1


def test_a_logger_with_no_console_handler_is_left_alone(setupLogs: SetupLogs) -> None:
    # ARRANGE
    # A HANDLER WITH ANOTHER NAME, SO THE TEST FAILS IF ANY HANDLER IS ATTACHED, NOT ONLY `consoleLog`
    bareLog = logging.getLogger("dy2284.bare")
    otherHandler = logging.StreamHandler(io.StringIO())
    otherHandler.set_name("otherLog")
    bareLog.addHandler(otherHandler)
    warningsHandlersBefore = list(logging.getLogger(toolkit.WARNINGS_LOGGER_NAME).handlers)

    # ACT
    try:
        toolkit.route_warnings_to_console(bareLog)
        warningsHandlersAfter = list(logging.getLogger(toolkit.WARNINGS_LOGGER_NAME).handlers)
    finally:
        bareLog.removeHandler(otherHandler)

    # ASSERT
    assert warningsHandlersAfter == warningsHandlersBefore


def test_a_spawned_worker_logger_has_no_terminal_handler_so_the_parent_alone_prints_a_warning() -> None:
    # ARRANGE
    import queue
    from logging.handlers import QueueHandler

    reducer: Any = import_module("soxspipe.commonutils.reducer")
    rootLogger = logging.getLogger()
    rootLevel = rootLogger.level
    rootHandlersBefore = list(rootLogger.handlers)
    wasCapturing = warnings.showwarning is logging._showwarning
    bareWorkerLog = logging.getLogger("dy2284.worker")

    # ACT
    try:
        workerLog = reducer._worker_log(bareWorkerLog, {"logging settings": {}}, queue.Queue())
        loggersToCheck = [workerLog, rootLogger, logging.getLogger(toolkit.WARNINGS_LOGGER_NAME)]
        handlerNames = [handler.get_name() for logger in loggersToCheck for handler in logger.handlers]
        addedToRoot = [handler for handler in rootLogger.handlers if handler not in rootHandlersBefore]
    finally:
        for handler in list(rootLogger.handlers):
            if handler not in rootHandlersBefore:
                rootLogger.removeHandler(handler)
        rootLogger.setLevel(rootLevel)
        logging.captureWarnings(wasCapturing)

    # ASSERT
    assert "consoleLog" not in handlerNames
    assert len(addedToRoot) == 1 and isinstance(addedToRoot[0], QueueHandler)


def test_the_command_line_routes_warnings_to_the_console_after_the_logging_setup(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, setupLogs: SetupLogs
) -> None:
    # ARRANGE
    cliLog = CliLogger()
    cliLog.handlers = setupLogs.log.handlers

    # ACT
    _run_cli(monkeypatch, tmp_path, ["soxspipe", "session", "list"], logger=cliLog)
    _raise_a_runtime_warning()

    # ASSERT
    assert setupLogs.terminal.getvalue().count(WARNING_TEXT) == 1
    assert setupLogs.file_text().count(WARNING_TEXT) == 1


def test_a_warning_queued_by_a_bulk_worker_prints_once_and_is_written_once_in_the_parent(setupLogs: SetupLogs) -> None:
    # ARRANGE
    toolkit.route_warnings_to_console(setupLogs.log)
    forwarder: Any = import_module("soxspipe.commonutils.reducer")._parent_log_forwarder()
    record = logging.makeLogRecord(
        {
            "name": toolkit.WARNINGS_LOGGER_NAME,
            "msg": f"RuntimeWarning: {WARNING_TEXT}",
            "levelno": logging.WARNING,
            "levelname": "WARNING",
        }
    )

    # ACT
    forwarder.handle(record)

    # ASSERT
    assert setupLogs.terminal.getvalue().count(WARNING_TEXT) == 1
    assert setupLogs.file_text().count(WARNING_TEXT) == 1
