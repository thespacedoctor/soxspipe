"""Tests that bulk reduction selects the spawn start method around its pool and restores it."""

from __future__ import annotations

import copy
import logging
import sqlite3
from importlib import import_module
from pathlib import Path
from typing import Any

import multiprocess
import pytest
import yaml

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


def test_run_recipe_bulk_selects_the_spawn_start_method_before_building_the_pool(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """Forked workers crash in Apple Accelerate, so the pool must start spawned workers."""
    # ARRANGE
    multiprocess.set_start_method("fork", force=True)
    seenInsidePool: list[str | None] = []

    def stub_multiprocess(**kwargs: Any) -> list[dict[str, Any]]:
        seenInsidePool.append(multiprocess.get_start_method(allow_none=True))
        return []

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", stub_multiprocess)

    # ACT
    _run_bulk(tmp_path, log)

    # ASSERT
    assert seenInsidePool == ["spawn"]


@pytest.mark.parametrize("originalMethod", [None, "fork"])
def test_run_recipe_bulk_restores_the_original_start_method_after_the_pool_returns(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    originalMethod: str | None,
) -> None:
    """The process-wide start method must not stay switched after a bulk reduction."""
    # ARRANGE
    multiprocess.set_start_method(originalMethod, force=True)
    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", lambda **kwargs: [])

    # ACT
    _run_bulk(tmp_path, log)

    # ASSERT
    assert multiprocess.get_start_method(allow_none=True) == originalMethod


def test_run_recipe_bulk_restores_the_original_start_method_when_the_pool_raises(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """A crashing pool must not leave the process on the spawn start method."""
    # ARRANGE
    multiprocess.set_start_method("fork", force=True)

    def exploding_multiprocess(**kwargs: Any) -> list[dict[str, Any]]:
        raise RuntimeError("pool exploded")

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", exploding_multiprocess)

    # ACT
    with pytest.raises(RuntimeError, match="pool exploded"):
        _run_bulk(tmp_path, log)

    # ASSERT
    assert multiprocess.get_start_method(allow_none=True) == "fork"


def test_run_recipe_bulk_targets_the_multiprocess_package_not_stdlib_multiprocessing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """`fundamentals` builds its pool from `multiprocess`, so the stdlib context must stay untouched."""
    # ARRANGE
    import multiprocessing

    stdlibMethodBefore = multiprocessing.get_start_method(allow_none=True)
    seenInsidePool: list[str | None] = []

    def stub_multiprocess(**kwargs: Any) -> list[dict[str, Any]]:
        seenInsidePool.append(multiprocessing.get_start_method(allow_none=True))
        return []

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", stub_multiprocess)

    # ACT
    _run_bulk(tmp_path, log)

    # ASSERT
    assert seenInsidePool == [stdlibMethodBefore]


def test_run_recipe_bulk_raises_before_starting_the_pool_when_logging_settings_are_missing(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """A spawned worker rebuilds its logger from `logging settings`, so the parent checks for them first."""
    # ARRANGE
    poolWasStarted = False

    def stub_multiprocess(**kwargs: Any) -> list[dict[str, Any]]:
        nonlocal poolWasStarted
        poolWasStarted = True
        return []

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", stub_multiprocess)

    # ACT
    with pytest.raises(ValueError, match="logging settings"):
        _run_bulk(tmp_path, log, settings={"workspace-root-dir": str(tmp_path)})

    # ASSERT
    assert poolWasStarted is False


class _StubLoggingSetup:
    """Stand in for `setup_dryx_logging`, recording the yaml it is given and returning an unregistered logger."""

    def __init__(self, logPath: Path) -> None:
        self.recordedYaml: list[dict[str, Any]] = []
        self.logPath = logPath
        self.returnedLog = logging.Logger("dy1532-probe")
        self.fileHandler = logging.FileHandler(logPath, delay=True)
        self.streamHandler = logging.StreamHandler()
        for handler in (self.fileHandler, self.streamHandler):
            handler.setLevel(logging.DEBUG)
            self.returnedLog.addHandler(handler)

    def __call__(self, yaml_file: str) -> logging.Logger:
        # THE TEMPORARY YAML IS DELETED WHEN THE CALLER RETURNS, SO READ IT NOW
        with open(yaml_file) as settingsFile:
            self.recordedYaml.append(yaml.safe_load(settingsFile))
        return self.returnedLog

    def close(self) -> None:
        for handler in (self.fileHandler, self.streamHandler):
            handler.close()


@pytest.fixture
def stub_logging_setup(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Any:
    """Replace the shared logger setup with a recording stub, and close its handlers afterwards."""
    stub = _StubLoggingSetup(tmp_path / "stub.log")
    monkeypatch.setattr(import_module("fundamentals.logs"), "setup_dryx_logging", stub)
    yield stub
    stub.close()


def test_worker_log_rebuilds_a_bare_logger_from_a_copy_of_the_logging_settings(
    tmp_path: Path,
    stub_logging_setup: _StubLoggingSetup,
) -> None:
    """A bare logger is rebuilt by the shared setup, muted on the console, from settings that stay untouched."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    settings = logging_settings(tmp_path / "shared.log")
    settingsBefore = copy.deepcopy(settings)

    # ACT
    workerLog = reducerModule._worker_log(logging.Logger("dy1532-bare"), settings)

    # ASSERT
    recordedConfig = stub_logging_setup.recordedYaml[0]["logging settings"]
    assert len(stub_logging_setup.recordedYaml) == 1
    assert recordedConfig["disable_existing_loggers"] is False
    assert recordedConfig["handlers"]["file"]["mode"] == "a"
    assert settings == settingsBefore
    assert workerLog is stub_logging_setup.returnedLog
    assert stub_logging_setup.streamHandler.level == logging.WARNING
    assert stub_logging_setup.fileHandler.level == logging.DEBUG


def test_worker_log_rebuilds_from_settings_that_have_no_file_handler(
    tmp_path: Path,
    stub_logging_setup: _StubLoggingSetup,
) -> None:
    """Settings with only a console handler are passed on without inventing a file handler."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    settings = logging_settings(tmp_path / "shared.log")
    del settings["logging settings"]["handlers"]["file"]
    settings["logging settings"]["root"]["handlers"] = ["console"]

    # ACT
    workerLog = reducerModule._worker_log(logging.Logger("dy1532-bare"), settings)

    # ASSERT
    recordedConfig = stub_logging_setup.recordedYaml[0]["logging settings"]
    assert "file" not in recordedConfig["handlers"]
    assert workerLog is stub_logging_setup.returnedLog


def test_worker_log_returns_a_working_logger_unchanged(
    tmp_path: Path,
    stub_logging_setup: _StubLoggingSetup,
) -> None:
    """A logger that has handlers and a `print` method needs no rebuild."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    workingLog = logging.Logger("dy1532-working")
    workingLog.addHandler(logging.NullHandler())
    workingLog.print = lambda message: None  # type: ignore[attr-defined]

    # ACT
    workerLog = reducerModule._worker_log(workingLog, logging_settings(tmp_path / "shared.log"))

    # ASSERT
    assert workerLog is workingLog
    assert stub_logging_setup.recordedYaml == []


def test_worker_log_returns_a_bare_logger_unchanged_when_settings_have_no_logging_block(
    stub_logging_setup: _StubLoggingSetup,
) -> None:
    """Without a `logging settings` block there is nothing to rebuild from."""
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    bareLog = logging.Logger("dy1532-bare")

    # ACT
    workerLog = reducerModule._worker_log(bareLog, {})

    # ASSERT
    assert workerLog is bareLog
    assert stub_logging_setup.recordedYaml == []
