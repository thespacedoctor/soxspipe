"""Integration tests that run bulk reduction workers in real spawned processes."""

from __future__ import annotations

import logging
import subprocess
import sys
import textwrap
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest
import yaml

from tests.factories import logging_settings

pytestmark = pytest.mark.integration

# SPAWNING A WORKER RE-IMPORTS THE PACKAGE, SO THE REAL-SPAWN TESTS ARE GIVEN HEADROOM
SPAWN_TIMEOUT_SECONDS = 120

# CHILD PROCESSES MUST START INSIDE THE REPOSITORY: FROM ELSEWHERE COVERAGE FINDS NO pyproject.toml,
# RECORDS STATEMENT-ONLY DATA, AND THE COMBINE STEP FAILS AFTER EVERY TEST HAS PASSED
REPOSITORY_ROOT = Path(__file__).resolve().parents[2]

SETUP_SCRIPT = textwrap.dedent(
    """
    import logging, sys, yaml
    from importlib import import_module

    reducer = import_module("soxspipe.commonutils.reducer")
    settings = yaml.safe_load(open(sys.argv[1]))
    """
)


def _run_real_pool(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    logPath: Path,
    commands: list[str],
) -> list[dict[str, Any]]:
    """Run the exact worker payload of `run_recipe_bulk` through a real spawned pool of one worker."""
    import sqlite3

    reducerModule = import_module("soxspipe.commonutils.reducer")
    realPool = reducerModule._run_in_spawn_pool
    results: list[dict[str, Any]] = []

    def run_real_pool(function: Any, inputDicts: list[dict[str, str]], poolSize: int, log: Any, **kwargs: Any) -> list:
        results.extend(realPool(function, inputDicts, 1, log, **kwargs))
        return []

    monkeypatch.chdir(REPOSITORY_ROOT)
    monkeypatch.setattr(reducerModule, "_run_in_spawn_pool", run_real_pool)
    # THE PARENT'S REAL LOGGER GAINS `print` FROM THE SETUP FUNCTION. THE CLASS IS PATCHED, NOT AN INSTANCE,
    # BECAUSE THE LOGGER IS PICKLED BY NAME AND A SPAWNED CHILD UNPICKLES A FRESH INSTANCE
    monkeypatch.setattr(logging.Logger, "print", lambda self, message: None, raising=False)
    connection = sqlite3.connect(":memory:")
    connection.execute("create table product_frames (sof text, status_night1 text, error_message text)")

    reducerModule.run_recipe_bulk(
        log=logging.getLogger("fundamentals.logs"),
        recipe="unsupported",
        sofList=[str(tmp_path / "sof" / f"{index}.sof") for index, _ in enumerate(commands)],
        commandList=commands,
        settings=logging_settings(logPath),
        overwrite=False,
        workspaceDirectory=str(tmp_path),
        conn=connection,
        sessionId="night1",
    )
    return results


def _run_in_fresh_python(tmp_path: Path, settings: dict[str, Any], body: str) -> subprocess.CompletedProcess:
    """Run a script body in a new interpreter, after loading `reducer` and `settings` from a yaml file."""
    settingsPath = tmp_path / "settings.yaml"
    settingsPath.write_text(yaml.safe_dump(settings))
    completed = subprocess.run(  # noqa: S603
        [sys.executable, "-c", SETUP_SCRIPT + textwrap.dedent(body), str(settingsPath)],
        capture_output=True,
        text=True,
        timeout=SPAWN_TIMEOUT_SECONDS,
        check=False,
        cwd=REPOSITORY_ROOT,
    )
    assert completed.returncode == 0, completed.stderr
    return completed


def test_spawned_pool_runs_the_bulk_worker_and_returns_its_result_with_a_working_log_file(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The worker payload survives a real spawn and the worker still logs to the settings log file.

    A `logging.Logger` pickles by name, so the unpickled copy in a spawned worker has no
    handlers and no `print` level. The worker must rebuild the logger from the settings.
    """
    # ARRANGE
    logPath = tmp_path / "worker.log"

    # ACT
    results = _run_real_pool(monkeypatch, tmp_path, logPath, ["only-command"])

    # ASSERT
    assert [result["status"] for result in results] == ["fail"]
    assert "Recipe failed with the following error" in logPath.read_text()


def test_reused_spawned_worker_logs_each_task_once_without_duplicate_handlers(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """One worker process handles both tasks, so a second rebuild would double every record."""
    # ARRANGE
    logPath = tmp_path / "reused.log"

    # ACT
    results = _run_real_pool(monkeypatch, tmp_path, logPath, ["first-command", "second-command"])

    # ASSERT
    logText = logPath.read_text()
    assert [result["status"] for result in results] == ["fail", "fail"]
    assert logText.count("Recipe Command: first-command") == 1
    assert logText.count("Recipe Command: second-command") == 1


def test_worker_logger_rebuilt_from_settings_has_print_level_and_writes_to_the_log_file(
    tmp_path: Path,
) -> None:
    """A fresh process given only a by-name logger gets handlers, the print level, and the log file."""
    # ARRANGE
    logPath = tmp_path / "rebuilt.log"

    # ACT
    _run_in_fresh_python(
        tmp_path,
        logging_settings(logPath),
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings)
        workerLog.print("printed from the worker")
        workerLog.error("errored from the worker")
        for handler in logging.getLogger().handlers:
            handler.flush()
        """,
    )

    # ASSERT
    assert "printed from the worker" in logPath.read_text()
    assert "errored from the worker" in logPath.read_text()


def test_worker_logger_rebuilt_from_settings_keeps_the_console_quiet_below_warning(
    tmp_path: Path,
) -> None:
    """Bulk workers are muted, so only warnings and errors may reach the terminal."""
    # ARRANGE
    settings = logging_settings(tmp_path / "quiet.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings)
        workerLog.print("printed from the worker")
        workerLog.error("errored from the worker")
        """,
    )

    # ASSERT
    assert "errored from the worker" in completed.stderr
    assert "printed from the worker" not in completed.stderr


def test_worker_logger_rebuilt_from_settings_appends_instead_of_truncating_the_shared_log(
    tmp_path: Path,
) -> None:
    """Several workers share one log file, so a rebuilt handler must not reopen it with `w+`."""
    # ARRANGE
    logPath = tmp_path / "shared.log"
    logPath.write_text("written by the parent\n")

    # ACT
    _run_in_fresh_python(
        tmp_path,
        logging_settings(logPath),
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings)
        workerLog.error("errored from the worker")
        for handler in logging.getLogger().handlers:
            handler.flush()
        """,
    )

    # ASSERT
    logText = logPath.read_text()
    assert "written by the parent" in logText
    assert "errored from the worker" in logText


def test_worker_logger_rebuilt_from_settings_leaves_other_loggers_enabled(
    tmp_path: Path,
) -> None:
    """Rebuilding must not silence loggers that already exist in the worker, such as astropy's."""
    # ARRANGE
    settings = logging_settings(tmp_path / "others.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        otherLog = logging.getLogger("astropy")
        reducer._worker_log(logging.getLogger("fundamentals.logs"), settings)
        print("OTHER DISABLED:", otherLog.disabled)
        """,
    )

    # ASSERT
    assert "OTHER DISABLED: False" in completed.stdout


def test_worker_logger_with_handlers_but_no_print_level_is_rebuilt(
    tmp_path: Path,
) -> None:
    """A logger that has handlers but lacks `print` would still fail recipes, so it must be rebuilt."""
    # ARRANGE
    settings = logging_settings(tmp_path / "noprint.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        logging.basicConfig()
        bareLog = logging.getLogger("fundamentals.logs")
        workerLog = reducer._worker_log(bareLog, settings)
        print("HAS PRINT:", hasattr(workerLog, "print"))
        """,
    )

    # ASSERT
    assert "HAS PRINT: True" in completed.stdout
