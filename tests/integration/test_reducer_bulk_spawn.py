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
    import logging, queue, sys, yaml
    from importlib import import_module

    reducer = import_module("soxspipe.commonutils.reducer")
    settings = yaml.safe_load(open(sys.argv[1]))
    logQueue = queue.Queue()
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
    # THE PARENT OWNS THE LOG FILE: THE WORKERS' RECORDS REACH IT THROUGH THE PARENT'S HANDLERS
    parentHandler = logging.FileHandler(logPath)
    parentHandler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
    logging.getLogger().addHandler(parentHandler)

    try:
        _call_run_recipe_bulk(reducerModule, tmp_path, logPath, commands, connection)
    finally:
        logging.getLogger().removeHandler(parentHandler)
        parentHandler.close()
    return results


def _call_run_recipe_bulk(
    reducerModule: Any, tmp_path: Path, logPath: Path, commands: list[str], connection: Any
) -> None:
    """Call `run_recipe_bulk` with one unsupported-recipe SOF per command."""
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
    """The worker payload survives a real spawn and the worker's records reach the parent's log file.

    A `logging.Logger` pickles by name, so the unpickled copy in a spawned worker has no
    handlers and no `print` level. The worker must give it the print level and send its records to the parent.
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
    """One worker process handles both tasks, so a second queue handler would double every record."""
    # ARRANGE
    logPath = tmp_path / "reused.log"

    # ACT
    results = _run_real_pool(monkeypatch, tmp_path, logPath, ["first-command", "second-command"])

    # ASSERT
    logText = logPath.read_text()
    assert [result["status"] for result in results] == ["fail", "fail"]
    assert logText.count("Recipe Command: first-command") == 1
    assert logText.count("Recipe Command: second-command") == 1


def test_worker_log_from_a_bare_logger_has_print_level_and_queues_its_records(
    tmp_path: Path,
) -> None:
    """A fresh process given only a by-name logger gets the print level and sends its records to the queue."""
    # ARRANGE
    settings = logging_settings(tmp_path / "bare.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings, logQueue)
        workerLog.print("printed from the worker")
        workerLog.error("errored from the worker")
        while not logQueue.empty():
            record = logQueue.get_nowait()
            print("QUEUED:", record.levelname, record.getMessage())
        """,
    )

    # ASSERT
    assert "QUEUED: PRINT printed from the worker" in completed.stdout
    assert "QUEUED: ERROR errored from the worker" in completed.stdout


def test_worker_log_writes_nothing_to_its_own_terminal(
    tmp_path: Path,
) -> None:
    """The parent shows the warnings and errors, so a worker with its own terminal handler would show them twice."""
    # ARRANGE
    settings = logging_settings(tmp_path / "quiet.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings, logQueue)
        workerLog.print("printed from the worker")
        workerLog.error("errored from the worker")
        """,
    )

    # ASSERT
    assert "from the worker" not in completed.stderr
    assert "from the worker" not in completed.stdout


def test_worker_log_leaves_the_shared_log_file_alone(
    tmp_path: Path,
) -> None:
    """Only the parent opens the shared log, so a worker that logs must not touch the file."""
    # ARRANGE
    logPath = tmp_path / "shared.log"
    logPath.write_text("written by the parent\n")

    # ACT
    _run_in_fresh_python(
        tmp_path,
        logging_settings(logPath),
        """
        workerLog = reducer._worker_log(logging.getLogger("fundamentals.logs"), settings, logQueue)
        workerLog.error("errored from the worker")
        for handler in logging.getLogger().handlers:
            handler.flush()
        """,
    )

    # ASSERT
    assert logPath.read_text() == "written by the parent\n"
    assert sorted(path.name for path in tmp_path.glob("shared.log*")) == ["shared.log"]


def test_worker_log_leaves_other_loggers_enabled(
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
        reducer._worker_log(logging.getLogger("fundamentals.logs"), settings, logQueue)
        print("OTHER DISABLED:", otherLog.disabled)
        """,
    )

    # ASSERT
    assert "OTHER DISABLED: False" in completed.stdout


def test_worker_log_gives_a_logger_that_has_handlers_but_no_print_level_the_print_level(
    tmp_path: Path,
) -> None:
    """A logger that has handlers but lacks `print` would still fail recipes, so it must be replaced."""
    # ARRANGE
    settings = logging_settings(tmp_path / "noprint.log")

    # ACT
    completed = _run_in_fresh_python(
        tmp_path,
        settings,
        """
        logging.basicConfig()
        bareLog = logging.getLogger("fundamentals.logs")
        workerLog = reducer._worker_log(bareLog, settings, logQueue)
        print("HAS PRINT:", hasattr(workerLog, "print"))
        """,
    )

    # ASSERT
    assert "HAS PRINT: True" in completed.stdout
