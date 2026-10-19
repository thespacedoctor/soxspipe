"""The bulk reduction pool survives a worker killed by a signal (DY-1615)."""

from __future__ import annotations

import contextlib
import multiprocessing
import os
import signal
import time
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from tests.factories import logging_settings

pytestmark = pytest.mark.unit

# A CRASHED WORKER USED TO LEAVE THE PARENT WAITING FOR TEN HOURS, SO ANY TEST THAT RETURNS IS FAST
MAX_POOL_SECONDS = 60
# THE POOL'S OWN DEADLINE IS THE BOUND ON EVERY TEST HERE: A REGRESSION FAILS AFTER THIS LONG INSTEAD OF HANGING CI
TEST_POOL_TIMEOUT_SECONDS = 90
# LONG ENOUGH FOR A SPAWNED WORKER TO START AND WRITE ITS PID FILE ON A LOADED CI RUNNER
SHORT_POOL_TIMEOUT_SECONDS = 15
PROCESS_EXIT_SECONDS = 10
CRASH_COMMAND = "crash"
SIBLING_COMMAND = "sibling"
SIBLING_RUN_SECONDS = 5
SIGNAL_WAIT_SECONDS = 30


@pytest.fixture(autouse=True)
def bounded_pool(monkeypatch: pytest.MonkeyPatch) -> int:
    """Give the pool a deadline well under the ten hours of production, so a hung pool fails the test.

    Return the production deadline.
    """
    reducerModule = import_module("soxspipe.commonutils.reducer")
    productionTimeout = reducerModule.POOL_TIMEOUT_SECONDS
    monkeypatch.setattr(reducerModule, "POOL_TIMEOUT_SECONDS", TEST_POOL_TIMEOUT_SECONDS)
    return productionTimeout


def reduce_or_crash(inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
    """Return a pass result, or kill this worker with SIGSEGV when the command asks for it."""
    if inputDict["command"] == CRASH_COMMAND:
        os.kill(os.getpid(), signal.SIGSEGV)
    return {"status": "pass", "sof": inputDict["sof"], "worker": os.getpid(), "error_message": None, **kwargs}


def return_an_unpicklable_result(inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
    """Return a result that cannot be sent back to the parent process."""
    return {"status": "pass", "sof": inputDict["sof"], "unpicklable": lambda: None}


def report_start_method(inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
    """Return the start method this worker process was created with."""
    return {"status": "pass", "sof": inputDict["sof"], "startMethod": multiprocessing.get_start_method()}


def sleep_after_recording_pid(inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
    """Write this worker's PID to the path in the input, then sleep far longer than any test waits."""
    Path(inputDict["pidPath"]).write_text(str(os.getpid()))
    time.sleep(600)
    return {"status": "pass", "sof": inputDict["sof"]}


def crash_after_sibling_starts(inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
    """Crash once the sibling task has started, or be the sibling: signal, run for a while, then pass."""
    signalPath = Path(inputDict["signalPath"])
    if inputDict["command"] == SIBLING_COMMAND:
        signalPath.write_text("started")
        time.sleep(SIBLING_RUN_SECONDS)
        return {"status": "pass", "sof": inputDict["sof"], "error_message": None}
    deadline = time.monotonic() + SIGNAL_WAIT_SECONDS
    while not signalPath.exists() and time.monotonic() < deadline:
        time.sleep(0.05)
    os.kill(os.getpid(), signal.SIGSEGV)
    return {"status": "pass", "sof": inputDict["sof"]}


def print_and_warn(
    inputDict: dict[str, Any], log: Any, settings: dict[str, Any], logQueue: Any, **kwargs: Any
) -> dict[str, Any]:
    """Print a plain line and log a warning, with the logger set up from the settings as the real worker does."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    workerLog = reducerModule._worker_log(log, settings, logQueue)
    print("PLAIN-PRINT-MARKER", flush=True)
    workerLog.warning("WARNING-MARKER")
    return {"status": "pass", "sof": inputDict["sof"], "error_message": None}


def _fail_on_unpickle() -> None:
    """Raise when a worker unpickles the object that names this function."""
    raise RuntimeError("the worker cannot unpickle this")


class unpicklable_in_worker:
    """A callable that pickles in the parent but raises when a spawned worker unpickles it."""

    def __reduce__(self) -> tuple[Any, tuple[()]]:
        return (_fail_on_unpickle, ())

    def __call__(self, inputDict: dict[str, Any], log: Any, **kwargs: Any) -> dict[str, Any]:
        return {"status": "pass", "sof": inputDict["sof"]}


def _is_alive(pid: int) -> bool:
    """Return whether a process with this PID exists."""
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _kill_recorded_worker(pidPath: Path) -> None:
    """Kill the worker whose PID a task recorded, if it is still there, so a failed test leaks nothing."""
    if not pidPath.exists():
        return
    with contextlib.suppress(ProcessLookupError):
        os.kill(int(pidPath.read_text()), signal.SIGKILL)


def _wait_for_file(path: Path) -> None:
    """Wait until a worker has written a file, or fail the test."""
    deadline = time.monotonic() + SIGNAL_WAIT_SECONDS
    while not path.exists():
        assert time.monotonic() < deadline, f"{path} was never written"
        time.sleep(0.05)


def _inputs(*commands: str) -> list[dict[str, str]]:
    """Return one input dictionary per command."""
    return [{"sof": f"/workspace/sof/{index}.sof", "command": command} for index, command in enumerate(commands)]


def test_a_worker_killed_by_a_signal_fails_only_its_own_sof_and_the_batch_finishes(log: Any) -> None:
    # ARRANGE: ONE WORKER RUNS THE SOFS IN TURN, SO THE CRASH FALLS ON THE THIRD AND THE LAST TWO NEVER STARTED
    reducerModule = import_module("soxspipe.commonutils.reducer")
    inputDicts = _inputs("ok", "ok", CRASH_COMMAND, "ok", "ok")
    startTime = time.monotonic()

    # ACT
    results = reducerModule._run_in_spawn_pool(reduce_or_crash, inputDicts, 1, log)

    # ASSERT
    assert time.monotonic() - startTime < MAX_POOL_SECONDS
    assert [result["sof"] for result in results] == [inputDict["sof"] for inputDict in inputDicts]
    assert [result["status"] for result in results] == ["pass", "pass", "fail", "pass", "pass"]
    assert results[2]["error_message"] == reducerModule.WORKER_CRASH_MESSAGE
    assert results[2]["qcTable"] is None
    assert any(level == "error" and inputDicts[2]["sof"] in message for level, message in log.messages)


def test_two_crashes_in_one_batch_each_fail_their_own_sof(log: Any) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    inputDicts = _inputs(CRASH_COMMAND, "ok", CRASH_COMMAND, "ok")

    # ACT
    results = reducerModule._run_in_spawn_pool(reduce_or_crash, inputDicts, 1, log)

    # ASSERT
    assert [result["status"] for result in results] == ["fail", "pass", "fail", "pass"]
    assert [result["sof"] for result in results] == [inputDict["sof"] for inputDict in inputDicts]


def test_results_keep_the_input_order_and_pass_keyword_arguments_to_the_function(log: Any) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    inputDicts = _inputs("ok", "ok", "ok", "ok", "ok", "ok")

    # ACT
    results = reducerModule._run_in_spawn_pool(reduce_or_crash, inputDicts, 3, log, recipe="mbias")

    # ASSERT
    assert [result["sof"] for result in results] == [inputDict["sof"] for inputDict in inputDicts]
    assert {result["recipe"] for result in results} == {"mbias"}


def test_the_workers_are_spawned_and_the_global_start_method_is_left_alone(log: Any) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    startMethodBefore = multiprocessing.get_start_method(allow_none=True)

    # ACT
    results = reducerModule._run_in_spawn_pool(report_start_method, _inputs("ok"), 1, log)

    # ASSERT
    assert results[0]["startMethod"] == "spawn"
    assert multiprocessing.get_start_method(allow_none=True) == startMethodBefore


def test_a_result_that_cannot_be_sent_back_fails_its_sof_with_a_string_message(log: Any) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    inputDicts = _inputs("ok")

    # ACT
    results = reducerModule._run_in_spawn_pool(return_an_unpicklable_result, inputDicts, 1, log)

    # ASSERT
    assert results[0]["status"] == "fail"
    assert isinstance(results[0]["error_message"], str)
    assert results[0]["error_message"]


def test_every_sof_fails_instead_of_looping_when_the_workers_cannot_start(log: Any) -> None:
    # ARRANGE: THE FUNCTION CANNOT BE UNPICKLED, SO EVERY WORKER DIES BEFORE IT STARTS A TASK
    reducerModule = import_module("soxspipe.commonutils.reducer")
    inputDicts = _inputs("ok", "ok")

    # ACT
    results = reducerModule._run_in_spawn_pool(unpicklable_in_worker(), inputDicts, 2, log)

    # ASSERT
    assert [result["status"] for result in results] == ["fail", "fail"]
    assert {result["error_message"] for result in results} == {reducerModule.WORKER_START_MESSAGE}


def test_an_empty_batch_returns_no_results_without_starting_a_pool(log: Any) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")

    # ACT
    results = reducerModule._run_in_spawn_pool(reduce_or_crash, [], 2, log)

    # ASSERT
    assert results == []


def test_a_sibling_task_in_flight_when_another_worker_crashes_still_passes(tmp_path: Path, log: Any) -> None:
    # ARRANGE: TWO WORKERS. THE CRASHER WAITS UNTIL THE SIBLING IS RUNNING, SO BOTH ARE IN FLIGHT AT THE CRASH
    reducerModule = import_module("soxspipe.commonutils.reducer")
    signalPath = str(tmp_path / "sibling-started")
    inputDicts = [
        {"sof": "/workspace/sof/crasher.sof", "command": CRASH_COMMAND, "signalPath": signalPath},
        {"sof": "/workspace/sof/sibling.sof", "command": SIBLING_COMMAND, "signalPath": signalPath},
    ]

    # ACT
    results = reducerModule._run_in_spawn_pool(crash_after_sibling_starts, inputDicts, 2, log)

    # ASSERT
    assert [result["sof"] for result in results] == [inputDict["sof"] for inputDict in inputDicts]
    assert [result["status"] for result in results] == ["fail", "pass"]
    assert results[0]["error_message"] == reducerModule.WORKER_CRASH_MESSAGE


def test_the_default_pool_deadline_is_ten_hours(bounded_pool: int) -> None:
    # ASSERT: THE FIXTURE RETURNS THE DEADLINE IT FOUND BEFORE IT SHORTENED IT FOR THE TEST
    assert bounded_pool == 36000


def test_a_sleeping_worker_is_stopped_and_its_sof_fails_when_the_pool_times_out(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: Any
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    monkeypatch.setattr(reducerModule, "POOL_TIMEOUT_SECONDS", SHORT_POOL_TIMEOUT_SECONDS)
    pidPath = tmp_path / "worker.pid"
    inputDicts = [{"sof": "/workspace/sof/slow.sof", "command": "ok", "pidPath": str(pidPath)}]
    startTime = time.monotonic()

    # ACT
    try:
        results = reducerModule._run_in_spawn_pool(sleep_after_recording_pid, inputDicts, 1, log)

        # ASSERT: BEFORE THE CLEANUP, WHICH WOULD HIDE A LEAKED WORKER
        assert pidPath.exists()
        assert not _is_alive(int(pidPath.read_text()))
    finally:
        _kill_recorded_worker(pidPath)
    assert time.monotonic() - startTime < MAX_POOL_SECONDS
    assert results[0]["status"] == "fail"
    assert results[0]["error_message"] == reducerModule.WORKER_TIMEOUT_MESSAGE


def test_an_interrupt_stops_the_workers_and_is_raised_again(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, log: Any
) -> None:
    # ARRANGE: THE INTERRUPT ARRIVES WHILE THE WORKER IS RUNNING
    reducerModule = import_module("soxspipe.commonutils.reducer")
    import concurrent.futures

    pidPath = tmp_path / "worker.pid"

    def interrupt_once_started(futures: Any, timeout: float | None = None) -> Any:
        _wait_for_file(pidPath)
        raise KeyboardInterrupt

    # THIS REACHES THE POOL BECAUSE _collect_round IMPORTS as_completed INSIDE THE FUNCTION, AT CALL TIME
    monkeypatch.setattr(concurrent.futures, "as_completed", interrupt_once_started)
    inputDicts = [{"sof": "/workspace/sof/slow.sof", "command": "ok", "pidPath": str(pidPath)}]

    # ACT
    try:
        with pytest.raises(KeyboardInterrupt):
            reducerModule._run_in_spawn_pool(sleep_after_recording_pid, inputDicts, 1, log)

        # ASSERT: BEFORE THE CLEANUP, WHICH WOULD HIDE A LEAKED WORKER
        assert pidPath.exists()
        assert not _is_alive(int(pidPath.read_text()))
    finally:
        _kill_recorded_worker(pidPath)


def _task(reducerModule: Any, log: Any, count: int) -> Any:
    """Return a pool task over ``count`` inputs, for tests that call the settling step directly."""
    return reducerModule._pool_task(reduce_or_crash, _inputs(*["ok"] * count), log, {}, "unused-marker-directory")


def _forbid_pool_rounds(monkeypatch: pytest.MonkeyPatch, reducerModule: Any) -> None:
    """Make any attempt to run a pool round fail the test."""

    def refuse_round(*args: Any, **kwargs: Any) -> Any:
        raise AssertionError("a pool round was started")

    monkeypatch.setattr(reducerModule, "_run_pool_round", refuse_round)


def test_in_flight_sofs_fail_as_timed_out_without_a_pool_when_the_deadline_has_passed(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    _forbid_pool_rounds(monkeypatch, reducerModule)
    task = _task(reducerModule, log, 3)

    # ACT
    settled = reducerModule._settle_in_flight(task, [0, 1, 2], time.monotonic() - 1)

    # ASSERT
    assert sorted(settled) == [0, 1, 2]
    assert {result["error_message"] for result in settled.values()} == {reducerModule.WORKER_TIMEOUT_MESSAGE}
    assert sum("deadline" in message for level, message in log.messages if level == "error") == 3


def test_a_lone_in_flight_sof_is_taken_to_be_the_crasher_without_a_second_run(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE
    reducerModule = import_module("soxspipe.commonutils.reducer")
    _forbid_pool_rounds(monkeypatch, reducerModule)
    task = _task(reducerModule, log, 4)

    # ACT
    settled = reducerModule._settle_in_flight(task, [3], time.monotonic() + TEST_POOL_TIMEOUT_SECONDS)

    # ASSERT
    assert settled[3]["error_message"] == reducerModule.WORKER_CRASH_MESSAGE


def test_an_isolated_sof_whose_pool_cannot_start_fails_with_the_start_message(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE: THE ISOLATION POOL BREAKS BEFORE THE SOF ENTERS ITS FUNCTION, SO IT IS NEITHER RESOLVED NOR IN FLIGHT
    reducerModule = import_module("soxspipe.commonutils.reducer")
    monkeypatch.setattr(
        reducerModule,
        "_run_pool_round",
        lambda task, indexes, poolSize, timeout: reducerModule._round_outcome({}, list(indexes), []),
    )
    task = _task(reducerModule, log, 2)

    # ACT
    settled = reducerModule._settle_in_flight(task, [0, 1], time.monotonic() + TEST_POOL_TIMEOUT_SECONDS)

    # ASSERT
    assert {result["error_message"] for result in settled.values()} == {reducerModule.WORKER_START_MESSAGE}


def test_a_worker_warning_reaches_the_terminal_while_a_plain_print_does_not(
    tmp_path: Path, log: Any, capfd: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE: THE PARENT'S CONSOLE HANDLER IS THE ONLY TERMINAL OUTPUT, AS IN A REAL RUN
    import logging
    import sys

    reducerModule = import_module("soxspipe.commonutils.reducer")
    settings = logging_settings(tmp_path / "worker.log")
    parentLog = logging.getLogger("fundamentals.logs")
    consoleHandler = logging.StreamHandler(sys.stdout)
    parentLog.addHandler(consoleHandler)

    # ACT
    try:
        with reducerModule._worker_log_queue(parentLog) as logQueue:
            reducerModule._run_in_spawn_pool(
                print_and_warn, _inputs("ok"), 1, parentLog, settings=settings, logQueue=logQueue
            )
    finally:
        parentLog.removeHandler(consoleHandler)

    # ASSERT
    terminalText = "".join(capfd.readouterr())
    assert "WARNING-MARKER" in terminalText
    assert "PLAIN-PRINT-MARKER" not in terminalText
