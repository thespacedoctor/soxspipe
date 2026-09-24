"""Characterization tests for the command-line dispatch boundary."""

from __future__ import annotations

import sys
from importlib import import_module
from pathlib import Path
from typing import Any, ClassVar

import pytest
from docopt import docopt as parse_docopt

pytestmark = pytest.mark.unit


class CliLogger:
    """Record logging calls made by the command-line adapter."""

    def __init__(self) -> None:
        self.messages: list[tuple[str, str]] = []

    def _record(self, level: str, message: object) -> None:
        self.messages.append((level, str(message)))

    def debug(self, message: object, **kwargs: object) -> None:
        self._record("debug", message)

    def info(self, message: object, **kwargs: object) -> None:
        self._record("info", message)

    def error(self, message: object, **kwargs: object) -> None:
        self._record("error", message)

    def print(self, message: object, **kwargs: object) -> None:
        self._record("print", message)


class RecordingOrganiser:
    """Record data-organizer commands without touching a real workspace."""

    calls: ClassVar[list[tuple[str, object]]] = []
    rawPaths: ClassVar[list[str]] = []
    currentSession: ClassVar[str | None] = None

    def __init__(self, **kwargs: object) -> None:
        self.calls.append(("init", kwargs))

    def session_list(self, **kwargs: object) -> tuple[str | None, list[str]]:
        self.calls.append(("session_list", kwargs))
        sessions = [self.currentSession] if self.currentSession else []
        return self.currentSession, sessions

    def prepare(self, **kwargs: object) -> None:
        self.calls.append(("prepare", kwargs))

    def session_switch(self, sessionId: str) -> None:
        self.calls.append(("session_switch", sessionId))

    def session_create(self, sessionId: str | None) -> str:
        self.calls.append(("session_create", sessionId))
        return sessionId or "generated"

    def list_obs(self) -> None:
        self.calls.append(("list_obs", None))

    def list_sofs(self) -> None:
        self.calls.append(("list_sofs", None))

    def list_raw(self, sofFile: str) -> tuple[list[str], Any]:
        import pandas as pd

        self.calls.append(("list_raw", sofFile))
        table = pd.DataFrame(
            [{"file": Path(path).name, "tag": "RAW"} for path in self.rawPaths],
            columns=["file", "tag"],
        )
        return list(self.rawPaths), table

    def close(self) -> None:
        self.calls.append(("close", None))


class RecordingDaemon:
    """Prevent daemonization while recording lifecycle commands."""

    calls: ClassVar[list[tuple[str, object]]] = []
    instances: ClassVar[list[RecordingDaemon]] = []

    def __init__(self, **kwargs: object) -> None:
        self.log = kwargs["log"]
        self.calls.append(("init", kwargs))
        self.instances.append(self)

    def start(self) -> None:
        self.calls.append(("start", None))

    def stop(self) -> None:
        self.calls.append(("stop", None))

    def status(self) -> None:
        self.calls.append(("status", None))


def _run_cli(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    argv: list[str],
    *,
    suppliedArguments: object = None,
    argumentOverrides: dict[str, object] | None = None,
    rawPaths: list[str] | None = None,
    currentSession: str | None = None,
) -> tuple[CliLogger, list[dict[str, object]]]:
    clUtils = import_module("soxspipe.cl_utils")
    commonutilsModule = import_module("soxspipe.commonutils")
    fundamentalsModule = import_module("fundamentals")
    parsedArguments = parse_docopt(clUtils.__doc__, argv=argv[1:])
    if argumentOverrides:
        parsedArguments = {**parsedArguments, **argumentOverrides}
    parserCalls: list[dict[str, object]] = []
    logger = CliLogger()

    def fake_docopt(docString: str, **kwargs: object) -> dict[str, object]:
        parserCalls.append(kwargs)
        return dict(parsedArguments)

    class RecordingTools:
        def __init__(self, **kwargs: object) -> None:
            self.arguments = kwargs["arguments"]

        def setup(self) -> tuple[dict[str, object], dict[str, object], CliLogger, None]:
            return (
                dict(parsedArguments),
                {"workspace-root-dir": str(tmp_path)},
                logger,
                parsedArguments.get("--dbConn"),
            )

    class FixedTimes:
        @staticmethod
        def get_now_sql_datetime() -> str:
            return "2026-09-10 00:00:00"

        @staticmethod
        def calculate_time_difference(start: str, end: str) -> str:
            return "0s"

    RecordingOrganiser.calls = []
    RecordingOrganiser.rawPaths = list(rawPaths or [])
    RecordingOrganiser.currentSession = currentSession
    RecordingDaemon.calls = []
    RecordingDaemon.instances = []
    monkeypatch.setattr(sys, "argv", argv.copy())
    monkeypatch.setattr(clUtils, "docopt", fake_docopt)
    monkeypatch.setattr(clUtils, "tools", RecordingTools)
    monkeypatch.setattr(clUtils, "times", FixedTimes)
    monkeypatch.setattr(commonutilsModule, "data_organiser", RecordingOrganiser)
    monkeypatch.setattr(fundamentalsModule, "daemonise", RecordingDaemon)
    monkeypatch.setattr(clUtils.readline, "set_completer_delims", lambda value: None)
    monkeypatch.setattr(clUtils.readline, "parse_and_bind", lambda value: None)
    monkeypatch.setattr(clUtils.readline, "set_completer", lambda value: None)
    monkeypatch.setattr(
        clUtils.pickle,
        "dump",
        lambda *args, **kwargs: pytest.fail("CLI attempted to write a pickle"),
    )
    originalChdir = clUtils.os.chdir

    def guarded_chdir(path: str) -> None:
        resolvedPath = (Path.cwd() / path).resolve()
        if not resolvedPath.is_relative_to(tmp_path.resolve()):
            pytest.fail(f"CLI attempted to leave the test workspace: {resolvedPath}")
        originalChdir(path)

    monkeypatch.setattr(clUtils.os, "chdir", guarded_chdir)
    monkeypatch.setattr(
        clUtils.time,
        "sleep",
        lambda seconds: pytest.fail("CLI attempted to sleep"),
    )

    clUtils.main(arguments=suppliedArguments)
    return logger, parserCalls


def test_tab_complete_returns_matching_path_then_sentinel(
    tmp_path: Path,
) -> None:
    """Expose filesystem completions followed by the readline sentinel."""
    clUtils = import_module("soxspipe.cl_utils")
    prefix = tmp_path / "calibration"
    (tmp_path / "calibration-file.fits").touch()

    assert clUtils.tab_complete(str(prefix), 0) == str(
        tmp_path / "calibration-file.fits"
    )
    assert clUtils.tab_complete(str(prefix), 1) is None


@pytest.mark.parametrize(
    ("argv", "recipeAttribute", "extraArguments"),
    [
        (["soxspipe", "mbias", "frames"], "soxs_mbias", {}),
        (["soxspipe", "mdark", "frames"], "soxs_mdark", {}),
        (
            ["soxspipe", "disp_solution", "frames", "--poly=1234"],
            "soxs_disp_solution",
            {"polyOrders": "1234"},
        ),
        (
            ["soxspipe", "order_centres", "frames", "--poly=12"],
            "soxs_order_centres",
            {"polyOrders": "12"},
        ),
        (["soxspipe", "mflat", "frames"], "soxs_mflat", {}),
        (
            ["soxspipe", "spat_solution", "frames", "--poly=123456"],
            "soxs_spatial_solution",
            {"polyOrders": "123456"},
        ),
        (["soxspipe", "stare_std", "frames"], "soxs_stare", {}),
        (["soxspipe", "nod_std", "frames"], "soxs_nod", {}),
        (["soxspipe", "offset", "frames"], "soxs_offset", {}),
    ],
)
def test_main_dispatches_recipe_commands_to_recording_adapters(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    argv: list[str],
    recipeAttribute: str,
    extraArguments: dict[str, object],
) -> None:
    """Pass CLI recipe inputs and flags to the matching recipe constructor."""
    recipesModule = import_module("soxspipe.recipes")
    constructorCalls: list[dict[str, object]] = []
    produceCalls = 0

    class RecordingRecipe:
        def __init__(self, **kwargs: object) -> None:
            constructorCalls.append(kwargs)

        def produce_product(self) -> str:
            nonlocal produceCalls
            produceCalls += 1
            return "product.fits"

    monkeypatch.setattr(recipesModule, recipeAttribute, RecordingRecipe)

    _run_cli(monkeypatch, tmp_path, argv)

    assert produceCalls == 1
    assert len(constructorCalls) == 1
    assert constructorCalls[0] == {
        "log": constructorCalls[0]["log"],
        "settings": {"workspace-root-dir": str(tmp_path)},
        "inputFrames": "frames",
        "verbose": False,
        "overwrite": False,
        "command": " ".join(argv),
        "debug": False,
        **extraArguments,
    }
    assert not any(
        call[0] in {"start", "stop", "status"} for call in RecordingDaemon.calls
    )


def test_main_ignores_supplied_arguments_and_reparses_process_argv(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Characterize the current ``arguments`` parameter behavior."""
    recipesModule = import_module("soxspipe.recipes")
    inputs: list[str] = []

    class RecordingRecipe:
        def __init__(self, **kwargs: object) -> None:
            inputs.append(str(kwargs["inputFrames"]))

        def produce_product(self) -> None:
            return None

    monkeypatch.setattr(recipesModule, "soxs_mbias", RecordingRecipe)

    _, parserCalls = _run_cli(
        monkeypatch,
        tmp_path,
        ["soxspipe", "mbias", "from-argv"],
        suppliedArguments={"mbias": False, "inputFrames": "from-parameter"},
    )

    assert inputs == ["from-argv"]
    assert parserCalls == [{}]


def test_main_reduce_all_converts_batch_and_forwards_flags(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Construct the reducer and forward normalized batch and multiprocessing flags."""
    commonutilsModule = import_module("soxspipe.commonutils")
    (tmp_path / "soxspipe.db").touch()
    initCalls: list[dict[str, object]] = []
    reduceCalls: list[dict[str, object]] = []

    class RecordingReducer:
        def __init__(self, **kwargs: object) -> None:
            initCalls.append(kwargs)

        def reduce(self, **kwargs: object) -> None:
            reduceCalls.append(kwargs)

    monkeypatch.setattr(commonutilsModule, "reducer", RecordingReducer)

    _run_cli(
        monkeypatch,
        tmp_path,
        ["soxspipe", "-qm", "reduce", "all"],
        argumentOverrides={"<batchSize>": "2", "--batch": True},
    )

    assert initCalls == [
        {
            "log": initCalls[0]["log"],
            "workspaceDirectory": ".",
            "reductionTarget": "all",
            "settings": {"workspace-root-dir": str(tmp_path)},
            "pathToSettings": "./soxspipe.yaml",
            "quitOnFail": True,
            "overwrite": False,
            "verbose": False,
            "refreshWorkspace": False,
        }
    ]
    assert reduceCalls == [{"batch": 2, "multiprocess": True}]


def test_main_reduce_returns_before_dispatch_when_workspace_is_unprepared(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Reject reduction when the workspace database is absent."""
    _run_cli(monkeypatch, tmp_path, ["soxspipe", "reduce", "all"])

    assert "Please run the 'soxspipe prep" in capsys.readouterr().out
    assert RecordingDaemon.calls == []


@pytest.mark.parametrize(
    ("argv", "expectedCall"),
    [
        (["soxspipe", "prep"], ("prepare", {"refresh": False})),
        (["soxspipe", "session", "ls"], ("session_list", {})),
        (["soxspipe", "session", "new", "night1"], ("session_create", "night1")),
        (["soxspipe", "session", "night1"], ("session_switch", "night1")),
        (["soxspipe", "list", "ob"], ("list_obs", None)),
        (["soxspipe", "list", "sof"], ("list_sofs", None)),
    ],
)
def test_main_dispatches_workspace_commands_to_data_organiser(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    argv: list[str],
    expectedCall: tuple[str, object],
) -> None:
    """Delegate workspace and session commands to the organizer adapter."""
    _run_cli(monkeypatch, tmp_path, argv)

    assert expectedCall in RecordingOrganiser.calls


@pytest.mark.parametrize("rawPaths", [[], ["raw-one.fits", "raw-two.fits"]])
def test_main_raw_reports_inventory_and_exports_available_files(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
    rawPaths: list[str],
) -> None:
    """Report an empty raw inventory or copy its files into the workspace export."""
    resolvedPaths: list[str] = []
    for rawName in rawPaths:
        rawPath = tmp_path / "raw" / rawName
        rawPath.parent.mkdir(exist_ok=True)
        rawPath.write_text(rawName)
        resolvedPaths.append(str(rawPath))

    _run_cli(
        monkeypatch,
        tmp_path,
        ["soxspipe", "raw", "sof", "science.sof"],
        rawPaths=resolvedPaths,
    )

    output = capsys.readouterr().out
    assert ("list_raw", "science.sof") in RecordingOrganiser.calls
    if resolvedPaths:
        assert "Exported 2 raw frames" in output
        assert (tmp_path / "exported" / "raw-one.fits").read_text() == "raw-one.fits"
        assert (tmp_path / "exported" / "raw-two.fits").read_text() == "raw-two.fits"
    else:
        assert "No raw frames found" in output
        assert not (tmp_path / "exported").exists()


def test_main_raw_rejects_frame_outside_workspace_raw_directory(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Do not export a database-supplied path outside the raw-frame directory."""
    sourcePath = tmp_path / "source" / "not-a-raw-frame.fits"
    sourcePath.parent.mkdir()
    sourcePath.write_text("private", encoding="utf-8")

    with pytest.raises(SystemExit) as error:
        _run_cli(
            monkeypatch,
            tmp_path,
            ["soxspipe", "raw", "sof", "science.sof"],
            rawPaths=[str(sourcePath)],
        )

    assert error.value.code == 1
    assert not (tmp_path / "exported" / sourcePath.name).exists()


def test_main_exits_cleanly_when_initial_session_lookup_is_unsafe(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Translate unsafe paths encountered before normal CLI setup into exit 1."""
    from soxspipe.commonutils.data_organiser import _UnsafePathError

    def reject_unsafe_workspace(self: object, **kwargs: object) -> None:
        raise _UnsafePathError("Unsafe sessions directory")

    monkeypatch.setattr(RecordingOrganiser, "__init__", reject_unsafe_workspace)

    with pytest.raises(SystemExit) as error:
        _run_cli(monkeypatch, tmp_path, ["soxspipe", "prep"])

    assert error.value.code == 1


def test_main_raw_rejects_export_directory_symlink_outside_workspace(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Do not export raw frames through a directory symlink outside the workspace."""
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    rawPath = tmp_path / "source" / "raw-one.fits"
    rawPath.parent.mkdir()
    rawPath.write_text("raw", encoding="utf-8")
    outsidePath = tmp_path / "outside"
    outsidePath.mkdir()
    (workspacePath / "exported").symlink_to(outsidePath, target_is_directory=True)

    with pytest.raises(SystemExit) as error:
        _run_cli(
            monkeypatch,
            tmp_path,
            ["soxspipe", "raw", "sof", "science.sof", str(workspacePath)],
            rawPaths=[str(rawPath)],
        )

    assert error.value.code == 1
    assert not (outsidePath / rawPath.name).exists()


def test_main_raw_rejects_destination_symlink_outside_workspace(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Do not overwrite or accept a raw-export symlink targeting another directory."""
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    rawPath = tmp_path / "source" / "raw-one.fits"
    rawPath.parent.mkdir()
    rawPath.write_text("raw", encoding="utf-8")
    exportPath = workspacePath / "exported"
    exportPath.mkdir()
    outsidePath = tmp_path / "outside.fits"
    outsidePath.write_text("outside", encoding="utf-8")
    (exportPath / rawPath.name).symlink_to(outsidePath)

    with pytest.raises(SystemExit) as error:
        _run_cli(
            monkeypatch,
            tmp_path,
            ["soxspipe", "raw", "sof", "science.sof", str(workspacePath)],
            rawPaths=[str(rawPath)],
        )

    assert error.value.code == 1
    assert outsidePath.read_text(encoding="utf-8") == "outside"


def test_main_changes_only_into_requested_test_workspace(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Honor an explicit workspace before normalizing CLI paths to the current directory."""
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()

    _run_cli(
        monkeypatch,
        tmp_path,
        ["soxspipe", "prep", str(workspacePath)],
    )

    assert Path.cwd() == workspacePath


def test_main_uses_current_session_settings_when_no_path_is_supplied(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Append the active session settings file to the process command."""
    recipesModule = import_module("soxspipe.recipes")

    class RecordingRecipe:
        def __init__(self, **kwargs: object) -> None:
            pass

        def produce_product(self) -> None:
            return None

    monkeypatch.setattr(recipesModule, "soxs_mbias", RecordingRecipe)

    _run_cli(
        monkeypatch,
        tmp_path,
        ["soxspipe", "mbias", "frames"],
        currentSession="night1",
    )

    assert sys.argv[-2:] == ["-s", "./sessions/night1/soxspipe.yaml"]


@pytest.mark.parametrize("daemonCommand", ["start", "stop", "status"])
def test_main_watch_dispatches_daemon_lifecycle_without_starting_action(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    daemonCommand: str,
) -> None:
    """Delegate watch lifecycle commands to the daemon adapter."""
    _run_cli(monkeypatch, tmp_path, ["soxspipe", "watch", daemonCommand])

    assert (daemonCommand, None) in RecordingDaemon.calls


def test_main_file_exists_failure_exits_successfully(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Preserve the CLI's successful exit status for existing products."""
    recipesModule = import_module("soxspipe.recipes")

    class ExistingRecipe:
        def __init__(self, **kwargs: object) -> None:
            pass

        def produce_product(self) -> None:
            raise FileExistsError("already exists")

    monkeypatch.setattr(recipesModule, "soxs_mbias", ExistingRecipe)

    with pytest.raises(SystemExit) as error:
        _run_cli(monkeypatch, tmp_path, ["soxspipe", "mbias", "frames"])

    assert error.value.code == 0


def test_main_logs_recipe_failure_then_finishes_dispatch(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Log unexpected recipe errors without re-raising them."""
    recipesModule = import_module("soxspipe.recipes")

    class FailingRecipe:
        def __init__(self, **kwargs: object) -> None:
            raise RuntimeError("recipe failed")

    monkeypatch.setattr(recipesModule, "soxs_mbias", FailingRecipe)

    logger, _ = _run_cli(monkeypatch, tmp_path, ["soxspipe", "mbias", "frames"])

    assert any(
        level == "error" and "recipe failed\nsoxspipe mbias frames" in message
        for level, message in logger.messages
    )
