"""Integration contracts for reducer selection and top-level control flow."""

from __future__ import annotations

import sqlite3
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytestmark = pytest.mark.integration


def _uninitialized_reducer(log: Any, databasePath: Path) -> object:
    reducerModule = import_module("soxspipe.commonutils.reducer")
    collection = reducerModule.reducer.__new__(reducerModule.reducer)
    collection.log = log
    collection.sessionDB = str(databasePath)
    collection.workspaceDirectory = str(databasePath.parent)
    collection.pathToSettings = False
    return collection


def _create_selection_database(databasePath: Path) -> None:
    with sqlite3.connect(databasePath) as connection:
        connection.execute(
            'CREATE TABLE raw_frame_sets (recipe_order INTEGER, complete INTEGER, recipe TEXT, sof TEXT, "eso seq arm" TEXT)'
        )
        connection.execute(
            "CREATE TABLE product_frames (sof TEXT, file TEXT, complete INTEGER)"
        )
        connection.execute("CREATE TABLE sof_map_base (sof TEXT, file TEXT)")


def test_select_all_filters_orders_limits_and_builds_commands(
    tmp_path: Path,
    log: Any,
) -> None:
    """Read processable SOFs from SQLite in stable recipe and filename order."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.executemany(
            "INSERT INTO raw_frame_sets VALUES (?, ?, ?, ?, ?)",
            [
                (2, 1, "mflat", "flat-b.sof", "UVB"),
                (1, 1, "mflat", "flat-a.sof", "UVB"),
                (1, 1, "mbias", "bias.sof", "UVB"),
                (0, 0, "mflat", "incomplete.sof", "UVB"),
                (0, 1, "mflat", "nir.sof", "NIR"),
            ],
        )
    collection = _uninitialized_reducer(log, databasePath)
    collection.pathToSettings = "/workspace/settings.yaml"

    selected = collection.select_sof_files_to_process(
        recipe="mflat",
        reductionTarget="all",
        batch=2,
        arm="UVB",
    )

    assert selected.to_dict("records") == [
        {
            "recipe": "mflat",
            "sof": "flat-a.sof",
            "command": "soxspipe mflat sof/flat-a.sof -s /workspace/settings.yaml",
        },
        {
            "recipe": "mflat",
            "sof": "flat-b.sof",
            "command": "soxspipe mflat sof/flat-b.sof -s /workspace/settings.yaml",
        },
    ]


def test_select_sof_includes_recursive_calibration_dependencies(
    tmp_path: Path,
    log: Any,
) -> None:
    """Resolve upstream SOFs through the real session tables."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.executemany(
            "INSERT INTO raw_frame_sets VALUES (?, ?, ?, ?, ?)",
            [
                (1, 1, "mbias", "bias.sof", "UVB"),
                (2, 1, "mflat", "flat.sof", "UVB"),
                (3, 1, "stare_obj", "science.sof", "UVB"),
            ],
        )
        connection.executemany(
            "INSERT INTO product_frames VALUES (?, ?, ?)",
            [
                ("bias.sof", "bias.fits", 1),
                ("flat.sof", "flat.fits", 1),
                ("science.sof", "science.fits", 1),
            ],
        )
        connection.executemany(
            "INSERT INTO sof_map_base VALUES (?, ?)",
            [
                ("flat.sof", "bias.fits"),
                ("science.sof", "flat.fits"),
            ],
        )
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="science.sof")

    assert selected["sof"].tolist() == ["bias.sof", "flat.sof", "science.sof"]
    assert selected["recipe"].tolist() == ["mbias", "mflat", "stare_obj"]
    assert selected["command"].tolist() == [
        "soxspipe mbias sof/bias.sof",
        "soxspipe mflat sof/flat.sof",
        "soxspipe stare sof/science.sof",
    ]


def test_select_empty_target_returns_stable_columns_and_warning(
    tmp_path: Path,
    log: Any,
) -> None:
    """Return an empty dispatch frame when a requested SOF is unavailable."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="missing.sof")

    assert selected.empty
    assert selected.columns.tolist() == ["recipe", "sof", "command"]
    assert (
        "warning",
        "The SOF file selected for processing is either missing or incomplete.",
    ) in log.messages


def test_select_unknown_reduction_target_preserves_unbound_failure(
    tmp_path: Path,
    log: Any,
) -> None:
    """Characterize the current failure for unsupported target forms."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    collection = _uninitialized_reducer(log, databasePath)

    with pytest.raises(UnboundLocalError):
        collection.select_sof_files_to_process(reductionTarget="ob-123")


def test_reduce_without_active_session_reports_preparation_requirement(
    tmp_path: Path,
    log: Any,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """Stop before workspace access when no reduction session exists."""
    collection = _uninitialized_reducer(log, tmp_path / "unused.db")
    collection.sessionId = None

    assert collection.reduce() is None
    assert "Please prepare this workspace" in capsys.readouterr().out


def test_reduce_serial_dispatches_rows_and_refreshes_after_failure(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """Run serial rows in order while preserving skip and failure behavior."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    commonutilsModule = import_module("soxspipe.commonutils")
    fundamentalsModule = import_module("fundamentals")
    organiserCalls: list[tuple[str, object]] = []
    runCalls: list[tuple[tuple[object, ...], dict[str, object]]] = []
    workspacePath = tmp_path / "workspace"
    collection = _uninitialized_reducer(log, workspacePath / "soxspipe.db")
    collection.sessionId = "night1"
    collection.sessionPath = str(workspacePath / "sessions" / "night1")
    collection.workspaceDirectory = str(workspacePath)
    collection.reductionTarget = "selected.sof"
    collection.recipeList = ["mbias"]
    collection.settings = {"workspace-root-dir": str(workspacePath)}
    collection.overwrite = True
    collection.verbose = True
    collection.quitOnFail = False
    collection.daemon = True
    selected = pd.DataFrame(
        [
            {"recipe": "stare_obj", "sof": "good.sof", "command": "good"},
            {"recipe": "mflat", "sof": "exists.sof", "command": "exists"},
            {"recipe": "nod_obj", "sof": "bad.sof", "command": "bad_obj"},
        ]
    )
    monkeypatch.setattr(
        collection,
        "select_sof_files_to_process",
        lambda **kwargs: selected.copy(),
    )

    class RecordingOrganiser:
        def __init__(self, **kwargs: object) -> None:
            organiserCalls.append(("init", kwargs))

        def session_refresh(self, **kwargs: object) -> bool:
            organiserCalls.append(("refresh", kwargs))
            return True

        def close(self) -> None:
            organiserCalls.append(("close", None))

    def fake_run_recipe(*args: object, **kwargs: object) -> None:
        runCalls.append((args, kwargs))
        sofPath = str(args[2])
        if sofPath.endswith("exists.sof"):
            raise FileExistsError("already exists")
        if sofPath.endswith("bad.sof"):
            raise RuntimeError("failed")

    class FixedTimes:
        @staticmethod
        def get_now_sql_datetime() -> str:
            return "2026-09-10 00:00:00"

        @staticmethod
        def calculate_time_difference(start: str, end: str) -> str:
            return "0s"

    monkeypatch.setattr(commonutilsModule, "data_organiser", RecordingOrganiser)
    monkeypatch.setattr(reducerModule, "run_recipe", fake_run_recipe)
    monkeypatch.setattr(fundamentalsModule, "times", FixedTimes)

    assert collection.reduce(batch=10, multiprocess=False) is None
    assert [Path(str(call[0][2])).name for call in runCalls] == [
        "good.sof",
        "exists.sof",
        "bad.sof",
    ]
    assert [call[0][1] for call in runCalls] == ["stare", "mflat", "nod"]
    assert collection.recipeList == [False]
    assert collection.overwrite is False
    assert ("refresh", {}) in organiserCalls
    assert sum(call == ("refresh", {"failure": None}) for call in organiserCalls) == 2
    assert any(
        level == "error" and "Recipe failed" in message
        for level, message in log.messages
    )


def test_reduce_multiprocess_forwards_selected_group(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
) -> None:
    """Pass selected SOFs and commands to the bulk reducer exactly once."""
    reducerModule = import_module("soxspipe.commonutils.reducer")
    commonutilsModule = import_module("soxspipe.commonutils")
    workspacePath = tmp_path / "workspace"
    workspacePath.mkdir()
    databasePath = workspacePath / "soxspipe.db"
    sqlite3.connect(databasePath).close()
    collection = _uninitialized_reducer(log, databasePath)
    collection.sessionId = "night1"
    collection.sessionPath = str(workspacePath / "sessions" / "night1")
    collection.workspaceDirectory = str(workspacePath)
    collection.reductionTarget = "all"
    collection.recipeList = ["mbias"]
    collection.settings = {"answer": 42}
    collection.overwrite = True
    collection.verbose = False
    collection.quitOnFail = False
    collection.daemon = False
    selectedFrames = [
        pd.DataFrame([{"recipe": "mbias", "sof": "bias.sof", "command": "run bias"}]),
        pd.DataFrame(columns=["recipe", "sof", "command"]),
    ]
    monkeypatch.setattr(
        collection,
        "select_sof_files_to_process",
        lambda **kwargs: selectedFrames.pop(0),
    )
    bulkCalls: list[dict[str, object]] = []

    class RecordingOrganiser:
        def __init__(self, **kwargs: object) -> None:
            pass

        def session_refresh(self, **kwargs: object) -> None:
            pass

        def get_incomplete_raw_frames_set(self) -> pd.DataFrame:
            return pd.DataFrame()

        def close(self) -> None:
            pass

    def fake_bulk(**kwargs: object) -> None:
        bulkCalls.append(kwargs)
        kwargs["conn"].close()

    monkeypatch.setattr(commonutilsModule, "data_organiser", RecordingOrganiser)
    monkeypatch.setattr(reducerModule, "run_recipe_bulk", fake_bulk)

    collection.reduce(batch=1, multiprocess=True)

    assert len(bulkCalls) == 1
    assert bulkCalls[0]["recipe"] == "mbias"
    assert bulkCalls[0]["sofList"] == [
        str(workspacePath / "sessions" / "night1" / "sof" / "bias.sof")
    ]
    assert bulkCalls[0]["commandList"] == ["run bias"]
    assert bulkCalls[0]["settings"] == {"answer": 42}
    assert bulkCalls[0]["overwrite"] is True
    assert bulkCalls[0]["sessionId"] == "night1"
