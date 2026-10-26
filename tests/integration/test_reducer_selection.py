"""Integration contracts for reducer selection and top-level control flow."""

from __future__ import annotations

import sqlite3
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

pytestmark = pytest.mark.integration

SESSION_ID = "20261022t120000"
STATUS_COLUMN = f"status_{SESSION_ID}"


def _uninitialized_reducer(log: Any, databasePath: Path) -> object:
    reducerModule = import_module("soxspipe.commonutils.reducer")
    collection = reducerModule.reducer.__new__(reducerModule.reducer)
    collection.log = log
    collection.sessionDB = str(databasePath)
    collection.workspaceDirectory = str(databasePath.parent)
    collection.pathToSettings = False
    collection.sessionId = SESSION_ID
    collection.overwrite = False
    return collection


def _create_selection_database(databasePath: Path) -> None:
    with sqlite3.connect(databasePath) as connection:
        connection.execute(
            'CREATE TABLE raw_frame_sets (recipe_order INTEGER, complete INTEGER, recipe TEXT, sof TEXT, "eso seq arm" TEXT)'
        )
        connection.execute(f"CREATE TABLE product_frames (sof TEXT, file TEXT, complete INTEGER, {STATUS_COLUMN} TEXT)")
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


def _insert_product_row(databasePath: Path, sofName: str | None, fileName: str, status: str | None) -> None:
    """Add one product file row of a SOF, with the given session status."""
    insertSql = f"INSERT INTO product_frames (sof, file, complete, {STATUS_COLUMN}) VALUES (?, ?, 1, ?)"  # noqa: S608
    with sqlite3.connect(databasePath) as connection:
        connection.execute(insertSql, (sofName, fileName, status))


def _seed_statuses(databasePath: Path, statuses: dict[str, str | None]) -> None:
    """Add a complete stare SOF to the selection database for each name, with the given session status."""
    for sofName, status in statuses.items():
        with sqlite3.connect(databasePath) as connection:
            connection.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'stare_obj', ?, 'VIS')", (sofName,))
        _insert_product_row(databasePath, sofName, f"{sofName}.fits", status)


def test_select_all_skips_a_sof_that_already_passed_this_session(tmp_path: Path, log: Any) -> None:
    """A per-frame stare SOF writes no predicted product, so the session status is what marks it as done (DY-1705)."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"done.sof": "pass"})
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected.empty


@pytest.mark.parametrize("status", ["fail", None])
def test_select_all_retries_a_sof_that_failed_or_never_ran(tmp_path: Path, log: Any, status: str | None) -> None:
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"retry.sof": status, "done.sof": "pass"})
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["retry.sof"]


def test_select_all_skips_a_passed_sof_even_when_another_passing_row_has_no_sof(tmp_path: Path, log: Any) -> None:
    """A NULL `sof` in `product_frames` must not make every candidate unknown, so nothing else is selected."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"done.sof": "pass", "retry.sof": "fail"})
    _insert_product_row(databasePath, None, "x.fits", "pass")
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["retry.sof"]


@pytest.mark.parametrize("otherStatus", ["fail", None])
def test_select_all_retries_a_sof_when_only_some_of_its_product_rows_passed(
    tmp_path: Path, log: Any, otherStatus: str | None
) -> None:
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"partial.sof": "pass", "done.sof": "pass"})
    _insert_product_row(databasePath, "partial.sof", "second.fits", otherStatus)
    _insert_product_row(databasePath, "done.sof", "second.fits", "pass")
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["partial.sof"]


def test_select_all_with_overwrite_selects_a_sof_that_already_passed(tmp_path: Path, log: Any) -> None:
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"done.sof": "pass"})
    collection = _uninitialized_reducer(log, databasePath)
    collection.overwrite = True

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["done.sof"]


def test_select_all_selects_a_set_that_has_no_product_row(tmp_path: Path, log: Any) -> None:
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'mbias', 'bias.sof', 'VIS')")
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["bias.sof"]


def test_select_sof_target_still_returns_a_sof_that_already_passed(tmp_path: Path, log: Any) -> None:
    """Naming a SOF on the command line is an explicit request, so the status filter is for `reduce all` only."""
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    _seed_statuses(databasePath, {"done.sof": "pass"})
    collection = _uninitialized_reducer(log, databasePath)

    selected = collection.select_sof_files_to_process(reductionTarget="done.sof")

    assert selected["sof"].tolist() == ["done.sof"]


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
            "INSERT INTO product_frames (sof, file, complete) VALUES (?, ?, ?)",
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


# ---------------------------------------------------------------------------
# DY-254 -- `recipe`, `arm` AND `reductionTarget` REACH THE DATABASE AS BOUND
# PARAMETERS, NOT INTERPOLATED SQL TEXT, IN `select_sof_files_to_process`.
# ---------------------------------------------------------------------------


def test_select_sof_files_binds_a_hostile_recipe_as_a_parameter_in_all_branch(
    tmp_path: Path,
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A hostile `recipe` filter reaches SQLite as a bound parameter, not as interpolated text.

    The method also re-filters its result in pandas by exact string match on
    `recipe` (see the `mask` step below the SQL call), which would make a
    row-count assertion vacuous here: an injected `OR '1'='1'` widens the SQL
    result, but the pandas mask then narrows it straight back down to nothing,
    since no real row's `recipe` column equals the hostile string verbatim.
    So this test inspects the query pandas actually sent instead.
    """
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'mflat', 'flat.sof', 'UVB')")
    collection = _uninitialized_reducer(log, databasePath)
    hostileRecipe = "mflat' OR '1'='1' -- "

    import pandas

    captured: dict[str, Any] = {}
    realReadSql = pandas.read_sql

    def _capturing_read_sql(sql: str, con: Any, params: Any = None, **kwargs: Any) -> pd.DataFrame:
        captured["sql"] = sql
        captured["params"] = params
        return realReadSql(sql, con=con, params=params, **kwargs)

    monkeypatch.setattr(pandas, "read_sql", _capturing_read_sql)

    collection.select_sof_files_to_process(recipe=hostileRecipe, reductionTarget="all")

    assert hostileRecipe not in captured["sql"]
    assert captured["params"]["recipe"] == hostileRecipe


def test_select_sof_files_binds_a_hostile_arm_as_a_parameter_in_all_branch(
    tmp_path: Path,
    log: Any,
) -> None:
    """A hostile `arm` filter matches no real row in a real database, rather than widening the query.

    Unlike `recipe`, `arm` is not re-filtered afterwards in pandas, so a real
    injected `OR '1'='1'` would be directly visible here as extra rows.
    """
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.executemany(
            "INSERT INTO raw_frame_sets VALUES (?, ?, ?, ?, ?)",
            [
                (1, 1, "mbias", "uvb.sof", "UVB"),
                (1, 1, "mbias", "vis.sof", "VIS"),
            ],
        )
    collection = _uninitialized_reducer(log, databasePath)
    # IF STILL INTERPOLATED, THIS CLOSES THE STRING LITERAL, INJECTS AN
    # ALWAYS-TRUE CONDITION, AND COMMENTS OUT THE TRAILING QUOTE, MATCHING
    # BOTH THE UVB AND VIS ROWS INSTEAD OF NEITHER.
    hostileArm = "UVB' OR '1'='1' -- "

    selected = collection.select_sof_files_to_process(reductionTarget="all", arm=hostileArm)

    assert selected.empty


def test_select_sof_files_binds_a_hostile_sof_reduction_target_as_a_parameter(
    tmp_path: Path,
    log: Any,
) -> None:
    """A hostile SOF-filename `reductionTarget` matches no real row, rather than dispatching an unrequested SOF.

    `reductionTarget` is the clearest attacker-reachable path into this
    method: it is a CLI argument (a SOF filename or `all`) forwarded from
    `soxspipe reduce`. This reproduces the recursive-subquery branch
    (`reductionTarget.split(".")[-1].lower() == "sof"`) against a real,
    file-backed database.
    """
    databasePath = tmp_path / "session.db"
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'mbias', 'target.sof', 'UVB')")
        connection.execute("INSERT INTO product_frames (sof, file, complete) VALUES ('target.sof', 'target.fits', 1)")
    collection = _uninitialized_reducer(log, databasePath)
    # BALANCED QUOTES, NO COMMENT MARKER: THE RECURSIVE QUERY-BUILDING LOOP
    # RE-EMBEDS THIS SAME STRING FOUR MORE TIMES, SO AN UNBALANCED PAYLOAD (E.G.
    # ONE RELYING ON `--` TO COMMENT OUT A TRAILING QUOTE) BREAKS EVERY LATER
    # NESTING LEVEL AND JUST RAISES A SYNTAX ERROR RATHER THAN DEMONSTRATING
    # THE INJECTION. IF STILL INTERPOLATED, THIS PAYLOAD MATCHES THE REAL
    # `target.sof` ROW REGARDLESS OF THE `and complete = 1` CONDITION, SO IT IS
    # DISPATCHED EVEN THOUGH IT WAS NEVER REQUESTED.
    hostileReductionTarget = "x' OR sof = 'target.sof' OR sof = 'z.sof"

    selected = collection.select_sof_files_to_process(reductionTarget=hostileReductionTarget)

    assert selected.empty


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
    assert any(level == "error" and "Recipe failed" in message for level, message in log.messages)


@pytest.mark.parametrize(
    ("poolArguments", "expectedPoolSize"),
    [({}, None), ({"poolSize": 2}, 2)],
)
def test_reduce_multiprocess_forwards_selected_group(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    log: Any,
    poolArguments: dict[str, int],
    expectedPoolSize: int | None,
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

        def get_incomplete_sets_report(self) -> tuple[pd.DataFrame, pd.DataFrame]:
            return pd.DataFrame(), pd.DataFrame()

        def close(self) -> None:
            pass

    def fake_bulk(**kwargs: object) -> None:
        bulkCalls.append(kwargs)
        kwargs["conn"].close()

    monkeypatch.setattr(commonutilsModule, "data_organiser", RecordingOrganiser)
    monkeypatch.setattr(reducerModule, "run_recipe_bulk", fake_bulk)

    collection.reduce(batch=1, multiprocess=True, **poolArguments)

    assert len(bulkCalls) == 1
    assert bulkCalls[0]["poolSize"] == expectedPoolSize
    assert bulkCalls[0]["recipe"] == "mbias"
    assert bulkCalls[0]["sofList"] == [str(workspacePath / "sessions" / "night1" / "sof" / "bias.sof")]
    assert bulkCalls[0]["commandList"] == ["run bias"]
    assert bulkCalls[0]["settings"] == {"answer": 42}
    assert bulkCalls[0]["overwrite"] is True
    assert bulkCalls[0]["sessionId"] == "night1"
