"""Every connection to the workspace database keeps SQLite's durable defaults (DY-1706).

`synchronous` stays `FULL` and the journal stays the rollback journal (`DELETE`), so an operating-system crash or a
power loss cannot corrupt `soxspipe.db`.
"""

from __future__ import annotations

import sqlite3
from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

from soxspipe.recipes import base_recipe
from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration

SYNCHRONOUS_FULL = 2
SESSION_ID = "base"


def _connection_settings(connection: sqlite3.Connection) -> tuple[int, str]:
    """Return the `synchronous` level and the journal mode, whatever the row factory of the connection is."""
    cursor = connection.cursor()
    cursor.row_factory = None
    synchronous = cursor.execute("PRAGMA synchronous").fetchone()[0]
    journalMode = cursor.execute("PRAGMA journal_mode").fetchone()[0]
    cursor.close()
    return synchronous, journalMode


def test_data_organiser_connection_uses_full_synchronous_and_the_rollback_journal(tmp_path: Path, log: Any) -> None:
    """The organiser's connection, as `prepare` and `session_create` open it, is `FULL` with a `DELETE` journal."""
    organiser = workspace_organiser(tmp_path, log=log)

    connection, _ = organiser._get_or_create_db_connection()

    assert _connection_settings(connection) == (SYNCHRONOUS_FULL, "delete")
    connection.close()


def test_recipe_session_database_connection_uses_full_synchronous_and_the_rollback_journal(
    tmp_path: Path, log: Any
) -> None:
    """The connection a recipe opens to record its status is `FULL` with a `DELETE` journal."""
    databasePath = tmp_path / "soxspipe.db"
    setup = sqlite3.connect(databasePath)
    setup.execute(f"create table product_frames (sof text, status_{SESSION_ID} text)")
    setup.execute("insert into product_frames values ('synthetic.sof', 'pass')")
    setup.commit()
    setup.close()
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.currentSession = SESSION_ID
    recipe.sofName = "synthetic"
    recipe.settings = {"workspace-root-dir": str(tmp_path)}
    recipe.conn = None

    recipe._open_session_database(home=str(tmp_path))

    assert _connection_settings(recipe.conn) == (SYNCHRONOUS_FULL, "delete")
    recipe.conn.close()


def test_reducer_selection_connection_does_not_set_synchronous(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Selecting SOFs to process opens a connection that never issues a `synchronous` pragma."""
    databasePath = tmp_path / "soxspipe.db"
    setup = sqlite3.connect(databasePath)
    setup.execute(
        "CREATE TABLE raw_frame_sets "
        '(recipe_order INTEGER, complete INTEGER, recipe TEXT, sof TEXT, "eso seq arm" TEXT)'
    )
    setup.execute(f"CREATE TABLE product_frames (sof TEXT, file TEXT, complete INTEGER, status_{SESSION_ID} TEXT)")
    setup.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'mbias', 'bias.sof', 'UVB')")
    setup.commit()
    setup.close()
    reducerModule = import_module("soxspipe.commonutils.reducer")
    collection = reducerModule.reducer.__new__(reducerModule.reducer)
    collection.log = log
    collection.sessionDB = str(databasePath)
    collection.workspaceDirectory = str(tmp_path)
    collection.pathToSettings = False
    collection.sessionId = SESSION_ID
    collection.overwrite = True
    statements: list[str] = []
    realConnect = sqlite3.connect

    def recording_connect(*args: Any, **kwargs: Any) -> sqlite3.Connection:
        connection = realConnect(*args, **kwargs)
        connection.set_trace_callback(statements.append)
        return connection

    monkeypatch.setattr(sqlite3, "connect", recording_connect)

    selected = collection.select_sof_files_to_process(reductionTarget="all")

    assert selected["sof"].tolist() == ["bias.sof"]
    assert [statement for statement in statements if "synchronous" in statement.lower()] == []
