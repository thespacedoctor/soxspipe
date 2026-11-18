"""The one place that opens a connection to the workspace database (DY-1706).

Every connection keeps SQLite's durable defaults, `synchronous = FULL` (2) and the rollback journal (`DELETE`), so an
operating-system crash or a power loss cannot corrupt `soxspipe.db`.
"""

from __future__ import annotations

import sqlite3
import threading
from pathlib import Path

import pytest

from soxspipe.commonutils.workspace_database import connect_workspace_database

pytestmark = pytest.mark.unit

SYNCHRONOUS_FULL = 2
BUSY_TIMEOUT_MILLISECONDS = 100000


@pytest.fixture(name="dbPath")
def db_path(tmp_path: Path) -> Path:
    path = tmp_path / "soxspipe.db"
    setup = sqlite3.connect(path)
    setup.execute("create table t (id integer primary key, value text)")
    setup.commit()
    setup.close()
    return path


def _pragma(connection: sqlite3.Connection, name: str) -> object:
    cursor = connection.cursor()
    cursor.row_factory = None
    value = cursor.execute(f"PRAGMA {name}").fetchone()[0]
    cursor.close()
    return value


def test_connection_keeps_synchronous_full_and_the_rollback_journal(dbPath: Path) -> None:
    """A new connection has `synchronous = 2` (FULL) and `journal_mode = delete`."""
    connection = connect_workspace_database(dbPath)

    assert _pragma(connection, "synchronous") == SYNCHRONOUS_FULL
    assert _pragma(connection, "journal_mode") == "delete"
    connection.close()


def test_connection_waits_up_to_the_busy_timeout_for_a_lock(dbPath: Path) -> None:
    """The busy timeout is 100000 ms, whatever the `timeout` argument of `sqlite3.connect` set."""
    connection = connect_workspace_database(dbPath)

    assert _pragma(connection, "busy_timeout") == BUSY_TIMEOUT_MILLISECONDS
    connection.close()


def test_connection_keeps_the_connect_timeout_as_its_busy_timeout_when_asked_not_to_apply_the_pragma(
    dbPath: Path,
) -> None:
    """A caller that must wait only `timeout` seconds on a lock can leave the `busy_timeout` pragma unset."""
    connection = connect_workspace_database(dbPath, timeout=2, applyBusyTimeout=False)

    assert _pragma(connection, "busy_timeout") == 2000
    connection.close()


def test_connection_commits_each_statement_at_once(dbPath: Path) -> None:
    """A statement is committed when it runs, so a second connection reads it with no `commit`."""
    writer = connect_workspace_database(dbPath)
    reader = sqlite3.connect(dbPath)

    writer.execute("insert into t (value) values ('a')")

    assert reader.execute("select count(*) from t").fetchone()[0] == 1
    writer.close()
    reader.close()


def test_connection_can_be_used_from_another_thread(dbPath: Path) -> None:
    """The connection is shared with worker threads, so SQLite's same-thread check is off."""
    connection = connect_workspace_database(dbPath)
    errors: list[Exception] = []

    def use_connection() -> None:
        try:
            connection.execute("select count(*) from t").fetchone()
        except sqlite3.ProgrammingError as error:
            errors.append(error)

    worker = threading.Thread(target=use_connection)
    worker.start()
    worker.join()

    assert errors == []
    connection.close()


def test_connection_applies_the_row_factory_it_is_given(dbPath: Path) -> None:
    """With a `rowFactory`, rows come back as that factory builds them."""
    connection = connect_workspace_database(dbPath, rowFactory=sqlite3.Row)

    row = connection.execute("select 1 as one").fetchone()

    assert row["one"] == 1
    connection.close()
