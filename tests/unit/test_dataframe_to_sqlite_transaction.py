"""Transaction contract for both `_dataframe_to_sqlite` copies (DY-1515).

Workspace connections open with `autocommit=True`, so a bare `executemany` commits every row on its own.
Each `to_sql` call must instead insert all its rows or none, on an autocommit connection and on a legacy one,
and a `replace=True` delete must live in the same transaction as the insert that follows it.
"""

from __future__ import annotations

import sqlite3
import time
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from soxspipe.commonutils.sql_identifiers import UnsafeSqlIdentifierError
from soxspipe.commonutils.sqlite_dataframe_writer import SqliteTransactionOpenError, _delete_rows
from soxspipe.recipes.base_recipe import base_recipe

pytestmark = pytest.mark.unit

ROW_COUNT = 6
ROWS_WRITTEN_BEFORE_FAILURE = 3
MAX_ATTEMPTS = 7

CONNECTION_MODES = ["autocommit", "legacy"]
WRITER_CLASSES = [base_recipe, data_organiser]


def _connect(databasePath: Path, mode: str, timeout: float = 5.0) -> sqlite3.Connection:
    """Open a connection the way the pipeline does (`autocommit`) or in the legacy transaction mode."""
    if mode == "autocommit":
        return sqlite3.connect(databasePath, timeout=timeout, autocommit=True)
    return sqlite3.connect(databasePath, timeout=timeout)


def _frame(prefix: str = "new", rowCount: int = ROW_COUNT) -> pd.DataFrame:
    """Build a dataframe of `rowCount` rows for the `rows` table."""
    return pd.DataFrame({"id": range(rowCount), "label": [f"{prefix}{i}" for i in range(rowCount)]})


def _labels(connection: sqlite3.Connection) -> list[str]:
    """Return the sorted labels currently committed in the `rows` table."""
    return sorted(row[0] for row in connection.execute("select label from rows"))


@pytest.fixture(name="databasePath")
def database_path_fixture(tmp_path: Path) -> Path:
    """Create a file database holding an empty `rows` table."""
    path = tmp_path / "workspace.db"
    setup = sqlite3.connect(path)
    setup.execute("create table rows (id integer, label text)")
    setup.commit()
    setup.close()
    return path


@pytest.fixture(params=CONNECTION_MODES)
def connection(request: pytest.FixtureRequest, databasePath: Path) -> Iterator[sqlite3.Connection]:
    """Yield the connection under test in each connection mode."""
    conn = _connect(databasePath, request.param)
    yield conn
    conn.close()


@pytest.fixture(params=WRITER_CLASSES, ids=lambda cls: cls.__name__)
def write(request: pytest.FixtureRequest, log: Any) -> Callable[..., None]:
    """Return a writer that calls the `_dataframe_to_sqlite` of each class on a given connection."""

    def _write(conn: sqlite3.Connection, dataframe: pd.DataFrame, replace: bool = False) -> None:
        owner = request.param.__new__(request.param)
        owner.log = log
        owner.conn = conn
        owner._dataframe_to_sqlite(dataframe, "rows", replace=replace)

    return _write


@pytest.fixture
def sleep_calls(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Record every `time.sleep` call instead of waiting."""
    calls: list[float] = []
    monkeypatch.setattr(time, "sleep", calls.append)
    return calls


@dataclass
class InsertFaults:
    """How many upcoming `to_sql` calls fail part-way, with what, and what to run just before they fail."""

    failures: int = 0
    attempts: int = 0
    error: BaseException = field(default_factory=lambda: sqlite3.OperationalError("database is locked"))
    probe: Callable[[], None] | None = None


@pytest.fixture
def failing_inserts(monkeypatch: pytest.MonkeyPatch) -> InsertFaults:
    """Make the next `failures` `to_sql` calls write part of their rows and then fail, as a lock timeout would."""
    faults = InsertFaults()
    originalToSql = pd.DataFrame.to_sql

    def _insert_part_then_fail(table: Any, cursor: Any, keys: Any, dataIter: Any) -> None:
        rows = list(dataIter)
        table._execute_insert(cursor, keys, iter(rows[:ROWS_WRITTEN_BEFORE_FAILURE]))
        if faults.probe:
            faults.probe()
        raise faults.error

    def flaky_to_sql(self: pd.DataFrame, *args: Any, **kwargs: Any) -> Any:
        faults.attempts += 1
        if faults.failures > 0:
            faults.failures -= 1
            kwargs["method"] = _insert_part_then_fail
        return originalToSql(self, *args, **kwargs)

    monkeypatch.setattr(pd.DataFrame, "to_sql", flaky_to_sql)
    return faults


def test_failed_insert_then_retry_leaves_each_row_exactly_once(connection, write, sleep_calls, failing_inserts) -> None:
    # ARRANGE
    failing_inserts.failures = 2

    # ACT
    write(connection, _frame())

    # ASSERT
    assert failing_inserts.attempts == 3
    assert _labels(connection) == sorted(f"new{i}" for i in range(ROW_COUNT))


def test_insert_that_fails_every_attempt_leaves_no_rows(connection, write, sleep_calls, failing_inserts) -> None:
    # ARRANGE
    failing_inserts.failures = MAX_ATTEMPTS

    # ACT
    with pytest.raises(sqlite3.OperationalError, match="database is locked"):
        write(connection, _frame())

    # ASSERT
    assert _labels(connection) == []
    assert connection.in_transaction is False


def test_successful_insert_is_committed_and_leaves_no_open_transaction(connection, write, databasePath) -> None:
    # ARRANGE
    reader = sqlite3.connect(databasePath)

    # ACT
    write(connection, _frame())

    # ASSERT
    assert connection.in_transaction is False
    assert len(_labels(reader)) == ROW_COUNT
    reader.close()


def test_replace_keeps_the_old_rows_when_the_insert_fails(connection, write, sleep_calls, failing_inserts) -> None:
    # ARRANGE
    write(connection, _frame("old", 2))
    failing_inserts.failures = MAX_ATTEMPTS

    # ACT
    with pytest.raises(sqlite3.OperationalError):
        write(connection, _frame(), replace=True)

    # ASSERT
    assert _labels(connection) == ["old0", "old1"]


def test_replace_that_succeeds_after_a_failed_attempt_holds_only_the_new_rows(
    connection, write, sleep_calls, failing_inserts
) -> None:
    # ARRANGE
    write(connection, _frame("old", 2))
    failing_inserts.failures = 1

    # ACT
    write(connection, _frame(), replace=True)

    # ASSERT
    assert _labels(connection) == sorted(f"new{i}" for i in range(ROW_COUNT))


def test_reader_never_sees_an_empty_table_while_replace_is_in_flight(
    connection, write, databasePath, sleep_calls, failing_inserts
) -> None:
    # ARRANGE
    write(connection, _frame("old", 2))
    reader = sqlite3.connect(databasePath, timeout=0)
    seenByReader: list[list[str]] = []
    failing_inserts.failures = 1
    # THE PROBE RUNS AFTER THE DELETE AND PART OF THE INSERT, BEFORE THE FAILURE
    failing_inserts.probe = lambda: seenByReader.append(_labels(reader))

    # ACT
    write(connection, _frame(), replace=True)

    # ASSERT
    assert seenByReader == [["old0", "old1"]]
    reader.close()


def test_insert_retries_after_another_connection_holds_the_write_lock(write, databasePath, monkeypatch) -> None:
    # ARRANGE
    lockHolder = sqlite3.connect(databasePath, autocommit=True)
    lockHolder.execute("BEGIN IMMEDIATE")
    impatient = sqlite3.connect(databasePath, autocommit=True, timeout=0)
    sleeps: list[float] = []

    def release_lock_on_first_sleep(seconds: float) -> None:
        sleeps.append(seconds)
        lockHolder.execute("COMMIT")

    monkeypatch.setattr(time, "sleep", release_lock_on_first_sleep)

    # ACT
    write(impatient, _frame())

    # ASSERT
    assert len(sleeps) == 1
    assert len(_labels(impatient)) == ROW_COUNT
    impatient.close()
    lockHolder.close()


def test_replace_on_a_table_that_does_not_exist_yet_creates_it(connection, log) -> None:
    # ARRANGE
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = connection

    # ACT
    owner._dataframe_to_sqlite(_frame(), "fresh_rows", replace=True)

    # ASSERT
    assert connection.execute("select count(*) from fresh_rows").fetchone() == (ROW_COUNT,)


def test_replace_does_not_swallow_a_delete_error_other_than_a_missing_table(
    connection, databasePath, log, sleep_calls
) -> None:
    # ARRANGE
    connection.execute("create view rows_view as select * from rows")
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = connection

    # ACT
    with pytest.raises(sqlite3.OperationalError, match="view"):
        owner._dataframe_to_sqlite(_frame(), "rows_view", replace=True)

    # ASSERT
    assert len(sleep_calls) == MAX_ATTEMPTS - 1
    assert _labels(connection) == []


def test_insert_inside_a_caller_owned_transaction_is_refused_before_anything_is_written(
    connection, write, databasePath, sleep_calls
) -> None:
    # ARRANGE
    connection.execute("BEGIN")
    connection.execute("insert into rows values (99, 'caller')")
    reader = sqlite3.connect(databasePath)

    # ACT
    with pytest.raises(SqliteTransactionOpenError, match="already inside a transaction"):
        write(connection, _frame())

    # ASSERT
    assert sleep_calls == []
    assert connection.in_transaction is True
    assert _labels(reader) == []
    connection.execute("COMMIT")
    assert _labels(reader) == ["caller"]
    reader.close()


def test_keyboard_interrupt_mid_write_rolls_back_and_is_not_retried(
    connection, write, sleep_calls, failing_inserts
) -> None:
    # ARRANGE
    write(connection, _frame("old", 2), replace=True)
    failing_inserts.error = KeyboardInterrupt()
    failing_inserts.failures = 1
    failing_inserts.attempts = 0

    # ACT
    with pytest.raises(KeyboardInterrupt):
        write(connection, _frame(), replace=True)

    # ASSERT
    assert failing_inserts.attempts == 1
    assert sleep_calls == []
    assert connection.in_transaction is False
    assert _labels(connection) == ["old0", "old1"]


def _qc_frame(sofName: str, value: str) -> pd.DataFrame:
    """Build two QC rows for one SOF."""
    return pd.DataFrame({"sof_name": [sofName, sofName], "qc_value": [value, value]})


@pytest.fixture
def qc_connection(connection) -> sqlite3.Connection:
    """Give the connection a `quality_control` table holding rows for two SOFs."""
    connection.execute("create table quality_control (sof_name text, qc_value text)")
    connection.executemany(
        "insert into quality_control values (?, ?)", [("a.sof", "old"), ("a.sof", "old"), ("b.sof", "keep")]
    )
    connection.commit()
    return connection


def _qc_rows(connection: sqlite3.Connection) -> list[tuple[str, str]]:
    """Return the sorted `(sof_name, qc_value)` rows of `quality_control`."""
    return sorted(connection.execute("select sof_name, qc_value from quality_control"))


def test_delete_where_in_replaces_only_the_named_rows_when_the_insert_succeeds(qc_connection, log) -> None:
    # ARRANGE
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = qc_connection

    # ACT
    owner._dataframe_to_sqlite(_qc_frame("a.sof", "new"), "quality_control", delete_where_in=("sof_name", ["a.sof"]))

    # ASSERT
    assert _qc_rows(qc_connection) == [("a.sof", "new"), ("a.sof", "new"), ("b.sof", "keep")]


def test_delete_where_in_keeps_the_old_rows_when_every_insert_attempt_fails(
    qc_connection, log, sleep_calls, failing_inserts
) -> None:
    # ARRANGE
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = qc_connection
    failing_inserts.failures = MAX_ATTEMPTS

    # ACT
    with pytest.raises(sqlite3.OperationalError):
        owner._dataframe_to_sqlite(
            _qc_frame("a.sof", "new"), "quality_control", delete_where_in=("sof_name", ["a.sof"])
        )

    # ASSERT
    assert _qc_rows(qc_connection) == [("a.sof", "old"), ("a.sof", "old"), ("b.sof", "keep")]


def test_delete_where_in_rejects_a_hostile_column_name(qc_connection, log) -> None:
    # ARRANGE
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = qc_connection

    # ACT / ASSERT
    with pytest.raises(UnsafeSqlIdentifierError):
        owner._dataframe_to_sqlite(
            _qc_frame("a.sof", "new"), "quality_control", delete_where_in=("sof_name; drop table quality_control", [])
        )
    assert len(_qc_rows(qc_connection)) == 3


def test_delete_that_names_a_different_missing_table_is_not_swallowed(connection, log, sleep_calls) -> None:
    # ARRANGE
    connection.execute("create trigger drop_guard before delete on rows begin select * from vanished; end")
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = connection

    # ACT / ASSERT
    with pytest.raises(sqlite3.OperationalError, match="vanished"):
        owner._dataframe_to_sqlite(_frame(), "rows", replace=True)


class CommitFaultConnection(sqlite3.Connection):
    """A connection whose next `COMMIT` statements raise the queued errors, once each."""

    commitErrors: list[BaseException]

    def execute(self, sql: str, *args: Any) -> sqlite3.Cursor:
        if sql == "COMMIT" and self.commitErrors:
            raise self.commitErrors.pop(0)
        return super().execute(sql, *args)


@pytest.fixture
def commit_fault_connection(databasePath: Path) -> Iterator[CommitFaultConnection]:
    """Yield an autocommit connection whose `COMMIT` can be made to fail."""
    conn = sqlite3.connect(databasePath, autocommit=True, factory=CommitFaultConnection)
    conn.commitErrors = []
    yield conn
    conn.close()


def test_commit_that_fails_once_is_rolled_back_and_the_retry_writes_the_rows_once(
    commit_fault_connection, write, databasePath, sleep_calls
) -> None:
    # ARRANGE
    write(commit_fault_connection, _frame("old", 2))
    commit_fault_connection.commitErrors = [sqlite3.OperationalError("database is locked")]

    # ACT
    write(commit_fault_connection, _frame(), replace=True)

    # ASSERT
    assert len(sleep_calls) == 1
    assert commit_fault_connection.in_transaction is False
    assert _labels(commit_fault_connection) == sorted(f"new{i}" for i in range(ROW_COUNT))


def test_keyboard_interrupt_during_commit_rolls_back_and_is_not_retried(
    commit_fault_connection, write, sleep_calls
) -> None:
    # ARRANGE
    write(commit_fault_connection, _frame("old", 2))
    commit_fault_connection.commitErrors = [KeyboardInterrupt()]

    # ACT
    with pytest.raises(KeyboardInterrupt):
        write(commit_fault_connection, _frame(), replace=True)

    # ASSERT
    assert sleep_calls == []
    assert commit_fault_connection.in_transaction is False
    assert _labels(commit_fault_connection) == ["old0", "old1"]


class _MissingTableConnection:
    """A stand-in connection whose every statement fails with a fixed `no such table` message."""

    def __init__(self, message: str) -> None:
        self.message = message

    def cursor(self) -> _MissingTableConnection:
        return self

    def execute(self, sql: str, params: Any = ()) -> None:
        raise sqlite3.OperationalError(self.message)

    def close(self) -> None:
        return None


def test_delete_that_names_a_table_sharing_the_requested_name_prefix_is_not_swallowed(log) -> None:
    # ARRANGE
    stub = _MissingTableConnection("no such table: quality_control_old")

    # ACT / ASSERT
    with pytest.raises(sqlite3.OperationalError, match="quality_control_old"):
        _delete_rows(stub, log, "quality_control")


def test_delete_of_the_requested_missing_table_is_treated_as_already_empty(log) -> None:
    # ARRANGE
    stub = _MissingTableConnection("no such table: quality_control")

    # ACT / ASSERT
    _delete_rows(stub, log, "quality_control")


def test_replace_combined_with_delete_where_in_is_rejected(connection, log, sleep_calls) -> None:
    # ARRANGE
    owner = base_recipe.__new__(base_recipe)
    owner.log = log
    owner.conn = connection

    # ACT / ASSERT
    with pytest.raises(ValueError, match="cannot be combined"):
        owner._dataframe_to_sqlite(_frame(), "rows", replace=True, delete_where_in=("label", ["x"]))
    assert sleep_calls == []
