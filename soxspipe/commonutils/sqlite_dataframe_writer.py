#!/usr/bin/env python
"""
*Write a dataframe to a SQLite table in one transaction, retrying a failed write*

Author
: David Young

Date Created
: October 19, 2026
"""

# NO IMPORTS FROM OTHER SOXSPIPE MODULES BEYOND THE LEAF `sql_identifiers`: THE
# PACKAGE `__init__` IMPORT ORDER IS LOAD-BEARING, AND THIS MODULE MUST STAY
# SAFE TO IMPORT FROM ANYWHERE.
import sqlite3
import time
from contextlib import closing, contextmanager

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

MAX_ATTEMPTS = 7
RETRY_DELAY_SECONDS = 1


class SqliteTransactionOpenError(RuntimeError):
    """*raised when a write starts on a connection that is already inside a transaction*"""


@contextmanager
def _write_transaction(conn):
    """*hold the database write lock for the length of one `with` block, then commit or roll back*

    `BEGIN IMMEDIATE` is an explicit statement, so it works on an `autocommit=True` connection (where each statement
    would otherwise commit on its own) and on a legacy connection alike. It also takes the write lock up front, so
    the busy timeout applies at once and a read lock is never upgraded part-way through. pandas calls
    `conn.commit()` when `to_sql` finishes: that is a no-op on an `autocommit=True` connection, and on a legacy one
    it commits this transaction early, so `conn.in_transaction` is checked before the final `COMMIT`.

    Any failure in the block or in the final `COMMIT` rolls back this transaction and propagates. A `BaseException`
    (such as `KeyboardInterrupt`) is never retried by the caller.

    **Key Arguments:**

    - ``conn`` -- the `sqlite3` connection to write on

    **Raises:**

    - ``SqliteTransactionOpenError`` if the connection is already inside a transaction. pandas would commit or roll
      back the caller's unrelated work as part of `to_sql`, so the write is refused before anything is written.
    """
    if conn.in_transaction:
        raise SqliteTransactionOpenError(
            "the connection is already inside a transaction, so the write would commit or roll back the caller's work"
        )
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield
        # PANDAS MAY ALREADY HAVE COMMITTED ON A LEGACY CONNECTION, SO COMMIT ONLY AN OPEN TRANSACTION
        if conn.in_transaction:
            conn.execute("COMMIT")
    except BaseException:
        # A FAILED COMMIT (BUSY, INTERRUPT) ALSO LANDS HERE, SO THE WRITE LOCK IS NEVER LEFT HELD
        if conn.in_transaction:
            conn.execute("ROLLBACK")
        raise


def _delete_rows(conn, log, tableName, whereColumn=None, whereValues=()):
    """*delete every row of a table, or the rows whose column value is in a list, treating a missing table as empty*"""
    # THE TABLE AND COLUMN NAMES HAVE PASSED `validate_sql_identifier` (A NAME CANNOT BE A BOUND PARAMETER);
    # EVERY VALUE IS BOUND.
    sqlQuery = f"delete from {tableName}"  # noqa: S608
    if whereColumn:
        sqlQuery += f" where {whereColumn} in ({', '.join(['?'] * len(whereValues))})"
    with closing(conn.cursor()) as c:
        try:
            c.execute(sqlQuery + ";", list(whereValues))
        except sqlite3.OperationalError as e:
            # THE QUERY NAMES THE TABLE UNQUALIFIED, SO SQLITE REPORTS EXACTLY `no such table: <name>`
            if str(e) != f"no such table: {tableName}":
                raise
            log.debug(f"_delete_rows: `{tableName}` does not exist yet, nothing to delete: {e}")


def write_dataframe_to_sqlite(
    conn,
    log,
    dataframe,
    tableName,
    replace=False,
    sentinels=("--",),
    method=None,
    maxAttempts=MAX_ATTEMPTS,
    deleteWhereIn=None,
):
    """*write a dataframe to a database table in one transaction, retrying a failed write*

    Each attempt inserts all the rows or none. With ``replace=True`` the delete shares the transaction, so a reader
    never sees the table empty and a failed attempt gives the old rows back before the next attempt starts.

    **Key Arguments:**

    - ``conn`` -- the `sqlite3` connection to write on (`autocommit=True` or legacy transaction control)
    - ``log`` -- logger
    - ``dataframe`` -- the dataframe containing the rows to insert
    - ``tableName`` -- the name of the database table to insert into
    - ``replace`` -- if True, delete the table's existing rows first; otherwise append. Default: False
    - ``sentinels`` -- values in the dataframe that are written as NULL. Default: ("--",)
    - ``method`` -- the `DataFrame.to_sql` insertion method. Default: None (one `executemany`)
    - ``maxAttempts`` -- how many times to try before re-raising the last error. Default: `MAX_ATTEMPTS`
    - ``deleteWhereIn`` -- a ``(column, values)`` pair: delete the rows whose column value is in ``values``, in the same
      transaction, before the insert. The column name is validated and the values are bound. Cannot be combined
      with ``replace``. Default: None

    **Raises:**

    - The error from the final attempt if every attempt fails.
    - ``SqliteTransactionOpenError`` at once, without retrying, if the connection is already inside a transaction.
    - ``ValueError`` if both ``replace`` and ``deleteWhereIn`` are given.
    - ``UnsafeSqlIdentifierError`` if ``tableName`` fails the safe-identifier grammar.
    """
    # A TABLE NAME CANNOT BE A BOUND PARAMETER, SO IT IS VALIDATED AGAINST
    # THE SAFE-IDENTIFIER GRAMMAR ONCE HERE, BEFORE ANY SQL INTERPOLATES IT.
    tableName = validate_sql_identifier(tableName, "table name")
    if replace and deleteWhereIn:
        raise ValueError("`replace` and `deleteWhereIn` cannot be combined: replace already deletes every row")
    whereColumn, whereValues = None, ()
    if deleteWhereIn:
        whereColumn = validate_sql_identifier(deleteWhereIn[0], "delete column name")
        whereValues = list(dict.fromkeys(deleteWhereIn[1]))

    for attempt in range(1, maxAttempts + 1):
        try:
            with _write_transaction(conn):
                if replace or whereColumn:
                    _delete_rows(conn, log, tableName, whereColumn, whereValues)
                dataframe.replace(list(sentinels), None).to_sql(
                    tableName, con=conn, index=False, if_exists="append", method=method
                )
            return
        except SqliteTransactionOpenError:
            raise
        except Exception:
            # RE-RAISE THE ORIGINAL ERROR AFTER THE FINAL ATTEMPT SO A LOST WRITE FAILS THE CALLER (DY-603)
            if attempt == maxAttempts:
                raise
            time.sleep(RETRY_DELAY_SECONDS)
