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


@contextmanager
def _write_transaction(conn):
    """*hold the database write lock for the length of one `with` block, then commit or roll back*

    `BEGIN IMMEDIATE` is an explicit statement, so it works on an `autocommit=True` connection (where each statement
    would otherwise commit on its own) and on a legacy connection alike. It also takes the write lock up front, so
    the busy timeout applies at once and a read lock is never upgraded part-way through. pandas calls
    `conn.commit()` when `to_sql` finishes: that is a no-op on an `autocommit=True` connection, and on a legacy one
    it commits this transaction early, so `conn.in_transaction` is checked before the final `COMMIT`.

    A connection that is already inside a transaction belongs to its caller, who keeps control of the commit.

    **Key Arguments:**

    - ``conn`` -- the `sqlite3` connection to write on
    """
    if conn.in_transaction:
        yield
        return
    conn.execute("BEGIN IMMEDIATE")
    try:
        yield
    except BaseException:
        if conn.in_transaction:
            conn.execute("ROLLBACK")
        raise
    if conn.in_transaction:
        conn.execute("COMMIT")


def _delete_existing_rows(conn, log, tableName):
    """*delete every row of a table, treating a table that does not exist yet as already empty*"""
    # THE TABLE NAME HAS PASSED `validate_sql_identifier`; A NAME CANNOT BE A BOUND PARAMETER.
    sqlQuery = f"delete from {tableName};"  # noqa: S608
    with closing(conn.cursor()) as c:
        try:
            c.execute(sqlQuery)
        except sqlite3.OperationalError as e:
            if "no such table" not in str(e):
                raise
            log.debug(f"_delete_existing_rows: `{tableName}` does not exist yet, nothing to delete: {e}")


def write_dataframe_to_sqlite(
    conn,
    log,
    dataframe,
    tableName,
    replace=False,
    sentinels=("--",),
    method=None,
    maxAttempts=MAX_ATTEMPTS,
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
    - ``maxAttempts`` -- how many times to try before re-raising the last error. Default: 7

    **Raises:**

    - The error from the final attempt if every attempt fails.
    - ``UnsafeSqlIdentifierError`` if ``tableName`` fails the safe-identifier grammar.
    """
    # A TABLE NAME CANNOT BE A BOUND PARAMETER, SO IT IS VALIDATED AGAINST
    # THE SAFE-IDENTIFIER GRAMMAR ONCE HERE, BEFORE ANY SQL INTERPOLATES IT.
    tableName = validate_sql_identifier(tableName, "table name")

    for attempt in range(1, maxAttempts + 1):
        try:
            with _write_transaction(conn):
                if replace:
                    _delete_existing_rows(conn, log, tableName)
                dataframe.replace(list(sentinels), None).to_sql(
                    tableName, con=conn, index=False, if_exists="append", method=method
                )
            return
        except Exception:
            # RE-RAISE THE ORIGINAL ERROR AFTER THE FINAL ATTEMPT SO A LOST WRITE FAILS THE CALLER (DY-603)
            if attempt == maxAttempts:
                raise
            time.sleep(RETRY_DELAY_SECONDS)
