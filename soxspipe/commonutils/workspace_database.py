"""Open a connection to the workspace database `soxspipe.db` (DY-1706).

Every read-write connection to the workspace database starts here, so the durability settings live in one place.
The read-only connection that the backup opens is the one exception. The connection keeps SQLite's default
`synchronous = FULL` and the rollback journal (`DELETE`): an operating-system crash or a power loss cannot corrupt
the database. Do not weaken the `synchronous` setting or switch to WAL here. WAL is
unsafe on network volumes, and a workspace can sit on one.
"""

# NO IMPORTS FROM OTHER SOXSPIPE MODULES: THE PACKAGE `__init__` IMPORT ORDER
# IS LOAD-BEARING, AND THIS MODULE MUST STAY SAFE TO IMPORT FROM ANYWHERE.
import sqlite3

# SECONDS A STATEMENT WAITS FOR ANOTHER CONNECTION'S LOCK, AS SET BY `sqlite3.connect`
DEFAULT_CONNECT_TIMEOUT_SECONDS = 300
# MILLISECONDS OF THE `busy_timeout` PRAGMA, WHICH REPLACES THE `timeout` OF `sqlite3.connect` ONCE SET
BUSY_TIMEOUT_MILLISECONDS = 100000


def connect_workspace_database(
    databasePath, timeout=DEFAULT_CONNECT_TIMEOUT_SECONDS, rowFactory=None, applyBusyTimeout=True
):
    """*open an autocommit connection to the workspace database that worker threads can share*

    **Key Arguments:**

    - ``databasePath`` -- path to the `soxspipe.db` file
    - ``timeout`` -- seconds to wait for a lock when connecting. Default *300*
    - ``rowFactory`` -- the `row_factory` for the connection, or *None* for tuples. Default *None*
    - ``applyBusyTimeout`` -- set the `busy_timeout` pragma to `BUSY_TIMEOUT_MILLISECONDS`. Pass *False* to keep the
      `timeout` seconds as the lock wait, for a caller that sets the pragma itself once its first query has run.
      Default *True*

    **Return:**

    - ``conn`` -- the open connection, with `synchronous = FULL` and the `DELETE` journal

    **Usage:**

    ```python
    from soxspipe.commonutils.workspace_database import connect_workspace_database

    conn = connect_workspace_database("/path/to/workspace/soxspipe.db")
    ```
    """
    conn = sqlite3.connect(
        databasePath,
        timeout=timeout,
        autocommit=True,
        check_same_thread=False,
    )
    if applyBusyTimeout:
        conn.execute(f"PRAGMA busy_timeout = {BUSY_TIMEOUT_MILLISECONDS}")
    if rowFactory is not None:
        conn.row_factory = rowFactory
    return conn
