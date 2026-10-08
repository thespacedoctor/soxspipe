#!/usr/bin/env python
"""
*Read the per-session pass/fail statuses out of a preserved workspace database and decide which can be restored*

Author
: David Young

Date Created
: September 30, 2026
"""

from __future__ import annotations

import logging
import shutil
import sqlite3
import tempfile
from collections.abc import Iterable, Iterator, Mapping
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

# NO IMPORTS FROM OTHER SOXSPIPE MODULES EXCEPT sql_identifiers: THE PACKAGE `__init__` IMPORT ORDER IS
# LOAD-BEARING, AND THIS MODULE MUST STAY SAFE TO IMPORT FROM ANYWHERE.

# THE ONLY STATUS VALUES THAT ARE WORTH CARRYING ACROSS A REBUILD (NULL MEANS "NOT REDUCED YET")
RESTORABLE_STATUSES = frozenset({"pass", "fail"})

# THE SQLITE FILES THAT MAKE UP ONE DATABASE: THE MAIN FILE AND ITS WRITE-AHEAD LOG AND SHARED-MEMORY SIDECARS
DATABASE_FILE_SUFFIXES = ("", "-wal", "-shm")

# NAME PREFIX OF THE TEMPORARY DIRECTORY A PRESERVED DATABASE IS COPIED TO FOR READING
SCRATCH_DIRECTORY_PREFIX = "soxspipe-backup-read-"

# THE LONGEST SLICE OF A VALUE FROM THE (UNTRUSTED) BACKUP THAT MAY REACH A LOG LINE OR AN EXCEPTION MESSAGE
MAX_LOGGED_VALUE_LENGTH = 80


@dataclass(frozen=True)
class session_status_snapshot:
    """The restorable statuses and frame sets of one session, as read from a preserved database.

    ``statuses`` maps each SOF name to ``pass`` or ``fail`` (NULL statuses are left out). ``frameSets`` maps each
    SOF name to the frozenset of `file` names it held, or is None when the backup has no `sof_map_<id>` table.
    """

    sessionId: str
    statuses: Mapping[str, str]
    frameSets: Mapping[str, frozenset[str]] | None

    def __post_init__(self) -> None:
        object.__setattr__(self, "statuses", MappingProxyType(dict(self.statuses)))
        if self.frameSets is not None:
            frozenSets = {sof: frozenset(files) for sof, files in self.frameSets.items()}
            object.__setattr__(self, "frameSets", MappingProxyType(frozenSets))


@dataclass(frozen=True)
class status_pass_result:
    """The outcome of classifying one session's candidate SOFs against the rebuilt database."""

    toRestore: Mapping[str, str]
    pending: frozenset[str]
    changed: frozenset[str]
    dropped: frozenset[str]

    def __post_init__(self) -> None:
        object.__setattr__(self, "toRestore", MappingProxyType(dict(self.toRestore)))


def bounded_repr(value: object) -> str:
    """*return the `repr` of a value from the backup, cut to `MAX_LOGGED_VALUE_LENGTH` characters*

    **Key Arguments:**

    - ``value`` -- the untrusted value to make safe to log

    **Return:**

    - ``text`` -- the `repr` of ``value``, with an ellipsis when it was cut
    """
    text = repr(value)
    if len(text) > MAX_LOGGED_VALUE_LENGTH:
        return text[:MAX_LOGGED_VALUE_LENGTH] + "..."
    return text


def format_session_summary(sessionId: str, restored: int, changed: int, dropped: int, backupPath: Path | str) -> str:
    """*return the one-line summary of what was restored for a session*

    **Key Arguments:**

    - ``sessionId`` -- the session the counts belong to
    - ``restored`` -- the number of statuses written back
    - ``changed`` -- the number of SOFs whose frames differ, so their status was not restored
    - ``dropped`` -- the number of SOFs whose product no longer exists
    - ``backupPath`` -- the preserved database, named only when some status was not restored

    **Return:**

    - ``summary`` -- the summary line
    """
    summary = (
        f"session '{sessionId}': {restored} statuses restored, {changed} changed (frames differ, requeued), "
        f"{dropped} dropped (product no longer exists)"
    )
    if changed + dropped > 0:
        summary += f"; the previous statuses are kept in {backupPath}"
    return summary


def classify_session_sofs(
    snapshot: session_status_snapshot,
    candidateSofs: Iterable[str],
    rebuiltSofs: Iterable[str],
    rebuiltFrameSets: Mapping[str, frozenset[str]],
) -> status_pass_result:
    """*sort candidate SOFs into restorable, pending, changed and dropped*

    A SOF is restorable when it holds exactly the frames it held in the backup. A SOF the rebuilt database no
    longer has is dropped. A SOF is changed when the backup has no frame sets to compare with, and pending when
    its frames differ now but may still match after a later pass re-matches its calibrations.

    **Key Arguments:**

    - ``snapshot`` -- the `session_status_snapshot` read from the backup
    - ``candidateSofs`` -- the SOF names to classify
    - ``rebuiltSofs`` -- the SOF names present in the rebuilt `product_frames`
    - ``rebuiltFrameSets`` -- SOF name to the frozenset of `file` names in the rebuilt `sof_map_<id>`

    **Return:**

    - ``result`` -- a `status_pass_result`
    """
    toRestore = {}
    pending = set()
    changed = set()
    dropped = set()
    for sof in candidateSofs:
        if sof not in rebuiltSofs:
            dropped.add(sof)
        elif snapshot.frameSets is None:
            changed.add(sof)
        elif snapshot.frameSets.get(sof, frozenset()) == rebuiltFrameSets.get(sof, frozenset()):
            toRestore[sof] = snapshot.statuses[sof]
        else:
            pending.add(sof)
    return status_pass_result(
        toRestore=toRestore,
        pending=frozenset(pending),
        changed=frozenset(changed),
        dropped=frozenset(dropped),
    )


@contextmanager
def open_backup_read_only(backupPath: Path | str) -> Iterator[sqlite3.Connection]:
    """*open a preserved database read-only through a scratch copy, after checking that SQLite can trust it*

    The database and any `-wal`/`-shm` sidecar files are copied to a temporary directory outside the backups
    directory, and the copy is opened. SQLite may create or rewrite sidecar files when it reads a WAL-mode
    database, so reading the original could change `backups/`; reading the copy cannot. Never opens with
    `immutable=1`. The scratch directory is removed when the block ends, including on error.

    **Key Arguments:**

    - ``backupPath`` -- path of the preserved database file

    **Return:**

    - ``connection`` -- a read-only `sqlite3.Connection` to the scratch copy, closed when the block ends

    **Raises:**

    - `sqlite3.DatabaseError` when the file fails `PRAGMA quick_check`
    - `sqlite3.OperationalError` when the file cannot be opened
    """
    sourcePath = Path(backupPath).resolve()
    if not sourcePath.is_file():
        raise sqlite3.OperationalError(f"unable to open the preserved database `{sourcePath}`")
    with tempfile.TemporaryDirectory(prefix=SCRATCH_DIRECTORY_PREFIX) as scratchDirectory:
        scratchPath = Path(scratchDirectory) / sourcePath.name
        for suffix in DATABASE_FILE_SUFFIXES:
            sourceFile = Path(f"{sourcePath}{suffix}")
            if sourceFile.exists():
                shutil.copyfile(sourceFile, f"{scratchPath}{suffix}")
        connection = sqlite3.connect(scratchPath.as_uri() + "?mode=ro", uri=True)
        try:
            check = connection.execute("PRAGMA quick_check;").fetchall()
            if check != [("ok",)]:
                raise sqlite3.DatabaseError(f"the preserved database failed its quick check: {bounded_repr(check)}")
            yield connection
        finally:
            connection.close()


def read_session_snapshots(
    backupPath: Path | str, sessionIds: Iterable[str], log: logging.Logger
) -> dict[str, session_status_snapshot]:
    """*read the status snapshot of every named session from a preserved database*

    **Key Arguments:**

    - ``backupPath`` -- path of the preserved database file
    - ``sessionIds`` -- the session IDs to read (already validated by the caller)
    - ``log`` -- logger

    **Return:**

    - ``snapshots`` -- session ID to `session_status_snapshot`, leaving out sessions the backup has no status column for

    **Raises:**

    - `sqlite3.Error` when the backup cannot be opened or fails its `quick_check`
    """
    snapshots = {}
    with open_backup_read_only(backupPath) as connection:
        for sessionId in sessionIds:
            snapshot = read_session_snapshot(connection, sessionId, log)
            if snapshot is not None:
                snapshots[sessionId] = snapshot
    return snapshots


def read_session_snapshot(
    connection: sqlite3.Connection, sessionId: str, log: logging.Logger
) -> session_status_snapshot | None:
    """*read one session's statuses and frame sets from an open preserved database*

    **Key Arguments:**

    - ``connection`` -- open connection to the preserved database
    - ``sessionId`` -- the session ID to read; callers validate it, and the names composed from it are checked
      again by `validate_sql_identifier`
    - ``log`` -- logger

    **Return:**

    - ``snapshot`` -- a `session_status_snapshot`, or None when the backup has no `status_<id>` column

    **Raises:**

    - `UnsafeSqlIdentifierError` when a name composed from ``sessionId`` is not a safe SQL identifier
    """
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    sofMapTable = validate_sql_identifier(f"sof_map_{sessionId}", "sof map table name")
    columns = {row[1] for row in connection.execute("PRAGMA table_info(product_frames);")}
    if statusColumn not in columns:
        log.debug(f"session_status_restore: the backup has no `{statusColumn}` column, so `{sessionId}` is skipped")
        return None

    frameSets = None
    if table_exists(connection, sofMapTable):
        frameSets = read_frame_sets(connection, sofMapTable)
    else:
        log.debug(f"session_status_restore: the backup has no `{sofMapTable}` table, so its frames are unknown")
    return session_status_snapshot(
        sessionId=sessionId,
        statuses=read_statuses(connection, statusColumn, log),
        frameSets=frameSets,
    )


def read_statuses(connection: sqlite3.Connection, statusColumn: str, log: logging.Logger) -> dict[str, str]:
    """*read the restorable status of each SOF from one `status_<id>` column*

    A SOF that has both a `pass` and a `fail` row is `fail`. Values other than `pass` and `fail` are skipped.

    **Key Arguments:**

    - ``connection`` -- open connection to the preserved database
    - ``statusColumn`` -- a validated `status_<id>` column name
    - ``log`` -- logger

    **Return:**

    - ``statuses`` -- SOF name to ``pass`` or ``fail``
    """
    statusColumn = validate_sql_identifier(statusColumn, "status column")
    # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier AND COMPOSED FROM A VALIDATED SESSION ID
    sqlQuery = f"SELECT DISTINCT sof, {statusColumn} FROM product_frames WHERE {statusColumn} IS NOT NULL;"  # noqa: S608
    statuses = {}
    for sof, status in connection.execute(sqlQuery):
        if status not in RESTORABLE_STATUSES:
            log.debug(
                f"session_status_restore: skipping status {bounded_repr(status)} of {bounded_repr(sof)}; "
                "it is not pass or fail"
            )
            continue
        if sof in statuses and statuses[sof] != status:
            log.debug(
                f"session_status_restore: {bounded_repr(sof)} has both pass and fail rows, so it is restored as fail"
            )
            statuses[sof] = "fail"
            continue
        statuses[sof] = status
    return statuses


def read_frame_sets(connection: sqlite3.Connection, sofMapTable: str) -> dict[str, frozenset[str]]:
    """*read the set of `file` names each SOF holds in one SOF-map table*

    **Key Arguments:**

    - ``connection`` -- open connection to the database
    - ``sofMapTable`` -- a validated `sof_map_<id>` table name

    **Return:**

    - ``frameSets`` -- SOF name to a frozenset of file names
    """
    sofMapTable = validate_sql_identifier(sofMapTable, "sof map table name")
    # TABLE NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier AND COMPOSED FROM A VALIDATED SESSION ID
    sqlQuery = f"SELECT sof, file FROM {sofMapTable};"  # noqa: S608
    files = {}
    for sof, file in connection.execute(sqlQuery):
        files.setdefault(sof, set()).add(file)
    return {sof: frozenset(names) for sof, names in files.items()}


def table_exists(connection: sqlite3.Connection, tableName: str) -> bool:
    """*check whether the database has a table of this name*

    **Key Arguments:**

    - ``connection`` -- open connection to the database
    - ``tableName`` -- the table name to look for

    **Return:**

    - ``exists`` -- True when the table exists
    """
    sqlQuery = "SELECT 1 FROM sqlite_master WHERE type = 'table' AND name = ?;"
    return connection.execute(sqlQuery, (tableName,)).fetchone() is not None
