"""Integration contracts for keeping quality-control history across database rebuilds (DY-59)."""

from __future__ import annotations

import os
import shutil
import sqlite3
import struct
from pathlib import Path

import pytest

from soxspipe.commonutils.data_organiser import DatabasePreservationError
from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration

QC_COLUMNS = (
    "soxspipe_recipe",
    "qc_name",
    "qc_value",
    "qc_unit",
    "qc_comment",
    "obs_date_utc",
    "reduction_date_utc",
    "sof_name",
    "qc_order",
    "qc_flag",
    "qc_value_min",
    "qc_value_max",
)

# ONE WHOLE-FRAME ROW THAT THE GUARDRAIL WILL FAIL (MASTER RON RANGE IS [0,10]) AND ONE
# PER-ORDER ROW (NOT WRITTEN TO ANY FITS HEADER) THAT THE GUARDRAIL NEVER TOUCHES
SEEDED_QC_ROWS = (
    (
        "soxs-mbias",
        "MASTER RON",
        "50.0",
        "electrons",
        "combined read noise",
        "2024-01-02T03:04:05",
        "2024-01-03T03:04:05",
        "bias_sof.sof",
        "-1",
        "pass",
        None,
        None,
    ),
    (
        "soxs-order-centre",
        "ORDER RES",
        "0.25",
        "pixels",
        "per-order residual",
        "2024-01-02T03:04:05",
        "2024-01-03T03:04:05",
        "order_sof.sof",
        "12",
        "pass",
        "0",
        "1",
    ),
)

CORRUPT_BYTES = b"this is not an sqlite database" * 64


def _insert_complete_raw_frame(connection) -> None:
    """Seed the minimum complete raw-frame row used by session workflows."""
    connection.execute(
        """INSERT INTO raw_frames (
            instrume, file, "mjd-obs", "mjd-date", "date-obs",
            "night start mjd", "eso dpr catg", "eso dpr type", "eso dpr tech",
            exptime, object, filepath
        ) VALUES ('SOXS', 'bias-1.fits', 60311.1, '2024-01-02',
            '2024-01-02T03:04:05', 60310, 'CALIB', 'BIAS', 'IMAGE', 0, 'BIAS',
            './raw/bias-1.fits')"""
    )


def _insert_bias_product(connection) -> None:
    """Seed the product row whose SOF the failing QC row belongs to."""
    connection.execute(
        """INSERT INTO product_frames (file, recipe, sof, status_base)
        VALUES ('master_bias.fits', 'mbias', 'bias_sof.sof', 'pass')"""
    )


def _seed_quality_control(connection) -> None:
    """Insert the known `SEEDED_QC_ROWS` into `quality_control`."""
    placeholders = ", ".join("?" for _ in QC_COLUMNS)
    columns = ", ".join(QC_COLUMNS)
    connection.executemany(
        f"INSERT INTO quality_control ({columns}) VALUES ({placeholders})",  # noqa: S608
        SEEDED_QC_ROWS,
    )


def _read_qc(connection, columns=QC_COLUMNS):
    selected = ", ".join(columns)
    return sorted(
        connection.execute(f"SELECT {selected} FROM quality_control").fetchall(),  # noqa: S608
        key=repr,
    )


def _backup_files(organiser):
    backupDir = Path(organiser.rootDir) / "backups"
    if not backupDir.is_dir():
        return []
    return sorted(backupDir.glob("*.db"))


def _write_corrupt_but_openable_quality_control(path) -> None:
    """Write a SQLite file with a `quality_control` table that opens but fails `PRAGMA quick_check`."""
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA page_size=4096")
    connection.execute("CREATE TABLE quality_control (qc_name text, sof_name text, qc_order text)")
    connection.executemany(
        "INSERT INTO quality_control VALUES (?, ?, '-1')",
        [(f"QC {i} " + "x" * 50, "a.sof") for i in range(500)],
    )
    connection.commit()
    connection.execute("VACUUM")
    connection.close()
    # CORRUPT THE FREELIST PAGE COUNT IN THE FILE HEADER (BYTES 36-39)
    with open(path, "r+b") as databaseFile:
        databaseFile.seek(36)
        databaseFile.write(struct.pack(">I", 999))


@pytest.fixture
def open_sqlite():
    """Open SQLite connections that are closed at teardown."""
    connections = []

    def _open(path):
        connection = sqlite3.connect(path)
        connections.append(connection)
        return connection

    yield _open
    for connection in connections:
        connection.close()


@pytest.fixture
def science_workspace(tmp_path, log, monkeypatch):
    """A workspace in a named `science` session, with seeded QC history and stubbed FITS harvesting."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_complete_raw_frame(connection)
    organiser.session_create("science")
    _seed_quality_control(connection)
    organiser.syncCalls = []

    def fake_sync_raw_frames():
        organiser.syncCalls.append(True)
        _insert_complete_raw_frame(organiser.conn)
        _insert_bias_product(organiser.conn)

    monkeypatch.setattr(organiser, "_select_instrument", lambda: None)
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_sync_raw_frames", fake_sync_raw_frames)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_flag_files_to_ignore", lambda: None)
    monkeypatch.setattr(organiser, "build_sof_files", lambda: None)
    monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
    yield organiser
    if getattr(organiser, "conn", None) is not None:
        organiser.conn.close()


@pytest.fixture
def corrupt_workspace(science_workspace, monkeypatch):
    """The science workspace with its database replaced by bytes SQLite cannot open."""
    science_workspace.conn.close()
    science_workspace.conn = None
    Path(science_workspace.rootDbPath).write_bytes(CORRUPT_BYTES)
    monkeypatch.setattr("time.sleep", lambda seconds: None)
    return science_workspace


# EXPLICIT REFRESH


def test_refresh_restores_quality_control_rows_in_a_named_session(science_workspace, capsys) -> None:
    # ARRANGE
    organiser = science_workspace
    assert organiser.sessionId == "science"
    keyColumns = QC_COLUMNS[:9]
    expectedRows = sorted((row[:9] for row in SEEDED_QC_ROWS), key=repr)

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert _read_qc(organiser.conn, keyColumns) == expectedRows
    backups = _backup_files(organiser)
    assert len(backups) == 1
    output = capsys.readouterr().out
    restoreLines = [line for line in output.splitlines() if str(backups[0]) in line]
    assert len(restoreLines) == 1
    assert "Restored 2 of 2 quality-control rows" in restoreLines[0]


def test_refresh_keeps_a_complete_copy_of_the_old_database(science_workspace, open_sqlite) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute("INSERT INTO quality_control (qc_name, sof_name) VALUES ('ONLY IN OLD DB', 'x.sof')")
    oldTables = {row[0] for row in organiser.conn.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert organiser.databaseBackupPath == str(backups[0])
    backup = open_sqlite(backups[0])
    assert backup.execute("PRAGMA integrity_check").fetchall() == [("ok",)]
    backupTables = {row[0] for row in backup.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
    assert backupTables == oldTables
    assert backup.execute("SELECT count(*) FROM quality_control").fetchone() == (3,)
    assert backup.execute("SELECT count(*) FROM raw_frames").fetchone() == (1,)
    assert Path(organiser.rootDbPath).is_file()


def test_refresh_closes_the_old_connection(science_workspace) -> None:
    # ARRANGE
    oldConnection = science_workspace.conn

    # ACT
    science_workspace.prepare(refresh=True, report=False)

    # ASSERT
    with pytest.raises(sqlite3.ProgrammingError):
        oldConnection.execute("SELECT 1")


def test_refresh_fsyncs_the_backup_before_deleting_the_original(science_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = science_workspace
    events = []
    realFsync = os.fsync
    realRemove = organiser._remove_database_files

    def recording_fsync(fd):
        events.append(("fsync", fd))
        return realFsync(fd)

    def recording_remove():
        events.append(("remove", None))
        return realRemove()

    monkeypatch.setattr("os.fsync", recording_fsync)
    monkeypatch.setattr(organiser, "_remove_database_files", recording_remove)

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    removeIndex = events.index(("remove", None))
    fsyncsBeforeRemove = [event for event in events[:removeIndex] if event[0] == "fsync"]
    assert len(fsyncsBeforeRemove) >= 2


def test_backup_is_one_self_contained_file_after_refresh_and_restore(science_workspace) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute("PRAGMA journal_mode=WAL")
    backupDir = Path(organiser.dbBackupsDir)

    # ACT
    organiser.prepare(refresh=True, report=False)
    afterRefresh = sorted(path.name for path in backupDir.iterdir())
    organiser._restore_quality_control_history(organiser.databaseBackupPath)
    afterRestore = sorted(path.name for path in backupDir.iterdir())

    # ASSERT
    assert afterRefresh == [Path(organiser.databaseBackupPath).name]
    assert afterRestore == afterRefresh


def test_refresh_keeps_rows_still_in_the_write_ahead_log(science_workspace, open_sqlite) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute("PRAGMA journal_mode=WAL")
    organiser.conn.execute("PRAGMA wal_autocheckpoint=0")
    organiser.conn.execute("INSERT INTO quality_control (qc_name, sof_name, qc_order) VALUES ('IN WAL', 'w.sof', '-1')")
    walFile = Path(organiser.rootDbPath + "-wal")
    assert walFile.is_file() and walFile.stat().st_size > 0

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    backup = open_sqlite(organiser.databaseBackupPath)
    query = "SELECT count(*) FROM quality_control WHERE qc_name = 'IN WAL'"
    assert backup.execute(query).fetchone() == (1,)
    assert organiser.conn.execute(query).fetchone() == (1,)


def test_guardrail_recomputes_status_from_restored_rows(science_workspace) -> None:
    # ACT
    science_workspace.prepare(refresh=True, report=False)

    # ASSERT
    conn = science_workspace.conn
    assert conn.execute(
        "SELECT qc_flag, qc_value_min, qc_value_max FROM quality_control WHERE qc_name = 'MASTER RON'"
    ).fetchone() == ("fail", "0", "10")
    assert conn.execute("SELECT status_base FROM product_frames WHERE sof = 'bias_sof.sof'").fetchone() == ("fail",)


def test_refused_refresh_raises_and_leaves_the_database_untouched(science_workspace, monkeypatch, open_sqlite) -> None:
    # ARRANGE
    organiser = science_workspace
    rootDb = Path(organiser.rootDbPath)

    def failing_snapshot(backupPath):
        raise OSError("disk full")

    monkeypatch.setattr(organiser, "_snapshot_database", failing_snapshot)

    # ACT
    with pytest.raises(DatabasePreservationError, match=str(rootDb)):
        organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert organiser.syncCalls == []
    assert open_sqlite(rootDb).execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
    assert _backup_files(organiser) == []


def test_locked_database_is_refused_not_copied_as_corrupt(science_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = science_workspace

    def locked_snapshot(backupPath):
        raise sqlite3.OperationalError("database is locked")

    monkeypatch.setattr(organiser, "_snapshot_database", locked_snapshot)

    # ACT
    with pytest.raises(DatabasePreservationError):
        organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert organiser.syncCalls == []
    assert _backup_files(organiser) == []
    assert Path(organiser.rootDbPath).is_file()


def test_explicit_refresh_of_an_unreadable_database_keeps_a_raw_copy(tmp_path, log, monkeypatch) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    rootDb = Path(organiser.rootDbPath)
    rootDb.write_bytes(b"stale")
    Path(f"{rootDb}-wal").write_bytes(b"stale wal")
    monkeypatch.setattr(organiser, "_select_instrument", lambda: None)
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: False)

    # ACT
    with pytest.raises(SystemExit):
        organiser.prepare(refresh=True)

    # ASSERT
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert "corrupt" in backups[0].name
    assert backups[0].read_bytes() == b"stale"
    assert Path(f"{backups[0]}-wal").read_bytes() == b"stale wal"
    assert not rootDb.exists()
    assert not Path(f"{rootDb}-wal").exists()


def test_refresh_without_an_existing_database_skips_the_backup(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)

    # ACT
    backupPath = organiser._preserve_database()

    # ASSERT
    assert backupPath is None
    assert _backup_files(organiser) == []


# AUTOMATIC REBUILD


def test_automatic_rebuild_keeps_the_failed_database_and_survives_unreadable_source(corrupt_workspace, capsys) -> None:
    # ARRANGE
    organiser = corrupt_workspace

    # ACT
    organiser.conn, reset = organiser._get_or_create_db_connection()

    # ASSERT
    assert reset is True
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert "corrupt" in backups[0].name
    assert backups[0].read_bytes() == CORRUPT_BYTES
    assert organiser.conn.execute("PRAGMA integrity_check").fetchall() == [("ok",)]
    assert organiser.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (0,)
    output = capsys.readouterr().out
    warnings = [line for line in output.splitlines() if "WARNING" in line and str(backups[0]) in line]
    assert len(warnings) == 1
    assert "failed to open" in warnings[0]


def test_automatic_rebuild_restores_rows_when_the_failed_file_is_readable(science_workspace, capsys) -> None:
    # ACT
    science_workspace.prepare(refresh=True, report=False, _failedToOpen=True)

    # ASSERT
    backups = _backup_files(science_workspace)
    assert len(backups) == 1
    assert "corrupt" in backups[0].name
    assert science_workspace.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
    restoreLines = [line for line in capsys.readouterr().out.splitlines() if "Restored 2 of 2" in line]
    assert len(restoreLines) == 1
    assert "failed to open" in restoreLines[0]


def test_automatic_rebuild_raises_and_keeps_the_database_when_it_cannot_be_copied(
    corrupt_workspace, monkeypatch
) -> None:
    # ARRANGE
    organiser = corrupt_workspace
    rootDb = Path(organiser.rootDbPath)

    def failing_copy(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr("shutil.copy2", failing_copy)

    # ACT / ASSERT
    with pytest.raises(DatabasePreservationError, match=str(rootDb)):
        organiser._get_or_create_db_connection()
    assert issubclass(DatabasePreservationError, sqlite3.DatabaseError)
    assert organiser.syncCalls == []
    assert rootDb.read_bytes() == CORRUPT_BYTES
    assert _backup_files(organiser) == []


def test_partial_raw_copy_is_removed_and_the_database_kept(science_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = science_workspace
    rootDb = Path(organiser.rootDbPath)
    walFile = Path(organiser.rootDbPath + "-wal")
    walFile.write_bytes(b"wal frames")
    realCopy = shutil.copy2

    def copy_main_then_fail(source, destination):
        if str(source).endswith("-wal"):
            raise OSError("disk full")
        return realCopy(source, destination)

    monkeypatch.setattr("shutil.copy2", copy_main_then_fail)

    # ACT
    with pytest.raises(DatabasePreservationError):
        organiser._preserve_database(failedToOpen=True)

    # ASSERT
    assert list(Path(organiser.dbBackupsDir).iterdir()) == []
    assert rootDb.is_file()
    assert walFile.read_bytes() == b"wal frames"


def test_raw_copy_with_the_wrong_size_is_refused(science_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = science_workspace

    def truncating_copy(source, destination):
        Path(destination).write_bytes(Path(source).read_bytes()[:10])

    monkeypatch.setattr("shutil.copy2", truncating_copy)

    # ACT
    with pytest.raises(DatabasePreservationError, match="size"):
        organiser._preserve_database(failedToOpen=True)

    # ASSERT
    assert list(Path(organiser.dbBackupsDir).iterdir()) == []
    assert Path(organiser.rootDbPath).is_file()


# RESTORE


def test_restore_is_idempotent(science_workspace, capsys) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.prepare(refresh=True, report=False)
    capsys.readouterr()

    # ACT
    restored = organiser._restore_quality_control_history(organiser.databaseBackupPath)

    # ASSERT
    assert restored == 0
    assert organiser.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
    assert "Restored 0 of 2 quality-control rows" in capsys.readouterr().out


def test_restore_is_idempotent_for_rows_with_a_null_order(science_workspace) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute(
        "INSERT INTO quality_control (qc_name, sof_name, qc_order) VALUES ('NULL ORDER', 'n.sof', NULL)"
    )
    organiser.prepare(refresh=True, report=False)

    # ACT
    restored = organiser._restore_quality_control_history(organiser.databaseBackupPath)

    # ASSERT
    assert restored == 0
    assert organiser.conn.execute(
        "SELECT count(*) FROM quality_control WHERE qc_name = 'NULL ORDER' AND qc_order IS NULL"
    ).fetchone() == (1,)


def test_restore_round_trips_every_column_unchanged(science_workspace, open_sqlite) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute(
        "INSERT INTO quality_control (soxspipe_recipe, qc_name, qc_value, sof_name, qc_order) "
        "VALUES ('soxs-mdark', 'SENTINEL', '-99.99', 'dark_sof.sof', '-1')"
    )
    organiser.prepare(refresh=True, report=False)
    preservedRows = _read_qc(open_sqlite(organiser.databaseBackupPath))
    organiser.conn.execute("DELETE FROM quality_control")

    # ACT
    restored = organiser._restore_quality_control_history(organiser.databaseBackupPath)

    # ASSERT
    assert restored == 3
    assert _read_qc(organiser.conn) == preservedRows
    assert organiser.conn.execute("SELECT qc_value FROM quality_control WHERE qc_name = 'SENTINEL'").fetchone() == (
        "-99.99",
    )


def test_restore_skips_columns_the_current_schema_does_not_have(science_workspace, tmp_path, log, open_sqlite) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.prepare(refresh=True, report=False)
    organiser.conn.execute("DELETE FROM quality_control")
    olderBackup = tmp_path / "older-schema.db"
    older = open_sqlite(olderBackup)
    older.execute("CREATE TABLE quality_control (qc_name text, sof_name text, qc_order text, legacy_col text)")
    older.execute("INSERT INTO quality_control VALUES ('OLD QC', 'old.sof', '-1', 'gone')")
    older.commit()

    # ACT
    restored = organiser._restore_quality_control_history(str(olderBackup))

    # ASSERT
    assert restored == 1
    assert organiser.conn.execute("SELECT qc_name, sof_name FROM quality_control").fetchall() == [("OLD QC", "old.sof")]
    assert "legacy_col" in " ".join(message for _, message in log.messages)


def test_restore_warns_and_returns_zero_when_the_source_cannot_be_read(science_workspace, tmp_path, capsys) -> None:
    # ARRANGE
    missingPath = str(tmp_path / "no-such-backup.db")

    # ACT
    restored = science_workspace._restore_quality_control_history(missingPath)

    # ASSERT
    assert restored == 0
    output = capsys.readouterr().out
    assert any("WARNING" in line and missingPath in line for line in output.splitlines())
    assert not Path(missingPath).exists()


def test_restore_skips_a_source_that_fails_its_quick_check(science_workspace, tmp_path, capsys) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute("DELETE FROM quality_control")
    damagedPath = tmp_path / "damaged.db"
    _write_corrupt_but_openable_quality_control(damagedPath)

    # ACT
    restored = organiser._restore_quality_control_history(str(damagedPath))

    # ASSERT
    assert restored == 0
    assert organiser.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (0,)
    output = capsys.readouterr().out
    assert any("WARNING" in line and str(damagedPath) in line for line in output.splitlines())


# WORKSPACE HOUSEKEEPING


def test_move_misc_files_leaves_the_backups_directory_in_place(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    backupDir = Path(organiser.rootDir) / "backups"
    backupDir.mkdir()
    backupFile = backupDir / "soxspipe_refresh_20260101T000000000000Z.db"
    backupFile.write_bytes(b"backup")

    # ACT
    organiser._move_misc_files()

    # ASSERT
    assert backupFile.read_bytes() == b"backup"


def test_snapshot_that_fails_its_quick_check_falls_back_to_a_raw_copy(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    _write_corrupt_but_openable_quality_control(organiser.rootDbPath)
    originalBytes = Path(organiser.rootDbPath).read_bytes()

    # ACT
    backupPath = organiser._preserve_database()

    # ASSERT
    assert "_corrupt_" in Path(backupPath).name
    assert Path(backupPath).read_bytes() == originalBytes
    assert [path.name for path in Path(organiser.dbBackupsDir).iterdir()] == [Path(backupPath).name]


def _sqlite_error(errorClass, message, errorCode):
    error = errorClass(message)
    error.sqlite_errorcode = errorCode
    return error


def test_automatic_rebuild_refuses_a_locked_database(science_workspace) -> None:
    # ARRANGE
    organiser = science_workspace
    lockedError = _sqlite_error(sqlite3.OperationalError, "database is locked", sqlite3.SQLITE_BUSY)

    # ACT
    with pytest.raises(DatabasePreservationError, match="locked"):
        organiser._rebuild_database_that_failed_to_open(lockedError)

    # ASSERT
    assert organiser.syncCalls == []
    assert _backup_files(organiser) == []
    assert Path(organiser.rootDbPath).is_file()


def test_automatic_rebuild_snapshots_a_database_that_is_not_known_to_be_unreadable(
    science_workspace, open_sqlite
) -> None:
    # ARRANGE
    organiser = science_workspace
    integrityError = sqlite3.DatabaseError("database integrity check failed: [('row 1 missing',)]")

    # ACT
    organiser._rebuild_database_that_failed_to_open(integrityError)

    # ASSERT
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert "_refresh_" in backups[0].name
    rebuilt = open_sqlite(organiser.rootDbPath)
    assert rebuilt.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
