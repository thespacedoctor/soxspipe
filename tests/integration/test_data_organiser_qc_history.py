"""Integration contracts for keeping quality-control history across database rebuilds (DY-59)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pytest

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
    placeholders = ", ".join("?" for _ in QC_COLUMNS)
    columns = ", ".join(QC_COLUMNS)
    connection.executemany(
        f"INSERT INTO quality_control ({columns}) VALUES ({placeholders})",  # noqa: S608
        SEEDED_QC_ROWS,
    )


def _read_qc(connection, columns=QC_COLUMNS):
    selected = ", ".join(columns)
    return sorted(
        connection.execute(f"SELECT {selected} FROM quality_control").fetchall()  # noqa: S608
    )


def _backup_files(organiser):
    backupDir = Path(organiser.rootDir) / "backups"
    if not backupDir.is_dir():
        return []
    return sorted(backupDir.glob("*.db"))


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

    def fake_sync_raw_frames():
        _insert_complete_raw_frame(organiser.conn)
        _insert_bias_product(organiser.conn)

    monkeypatch.setattr(organiser, "_select_instrument", lambda: None)
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_sync_raw_frames", fake_sync_raw_frames)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_flag_files_to_ignore", lambda: None)
    monkeypatch.setattr(organiser, "build_sof_files", lambda: None)
    monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
    return organiser


def test_refresh_restores_quality_control_rows_in_a_named_session(science_workspace, capsys) -> None:
    # ARRANGE
    organiser = science_workspace
    assert organiser.sessionId == "science"
    keyColumns = QC_COLUMNS[:9]
    expectedRows = sorted(row[:9] for row in SEEDED_QC_ROWS)

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


def test_refresh_keeps_a_complete_copy_of_the_old_database(science_workspace) -> None:
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
    backup = sqlite3.connect(backups[0])
    try:
        assert backup.execute("PRAGMA integrity_check").fetchall() == [("ok",)]
        backupTables = {row[0] for row in backup.execute("SELECT name FROM sqlite_master WHERE type = 'table'")}
        assert backupTables == oldTables
        assert backup.execute("SELECT count(*) FROM quality_control").fetchone() == (3,)
        assert backup.execute("SELECT count(*) FROM raw_frames").fetchone() == (1,)
    finally:
        backup.close()
    assert Path(organiser.rootDbPath).is_file()


def test_guardrail_recomputes_status_from_restored_rows(science_workspace) -> None:
    # ACT
    science_workspace.prepare(refresh=True, report=False)

    # ASSERT
    conn = science_workspace.conn
    assert conn.execute(
        "SELECT qc_flag, qc_value_min, qc_value_max FROM quality_control WHERE qc_name = 'MASTER RON'"
    ).fetchone() == ("fail", "0", "10")
    assert conn.execute("SELECT status_base FROM product_frames WHERE sof = 'bias_sof.sof'").fetchone() == ("fail",)


def test_restore_is_idempotent(science_workspace, capsys) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.prepare(refresh=True, report=False)
    backupPath = organiser.databaseBackupPath
    capsys.readouterr()

    # ACT
    restored = organiser._restore_quality_control_history(backupPath)

    # ASSERT
    assert restored == 0
    assert organiser.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
    assert "Restored 0 of 2 quality-control rows" in capsys.readouterr().out


def test_restore_round_trips_every_column_unchanged(science_workspace) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.execute(
        "INSERT INTO quality_control (soxspipe_recipe, qc_name, qc_value, sof_name, qc_order) "
        "VALUES ('soxs-mdark', 'SENTINEL', '-99.99', 'dark_sof.sof', '-1')"
    )
    organiser.prepare(refresh=True, report=False)
    backup = sqlite3.connect(organiser.databaseBackupPath)
    try:
        preservedRows = _read_qc(backup)
    finally:
        backup.close()
    organiser.conn.execute("DELETE FROM quality_control")

    # ACT
    restored = organiser._restore_quality_control_history(organiser.databaseBackupPath)

    # ASSERT
    assert restored == 3
    assert _read_qc(organiser.conn) == preservedRows
    assert organiser.conn.execute("SELECT qc_value FROM quality_control WHERE qc_name = 'SENTINEL'").fetchone() == (
        "-99.99",
    )


def test_restore_skips_columns_the_current_schema_does_not_have(science_workspace, tmp_path, log) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.prepare(refresh=True, report=False)
    organiser.conn.execute("DELETE FROM quality_control")
    olderBackup = tmp_path / "older-schema.db"
    older = sqlite3.connect(olderBackup)
    try:
        older.execute("CREATE TABLE quality_control (qc_name text, sof_name text, qc_order text, legacy_col text)")
        older.execute("INSERT INTO quality_control VALUES ('OLD QC', 'old.sof', '-1', 'gone')")
        older.commit()
    finally:
        older.close()

    # ACT
    restored = organiser._restore_quality_control_history(str(olderBackup))

    # ASSERT
    assert restored == 1
    assert organiser.conn.execute("SELECT qc_name, sof_name FROM quality_control").fetchall() == [("OLD QC", "old.sof")]
    assert "legacy_col" in " ".join(message for _, message in log.messages)


def test_automatic_rebuild_moves_the_failed_database_and_survives_unreadable_source(
    science_workspace, monkeypatch, capsys
) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.close()
    organiser.conn = None
    corruptBytes = b"this is not an sqlite database" * 64
    Path(organiser.rootDbPath).write_bytes(corruptBytes)
    monkeypatch.setattr("time.sleep", lambda seconds: None)

    # ACT
    connection, reset = organiser._get_or_create_db_connection()

    # ASSERT
    assert reset is True
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert "corrupt" in backups[0].name
    assert backups[0].read_bytes() == corruptBytes
    assert connection.execute("PRAGMA integrity_check").fetchall() == [("ok",)]
    assert connection.execute("SELECT count(*) FROM quality_control").fetchone() == (0,)
    output = capsys.readouterr().out
    warnings = [line for line in output.splitlines() if "WARNING" in line and str(backups[0]) in line]
    assert len(warnings) == 1
    connection.close()


def test_automatic_rebuild_restores_rows_when_the_moved_file_is_readable(science_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = science_workspace
    monkeypatch.setattr(organiser, "_databaseFailedToOpen", True)

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    backups = _backup_files(organiser)
    assert len(backups) == 1
    assert "corrupt" in backups[0].name
    assert organiser.conn.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)


def test_refresh_does_not_delete_the_database_when_it_cannot_be_preserved(
    science_workspace, monkeypatch, capsys
) -> None:
    # ARRANGE
    organiser = science_workspace
    rootDb = Path(organiser.rootDbPath)

    def failing_snapshot(backupPath):
        raise OSError("disk full")

    monkeypatch.setattr(organiser, "_snapshot_database", failing_snapshot)

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert rootDb.is_file()
    check = sqlite3.connect(rootDb)
    try:
        assert check.execute("SELECT count(*) FROM quality_control").fetchone() == (2,)
    finally:
        check.close()
    output = capsys.readouterr().out
    assert any("WARNING" in line and str(rootDb) in line for line in output.splitlines())
    assert _backup_files(organiser) == []


def test_automatic_rebuild_raises_and_keeps_the_database_when_it_cannot_be_copied(
    science_workspace, monkeypatch
) -> None:
    # ARRANGE
    organiser = science_workspace
    organiser.conn.close()
    organiser.conn = None
    corruptBytes = b"this is not an sqlite database" * 64
    rootDb = Path(organiser.rootDbPath)
    rootDb.write_bytes(corruptBytes)
    monkeypatch.setattr("time.sleep", lambda seconds: None)

    def failing_copy(*args, **kwargs):
        raise OSError("disk full")

    monkeypatch.setattr("shutil.copy2", failing_copy)

    # ACT / ASSERT
    with pytest.raises(sqlite3.DatabaseError, match="left in place"):
        organiser._get_or_create_db_connection()
    assert rootDb.read_bytes() == corruptBytes
    assert _backup_files(organiser) == []


def test_partial_raw_copy_is_removed_and_the_database_kept(science_workspace, monkeypatch) -> None:
    # ARRANGE
    import shutil

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
    backupPath = organiser._preserve_database(failedToOpen=True)

    # ASSERT
    assert backupPath is None
    assert list(Path(organiser.dbBackupsDir).iterdir()) == []
    assert rootDb.is_file()
    assert walFile.read_bytes() == b"wal frames"


def test_refresh_without_an_existing_database_skips_the_backup(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)

    # ACT
    backupPath = organiser._preserve_database()

    # ASSERT
    assert backupPath is None
    assert _backup_files(organiser) == []


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
