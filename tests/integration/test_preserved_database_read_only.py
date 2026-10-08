"""Integration contracts for reading preserved databases without changing `backups/` (DY-273)."""

from __future__ import annotations

import hashlib
import sqlite3
import tempfile
from pathlib import Path

import pytest

from soxspipe.commonutils.session_status_restore import read_session_snapshots
from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration

SESSION_ID = "science"

QC_ROW = ("soxs-mbias", "MASTER RON", "2.5", "bias_sof.sof", "-1")

STATUS_ROW = ("master_bias.fits", "bias_sof.sof", "pass")


def _write_wal_database(path: Path) -> sqlite3.Connection:
    """Write a WAL-mode database with QC and status rows, and return the writer, still holding the WAL open.

    Auto-checkpointing is off, so the rows stay in the `-wal` file until the writer closes.
    """
    writer = sqlite3.connect(path)
    writer.execute("PRAGMA journal_mode=WAL;")
    writer.execute("PRAGMA wal_autocheckpoint=0;")
    writer.execute(
        "CREATE TABLE quality_control (soxspipe_recipe TEXT, qc_name TEXT, qc_value TEXT, sof_name TEXT, qc_order TEXT)"
    )
    writer.execute(f"CREATE TABLE product_frames (file TEXT, sof TEXT, status_{SESSION_ID} TEXT)")
    writer.execute(f"CREATE TABLE sof_map_{SESSION_ID} (sof TEXT, file TEXT)")
    writer.execute("INSERT INTO quality_control VALUES (?, ?, ?, ?, ?)", QC_ROW)
    writer.execute("INSERT INTO product_frames VALUES (?, ?, ?)", STATUS_ROW)
    writer.execute(f"INSERT INTO sof_map_{SESSION_ID} VALUES ('bias_sof.sof', 'bias-1.fits')")  # noqa: S608
    writer.commit()
    return writer


def _raw_copy(organiser, sourcePath: Path) -> Path:
    """Copy ``sourcePath`` and its sidecars into the backups directory the way a raw `_corrupt_` backup is made."""
    backupsDir = Path(organiser.dbBackupsDir)
    backupsDir.mkdir(parents=True, exist_ok=True)
    backupPath = backupsDir / "soxspipe_corrupt_20260101T000000.db"
    organiser.rootDbPath = str(sourcePath)
    organiser._copy_database_files(str(backupPath))
    return backupPath


def _fingerprint(directory: Path) -> list[tuple[str, str]]:
    return sorted((path.name, hashlib.sha256(path.read_bytes()).hexdigest()) for path in directory.iterdir())


def _read_qc(organiser, backupPath, log):
    del log
    return [tuple(row) for row in organiser._read_preserved_quality_control(backupPath).itertuples(index=False)]


def _read_statuses(organiser, backupPath, log):
    del organiser
    snapshots = read_session_snapshots(backupPath, [SESSION_ID], log)
    return dict(snapshots[SESSION_ID].statuses)


READERS = {
    "quality_control": (_read_qc, [QC_ROW]),
    "session_statuses": (_read_statuses, {"bias_sof.sof": "pass"}),
}


@pytest.fixture
def organiser(tmp_path, log, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    return workspace_organiser(tmp_path / "workspace", log=log)


@pytest.fixture
def scratch_root(tmp_path, monkeypatch) -> Path:
    """Send every temporary directory to a known folder, so a leaked scratch directory is visible."""
    root = tmp_path / "scratch-root"
    root.mkdir()
    monkeypatch.setattr(tempfile, "tempdir", str(root))
    return root


@pytest.mark.parametrize("readerName", sorted(READERS))
def test_reading_a_raw_copy_without_sidecars_leaves_the_backups_directory_unchanged(
    organiser, tmp_path, log, readerName
) -> None:
    # ARRANGE
    writer = _write_wal_database(tmp_path / "source.db")
    writer.close()
    backupPath = _raw_copy(organiser, tmp_path / "source.db")
    assert [path.name for path in backupPath.parent.iterdir()] == [backupPath.name]
    before = _fingerprint(backupPath.parent)
    read, expected = READERS[readerName]

    # ACT
    result = read(organiser, backupPath, log)

    # ASSERT
    assert _fingerprint(backupPath.parent) == before
    assert result == expected


@pytest.mark.parametrize("readerName", sorted(READERS))
def test_reading_a_raw_copy_with_sidecars_leaves_the_backups_directory_unchanged(
    organiser, tmp_path, log, readerName
) -> None:
    # ARRANGE
    writer = _write_wal_database(tmp_path / "source.db")
    try:
        backupPath = _raw_copy(organiser, tmp_path / "source.db")
    finally:
        writer.close()
    assert sorted(path.name for path in backupPath.parent.iterdir()) == [
        backupPath.name,
        f"{backupPath.name}-shm",
        f"{backupPath.name}-wal",
    ]
    before = _fingerprint(backupPath.parent)
    read, expected = READERS[readerName]

    # ACT
    result = read(organiser, backupPath, log)

    # ASSERT
    assert _fingerprint(backupPath.parent) == before
    assert result == expected


@pytest.mark.parametrize("readerName", sorted(READERS))
def test_reading_a_preserved_database_leaves_no_scratch_directory_behind(
    organiser, tmp_path, scratch_root, log, readerName
) -> None:
    # ARRANGE
    writer = _write_wal_database(tmp_path / "source.db")
    try:
        backupPath = _raw_copy(organiser, tmp_path / "source.db")
    finally:
        writer.close()
    read, _ = READERS[readerName]

    # ACT
    read(organiser, backupPath, log)

    # ASSERT
    assert list(scratch_root.iterdir()) == []
