"""Integration contracts for restoring per-session pass/fail statuses across a database rebuild (DY-218)."""

from __future__ import annotations

import hashlib
import sqlite3
from contextlib import closing
from pathlib import Path

import pandas as pd
import pytest

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
from tests.factories import raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

INSTRUMENT = "SOXS"

MBIAS_ONE = "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"
MBIAS_TWO = "20240103T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"
DISP_SOLUTION = "20240102T040405_VIS_1X1_FAST_DISP_SOLUTION_10_0S_SOXS.sof"


def _raw_frames_columns(connection):
    return {row[1] for row in connection.execute("PRAGMA table_info(raw_frames);")}


def _night_two_biases(nightOne):
    """Two biases one day after `nightOne`, sharing every grouping key."""
    frames = nightOne.copy()
    frames["file"] = ["bias-3.fits", "bias-4.fits"]
    frames["filepath"] = [f"./raw/2024-01-02/{name}" for name in frames["file"]]
    frames["date-obs"] = ["2024-01-03T03:04:05.000", "2024-01-03T03:04:06.000"]
    frames["mjd-obs"] = frames["mjd-obs"] + 1
    frames["eso tpl start"] = "2024-01-03T03:04:00"
    frames["night start date"] = "2024-01-02"
    frames["night start mjd"] = frames["night start mjd"] + 1
    frames["set_first_file"] = "bias-3.fits"
    return frames


def _pinhole_arc(nightOne):
    """One VIS pinhole arc taken between the two bias nights."""
    frame = nightOne.iloc[[0]].copy()
    frame["file"] = "arc-1.fits"
    frame["filepath"] = "./raw/2024-01-01/arc-1.fits"
    frame["eso dpr type"] = "LAMP,WAVE"
    frame["eso dpr tech"] = "ECHELLE,PINHOLE"
    frame["eso tpl nexp"] = 1
    frame["eso tpl name"] = "SOXS_arc"
    frame["exptime"] = 10.0
    frame["date-obs"] = "2024-01-02T04:04:05.000"
    frame["mjd-obs"] = frame["mjd-obs"] + 0.04
    frame["eso tpl start"] = "2024-01-02T04:04:00"
    frame["set_first_file"] = "arc-1.fits"
    return frame


def _insert_frames(connection, frames):
    """Insert dataframe rows into `raw_frames`, keeping only the columns the table has."""
    known = frames.loc[:, frames.columns.isin(_raw_frames_columns(connection))]
    known.to_sql("raw_frames", con=connection, index=False, if_exists="append")


def _replace_second_night_one_bias(nightOne):
    """Swap `bias-2` for a later `bias-5` that groups with `bias-1`, as when a re-download changes a set."""
    frames = nightOne.copy()
    frames.loc[frames["file"] == "bias-2.fits", ["file", "filepath", "date-obs"]] = [
        "bias-5.fits",
        "./raw/2024-01-01/bias-5.fits",
        "2024-01-02T03:04:09.000",
    ]
    return frames


def _insert_all_frames(connection, replaceSecondBias=False):
    nightOne = raw_frame_table()
    frames = [nightOne, _night_two_biases(nightOne), _pinhole_arc(nightOne)]
    if replaceSecondBias:
        frames[0] = _replace_second_night_one_bias(nightOne)
    _insert_frames(connection, pd.concat(frames))


def _prepared_workspace(tmp_path, log, monkeypatch, sessionIds):
    """Prepare a workspace with real SOF building, creating ``sessionIds`` in order (the last is current)."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_all_frames(connection)
    for sessionId in sessionIds:
        organiser.session_create(sessionId)

    organiser.replaceSecondBias = False

    def fake_sync_raw_frames():
        if not organiser.conn.execute("SELECT count(*) FROM raw_frames").fetchone()[0]:
            _insert_all_frames(organiser.conn, organiser.replaceSecondBias)
        organiser._select_instrument(inst=INSTRUMENT)

    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_sync_raw_frames", fake_sync_raw_frames)
    monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
    organiser.prepare(report=False)
    return organiser


@pytest.fixture
def two_session_workspace(tmp_path, log, monkeypatch):
    """A prepared workspace with an `archive` and a current `science` session and real SOF building."""
    organiser = _prepared_workspace(tmp_path, log, monkeypatch, ["archive", "science"])
    yield organiser
    if getattr(organiser, "conn", None) is not None:
        organiser.conn.close()


@pytest.fixture
def base_session_workspace(tmp_path, log, monkeypatch):
    """A prepared workspace whose only and current session is `base`."""
    organiser = _prepared_workspace(tmp_path, log, monkeypatch, [])
    yield organiser
    if getattr(organiser, "conn", None) is not None:
        organiser.conn.close()


def _seed_statuses(organiser, sessionId, statuses):
    """Write ``statuses`` (SOF name to status) into one session's `status_<id>` column."""
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    for sof, status in statuses.items():
        # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE (DY-254)
        organiser.conn.execute(
            f"UPDATE product_frames SET {statusColumn} = ? WHERE sof = ?",  # noqa: S608
            (status, sof),
        )


def _seed_archive_session(organiser, statuses):
    """Give the non-current `archive` session the current SOF map plus ``statuses``."""
    organiser.conn.execute("INSERT INTO sof_map_archive SELECT * FROM sof_map_science")
    _seed_statuses(organiser, "archive", statuses)


def _read_statuses(organiser, sessionId):
    """Return SOF name to status for one session, including NULL."""
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE (DY-254)
    rows = organiser.conn.execute(
        f"SELECT DISTINCT sof, {statusColumn} FROM product_frames"  # noqa: S608
    ).fetchall()
    return dict(rows)


def _summary_lines(capsys):
    return [line for line in capsys.readouterr().out.splitlines() if line.startswith("session '")]


def _summary_lines_in(output):
    return [line for line in output.splitlines() if line.startswith("session '")]


def _debug_messages(log):
    return [message for level, message in log.messages if level == "debug"]


# RESTORING STATUSES


def test_refresh_restores_pass_fail_and_null_statuses_in_every_session(two_session_workspace) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass", MBIAS_TWO: "fail"})
    _seed_archive_session(organiser, {MBIAS_ONE: "fail", DISP_SOLUTION: "pass"})

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: "pass", MBIAS_TWO: "fail", DISP_SOLUTION: None}
    assert _read_statuses(organiser, "archive") == {MBIAS_ONE: "fail", MBIAS_TWO: None, DISP_SOLUTION: "pass"}


def test_restored_fail_with_in_range_quality_control_rows_ends_pass(two_session_workspace) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "fail", MBIAS_TWO: "fail"})
    organiser.conn.execute(
        "INSERT INTO quality_control (soxspipe_recipe, qc_name, qc_value, sof_name, qc_order, qc_flag) "
        "VALUES ('soxs-mbias', 'MASTER RON', '5.0', ?, '-1', 'pass')",
        (MBIAS_ONE,),
    )

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: "pass", MBIAS_TWO: "fail", DISP_SOLUTION: None}


def _calibration_products(organiser, sof):
    """Return the master-bias products held in one SOF of the current session's map."""
    rows = organiser.conn.execute("SELECT file FROM sof_map_science WHERE sof = ? AND file LIKE '%MBIAS%'", (sof,))
    return {row[0] for row in rows}


def test_refresh_restores_a_downstream_status_whose_calibration_was_rematched_after_a_failure(
    two_session_workspace,
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    firstBias = MBIAS_ONE.replace(".sof", ".fits")
    secondBias = MBIAS_TWO.replace(".sof", ".fits")
    assert _calibration_products(organiser, DISP_SOLUTION) == {firstBias}
    _seed_statuses(organiser, "science", {MBIAS_ONE: "fail"})
    organiser.session_refresh(failure=None)
    assert _calibration_products(organiser, DISP_SOLUTION) == {secondBias}
    _seed_statuses(organiser, "science", {DISP_SOLUTION: "pass"})

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: "fail", MBIAS_TWO: None, DISP_SOLUTION: "pass"}
    assert _calibration_products(organiser, DISP_SOLUTION) == {secondBias}


def test_pass_limit_counts_the_sofs_still_differing_as_changed(two_session_workspace, log, capsys) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "fail"})
    organiser.session_refresh(failure=None)
    _seed_statuses(organiser, "science", {DISP_SOLUTION: "pass"})
    organiser._STATUS_RESTORE_MAX_PASSES = 1
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    (scienceLine,) = [line for line in _summary_lines(capsys) if line.startswith("session 'science'")]
    assert scienceLine.startswith("session 'science': 1 statuses restored, 1 changed (frames differ, requeued)")
    assert _read_statuses(organiser, "science")[DISP_SOLUTION] is None
    assert any("pass limit" in message for message in _debug_messages(log))


# QC FAILURES MARK THE CURRENT SESSION (DY-272)


def _insert_master_ron(organiser, sof, value):
    """Give ``sof`` a whole-frame MASTER RON QC row; the default acceptable range is [0, 10]."""
    organiser.conn.execute(
        "INSERT INTO quality_control (soxspipe_recipe, qc_name, qc_value, sof_name, qc_order, qc_flag) "
        "VALUES ('soxs-mbias', 'MASTER RON', ?, ?, '-1', 'pass')",
        (value, sof),
    )


def test_out_of_range_qc_fails_the_sof_in_the_named_session_and_leaves_status_base_alone(
    two_session_workspace,
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass"})
    _seed_statuses(organiser, "base", {MBIAS_ONE: "pass"})
    _insert_master_ron(organiser, MBIAS_ONE, "50.0")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert _read_statuses(organiser, "science")[MBIAS_ONE] == "fail"
    assert _read_statuses(organiser, "base")[MBIAS_ONE] == "pass"


def test_out_of_range_qc_fails_the_sof_in_the_base_session(base_session_workspace) -> None:
    # ARRANGE
    organiser = base_session_workspace
    _seed_statuses(organiser, "base", {MBIAS_ONE: "pass"})
    _insert_master_ron(organiser, MBIAS_ONE, "50.0")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert organiser.sessionId == "base"
    assert _read_statuses(organiser, "base")[MBIAS_ONE] == "fail"


@pytest.mark.parametrize(
    ("workspaceFixture", "sessionId"),
    [("base_session_workspace", "base"), ("two_session_workspace", "science")],
)
def test_in_range_qc_ends_pass_in_the_current_session(request, workspaceFixture, sessionId) -> None:
    # ARRANGE
    organiser = request.getfixturevalue(workspaceFixture)
    _seed_statuses(organiser, sessionId, {MBIAS_ONE: "fail"})
    _insert_master_ron(organiser, MBIAS_ONE, "5.0")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert _read_statuses(organiser, sessionId)[MBIAS_ONE] == "pass"


# REPORTING


def test_refresh_prints_one_summary_line_per_session(two_session_workspace, log, capsys) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass", MBIAS_TWO: "fail"})
    _seed_archive_session(organiser, {MBIAS_ONE: "fail", DISP_SOLUTION: "pass"})
    (Path(organiser.sessionsDir) / "late").mkdir()
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    tail = "0 changed (frames differ, requeued), 0 dropped (product no longer exists)"
    assert _summary_lines(capsys) == [
        f"session 'archive': 2 statuses restored, {tail}",
        f"session 'science': 2 statuses restored, {tail}",
    ]
    assert any("status_late" in message for message in _debug_messages(log))


def test_sof_only_in_the_backup_is_dropped_and_counted_with_the_backup_path(two_session_workspace, capsys) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass"})
    organiser.conn.execute(
        "INSERT INTO product_frames (file, recipe, sof, status_science) "
        "VALUES ('gone.fits', 'mbias', 'GONE.sof', 'fail')"
    )
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    (scienceLine,) = [line for line in _summary_lines(capsys) if line.startswith("session 'science'")]
    assert scienceLine == (
        "session 'science': 1 statuses restored, 0 changed (frames differ, requeued), "
        f"1 dropped (product no longer exists); the previous statuses are kept in {organiser.databaseBackupPath}"
    )
    assert "GONE.sof" not in _read_statuses(organiser, "science")
    with closing(sqlite3.connect(organiser.databaseBackupPath)) as backup:
        goneStatus = backup.execute("SELECT status_science FROM product_frames WHERE sof = 'GONE.sof'").fetchone()
    assert goneStatus == ("fail",)


def test_sof_with_same_name_but_different_frames_is_not_restored_and_counted_changed(
    two_session_workspace, capsys
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "fail", MBIAS_TWO: "pass"})
    organiser.replaceSecondBias = True
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    (scienceLine,) = [line for line in _summary_lines(capsys) if line.startswith("session 'science'")]
    assert scienceLine == (
        "session 'science': 1 statuses restored, 1 changed (frames differ, requeued), "
        f"0 dropped (product no longer exists); the previous statuses are kept in {organiser.databaseBackupPath}"
    )
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: None, MBIAS_TWO: "pass", DISP_SOLUTION: None}


# BACKUPS THAT CANNOT BE FULLY USED


def _alter_backup_after_it_is_written(organiser, monkeypatch, alter):
    """Make `_preserve_database` hand ``alter`` the backup path right after the real backup is written."""
    realPreserve = organiser._preserve_database

    def altering_preserve(failedToOpen=False):
        backupPath = realPreserve(failedToOpen=failedToOpen)
        alter(Path(backupPath))
        return backupPath

    monkeypatch.setattr(organiser, "_preserve_database", altering_preserve)


def _drop_science_sof_map(backupPath):
    with closing(sqlite3.connect(backupPath)) as backup:
        backup.execute("DROP VIEW sof_map")
        backup.execute("DROP TABLE sof_map_science")
        backup.commit()


def test_refresh_with_a_corrupt_backup_completes_with_a_warning_and_empty_statuses(
    two_session_workspace, monkeypatch, capsys, log
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass"})
    _alter_backup_after_it_is_written(organiser, monkeypatch, lambda path: path.write_bytes(b"not a database" * 64))
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    output = capsys.readouterr().out
    warnings = [line for line in output.splitlines() if line.startswith("WARNING: could not restore session statuses")]
    assert len(warnings) == 1
    assert organiser.databaseBackupPath in warnings[0]
    assert warnings[0].endswith("Status columns are left empty; the preserved file has been kept.")
    assert _summary_lines_in(output) == []
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: None, MBIAS_TWO: None, DISP_SOLUTION: None}
    assert any("session statuses" in message for level, message in log.messages if level == "warning")


def test_restoring_from_a_missing_backup_warns_and_returns_no_counts(two_session_workspace, tmp_path, capsys) -> None:
    # ARRANGE
    capsys.readouterr()

    # ACT
    counts = two_session_workspace._restore_session_statuses(tmp_path / "missing.db")

    # ASSERT
    assert counts == {}
    assert "WARNING: could not restore session statuses" in capsys.readouterr().out


def test_backup_without_a_session_sof_map_counts_that_sessions_sofs_as_changed(
    two_session_workspace, monkeypatch, capsys
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass", MBIAS_TWO: "fail"})
    _alter_backup_after_it_is_written(organiser, monkeypatch, _drop_science_sof_map)
    capsys.readouterr()

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    (scienceLine,) = [line for line in _summary_lines(capsys) if line.startswith("session 'science'")]
    assert scienceLine.startswith("session 'science': 0 statuses restored, 2 changed (frames differ, requeued)")
    assert _read_statuses(organiser, "science") == {MBIAS_ONE: None, MBIAS_TWO: None, DISP_SOLUTION: None}


def test_status_restore_leaves_the_backups_directory_unchanged(two_session_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = two_session_workspace
    _seed_statuses(organiser, "science", {MBIAS_ONE: "pass", MBIAS_TWO: "fail"})
    seen = {}
    realRestore = organiser._restore_session_statuses

    def snapshot_backups():
        backupDir = Path(organiser.dbBackupsDir)
        return sorted((path.name, hashlib.sha256(path.read_bytes()).hexdigest()) for path in backupDir.iterdir())

    def recording_restore(backupPath):
        seen["before"] = snapshot_backups()
        counts = realRestore(backupPath)
        seen["after"] = snapshot_backups()
        return counts

    monkeypatch.setattr(organiser, "_restore_session_statuses", recording_restore)

    # ACT
    organiser.prepare(refresh=True, report=False)

    # ASSERT
    assert len(seen["before"]) == 1
    assert seen["after"] == seen["before"]


# CURSOR HYGIENE


class _SpyConnection:
    """Stand-in connection whose single cursor records `close` and can fail on `execute`."""

    def __init__(self, failOnExecute):
        self.failOnExecute = failOnExecute
        self.cursorClosed = False

    def cursor(self):
        return self

    def execute(self, sqlQuery, sqlParams=()):
        if self.failOnExecute:
            raise sqlite3.OperationalError("simulated failure")

    def commit(self):
        pass

    def close(self):
        self.cursorClosed = True


@pytest.mark.parametrize("failOnExecute", [False, True])
def test_qc_acceptable_range_pass_closes_its_cursor_even_when_a_statement_fails(tmp_path, log, failOnExecute) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "science"
    organiser.settings = {}
    organiser.conn = _SpyConnection(failOnExecute)

    # ACT
    if failOnExecute:
        with pytest.raises(sqlite3.OperationalError):
            organiser._apply_qc_acceptable_ranges()
    else:
        organiser._apply_qc_acceptable_ranges()

    # ASSERT
    assert organiser.conn.cursorClosed
