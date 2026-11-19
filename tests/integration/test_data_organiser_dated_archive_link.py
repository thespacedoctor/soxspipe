"""Integration contracts for a workspace-root link to a frame in a dated folder of an external raw archive (DY-2260)."""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from tests.factories import harvestable_raw_fits, workspace_layout, workspace_organiser

pytestmark = pytest.mark.integration

NIGHT = "2024-01-02"
# A ROW WITH A VALUE FOR EVERY NOT NULL COLUMN OF `raw_frames`, SHARING THE FILEPATH OF THE ARCHIVED FRAME
OTHER_ROW = {
    "uuid": 999999,
    "file": "x.fits",
    "mjd-obs": 60311.75,
    "mjd-date": NIGHT,
    "date-obs": f"{NIGHT}T06:00:00",
    "night start mjd": 60311,
    "eso dpr catg": "CALIB",
    "eso dpr type": "BIAS",
    "eso dpr tech": "OTHER,TECH",
    "exptime": 0.0,
    "object": "BIAS",
    "filepath": f"./raw/{NIGHT}/x.fits",
    "ignore": 0,
    "simulation": 0,
    "processed": 0,
}


def _archive_workspace(tmp_path: Path, log: Any) -> tuple[data_organiser, Path, Path]:
    """Return an organiser whose raw directory is an external archive, the workspace root and the archived frame.

    The workspace `raw/` folder is not the archive, so the frame is not where its `./raw/<date>/` path would file it.
    """
    layout = workspace_layout(tmp_path)
    archive = (tmp_path / "archive" / "raw").resolve()
    (archive / NIGHT).mkdir(parents=True)
    frame = harvestable_raw_fits(archive / NIGHT / "x.fits")
    (layout.root / "x.fits").symlink_to(frame)
    organiser = data_organiser(log=log, rootDir=str(layout.root), dbConnect=False)
    organiser.rawDir = str(archive)
    return organiser, layout.root, frame


def _raw_frame_paths(organiser: data_organiser) -> list[str]:
    return [row[0] for row in organiser.conn.execute("SELECT filepath FROM raw_frames").fetchall()]


def test_sync_raw_frames_deletes_the_row_of_a_dated_archive_frame_it_leaves_in_the_archive(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    organiser, _, frame = _archive_workspace(tmp_path, log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser._select_instrument()

    # ACT
    try:
        organiser._sync_raw_frames(skipSqlSync=True)
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []
    assert frame.is_file()


def test_sync_raw_frames_keeps_the_row_of_another_frame_that_shares_the_deleted_filepath(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: A KNOWN ROW AT THE SAME FILEPATH WITH A DIFFERENT KEY (`file`, `eso dpr tech`)
    organiser, _, _ = _archive_workspace(tmp_path, log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser._select_instrument()
    pd.DataFrame([OTHER_ROW]).to_sql("raw_frames", organiser.conn, index=False, if_exists="append")

    # ACT
    try:
        organiser._sync_raw_frames(skipSqlSync=True)
        techs = [row[0] for row in organiser.conn.execute('SELECT "eso dpr tech" FROM raw_frames').fetchall()]
    finally:
        organiser.close()

    # ASSERT
    assert techs == ["OTHER,TECH"]


def test_prepare_does_not_read_the_session_id_when_a_dated_archive_link_leaves_a_row_behind(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE
    organiser, root, frame = _archive_workspace(tmp_path, log)

    # ACT
    try:
        organiser.prepare(report=False)
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []
    assert frame.is_file()
    assert not (root / "x.fits").exists()


def test_sync_sql_table_to_directory_removes_a_missing_frame_row_before_any_session_exists(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: A ROW WHOSE FILE IS GONE, IN A WORKSPACE WITH NO SESSION YET
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES (?)", (f"./raw/{NIGHT}/gone.fits",))

    # ACT
    try:
        organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []


def test_sync_sql_table_to_directory_clears_the_sof_map_of_a_missing_frame_once_a_session_exists(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: A SESSION WHOSE SOF MAP HOLDS A SET THAT USES THE FRAME THAT IS GONE
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES (?)", (f"./raw/{NIGHT}/gone.fits",))
    organiser.conn.execute("CREATE TABLE sof_map_base (sof TEXT, filepath TEXT)")
    organiser.conn.executemany(
        "INSERT INTO sof_map_base VALUES (?, ?)",
        [("set.sof", f"./raw/{NIGHT}/gone.fits"), ("set.sof", f"./raw/{NIGHT}/kept.fits")],
    )

    # ACT
    try:
        organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")
        remainingRows = organiser.conn.execute("SELECT COUNT(*) FROM sof_map_base").fetchone()
    finally:
        organiser.close()

    # ASSERT
    assert remainingRows == (0,)


def test_sync_sql_table_to_directory_skips_the_sof_map_when_the_session_has_no_table(log: Any, tmp_path: Path) -> None:
    # ARRANGE: A STALE `.sessionid` NAMES A SESSION WHOSE DIRECTORY, AND SO WHOSE SOF MAP TABLE, DOES NOT EXIST
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "stale"
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES (?)", (f"./raw/{NIGHT}/gone.fits",))

    # ACT
    try:
        organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []


def test_sync_raw_frames_leaves_no_row_for_a_dated_archive_frame_whose_dpr_tech_is_missing(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE FRAME'S TECH READS `--`; THE COLUMN IS NOT NULL ON CONFLICT IGNORE, SO NO ROW IS STORED FOR IT
    organiser, _, _ = _archive_workspace(tmp_path, log)
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser._select_instrument()
    scan = organiser._create_directory_table

    def scan_without_tech(**kwargs: Any) -> tuple[Any, Any, Any]:
        rawFrames, fitsPaths, remainingFiles = scan(**kwargs)
        if rawFrames is not None:
            rawFrames["eso dpr tech"] = "--"
        return rawFrames, fitsPaths, remainingFiles

    monkeypatch.setattr(organiser, "_create_directory_table", scan_without_tech)

    # ACT
    try:
        organiser._sync_raw_frames(skipSqlSync=True)
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []


def test_sync_raw_frames_deletes_the_archive_row_and_holds_a_root_frame_whose_destination_is_taken(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: ONE ROOT LINK TO A DATED ARCHIVE FRAME, AND ONE ROOT FRAME WHOSE `raw/<date>/` DESTINATION IS TAKEN
    organiser, root, archivedFrame = _archive_workspace(tmp_path, log)
    heldFrame = harvestable_raw_fits(root / "y.fits")
    takenDestination = root / "raw" / NIGHT / "y.fits"
    takenDestination.parent.mkdir(parents=True)
    takenDestination.write_bytes(b"already here")
    organiser.conn, _ = organiser._get_or_create_db_connection()
    organiser._select_instrument()

    # ACT
    try:
        organiser._sync_raw_frames(skipSqlSync=True)
        rowPaths = _raw_frame_paths(organiser)
    finally:
        organiser.close()

    # ASSERT
    assert rowPaths == []
    assert archivedFrame.is_file()
    assert not (root / "x.fits").exists()
    assert heldFrame.is_file()
    assert takenDestination.read_bytes() == b"already here"
