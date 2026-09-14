"""Integration contracts for isolated data-organiser workspaces."""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.data_organiser import (
    _UnsafePathError,
    _harvest_fits_headers,
)
from tests.factories import raw_fits, workspace_organiser

pytestmark = pytest.mark.integration


def _organiser_header_keywords(log) -> tuple[object, list[str]]:
    kw = keyword_lookup(log=log, instrument="SOXS").get
    aliases = [
        "MJDOBS",
        "DATE_OBS",
        "SEQ_ARM",
        "DPR_CATG",
        "DPR_TECH",
        "DPR_TYPE",
        "PRO_CATG",
        "PRO_TECH",
        "PRO_TYPE",
        "EXPTIME",
        "EXPTIME2",
        "WIN_BINX",
        "WIN_BINY",
        "DET_READ_SPEED",
        "TPL_ID",
        "ACFW_ID",
        "DET",
        "ABSROT",
    ]
    return kw, ["file", *(kw(alias).lower() for alias in aliases)]


def _write_harvestable_fits(destination: Path, *, includeDprType: bool = True) -> Path:
    raw_fits(destination)
    with fits.open(destination, mode="update") as hdus:
        header = hdus[0].header
        header["MJD-OBS"] = 60311.75
        header["ESO DET3 EXPO TIME"] = 8.0
        header["ESO DET BINX"] = 1
        header["ESO DET BINY"] = 2
        header["ESO TPL ID"] = "SOXS_cal_bias"
        header["ESO INS ACFW ID"] = "g"
        header["ESO INS VISE NAME"] = "SLIT_1.0"
        header["ESO DET3 CAM NAME"] = "VIS"
        header["ESO ADA ABSROT END"] = 12.5
        if not includeDprType:
            del header["ESO DPR TYPE"]
    return destination


def test_header_harvest_filters_incomplete_fits_and_derives_night_metadata(
    tmp_path, log, capsys
) -> None:
    """Use real FITS headers and ignore frames missing mandatory DPR metadata."""
    validPath = _write_harvestable_fits(tmp_path / "valid.fits")
    invalidPath = _write_harvestable_fits(
        tmp_path / "invalid.fits", includeDprType=False
    )
    kw, keywords = _organiser_header_keywords(log)

    harvested = _harvest_fits_headers(
        batch=[str(validPath), str(invalidPath)],
        log=log,
        pathToDirectory=str(tmp_path),
        keywords=keywords,
        filterKeys=[
            "mjd-obs",
            "eso dpr catg",
            "eso dpr tech",
            "eso dpr type",
            "exptime",
        ],
        instrument="SOXS",
        kw=kw,
    )

    assert len(harvested.index) == 1
    assert Path(harvested.loc[0, "file"]).name == "valid.fits"
    assert harvested.loc[0, "night start date"] == "2024-01-02"
    assert harvested.loc[0, "mjd-date"] == "2024-01-02"
    assert harvested.loc[0, "binning"] == "1x2"
    assert harvested.loc[0, "rospeed"] == 2
    assert harvested.loc[0, "filepath"] == "--"
    assert "missing DPR keywords" in capsys.readouterr().out


def test_directory_table_returns_empty_contract_when_no_fits_exist(
    tmp_path, log
) -> None:
    """Return the established all-None result for an empty workspace scan."""
    organiser = workspace_organiser(tmp_path, log=log)

    result = organiser._create_directory_table(
        str(Path(organiser.rootDir)), organiser.filterKeywords
    )

    assert result == (None, None, None)


def test_directory_table_rejects_mixed_instruments(tmp_path, log) -> None:
    """Fail before indexing when one directory mixes SOXS and XSH frames."""
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)
    raw_fits(rootPath / "soxs.fits", instrument="soxs")
    raw_fits(rootPath / "xsh.fits", instrument="xsh")

    with pytest.raises(AssertionError):
        organiser._create_directory_table(str(rootPath), organiser.filterKeywords)

    assert any(
        level == "error" and "mix of instruments" in message
        for level, message in log.messages
    )


def test_fits_detection_checks_workspace_root_and_raw_tree(tmp_path, log) -> None:
    """Detect standard lowercase FITS files in either supported location."""
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)

    assert organiser._fits_files_exist() is False
    raw_fits(rootPath / "root.fits")
    assert organiser._fits_files_exist() is True
    (rootPath / "root.fits").unlink()
    nestedPath = Path(organiser.rawDir) / "2024-01-01"
    nestedPath.mkdir()
    raw_fits(nestedPath / "nested.fits")
    assert organiser._fits_files_exist() is True


def test_prepare_refresh_removes_only_workspace_database_before_empty_inventory_exit(
    tmp_path, log, monkeypatch
) -> None:
    """A requested refresh clears local SQLite sidecars before reporting no data."""
    organiser = workspace_organiser(tmp_path, log=log)
    databasePaths = [
        Path(organiser.rootDbPath),
        Path(f"{organiser.rootDbPath}-shm"),
        Path(f"{organiser.rootDbPath}-wal"),
    ]
    for databasePath in databasePaths:
        databasePath.write_text("stale", encoding="utf-8")
    monkeypatch.setattr(organiser, "_select_instrument", lambda: None)
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: False)

    with pytest.raises(SystemExit):
        organiser.prepare(refresh=True)

    assert not any(databasePath.exists() for databasePath in databasePaths)


def test_instrument_selection_normalises_xshooter_and_loads_sof_map(
    tmp_path, log
) -> None:
    """Select the XSH configuration from the instrument stored in SQLite."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (instrume TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES ('XSHOOTER')")

    organiser._select_instrument()

    assert organiser.instrument == "XSH"
    assert organiser.kw("SEQ_ARM") == "ESO SEQ ARM"
    assert organiser.sofMapLookup["bias_frame"]["recipe"] == "mbias"


def test_session_list_and_switch_retarget_workspace_symlinks(tmp_path, log) -> None:
    """List sessions deterministically and retarget all root assets on switch."""
    organiser = workspace_organiser(tmp_path, log=log)
    sessionsPath = Path(organiser.sessionsDir)
    for sessionName in ("base", "science"):
        sessionPath = sessionsPath / sessionName
        sessionPath.mkdir()
        for directoryName in ("reduced", "qc", "sof"):
            (sessionPath / directoryName).mkdir()
        for fileName in ("soxspipe.yaml", "soxspipe.log"):
            (sessionPath / fileName).write_text(sessionName, encoding="utf-8")
    Path(organiser.sessionIdFile).write_text("base", encoding="utf-8")
    organiser.sessionPath = str(sessionsPath / "base")
    organiser._symlink_session_assets_to_workspace_root()

    currentSession, allSessions = organiser.session_list(silent=True)
    organiser.session_switch("science")

    assert currentSession == "base"
    assert allSessions == ["base", "science"]
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "science"
    for assetName in ("reduced", "qc", "sof", "soxspipe.yaml", "soxspipe.log"):
        rootAsset = Path(organiser.rootDir) / assetName
        assert rootAsset.is_symlink()
        assert Path(os.readlink(rootAsset)) == sessionsPath / "science" / assetName


def test_session_list_and_switch_preserve_missing_and_active_session_contracts(
    tmp_path, log, capsys
) -> None:
    """Session commands report absent, unknown, and already-active sessions safely."""
    organiser = workspace_organiser(tmp_path, log=log)

    assert organiser.session_list(silent=True) == (None, None)

    sessionPath = Path(organiser.sessionsDir) / "base"
    sessionPath.mkdir()
    Path(organiser.sessionIdFile).write_text("base", encoding="utf-8")
    organiser.sessionPath = str(sessionPath)

    assert organiser.session_switch("base") is None
    assert organiser.session_switch("missing") is None

    output = capsys.readouterr().out
    assert "already in use" in output
    assert "There is no session" in output


def test_move_misc_files_keeps_workspace_metadata_and_moves_unrelated_files(
    tmp_path, log
) -> None:
    """Workspace preparation archives unrelated root files without moving metadata."""
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)
    (rootPath / "notes.txt").write_text("move", encoding="utf-8")
    (rootPath / "README.txt").write_text("keep", encoding="utf-8")
    (rootPath / "soxspipe-helper.txt").write_text("keep", encoding="utf-8")
    (rootPath / "settings.yaml").write_text("keep", encoding="utf-8")

    organiser._move_misc_files()

    assert (Path(organiser.miscDir) / "notes.txt").read_text(encoding="utf-8") == "move"
    assert (rootPath / "README.txt").is_file()
    assert (rootPath / "soxspipe-helper.txt").is_file()
    assert (rootPath / "settings.yaml").is_file()


@pytest.mark.parametrize("failure", [True, False])
def test_session_refresh_rebuilds_the_active_session_inventory(
    tmp_path, log, monkeypatch, failure: bool
) -> None:
    """Refreshing an active session rebuilds its SOFs for either status transition."""
    organiser = workspace_organiser(tmp_path, log=log)
    sessionPath = Path(organiser.sessionsDir) / "base"
    sessionPath.mkdir()
    Path(organiser.sessionIdFile).write_text("base", encoding="utf-8")
    connection = sqlite3.connect(":memory:")
    calls: list[str] = []

    monkeypatch.setattr(
        organiser,
        "_get_or_create_db_connection",
        lambda: (connection, False),
    )
    monkeypatch.setattr(organiser, "_select_instrument", lambda: calls.append("instrument"))
    monkeypatch.setattr(organiser, "build_sof_files", lambda: calls.append("sofs"))

    reset = organiser.session_refresh(silent=True, failure=failure)

    assert reset is False
    assert organiser.sessionId == "base"
    assert organiser.sessionPath == str(sessionPath)
    assert calls == ["instrument", "sofs"]


def test_write_sof_files_writes_only_complete_new_inventories(tmp_path, log) -> None:
    """Write complete database inventories while preserving existing SOFs."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.sessionPath = str(Path(organiser.sessionsDir) / "base")
    Path(organiser.sessionPath).mkdir()
    organiser.conn = sqlite3.connect(":memory:")
    pd.DataFrame(
        [
            {
                "sof": "new.sof",
                "filepath": "./raw/bias-1.fits",
                "tag": "BIAS_VIS",
                "complete": 1,
            },
            {
                "sof": "new.sof",
                "filepath": "./raw/bias-2.fits",
                "tag": "BIAS_VIS",
                "complete": 1,
            },
            {
                "sof": "incomplete.sof",
                "filepath": "./raw/dark.fits",
                "tag": "DARK_VIS",
                "complete": 0,
            },
            {
                "sof": "existing.sof",
                "filepath": "./raw/old.fits",
                "tag": "BIAS_VIS",
                "complete": 1,
            },
        ]
    ).to_sql("sof_map_base", organiser.conn, index=False)
    sofPath = Path(organiser.sessionPath) / "sof"
    sofPath.mkdir()
    existingPath = sofPath / "existing.sof"
    existingPath.write_text("maintainer content\n", encoding="utf-8")

    organiser._write_sof_files()

    rows = [line.split() for line in (sofPath / "new.sof").read_text().splitlines()]
    assert rows == [
        ["./raw/bias-1.fits", "BIAS_VIS"],
        ["./raw/bias-2.fits", "BIAS_VIS"],
    ]
    assert not (sofPath / "incomplete.sof").exists()
    assert existingPath.read_text(encoding="utf-8") == "maintainer content\n"


def test_directory_sync_removes_stale_database_and_sof_rows(tmp_path, log) -> None:
    """Drop database inventory for a raw file that no longer exists on disk."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES ('./raw/missing.fits')")
    organiser.conn.execute("CREATE TABLE sof_map_base (sof TEXT, filepath TEXT)")
    organiser.conn.execute(
        "INSERT INTO sof_map_base VALUES ('missing.sof', './raw/missing.fits')"
    )

    organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")

    assert organiser.conn.execute("SELECT * FROM raw_frames").fetchall() == []
    assert organiser.conn.execute("SELECT * FROM sof_map_base").fetchall() == []


def test_directory_sync_removes_a_stale_absolute_database_path(tmp_path, log) -> None:
    """Remove stale legacy absolute paths from both database inventories."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")
    stalePath = str(Path(organiser.rawDir) / "missing.fits")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES (?)", (stalePath,))
    organiser.conn.execute("CREATE TABLE sof_map_base (sof TEXT, filepath TEXT)")
    organiser.conn.execute(
        "INSERT INTO sof_map_base VALUES (?, ?)", ("missing.sof", stalePath)
    )

    organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")

    assert organiser.conn.execute("SELECT * FROM raw_frames").fetchall() == []
    assert organiser.conn.execute("SELECT * FROM sof_map_base").fetchall() == []


def test_directory_sync_rejects_database_paths_outside_the_workspace(tmp_path, log) -> None:
    """Reject a database path that escapes the workspace ownership boundary."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute("CREATE TABLE raw_frames (filepath TEXT)")
    organiser.conn.execute("INSERT INTO raw_frames VALUES ('../../outside.fits')")

    with pytest.raises(_UnsafePathError, match="database filepath"):
        organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames")


def test_directory_sync_rejects_an_unknown_sql_table_name(tmp_path, log) -> None:
    """Reject table names before they can be interpolated into a query."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")

    with pytest.raises(ValueError, match="raw_frames"):
        organiser._sync_sql_table_to_directory(organiser.rawDir, "raw_frames; DROP TABLE raw_frames")


def test_raw_frame_sync_indexes_and_moves_a_new_root_frame(tmp_path, log, monkeypatch) -> None:
    """A discovered frame is indexed once and moved into its dated raw directory."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    sourcePath = Path(organiser.rootDir) / "new-frame.fits"
    sourcePath.write_text("synthetic frame", encoding="utf-8")
    destinationPath = Path(organiser.rawDir) / "2024-01-02" / sourcePath.name
    frames = pd.DataFrame(
        [
            {
                "file": sourcePath.name,
                "eso dpr tech": "ECHELLE,SLIT,STARE",
                "mjd-date": "2024-01-02",
                "filepath": str(destinationPath),
            }
        ]
    )
    pd.DataFrame(columns=frames.columns).to_sql("raw_frames", organiser.conn, index=False)
    scans = iter([(frames.copy(), [str(sourcePath)], 1), (None, None, 0)])

    monkeypatch.setattr(
        organiser,
        "_create_directory_table",
        lambda **_: next(scans),
    )
    monkeypatch.setattr(
        organiser,
        "_populate_raw_frames_extra_columns",
        lambda table: table.copy(),
    )

    organiser._sync_raw_frames(skipSqlSync=True)

    assert not sourcePath.exists()
    assert destinationPath.read_text(encoding="utf-8") == "synthetic frame"
    indexed = pd.read_sql("SELECT file, filepath FROM raw_frames", organiser.conn)
    assert indexed.to_dict("records") == [
        {"file": "new-frame.fits", "filepath": str(destinationPath)}
    ]


def test_prepare_indexes_a_raw_frame_and_creates_a_base_session(
    tmp_path, log, monkeypatch
) -> None:
    """Prepare a disposable workspace through its public setup workflow."""
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)
    rawPath = _write_harvestable_fits(rootPath / "bias.fits")
    monkeypatch.chdir(rootPath)

    organiser.prepare(report=False)

    assert not rawPath.exists()
    indexedFrames = pd.read_sql("SELECT file, filepath FROM raw_frames", organiser.conn)
    assert indexedFrames.to_dict("records") == [
        {"file": "bias.fits", "filepath": "./raw/2024-01-02/bias.fits"}
    ]
    assert (rootPath / "sessions" / "base" / "soxspipe.yaml").is_file()
    assert (rootPath / "sof").is_symlink()


def test_prepare_report_describes_the_completed_isolated_workspace(
    tmp_path, log, monkeypatch, capsys
) -> None:
    """A second public preparation run reports the established workspace inventory."""
    organiser = workspace_organiser(tmp_path, log=log)
    rootPath = Path(organiser.rootDir)
    _write_harvestable_fits(rootPath / "bias.fits")
    monkeypatch.chdir(rootPath)

    organiser.prepare(report=False)
    organiser.prepare(report=True)

    output = capsys.readouterr().out
    assert "WORKSPACE FOR HAS BEEN PREPARED" in output
    assert "`misc/`: a lost-and-found archive" in output
    assert "`sessions/`: directory of data-reduction sessions" in output
