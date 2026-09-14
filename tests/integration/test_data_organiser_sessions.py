"""Integration contracts for public data-organiser session workflows."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.factories import raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration


def test_session_create_initialises_sqlite_view_assets_and_workspace_links(
    tmp_path, log
) -> None:
    """Create an isolated session with its database view and root-facing assets."""
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    connection.execute(
        """INSERT INTO raw_frames (
            instrume, file, "mjd-obs", "mjd-date", "date-obs",
            "night start mjd", "eso dpr catg", "eso dpr type", "eso dpr tech",
            exptime, object, filepath
        ) VALUES ('SOXS', 'bias-1.fits', 60311.1, '2024-01-02',
            '2024-01-02T03:04:05', 60310, 'CALIB', 'BIAS', 'IMAGE', 0, 'BIAS',
            './raw/bias-1.fits')"""
    )
    sessionPath = Path(organiser.sessionsDir) / "science"
    sessionPath.mkdir()
    (sessionPath / "soxspipe.yaml").write_text("# synthetic settings\n")

    sessionId = organiser.session_create("science")

    assert sessionId == "science"
    assert (sessionPath / "soxspipe.yaml").is_file()
    assert {(sessionPath / name).is_dir() for name in ("sof", "qc", "reduced")} == {
        True
    }
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "science"
    assert connection.execute(
        "SELECT name FROM sqlite_master WHERE type = 'table' AND name = 'sof_map_science'"
    ).fetchone() == ("sof_map_science",)
    assert connection.execute(
        "SELECT name FROM sqlite_master WHERE type = 'view' AND name = 'sof_map'"
    ).fetchone() == ("sof_map",)
    assert connection.execute("SELECT * FROM sof_map").fetchall() == []
    assert connection.execute("PRAGMA table_info(product_frames)").fetchall()
    assert "status_science" in {
        row[1] for row in connection.execute("PRAGMA table_info(product_frames)")
    }
    for assetName in ("reduced", "qc", "sof", "soxspipe.yaml"):
        rootAsset = Path(organiser.rootDir) / assetName
        assert rootAsset.is_symlink()
        assert rootAsset.resolve() == sessionPath / assetName


def test_session_build_sof_files_creates_complete_bias_inventory(tmp_path, log) -> None:
    """Build a session's master-bias product and its usable SOF from raw rows."""
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    rawColumns = {
        row[1] for row in connection.execute("PRAGMA table_info(raw_frames)")
    }
    rawFrames = raw_frame_table()
    rawFrames.loc[:, rawFrames.columns.isin(rawColumns)].to_sql(
        "raw_frames", connection, if_exists="append", index=False
    )
    sessionPath = Path(organiser.sessionsDir) / "science"
    sessionPath.mkdir()
    (sessionPath / "soxspipe.yaml").write_text("# synthetic settings\n")
    organiser.session_create("science")
    organiser._select_instrument()

    organiser.build_sof_files()

    productRows = connection.execute(
        "SELECT recipe, complete FROM product_frames"
    ).fetchall()
    sofRows = connection.execute(
        "SELECT file, tag, complete FROM sof_map_science ORDER BY file"
    ).fetchall()
    sofPath = sessionPath / "sof" / "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"

    assert productRows == [("mbias", 1)]
    assert sofRows == [
        ("bias-1.fits", "BIAS_VIS", 1),
        ("bias-2.fits", "BIAS_VIS", 1),
    ]
    assert [line.split() for line in sofPath.read_text().splitlines()] == [
        ["./raw/2024-01-01/bias-1.fits", "BIAS_VIS"],
        ["./raw/2024-01-01/bias-2.fits", "BIAS_VIS"],
    ]
