"""Integration contracts for public data-organiser session workflows."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.factories import raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration


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


def test_session_create_initialises_sqlite_view_assets_and_workspace_links(
    tmp_path, log
) -> None:
    """Create an isolated session with its database view and root-facing assets."""
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_complete_raw_frame(connection)
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


def test_prepare_refresh_restores_every_session_database_object(
    tmp_path, log, monkeypatch
) -> None:
    """Rebuild a named-session workspace without leaving its schema unusable."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_complete_raw_frame(connection)

    organiser.session_create("science")
    organiser.session_create("archive")
    Path(organiser.sessionIdFile).write_text("science", encoding="utf-8")
    organiser.sessionId = "science"
    organiser.sessionPath = str(Path(organiser.sessionsDir) / "science")
    (Path(organiser.sessionsDir) / "not-a-session").mkdir()
    (Path(organiser.sessionsDir) / "README.txt").write_text(
        "not a session", encoding="utf-8"
    )

    monkeypatch.setattr(organiser, "_select_instrument", lambda: None)
    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(
        organiser, "_sync_raw_frames", lambda: _insert_complete_raw_frame(organiser.conn)
    )
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_flag_files_to_ignore", lambda: None)
    monkeypatch.setattr(organiser, "build_sof_files", lambda **_: None)

    organiser.prepare(refresh=True, report=False)

    columns = {
        row[1] for row in organiser.conn.execute("PRAGMA table_info(product_frames)")
    }
    assert {"status_science", "status_archive"} <= columns
    tables = {
        row[0]
        for row in organiser.conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }
    assert {"sof_map_science", "sof_map_archive"} <= tables
    viewSql = organiser.conn.execute(
        "SELECT sql FROM sqlite_master WHERE type = 'view' AND name = 'sof_map'"
    ).fetchone()[0]
    assert "sof_map_science" in viewSql
    assert "Skipping invalid session directory" in " ".join(
        message for _, message in log.messages
    )

    organiser.prepare(report=False)


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


def _sof_map_view_sql(connection) -> str:
    """Return the stored definition of the shared `sof_map` view."""
    return connection.execute("SELECT sql FROM sqlite_master WHERE type = 'view' AND name = 'sof_map'").fetchone()[0]


@pytest.fixture
def two_session_workspace(tmp_path, log):
    """Create session `alpha` and then session `beta`, which leaves `beta` active."""
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_complete_raw_frame(connection)
    organiser.session_create("alpha")
    organiser.session_create("beta")
    return organiser


def test_session_switch_points_the_sof_map_view_at_the_session_switched_to(two_session_workspace, capsys) -> None:
    # ARRANGE
    organiser = two_session_workspace
    assert _sof_map_view_sql(organiser.conn) == "CREATE VIEW sof_map as select * from sof_map_beta"
    organiser.conn.execute("INSERT INTO sof_map_alpha (filepath, tag, sof) VALUES ('a.fits', 'BIAS', 'alpha.sof')")

    # ACT
    organiser.session_switch("alpha")

    # ASSERT
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "alpha"
    assert _sof_map_view_sql(organiser.conn) == "CREATE VIEW sof_map as select * from sof_map_alpha"
    assert organiser.conn.execute("SELECT sof FROM sof_map").fetchall() == [("alpha.sof",)]
    assert "Session successfully switched to 'alpha'." in capsys.readouterr().out


def test_session_switch_back_points_the_sof_map_view_at_the_original_session(two_session_workspace) -> None:
    # ARRANGE
    organiser = two_session_workspace
    organiser.session_switch("alpha")

    # ACT
    organiser.session_switch("beta")

    # ASSERT
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "beta"
    assert _sof_map_view_sql(organiser.conn) == "CREATE VIEW sof_map as select * from sof_map_beta"


@pytest.mark.parametrize(
    ("sessionId", "expectedMessage"),
    [
        ("beta", "Session 'beta' is already in use."),
        ("gamma", "There is no session with the ID 'gamma'. List existing sessions with `soxspipe session ls`."),
    ],
)
def test_refused_session_switch_leaves_the_sof_map_view_unchanged(
    two_session_workspace, capsys, sessionId, expectedMessage
) -> None:
    # ARRANGE
    organiser = two_session_workspace
    viewBefore = _sof_map_view_sql(organiser.conn)
    capsys.readouterr()

    # ACT
    result = organiser.session_switch(sessionId)

    # ASSERT
    assert result is None
    assert _sof_map_view_sql(organiser.conn) == viewBefore
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "beta"
    assert capsys.readouterr().out.strip() == expectedMessage


def test_the_sof_map_view_is_created_in_exactly_one_place_in_the_package() -> None:
    # ARRANGE
    packageDir = Path(__file__).resolve().parents[2] / "soxspipe"

    # ACT
    creators = [
        f"{path.relative_to(packageDir)}:{lineNumber}"
        for path in sorted(packageDir.rglob("*.py"))
        if "tests" not in path.relative_to(packageDir).parts
        for lineNumber, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if "create view sof_map" in line.lower()
    ]

    # ASSERT
    assert len(creators) == 1, creators
