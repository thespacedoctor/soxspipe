"""Integration contracts for the data organiser's temporary SQLite state."""

from __future__ import annotations

import os
import sqlite3
import struct

import pandas as pd
import pytest

from soxspipe.recipes import base_recipe
from tests.factories import (
    product_table,
    qc_table,
    raw_frame_table,
    raw_group_table,
    workspace_organiser,
)

pytestmark = pytest.mark.integration


def _write_corrupt_but_openable_database(path: str | os.PathLike[str]) -> None:
    """Write a real SQLite file whose `PRAGMA integrity_check` returns non-`ok`
    rows without SQLite refusing to open it or raising while executing the
    check.
    """
    connection = sqlite3.connect(path)
    connection.execute("PRAGMA page_size=4096")
    connection.execute("CREATE TABLE t (id INTEGER PRIMARY KEY, val TEXT)")
    connection.executemany(
        "INSERT INTO t (val) VALUES (?)",
        [(f"row-{i}-" + "x" * 50,) for i in range(500)],
    )
    connection.commit()
    connection.execute("VACUUM")
    connection.close()

    # CORRUPT THE "NUMBER OF FREELIST PAGES" FIELD IN THE FILE HEADER (BYTES
    # 36-39, BIG-ENDIAN) SO INTEGRITY_CHECK REPORTS A MISMATCH AS DATA ROWS
    # RATHER THAN RAISING WHILE READING THE FILE.
    with open(path, "r+b") as databaseFile:
        databaseFile.seek(36)
        databaseFile.write(struct.pack(">I", 999))


def test_database_connection_creates_and_reuses_valid_workspace_database(
    tmp_path, log
) -> None:
    """Copy the schema template once and keep a healthy connection reusable."""
    organiser = workspace_organiser(tmp_path, log=log)

    connection, wasReset = organiser._get_or_create_db_connection()
    organiser.conn = connection
    reusedConnection, wasReusedReset = organiser._get_or_create_db_connection()

    assert wasReset is False
    assert wasReusedReset is False
    assert reusedConnection is connection
    assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    tableNames = {
        row[0]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }
    assert {"raw_frames", "product_frames", "z_sof_map"} <= tableNames
    connection.close()


def test_database_connection_rebuilds_a_corrupt_but_openable_database(
    tmp_path, log, monkeypatch
) -> None:
    """Detect a corrupt-but-openable database via `PRAGMA integrity_check`
    and rebuild it through the same recovery path used for a connection
    failure, instead of silently reusing it."""
    monkeypatch.setattr("time.sleep", lambda seconds: None)
    organiser = workspace_organiser(tmp_path, log=log)
    _write_corrupt_but_openable_database(organiser.rootDbPath)

    # CONFIRM THE FIXTURE'S OWN PRECONDITION: THE FILE MUST STAY OPENABLE
    # AND MUST MAKE `PRAGMA INTEGRITY_CHECK` REPORT ROWS OTHER THAN
    # `[("OK",)]`. WITHOUT THIS, A SQLITE BUILD THAT INSTEAD RAISED ON
    # OPEN/CHECK WOULD STILL ROUTE INTO THE PRE-EXISTING BROAD
    # `EXCEPT EXCEPTION:` RECOVERY PATH, AND THIS TEST WOULD PASS WITHOUT
    # EVER EXERCISING THE FETCHALL-INSPECTION CODE IT IS MEANT TO PIN.
    precondingConnection = sqlite3.connect(organiser.rootDbPath)
    precondingRows = precondingConnection.execute(
        "PRAGMA integrity_check"
    ).fetchall()
    precondingConnection.close()
    assert precondingRows != [("ok",)]

    # `PREPARE(REFRESH=TRUE)` IS THE EXISTING RECOVERY PATH THIS FIX MUST
    # ROUTE INTO; IT ALSO DOES A FULL WORKSPACE SCAN (FITS INDEXING,
    # SESSION BOOTSTRAP, AND SO ON) THAT IS UNRELATED TO THIS DEFECT AND
    # REQUIRES A MUCH BIGGER FIXTURE, SO REPLACE IT WITH A STAND-IN THAT
    # PERFORMS THE SAME FILE-REPLACEMENT STEP `PREPARE(REFRESH=TRUE)`
    # ITSELF DOES, THEN RE-RUNS THE REAL, UNMOCKED `_GET_OR_CREATE_DB_CONNECTION`
    # TO REBUILD THE CONNECTION FROM THE SCHEMA TEMPLATE - EXACTLY WHAT
    # `PREPARE` DOES INTERNALLY.
    prepareRefreshValues: tuple[bool, ...] = ()

    def _rebuild_from_template(*, refresh: bool = False, report: bool = True) -> None:
        """Stand in for `prepare(refresh=True)`'s file-replacement step, then
        rebuild the connection via the real `_get_or_create_db_connection`.
        """
        nonlocal prepareRefreshValues
        prepareRefreshValues = (*prepareRefreshValues, refresh)
        if refresh and os.path.exists(organiser.rootDbPath):
            os.remove(organiser.rootDbPath)
        organiser.conn, _ = organiser._get_or_create_db_connection()

    monkeypatch.setattr(organiser, "prepare", _rebuild_from_template)

    connection = None
    try:
        connection, wasReset = organiser._get_or_create_db_connection()

        assert prepareRefreshValues == (True,)
        assert wasReset is True
        assert connection.execute("PRAGMA integrity_check").fetchone() == ("ok",)
        tableNames = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type = 'table'"
            )
        }
        assert {"raw_frames", "product_frames", "z_sof_map"} <= tableNames
    finally:
        # `_REBUILD_FROM_TEMPLATE` (THE `PREPARE` STAND-IN) SETS
        # `ORGANISER.CONN` VIA A NESTED CALL TO
        # `_GET_OR_CREATE_DB_CONNECTION`, WHICH IS A GENUINELY SEPARATE
        # CONNECTION OBJECT FROM THE ONE THE OUTER CALL RETURNS - CLOSE
        # BOTH ON EVERY EXIT PATH.
        if connection is not None:
            connection.close()
        if organiser.conn is not None and organiser.conn is not connection:
            organiser.conn.close()


def test_dataframe_write_normalises_sentinels_and_replaces_existing_rows(
    tmp_path, log
) -> None:
    """Persist missing-value sentinels as SQL NULL and honor replacement."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    firstRows = pd.DataFrame(
        {"file": ["old.fits", "missing.fits"], "value": [1.0, -99.99]}
    )
    replacement = pd.DataFrame({"file": ["new.fits"], "value": [2.0]})

    organiser._dataframe_to_sqlite(firstRows, "safe_rows")
    initialRows = organiser.conn.execute(
        "SELECT file, value FROM safe_rows ORDER BY file"
    ).fetchall()
    organiser._dataframe_to_sqlite(replacement, "safe_rows", replace=True)

    rows = organiser.conn.execute(
        "SELECT file, value FROM safe_rows ORDER BY file"
    ).fetchall()
    assert initialRows == [("missing.fits", None), ("old.fits", "1.0")]
    assert rows == [("new.fits", "2.0")]


def test_recipe_report_replaces_only_current_sof_qc_rows_in_workspace_database(
    tmp_path, log
) -> None:
    """Persist current-SOF QC rows without deleting another reduction's rows."""
    databasePath = tmp_path / "workspace.db"
    connection = sqlite3.connect(databasePath)
    currentQc = qc_table().assign(qc_order="-1", qc_flag="pass", sof_name="current.sof")
    databaseColumns = [column for column in currentQc.columns if column != "to_header"]
    existingRows = pd.concat(
        [
            currentQc.assign(qc_value=999.0),
            currentQc.assign(sof_name="other.sof", qc_value=7.5),
        ],
        ignore_index=True,
    )
    existingRows[databaseColumns].to_sql(
        "quality_control", connection, index=False
    )
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.verbose = True
    recipe.conn = connection
    recipe.sofName = "current"
    recipe.dateObs = "2024-01-02T03:04:05.678"
    recipe.recipeName = "soxs-mbias"
    recipe.inst = "XSHOOTER"
    recipe.recipeSettings = {}
    recipe.qc = qc_table()
    recipe.products = product_table()

    reported = recipe.report_output()

    rows = connection.execute(
        "SELECT sof_name, qc_name, qc_value FROM quality_control ORDER BY sof_name"
    ).fetchall()
    assert rows == [("current.sof", "RON", 3.2), ("other.sof", "RON", 7.5)]
    assert reported["sof_name"].tolist() == ["current.sof"]
    connection.close()


def test_raw_groups_add_only_selected_members_to_sof_map(tmp_path, log) -> None:
    """Explode selected group members and flag their raw rows as processed."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute(
        "CREATE TABLE raw_frames (file TEXT PRIMARY KEY, processed INTEGER)"
    )
    organiser.conn.executemany(
        "INSERT INTO raw_frames VALUES (?, 0)",
        [("bias-1.fits",), ("bias-2.fits",), ("other.fits",)],
    )
    rawGroups = pd.concat(
        [
            raw_group_table(),
            raw_group_table().assign(
                sof="UNSELECTED.sof",
                filepaths=[["./raw/2024-01-01/other.fits"]],
            ),
        ],
        ignore_index=True,
    )
    selectedSof = rawGroups.loc[0, "sof"]

    organiser.raw_frames_to_sof_map(rawGroups, [selectedSof])

    sofRows = organiser.conn.execute(
        "SELECT file, tag, sof, filepath, complete FROM sof_map_base ORDER BY file"
    ).fetchall()
    assert sofRows == [
        (
            "bias-1.fits",
            "BIAS_VIS",
            selectedSof,
            "./raw/2024-01-01/bias-1.fits",
            1,
        ),
        (
            "bias-2.fits",
            "BIAS_VIS",
            selectedSof,
            "./raw/2024-01-01/bias-2.fits",
            1,
        ),
    ]
    processedRows = organiser.conn.execute(
        "SELECT file, processed FROM raw_frames ORDER BY file"
    ).fetchall()
    assert processedRows == [
        ("bias-1.fits", 1),
        ("bias-2.fits", 1),
        ("other.fits", 0),
    ]


def test_raw_groups_skip_sof_mapping_when_no_groups_exist(tmp_path, log) -> None:
    """Leave SQLite untouched when there are no raw groups to map."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.sessionId = "base"
    organiser.conn = sqlite3.connect(":memory:")

    result = organiser.raw_frames_to_sof_map(pd.DataFrame(), [])

    assert result is None
    assert (
        organiser.conn.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        ).fetchall()
        == []
    )


def test_raw_frame_inventory_groups_complete_bias_set_and_names_sof(
    tmp_path, log
) -> None:
    """Read valid raw inventory and build one complete, deterministically named set."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.instrument = "SOXS"
    organiser.conn = sqlite3.connect(":memory:")
    raw_frame_table().to_sql("raw_frames_valid", organiser.conn, index=False)

    rawFrames, rawGroups = organiser.get_raw_frames_and_groups(
        ttype="BIAS",
        recipe="mbias",
        recipeOrder=0,
        filterName="bias",
        unprocessedOnly=True,
    )

    assert rawFrames["file"].tolist() == ["bias-1.fits", "bias-2.fits"]
    assert len(rawGroups.index) == 1
    assert rawGroups.loc[0, "counts"] == 2
    assert rawGroups.loc[0, "complete"] == 1
    assert rawGroups.loc[0, "recipe_order"] == 0
    assert rawGroups.loc[0, "sof"] == ("20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof")


def test_raw_frame_inventory_binds_a_hostile_arm_as_a_parameter(tmp_path, log) -> None:
    """A hostile `arm` filter matches no real row in a real database, rather than widening the query.

    DY-254: `ttype`, `arm` and `tech` were interpolated directly into the
    `raw_frames_valid` filter text in `get_raw_frames_and_groups`.
    """
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.instrument = "SOXS"
    organiser.conn = sqlite3.connect(":memory:")
    frames = raw_frame_table()
    frames.loc[frames.index[1], "eso seq arm"] = "NIR"
    frames.to_sql("raw_frames_valid", organiser.conn, index=False)

    # IF STILL INTERPOLATED, THIS CLOSES THE STRING LITERAL, INJECTS AN
    # ALWAYS-TRUE CONDITION, AND COMMENTS OUT THE TRAILING QUOTE, MATCHING
    # BOTH THE VIS AND NIR ROWS INSTEAD OF NEITHER.
    hostileArm = "VIS' OR '1'='1' -- "

    rawFrames, rawGroups = organiser.get_raw_frames_and_groups(
        arm=hostileArm,
        recipe="mbias",
        recipeOrder=0,
        filterName="bias",
        unprocessedOnly=True,
    )

    assert rawFrames.empty
    assert rawGroups.empty
