"""Public data-organiser inventory listing contracts."""

from __future__ import annotations

import sqlite3

import pytest

from tests.factories import workspace_organiser

pytestmark = pytest.mark.integration


def _organiser_with_inventory(tmp_path, log):
    """Return an organiser backed by the minimal listing schema."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute(
        '''CREATE TABLE raw_frame_sets (
            "eso seq arm" TEXT, "night start date" TEXT, "night start mjd" REAL,
            "eso obs name" TEXT, "eso obs id" TEXT, "mjd-obs" REAL, slit TEXT,
            sof TEXT, complete INTEGER, recipe TEXT, "eso dpr type" TEXT
        )'''
    )
    organiser.conn.execute(
        '''INSERT INTO raw_frame_sets VALUES (
            'VIS', '2024-01-02', 60311.0, 'Synthetic target', '42', 60311.2,
            'SLIT1.0', 'science.sof', 1, 'stare_obj', 'OBJECT')'''
    )
    organiser.conn.execute(
        '''CREATE TABLE product_frames (sof TEXT, complete INTEGER, file TEXT)'''
    )
    organiser.conn.executemany(
        "INSERT INTO product_frames VALUES (?, ?, ?)",
        [
            ("science.sof", 1, "science-product.fits"),
            ("calibration.sof", 1, "master.fits"),
        ],
    )
    organiser.conn.execute(
        '''CREATE TABLE sof_map_base (sof TEXT, file TEXT, filepath TEXT)'''
    )
    organiser.conn.executemany(
        "INSERT INTO sof_map_base VALUES (?, ?, ?)",
        [
            ("science.sof", "master.fits", "./reduced/master.fits"),
            ("calibration.sof", "master.fits", "./raw/2024-01-02/raw-1.fits"),
            ("calibration.sof", "raw-1.fits", "./raw/2024-01-02/raw-1.fits"),
        ],
    )
    organiser.conn.execute(
        "CREATE VIEW sof_map AS SELECT * FROM sof_map_base"
    )
    return organiser


def test_list_obs_and_sofs_return_sorted_complete_science_inventory(tmp_path, log) -> None:
    """Public listing methods return only complete science rows and their fields."""
    organiser = _organiser_with_inventory(tmp_path, log)

    observations = organiser.list_obs()
    sofs = organiser.list_sofs()

    assert observations.to_dict(orient="records") == [
        {
            "eso seq arm": "VIS",
            "night start date": "2024-01-02",
            "night start mjd": 60311.0,
            "eso obs name": "Synthetic target",
            "eso obs id": "42",
        }
    ]
    assert sofs.to_dict(orient="records") == [
        {
            "eso seq arm": "VIS",
            "recipe": "STARE",
            "night start date": "2024-01-02",
            "mjd-obs": 60311.2,
            "eso obs name": "Synthetic target",
            "eso obs id": "42",
            "slit": "SLIT1.0",
            "sof": "science.sof",
        }
    ]


def test_list_raw_recursively_traces_raw_members_and_deduplicates_paths(tmp_path, log) -> None:
    """Raw inventory follows dependent SOFs and returns each raw file once."""
    organiser = _organiser_with_inventory(tmp_path, log)

    paths, table = organiser.list_raw("science.sof")

    assert paths == ["./raw/2024-01-02/raw-1.fits"]
    assert table[["sof", "file", "filepath"]].to_dict(orient="records") == [
        {
            "sof": "calibration.sof",
            "file": "master.fits",
            "filepath": "./raw/2024-01-02/raw-1.fits",
        },
        {
            "sof": "calibration.sof",
            "file": "raw-1.fits",
            "filepath": "./raw/2024-01-02/raw-1.fits",
        },
    ]


def test_list_raw_binds_a_hostile_sof_file_as_a_parameter(tmp_path, log) -> None:
    """An injection-shaped `sofFile` matches no real row rather than every row.

    `list_raw` interpolated `sofFile` into the base query text. A value ending
    `' OR '1'='1' --` closed the string literal early, matched every row in
    `product_frames`, and propagated through the five-level recursive nesting.
    Bound as a parameter, it matches nothing.
    """
    organiser = _organiser_with_inventory(tmp_path, log)
    hostileSofFile = "science.sof' OR '1'='1' --"

    paths, table = organiser.list_raw(hostileSofFile)

    assert paths == []
    assert table.empty
