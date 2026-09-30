"""Unit contracts for grouping and product prediction in the data organiser."""

from __future__ import annotations

import sqlite3

import pandas as pd
import pytest

from tests.factories import raw_group_table, workspace_organiser

pytestmark = pytest.mark.unit


def test_group_raw_frames_preserves_members_and_uses_earliest_timestamp(
    tmp_path, log
) -> None:
    """Group equivalent frames without losing their ordered source paths."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.filterKeywordsExtras = ["mjd-obs", "vis temp c"]
    rawFrames = pd.DataFrame(
        {
            "eso seq arm": ["VIS", "VIS"],
            "eso dpr type": ["BIAS", "BIAS"],
            "file": ["bias-2.fits", "bias-1.fits"],
            "filepath": ["./raw/bias-2.fits", "./raw/bias-1.fits"],
            "date-obs": [
                "2024-01-02T03:05:06.900",
                "2024-01-02T03:04:05.100",
            ],
            "mjd-obs": [60311.12855, 60311.12784],
            "vis temp c": [10.0, 12.0],
        }
    )

    grouped = organiser._group_raw_frames(
        rawFrames,
        ["eso seq arm", "eso dpr type", "set_first_file"],
    )

    assert grouped.loc[0, "counts"] == 2
    assert grouped.loc[0, "date-obs"] == "20240102T030405"
    assert grouped.loc[0, "filepaths"] == [
        "./raw/bias-2.fits",
        "./raw/bias-1.fits",
    ]
    assert grouped.loc[0, "mjd-obs"] == pytest.approx(60311.128195)
    assert grouped.loc[0, "vis temp c"] == pytest.approx(11.0)


def test_group_raw_frames_can_omit_member_paths_and_start_timestamp(
    tmp_path, log
) -> None:
    """Allow callers to request aggregation metadata without member details."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.filterKeywordsExtras = ["mjd-obs"]
    rawFrames = pd.DataFrame(
        {
            "eso seq arm": ["VIS", "VIS"],
            "file": ["bias-1.fits", "bias-2.fits"],
            "filepath": ["one", "two"],
            "date-obs": ["2024-01-02T03:04:05", "2024-01-02T03:05:05"],
            "mjd-obs": [60311.1, 60311.2],
        }
    )

    grouped = organiser._group_raw_frames(
        rawFrames,
        ["eso seq arm"],
        addFilepaths=False,
        addStartDate=False,
    )

    assert grouped.to_dict(orient="records") == [
        {"eso seq arm": "VIS", "mjd-obs": pytest.approx(60311.15), "counts": 2}
    ]


def test_predict_product_frames_builds_stable_reduced_fits_path(tmp_path, log) -> None:
    """Translate a complete raw group into its predicted FITS product row."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    rawGroups = raw_group_table()
    originalGroups = rawGroups.copy(deep=True)

    result = organiser.predict_product_frames(
        [
            {
                "fits image": {
                    "eso pro type": "REDUCED",
                    "eso pro tech": "IMAGE",
                    "eso pro catg": "MASTER_BIAS",
                    "replace": [{"from": "MBIAS", "to": "MASTER_BIAS"}],
                }
            }
        ],
        rawGroups,
        "mbias",
    )

    product = pd.read_sql("SELECT * FROM product_frames", organiser.conn).iloc[0]
    assert result == 1
    assert product["file"] == ("20240102T030405_VIS_1X1_FAST_MASTER_BIAS_SOXS.fits")
    assert product["filepath"] == (
        "./reduced/2024-01-01/soxs-mbias/"
        "20240102T030405_VIS_1X1_FAST_MASTER_BIAS_SOXS.fits"
    )
    assert product["eso pro catg"] == "MASTER_BIAS_VIS"
    pd.testing.assert_frame_equal(rawGroups, originalGroups)


def test_predict_product_frames_reports_existing_incomplete_products(
    tmp_path, log
) -> None:
    """Report outstanding products when no new raw groups are available."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute(
        "CREATE TABLE product_frames (recipe TEXT, complete INTEGER)"
    )
    organiser.conn.executemany(
        "INSERT INTO product_frames VALUES (?, ?)",
        [("mflat", 0), ("mflat", 1), ("mbias", 0)],
    )

    incompleteCount = organiser.predict_product_frames([], pd.DataFrame(), "mflat")

    assert incompleteCount == 1


def test_predict_product_frames_binds_a_hostile_recipe_as_a_parameter(
    tmp_path, log
) -> None:
    """A hostile recipe name matches no real row in a real database, rather than counting every row.

    DY-254: `recipe` was interpolated directly into the `count(*)` query text.
    """
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    organiser.conn.execute(
        "CREATE TABLE product_frames (recipe TEXT, complete INTEGER)"
    )
    organiser.conn.executemany(
        "INSERT INTO product_frames VALUES (?, ?)",
        [("mflat", 0), ("mflat", 1), ("mbias", 0)],
    )
    # IF STILL INTERPOLATED, THE COMMENT MARKER `--` STRIPS THE TRAILING
    # `and complete< 1` CONDITION AND THE `OR '1'='1'` MATCHES EVERY ROW,
    # COUNTING ALL 3 ROWS INSTEAD OF THE 0 WHOSE `recipe` REALLY EQUALS THE
    # HOSTILE STRING.
    hostileRecipe = "mflat' OR '1'='1' -- "

    incompleteCount = organiser.predict_product_frames(
        [], pd.DataFrame(), hostileRecipe
    )

    assert incompleteCount == 0
