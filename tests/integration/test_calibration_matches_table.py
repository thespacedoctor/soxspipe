"""Integration tests for `refresh_calibration_matches` against a copy of the template `soxspipe.db`.

The equivalence test fills `product_frames` with a seeded mix of products and checks that each new `cal_<type>`
view, read after a refresh, holds the same `sof -> file` pairs as the pre-DY-1513 view SQL kept in
`data/legacy_calibration_views.sql`. Every SOF has its own `mjd-obs` and writes as many products as a real SOF of
its recipe, so the old views' unordered tie-breaks never decide a match.
"""

from __future__ import annotations

import random
import shutil
import sqlite3
from pathlib import Path

import pytest

import soxspipe
from soxspipe.commonutils.calibration_matching import CALIBRATION_RULES, refresh_calibration_matches

pytestmark = pytest.mark.integration

TEMPLATE_DB = Path(soxspipe.__file__).parent / "resources" / "soxspipe.db"
LEGACY_VIEWS = Path(__file__).parent / "data" / "legacy_calibration_views.sql"

RECIPE_ORDER = {
    "mbias": 0,
    "mdark": 1,
    "disp_solution": 2,
    "order_centres": 3,
    "mflat": 4,
    "spat_solution": 5,
    "stare_std": 6,
    "nod_std": 7,
    "offset_std": 8,
    "stare_obj": 9,
    "nod_obj": 10,
    "offset_obj": 11,
}

SEEDED_SOF_COUNT = 400

PRODUCT_TAGS = {
    "mflat": ["MASTER_FLAT", "ORDER_TAB"],
    "spat_solution": ["DISP_TAB", "DISP_IMAGE"],
    "stare_std": ["OBJECT_TAB", "RESP_TAB"],
    "nod_std": ["OBJECT_TAB", "RESP_TAB"],
    "offset_std": ["OBJECT_TAB", "RESP_TAB"],
}

PRODUCT_COLUMNS = [
    "file",
    "sof",
    "recipe",
    "recipe_order",
    "eso seq arm",
    "mjd-obs",
    "eso pro catg",
    "exptime",
    "binning",
    "gain",
    "rospeed",
    "slit",
    "status",
    "complete",
    "filepath",
]


@pytest.fixture
def conn(tmp_path):
    """Open an autocommit connection to a fresh copy of the template database."""
    dbPath = tmp_path / "soxspipe.db"
    shutil.copyfile(TEMPLATE_DB, dbPath)
    connection = sqlite3.connect(dbPath, autocommit=True)
    yield connection
    connection.close()


def insert_products(connection, rows):
    """Insert `product_frames` rows given as dicts keyed by `PRODUCT_COLUMNS`."""
    columns = ", ".join(f'"{column}"' for column in PRODUCT_COLUMNS)
    placeholders = ", ".join("?" for _ in PRODUCT_COLUMNS)
    connection.executemany(
        f"INSERT INTO product_frames ({columns}) VALUES ({placeholders})",  # noqa: S608
        [tuple(row.get(column) for column in PRODUCT_COLUMNS) for row in rows],
    )


def product_row(file, recipe, mjd, **overrides):
    """Return one VIS product row; `sof` defaults to one SOF per product."""
    row = {
        "file": file,
        "sof": f"{file}.sof",
        "recipe": recipe,
        "recipe_order": RECIPE_ORDER[recipe],
        "eso seq arm": "VIS",
        "mjd-obs": mjd,
        "eso pro catg": f"{recipe.upper()}_VIS",
        "exptime": 0.0,
        "binning": "1x1",
        "gain": 1.2,
        "rospeed": "1",
        "slit": "SLIT1.0",
        "status": "pass",
        "complete": 1,
        "filepath": f"./product/{file}",
    }
    row.update(overrides)
    return row


def pairs(connection, view):
    """Return the set of `(sof, file)` pairs a view holds."""
    return set(connection.execute(f"SELECT sof, file FROM {view}"))  # noqa: S608


def seeded_products(seed):
    """Return a seeded mix of upstream and downstream products covering every rule's hard matches."""
    rng = random.Random(seed)  # noqa: S311 -- SEEDED TEST DATA, NOT CRYPTOGRAPHY
    mjds = rng.sample(range(1_000_000), SEEDED_SOF_COUNT)
    rows = []

    def next_mjd():
        return 60000 + mjds.pop() / 10_000

    def mode(arm):
        if arm == "NIR":
            return {"eso seq arm": "NIR", "binning": None, "gain": None, "rospeed": "3"}
        return {
            "eso seq arm": arm,
            "binning": rng.choice(["1x1", "2x2"]),
            "gain": rng.choice([1.2, 2.4]),
            "rospeed": rng.choice(["1", "2"]),
        }

    for index in range(SEEDED_SOF_COUNT):
        recipe = rng.choice(list(RECIPE_ORDER))
        arm = rng.choice(["VIS", "NIR"])
        isDownstream = rng.random() < 0.4
        sof = f"sof-{index}.sof"
        base = {
            "recipe": recipe,
            "recipe_order": RECIPE_ORDER[recipe],
            "sof": sof,
            "exptime": rng.choice([0.0, 10.0, 300.0, 600.0]),
            "slit": rng.choice(["SLIT1.0", "SLIT1.5", "SLIT5.0"]),
            "status": None if isDownstream else rng.choice(["pass", "pass", None, "fail"]),
            "complete": rng.choice([0, -1]) if isDownstream else rng.choice([1, 1, 1, 0]),
            # EVERY PRODUCT OF ONE SOF SHARES ITS `mjd-obs`, AS IN A REAL WORKSPACE
            "mjd-obs": next_mjd(),
            **mode(arm),
        }
        # REAL SOFS OF THESE RECIPES WRITE TWO PRODUCTS; THE OLD VIEWS' k = 2 AND `ORDER%` FILTER RELY ON IT
        tags = PRODUCT_TAGS.get(recipe, [recipe.upper()])
        if isDownstream and len(tags) == 1 and rng.random() < 0.3:
            tags = [*tags, f"{tags[0]}_EXTRA"]
        for part, tag in enumerate(tags):
            file = f"p-{index}-{part}.fits"
            rows.append({**base, "file": file, "eso pro catg": f"{tag}_{arm}", "filepath": f"./product/{file}"})
    return rows


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_new_views_agree_with_the_legacy_view_sql(conn, seed):
    # ARRANGE
    insert_products(conn, seeded_products(seed))
    conn.executescript(LEGACY_VIEWS.read_text())

    # ACT
    refresh_calibration_matches(conn)

    # ASSERT
    for calType in CALIBRATION_RULES:
        legacy = pairs(conn, f"legacy_cal_{calType}")
        assert pairs(conn, f"cal_{calType}") == legacy, calType


def test_the_seeded_products_exercise_every_calibration_type(conn):
    # ARRANGE
    insert_products(conn, seeded_products(1))

    # ACT
    refresh_calibration_matches(conn)

    # ASSERT -- THE EQUIVALENCE TEST WOULD BE VACUOUS FOR A TYPE THAT NEVER MATCHES
    counts = dict(conn.execute("SELECT cal_type, COUNT(*) FROM calibration_matches GROUP BY cal_type"))
    assert set(counts) == set(CALIBRATION_RULES)


def test_a_scoped_refresh_replaces_only_that_recipe_and_arm(conn):
    # ARRANGE
    insert_products(
        conn,
        [
            product_row("bias-old", "mbias", 100.0),
            product_row("dsol-vis", "disp_solution", 105.0, status=None, complete=0),
            product_row("mflat-vis", "mflat", 105.0, status=None, complete=0),
        ],
    )
    refresh_calibration_matches(conn, ["mbias"])
    insert_products(conn, [product_row("bias-new", "mbias", 104.9)])

    # ACT
    refresh_calibration_matches(conn, ["mbias"], recipe="disp_solution", arm="VIS")

    # ASSERT
    assert pairs(conn, "cal_mbias") == {("dsol-vis.sof", "bias-new"), ("mflat-vis.sof", "bias-old")}


def test_a_refresh_drops_matches_for_sofs_that_are_now_complete(conn):
    # ARRANGE
    insert_products(
        conn,
        [product_row("bias", "mbias", 100.0), product_row("dsol", "disp_solution", 105.0, status=None, complete=0)],
    )
    refresh_calibration_matches(conn, ["mbias"])
    conn.execute("UPDATE product_frames SET complete = 1 WHERE recipe = 'disp_solution'")

    # ACT
    refresh_calibration_matches(conn, ["mbias"])

    # ASSERT
    assert pairs(conn, "cal_mbias") == set()


def test_the_views_expose_the_columns_the_sof_build_reads(conn):
    # ARRANGE
    insert_products(
        conn,
        [product_row("bias", "mbias", 100.0), product_row("dsol", "disp_solution", 105.0, status=None, complete=-1)],
    )

    # ACT
    refresh_calibration_matches(conn, ["mbias"])
    row = conn.execute("SELECT sof, file, upstream_tag, filepath, upstream_status, complete FROM cal_mbias").fetchall()

    # ASSERT
    assert row == [("dsol.sof", "bias", "MBIAS_VIS", "./product/bias", "pass", -1)]


def test_a_refresh_inside_a_legacy_transaction_leaves_it_to_the_caller(tmp_path):
    # ARRANGE
    dbPath = tmp_path / "soxspipe.db"
    shutil.copyfile(TEMPLATE_DB, dbPath)
    connection = sqlite3.connect(dbPath)
    insert_products(
        connection,
        [product_row("bias", "mbias", 100.0), product_row("dsol", "disp_solution", 105.0, status=None, complete=0)],
    )

    # ACT
    refresh_calibration_matches(connection, ["mbias"])
    connection.rollback()

    # ASSERT
    assert pairs(connection, "cal_mbias") == set()
    connection.close()


def test_an_unknown_calibration_type_is_rejected(conn):
    # ACT / ASSERT
    with pytest.raises(ValueError, match="bogus"):
        refresh_calibration_matches(conn, ["bogus"])


def test_a_single_calibration_type_may_be_named_as_a_string(conn):
    # ARRANGE
    insert_products(
        conn,
        [product_row("bias", "mbias", 100.0), product_row("dsol", "disp_solution", 105.0, status=None, complete=0)],
    )

    # ACT
    refresh_calibration_matches(conn, "mbias")

    # ASSERT
    assert pairs(conn, "cal_mbias") == {("dsol.sof", "bias")}
