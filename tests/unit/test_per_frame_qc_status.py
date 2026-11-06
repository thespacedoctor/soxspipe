"""Per-frame stare QCs are reported under `<sof>_F<N>.sof` and still decide the SOF's status (DY-1611)."""

from __future__ import annotations

import sqlite3

import pytest

from soxspipe.commonutils.data_organiser import data_organiser
from soxspipe.recipes.soxs_nod import CYCLE_INSPECTION_SUFFIX

pytestmark = pytest.mark.unit

PARENT = "20240102T030500_VIS_STARE_OBJ.sof"
OTHER = "20240102T030500_VIS_STARE_OBJ_F1_EXTRA.sof"
SETTINGS = {"soxs-stare": {"qc-acceptable-ranges": {"snr-median": [5.0, 1000.0]}}}


def _status_after_range_check(qcRows, startStatus, storedFlag="pass"):
    """Apply the QC-range queries to one parent SOF and return each product's session status."""
    organiser = data_organiser.__new__(data_organiser)
    organiser.sessionId = "base"
    organiser.settings = SETTINGS
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "create table quality_control (soxspipe_recipe text, qc_name text, qc_order text, qc_value_min real, "
        "qc_value_max real, qc_flag text, qc_value text, sof_name text)"
    )
    connection.executemany(
        "insert into quality_control values ('soxs-stare', 'SNR MEDIAN', '-1', null, null, ?, ?, ?)",
        [(storedFlag, *row) for row in qcRows],
    )
    connection.execute("create table product_frames (sof text, status_base text)")
    connection.executemany("insert into product_frames values (?, ?)", [(PARENT, startStatus), (OTHER, startStatus)])
    for sqlQuery, sqlParams in organiser._qc_acceptable_range_queries():
        connection.execute(sqlQuery, sqlParams)
    statuses = dict(connection.execute("select sof, status_base from product_frames").fetchall())
    connection.close()
    return statuses


def test_a_failed_frame_qc_fails_the_parent_sof() -> None:
    # ARRANGE
    qcRows = [("20.0", "20240102T030500_VIS_STARE_OBJ_F1.sof"), ("2.0", "20240102T030500_VIS_STARE_OBJ_F2.sof")]

    # ACT
    statuses = _status_after_range_check(qcRows, startStatus="pass")

    # ASSERT
    assert statuses[PARENT] == "fail"
    assert statuses[OTHER] == "pass"


def test_a_parent_sof_whose_frames_now_all_pass_is_passed_again() -> None:
    # ARRANGE
    qcRows = [("20.0", "20240102T030500_VIS_STARE_OBJ_F1.sof"), ("30.0", "20240102T030500_VIS_STARE_OBJ_F2.sof")]

    # ACT
    statuses = _status_after_range_check(qcRows, startStatus="fail", storedFlag="fail")

    # ASSERT
    assert statuses[PARENT] == "pass"
    assert statuses[OTHER] == "fail"


def test_a_failed_parent_sof_with_only_passing_frame_qcs_stays_failed() -> None:
    """A `fail` that no failing frame QC explains is an unfinished or crashed run, not a QC-range failure (DY-1698)."""
    # ARRANGE
    qcRows = [("20.0", "20240102T030500_VIS_STARE_OBJ_F1.sof"), ("30.0", "20240102T030500_VIS_STARE_OBJ_F2.sof")]

    # ACT
    statuses = _status_after_range_check(qcRows, startStatus="fail", storedFlag="pass")

    # ASSERT
    assert statuses[PARENT] == "fail"


CYCLES = "20240102T030500_VIS_NOD_STD_CYCLES.sof"
CYCLES_PARENT = "20240102T030500_VIS_NOD_STD.sof"


def _status_of_cycles_parent(qcRows, startStatus, extraSofs=(), storedFlag="pass"):
    """Apply the QC-range queries and return the session status of the nod standard SOF and any extra SOFs."""
    organiser = data_organiser.__new__(data_organiser)
    organiser.sessionId = "base"
    organiser.settings = {"soxs-nod-std": {"qc-acceptable-ranges": {"snr-median": [5.0, 1000.0]}}}
    connection = sqlite3.connect(":memory:")
    connection.execute(
        "create table quality_control (soxspipe_recipe text, qc_name text, qc_order text, qc_value_min real, "
        "qc_value_max real, qc_flag text, qc_value text, sof_name text)"
    )
    connection.executemany(
        "insert into quality_control values ('soxs-nod-std', 'SNR MEDIAN', '-1', null, null, ?, ?, ?)",
        [(storedFlag, *row) for row in qcRows],
    )
    connection.execute("create table product_frames (sof text, status_base text)")
    sofs = [CYCLES_PARENT, *extraSofs]
    connection.executemany("insert into product_frames values (?, ?)", [(sof, startStatus) for sof in sofs])
    for sqlQuery, sqlParams in organiser._qc_acceptable_range_queries():
        connection.execute(sqlQuery, sqlParams)
    statuses = dict(connection.execute("select sof, status_base from product_frames").fetchall())
    connection.close()
    return statuses


def test_a_failed_cycles_qc_fails_the_parent_sof() -> None:
    """The inspection pass of a `--no-stack` standard reports under `<sof>_CYCLES.sof` and still fails its SOF."""
    # ARRANGE
    qcRows = [("2.0", CYCLES)]

    # ACT
    statuses = _status_of_cycles_parent(qcRows, startStatus="pass")

    # ASSERT
    assert statuses[CYCLES_PARENT] == "fail"


def test_a_passing_cycles_qc_passes_a_failed_parent_sof_again() -> None:
    # ARRANGE
    qcRows = [("30.0", CYCLES)]

    # ACT
    statuses = _status_of_cycles_parent(qcRows, startStatus="fail", storedFlag="fail")

    # ASSERT
    assert statuses[CYCLES_PARENT] == "pass"


@pytest.mark.parametrize(
    "lookalike",
    [
        "20240102T030500_VIS_NOD_STD_CYCLESX.sof",
        "20240102T030500_VIS_NOD_STD_cycles.sof",
        "20240102T030500_VIS_NOD_STD_CYCLE.sof",
        "20240102T030500_VIS_NOD_STD_CYCLES1.sof",
    ],
)
def test_a_name_that_only_resembles_the_cycles_suffix_does_not_fail_a_sof(lookalike: str) -> None:
    # ARRANGE
    parent = lookalike.replace(".sof", "").rsplit("_", 1)[0] + ".sof"
    qcRows = [("2.0", lookalike)]

    # ACT
    statuses = _status_of_cycles_parent(qcRows, startStatus="pass", extraSofs=[parent])

    # ASSERT
    assert statuses[parent] == "pass"
    assert statuses[CYCLES_PARENT] == "pass"


def test_the_parent_rule_maps_a_name_built_from_the_recipe_suffix() -> None:
    """The SQL repeats the recipe's suffix as a literal, so a name built from the constant must still map."""
    # ARRANGE
    cyclesName = CYCLES_PARENT.removesuffix(".sof") + CYCLE_INSPECTION_SUFFIX + ".sof"

    # ACT
    statuses = _status_of_cycles_parent([("2.0", cyclesName)], startStatus="pass")

    # ASSERT
    assert statuses[CYCLES_PARENT] == "fail"
