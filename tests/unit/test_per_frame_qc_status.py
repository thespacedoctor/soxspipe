"""Per-frame stare QCs are reported under `<sof>_F<N>.sof` and still decide the SOF's status (DY-1611)."""

from __future__ import annotations

import sqlite3

import pytest

from soxspipe.commonutils.data_organiser import data_organiser

pytestmark = pytest.mark.unit

PARENT = "20240102T030500_VIS_STARE_OBJ.sof"
OTHER = "20240102T030500_VIS_STARE_OBJ_F1_EXTRA.sof"
SETTINGS = {"soxs-stare": {"qc-acceptable-ranges": {"snr-median": [5.0, 1000.0]}}}


def _status_after_range_check(qcRows, startStatus):
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
        "insert into quality_control values ('soxs-stare', 'SNR MEDIAN', '-1', null, null, 'pass', ?, ?)",
        qcRows,
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
    statuses = _status_after_range_check(qcRows, startStatus="fail")

    # ASSERT
    assert statuses[PARENT] == "pass"
    assert statuses[OTHER] == "fail"
