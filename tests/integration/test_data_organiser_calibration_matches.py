"""Integration contracts for how the data organiser keeps `calibration_matches` current (DY-1513)."""

from __future__ import annotations

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import OutdatedWorkspaceError
from tests.factories import raw_flat_run_table, raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

INSTRUMENT = "SOXS"
MFLAT = "20240102T060411_VIS_1X1_FAST_MFLAT_1_0_10_0S_SOXS.sof"


def _insert_frames(connection, frames):
    """Insert dataframe rows into `raw_frames`, keeping only the columns the table has."""
    columns = {row[1] for row in connection.execute("PRAGMA table_info(raw_frames);")}
    frames.loc[:, frames.columns.isin(columns)].to_sql("raw_frames", con=connection, index=False, if_exists="append")


def _pinhole_frame(biases, name, dprType, tplStart, obsTime):
    """One VIS pinhole exposure that is its own template run."""
    frame = biases.iloc[[0]].copy()
    frame["file"] = f"{name}.fits"
    frame["filepath"] = f"./raw/2024-01-01/{name}.fits"
    frame["eso dpr type"] = dprType
    frame["eso dpr tech"] = "ECHELLE,PINHOLE"
    frame["eso tpl nexp"] = 1
    frame["eso tpl name"] = "SOXS_orderdef" if dprType == "LAMP,FLAT" else "SOXS_arc"
    frame["exptime"] = 10.0
    frame["date-obs"] = f"{obsTime}.000"
    frame["mjd-obs"] = frame["mjd-obs"] + 0.04
    frame["eso tpl start"] = tplStart
    frame["set_first_file"] = f"{name}.fits"
    return frame


def _calibration_chain_frames():
    """Biases, a pinhole arc, an order-centre flat and a full slit-flat run."""
    biases = raw_frame_table()
    arc = _pinhole_frame(biases, "arc-1", "LAMP,WAVE", "2024-01-02T04:04:00", "2024-01-02T04:04:05")
    orderCentreFlat = _pinhole_frame(biases, "ocflat-1", "LAMP,FLAT", "2024-01-02T05:04:00", "2024-01-02T05:04:05")
    return pd.concat([biases, arc, orderCentreFlat, raw_flat_run_table(range(1, 11))])


def _workspace(tmp_path, log, monkeypatch, frames):
    """Return an organiser on a fresh workspace holding ``frames``, ready to `prepare`."""
    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection
    _insert_frames(connection, frames)
    organiser.session_create("science")

    monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
    monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
    monkeypatch.setattr(organiser, "_sync_raw_frames", lambda: organiser._select_instrument(inst=INSTRUMENT))
    return organiser


@pytest.fixture
def new_workspace(tmp_path, log, monkeypatch):
    """A fresh workspace holding the calibration chain."""
    organiser = _workspace(tmp_path, log, monkeypatch, _calibration_chain_frames())
    yield organiser
    organiser.conn.close()


@pytest.fixture
def flats_only_workspace(tmp_path, log, monkeypatch):
    """A fresh workspace holding a slit-flat run but none of the calibrations it needs."""
    organiser = _workspace(tmp_path, log, monkeypatch, raw_flat_run_table(range(1, 11)))
    yield organiser
    organiser.conn.close()


def _sof_members(organiser, sof):
    return {row[0] for row in organiser.conn.execute("SELECT file FROM sof_map WHERE sof = ?", (sof,))}


def test_prepare_still_builds_the_mflat_sof_from_its_nearest_calibrations(new_workspace):
    # ACT
    new_workspace.prepare(report=False)

    # ASSERT
    members = _sof_members(new_workspace, MFLAT)
    assert {f"flat-{expno}.fits" for expno in range(1, 11)} <= members
    tags = {row[0] for row in new_workspace.conn.execute("SELECT tag FROM sof_map WHERE sof = ?", (MFLAT,))}
    assert {"MASTER_BIAS_VIS", "ORDER_TAB_VIS"} <= tags


def test_prepare_leaves_no_calibration_match_for_a_sof_it_completed(new_workspace):
    # ACT
    new_workspace.prepare(report=False)

    # ASSERT
    completed = {row[0] for row in new_workspace.conn.execute("SELECT sof FROM product_frames WHERE complete = 1")}
    matched = {row[0] for row in new_workspace.conn.execute("SELECT DISTINCT sof FROM calibration_matches")}
    assert MFLAT in completed
    assert not completed & matched


def test_prepare_refuses_a_workspace_built_from_an_older_template(new_workspace):
    # ARRANGE
    new_workspace.conn.execute("DROP TABLE calibration_matches")

    # ACT / ASSERT
    with pytest.raises(OutdatedWorkspaceError, match="soxspipe prep --refresh"):
        new_workspace.prepare(report=False)


def test_a_session_sof_map_is_indexed_on_sof_and_filepath(new_workspace):
    # ACT
    indexedColumns = set()
    for index in new_workspace.conn.execute("PRAGMA index_list(sof_map_science)").fetchall():
        indexedColumns |= {row[2] for row in new_workspace.conn.execute(f"PRAGMA index_info('{index[1]}')")}

    # ASSERT
    assert {"sof", "filepath"} <= indexedColumns


def test_a_sof_left_mid_build_by_an_interrupted_prep_is_not_promoted_to_complete(flats_only_workspace):
    # ARRANGE -- A PREP KILLED BETWEEN MARKING A SOF -1 AND PROMOTING IT LEAVES THE -1 BEHIND
    flats_only_workspace.prepare(report=False)
    flats_only_workspace.conn.execute("UPDATE product_frames SET complete = -1 WHERE sof = ?", (MFLAT,))

    # ACT
    flats_only_workspace.prepare(report=False)

    # ASSERT -- NO MASTER BIAS OR ORDER CENTRES EXIST, SO THE FLAT SOF CANNOT BE COMPLETE
    completes = {
        row[0]
        for row in flats_only_workspace.conn.execute("SELECT complete FROM product_frames WHERE sof = ?", (MFLAT,))
    }
    assert completes == {0}
