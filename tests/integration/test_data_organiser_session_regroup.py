"""Integration contracts for a session regrouping every raw frame in the workspace (DY-274)."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import OutdatedWorkspaceError
from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
from tests.factories import raw_flat_run_table, raw_frame_table, workspace_organiser
from tests.integration.test_data_organiser_late_frames import (
    _calibration_chain_frames,
    _insert_frames,
    _pinhole_frame,
)

pytestmark = pytest.mark.integration

INSTRUMENT = "SOXS"
MBIAS = "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"
MFLAT = "20240102T060411_VIS_1X1_FAST_MFLAT_1_0_10_0S_SOXS.sof"
FIRST_HALF = [1, 2, 3, 4, 5]
SECOND_HALF = [6, 7, 8, 9, 10]


@pytest.fixture
def prepared_workspace(tmp_path, log, monkeypatch):
    """Build a workspace holding ``frames``, create session `alpha` and run a plain `prep` in it."""
    organisers = []

    def build(frames):
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        organiser = workspace_organiser(tmp_path, log=log)
        connection, _ = organiser._get_or_create_db_connection()
        organiser.conn = connection
        _insert_frames(connection, frames)
        organiser.session_create("alpha")

        monkeypatch.setattr(organiser, "_fits_files_exist", lambda: True)
        monkeypatch.setattr(organiser, "_move_misc_files", lambda: None)
        monkeypatch.setattr(organiser, "_sync_raw_frames", lambda: organiser._select_instrument(inst=INSTRUMENT))
        monkeypatch.setattr(organiser, "get_incomplete_sets_report", lambda: (None, None))
        organiser.prepare(report=False)
        organisers.append(organiser)
        return organiser

    yield build
    for organiser in organisers:
        if getattr(organiser, "conn", None) is not None:
            organiser.conn.close()


def _sof_map_rows(organiser, sessionId):
    """Return every row of one session's SOF map in a stable order."""
    # TABLE NAME IS A LITERAL PREFIX PLUS A SESSION ID THAT THE TEST CONTROLS
    sqlQuery = f"SELECT file, tag, sof, filepath, complete FROM sof_map_{sessionId} ORDER BY sof, file"  # noqa: S608
    return organiser.conn.execute(sqlQuery).fetchall()


def _sof_files(organiser, sessionId):
    """Return the name and text of every `.sof` file of one session."""
    sofDir = Path(organiser.sessionsDir) / sessionId / "sof"
    return {path.name: path.read_text() for path in sorted(sofDir.glob("*.sof"))}


# A NEW SESSION GROUPS FRAMES THAT ANOTHER SESSION ALREADY GROUPED


def test_session_create_gives_the_new_session_the_sofs_another_session_grouped(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    alphaRows = _sof_map_rows(organiser, "alpha")
    alphaFiles = _sof_files(organiser, "alpha")
    assert {row[2] for row in alphaRows} == {MBIAS}

    # ACT
    organiser.session_create("beta")

    # ASSERT
    assert _sof_map_rows(organiser, "beta") == alphaRows
    assert _sof_files(organiser, "beta") == alphaFiles


def _second_bias_pair():
    """Two more biases, taken a night later, that form a second master-bias set."""
    frames = raw_frame_table()
    frames["file"] = ["bias-3.fits", "bias-4.fits"]
    frames["filepath"] = ["./raw/2024-01-01/bias-3.fits", "./raw/2024-01-01/bias-4.fits"]
    frames["date-obs"] = ["2024-01-03T03:04:05.000", "2024-01-03T03:04:06.000"]
    frames["mjd-obs"] = frames["mjd-obs"] + 1
    frames["eso tpl start"] = "2024-01-03T03:04:00"
    frames["night start date"] = "2024-01-02"
    frames["night start mjd"] = 60311
    frames["set_first_file"] = "bias-3.fits"
    return frames


def test_session_switch_gives_the_session_frames_that_were_grouped_while_another_session_was_active(
    prepared_workspace,
) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.session_create("beta")
    organiser.session_switch("alpha")
    _insert_frames(organiser.conn, _second_bias_pair())
    organiser.prepare(report=False)
    laterMbias = "20240103T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"
    assert {row[2] for row in _sof_map_rows(organiser, "alpha")} == {MBIAS, laterMbias}
    assert {row[2] for row in _sof_map_rows(organiser, "beta")} == {MBIAS}

    # ACT
    organiser.session_switch("beta")

    # ASSERT
    assert _sof_map_rows(organiser, "beta") == _sof_map_rows(organiser, "alpha")
    assert set(_sof_files(organiser, "beta")) == {MBIAS, laterMbias}


def _seed_status(organiser, sessionId, sof, status):
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE
    organiser.conn.execute(
        f"UPDATE product_frames SET {statusColumn} = ? WHERE sof = ?",  # noqa: S608
        (status, sof),
    )


def _statuses(organiser, sof):
    """Return the distinct (status_alpha, status_beta, complete) tuples of one SOF."""
    return set(
        organiser.conn.execute(
            "SELECT DISTINCT status_alpha, status_beta, complete FROM product_frames WHERE sof = ?", (sof,)
        ).fetchall()
    )


# A REGROUP KEEPS WHAT A SET ALREADY HAS


def test_prep_twice_keeps_every_sof_name_and_status_in_every_session(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.session_create("beta")
    _seed_status(organiser, "beta", MBIAS, "fail")
    organiser.session_switch("alpha")
    _seed_status(organiser, "alpha", MBIAS, "pass")
    mapBefore = _sof_map_rows(organiser, "alpha")
    filesBefore = _sof_files(organiser, "alpha")

    # ACT
    organiser.prepare(report=False)
    organiser.prepare(report=False)
    organiser.session_switch("beta")
    organiser.session_switch("alpha")

    # ASSERT
    assert _statuses(organiser, MBIAS) == {("pass", "fail", 1)}
    assert _sof_map_rows(organiser, "alpha") == mapBefore
    assert _sof_files(organiser, "alpha") == filesBefore


def test_regrouping_two_sessions_over_the_same_frames_adds_no_product_rows(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    productsBefore = organiser.conn.execute("SELECT file, recipe, sof FROM product_frames ORDER BY file").fetchall()

    # ACT
    organiser.session_create("beta")
    organiser.prepare(report=False)
    organiser.session_switch("alpha")

    # ASSERT
    productsAfter = organiser.conn.execute("SELECT file, recipe, sof FROM product_frames ORDER BY file").fetchall()
    duplicates = organiser.conn.execute(
        "SELECT file, recipe FROM product_frames GROUP BY file, recipe HAVING COUNT(*) > 1"
    ).fetchall()
    assert productsBefore
    assert productsAfter == productsBefore
    assert duplicates == []


# A LATE FRAME THAT REJOINS A SET IN ONE SESSION REACHES THE OTHER SESSIONS WHEN THEY REGROUP


def _half_flat_workspace(prepared_workspace):
    """Return a workspace whose active session `alpha` has grouped half of a flat run, plus an idle session `beta`."""
    frames = pd.concat([_calibration_chain_frames(), raw_flat_run_table(FIRST_HALF)])
    organiser = prepared_workspace(frames)
    organiser.session_create("beta")
    organiser.session_switch("alpha")
    _seed_status(organiser, "alpha", MFLAT, "pass")
    _seed_status(organiser, "beta", MFLAT, "pass")
    return organiser


def _late_prepare(organiser):
    """Add the second half of the flat run to the workspace, then run a plain `prep` in the active session."""
    _insert_frames(organiser.conn, raw_flat_run_table(SECOND_HALF))
    organiser.prepare(report=False)


def _mflat_members(organiser, sessionId):
    # TABLE NAME IS A LITERAL PREFIX PLUS A SESSION ID THAT THE TEST CONTROLS
    sqlQuery = f"SELECT file FROM sof_map_{sessionId} WHERE sof = ? AND file LIKE 'flat-%'"  # noqa: S608
    return {row[0] for row in organiser.conn.execute(sqlQuery, (MFLAT,))}


def test_a_late_frame_rejoined_in_one_session_leaves_the_other_sessions_untouched(prepared_workspace) -> None:
    # ARRANGE
    organiser = _half_flat_workspace(prepared_workspace)
    betaMap = _sof_map_rows(organiser, "beta")
    betaFiles = _sof_files(organiser, "beta")

    # ACT
    _late_prepare(organiser)

    # ASSERT
    assert _mflat_members(organiser, "alpha") == {f"flat-{n}.fits" for n in FIRST_HALF + SECOND_HALF}
    assert _sof_map_rows(organiser, "beta") == betaMap
    assert _sof_files(organiser, "beta") == betaFiles
    assert {row[1] for row in _statuses(organiser, MFLAT)} == {"pass"}


def test_switching_to_a_session_gives_it_the_late_frame_and_queues_the_set_again(prepared_workspace) -> None:
    # ARRANGE
    organiser = _half_flat_workspace(prepared_workspace)
    _late_prepare(organiser)
    _seed_status(organiser, "alpha", MFLAT, "pass")
    alphaMap = _sof_map_rows(organiser, "alpha")

    # ACT
    organiser.session_switch("beta")

    # ASSERT
    wholeRun = {f"flat-{n}.fits" for n in FIRST_HALF + SECOND_HALF}
    sofText = _sof_files(organiser, "beta")[MFLAT]
    assert _mflat_members(organiser, "beta") == wholeRun
    assert all(f"flat-{n}.fits" in sofText for n in FIRST_HALF + SECOND_HALF)
    assert _statuses(organiser, MFLAT) == {("pass", None, 1)}
    assert _sof_map_rows(organiser, "alpha") == alphaMap


# TWO TEMPLATE RUNS THAT SHARE EVERY GROUPING KEY STAY TWO SOFS IN EVERY SESSION


def test_session_create_keeps_two_template_runs_with_the_same_keys_as_two_sofs(prepared_workspace) -> None:
    # ARRANGE: THE SECOND RUN ARRIVES AFTER THE FIRST WAS GROUPED, SO IT FORMS ITS OWN SOF
    frames = pd.concat([_calibration_chain_frames(), raw_flat_run_table(FIRST_HALF + SECOND_HALF)])
    organiser = prepared_workspace(frames)
    _insert_frames(
        organiser.conn, raw_flat_run_table(FIRST_HALF + SECOND_HALF, namePrefix="later", tplStart="2024-01-02T09:04:00")
    )
    organiser.prepare(report=False)
    alphaRows = _sof_map_rows(organiser, "alpha")
    assert len({row[2] for row in alphaRows if "MFLAT" in row[2]}) == 2

    # ACT
    organiser.session_create("beta")

    # ASSERT
    assert _sof_map_rows(organiser, "beta") == alphaRows


# A SET THAT GAINED AN EARLIER FRAME KEEPS ITS NAME AFTER THE IN-MEMORY PIN IS GONE, AND IN NEW SESSIONS

MFLAT_NAMED_AFTER_EXPNO_6 = "20240102T060416_VIS_1X1_FAST_MFLAT_1_0_10_0S_SOXS.sof"


def _workspace_with_early_late_frames(prepared_workspace):
    """Return a workspace whose active session `alpha` grouped the second half of a run, then gained the first half."""
    frames = pd.concat([_calibration_chain_frames(), raw_flat_run_table(SECOND_HALF)])
    organiser = prepared_workspace(frames)
    _insert_frames(organiser.conn, raw_flat_run_table(FIRST_HALF))
    organiser.prepare(report=False)
    return organiser


def _flat_sofs(organiser, sessionId):
    return {row[2] for row in _sof_map_rows(organiser, sessionId) if "MFLAT" in row[2]}


def test_prep_in_a_new_process_keeps_the_name_of_a_set_that_gained_an_earlier_frame(prepared_workspace) -> None:
    # ARRANGE
    organiser = _workspace_with_early_late_frames(prepared_workspace)
    assert _flat_sofs(organiser, "alpha") == {MFLAT_NAMED_AFTER_EXPNO_6}
    organiser.pinnedSofNames = {}

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    productSofs = organiser.conn.execute("SELECT DISTINCT sof FROM product_frames WHERE recipe = 'mflat'").fetchall()
    assert _flat_sofs(organiser, "alpha") == {MFLAT_NAMED_AFTER_EXPNO_6}
    assert productSofs == [(MFLAT_NAMED_AFTER_EXPNO_6,)]
    assert set(_sof_files(organiser, "alpha")) == {row[2] for row in _sof_map_rows(organiser, "alpha")}


def test_a_session_created_later_names_a_set_that_gained_an_earlier_frame_as_the_other_session_does(
    prepared_workspace,
) -> None:
    # ARRANGE
    organiser = _workspace_with_early_late_frames(prepared_workspace)
    organiser.pinnedSofNames = {}

    # ACT
    organiser.session_create("beta")

    # ASSERT
    assert _flat_sofs(organiser, "beta") == {MFLAT_NAMED_AFTER_EXPNO_6}
    assert _sof_map_rows(organiser, "beta") == _sof_map_rows(organiser, "alpha")


# A FAILED SET STAYS FAILED THROUGH A REGROUP


def test_a_failed_set_keeps_its_status_and_its_map_rows_through_a_regroup(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.session_create("beta")
    organiser.session_switch("alpha")
    _seed_status(organiser, "alpha", MBIAS, "fail")
    mapBefore = _sof_map_rows(organiser, "alpha")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert _statuses(organiser, MBIAS) == {("fail", None, 1)}
    assert _sof_map_rows(organiser, "alpha") == mapBefore


# TWO KNOWN PINHOLE SETS THAT SHARE EVERY KEY STAY TWO SOFS


def test_two_pinhole_arcs_that_share_every_key_stay_two_sofs_in_every_session(prepared_workspace) -> None:
    # ARRANGE: THE SECOND ARC ARRIVES AFTER THE FIRST WAS GROUPED AND HAS THE SAME KEYS, TEMPLATE START AND `mjd-obs`
    organiser = prepared_workspace(_calibration_chain_frames())
    secondArc = _pinhole_frame(raw_frame_table(), "arc-2", "LAMP,WAVE", "2024-01-02T04:04:00", "2024-01-02T04:04:09")
    _insert_frames(organiser.conn, secondArc)
    organiser.prepare(report=False)
    alphaArcSofs = {row[2] for row in _sof_map_rows(organiser, "alpha") if "DISP_SOLUTION" in row[2]}
    assert len(alphaArcSofs) == 2

    # ACT
    organiser.prepare(report=False)
    organiser.session_create("beta")

    # ASSERT
    assert {row[2] for row in _sof_map_rows(organiser, "alpha") if "DISP_SOLUTION" in row[2]} == alphaArcSofs
    assert {row[2] for row in _sof_map_rows(organiser, "beta") if "DISP_SOLUTION" in row[2]} == alphaArcSofs


# AN OUT-OF-DATE WORKSPACE IS REFUSED BEFORE ANYTHING CHANGES


def _workspace_state(organiser):
    """Return what `session_create` and `session_switch` change: session directories, ID file, tables, view, links."""
    sessionsDir = Path(organiser.sessionsDir)
    rootDir = Path(organiser.rootDir)
    return {
        "sessions": sorted(path.name for path in sessionsDir.iterdir()),
        "sessionId": Path(organiser.sessionIdFile).read_text(encoding="utf-8"),
        "tables": sorted(row[0] for row in organiser.conn.execute("SELECT name FROM sqlite_master")),
        "columns": sorted(row[1] for row in organiser.conn.execute("PRAGMA table_info(product_frames)")),
        "view": organiser.conn.execute("SELECT sql FROM sqlite_master WHERE name = 'sof_map'").fetchone(),
        "links": {name: Path(rootDir / name).resolve() for name in ("sof", "qc", "reduced", "soxspipe.yaml")},
        "attributes": (organiser.sessionId, organiser.sessionPath),
    }


def test_session_create_refuses_an_outdated_workspace_before_it_changes_anything(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.conn.execute("DROP TABLE calibration_matches")
    stateBefore = _workspace_state(organiser)

    # ACT
    with pytest.raises(OutdatedWorkspaceError):
        organiser.session_create("beta")

    # ASSERT
    assert _workspace_state(organiser) == stateBefore


def test_session_switch_refuses_an_outdated_workspace_before_it_changes_anything(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.session_create("beta")
    organiser.session_switch("alpha")
    organiser.conn.execute("DROP TABLE calibration_matches")
    stateBefore = _workspace_state(organiser)

    # ACT
    with pytest.raises(OutdatedWorkspaceError):
        organiser.session_switch("beta")

    # ASSERT
    assert _workspace_state(organiser) == stateBefore


# A WORKSPACE WITH NO RAW FRAME CAN STILL GET A SESSION


def test_session_create_works_in_a_workspace_with_no_raw_frames(tmp_path, log) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    connection, _ = organiser._get_or_create_db_connection()
    organiser.conn = connection

    # ACT
    sessionId = organiser.session_create("alpha")

    # ASSERT
    assert sessionId == "alpha"
    assert (Path(organiser.sessionsDir) / "alpha" / "soxspipe.yaml").is_file()
    assert Path(organiser.sessionIdFile).read_text(encoding="utf-8") == "alpha"


# A SYMLINK IN `sof/` IS NEVER FOLLOWED WHEN A STALE SOF FILE IS DELETED


def test_a_stale_sof_symlink_is_skipped_and_its_target_survives(prepared_workspace) -> None:
    # ARRANGE
    frames = pd.concat([_calibration_chain_frames(), raw_flat_run_table(FIRST_HALF)])
    organiser = prepared_workspace(frames)
    sofDir = Path(organiser.sessionsDir) / "alpha" / "sof"
    target = sofDir / "keep.txt"
    target.write_text("keep")
    (sofDir / MFLAT).unlink()
    (sofDir / MFLAT).symlink_to(target)

    # ACT
    _late_prepare(organiser)

    # ASSERT
    warnings = [message for level, message in organiser.log.messages if level == "warning"]
    assert target.read_text() == "keep"
    assert (sofDir / MFLAT).is_symlink()
    assert any(MFLAT in message and "symlink" in message for message in warnings)


# THE FIRST STEPS OF A REGROUP CHANGE THE DATABASE IN ONE TRANSACTION


class _FailsOnRequeue:
    """Connection proxy that behaves like the real one except that the product requeue statement fails."""

    def __init__(self, connection):
        self._connection = connection

    def __getattr__(self, name):
        return getattr(self._connection, name)

    def execute(self, sqlQuery, *args):
        if "SET complete = 0 WHERE complete != 0" in sqlQuery:
            raise sqlite3.OperationalError("simulated failure")
        return self._connection.execute(sqlQuery, *args)


def test_a_failed_regroup_start_leaves_processed_and_product_flags_as_they_were(prepared_workspace) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiser.conn.execute("UPDATE raw_frames SET processed = 5")
    realConnection = organiser.conn
    organiser.conn = _FailsOnRequeue(realConnection)

    # ACT
    with pytest.raises(sqlite3.OperationalError):
        organiser.build_sof_files(rejoinLateFrames=True)

    # ASSERT
    organiser.conn = realConnection
    assert organiser.conn.execute("SELECT DISTINCT processed FROM raw_frames").fetchall() == [(5,)]
    assert not organiser.conn.in_transaction


# GROUPING WORK IS NOT REPEATED


def test_the_known_set_labels_are_built_once_per_recipe_in_a_regroup(prepared_workspace, monkeypatch) -> None:
    # ARRANGE
    organiser = prepared_workspace(raw_frame_table())
    organiserModule = sys.modules["soxspipe.commonutils.data_organiser"]
    realLabels = organiserModule.partition_labels
    recipes = []

    def counting_labels(pins, recipe, excludedSofs=()):
        recipes.append(recipe)
        return realLabels(pins, recipe, excludedSofs)

    monkeypatch.setattr(organiserModule, "partition_labels", counting_labels)

    # ACT
    organiser.build_sof_files(rejoinLateFrames=True)

    # ASSERT
    assert recipes
    assert len(recipes) == len(set(recipes))
