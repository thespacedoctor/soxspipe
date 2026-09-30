"""Integration contracts for a raw frame that arrives after its set was grouped by `prep` (DY-263)."""

from __future__ import annotations

import sqlite3
import sys
from pathlib import Path

import pandas as pd
import pytest

from soxspipe.commonutils.data_organiser import _validate_owned_path
from soxspipe.commonutils.late_frame_rejoin import find_rejoined_sofs, reset_rejoined_sofs
from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
from tests.factories import raw_flat_run_table, raw_frame_table, workspace_organiser

pytestmark = pytest.mark.integration

INSTRUMENT = "SOXS"

FIRST_HALF = [1, 2, 3, 4, 5]
SECOND_HALF = [6, 7, 8, 9, 10]
WHOLE_RUN = FIRST_HALF + SECOND_HALF

# THE SOF NAME IS DERIVED FROM THE EARLIEST FRAME OF THE SET THAT WAS GROUPED
MFLAT_FROM_EXPNO_1 = "20240102T060411_VIS_1X1_FAST_MFLAT_1_0_10_0S_SOXS.sof"
MBIAS = "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof"
DISP_SOLUTION = "20240102T040405_VIS_1X1_FAST_DISP_SOLUTION_10_0S_SOXS.sof"
MFLAT_FROM_EXPNO_6 = "20240102T060416_VIS_1X1_FAST_MFLAT_1_0_10_0S_SOXS.sof"

PRODUCT_TYPES_OF_MFLAT = 2


def _raw_frames_columns(connection):
    return {row[1] for row in connection.execute("PRAGMA table_info(raw_frames);")}


def _insert_frames(connection, frames):
    """Insert dataframe rows into `raw_frames`, keeping only the columns the table has."""
    known = frames.loc[:, frames.columns.isin(_raw_frames_columns(connection))]
    known.to_sql("raw_frames", con=connection, index=False, if_exists="append")


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
    """Two biases, a pinhole arc and an order-centre flat: everything a slit-flat SOF needs upstream."""
    biases = raw_frame_table()
    arc = _pinhole_frame(biases, "arc-1", "LAMP,WAVE", "2024-01-02T04:04:00", "2024-01-02T04:04:05")
    orderCentreFlat = _pinhole_frame(biases, "ocflat-1", "LAMP,FLAT", "2024-01-02T05:04:00", "2024-01-02T05:04:05")
    return pd.concat([biases, arc, orderCentreFlat])


@pytest.fixture
def workspace_with_frames(tmp_path, log, monkeypatch):
    """Build a prepared two-session workspace holding the calibration chain plus any flat exposures given."""
    organisers = []

    def build(flatFrames):
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        organiser = workspace_organiser(tmp_path, log=log)
        connection, _ = organiser._get_or_create_db_connection()
        organiser.conn = connection
        _insert_frames(connection, _calibration_chain_frames())
        _insert_frames(connection, flatFrames)
        organiser.session_create("archive")
        organiser.session_create("science")

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


def _late_prepare(organiser, frames):
    """Add ``frames`` to the workspace after its first `prep`, then run a plain `prep` again."""
    _insert_frames(organiser.conn, frames)
    organiser.conn.commit()
    organiser.prepare(report=False)


def _members(organiser, sof):
    """Return the raw frames of one SOF, leaving out the calibration products the SOF map also lists."""
    rows = organiser.conn.execute(
        "SELECT filepath FROM sof_map_science WHERE sof = ? AND filepath LIKE './raw/%'", (sof,)
    ).fetchall()
    return {Path(row[0]).name for row in rows}


def _mflat_sofs(organiser, table):
    rows = organiser.conn.execute(f"SELECT DISTINCT sof FROM {table} WHERE recipe = 'mflat'").fetchall()  # noqa: S608
    return {row[0] for row in rows}


def _seed_status(organiser, sessionId, sof, status):
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE (DY-254)
    organiser.conn.execute(
        f"UPDATE product_frames SET {statusColumn} = ? WHERE sof = ?",  # noqa: S608
        (status, sof),
    )
    organiser.conn.commit()


def _status_row(organiser, sof):
    """Return the distinct (status, status_science, status_archive, complete) tuples of one SOF."""
    return set(
        organiser.conn.execute(
            "SELECT DISTINCT status, status_science, status_archive, complete FROM product_frames WHERE sof = ?",
            (sof,),
        ).fetchall()
    )


def _flat_flats(indexes, **kwargs):
    return raw_flat_run_table(indexes, **kwargs)


# A LATE FRAME REJOINS ITS SET


def test_late_half_set_joins_the_existing_sof_map_entry(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in FIRST_HALF}

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in WHOLE_RUN}


def test_late_half_set_is_written_into_the_existing_sof_file(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    sofPath = Path(organiser.sessionPath) / "sof" / MFLAT_FROM_EXPNO_1

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    listed = sofPath.read_text()
    assert all(f"flat-{n}.fits" in listed for n in WHOLE_RUN)


def test_late_half_set_leaves_one_mflat_sof_holding_every_flat(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    counts = organiser.conn.execute("SELECT DISTINCT counts FROM raw_frame_sets WHERE recipe = 'mflat'").fetchall()
    assert _mflat_sofs(organiser, "product_frames") == {MFLAT_FROM_EXPNO_1}
    assert _mflat_sofs(organiser, "raw_frame_sets") == {MFLAT_FROM_EXPNO_1}
    assert counts == [(10,)]


def test_rejoined_set_is_queued_again_in_every_session(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")
    _seed_status(organiser, "archive", MFLAT_FROM_EXPNO_1, "fail")

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {(None, None, None, 1)}


def test_early_late_frames_keep_the_name_of_the_set_they_join(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(SECOND_HALF))

    # ACT
    # THE REJOINED SET IS NOT FINISHED IN THE FIRST `build_sof_files` PASS, SO THE SECOND PASS MUST KEEP THE PIN
    _late_prepare(organiser, _flat_flats(FIRST_HALF))

    # ASSERT
    productRows = organiser.conn.execute(
        "SELECT COUNT(*) FROM product_frames WHERE recipe = 'mflat' GROUP BY \"eso pro catg\""
    ).fetchall()
    assert _mflat_sofs(organiser, "product_frames") == {MFLAT_FROM_EXPNO_6}
    assert _members(organiser, MFLAT_FROM_EXPNO_6) == {f"flat-{n}.fits" for n in WHOLE_RUN}
    assert productRows == [(1,)] * PRODUCT_TYPES_OF_MFLAT


# A LATE FRAME THAT DOES NOT BELONG TO THE SET STAYS OUT OF IT


@pytest.mark.parametrize(
    "lateRun",
    [
        pytest.param({"namePrefix": "night2", "nightOffset": 1, "tplStart": "2024-01-03T06:04:00"}, id="other-night"),
        pytest.param({"namePrefix": "slit2", "slit": "5.0", "tplStart": "2024-01-02T07:04:00"}, id="other-slit"),
    ],
)
def test_second_full_run_with_other_keys_forms_its_own_sof(workspace_with_frames, lateRun) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(WHOLE_RUN))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")

    # ACT
    _late_prepare(organiser, _flat_flats(WHOLE_RUN, **lateRun))

    # ASSERT
    assert len(_mflat_sofs(organiser, "product_frames")) == 2
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in WHOLE_RUN}
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {("pass", "pass", None, 1)}


def test_full_run_on_the_same_night_with_a_new_template_start_forms_its_own_sof(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(WHOLE_RUN))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")

    # ACT
    _late_prepare(organiser, _flat_flats(WHOLE_RUN, namePrefix="later", tplStart="2024-01-02T09:04:00"))

    # ASSERT
    assert len(_mflat_sofs(organiser, "product_frames")) == 2
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in WHOLE_RUN}
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {("pass", "pass", None, 1)}


# WHAT THE RESET DOES


def _grouping_keys(organiser):
    return [key for key in organiser.filterKeywords if key not in organiser.proKeywords]


def _product_paths(organiser, sof):
    rows = organiser.conn.execute(
        "SELECT filepath FROM product_frames WHERE sof = ? AND file != 'XXXX' ORDER BY filepath", (sof,)
    ).fetchall()
    return [Path(organiser.rootDir) / row[0] for row in rows]


def _touch(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("stale")


def _insert_late_frames(organiser, frames):
    """Insert ``frames`` and give every flat the set_size `build_sof_files` computes before the rejoin step runs."""
    _insert_frames(organiser.conn, frames)
    organiser.conn.execute("UPDATE raw_frames SET set_size = ? WHERE file LIKE 'flat-%'", (len(WHOLE_RUN),))


def _reset_with_late_frames(organiser, frames):
    """Insert ``frames`` late, then run only the reset step and return what it reset and the stale paths it found."""
    _insert_late_frames(organiser, frames)
    rejoined = find_rejoined_sofs(organiser.conn, "sof_map_science", _grouping_keys(organiser))
    stalePaths = reset_rejoined_sofs(
        organiser.conn, "sof_map_science", rejoined, organiser.rootDir, _validate_owned_path, organiser.log
    )
    return rejoined, stalePaths


def _warnings(organiser):
    return [message for level, message in organiser.log.messages if level == "warning"]


def test_reset_clears_every_status_and_frees_the_frames_but_keeps_the_product_rows(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")
    _seed_status(organiser, "archive", MFLAT_FROM_EXPNO_1, "fail")
    organiser.conn.execute("UPDATE product_frames SET status = 'pass' WHERE sof = ?", (MFLAT_FROM_EXPNO_1,))
    rowsBefore = organiser.conn.execute(
        "SELECT COUNT(*) FROM product_frames WHERE sof = ?", (MFLAT_FROM_EXPNO_1,)
    ).fetchone()

    # ACT
    rejoined, _ = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    rowsAfter = organiser.conn.execute(
        "SELECT COUNT(*) FROM product_frames WHERE sof = ?", (MFLAT_FROM_EXPNO_1,)
    ).fetchone()
    processed = organiser.conn.execute("SELECT DISTINCT processed FROM raw_frames WHERE file LIKE 'flat-%'").fetchall()
    assert list(rejoined) == [MFLAT_FROM_EXPNO_1]
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {(None, None, None, 0)}
    assert rowsAfter == rowsBefore
    assert processed == [(0,)]
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == set()
    assert _mflat_sofs(organiser, "raw_frame_sets") == set()


def test_reset_clears_set_first_file_on_the_kept_product_rows(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    organiser.conn.execute(
        "UPDATE product_frames SET set_first_file = 'stale.fits' WHERE sof = ?", (MFLAT_FROM_EXPNO_1,)
    )

    # ACT
    _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    firstFiles = organiser.conn.execute(
        "SELECT DISTINCT set_first_file FROM product_frames WHERE sof = ?", (MFLAT_FROM_EXPNO_1,)
    ).fetchall()
    assert firstFiles == [(None,)]


def test_reset_returns_the_stale_product_files_and_error_logs_without_deleting_them(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    productPaths = _product_paths(organiser, MFLAT_FROM_EXPNO_1)
    errorLogs = [path.with_name(path.stem + "_ERROR.log") for path in productPaths]
    for path in productPaths + errorLogs:
        _touch(path)

    # ACT
    _, stalePaths = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    expected = {path.resolve() for path in productPaths + errorLogs}
    assert productPaths and set(stalePaths) == expected
    assert all(path.exists() for path in productPaths + errorLogs)


def test_prepare_deletes_the_stale_product_files_and_error_logs_of_a_rejoined_set(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    productPaths = _product_paths(organiser, MFLAT_FROM_EXPNO_1)
    errorLogs = [path.with_name(path.stem + "_ERROR.log") for path in productPaths]
    for path in productPaths + errorLogs:
        _touch(path)

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert productPaths and not any(path.exists() for path in productPaths + errorLogs)


def test_reset_does_not_return_a_product_path_outside_the_workspace_and_warns(workspace_with_frames, tmp_path) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    outsideFile = tmp_path / "outside" / "MASTER_FLAT.fits"
    _touch(outsideFile)
    organiser.conn.execute(
        "UPDATE product_frames SET filepath = ? WHERE sof = ?", ("../outside/MASTER_FLAT.fits", MFLAT_FROM_EXPNO_1)
    )

    # ACT
    _, stalePaths = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert stalePaths == []
    assert outsideFile.exists()
    assert any("MASTER_FLAT.fits" in message for message in _warnings(organiser))


def test_prepare_leaves_a_product_path_outside_the_workspace_alone(workspace_with_frames, tmp_path) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    outsideFile = tmp_path / "outside" / "MASTER_FLAT.fits"
    _touch(outsideFile)
    organiser.conn.execute(
        "UPDATE product_frames SET filepath = ? WHERE sof = ?", ("../outside/MASTER_FLAT.fits", MFLAT_FROM_EXPNO_1)
    )

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert outsideFile.exists()


def _point_products_at(organiser, filepath):
    organiser.conn.execute("UPDATE product_frames SET filepath = ? WHERE sof = ?", (filepath, MFLAT_FROM_EXPNO_1))


@pytest.mark.parametrize("filepath", ["./raw/r.fits", "./reduced/../raw/r.fits", "./soxspipe.db"])
def test_reset_does_not_return_a_product_path_inside_the_workspace_but_outside_reduced(
    workspace_with_frames, filepath
) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _touch(Path(organiser.rootDir) / "raw" / "r.fits")
    _point_products_at(organiser, filepath)

    # ACT
    _, stalePaths = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert stalePaths == []
    assert any(Path(filepath).name in message for message in _warnings(organiser))


def test_prepare_leaves_a_raw_file_named_as_a_product_path_alone(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    rawFile = Path(organiser.rootDir) / "raw" / "r.fits"
    _touch(rawFile)
    _point_products_at(organiser, "./raw/r.fits")

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert rawFile.exists()


def test_prepare_does_not_follow_a_product_symlink_to_a_raw_file(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    rawFile = Path(organiser.rootDir) / "raw" / "r.fits"
    _touch(rawFile)
    link = Path(organiser.rootDir) / "reduced" / "2024-01-01" / "soxs-mflat" / "link.fits"
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(rawFile)
    _point_products_at(organiser, "./reduced/2024-01-01/soxs-mflat/link.fits")

    # ACT
    _late_prepare(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert rawFile.exists()
    assert link.is_symlink()


def test_reset_does_not_return_a_product_symlink_and_warns(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    rawFile = Path(organiser.rootDir) / "raw" / "r.fits"
    _touch(rawFile)
    link = Path(organiser.rootDir) / "reduced" / "2024-01-01" / "soxs-mflat" / "link.fits"
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(rawFile)
    _point_products_at(organiser, "./reduced/2024-01-01/soxs-mflat/link.fits")

    # ACT
    _, stalePaths = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert [path.name for path in stalePaths] == ["link_ERROR.log"] * PRODUCT_TYPES_OF_MFLAT
    assert any("link.fits" in message for message in _warnings(organiser))


class _CommitFails:
    """Connection proxy that behaves like the real one except that `commit` fails."""

    def __init__(self, connection):
        self._connection = connection

    def __getattr__(self, name):
        return getattr(self._connection, name)

    def commit(self):
        raise sqlite3.OperationalError("simulated commit failure")


def test_stale_files_survive_when_the_commit_of_the_reset_fails(workspace_with_frames, monkeypatch) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    productPaths = _product_paths(organiser, MFLAT_FROM_EXPNO_1)
    errorLogs = [path.with_name(path.stem + "_ERROR.log") for path in productPaths]
    for path in productPaths + errorLogs:
        _touch(path)
    _insert_late_frames(organiser, _flat_flats(SECOND_HALF))
    dataOrganiserModule = sys.modules["soxspipe.commonutils.data_organiser"]
    realReset = dataOrganiserModule.reset_rejoined_sofs

    def reset_then_break_commit(*args, **kwargs):
        stalePaths = realReset(*args, **kwargs)
        organiser.conn = _CommitFails(organiser.conn)
        return stalePaths

    monkeypatch.setattr(dataOrganiserModule, "reset_rejoined_sofs", reset_then_break_commit)

    # ACT
    with pytest.raises(sqlite3.OperationalError):
        organiser._rejoin_late_frames("sof_map_science")

    # ASSERT
    assert productPaths and all(path.exists() for path in productPaths + errorLogs)


def test_reset_tolerates_product_files_that_are_already_gone(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))

    # ACT
    rejoined, _ = _reset_with_late_frames(organiser, _flat_flats(SECOND_HALF))

    # ASSERT
    assert list(rejoined) == [MFLAT_FROM_EXPNO_1]
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {(None, None, None, 0)}


# FRAMES THAT MUST NOT RESET A SET


def test_a_rejoined_set_keeps_a_status_it_earns_after_the_rejoin(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _late_prepare(organiser, _flat_flats(SECOND_HALF))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {("pass", "pass", None, 1)}


def test_a_third_bias_that_breaks_the_set_size_resets_nothing(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", MBIAS, "pass")
    thirdBias = raw_frame_table().iloc[[0]].copy()
    thirdBias["file"] = "bias-3.fits"
    thirdBias["filepath"] = "./raw/2024-01-01/bias-3.fits"
    thirdBias["eso tpl expno"] = 3
    thirdBias["date-obs"] = "2024-01-02T03:04:07.000"

    # ACT
    _late_prepare(organiser, thirdBias)

    # ASSERT
    assert _status_row(organiser, MBIAS) == {("pass", "pass", None, 1)}


def test_a_late_frame_without_a_template_start_never_joins_a_set(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")
    lateFlats = _flat_flats(SECOND_HALF)
    lateFlats["eso tpl start"] = None

    # ACT
    _late_prepare(organiser, lateFlats)

    # ASSERT
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in FIRST_HALF}
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {("pass", "pass", None, 1)}


def test_a_late_pinhole_arc_sharing_a_template_start_resets_no_set(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", DISP_SOLUTION, "pass")
    lateArc = _pinhole_frame(raw_frame_table(), "arc-2", "LAMP,WAVE", "2024-01-02T04:04:00", "2024-01-02T04:04:09")

    # ACT
    _late_prepare(organiser, lateArc)

    # ASSERT
    assert _status_row(organiser, DISP_SOLUTION) == {("pass", "pass", None, 1)}


# WHERE THE REJOIN RUNS AND WHAT IT LEAVES BEHIND


def _sox_named_flats(expnos):
    """Flats named like ESO files, which the first-file lookup needs, with an `absrot` that identifies the exposure."""
    frames = _flat_flats(expnos, namePrefix="SOXS.flat")
    frames["absrot"] = frames["eso tpl expno"].astype(float)
    return frames


def test_rejoined_set_takes_its_first_file_and_absrot_from_the_new_first_frame(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_sox_named_flats(SECOND_HALF))

    # ACT
    _late_prepare(organiser, _sox_named_flats(FIRST_HALF))

    # ASSERT
    rows = organiser.conn.execute(
        "SELECT DISTINCT set_first_file, absrot FROM product_frames WHERE sof = ?", (MFLAT_FROM_EXPNO_6,)
    ).fetchall()
    assert rows == [("SOXS.flat-1.fits", 1.0)]


def test_a_direct_build_of_sof_files_does_not_rejoin_late_frames(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    _seed_status(organiser, "science", MFLAT_FROM_EXPNO_1, "pass")
    _insert_frames(organiser.conn, _flat_flats(SECOND_HALF))

    # ACT
    organiser.build_sof_files()

    # ASSERT
    assert _members(organiser, MFLAT_FROM_EXPNO_1) == {f"flat-{n}.fits" for n in FIRST_HALF}
    assert _status_row(organiser, MFLAT_FROM_EXPNO_1) == {("pass", "pass", None, 1)}


def test_a_compromised_sof_name_cannot_delete_a_file_outside_the_sof_directory(workspace_with_frames) -> None:
    # ARRANGE
    organiser = workspace_with_frames(_flat_flats(FIRST_HALF))
    outsideFile = Path(organiser.sessionPath).parent / "outside.sof"
    _touch(outsideFile)
    evilSof = "../../outside.sof"
    organiser.conn.execute(
        "INSERT INTO sof_map_science (file, tag, sof, filepath, complete) "
        "VALUES ('x.fits', 'T', ?, './reduced/x.fits', 1)",
        (evilSof,),
    )
    organiser.conn.execute(
        "INSERT INTO product_frames (file, sof, recipe, filepath, complete) "
        "VALUES ('x.fits', ?, 'mflat', './reduced/x.fits', 0)",
        (evilSof,),
    )

    # ACT
    organiser.build_sof_files()

    # ASSERT
    assert outsideFile.exists()
    assert any("outside.sof" in message for message in _warnings(organiser))
