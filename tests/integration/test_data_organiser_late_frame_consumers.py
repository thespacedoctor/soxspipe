"""Integration contracts for re-queuing the products built from a set that a late raw frame rejoined (DY-264)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
from tests.factories import raw_flat_run_table, raw_frame_table, workspace_organiser
from tests.integration.test_data_organiser_late_frames import (
    FIRST_HALF,
    INSTRUMENT,
    MBIAS,
    MFLAT_FROM_EXPNO_1,
    SECOND_HALF,
    _calibration_chain_frames,
    _insert_frames,
    _late_prepare,
    _pinhole_frame,
)

pytestmark = pytest.mark.integration

DISP_SOLUTION = "20240102T040405_VIS_1X1_FAST_DISP_SOLUTION_10_0S_SOXS.sof"
ORDER_CENTRES = "20240102T050405_VIS_1X1_FAST_OLOC_10_0S_SOXS.sof"
SPAT_SOLUTION = "20240102T070405_VIS_1X1_FAST_SPAT_SOLUTION_10_0S_SOXS.sof"
STARE_STD = "20240102T080100_VIS_1X1_FAST_STARE_STD_FLUX_1_0_60_0S_SOXS.sof"

# THE MFLAT PRODUCTS A SET OF CONSUMERS TAKES, AT TWO LEVELS: THE SPATIAL SOLUTION TAKES THE MFLAT, AND THE STANDARD
# STAR TAKES THE MFLAT AND THE SPATIAL SOLUTION
CONSUMERS = (SPAT_SOLUTION, STARE_STD)
UNRELATED = (MBIAS, DISP_SOLUTION, ORDER_CENTRES)
CHAIN = (MFLAT_FROM_EXPNO_1, *CONSUMERS)


def _spatial_solution_arc():
    """A multi-pinhole arc: the input of the spatial-solution recipe."""
    frame = _pinhole_frame(raw_frame_table(), "spat-1", "LAMP,WAVE", "2024-01-02T07:04:00", "2024-01-02T07:04:05")
    frame["eso dpr tech"] = "ECHELLE,MULTI-PINHOLE"
    return frame


def _flux_standard_stare(count=2):
    """One OB of flux-standard stare exposures: the input of the stare_std recipe."""
    template = raw_frame_table().iloc[0].to_dict()
    rows = []
    for index in range(count):
        name = f"std-{index + 1}.fits"
        rows.append(
            {
                **template,
                "eso dpr catg": "SCIENCE",
                "eso dpr tech": "ECHELLE,SLIT,STARE",
                "eso dpr type": "STD,FLUX",
                "eso obs id": "201",
                # A STANDARD STAR THAT IS NOT IN THE STATIC LIBRARY IS FLAGGED TO IGNORE
                "eso obs name": "GD71",
                "eso obs targ name": "GD71",
                "exptime": 60.0,
                "slit": "1.0",
                "object": "STD",
                "template": "SOXS_slt_cal_Stare",
                "eso tpl name": "SOXS_slt_cal_Stare",
                "eso tpl nexp": count,
                "eso tpl expno": index + 1,
                "eso tpl start": "2024-01-02T08:00:00",
                "file": name,
                "filepath": f"./raw/2024-01-01/{name}",
                "date-obs": f"2024-01-02T08:01:{index:02d}.000",
                "mjd-obs": 60311.34 + index * 0.0001,
                "set_first_file": "std-1.fits",
            }
        )
    return pd.DataFrame(rows)


@pytest.fixture
def chain_workspace(tmp_path, log, monkeypatch):
    """Build a prepared workspace: calibrations, half a slit-flat run, a spatial solution and a standard star."""
    organisers = []

    def build():
        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        organiser = workspace_organiser(tmp_path, log=log)
        connection, _ = organiser._get_or_create_db_connection()
        organiser.conn = connection
        for frames in (
            _calibration_chain_frames(),
            raw_flat_run_table(FIRST_HALF),
            _spatial_solution_arc(),
            _flux_standard_stare(),
        ):
            _insert_frames(connection, frames)
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


def _seed_status(organiser, sessionId, sofs, status):
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    for sof in sofs:
        # COLUMN NAME CANNOT BE BOUND; CHECKED BY validate_sql_identifier ABOVE
        organiser.conn.execute(
            f"UPDATE product_frames SET {statusColumn} = ? WHERE sof = ?",  # noqa: S608
            (status, sof),
        )
    organiser.conn.commit()


def _statuses(organiser, sof):
    """Return the distinct (status_science, status_archive) pairs of one SOF."""
    return set(
        organiser.conn.execute(
            "SELECT DISTINCT status_science, status_archive FROM product_frames WHERE sof = ?", (sof,)
        )
    )


def _product_paths(organiser, sof):
    rows = organiser.conn.execute(
        "SELECT filepath FROM product_frames WHERE sof = ? AND file != 'XXXX' ORDER BY filepath", (sof,)
    ).fetchall()
    return [Path(organiser.rootDir) / row[0] for row in rows]


def _touch_products(organiser, sofs):
    """Write a stale product and a stale `_ERROR.log` for every SOF; return every path written."""
    written = []
    for sof in sofs:
        for product in _product_paths(organiser, sof):
            for path in (product, product.with_name(product.stem + "_ERROR.log")):
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("stale")
                written.append(path)
    return written


def _stale_files(organiser, sofs):
    return [path for sof in sofs for path in _touch_products(organiser, [sof])]


def test_a_late_frame_queues_the_rejoined_set_and_both_levels_of_its_consumers(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "science", CHAIN, "pass")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_statuses(organiser, sof) for sof in CHAIN] == [{(None, None)}] * len(CHAIN)


def test_a_late_frame_deletes_the_stale_products_and_error_logs_of_the_consumers(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    stale = _stale_files(organiser, CONSUMERS)

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert stale and not any(path.exists() for path in stale)


def test_a_consumer_that_failed_is_queued_again(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "science", [SPAT_SOLUTION], "pass")
    _seed_status(organiser, "science", [STARE_STD], "fail")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert _statuses(organiser, STARE_STD) == {(None, None)}


def test_a_consumer_keeps_its_sof_name_and_its_raw_frames(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    framesBefore = organiser.conn.execute(
        "SELECT DISTINCT sof, file FROM sof_map_science WHERE file LIKE 'std-%' OR file LIKE 'spat-%'"
    ).fetchall()

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    framesAfter = organiser.conn.execute(
        "SELECT DISTINCT sof, file FROM sof_map_science WHERE file LIKE 'std-%' OR file LIKE 'spat-%'"
    ).fetchall()
    assert framesBefore and sorted(framesAfter) == sorted(framesBefore)


def test_a_product_that_does_not_use_the_rejoined_set_keeps_its_status_and_its_files(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "science", UNRELATED, "pass")
    products = _stale_files(organiser, UNRELATED)

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_statuses(organiser, sof) for sof in UNRELATED] == [{("pass", None)}] * len(UNRELATED)
    assert products and all(path.exists() for path in products)


def test_the_consumers_of_another_session_keep_their_status(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "archive", CONSUMERS, "pass")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_statuses(organiser, sof) for sof in CONSUMERS] == [{(None, "pass")}] * len(CONSUMERS)


def test_a_consumer_product_outside_the_reduced_directory_is_not_deleted(chain_workspace, tmp_path) -> None:
    # ARRANGE
    organiser = chain_workspace()
    outsideFile = tmp_path / "outside" / "SSOL.fits"
    outsideFile.parent.mkdir(parents=True)
    outsideFile.write_text("precious")
    organiser.conn.execute(
        "UPDATE product_frames SET filepath = ? WHERE sof = ? AND filepath LIKE '%SSOL_10_0S_SOXS.fits'",
        ("../outside/SSOL.fits", SPAT_SOLUTION),
    )

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert outsideFile.read_text() == "precious"


def test_a_consumer_product_that_is_a_symlink_is_not_deleted(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    target = Path(organiser.rootDir) / "keep.fits"
    target.write_text("precious")
    link = _product_paths(organiser, STARE_STD)[0]
    link.parent.mkdir(parents=True, exist_ok=True)
    link.symlink_to(target)

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert link.is_symlink() and target.read_text() == "precious"


def test_a_prepare_with_no_late_frame_leaves_every_status_and_product_alone(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "science", CHAIN + UNRELATED, "pass")
    products = _stale_files(organiser, CHAIN + UNRELATED)

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert [_statuses(organiser, sof) for sof in CHAIN + UNRELATED] == [{("pass", None)}] * len(CHAIN + UNRELATED)
    assert all(path.exists() for path in products)
