"""Integration contracts for re-queuing the products built from a set that a late raw frame rejoined (DY-264)."""

from __future__ import annotations

from pathlib import Path

import pytest

from tests.factories import flux_standard_stare, raw_flat_run_table, spatial_solution_arc
from tests.integration.test_data_organiser_late_frames import (
    DISP_SOLUTION,
    FIRST_HALF,
    MBIAS,
    MFLAT_FROM_EXPNO_1,
    SECOND_HALF,
    _late_prepare,
    _product_paths,
    _seed_status,
    _status_row,
    _touch,
)

pytestmark = pytest.mark.integration

ORDER_CENTRES = "20240102T050405_VIS_1X1_FAST_OLOC_10_0S_SOXS.sof"
SPAT_SOLUTION = "20240102T070405_VIS_1X1_FAST_SPAT_SOLUTION_10_0S_SOXS.sof"
STARE_STD = "20240102T080100_VIS_1X1_FAST_STARE_STD_FLUX_1_0_60_0S_SOXS.sof"

# THE SPATIAL SOLUTION TAKES THE MASTER FLAT (LEVEL 1); THE STANDARD STAR TAKES THE SPATIAL SOLUTION (LEVEL 2)
CONSUMERS = (SPAT_SOLUTION, STARE_STD)
UNRELATED = (MBIAS, DISP_SOLUTION, ORDER_CENTRES)
CHAIN = (MFLAT_FROM_EXPNO_1, *CONSUMERS)

QUEUED = (None, None, None, 1)
PASSED = ("pass", "pass", None, 1)


@pytest.fixture
def chain_workspace(workspace_with_frames):
    """Build a prepared workspace whose standard star takes the spatial solution but no master-flat product.

    The recipes the data organiser knows give a standard star both the master flat and the spatial solution. The
    standard star's rows for the master flat are removed from the session's SOF map, so only a walk through the
    spatial solution can reach it from the master flat.
    """

    def build():
        organiser = workspace_with_frames(
            raw_flat_run_table(FIRST_HALF), extraFrames=[spatial_solution_arc(), flux_standard_stare()]
        )
        organiser.conn.execute(
            "DELETE FROM sof_map_science WHERE sof = ? AND filepath LIKE '%/soxs-mflat/%'", (STARE_STD,)
        )
        organiser.conn.commit()
        return organiser

    return build


def _seed_all(organiser, sofs, status, sessionId="science"):
    for sof in sofs:
        _seed_status(organiser, sessionId, sof, status)


def _touch_products(organiser, sofs):
    """Write a stale product and a stale `_ERROR.log` for every SOF, and return every path written."""
    written = []
    for sof in sofs:
        for product in _product_paths(organiser, sof):
            for path in (product, product.with_name(product.stem + "_ERROR.log")):
                _touch(path)
                written.append(path)
    return written


def test_a_late_frame_queues_the_rejoined_set_and_both_levels_of_its_consumers(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_all(organiser, CHAIN, "pass")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_status_row(organiser, sof) for sof in CHAIN] == [{QUEUED}] * len(CHAIN)


def test_a_late_frame_deletes_the_stale_products_and_error_logs_of_the_consumers(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    stale = _touch_products(organiser, CONSUMERS)

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert stale and not any(path.exists() for path in stale)


def test_a_consumer_that_failed_is_queued_again(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_status(organiser, "science", SPAT_SOLUTION, "pass")
    _seed_status(organiser, "science", STARE_STD, "fail")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert _status_row(organiser, STARE_STD) == {QUEUED}


def test_a_consumer_of_a_failed_consumer_is_queued_again_with_its_stale_product_deleted(chain_workspace) -> None:
    # ARRANGE: THE FAILED SPATIAL SOLUTION MAKES THE ORGANISER DROP THE STANDARD STAR FROM THE SOF MAP BEFORE THE REJOIN
    organiser = chain_workspace()
    _seed_status(organiser, "science", SPAT_SOLUTION, "fail")
    _seed_status(organiser, "science", STARE_STD, "pass")
    stale = _touch_products(organiser, [STARE_STD])

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert _status_row(organiser, STARE_STD) == {QUEUED}
    assert stale and not any(path.exists() for path in stale)


def test_a_consumer_keeps_its_sof_name_and_its_raw_frames(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    query = "SELECT DISTINCT sof, file FROM sof_map_science WHERE file LIKE 'std-%' OR file LIKE 'spat-%'"
    framesBefore = organiser.conn.execute(query).fetchall()

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    framesAfter = organiser.conn.execute(query).fetchall()
    assert framesBefore and sorted(framesAfter) == sorted(framesBefore)


def test_a_product_that_does_not_use_the_rejoined_set_keeps_its_status_and_its_files(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_all(organiser, UNRELATED, "pass")
    products = _touch_products(organiser, UNRELATED)

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_status_row(organiser, sof) for sof in UNRELATED] == [{PASSED}] * len(UNRELATED)
    assert products and all(path.exists() for path in products)


def test_the_consumers_of_another_session_keep_their_status(chain_workspace) -> None:
    # ARRANGE
    organiser = chain_workspace()
    _seed_all(organiser, CONSUMERS, "pass", sessionId="archive")

    # ACT
    _late_prepare(organiser, raw_flat_run_table(SECOND_HALF))

    # ASSERT
    assert [_status_row(organiser, sof) for sof in CONSUMERS] == [{(None, None, "pass", 1)}] * len(CONSUMERS)


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
    sofs = CHAIN + UNRELATED
    _seed_all(organiser, sofs, "pass")
    products = _touch_products(organiser, sofs)

    # ACT
    organiser.prepare(report=False)

    # ASSERT
    assert [_status_row(organiser, sof) for sof in sofs] == [{PASSED}] * len(sofs)
    assert all(path.exists() for path in products)
