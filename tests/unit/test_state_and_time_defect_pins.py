"""Strict-xfail pins for the state and persistence defects found by the DY-37 defect hunt.

Every test here pins a defect that has been found but not yet filed or fixed. Each test
asserts the CORRECT behaviour, so it fails today and is marked ``xfail(strict=True)``. When
someone fixes the defect the test XPASSes, strict mode turns that into a suite failure, and
the fixer must delete the ``xfail`` marker. If another test file pins the old (wrong)
behaviour as a characterization, the fixer must update or delete that pin in the same change.

The ``reason`` of each marker reads ``<Linear issue> (DY-37 <ID>): <defect>``, naming the
issue that tracks the defect and the defect's ID in the DY-37 findings table.

Each marker names the exception the defect produces today (``AssertionError`` for a wrong
value), so a setup, import or fixture error fails the suite instead of hiding as an xfail.
"""

from __future__ import annotations

import contextlib
import copy
import logging
import shutil
import sqlite3
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.nddata import CCDData

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.data_organiser import data_organiser as real_data_organiser
from soxspipe.commonutils.reducer import reducer
from soxspipe.recipes import base_recipe
from tests.factories import harvestable_raw_fits, workspace_organiser

pytestmark = pytest.mark.unit

# THE SESSION NAME USED BY THE DATABASE FIXTURES BELOW; ITS STATUS COLUMN IS `status_base`
SESSION = "base"


class _IsolatedOrganiser:
    """Session lookup boundary with no external workspace state."""

    def __init__(self, **_: object) -> None:
        return None

    def session_list(self, *, silent: bool) -> tuple[bool, list[str]]:
        return False, []

    def close(self) -> None:
        return None


def _isolate_constructor(monkeypatch: pytest.MonkeyPatch, tmpPath: Path) -> None:
    """Replace every boundary `base_recipe.__init__` reaches outside its settings."""
    monkeypatch.setattr(commonutils, "data_organiser", _IsolatedOrganiser)
    monkeypatch.setattr(base_recipe, "_advanced_settings", lambda self: {})
    monkeypatch.setattr(toolkit, "start_night_date_from_frames", lambda **_: ("2024-01-02", "VIS"))
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmpPath / "qc"), str(tmpPath / "products")))


def _recipe_from(settings: dict[str, Any], recipeName: str, arm: str, log: Any) -> base_recipe:
    """Construct a real recipe from a shared settings dictionary, as the serial reducer does."""
    recipe = base_recipe(log=log, settings=settings, inputFrames=[], recipeName=recipeName, turnOffMP=True)
    recipe.arm = arm
    return recipe


# ---------------------------------------------------------------------------
# S1: RECIPE SETTINGS LEAK BETWEEN RECIPES THAT SHARE ONE SETTINGS DICTIONARY
# ---------------------------------------------------------------------------


def test_a_nir_recipe_does_not_hand_its_arm_settings_to_a_later_vis_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _isolate_constructor(monkeypatch, tmp_path)
    settings = {
        "workspace-root-dir": str(tmp_path),
        "instrument": "soxs",
        "soxs-stare-std": {
            "vis": {"lacosmic-sigma": 25},
            "nir": {"lacosmic-sigma": 25, "qc-acceptable-ranges": {"detector-temp": [0, 50]}},
        },
    }
    snapshot = copy.deepcopy(settings)
    nir = _recipe_from(settings, "soxs-stare-std", "NIR", log)
    if "qc-acceptable-ranges" not in nir.get_recipe_settings():
        pytest.fail("the NIR recipe did not read its own arm's qc-acceptable-ranges")

    # ACT
    vis = _recipe_from(settings, "soxs-stare-std", "VIS", log)
    visSettings = vis.get_recipe_settings()

    # ASSERT
    assert "qc-acceptable-ranges" not in visSettings, "the VIS recipe inherited the NIR-only QC range"
    assert settings == snapshot, "the caller's settings dictionary was changed"


def test_a_recipe_turning_sky_subtraction_off_does_not_turn_it_off_for_the_next_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _isolate_constructor(monkeypatch, tmp_path)
    settings = {
        "workspace-root-dir": str(tmp_path),
        "instrument": "soxs",
        "soxs-stare": {"sky-subtraction": {"subtract_sky": True}, "vis": {}},
    }
    first = _recipe_from(settings, "soxs-stare", "VIS", log)
    first.recipeSettings = first.get_recipe_settings()

    # ACT: WHAT soxs_stare.py DOES ON ITS NO-SKY-MODEL FALLBACK
    first.recipeSettings["sky-subtraction"]["subtract_sky"] = False
    second = _recipe_from(settings, "soxs-stare", "VIS", log)

    # ASSERT
    assert second.get_recipe_settings()["sky-subtraction"]["subtract_sky"] is True


# ---------------------------------------------------------------------------
# W3: SAVED PREPARED FRAMES LAND IN THE WORKSPACE ROOT AND PREP INDEXES THEM AS RAW
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1703 (DY-37 W3): a saved _pre frame is written to the workspace root and the next prep indexes it as raw"
    ),
)
def test_a_saved_prepared_frame_is_not_indexed_as_a_raw_frame(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    harvestable_raw_fits(root / "frame1.fits")
    monkeypatch.chdir(root)
    try:
        organiser.prepare(report=False)
    finally:
        organiser.close()
    rawPath = next(root.rglob("frame1.fits"))
    with fits.open(rawPath) as hdus:
        frame = CCDData(hdus[0].data.astype("float32"), unit="adu", header=hdus[0].header.copy())
    frame.mask = np.zeros(frame.data.shape, dtype=bool)
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame(columns=["qc_name", "qc_value", "qc_comment", "qc_unit", "to_header", "qc_order"])
    recipe.workspaceRootPath = str(root)
    recipe.outDir = str(root / "tmp" / "scratch")
    recipe.sofName = None
    recipe.settings = {}

    # WHAT `save-intermediate-products: True` DOES FOR EVERY PREPARED FRAME
    prePath = Path(recipe._write_prepared_frame(frame, str(rawPath), save=True))
    if not prePath.is_file():
        pytest.fail(f"_write_prepared_frame(save=True) wrote no file at {prePath}")

    # ACT
    again = real_data_organiser(log=log, rootDir=str(root))
    try:
        again.prepare(report=False)
        indexed = pd.read_sql("select file from raw_frames", again.conn)["file"].tolist()
    finally:
        again.close()

    # ASSERT
    assert prePath.name not in indexed, f"the saved prepared frame {prePath.name} was indexed as raw"
    assert indexed == ["frame1.fits"]


# ---------------------------------------------------------------------------
# W1: THE PYTHON API PREPARES RELATIVE TO THE PROCESS CWD, NOT THE WORKSPACE
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason="DY-1702 (DY-37 W1): data_organiser.prepare builds ./raw paths relative to the cwd, not rootDir",
)
def test_prepare_through_the_api_files_frames_inside_the_workspace_whatever_the_cwd(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    organiser = workspace_organiser(tmp_path, log=log)
    root = Path(organiser.rootDir)
    source = tmp_path / "source"
    source.mkdir()
    harvestable_raw_fits(source / "frame1.fits")
    shutil.copy(source / "frame1.fits", root / "frame1.fits")
    monkeypatch.chdir(source)

    # ACT: TODAY PREPARE MOVES THE FRAME INTO <cwd>/raw/ AND THEN CRASHES WITH THIS AttributeError; ANY OTHER
    # EXCEPTION IS NOT THE DEFECT AND FAILS THE SUITE
    crash = None
    try:
        organiser.prepare(report=False)
    except AttributeError as error:
        crash = error
    finally:
        organiser.close()

    # ASSERT
    assert not (source / "raw").exists(), f"prepare filed frames under the cwd (crash: {crash!r})"
    assert list((root / "raw").rglob("frame1.fits")), f"the frame is not in the workspace raw/ (crash: {crash!r})"
    assert crash is None, f"prepare filed the frame but still crashed: {crash!r}"


# ---------------------------------------------------------------------------
# I1: AN INTERRUPTED PRODUCT WRITE LEAVES A PARTIAL FILE THAT EVERY RERUN TREATS AS DONE
# ---------------------------------------------------------------------------

I1_SOF_NAME = "20240102T030405_VIS_1X1_1_MBIAS_SOXS"


def _product_recipe(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, truncate: bool) -> tuple[base_recipe, Path]:
    """Write a product FITS at the predicted path (optionally cut in half) and return a recipe aimed at it."""
    productPath = tmp_path / "reduced" / "2024-01-01" / "soxs-mbias" / f"{I1_SOF_NAME}.fits"
    productPath.parent.mkdir(parents=True)
    fits.HDUList([fits.PrimaryHDU(np.zeros((256, 256), dtype="float32"))]).writeto(productPath, checksum=True)
    if truncate:
        fullSize = productPath.stat().st_size
        with open(productPath, "r+b") as handle:
            handle.truncate(fullSize // 2)
    monkeypatch.setattr(toolkit, "predict_product_path", lambda *_: (str(productPath), "2024-01-01"))
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)
    recipe = base_recipe.__new__(base_recipe)
    recipe.recipeName = "soxs-mbias"
    recipe.settings = {}
    return recipe, productPath


def test_a_complete_product_still_blocks_a_rerun_without_overwrite(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Control for the I1 pin: a valid product must keep blocking a rerun, before and after the fix."""
    # ARRANGE
    recipe, _ = _product_recipe(tmp_path, monkeypatch, truncate=False)

    # ACT / ASSERT
    with pytest.raises(FileExistsError):
        recipe._resolve_product_path(log, str(tmp_path / "sof" / f"{I1_SOF_NAME}.sof"), False, False)


@pytest.mark.xfail(
    strict=True,
    raises=FileExistsError,
    reason=(
        "DY-1697 (DY-37 I1): a truncated product from an interrupted run is skipped as 'already exists' on every rerun"
    ),
)
def test_a_truncated_product_from_an_interrupted_run_is_not_treated_as_an_existing_product(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, productPath = _product_recipe(tmp_path, monkeypatch, truncate=True)

    # ACT: A HALF-WRITTEN FILE MUST NOT BLOCK THE RERUN AS IF IT WERE THE PRODUCT
    recipe._resolve_product_path(log, str(tmp_path / "sof" / f"{I1_SOF_NAME}.sof"), False, False)

    # ASSERT: THE RERUN PROCEEDS AND WILL WRITE THE PRODUCT AT THE SAME PATH
    assert recipe.productPath == str(productPath)


# ---------------------------------------------------------------------------
# I2: A RECIPE INTERRUPTED AFTER WRITING ITS QC ROWS IS MARKED 'pass' BY THE NEXT PREP
# ---------------------------------------------------------------------------


def _session_database(path: Path) -> None:
    """Create the two tables the status bookkeeping reads, holding one SOF."""
    with contextlib.closing(sqlite3.connect(path)) as conn, conn:
        conn.execute("create table product_frames (sof text, error_message text, status_base text)")
        conn.execute(
            "create table quality_control (soxspipe_recipe text, sof_name text, qc_name text, qc_value text, "
            "qc_value_min real, qc_value_max real, qc_flag text, qc_order text)"
        )
        conn.execute("insert into product_frames (sof) values ('A.sof')")


def test_a_recipe_that_never_reached_clean_up_is_not_marked_pass_by_the_qc_range_reconciliation(
    log: Any, tmp_path: Path
) -> None:
    # ARRANGE: THE REAL CONSTRUCTOR BOOKKEEPING MARKS THE SOF BEFORE THE RECIPE RUNS
    _session_database(tmp_path / "soxspipe.db")
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.currentSession = SESSION
    recipe.sofName = "A"
    recipe.settings = {"workspace-root-dir": str(tmp_path)}
    recipe.conn = None
    recipe._open_session_database(home=str(tmp_path))
    # THE RECIPE WRITES ITS (PASSING) QC ROWS, THEN THE PROCESS IS KILLED BEFORE clean_up
    recipe.conn.execute(
        "insert into quality_control values ('soxs-mflat', 'A.sof', 'FLUX MEDIAN', '1.0', null, null, 'pass', '-1')"
    )
    recipe.conn.close()

    # ACT: WHAT THE NEXT `soxspipe prep` RUNS
    organiser = SimpleNamespace(sessionId=SESSION, settings={})
    with contextlib.closing(sqlite3.connect(tmp_path / "soxspipe.db")) as conn, conn:
        for sqlQuery, sqlParams in real_data_organiser._qc_acceptable_range_queries(organiser):
            conn.execute(sqlQuery, sqlParams)
        status = conn.execute("select status_base from product_frames where sof = 'A.sof'").fetchone()[0]

    # ASSERT
    assert status != "pass"


def test_a_recipe_that_finished_failing_on_qc_range_is_passed_when_the_ranges_now_accept_its_value(
    tmp_path: Path,
) -> None:
    # ARRANGE: A FINISHED RUN LEFT `fail` AND A FAILING QC ROW (VALUE 5.0, THEN-RANGE 0 TO 1)
    _session_database(tmp_path / "soxspipe.db")
    with contextlib.closing(sqlite3.connect(tmp_path / "soxspipe.db")) as conn, conn:
        conn.execute("update product_frames set status_base = 'fail' where sof = 'A.sof'")
        conn.execute(
            "insert into quality_control (soxspipe_recipe, sof_name, qc_name, qc_value, qc_value_min, qc_value_max, "
            "qc_flag, qc_order) values ('soxs-mflat', 'A.sof', 'FLUX MEDIAN', '5.0', 0, 1, 'fail', '-1')"
        )
    widenedRanges = {"soxs-mflat": {"qc-acceptable-ranges": {"flux-median": [0, 10]}}}

    # ACT: THE NEXT `soxspipe prep` RUNS WITH THE WIDENED RANGE
    organiser = SimpleNamespace(sessionId=SESSION, settings=widenedRanges)
    with contextlib.closing(sqlite3.connect(tmp_path / "soxspipe.db")) as conn, conn:
        for sqlQuery, sqlParams in real_data_organiser._qc_acceptable_range_queries(organiser):
            conn.execute(sqlQuery, sqlParams)
        status = conn.execute("select status_base from product_frames where sof = 'A.sof'").fetchone()[0]

    # ASSERT
    assert status == "pass"


# ---------------------------------------------------------------------------
# I3: A PREVIOUSLY FAILED RECIPE THAT NOW PASSES NEVER REFRESHES THE DOWNSTREAM SOFS
# ---------------------------------------------------------------------------


class _RecordingCursor:
    def execute(self, *_: object) -> None:
        return None

    def close(self) -> None:
        return None


class _RecordingConnection:
    def cursor(self) -> _RecordingCursor:
        return _RecordingCursor()

    def close(self) -> None:
        return None


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1699 (DY-37 I3): clean_up's 'previously failed, now passing' branch can never run (passToFail needs 'pass')"
    ),
)
def test_a_previously_failed_recipe_that_now_passes_refreshes_the_session(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["FLUX MEDIAN"]})
    recipe.status = "fail"
    recipe.currentSession = SESSION
    recipe.sofName = "A"
    recipe.workspaceRootPath = str(tmp_path)
    recipe.outDir = str(tmp_path / "intermediate")
    Path(recipe.outDir).mkdir()
    recipe.conn = _RecordingConnection()
    refreshes: list[object] = []

    class RecordingOrganiser:
        def __init__(self, **_: object) -> None:
            return None

        def session_refresh(self, *, failure: object = True) -> None:
            refreshes.append(failure)

        def close(self) -> None:
            return None

    monkeypatch.setattr(commonutils, "data_organiser", RecordingOrganiser)

    # ACT
    recipe.clean_up()

    # ASSERT: EXACTLY ONE REFRESH, AND NOT A FAILURE REFRESH
    assert len(refreshes) == 1, f"expected one session refresh, got {refreshes!r}"
    assert refreshes[0] is not True, "the refresh was issued as a failure refresh"


# ---------------------------------------------------------------------------
# I4: THE SERIAL REDUCER NEVER CLEARS SCRATCH LEFT IN <workspace>/tmp/ BY AN INTERRUPTED RECIPE
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1700 (DY-37 I4): only the multiprocess reducer removes <workspace>/tmp/, so interrupted scratch piles up; "
        "related to DY-36 L1 (raise path, tests/unit/test_resource_leak_pins.py); file together, fix must clear "
        "both pins"
    ),
)
def test_the_serial_reducer_removes_scratch_left_by_an_interrupted_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Pin I4: the serial reducer must clear scratch an interrupted recipe left in <workspace>/tmp/.

    Related to DY-36 L1 (raise path, tests/unit/test_resource_leak_pins.py); file together, fix must clear both pins.
    """
    # ARRANGE
    orphan = tmp_path / "tmp" / "tmpinterrupted"
    orphan.mkdir(parents=True)
    (orphan / "SOXS.2024-01-02T03:04:05.678_pre.fits").write_bytes(b"\0" * 2880)

    class StubOrganiser:
        def __init__(self, **_: object) -> None:
            return None

        def session_refresh(self, **_: object) -> None:
            return None

        def close(self) -> None:
            return None

    monkeypatch.setattr(commonutils, "data_organiser", StubOrganiser)
    runner = reducer.__new__(reducer)
    runner.log = log
    runner.settings = {}
    runner.workspaceDirectory = str(tmp_path)
    runner.reductionTarget = "20240102T030405_VIS_1X1_1_MBIAS_SOXS.sof"
    runner.sessionId = SESSION
    runner.sessionPath = str(tmp_path / "sessions" / SESSION)
    runner.sessionDB = str(tmp_path / "soxspipe.db")
    runner.overwrite = False
    runner.quitOnFail = False
    runner.daemon = True
    runner.verbose = False
    runner.pathToSettings = False
    monkeypatch.setattr(
        runner, "select_sof_files_to_process", lambda **_: pd.DataFrame(columns=["recipe", "sof", "command"])
    )

    # ACT
    runner.reduce()

    # ASSERT
    assert not orphan.exists()


# ---------------------------------------------------------------------------
# S2: A FINISHED RECIPE'S LOG HANDLERS STAY ON THE SHARED LOGGER
# ---------------------------------------------------------------------------


@pytest.mark.xfail(
    strict=True,
    raises=AssertionError,
    reason=(
        "DY-1704 (DY-37 S2): clean_up leaves the recipe log handlers on the shared logger; later errors "
        "hit its _ERROR.log"
    ),
)
def test_a_later_error_does_not_land_in_a_finished_recipes_error_log(tmp_path: Path) -> None:
    """Pin S2: clean_up must detach and close the recipe's file handlers.

    The pin expects ``clean_up`` itself to detach them. If the fix lands in the reducer instead, move the pin there.
    """
    # ARRANGE: AN UNREGISTERED LOGGER, SO NOTHING IS LEFT IN THE GLOBAL LOGGER REGISTRY
    sharedLog = logging.Logger("dy37-s2", level=logging.DEBUG)
    productPath = tmp_path / "A" / "A.fits"
    productPath.parent.mkdir()
    try:
        recipe = base_recipe.__new__(base_recipe)
        recipe.log = toolkit.add_recipe_logger(sharedLog, str(productPath))
        recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
        recipe.status = None
        recipe.conn = None
        recipe.outDir = str(tmp_path / "intermediate")
        Path(recipe.outDir).mkdir()
        recipeHandlers = list(recipe.log.handlers)
        if not any(isinstance(handler, logging.FileHandler) for handler in recipeHandlers):
            pytest.fail("add_recipe_logger attached no file handlers")
        recipe.clean_up()

        # ACT: THE REDUCER LOGS THE NEXT RECIPE'S FAILURE ON THE SAME LOGGER
        sharedLog.error("Recipe failed with the following error: a later recipe")
        stillAttached = [handler.get_name() for handler in recipeHandlers if handler in sharedLog.handlers]
        stillOpen = [
            handler.get_name()
            for handler in recipeHandlers
            if isinstance(handler, logging.FileHandler) and handler.stream is not None
        ]
    finally:
        for handler in list(sharedLog.handlers):
            handler.close()
            sharedLog.removeHandler(handler)

    # ASSERT
    assert not stillAttached, f"clean_up left the recipe handlers on the shared logger: {stillAttached}"
    assert not stillOpen, f"clean_up left the recipe file handlers open: {stillOpen}"
    assert not (tmp_path / "A" / "A_ERROR.log").exists()
