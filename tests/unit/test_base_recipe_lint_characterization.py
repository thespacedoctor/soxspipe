"""Characterization tests pinning three lint-flagged expressions in `base_recipe`.

DY-88 commit 5 fixes three Ruff findings in this module: an "unused" local
(F841), two f-string-built SQL statements slated to move to bound parameters,
and a percent-format string (UP031). Each finding looks like a mechanical
rewrite, but the first two are not dead code -- deleting or reshaping them
would change what the recipe does. These tests pin today's behaviour so
commit 5 can be checked against them before touching the expressions.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes import base_recipe
from tests.factories import qc_table, synthetic_ccd

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# 1. `qc_median_flux_level` -- THE F841 `fluxRange` ASSIGNMENT IS A DISGUISED
# HEADER-LOOKUP GUARD, NOT DEAD CODE. DELETING IT WOULD SILENCE THE KeyError
# BELOW WHEN A FRAME'S HEADER HAS NO EXPTIME KEYWORD.
# ---------------------------------------------------------------------------


def _qc_recipe(log: Any) -> base_recipe:
    """Build a recipe carrying only the attributes `qc_median_flux_level` reads."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.kw = lambda keyword: keyword
    return recipe


def test_qc_median_flux_level_raises_keyerror_when_exptime_header_missing(
    log: Any,
) -> None:
    """The unread `fluxRange` line still runs the EXPTIME header lookup, which raises."""
    # ARRANGE
    recipe = _qc_recipe(log)
    frame = synthetic_ccd(shape=(2, 2), prepared=True)
    del frame.header["EXPTIME"]

    # ACT / ASSERT
    with pytest.raises(KeyError):
        recipe.qc_median_flux_level(frame, frameType="MDARK", frameName="synthetic dark")


def test_qc_median_flux_level_returns_the_unmasked_median_when_exptime_is_present(
    log: Any,
) -> None:
    """With EXPTIME present, the method proceeds past the guard and returns today's median."""
    # ARRANGE
    recipe = _qc_recipe(log)
    frame = synthetic_ccd(
        shape=(2, 2),
        prepared=True,
        headerOverrides={"EXPTIME": 10.0},
    )
    frame.data = np.array([[2.0, 4.0], [6.0, 8.0]])
    frame.mask = np.zeros((2, 2), dtype=bool)

    # ACT
    medianFlux = recipe.qc_median_flux_level(
        frame,
        frameType="MDARK",
        frameName="synthetic dark",
    )

    # ASSERT
    assert medianFlux == pytest.approx(5.0)
    metric = recipe.qc.iloc[-1]
    assert metric["qc_name"] == "MDARK MEDIAN"
    assert metric["qc_value"] == pytest.approx(5.0)
    assert metric["qc_comment"] == "[e-] Median flux level of synthetic dark"


# ---------------------------------------------------------------------------
# 2. `clean_up` -- THE EMBEDDED SQL STRINGS BEFORE THEY MOVE TO BOUND `?`
# PARAMETERS. PINS THE EXACT RENDERED QUERY TEXT FOR BOTH BRANCHES.
# ---------------------------------------------------------------------------


class _RecordingCursor:
    """A stub cursor that records every `execute()` call verbatim."""

    def __init__(self, executedQueries: list[tuple[str, tuple[Any, ...] | None]]) -> None:
        self.executedQueries = executedQueries
        self.closed = False

    def execute(self, sqlQuery: str, params: tuple[Any, ...] | None = None) -> None:
        self.executedQueries.append((sqlQuery, params))

    def close(self) -> None:
        self.closed = True


class _RecordingConnection:
    """A stub DB connection that hands out recording cursors.

    Stands in for `sqlite3.Connection`, which returns a distinct `Cursor`
    from each `cursor()` call. Keeping them distinct means closing a cursor
    does not close the connection, so the close order the method under test
    uses stays visible.
    """

    def __init__(self) -> None:
        self.executedQueries: list[tuple[str, tuple[Any, ...] | None]] = []
        self.cursors: list[_RecordingCursor] = []
        self.closed = False

    def cursor(self) -> _RecordingCursor:
        cursor = _RecordingCursor(self.executedQueries)
        self.cursors.append(cursor)
        return cursor

    def close(self) -> None:
        self.closed = True


def _clean_up_recipe(
    log: Any,
    tmp_path: Path,
    *,
    status: str,
) -> tuple[base_recipe, _RecordingConnection]:
    """Build a recipe carrying only the attributes `clean_up` reads, plus its stub connection."""
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    # A PASSING QC ROW, SO `clean_up` NEITHER FLIPS `forceFail` TO `True`
    # NOR TAKES THE `data_organiser` SESSION-REFRESH BRANCH.
    recipe.qc = pd.DataFrame({"qc_flag": ["pass"], "qc_name": ["RON"]})
    recipe.status = status
    recipe.currentSession = "base"
    recipe.sofName = "20240102T030405_VIS_1X1_FAST_MBIAS_SOXS"
    conn = _RecordingConnection()
    recipe.conn = conn
    # A DIRECTORY THAT DOES NOT EXIST, SO `shutil.rmtree` FAILS AND THE
    # METHOD SWALLOWS THE ERROR, MATCHING PRODUCTION WHEN THE SCRATCH
    # DIRECTORY IS ALREADY GONE. DELIBERATELY NOT NAMED "cache".
    recipe.outDir = str(tmp_path / "already-removed")
    return recipe, conn


def test_clean_up_sets_products_to_pass_with_todays_interpolated_sql(
    log: Any,
    tmp_path: Path,
) -> None:
    """A clean pass renders `status_<session> = 'pass'` with the sof name inlined into the SQL text."""
    # ARRANGE
    recipe, conn = _clean_up_recipe(log, tmp_path, status="pass")

    # ACT
    recipe.clean_up(forceFail=False)

    # ASSERT
    assert len(conn.executedQueries) == 1
    sqlQuery, params = conn.executedQueries[0]
    assert sqlQuery == (
        "update product_frames set status_base = 'pass' where sof = '20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof'"
    )
    # NO BOUND PARAMETERS TODAY -- EVERY VALUE IS ALREADY INSIDE `sqlQuery`.
    assert params is None
    # THE CONNECTION OUTLIVES THE CURSOR: `clean_up` CLOSES THE CURSOR IT OPENED,
    # THEN THE CONNECTION.
    assert conn.cursors[0].closed
    assert conn.closed


def test_clean_up_records_the_force_fail_message_with_todays_interpolated_sql(
    log: Any,
    tmp_path: Path,
) -> None:
    """A string `forceFail` renders `error_message = '<message>'` with the message inlined into the SQL text."""
    # ARRANGE
    recipe, conn = _clean_up_recipe(log, tmp_path, status="fail")

    # ACT
    recipe.clean_up(forceFail="boom: something went wrong")

    # ASSERT
    assert len(conn.executedQueries) == 1
    sqlQuery, params = conn.executedQueries[0]
    assert sqlQuery == (
        "update product_frames set error_message = 'boom: something went wrong' "
        "where sof = '20240102T030405_VIS_1X1_FAST_MBIAS_SOXS.sof'"
    )
    assert params is None
    assert conn.cursors[0].closed
    assert conn.closed


# ---------------------------------------------------------------------------
# 3. `_bad_pixel_mask` -- THE UP031 PERCENT-FORMAT MESSAGE. PINS THE FULLY
# RENDERED SENTENCE RAISED WHEN THE BITMAP PATH DOES NOT EXIST.
# ---------------------------------------------------------------------------


def test_bad_pixel_mask_raises_oserror_naming_the_missing_bitmap_path(
    log: Any,
    tmp_path: Path,
) -> None:
    """The missing bitmap path is rendered into the raised `OSError`'s message, verbatim."""
    # ARRANGE
    calibrationRootPath = tmp_path / "calibrations"
    calibrationRootPath.mkdir(parents=True, exist_ok=True)
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.arm = "VIS"
    recipe.detectorParams = {"bad-pixel map": {"1x1": "absent.fits"}}
    recipe.calibrationRootPath = str(calibrationRootPath)
    frame = CCDData(np.zeros((4, 4)), unit=u.adu)

    expectedPath = f"{calibrationRootPath}/absent.fits"
    expectedMessage = f"the path to the bitMapPath {expectedPath} does not exist on this machine"

    # ACT / ASSERT
    with pytest.raises(OSError) as excInfo:
        recipe._bad_pixel_mask(frame)

    assert str(excInfo.value) == expectedMessage
    # THE METHOD WRITES THE DUMMY MAP BEFORE RAISING, SO A LATER RUN WOULD
    # FIND IT ALREADY IN PLACE.
    assert (calibrationRootPath / "absent.fits").exists()


# ---------------------------------------------------------------------------
# 4. `_dataframe_to_sqlite` -- THE OTHER SUPPRESSED SQL SITE. ITS `replace`
# BRANCH INTERPOLATES THE TABLE NAME, WHICH NO BOUND PARAMETER CAN CARRY.
# ---------------------------------------------------------------------------


def test_dataframe_to_sqlite_with_replace_clears_the_table_before_inserting(
    log: Any,
) -> None:
    """`replace=True` deletes the table's existing rows, then writes the new ones."""
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.conn = sqlite3.connect(":memory:")
    recipe.conn.execute("create table quality_control (qc_name text, qc_value text)")
    recipe.conn.execute("insert into quality_control values ('OLD', '1')")
    recipe.conn.commit()
    dataframe = pd.DataFrame({"qc_name": ["RON"], "qc_value": ["--"]})

    # ACT
    recipe._dataframe_to_sqlite(dataframe, "quality_control", replace=True)

    # ASSERT
    rows = recipe.conn.execute("select qc_name, qc_value from quality_control").fetchall()
    # THE OLD ROW IS GONE, AND THE "--" PLACEHOLDER IS WRITTEN AS NULL.
    assert rows == [("RON", None)]
    recipe.conn.close()
