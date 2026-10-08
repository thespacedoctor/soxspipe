"""A recipe constructed from a list of frames or a directory instead of a `.sof` file (DY-90).

`base_recipe` documents three ways in: a `.sof` path, a list of frame paths and
a directory. Only the `.sof` route used to set the start-of-night date, so the
other two raised `AttributeError` in the constructor. These tests drive the
real constructor over real FITS files written by the factories. Only the
workspace session lookup and the calibrations path are replaced.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pytest
from astropy.time import Time

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes import base_recipe
from tests.factories import pipeline_settings, product_table, qc_table, raw_fits

pytestmark = pytest.mark.unit

TEST_SESSION = "base"


def _mjd(isoTime: str) -> float:
    """Return the MJD of a UTC ISO timestamp."""
    return float(Time(isoTime, scale="utc").mjd)


def _isolate(monkeypatch: pytest.MonkeyPatch, *, currentSession: str | bool = False) -> None:
    """Replace the session lookup, the calibrations path and the recipe logger."""

    class IsolatedOrganiser:
        def __init__(self, **_: object) -> None:
            return None

        def session_list(self, *, silent: bool) -> tuple[str | bool, list[str]]:
            return currentSession, []

        def close(self) -> None:
            return None

    monkeypatch.setattr(commonutils, "data_organiser", IsolatedOrganiser)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)


def _write_frame(directory: Path, name: str, isoTime: str | None, *, arm: str = "VIS", seed: int = 1) -> Path:
    """Write one raw frame observed at `isoTime`, or with no MJD-OBS when it is None."""
    directory.mkdir(parents=True, exist_ok=True)
    overrides: dict[str, object] = {"SEQ_ARM": arm}
    if isoTime is not None:
        overrides["MJD-OBS"] = _mjd(isoTime)
    return raw_fits(directory / name, seed=seed, headerOverrides=overrides)


def _recipe(log: Any, tmp_path: Path, inputFrames: Any, *, turnOffMP: bool = True) -> base_recipe:
    return base_recipe(
        log=log,
        settings=pipeline_settings(tmp_path / "workspace"),
        inputFrames=inputFrames,
        recipeName="soxs-mbias",
        turnOffMP=turnOffMP,
    )


def test_a_frame_list_takes_the_night_date_from_the_earliest_frame(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The earliest MJD-OBS minus the 15 hour offset, not the first frame in the list, sets the night."""
    # ARRANGE: THE LATER FRAME IS LISTED FIRST AND FALLS AFTER THE 15:00 UTC NIGHT BOUNDARY
    _isolate(monkeypatch)
    later = _write_frame(tmp_path / "raw", "later.fits", "2024-01-02T15:30:00", seed=1)
    earlier = _write_frame(tmp_path / "raw", "earlier.fits", "2024-01-02T14:30:00", seed=2)

    # ACT
    recipe = _recipe(log, tmp_path, [str(later), str(earlier)])

    # ASSERT
    assert recipe.startNightDate == "2024-01-01"


def test_a_frame_list_that_straddles_midnight_utc_belongs_to_the_night_that_began_the_day_before(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Frames either side of midnight UTC are one night, dated by the earliest."""
    # ARRANGE
    _isolate(monkeypatch)
    beforeMidnight = _write_frame(tmp_path / "raw", "a.fits", "2024-01-01T23:50:00", seed=1)
    afterMidnight = _write_frame(tmp_path / "raw", "b.fits", "2024-01-02T00:10:00", seed=2)

    # ACT
    recipe = _recipe(log, tmp_path, [str(afterMidnight), str(beforeMidnight)])

    # ASSERT
    assert recipe.startNightDate == "2024-01-01"


def test_a_directory_of_frames_gives_the_same_night_date_as_the_equivalent_list(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The list and directory routes share one date."""
    # ARRANGE
    _isolate(monkeypatch)
    paths = [
        _write_frame(tmp_path / "raw", "a.fits", "2024-01-02T14:30:00", seed=1),
        _write_frame(tmp_path / "raw", "b.fits", "2024-01-02T15:30:00", seed=2),
    ]

    # ACT
    fromList = _recipe(log, tmp_path, [str(path) for path in paths])
    fromDirectory = _recipe(log, tmp_path, str(tmp_path / "raw"))

    # ASSERT
    assert fromDirectory.startNightDate == fromList.startNightDate == "2024-01-01"


def test_the_list_route_and_the_sof_route_date_the_same_instant_alike(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """One implementation of the 15 hour convention serves both routes."""
    # ARRANGE
    _isolate(monkeypatch)
    frame = _write_frame(tmp_path / "raw", "a.fits", "2024-01-02T00:12:34")
    sofRecipe = _recipe(log, tmp_path, str(tmp_path / "20240102T001234_BIAS_VIS.sof"))

    # ACT
    listRecipe = _recipe(log, tmp_path, [str(frame)])

    # ASSERT
    assert listRecipe.startNightDate == sofRecipe.startNightDate == "2024-01-01"


def test_a_frame_without_mjd_obs_is_ignored_while_another_frame_has_one(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Only frames with a readable MJD-OBS date the night."""
    # ARRANGE
    _isolate(monkeypatch)
    undated = _write_frame(tmp_path / "raw", "undated.fits", None, seed=1)
    dated = _write_frame(tmp_path / "raw", "dated.fits", "2024-01-02T20:00:00", seed=2)

    # ACT
    recipe = _recipe(log, tmp_path, [str(undated), str(dated)])

    # ASSERT
    assert recipe.startNightDate == "2024-01-02"


def test_an_empty_frame_list_fails_early_naming_the_missing_frames(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """No frames is a clear `ValueError`, not an `AttributeError` from the night date."""
    # ARRANGE
    _isolate(monkeypatch)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no input frames"):
        _recipe(log, tmp_path, [])


def test_the_default_input_frames_fail_early_naming_the_missing_frames(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A recipe built with no `inputFrames` at all fails the same clear way."""
    # ARRANGE
    _isolate(monkeypatch)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no input frames"):
        base_recipe(
            log=log,
            settings=pipeline_settings(tmp_path / "workspace"),
            recipeName="soxs-mbias",
            turnOffMP=True,
        )


def test_frames_without_mjd_obs_fail_early_naming_the_keyword(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """When no frame can date the night, the error names MJD-OBS."""
    # ARRANGE
    _isolate(monkeypatch)
    undated = _write_frame(tmp_path / "raw", "undated.fits", None)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="MJD-OBS"):
        _recipe(log, tmp_path, [str(undated)])


def test_an_empty_directory_fails_early_naming_the_missing_frames(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A directory without FITS frames is the same clear error."""
    # ARRANGE
    _isolate(monkeypatch)
    (tmp_path / "raw").mkdir()

    # ACT / ASSERT
    with pytest.raises(ValueError, match="no input frames"):
        _recipe(log, tmp_path, str(tmp_path / "raw"))


def test_a_path_that_is_neither_a_sof_file_nor_a_directory_is_rejected_by_name(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """A string that names no directory is a `TypeError` that quotes the path."""
    # ARRANGE
    _isolate(monkeypatch)
    missing = str(tmp_path / "no_such_directory")

    # ACT / ASSERT
    with pytest.raises(TypeError, match="no_such_directory"):
        _recipe(log, tmp_path, missing)


def test_the_list_route_names_the_recipe_after_recipe_arm_and_night(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """`sofName` is `<recipe>_<ARM>_<night>`, and no product path is predicted."""
    # ARRANGE
    _isolate(monkeypatch)
    frame = _write_frame(tmp_path / "raw", "a.fits", "2024-01-02T20:00:00", arm="NIR")

    # ACT
    recipe = _recipe(log, tmp_path, [str(frame)])

    # ASSERT
    assert recipe.sofName == "soxs-mbias_NIR_2024-01-02"
    assert recipe.productPath is False


def test_the_list_route_creates_the_qc_and_product_directories_under_the_night(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """QC and products land under `<workspace>/{qc,reduced}/<night>/<recipe>/`."""
    # ARRANGE
    _isolate(monkeypatch)
    frame = _write_frame(tmp_path / "raw", "a.fits", "2024-01-02T20:00:00")
    workspace = tmp_path / "workspace"

    # ACT
    recipe = _recipe(log, tmp_path, [str(frame)])

    # ASSERT
    assert Path(recipe.qcDir) == workspace / "qc" / "2024-01-02" / "soxs-mbias"
    assert Path(recipe.productDir) == workspace / "reduced" / "2024-01-02" / "soxs-mbias"


def test_two_arms_on_one_night_keep_their_own_qc_rows(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """The QC write deletes by `sof_name`, so the arm in the name stops one arm erasing the other."""
    # ARRANGE
    _isolate(monkeypatch, currentSession=TEST_SESSION)
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    database = sqlite3.connect(str(workspace / "soxspipe.db"))
    database.execute(f"create table product_frames (sof text, status_{TEST_SESSION} text)")
    seedRow = qc_table().assign(qc_order="-1", qc_flag="pass", sof_name="other.sof").drop(columns=["to_header"])
    seedRow.to_sql("quality_control", database, index=False)
    database.commit()
    database.close()

    # ACT
    for arm in ("UVB", "VIS"):
        frame = _write_frame(tmp_path / f"raw_{arm}", "a.fits", "2024-01-02T20:00:00", arm=arm)
        recipe = _recipe(log, tmp_path, [str(frame)], turnOffMP=False)
        recipe.inst = "XSHOOTER"
        recipe.recipeSettings = {}
        recipe.dateObs = "2024-01-02T20:00:00"
        recipe.verbose = True
        recipe.qc = qc_table()
        recipe.products = product_table()
        recipe.report_output()
        recipe.conn.close()

    # ASSERT
    stored = sqlite3.connect(str(workspace / "soxspipe.db"))
    storedSofNames = sorted(row[0] for row in stored.execute("select sof_name from quality_control").fetchall())
    stored.close()
    assert storedSofNames == [
        "other.sof",
        "soxs-mbias_UVB_2024-01-02.sof",
        "soxs-mbias_VIS_2024-01-02.sof",
    ]
