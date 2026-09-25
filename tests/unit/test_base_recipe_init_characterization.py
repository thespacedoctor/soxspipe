"""Characterization of `base_recipe.__init__`.

These tests pin what the constructor does today, before it is split. It is 205
lines long and the lines the offline suite never reaches are the guards and the
session bookkeeping: the standard-star recipe renaming, the two refusals to
overwrite an existing product, and the workspace database block that reads a
recipe's previous status and marks it failed before the recipe runs.

Every test isolates the constructor from the user's workspace: the session
lookup, the product-path prediction and the recipe logger are all replaced, and
the database is a file inside `tmp_path`.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

# THE SESSION NAME THE FAKE ORGANISER REPORTS. THE CONSTRUCTOR BUILDS A COLUMN
# NAME OUT OF IT, SO THE DATABASE FIXTURE AND THE ASSERTIONS SHARE IT.
TEST_SESSION = "base"


def _isolate(
    monkeypatch: pytest.MonkeyPatch,
    tmpPath: Path,
    *,
    productPath: Path,
    currentSession: str | bool = False,
) -> None:
    """Replace every boundary `__init__` reaches outside its own settings."""

    class IsolatedOrganiser:
        """Session lookup boundary with no external workspace state."""

        def __init__(self, **_: object) -> None:
            return None

        def session_list(self, *, silent: bool) -> tuple[str | bool, list[str]]:
            assert silent is True
            return currentSession, []

        def close(self) -> None:
            return None

    monkeypatch.setattr(commonutils, "data_organiser", IsolatedOrganiser)
    monkeypatch.setattr(
        toolkit,
        "predict_product_path",
        lambda *_: (str(productPath), "2024-01-02"),
    )
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmpPath / "qc"), str(tmpPath / "products")),
    )


def _settings(tmpPath: Path, **overrides: Any) -> dict[str, Any]:
    """Return the minimum settings the constructor reads."""
    return {
        "workspace-root-dir": str(tmpPath),
        "instrument": "soxs",
        **overrides,
    }


def test_a_standard_star_set_of_files_switches_the_recipe_to_its_std_variant(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`_STD_` in the set-of-files name renames the nodding, staring and offset recipes."""
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "reduced" / "std.fits")

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "2024_STD_nod.sof"),
        recipeName="soxs-nod",
        turnOffMP=True,
    )

    # ASSERT
    assert recipe.recipeName == "soxs-nod-std"


def test_a_previous_failure_refuses_to_run_again_without_the_overwrite_flag(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """An `_ERROR.log` beside the product raises `FileExistsError` and prints why."""
    # ARRANGE
    productPath = tmp_path / "reduced" / "synthetic.fits"
    productPath.parent.mkdir(parents=True, exist_ok=True)
    (productPath.parent / "synthetic_ERROR.log").write_text("previous failure\n")
    _isolate(monkeypatch, tmp_path, productPath=productPath)

    # ACT / ASSERT
    with pytest.raises(FileExistsError, match="previously failed"):
        base_recipe(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=str(tmp_path / "synthetic.sof"),
            recipeName="soxs-mbias",
            verbose=True,
            turnOffMP=True,
        )

    assert "run the recipe command with the overwrite flag (-x)" in capsys.readouterr().out


def test_an_existing_product_refuses_to_run_again_without_the_overwrite_flag(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """An existing product raises `FileExistsError` and prints the product name."""
    # ARRANGE
    productPath = tmp_path / "reduced" / "synthetic.fits"
    productPath.parent.mkdir(parents=True, exist_ok=True)
    productPath.write_bytes(b"")
    _isolate(monkeypatch, tmp_path, productPath=productPath)

    # ACT / ASSERT
    with pytest.raises(FileExistsError, match="product of this recipe already exists"):
        base_recipe(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=str(tmp_path / "synthetic.sof"),
            recipeName="soxs-mbias",
            verbose=True,
            turnOffMP=True,
        )

    assert "rerun the pipeline command with the overwrite flag (-x)" in capsys.readouterr().out


def test_the_overwrite_flag_runs_over_an_existing_product(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """With `overwrite` set, an existing product and error log are both ignored."""
    # ARRANGE
    productPath = tmp_path / "reduced" / "synthetic.fits"
    productPath.parent.mkdir(parents=True, exist_ok=True)
    productPath.write_bytes(b"")
    (productPath.parent / "synthetic_ERROR.log").write_text("previous failure\n")
    _isolate(monkeypatch, tmp_path, productPath=productPath)

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        overwrite=True,
        turnOffMP=True,
    )

    # ASSERT
    assert recipe.sofName == "synthetic"


def test_input_that_is_not_a_set_of_files_fails_on_the_unset_night_date(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A frame list, or no frames at all, cannot construct a recipe today.

    `self.startNightDate` is set only on the set-of-files branch, and the
    constructor reads it unconditionally when it sets up the QC and product
    directories. Pinned as it stands: every route into `base_recipe` other
    than a `.sof` path raises `AttributeError`. Reported as a defect.
    """
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "unused.fits")

    # ACT / ASSERT
    with pytest.raises(AttributeError, match="startNightDate"):
        base_recipe(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=["raw_one.fits", "raw_two.fits"],
            recipeName="soxs-mbias",
            turnOffMP=True,
        )


def test_the_session_database_records_the_recipe_as_failed_before_it_runs(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """With a session and a set-of-files name, the constructor reads and then fails the status.

    The previous status is carried on the recipe so a completed recipe can be
    restored, and the stored status is set to `fail` immediately, so a crashed
    recipe leaves a failure behind rather than a stale pass.
    """
    # ARRANGE
    databasePath = tmp_path / "soxspipe.db"
    connection = sqlite3.connect(str(databasePath))
    connection.execute(f"create table product_frames (sof text, status_{TEST_SESSION} text)")
    connection.execute("insert into product_frames values ('synthetic.sof', 'pass')")
    connection.commit()
    connection.close()
    _isolate(
        monkeypatch,
        tmp_path,
        productPath=tmp_path / "reduced" / "synthetic.fits",
        currentSession=TEST_SESSION,
    )

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
    )

    # ASSERT
    assert recipe.currentSession == TEST_SESSION
    assert recipe.status == "pass"
    assert recipe.sessionDb == str(databasePath)

    stored = sqlite3.connect(str(databasePath))
    # THE COLUMN NAME COMES FROM A CONSTANT IN THIS FILE, NOT FROM INPUT, AND
    # SQLITE CANNOT PARAMETERISE A COLUMN NAME.
    storedStatus = stored.execute(
        f"select status_{TEST_SESSION} from product_frames where sof = 'synthetic.sof'"  # noqa: S608
    ).fetchone()[0]
    stored.close()
    assert storedStatus == "fail"
    recipe.conn.close()


def test_a_set_of_files_missing_from_the_database_leaves_the_status_unset(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A status lookup that matches no row warns and carries `None`."""
    # ARRANGE
    databasePath = tmp_path / "soxspipe.db"
    connection = sqlite3.connect(str(databasePath))
    connection.execute(f"create table product_frames (sof text, status_{TEST_SESSION} text)")
    connection.commit()
    connection.close()
    _isolate(
        monkeypatch,
        tmp_path,
        productPath=tmp_path / "reduced" / "synthetic.fits",
        currentSession=TEST_SESSION,
    )

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
    )

    # ASSERT
    assert recipe.status is None
    assert any("`self.status = c.fetchone()['status']` failed" in str(message) for message in log.messages)
    recipe.conn.close()


def test_multiprocessing_turned_off_opens_no_database_connection(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`turnOffMP` skips the session database entirely, session or not."""
    # ARRANGE
    _isolate(
        monkeypatch,
        tmp_path,
        productPath=tmp_path / "reduced" / "synthetic.fits",
        currentSession=TEST_SESSION,
    )

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        turnOffMP=True,
    )

    # ASSERT
    assert recipe.conn is None
    assert recipe.status is None


def test_user_settings_override_the_packaged_advanced_settings(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The merged settings take the user's value for any key both sides define."""
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "unused.fits")

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path, **{"soxs-mbias": "user value"}),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        turnOffMP=True,
    )

    # ASSERT
    assert recipe.settings["soxs-mbias"] == "user value"
    # A KEY ONLY THE PACKAGED FILE DEFINES SURVIVES THE MERGE, WHICH IS WHAT
    # PROVES THE PACKAGED FILE WAS FOUND AND READ AT ALL.
    assert len(recipe.settings) > len(_settings(tmp_path)) + 1


def test_the_constructor_assigns_the_empty_qc_and_product_tables(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A constructed recipe carries both empty tables, in the column order rows append into.

    The literal moved into `_empty_qc_and_product_tables`, so this asserts the
    wiring as well as the order: a constructor that stopped calling the helper,
    or stopped assigning its results, fails here.
    """
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "reduced" / "synthetic.fits")

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        turnOffMP=True,
    )

    # ASSERT
    assert list(recipe.qc.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_unit",
        "qc_order",
        "qc_comment",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]
    assert list(recipe.products.columns) == [
        "soxspipe_recipe",
        "product_label",
        "file_name",
        "file_type",
        "obs_date_utc",
        "reduction_date_utc",
        "file_path",
        "label",
    ]
    assert recipe.qc.empty
    assert recipe.products.empty


def test_the_scratch_directory_sits_under_the_workspace_tmp_directory(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Each recipe gets its own randomly named scratch directory under `tmp`."""
    # ARRANGE
    _isolate(monkeypatch, tmp_path, productPath=tmp_path / "unused.fits")

    # ACT
    recipe = base_recipe(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(tmp_path / "synthetic.sof"),
        recipeName="soxs-mbias",
        turnOffMP=True,
    )

    # ASSERT
    assert Path(recipe.outDir).parent == tmp_path / "tmp"
    assert Path(recipe.outDir).name.isdigit()
