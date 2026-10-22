"""`reduce all` must not reduce a per-frame stare SOF again once it has passed (DY-1705).

A SOF reduced one frame at a time writes `<sof>_F<N>_...` products and never the `<sof>_EXTRACTED_MERGED.fits` that
`toolkit.predict_product_path` predicts, so the product-exists guard in `base_recipe` cannot see it. The session
status in `product_frames` is what marks the SOF as done.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pytest

from soxspipe.recipes.base_recipe import base_recipe
from tests.factories.dataframes import qc_table
from tests.integration.test_reducer_selection import (
    SESSION_ID,
    STATUS_COLUMN,
    _create_selection_database,
    _uninitialized_reducer,
)
from tests.integration.test_soxs_stare_stacking import _recipe

pytestmark = pytest.mark.integration

SOF_NAME = "OBJECT_VIS"


def _workspace_database(databasePath: Path) -> None:
    """Create the selection tables with the stare SOF as a complete set that has not been reduced."""
    _create_selection_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        connection.execute("CREATE TABLE quality_control (qc_name TEXT, sof_name TEXT)")
        connection.execute("ALTER TABLE product_frames ADD COLUMN error_message TEXT")
        connection.execute("INSERT INTO raw_frame_sets VALUES (1, 1, 'stare_obj', ?, 'VIS')", (f"{SOF_NAME}.sof",))
        connection.execute(
            "INSERT INTO product_frames (sof, file, complete) VALUES (?, ?, 1)", (f"{SOF_NAME}.sof", "unused.fits")
        )


def _reduce_per_frame(log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, databasePath: Path, **recipeArgs: Any):
    """Run the stare recipe frame by frame over the workspace database, with its real `clean_up`."""
    recipe, captured = _recipe(log, tmpPath, monkeypatch, stack=False, **recipeArgs)
    recipe.conn = sqlite3.connect(databasePath, autocommit=True)
    recipe.currentSession = SESSION_ID
    recipe.status = None
    recipe.outDir = str(tmpPath / "recipe-out")
    passingQc = qc_table().assign(qc_flag="pass")
    monkeypatch.setattr(recipe, "report_output", lambda: passingQc)
    monkeypatch.setattr(recipe, "clean_up", lambda **kwargs: base_recipe.clean_up(recipe, **kwargs))
    recipe.produce_product()
    return captured


def test_a_per_frame_stare_sof_that_passed_is_not_selected_by_reduce_all_again(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    databasePath = tmp_path / "session.db"
    _workspace_database(databasePath)
    collection = _uninitialized_reducer(log, databasePath)
    assert collection.select_sof_files_to_process(reductionTarget="all")["sof"].tolist() == [f"{SOF_NAME}.sof"]

    # ACT
    captured = _reduce_per_frame(log, tmp_path, monkeypatch, databasePath)
    selected = collection.select_sof_files_to_process(reductionTarget="all")

    # ASSERT
    assert [call["sofName"] for call in captured["extractor"]] == [f"{SOF_NAME}_F{n}" for n in (1, 2, 3)]
    assert selected.empty


def test_a_per_frame_stare_sof_that_failed_is_selected_by_reduce_all_again(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    databasePath = tmp_path / "session.db"
    _workspace_database(databasePath)
    with sqlite3.connect(databasePath) as connection:
        # THE RECIPE MARKS ITS SOF AS FAILED WHEN IT STARTS, AND ONLY CLEAN-UP CAN TURN THAT INTO A PASS
        connection.execute(f"update product_frames set {STATUS_COLUMN} = 'fail'")  # noqa: S608
    collection = _uninitialized_reducer(log, databasePath)
    recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=False)
    recipe.conn = sqlite3.connect(databasePath, autocommit=True)
    recipe.currentSession = SESSION_ID
    recipe.status = None
    recipe.outDir = str(tmp_path / "recipe-out")
    recipe.qc = qc_table().assign(qc_flag="pass")

    # ACT
    base_recipe.clean_up(recipe, forceFail="frame 2 has no trace")
    selected = collection.select_sof_files_to_process(reductionTarget="all")

    # ASSERT
    with sqlite3.connect(databasePath) as connection:
        status = connection.execute(f"select {STATUS_COLUMN} from product_frames").fetchone()[0]  # noqa: S608
    assert status == "fail"
    assert selected["sof"].tolist() == [f"{SOF_NAME}.sof"]
