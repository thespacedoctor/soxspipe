"""Contracts for the `stack_before_extraction` setting of `soxs_stare` (DY-1611).

With the setting on (the default), every object frame of the SOF is stacked and extracted once. With it off, each
object frame is reduced and extracted on its own, with its own products, and its header lists only its own raw
frame.
"""

from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from tests.integration.test_observing_recipe_orchestration import _route, _write_prepared_frame
from tests.integration.test_soxs_stare_reduction_characterization import (
    STARE_TECHNIQUE,
    _patch_stare,
    _stare_recipe,
)

pytestmark = pytest.mark.integration

OBJECT_SEEDS = (601, 602, 603)


def _object_routes(tmpPath: Path, dprType: str = "OBJECT") -> dict[frozenset[tuple[str, object]], list[str]]:
    """Route three prepared object frames, named like prepared raw frames."""
    paths = [
        str(
            _write_prepared_frame(
                tmpPath / f"raw-object-{index}_pre.fits",
                seed=seed,
                headerOverrides={"DPR_TYPE": dprType, "DPR_TECH": STARE_TECHNIQUE},
            )
        )
        for index, seed in enumerate(OBJECT_SEEDS, start=1)
    ]
    # THE HARNESS ALWAYS ROUTES ONE OBJECT FRAME; A STANDARD SOF HOLDS NONE
    return {
        _route(DPR_TYPE="OBJECT", DPR_TECH=STARE_TECHNIQUE): [],
        _route(DPR_TYPE=dprType, DPR_TECH=STARE_TECHNIQUE): paths,
    }


def _recipe(log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, stack: bool | None, dprType: str = "OBJECT"):
    recipeSettings: dict[str, object] = {"use_flat": True, "sky-subtraction": {"subtract_sky": False}}
    if stack is not None:
        recipeSettings["stack_before_extraction"] = stack
    recipe = _stare_recipe(log, tmpPath, extraRoutes=_object_routes(tmpPath, dprType), recipeSettings=recipeSettings)
    captured = _patch_stare(recipe, monkeypatch, tmpPath)
    captured["keywords"] = []
    monkeypatch.setattr(
        recipe, "update_fits_keywords", lambda frame, rawFrames=False: captured["keywords"].append(rawFrames)
    )
    return recipe, captured


@pytest.mark.parametrize("stack", [True, None], ids=["setting-on", "setting-absent"])
def test_stacking_extracts_all_object_frames_once(log: Any, tmp_path: Path, monkeypatch, stack) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["stack"]) == 1
    assert len(captured["stack"][0]["frames"]) == 3
    assert len(captured["extractor"]) == 1
    assert captured["extractor"][0]["sofName"] == "OBJECT_VIS"
    # NO FILTER: THE HEADER LISTS EVERY RAW FRAME OF THE SOF
    assert captured["keywords"] == [False]


def test_without_stacking_each_object_frame_is_extracted_on_its_own(log: Any, tmp_path: Path, monkeypatch) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [len(call["frames"]) for call in captured["stack"]] == [1, 1, 1]
    assert [call["sofName"] for call in captured["extractor"]] == ["OBJECT_VIS_F1", "OBJECT_VIS_F2", "OBJECT_VIS_F3"]
    assert captured["keywords"] == [["raw-object-1.fits"], ["raw-object-2.fits"], ["raw-object-3.fits"]]
    assert recipe.sofName == "OBJECT_VIS"
    assert recipe.filenameTemplate == "OBJECT_VIS.fits"


def test_without_stacking_each_frame_gets_its_own_product_names(log: Any, tmp_path: Path, monkeypatch) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)
    templates = []
    monkeypatch.setattr(
        recipe, "_plot_merged_spectrum_qcs", lambda *args, **kwargs: templates.append(recipe.filenameTemplate)
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert templates == ["OBJECT_VIS_F1.fits", "OBJECT_VIS_F2.fits", "OBJECT_VIS_F3.fits"]


def test_without_stacking_the_qc_of_every_frame_is_reported(log: Any, tmp_path: Path, monkeypatch) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)
    reported = []

    def fake_report():
        reported.append(recipe.sofName)
        return pd.DataFrame({"qc_name": ["FLUX"], "sof_name": [recipe.sofName + ".sof"], "qc_flag": ["pass"]})

    monkeypatch.setattr(recipe, "report_output", fake_report)

    # ACT
    _, qcTable = recipe.produce_product()

    # ASSERT
    assert reported == ["OBJECT_VIS_F1", "OBJECT_VIS_F2", "OBJECT_VIS_F3"]
    assert sorted(set(qcTable["sof_name"])) == ["OBJECT_VIS_F1.sof", "OBJECT_VIS_F2.sof", "OBJECT_VIS_F3.sof"]


def test_a_flux_standard_is_stacked_even_without_stacking(log: Any, tmp_path: Path, monkeypatch) -> None:
    """One response curve per SOF: the reducer predicts exactly one, and later SOFs read it."""
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False, dprType="STD,FLUX")
    monkeypatch.setattr(recipe, "_generate_response_curve", lambda *args, **kwargs: False)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["stack"]) == 1
    assert len(captured["stack"][0]["frames"]) == 3
    assert len(captured["extractor"]) == 1


def _qc_database(sofNames):
    connection = sqlite3.connect(":memory:")
    connection.execute("create table quality_control (qc_name text, sof_name text)")
    connection.executemany("insert into quality_control values ('SNR MEDIAN', ?)", [(name,) for name in sofNames])
    return connection


def _remaining_qc_sofs(connection):
    return sorted(row[0] for row in connection.execute("select sof_name from quality_control"))


@pytest.mark.parametrize(
    ("stack", "remaining"),
    [
        (True, ["OBJECT_VIS.sof", "OBJECT_VIS_FX.sof", "OTHER_F1.sof"]),
        (False, ["OBJECT_VIS_F1.sof", "OBJECT_VIS_F2.sof", "OBJECT_VIS_FX.sof", "OTHER_F1.sof"]),
    ],
    ids=["stacked-drops-frame-rows", "per-frame-drops-stack-rows"],
)
def test_qc_rows_left_by_the_other_mode_are_removed(log: Any, tmp_path: Path, monkeypatch, stack, remaining) -> None:
    """A frame's QC row decides its SOF's status, so a stale one from an earlier run in the other mode must go."""
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=stack)
    recipe.conn = _qc_database(
        ["OBJECT_VIS.sof", "OBJECT_VIS_F1.sof", "OBJECT_VIS_F2.sof", "OBJECT_VIS_FX.sof", "OTHER_F1.sof"]
    )
    connection = recipe.conn

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _remaining_qc_sofs(connection) == remaining


def test_a_frame_that_raises_leaves_the_sof_names_restored(log: Any, tmp_path: Path, monkeypatch) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=False)

    def failing_reduction(*args, **kwargs):
        raise RuntimeError("extraction failed")

    monkeypatch.setattr(recipe, "_reduce_object_frames", failing_reduction)

    # ACT
    with pytest.raises(RuntimeError, match="extraction failed"):
        recipe.produce_product()

    # ASSERT
    assert recipe.sofName == "OBJECT_VIS"
    assert recipe.filenameTemplate == "OBJECT_VIS.fits"


@pytest.mark.parametrize(
    ("frameFailures", "expected"),
    [
        ([True, "no trace", False], "no trace"),
        ([True, "first message", "second message"], "first message"),
        ([False, True, False], True),
        ([False, False, False], False),
    ],
    ids=["string-beats-earlier-true", "first-string-wins", "true-kept", "no-failure"],
)
def test_without_stacking_a_failure_message_is_not_lost_to_a_bare_failure(
    log: Any, tmp_path: Path, monkeypatch, frameFailures, expected
) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)
    failures = iter(frameFailures)
    monkeypatch.setattr(recipe, "_reduce_object_frames", lambda *args, **kwargs: ("product.fits", next(failures)))

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["cleanUp"] == [{"forceFail": expected}]


def test_non_stare_fallback_frames_are_always_stacked(log: Any, tmp_path: Path, monkeypatch) -> None:
    """Only stare frames are extracted one by one; the PAE flat-lamp fallback is stacked as before."""
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)
    monkeypatch.setattr(
        recipe, "_read_stare_object_frames", lambda kw: ([_flat_lamp(seed) for seed in OBJECT_SEEDS], ["a", "b", "c"])
    )

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(captured["stack"]) == 1
    assert len(captured["stack"][0]["frames"]) == 3


def _flat_lamp(seed):
    from tests.factories import synthetic_ccd

    return synthetic_ccd(
        seed=seed, prepared=True, headerOverrides={"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,PINHOLE"}
    )


def test_a_sky_fallback_on_one_frame_does_not_carry_to_the_next(log: Any, tmp_path: Path, monkeypatch) -> None:
    """The sky-model fallback edits the recipe settings; each frame must start from the settings the run began with."""
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=False)
    recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = True
    seen = []

    def reduce_and_fall_back(frames, rawFrames=False):
        seen.append((recipe.subtractSky, recipe.recipeSettings["sky-subtraction"]["subtract_sky"]))
        recipe.subtractSky = False
        recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = False
        return "product.fits", False

    monkeypatch.setattr(recipe, "_reduce_object_frames", reduce_and_fall_back)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert seen == [(True, True), (True, True), (True, True)]
