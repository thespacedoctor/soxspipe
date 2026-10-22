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


def _recipe(
    log: Any,
    tmpPath: Path,
    monkeypatch: pytest.MonkeyPatch,
    stack: bool | None,
    dprType: str = "OBJECT",
    noStack: bool = False,
):
    recipeSettings: dict[str, object] = {"use_flat": True, "sky-subtraction": {"subtract_sky": False}}
    if stack is not None:
        recipeSettings["stack_before_extraction"] = stack
    recipe = _stare_recipe(
        log,
        tmpPath,
        extraRoutes=_object_routes(tmpPath, dprType),
        recipeSettings=recipeSettings,
        settingsOverrides={"no-stack": True} if noStack else None,
    )
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
    # THE PRODUCTION CONNECTION IS AUTOCOMMIT
    connection = sqlite3.connect(":memory:", autocommit=True)
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


def test_with_no_stack_a_science_sof_is_extracted_frame_by_frame_even_when_the_setting_is_on(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=True, noStack=True)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [call["sofName"] for call in captured["extractor"]] == ["OBJECT_VIS_F1", "OBJECT_VIS_F2", "OBJECT_VIS_F3"]


def _standard_with_no_stack(log: Any, tmp_path: Path, monkeypatch, *, noStack: bool = True, responseFailure=False):
    """Return a three-frame flux standard whose response stub records the SOF name it is built under."""
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False, dprType="STD,FLUX", noStack=noStack)
    captured["responses"] = []
    captured["products"] = []
    failures = iter(responseFailure if isinstance(responseFailure, list) else [responseFailure] * 4)

    def fake_response(*args, **kwargs):
        captured["responses"].append(recipe.sofName)
        return next(failures)

    monkeypatch.setattr(recipe, "_generate_response_curve", fake_response)
    reduce_object_frames = recipe._reduce_object_frames

    def recording_reduction(*args, **kwargs):
        _, failure = reduce_object_frames(*args, **kwargs)
        path = f"product-{recipe.sofName}.fits"
        captured["products"].append(path)
        return path, failure

    monkeypatch.setattr(recipe, "_reduce_object_frames", recording_reduction)
    return recipe, captured


def test_with_no_stack_a_flux_standard_is_stacked_first_then_extracted_frame_by_frame(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """The stacked pass writes the SOF-level response later science SOFs read; the frames follow with their own."""
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [len(call["frames"]) for call in captured["stack"]] == [3, 1, 1, 1]
    expected = ["OBJECT_VIS", "OBJECT_VIS_F1", "OBJECT_VIS_F2", "OBJECT_VIS_F3"]
    assert [call["sofName"] for call in captured["extractor"]] == expected
    assert captured["responses"] == expected
    assert recipe.sofName == "OBJECT_VIS"
    assert recipe.filenameTemplate == "OBJECT_VIS.fits"


def test_with_no_stack_a_flux_standard_returns_the_stacked_product_path(log: Any, tmp_path: Path, monkeypatch) -> None:
    """The reducer predicts the SOF-level product, so that is the path it gets back."""
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch)

    # ACT
    productPath, _ = recipe.produce_product()

    # ASSERT
    assert productPath == "product-OBJECT_VIS.fits"
    assert captured["products"][-1] == "product-OBJECT_VIS_F3.fits"


def test_with_no_stack_a_flux_standard_reports_the_qcs_of_both_passes_and_cleans_up_once(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch)

    def fake_report():
        return pd.DataFrame({"qc_name": ["FLUX"], "sof_name": [recipe.sofName + ".sof"], "qc_flag": ["pass"]})

    monkeypatch.setattr(recipe, "report_output", fake_report)

    # ACT
    _, qcTable = recipe.produce_product()

    # ASSERT
    assert sorted(qcTable["sof_name"]) == [
        "OBJECT_VIS.sof",
        "OBJECT_VIS_F1.sof",
        "OBJECT_VIS_F2.sof",
        "OBJECT_VIS_F3.sof",
    ]
    assert captured["cleanUp"] == [{"forceFail": False}]


def test_with_no_stack_a_flux_standard_lists_its_own_raw_frame_in_each_frame_product(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["keywords"] == [False, ["raw-object-1.fits"], ["raw-object-2.fits"], ["raw-object-3.fits"]]


@pytest.mark.parametrize(
    ("responseFailure", "expected"),
    [
        ([False, True, False, False], True),
        ([False, False, False, "no trace"], "no trace"),
        ([True, False, "later message", False], "later message"),
        ([False, False, False, False], False),
    ],
    ids=["frame-fails", "frame-message", "message-beats-stacked-true", "all-pass"],
)
def test_with_no_stack_a_failed_pass_fails_the_standard_sof(
    log: Any, tmp_path: Path, monkeypatch, responseFailure, expected
) -> None:
    """House rule: one failed pass, stacked or per frame, fails the SOF."""
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch, responseFailure=responseFailure)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert captured["cleanUp"] == [{"forceFail": expected}]


def test_with_no_stack_stale_frame_qc_rows_beyond_the_new_frame_count_are_removed(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """An earlier run with more frames left `_F5` rows; through the parent rule they would fail this SOF."""
    # ARRANGE
    recipe, _ = _standard_with_no_stack(log, tmp_path, monkeypatch)
    names = ["OBJECT_VIS.sof", "OBJECT_VIS_F1.sof", "OBJECT_VIS_F5.sof", "OBJECT_VIS_FX.sof", "OTHER_F1.sof"]
    recipe.conn = _qc_database(names)
    connection = recipe.conn

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _remaining_qc_sofs(connection) == ["OBJECT_VIS.sof", "OBJECT_VIS_FX.sof", "OTHER_F1.sof"]


def test_with_no_stack_each_pass_of_a_flux_standard_starts_from_the_settings_the_run_began_with(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """The sky-model fallback edits the recipe settings; neither the stacked pass nor a frame may leak it on."""
    # ARRANGE
    recipe, _ = _standard_with_no_stack(log, tmp_path, monkeypatch)
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
    assert seen == [(True, True)] * 4


def test_without_no_stack_a_flux_standard_is_stacked_only_and_the_log_says_so(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """A settings-file `False` alone keeps today's behaviour for standards."""
    # ARRANGE
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch, noStack=False)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [len(call["frames"]) for call in captured["stack"]] == [3]
    assert captured["responses"] == ["OBJECT_VIS"]
    assert captured["cleanUp"] == [{"forceFail": False}]
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("stack_before_extraction" in message and "flux-standard" in message for message in warnings)


@pytest.mark.parametrize(
    ("noStack", "isStandard", "readsStareFrames", "frameCount", "expected"),
    [
        (True, True, True, 3, True),
        (True, True, True, 2, True),
        (False, True, True, 3, False),
        (True, False, True, 3, False),
        (True, True, False, 3, False),
        (True, True, True, 1, False),
    ],
    ids=["all-met", "two-frames", "no-override", "science", "not-stare-frames", "single-frame"],
)
def test_also_extract_each_standard_frame_needs_the_override_a_standard_and_two_stare_frames(
    log: Any, tmp_path: Path, monkeypatch, noStack, isStandard, readsStareFrames, frameCount, expected
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=False, noStack=noStack)
    recipe.generateReponseCurve = isStandard
    recipe.readStareFrames = readsStareFrames

    # ACT
    result = recipe._also_extract_each_standard_frame([object()] * frameCount)

    # ASSERT
    assert result is expected


def test_with_no_stack_a_stacked_pass_that_raises_leaves_the_settings_restored(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """The sky-model fallback edits the recipe settings; a pass that fails after editing must not leak them."""
    # ARRANGE
    recipe, _ = _standard_with_no_stack(log, tmp_path, monkeypatch)
    recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = True
    recipe.subtractSky = True

    def edit_then_fail(frames, rawFrames=False):
        recipe.subtractSky = False
        recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = False
        raise RuntimeError("stacked pass failed")

    monkeypatch.setattr(recipe, "_reduce_object_frames", edit_then_fail)

    # ACT
    with pytest.raises(RuntimeError, match="stacked pass failed"):
        recipe.produce_product()

    # ASSERT
    assert recipe.subtractSky is True
    assert recipe.recipeSettings["sky-subtraction"]["subtract_sky"] is True


def test_with_no_stack_the_response_science_sofs_read_is_the_stacked_one_never_a_frame_response(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    """The data organiser predicts one response name per standard SOF; a `_F<N>` response must never match it."""
    # ARRANGE
    import yaml

    import soxspipe

    sofMapPath = Path(soxspipe.__file__).parent / "resources" / "soxs_sof_map.yaml"
    products = yaml.safe_load(sofMapPath.read_text())["flux_standard_stare"]["products"]
    rule = next(
        item["replace"][0] for product in products for item in product.values() if item["eso pro catg"] == "RESP_TAB"
    )
    recipe, captured = _standard_with_no_stack(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    predicted = f"{recipe.sofName}.fits".replace(rule["from"], rule["to"])
    written = [f"{sofName}_RESP.fits" for sofName in captured["responses"]]
    assert written[0] == predicted
    assert len(written) == len(OBJECT_SEEDS) + 1
    assert predicted not in written[1:]
    assert len(set(written)) == len(written)


@pytest.mark.parametrize("mode", ["stacked-then-frames", "frames-only"])
def test_the_settings_of_the_run_are_unchanged_after_a_sky_fallback(
    log: Any, tmp_path: Path, monkeypatch, mode: str
) -> None:
    """`recipeSettings` is `settings[recipeName]`; the sky fallback edits it in place and must not reach `settings`."""
    # ARRANGE
    if mode == "stacked-then-frames":
        recipe, _ = _standard_with_no_stack(log, tmp_path, monkeypatch)
    else:
        recipe, _ = _recipe(log, tmp_path, monkeypatch, stack=False)
    recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = True
    recipe.settings[recipe.recipeName] = recipe.recipeSettings
    original = recipe.recipeSettings

    def reduce_and_fall_back(frames, rawFrames=False):
        recipe.subtractSky = False
        recipe.recipeSettings["sky-subtraction"]["subtract_sky"] = False
        return "product.fits", False

    monkeypatch.setattr(recipe, "_reduce_object_frames", reduce_and_fall_back)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.settings[recipe.recipeName]["sky-subtraction"]["subtract_sky"] is True
    assert recipe.recipeSettings is original


def test_without_stacking_the_cosmic_ray_qc_is_reported_under_the_sof_name(
    log: Any, tmp_path: Path, monkeypatch
) -> None:
    # ARRANGE
    recipe, captured = _recipe(log, tmp_path, monkeypatch, stack=False)
    monkeypatch.delattr(recipe, "report_output")
    recipe.verbose = False
    recipe.conn = None
    recipe.inst = "XSHOOTER"
    recipe.cosmicRayProduct = {
        "sofName": recipe.sofName,
        "fileName": f"{recipe.sofName}_LACOSMIC_MASK.fits",
        "filePath": str(tmp_path / "mask.fits"),
        "dateObs": None,
        "count": 7,
    }

    # ACT
    _, qcTable = recipe.produce_product()

    # ASSERT
    cosmicRows = qcTable[qcTable["qc_name"] == "COSMIC RAY PIX NUM"]
    assert cosmicRows["qc_value"].tolist() == [7]
    assert cosmicRows["sof_name"].tolist() == ["OBJECT_VIS.sof"]
