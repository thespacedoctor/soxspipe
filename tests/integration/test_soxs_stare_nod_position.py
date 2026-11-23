"""One nod position of a nodding SOF is extracted under its own names, with its own copy of the settings (DY-2290).

`_extract_position` runs the stare tail's extraction half on the stack of one nod position. The collaborators are the
recorders of the stare characterization tests.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import position_stack, stack_extraction
from tests.factories import synthetic_ccd
from tests.integration.test_soxs_stare_nod_mode import MARKER, SOF_NAME, _sky_result, _tail_inputs, _tail_recipe

pytestmark = pytest.mark.integration

MERGED_LABEL = "EXTRACTED_MERGED_TABLE"


def _stack(position: str, seed: int = 3) -> position_stack:
    """Return the stack of one nod position: its raw frame, its two stacks and an uncovered mask of two columns."""
    notFlattened = synthetic_ccd(seed=seed, prepared=True)
    flattened = synthetic_ccd(seed=seed + 1, prepared=True)
    uncovered = np.zeros(notFlattened.data.shape, dtype=bool)
    uncovered[:, -2:] = True
    return position_stack(position, [synthetic_ccd(seed=seed + 2, prepared=True)], notFlattened, flattened, uncovered)


def test_a_position_is_extracted_under_its_own_names(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe._extract_position(_stack("A"), _tail_inputs(), False)
    recipe._extract_position(_stack("B"), _tail_inputs(), False)

    # ASSERT: THE POSITION GOES BEFORE THE MARKER, AS THE FRAME NUMBER DOES IN A NO-STACK RUN
    assert [call["sofName"] for call in captured["sky"]] == [f"{SOF_NAME}_A{MARKER}", f"{SOF_NAME}_B{MARKER}"]
    assert [call["sofName"] for call in captured["extractor"]] == [f"{SOF_NAME}_A{MARKER}", f"{SOF_NAME}_B{MARKER}"]
    assert [call["filename"] for call in captured["write"]] == [
        f"{SOF_NAME}_A{MARKER}_SKYSUB.fits",
        f"{SOF_NAME}_A{MARKER}_SKYMODEL.fits",
        f"{SOF_NAME}_A{MARKER}_SKYSUB_RESIDUALS.fits",
        f"{SOF_NAME}_B{MARKER}_SKYSUB.fits",
        f"{SOF_NAME}_B{MARKER}_SKYMODEL.fits",
        f"{SOF_NAME}_B{MARKER}_SKYSUB_RESIDUALS.fits",
    ]


def test_the_sof_name_and_filename_template_are_restored_after_a_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _tail_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe._extract_position(_stack("A"), _tail_inputs(), False)

    # ASSERT
    assert recipe.sofName == SOF_NAME
    assert recipe.filenameTemplate == f"{SOF_NAME}{MARKER}.fits"


def test_the_sof_name_and_filename_template_are_restored_when_a_position_raises(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _, _ = _tail_recipe(log, tmp_path, monkeypatch)

    def fail(*_: Any) -> None:
        raise RuntimeError("the sky fit failed")

    monkeypatch.setattr(recipe, "_subtract_sky", fail)

    # ACT / ASSERT
    with pytest.raises(RuntimeError, match="the sky fit failed"):
        recipe._extract_position(_stack("A"), _tail_inputs(), False)
    assert recipe.sofName == SOF_NAME
    assert recipe.filenameTemplate == f"{SOF_NAME}{MARKER}.fits"


def test_a_position_is_extracted_from_its_own_stacks_and_raw_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)
    stack = _stack("B")
    rawFrames = ["b1.fits", "b2.fits"]
    seen: dict[str, Any] = {}
    monkeypatch.setattr(
        recipe,
        "_slit_illumination_exclude_mask",
        lambda combined, flat, twoDMap, objectFrames: seen.update(frames=objectFrames),
    )
    keywordCalls: list[Any] = []
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **kwargs: keywordCalls.append(kwargs["rawFrames"]))

    # ACT
    extraction = recipe._extract_position(stack, _tail_inputs(), rawFrames)

    # ASSERT: THE FLAT-FIELDED STACK IS FITTED AND EXTRACTED, AND THE POSITION'S UNCOVERED PIXELS ARE LEFT OUT
    assert isinstance(extraction, stack_extraction)
    assert captured["sky"][0]["objectFrame"] is stack.flattened
    np.testing.assert_array_equal(captured["sky"][0]["uncoveredMask"], stack.uncovered)
    assert extraction.combinedObject is stack.flattened
    assert seen["frames"] is stack.objectFrames
    assert keywordCalls == [rawFrames]


def test_the_sky_fallback_of_one_position_does_not_leak_into_the_next(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE SKY FIT OF POSITION A GIVES NO MODEL, WHICH SWITCHES SKY SUBTRACTION OFF IN ITS SETTINGS COPY
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)
    results = iter([(None, None, None), _sky_result()])
    stareModule = import_module("soxspipe.recipes.soxs_stare")

    class Sky:
        def __init__(self, **kwargs: object) -> None:
            captured["sky"].append(kwargs)

        def subtract(self) -> tuple[Any, ...]:
            return (*next(results), recipe.qc, recipe.products)

    monkeypatch.setattr(stareModule, "subtract_sky", Sky)
    originalSettings = recipe.recipeSettings

    # ACT
    recipe._extract_position(_stack("A"), _tail_inputs(), False)
    afterA = (recipe.subtractSky, recipe.recipeSettings is originalSettings, dict(originalSettings))
    recipe._extract_position(_stack("B"), _tail_inputs(), False)

    # ASSERT: A RAN WITHOUT SKY AND ITS COPY SAW THE FALLBACK. THE ORIGINAL SETTINGS AND THE SKY FLAG ARE UNTOUCHED,
    # SO B FITS ITS SKY
    extractorSettings = [call["recipeSettings"] for call in captured["extractor"]]
    assert extractorSettings[0]["horne-extraction-slit-length"] == 30.0
    assert "horne-extraction-slit-length" not in extractorSettings[1]
    assert extractorSettings[1]["sky-subtraction"]["subtract_sky"] is True
    assert afterA[0] is True
    assert afterA[1] is True
    assert "horne-extraction-slit-length" not in afterA[2]
    assert originalSettings["sky-subtraction"]["subtract_sky"] is True
    assert recipe.subtractSky is True


def _recipe_with_merged_row(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> tuple[soxs_stare, dict[str, list[dict[str, object]]]]:
    """Return a tail recipe whose extraction adds the merged-table product row, as the Horne extraction does."""
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)

    class Extractor:
        slitEdgeOrders: list = []
        extractionNotes: list = []
        flatPassHandoff = None
        traceFallback = False

        def __init__(self, **kwargs: object) -> None:
            captured["extractor"].append(kwargs)

        def extract(self) -> tuple[Any, pd.DataFrame, pd.DataFrame, dict[str, int], str]:
            row = dict.fromkeys(recipe.products.columns, "x")
            row["product_label"] = MERGED_LABEL
            row["file_name"] = f"{captured['extractor'][-1]['sofName']}_EXTRACTED_MERGED_A.fits"
            products = pd.concat([recipe.products, pd.DataFrame([row])], ignore_index=True)
            return recipe.qc, products, pd.DataFrame({"WAVE": [500.0]}), {"10-11": 11}, "extracted.fits"

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", Extractor)
    return recipe, captured


def test_the_merged_table_row_of_each_position_is_labelled_with_the_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe_with_merged_row(log, tmp_path, monkeypatch)

    # ACT
    recipe._extract_position(_stack("A"), _tail_inputs(), False)
    recipe._extract_position(_stack("B"), _tail_inputs(), False)

    # ASSERT: THE LABEL OF THE COMBINED SPECTRUM STAYS FREE FOR THE FINAL PRODUCT
    merged = recipe.products[recipe.products["product_label"].str.startswith(MERGED_LABEL)]
    assert merged["product_label"].tolist() == [f"{MERGED_LABEL}_A", f"{MERGED_LABEL}_B"]
    assert MERGED_LABEL not in set(recipe.products["product_label"])


def test_a_product_row_of_an_earlier_step_keeps_its_label(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe_with_merged_row(log, tmp_path, monkeypatch)
    earlier = dict.fromkeys(recipe.products.columns, "x")
    earlier["product_label"] = MERGED_LABEL
    recipe.products = pd.concat([recipe.products, pd.DataFrame([earlier])], ignore_index=True)

    # ACT
    recipe._extract_position(_stack("A"), _tail_inputs(), False)

    # ASSERT
    labels = recipe.products["product_label"]
    assert labels[labels.str.startswith(MERGED_LABEL)].tolist() == [MERGED_LABEL, f"{MERGED_LABEL}_A"]
