"""Public construction and validation contracts for pipeline recipes."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

from soxspipe.recipes import (
    base_recipe,
    soxs_disp_solution,
    soxs_mbias,
    soxs_mdark,
    soxs_mflat,
    soxs_nod,
    soxs_offset,
    soxs_order_centres,
    soxs_spatial_solution,
    soxs_stare,
    soxs_straighten,
)
from tests.factories import pipeline_settings

pytestmark = pytest.mark.unit


@dataclass
class FakeFrameCollection:
    """Small image-collection stand-in for constructor and validation seams."""

    calls: list[tuple[str, object]] = field(default_factory=list)
    summary: str = "synthetic frame inventory"
    files: list[str] = field(default_factory=list)
    exposureTimes: list[float] = field(default_factory=lambda: [10.0])

    def sort(self, keys: list[str]) -> None:
        self.calls.append(("sort", tuple(keys)))

    def filter(self, **kwargs: object) -> FakeFrameCollection:
        self.calls.append(("filter", kwargs))
        return FakeFrameCollection()

    def values(self, *, keyword: str, unique: bool) -> list[float]:
        self.calls.append(("values", (keyword, unique)))
        return list(self.exposureTimes)


RECIPE_NAMES = (
    (soxs_mbias, "soxs-mbias"),
    (soxs_mdark, "soxs-mdark"),
    (soxs_mflat, "soxs-mflat"),
    (soxs_disp_solution, "soxs-disp-solution"),
    (soxs_order_centres, "soxs-order-centres"),
    (soxs_spatial_solution, "soxs-spat-solution"),
    (soxs_stare, "soxs-stare"),
    (soxs_nod, "soxs-nod"),
    (soxs_offset, "soxs-offset"),
    (soxs_straighten, "soxs-straighten"),
)


@pytest.mark.parametrize(("recipeClass", "recipeName"), RECIPE_NAMES)
def test_constructor_runs_shared_input_workflow_and_forwards_public_options(
    recipeClass: type[base_recipe],
    recipeName: str,
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every public constructor loads, verifies, sorts, and prepares its input."""
    settings = pipeline_settings(tmp_path, overrides={"data-extension": 0})
    frames = FakeFrameCollection()
    baseCalls: list[dict[str, object]] = []
    workflowCalls: list[tuple[str, object]] = []
    frames.sort = lambda keys: workflowCalls.append(("sort", tuple(keys)))

    def fake_base_init(self: base_recipe, **kwargs: object) -> None:
        baseCalls.append(dict(kwargs))
        self.log = kwargs["log"]
        self.settings = kwargs["settings"]
        self.inputFrames = kwargs["inputFrames"]
        self.recipeName = kwargs["recipeName"]
        self.sofName = "synthetic"

    class FakeSetOfFiles:
        def __init__(self, **kwargs: object) -> None:
            workflowCalls.append(("load", kwargs["inputFrames"]))

        def get(self) -> tuple[FakeFrameCollection, dict[str, object]]:
            return frames, {}

    def fake_verify(self: base_recipe) -> None:
        workflowCalls.append(("verify", self.recipeName))

    def fake_prepare(self: base_recipe, *, save: bool) -> FakeFrameCollection:
        workflowCalls.append(("prepare", save))
        return frames

    monkeypatch.setattr(base_recipe, "__init__", fake_base_init)
    monkeypatch.setattr(base_recipe, "get_recipe_settings", lambda self: {})
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare)
    monkeypatch.setattr(recipeClass, "verify_input_frames", fake_verify)
    monkeypatch.setattr(
        "soxspipe.commonutils.set_of_files.set_of_files",
        FakeSetOfFiles,
    )

    recipeClass(
        log=log,
        settings=settings,
        inputFrames=["raw.fits"],
        verbose=True,
        overwrite=True,
        command="recipe command",
        debug=True,
        turnOffMP=True,
    )

    assert baseCalls == [
        {
            "log": log,
            "settings": settings,
            "inputFrames": ["raw.fits"],
            "overwrite": True,
            "recipeName": recipeName,
            "command": "recipe command",
            "debug": True,
            "verbose": True,
            "turnOffMP": True,
        }
    ]
    expectedWorkflow = [
        ("load", ["raw.fits"]),
        ("verify", recipeName),
        ("sort", ("MJD-OBS",)),
        ("prepare", False),
    ]
    if recipeClass is soxs_offset:
        expectedWorkflow = expectedWorkflow * 2
    assert workflowCalls == expectedWorkflow


VALID_INPUTS = (
    (soxs_mbias, "soxs-mbias", ["BIAS"], ["IMAGE"], [], {}),
    (soxs_mdark, "soxs-mdark", ["DARK"], ["IMAGE"], [], {}),
    (
        soxs_mflat,
        "soxs-mflat",
        ["LAMP,FLAT"],
        ["ECHELLE,SLIT"],
        ["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        {},
    ),
    (
        soxs_disp_solution,
        "soxs-disp-solution",
        ["LAMP,FMTCHK"],
        ["ECHELLE,PINHOLE"],
        ["MASTER_BIAS_VIS"],
        {},
    ),
    (
        soxs_order_centres,
        "soxs-order-centres",
        ["FLAT,LAMP"],
        ["ECHELLE,PINHOLE"],
        ["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
        {},
    ),
    (
        soxs_spatial_solution,
        "soxs-spat-solution",
        ["LAMP,WAVE"],
        ["ECHELLE,MULTI-PINHOLE"],
        ["MASTER_BIAS_VIS", "ORDER_TAB_VIS", "DISP_TAB_VIS"],
        {},
    ),
    (
        soxs_stare,
        "soxs-stare",
        ["OBJECT"],
        ["ECHELLE,SLIT,STARE"],
        [
            "MASTER_BIAS_VIS",
            "MASTER_FLAT_VIS",
            "ORDER_TAB_VIS",
            "DISP_TAB_VIS",
            "DISP_IMAGE_VIS",
        ],
        {},
    ),
    (
        soxs_nod,
        "soxs-nod",
        ["OBJECT"],
        ["ECHELLE,SLIT,NODDING"],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS"],
        {},
    ),
    (
        soxs_offset,
        "soxs-offset",
        ["OBJECT"],
        ["ECHELLE,SLIT,OFFSET"],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS"],
        {},
    ),
    (
        soxs_straighten,
        "soxs-straighten",
        ["OBJECT"],
        ["ECHELLE,SLIT"],
        [],
        {"VIS": {"2D_MAP": "dispersion-map.fits"}},
    ),
)


def _validation_recipe(
    recipeClass: type[base_recipe],
    log: Any,
    recipeName: str,
    imageTypes: list[str],
    imageTechniques: list[str],
    imageCategories: list[str],
    supplementaryInput: dict[str, object],
) -> base_recipe:
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.recipeName = recipeName
    recipe.settings = {}
    recipe.inputFrames = FakeFrameCollection()
    recipe.supplementaryInput = supplementaryInput
    recipe._verify_input_frames_basics = lambda: (
        list(imageTypes),
        list(imageTechniques),
        list(imageCategories),
    )
    return recipe


@pytest.mark.parametrize(
    (
        "recipeClass",
        "recipeName",
        "imageTypes",
        "imageTechniques",
        "imageCategories",
        "supplementaryInput",
    ),
    VALID_INPUTS,
)
def test_verify_input_frames_accepts_complete_recipe_inventory(
    recipeClass: type[base_recipe],
    recipeName: str,
    imageTypes: list[str],
    imageTechniques: list[str],
    imageCategories: list[str],
    supplementaryInput: dict[str, object],
    log: Any,
) -> None:
    """Every recipe accepts its minimal complete VIS input inventory."""
    recipe = _validation_recipe(
        recipeClass,
        log,
        recipeName,
        imageTypes,
        imageTechniques,
        imageCategories,
        supplementaryInput,
    )

    result = recipe.verify_input_frames()

    assert result is None
    assert recipe.imageType == imageTypes[0]


INVALID_INPUTS = (
    (*VALID_INPUTS[0], "not BIAS frames"),
    (*VALID_INPUTS[1], "not DARK frames"),
    (*VALID_INPUTS[2], "need to be flat-lamp frames"),
    (*VALID_INPUTS[3], "need to be single pinhole lamp"),
    (*VALID_INPUTS[4], "need to be single pinhole flat-lamp"),
    (*VALID_INPUTS[5], "Found a UNEXPECTED frame"),
    (*VALID_INPUTS[6], "soxspipe stare need to be an object frame"),
    (*VALID_INPUTS[7], "Found a UNEXPECTED file"),
    (*VALID_INPUTS[8], "soxspipe nod need to be an object/std nodding frames"),
    (*VALID_INPUTS[9], "Need a full dispersion/spatial solution for VIS"),
)


@pytest.mark.parametrize(
    (
        "recipeClass",
        "recipeName",
        "imageTypes",
        "imageTechniques",
        "imageCategories",
        "supplementaryInput",
        "errorFragment",
    ),
    INVALID_INPUTS,
)
def test_verify_input_frames_rejects_invalid_recipe_inventory(
    recipeClass: type[base_recipe],
    recipeName: str,
    imageTypes: list[str],
    imageTechniques: list[str],
    imageCategories: list[str],
    supplementaryInput: dict[str, object],
    errorFragment: str,
    log: Any,
) -> None:
    """Every recipe reports a meaningful type error for an invalid inventory."""
    if recipeClass is soxs_straighten:
        invalidTypes = imageTypes
        invalidSupplementary: dict[str, object] = {}
    else:
        invalidTypes = ["UNEXPECTED"]
        invalidSupplementary = supplementaryInput
    recipe = _validation_recipe(
        recipeClass,
        log,
        recipeName,
        invalidTypes,
        imageTechniques,
        imageCategories,
        invalidSupplementary,
    )

    with pytest.raises(TypeError, match=errorFragment):
        recipe.verify_input_frames()


MISSING_CALIBRATION_INPUTS = (
    (*VALID_INPUTS[2][:4], ["ORDER_TAB_VIS"], {}, "master-bias frame"),
    (*VALID_INPUTS[3][:4], [], {}, "a master-bias"),
    (*VALID_INPUTS[4][:4], ["MASTER_BIAS_VIS"], {}, "dispersion solution table"),
    (
        *VALID_INPUTS[5][:4],
        ["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
        {},
        "order location table",
    ),
    (
        *VALID_INPUTS[6][:4],
        ["MASTER_BIAS_VIS"],
        {},
        "missing a DISP_TAB_VIS frame",
    ),
    (
        *VALID_INPUTS[7][:4],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS"],
        {},
        "missing a DISP_IMAGE_VIS frame",
    ),
    (
        *VALID_INPUTS[8][:4],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS"],
        {},
        "missing a DISP_IMAGE_VIS frame",
    ),
    (*VALID_INPUTS[9][:5], {}, "Need a full dispersion/spatial solution for VIS"),
)


@pytest.mark.parametrize(
    (
        "recipeClass",
        "recipeName",
        "imageTypes",
        "imageTechniques",
        "imageCategories",
        "supplementaryInput",
        "errorFragment",
    ),
    MISSING_CALIBRATION_INPUTS,
)
def test_verify_input_frames_rejects_missing_calibration(
    recipeClass: type[base_recipe],
    recipeName: str,
    imageTypes: list[str],
    imageTechniques: list[str],
    imageCategories: list[str],
    supplementaryInput: dict[str, object],
    errorFragment: str,
    log: Any,
) -> None:
    """Recipes report which required calibration is absent."""
    recipe = _validation_recipe(
        recipeClass,
        log,
        recipeName,
        imageTypes,
        imageTechniques,
        imageCategories,
        supplementaryInput,
    )

    with pytest.raises(TypeError, match=errorFragment):
        recipe.verify_input_frames()


INVALID_TECHNIQUE_INPUTS = (
    (
        soxs_mflat,
        "soxs-mflat",
        ["LAMP,FLAT"],
        ["UNEXPECTED"],
        ["ORDER_TAB_NIR"],
        "flat-lamp on and lamp off frames for NIR",
    ),
    (
        soxs_disp_solution,
        "soxs-disp-solution",
        ["LAMP,FMTCHK"],
        ["UNEXPECTED"],
        [],
        "single pinhole lamp on and lamp off frames for NIR",
    ),
    (
        soxs_order_centres,
        "soxs-order-centres",
        ["FLAT,LAMP"],
        ["UNEXPECTED"],
        ["DISP_TAB_NIR"],
        "single pinhole flat-lamp on and lamp off frames",
    ),
    (
        soxs_spatial_solution,
        "soxs-spat-solution",
        ["LAMP,WAVE"],
        ["UNEXPECTED"],
        ["ORDER_TAB_NIR", "DISP_TAB_NIR"],
        "Found a UNEXPECTED file",
    ),
    (
        soxs_stare,
        "soxs-stare",
        ["OBJECT"],
        ["UNEXPECTED"],
        [],
        "missing a UNEXPECTED frame",
    ),
    (
        soxs_nod,
        "soxs-nod",
        ["OBJECT"],
        ["UNEXPECTED"],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS"],
        "soxspipe nod need to be an object/std nodding frames",
    ),
    (
        soxs_offset,
        "soxs-offset",
        ["OBJECT"],
        ["UNEXPECTED"],
        ["ORDER_TAB_VIS", "DISP_TAB_VIS", "DISP_IMAGE_VIS"],
        "soxspipe offset need to be an object/std offset frames",
    ),
)


@pytest.mark.parametrize(
    (
        "recipeClass",
        "recipeName",
        "imageTypes",
        "imageTechniques",
        "imageCategories",
        "errorFragment",
    ),
    INVALID_TECHNIQUE_INPUTS,
)
def test_verify_input_frames_rejects_invalid_technique(
    recipeClass: type[base_recipe],
    recipeName: str,
    imageTypes: list[str],
    imageTechniques: list[str],
    imageCategories: list[str],
    errorFragment: str,
    log: Any,
) -> None:
    """Recipes that define technique contracts reject unsupported techniques."""
    recipe = _validation_recipe(
        recipeClass,
        log,
        recipeName,
        imageTypes,
        imageTechniques,
        imageCategories,
        {},
    )
    if recipeClass not in (soxs_nod, soxs_offset):
        recipe.arm = "NIR"

    with pytest.raises(TypeError, match=errorFragment):
        recipe.verify_input_frames()
