"""Characterization of `soxs_spatial_solution.__init__` and `soxs_spatial_solution.verify_input_frames`.

These tests pin what the recipe does today, before its constructor and its
frame verification are split. They are the recipe-layer construction seam the
tier 3 spec asks for: the orchestration tests build a recipe with `__new__`,
which skips `__init__` entirely, so the constructor this ticket edits is code
the offline suite otherwise never reaches.

The constructor tests call the real `soxs_spatial_solution.__init__`.
Everything the recipe inherits is replaced at its own boundary -- the session
lookup, the product-path prediction, the recipe logger, the basic frame
verification and the frame preparation -- so what runs for real is the
recipe's own constructor and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call. Each rejection branch's exact
rendered message is pinned here, so the line-length rewrites of the long
message strings can be shown to render byte-identical strings.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_spatial_solution import soxs_spatial_solution

pytestmark = pytest.mark.unit

# THE RENDERED REJECTION MESSAGES. THE "Found a" MESSAGES NAME THE LAST
# OFFENDING VALUE THE LOOP SAW, NOT THE FIRST.
NIR_TAIL = (
    "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames, "
    "a first-guess dispersion solution table and an order location table for NIR. "
    "Can optionally supply a master-flat for NIR."
)
NIR_ERROR = NIR_TAIL
UVB_VIS_ERROR = (
    "Input frames for soxspipe spatial_solution need to be LAMP,WAVE, a master-bias, a first-guess "
    "dispersion solution table and an order location table. Can optionally supply a master-flat "
    "and/or master-dark for UVB/VIS."
)


def _nir_type_error(found: str) -> str:
    return (
        f"Found a {found} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE. "
        "Can optionally supply a master-flat for NIR."
    )


def _nir_tech_error(found: str) -> str:
    return f"Found a {found} file. {NIR_TAIL}"


def _uvb_vis_type_error(found: str) -> str:
    return (
        f"Found a {found} frame. Input frames for soxspipe spatial_solution need to be LAMP,WAVE and a "
        "master-bias, a first-guess dispersion solution table and an order location table. Can "
        "optionally supply a master-flat and/or master-dark for UVB/VIS."
    )


# A VALID NIR SET: THE LAMP-ON MULTI-PINHOLE FRAME AND ITS LAMP-OFF PARTNER.
VALID_NIR_TECH = ["ECHELLE,MULTI-PINHOLE", "IMAGE"]
NIR_CATALOGUE = ["ORDER_TAB_NIR", "DISP_TAB_NIR"]
VIS_CATALOGUE = ["MASTER_BIAS_VIS", "ORDER_TAB_VIS", "DISP_TAB_VIS"]

# THE PACKAGE `__init__` BINDS THIS NAME TO THE RECIPE CLASS, SO THE MODULE
# ITSELF HAS TO BE IMPORTED BY PATH.
SPATIAL_MODULE = import_module("soxspipe.recipes.soxs_spatial_solution")


class FrameInventory:
    """Inventory exposing the `ImageFileCollection` seams the recipe uses."""

    def __init__(self, *, calls: list[str] | None = None) -> None:
        self.calls = calls
        self.sortedBy: list[str] | None = None
        self.summary = "SUMMARY TABLE"

    def sort(self, keywords: list[str]) -> None:
        self.sortedBy = list(keywords)
        if self.calls is not None:
            self.calls.append("sort")


def _recipe_with_inventory(log: Any) -> soxs_spatial_solution:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    recipe.log = log
    recipe.inputFrames = FrameInventory()
    recipe.kw = lambda keyword: keyword
    return recipe


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    *,
    imageTypes: list[str],
    imageTech: list[str],
    imageCat: list[str],
    arm: str,
    inst: str = "SOXS",
    calls: list[str] | None = None,
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real `_verify_input_frames_basics` sets `self.arm` and `self.inst`
    before it returns, and `verify_input_frames` reads both, so the stub sets
    them too.
    """

    def fake_basics(self: soxs_spatial_solution) -> tuple[list[str], list[str], list[str]]:
        if calls is not None:
            calls.append("verify")
        self.arm = arm
        self.inst = inst
        return list(imageTypes), list(imageTech), list(imageCat)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


def _isolate(
    monkeypatch: pytest.MonkeyPatch,
    tmpPath: Path,
    *,
    productPath: Path,
) -> None:
    """Replace every boundary the inherited constructor reaches outside its settings."""

    class IsolatedOrganiser:
        """Session lookup boundary with no external workspace state."""

        def __init__(self, **_: object) -> None:
            return None

        def session_list(self, *, silent: bool) -> tuple[str | bool, list[str]]:
            assert silent is True
            return False, []

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
    """Return the minimum settings the recipe constructor reads."""
    return {
        "workspace-root-dir": str(tmpPath),
        "instrument": "soxs",
        "data-extension": 0,
        "save-intermediate-products": False,
        **overrides,
    }


def _construct(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmpPath: Path,
    *,
    calls: list[str],
    verbose: bool = False,
    polyOrders: Any = False,
    create2DMap: bool = True,
    debug: bool = False,
    settings: dict[str, Any] | None = None,
) -> tuple[soxs_spatial_solution, dict[str, Any], FrameInventory]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, the keyword arguments `set_of_files` received,
      and the inventory the stubbed `set_of_files` handed back.
    """
    _isolate(monkeypatch, tmpPath, productPath=tmpPath / "reduced" / "spatial_map.fits")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["ECHELLE,MULTI-PINHOLE"],
        imageCat=VIS_CATALOGUE,
        arm="VIS",
        calls=calls,
    )

    inventory = FrameInventory(calls=calls)
    sofArguments: dict[str, Any] = {}
    preparedFrames = object()

    class StubSetOfFiles:
        """Set-of-files boundary returning the inventory under test."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.update(kwargs)
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, str]:
            return inventory, "SUPPLEMENTARY"

    def fake_prepare_frames(self: object, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        return preparedFrames

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)

    recipe = soxs_spatial_solution(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_spatial_solution.sof"),
        verbose=verbose,
        polyOrders=polyOrders,
        create2DMap=create2DMap,
        debug=debug,
        turnOffMP=True,
    )
    recipe.preparedFramesSentinel = preparedFrames
    return recipe, sofArguments, inventory


def _raised_message(recipe: soxs_spatial_solution) -> str:
    """Run the verification, which must reject, and return the rendered message."""
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()
    return str(raised.value)


# ---------------------------------------------------------------------------
# THE CONSTRUCTION SEAM
# ---------------------------------------------------------------------------


def test_the_constructor_establishes_the_attributes_the_reduction_reads(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A real construction leaves the prepared frames, the supplement and the image type."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert recipe.inputFrames is recipe.preparedFramesSentinel
    assert recipe.supplementaryInput == "SUPPLEMENTARY"
    assert recipe.imageType == "LAMP,WAVE"
    assert recipe.recipeName == "soxs-spat-solution"
    assert recipe.verbose is False
    assert recipe.polyOrders is False
    assert recipe.create2DMap is True
    assert recipe.debug is False
    assert recipe.settings["data-extension"] == 0


def test_the_constructor_keeps_the_map_and_debug_switches(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`create2DMap` and `debug` are stored as given, for the reduction to read."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, create2DMap=False, debug=True)

    # ASSERT
    assert recipe.create2DMap is False
    assert recipe.debug is True


def test_the_constructor_verifies_then_sorts_then_prepares(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The set of files is read, then verified, then sorted, then prepared."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _, _, inventory = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert calls == ["set_of_files", "verify", "sort", "prepare_frames(save=False)"]
    assert inventory.sortedBy == ["MJD-OBS"]


def test_the_set_of_files_receives_the_configured_data_extension(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`set_of_files` is built from the recipe's own log, settings, frames and extension."""
    # ARRANGE
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_spatial_solution.sof")

    # ACT
    recipe, sofArguments, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        settings=_settings(tmp_path, **{"data-extension": 1}),
    )

    # ASSERT
    assert sofArguments["ext"] == 1
    assert sofArguments["inputFrames"] == sofPath
    assert sofArguments["log"] is recipe.log
    assert sofArguments["settings"]["data-extension"] == 1


def test_saving_intermediate_products_reaches_the_frame_preparation(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The `save-intermediate-products` setting is the `save` argument, nothing else."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        settings=_settings(tmp_path, **{"save-intermediate-products": True}),
    )

    # ASSERT
    assert calls == ["set_of_files", "verify", "sort", "prepare_frames(save=True)"]


def test_a_verbose_construction_prints_the_raw_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Verbose mode prints the inventory summary between sorting and preparation."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls, verbose=True)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed[-4:] == [
        "# VERIFYING INPUT FRAMES - ALL GOOD",
        "# RAW INPUT FRAMES - SUMMARY",
        "SUMMARY TABLE",
        "\n",
    ]


def test_a_quiet_construction_prints_no_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Without verbose mode the last print is the verification all-clear."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed[-2:] == ["# VERIFYING INPUT FRAMES", "# VERIFYING INPUT FRAMES - ALL GOOD"]
    assert "SUMMARY TABLE" not in printed


# ---------------------------------------------------------------------------
# THE `polyOrders` ARGUMENT, WHICH ONLY THE CONSTRUCTOR VALIDATES
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("polyOrders", "expected"),
    [(345435, 345435), ("345435", 345435), (345435.9, 345435), (7, 7), (3454, 3454)],
)
def test_an_integer_like_poly_orders_is_kept_as_an_integer(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    polyOrders: Any,
    expected: int,
) -> None:
    """`int()` is the only check: floats truncate and the digit count is never tested."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=polyOrders)

    # ASSERT
    assert recipe.polyOrders == expected
    assert isinstance(recipe.polyOrders, int)


@pytest.mark.parametrize("polyOrders", ["not-a-number", [3, 4, 5, 4, 3, 5]])
def test_a_poly_orders_int_cannot_take_is_rejected_with_a_type_error(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    polyOrders: Any,
) -> None:
    """A `ValueError` or `TypeError` from `int()` is logged, then rejected."""
    # ARRANGE
    calls: list[str] = []

    # ACT / ASSERT
    with pytest.raises(TypeError, match="^THE poly VALUE NEEDS TO BE A 6 DIGIT INTEGER$"):
        _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=polyOrders)

    # THE FAILED COERCION IS SWALLOWED AND RECORDED BEFORE THE REJECTION, AND
    # THE REJECTION COMES BEFORE THE SET OF FILES IS READ.
    debugged = [message for level, message in log.messages if level == "debug"]
    assert any("`self.polyOrders = int(self.polyOrders)` failed" in message for message in debugged)
    assert calls == []


# ---------------------------------------------------------------------------
# `verify_input_frames`, NIR
# ---------------------------------------------------------------------------


def test_nir_rejects_an_unlisted_image_type_by_the_last_one_found(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Every type is checked, and the message names the last offender."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS", "LAMP,WAVE", "DARK"],
        imageTech=VALID_NIR_TECH,
        imageCat=NIR_CATALOGUE,
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == _nir_type_error("DARK")


def test_nir_rejects_an_unlisted_technique_by_the_last_one_found(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The technique check runs after the type check, over every technique."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["ECHELLE,MULTI-PINHOLE", "MOS", "IFU"],
        imageCat=NIR_CATALOGUE,
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == _nir_tech_error("IFU")


@pytest.mark.parametrize("imageTypes", [["LAMP,FLAT"], ["FLAT,LAMP", "LAMP,FLAT"]])
def test_nir_rejects_a_set_with_no_arc_lamp(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageTypes: list[str],
) -> None:
    """Flat lamps are allowed, but one of the two arc-lamp spellings must be present."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTech=VALID_NIR_TECH,
        imageCat=NIR_CATALOGUE,
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_ERROR


def test_nir_rejects_a_set_with_no_multi_pinhole_frame(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Listed techniques pass one by one, but the multi-pinhole technique is required."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["WAVE,LAMP"],
        imageTech=["ECHELLE,SLIT", "IMAGE", "ECHELLE,PINHOLE"],
        imageCat=NIR_CATALOGUE,
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_ERROR


@pytest.mark.parametrize("imageCat", [["ORDER_TAB_NIR"], ["DISP_TAB_NIR"], ["ORDER_TAB_VIS", "DISP_TAB_VIS"]])
def test_nir_rejects_a_missing_order_or_dispersion_table_for_its_own_arm(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageCat: list[str],
) -> None:
    """Both tables must be catalogued for the NIR arm."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=VALID_NIR_TECH,
        imageCat=imageCat,
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_ERROR


@pytest.mark.parametrize(
    ("imageTypes", "imageTech"),
    [
        (["LAMP,WAVE"], ["ECHELLE,MULTI-PINHOLE"]),
        (["FLAT,LAMP", "WAVE,LAMP"], ["ECHELLE,MULTI-PINHOLE", "IMAGE", "ECHELLE,SLIT", "ECHELLE,PINHOLE"]),
    ],
)
def test_nir_accepts_an_arc_set_and_records_its_first_image_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageTypes: list[str],
    imageTech: list[str],
) -> None:
    """A valid NIR set passes silently and records the first type, which may be a flat."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTech=imageTech,
        imageCat=NIR_CATALOGUE,
        arm="NIR",
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageTypes[0]


# ---------------------------------------------------------------------------
# `verify_input_frames`, UVB AND VIS
# ---------------------------------------------------------------------------


def test_uvb_vis_rejects_an_unlisted_image_type_by_the_last_one_found(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`FLAT,LAMP` is allowed on NIR only, and the message names the last offender."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS", "LAMP,WAVE", "FLAT,LAMP"],
        imageTech=["ECHELLE,MULTI-PINHOLE"],
        imageCat=VIS_CATALOGUE,
        arm="VIS",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == _uvb_vis_type_error("FLAT,LAMP")


@pytest.mark.parametrize(
    "imageCat",
    [
        ["ORDER_TAB_UVB", "DISP_TAB_UVB"],
        ["MASTER_BIAS_UVB", "DISP_TAB_UVB"],
        ["MASTER_BIAS_UVB", "ORDER_TAB_UVB"],
        ["MASTER_BIAS_VIS", "ORDER_TAB_VIS", "DISP_TAB_VIS"],
    ],
)
def test_uvb_vis_rejects_a_missing_bias_order_or_dispersion_table(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageCat: list[str],
) -> None:
    """All three calibrations must be catalogued for the arm under reduction."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["ECHELLE,MULTI-PINHOLE"],
        imageCat=imageCat,
        arm="UVB",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == UVB_VIS_ERROR


@pytest.mark.parametrize("imageTypes", [["LAMP,WAVE"], ["WAVE,LAMP"], ["LAMP,FLAT"], ["LAMP,FLAT", "LAMP,WAVE"]])
def test_uvb_vis_accepts_listed_types_without_requiring_an_arc_or_a_technique(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageTypes: list[str],
) -> None:
    """Unlike NIR, UVB/VIS checks neither for an arc lamp nor for any technique."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTech=["BIAS-LIKE-TECH"],
        imageCat=VIS_CATALOGUE,
        arm="VIS",
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageTypes[0]


def test_the_image_type_check_runs_before_the_catalogue_check(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With both faults present, the image-type message is the one raised."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(monkeypatch, imageTypes=["BIAS"], imageTech=[], imageCat=[], arm="VIS")

    # ACT / ASSERT
    assert _raised_message(recipe) == _uvb_vis_type_error("BIAS")


def test_a_rejection_prints_the_inventory_summary_and_leaves_the_image_type_unset(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The error path prints the summary, raises, and writes no `imageType`."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(monkeypatch, imageTypes=["BIAS"], imageTech=[], imageCat=VIS_CATALOGUE, arm="VIS")

    # ACT
    _raised_message(recipe)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# VERIFYING INPUT FRAMES - **ERROR**\n", "SUMMARY TABLE", ""]
    assert not hasattr(recipe, "imageType")


def test_verification_logs_its_own_entry_and_exit(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The passing path emits exactly the method's own debug pair."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(monkeypatch, imageTypes=["LAMP,WAVE"], imageTech=[], imageCat=VIS_CATALOGUE, arm="VIS")

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert debugged == [
        "starting the ``verify_input_frames`` method",
        "completed the ``verify_input_frames`` method",
    ]


# ---------------------------------------------------------------------------
# THE TUNING WORKER
# ---------------------------------------------------------------------------


def test_the_tuning_worker_writes_all_three_degree_pairs_then_fails_on_an_undefined_self(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """DY-61: the worker reads `self.debug` at module level and raises `NameError`.

    The three degree writes land first, and the dispersion map is never built.
    """
    # ARRANGE
    recipeSettings: dict[str, Any] = {}
    built: list[object] = []
    monkeypatch.setattr(commonutils, "create_dispersion_map", lambda **kwargs: built.append(kwargs))

    # ACT / ASSERT
    with pytest.raises(NameError, match="self"):
        SPATIAL_MODULE.parameterTuning(
            (2, 3, 4, 5, 1, 2),
            log=log,
            recipeSettings=recipeSettings,
            settings={},
            multiPinholeFrame=None,
            disp_map_table="DISP_TAB",
            order_table="ORDER_TAB",
            qc=None,
            products=None,
            sofName="sof",
            lineDetectionTable=None,
        )

    assert recipeSettings == {"order-deg": [2, 3], "wavelength-deg": [4, 5], "slit-deg": [1, 2]}
    assert built == []
