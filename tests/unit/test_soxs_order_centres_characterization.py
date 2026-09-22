"""Characterization of `soxs_order_centres.__init__` and `soxs_order_centres.verify_input_frames`.

These tests pin what the recipe does today, before its constructor and its
frame verification are split. They are the recipe-layer construction seam the
tier 3 spec asks for: the orchestration tests build a recipe with `__new__`,
which skips `__init__` entirely, so the constructor this ticket edits is code
the offline suite otherwise never reaches.

The constructor tests call the real `soxs_order_centres.__init__`. Everything
the recipe inherits is replaced at its own boundary -- the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation -- so what runs for real is the recipe's own constructor
and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call. Each rejection branch's exact
rendered message is pinned here, and those pins are the evidence that the
`F507`/`UP031` rewrite of the `% locals()` sites renders a byte-identical
string. Two of the pinned messages are wrong today, and are pinned as they are.
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
from soxspipe.recipes.soxs_order_centres import soxs_order_centres

pytestmark = pytest.mark.unit

# THE FOUR RENDERED REJECTION MESSAGES. THE TWO NIR MESSAGES DIFFER BY ONE WORD,
# "and", AND THE FIRST UVB/VIS MESSAGE ENDS WITH A LITERAL, UNINTERPOLATED
# `{i}`, BECAUSE A `% locals()` SUFFIX DOES NOT FILL A BRACE PLACEHOLDER.
NIR_TYPE_ERROR = (
    "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off "
    "frames and a first-guess dispersion solution table for NIR"
)
NIR_ERROR = (
    "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off "
    "frames a first-guess dispersion solution table for NIR"
)
UVB_VIS_TYPE_ERROR = (
    "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias frame, "
    "a first-guess dispersion solution table and possibly a master dark for UVB/VIS. Found {i}"
)
UVB_VIS_ERROR = (
    "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias frame, "
    "a first-guess dispersion solution table and possibly a master dark for UVB/VIS."
)

# THE TECHNIQUES A VALID NIR SET CARRIES: THE PINHOLE FRAME AND THE LAMP-OFF FRAME.
VALID_NIR_TECH = ["ECHELLE,PINHOLE", "IMAGE"]

# THE PACKAGE `__init__` BINDS THIS NAME TO THE RECIPE CLASS, SO THE MODULE
# ITSELF HAS TO BE IMPORTED BY PATH.
ORDER_MODULE = import_module("soxspipe.recipes.soxs_order_centres")


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


def _recipe_with_inventory(log: Any) -> soxs_order_centres:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_order_centres.__new__(soxs_order_centres)
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

    def fake_basics(self: soxs_order_centres) -> tuple[list[str], list[str], list[str]]:
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
    settings: dict[str, Any] | None = None,
) -> tuple[soxs_order_centres, dict[str, Any], FrameInventory]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, the keyword arguments `set_of_files` received,
      and the inventory the stubbed `set_of_files` handed back.
    """
    _isolate(monkeypatch, tmpPath, productPath=tmpPath / "reduced" / "order_table.fits")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
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

    recipe = soxs_order_centres(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_order_centres.sof"),
        verbose=verbose,
        polyOrders=polyOrders,
        turnOffMP=True,
    )
    recipe.preparedFramesSentinel = preparedFrames
    return recipe, sofArguments, inventory


def _raised_message(recipe: soxs_order_centres) -> str:
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
    assert recipe.imageType == "LAMP,FLAT"
    assert recipe.recipeName == "soxs-order-centres"
    assert recipe.verbose is False
    assert recipe.polyOrders is False
    assert recipe.settings["data-extension"] == 0


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
    sofPath = str(tmp_path / "2024-01-02_order_centres.sof")

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
    [(34, 34), ("34", 34), (34.9, 34), (7, 7), (3454, 3454)],
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


@pytest.mark.parametrize("polyOrders", ["not-a-number", [3, 4]])
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
    with pytest.raises(TypeError, match="^THE poly VALUE NEEDS TO BE A 2 DIGIT INTEGER$"):
        _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=polyOrders)

    # THE FAILED COERCION IS SWALLOWED AND RECORDED BEFORE THE REJECTION, AND
    # THE REJECTION COMES BEFORE THE SET OF FILES IS READ.
    debugged = [message for level, message in log.messages if level == "debug"]
    assert any("`self.polyOrders = int(self.polyOrders)` failed" in message for message in debugged)
    assert calls == []


# ---------------------------------------------------------------------------
# `verify_input_frames`, NIR
# ---------------------------------------------------------------------------


def test_nir_mixed_image_types_raise_the_image_type_message_not_the_mix_message(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The mix message is bound to a misspelt name, so the next check rejects instead.

    The joined string's first character is then tested for the lamp type,
    which cannot match, so the set is still rejected -- with the wrong message.
    """
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["FLAT,LAMP", "LAMP,FLAT"],
        imageTech=VALID_NIR_TECH,
        imageCat=["DISP_TAB_NIR"],
        arm="NIR",
    )

    # ACT
    message = _raised_message(recipe)

    # ASSERT
    assert message == NIR_TYPE_ERROR
    assert "mix" not in message


@pytest.mark.parametrize(
    ("inst", "imageType"),
    [("SOXS", "BIAS"), ("SOXS", "LAMP,ORDERDEF"), ("XSH", "FLAT,LAMP"), ("XSH", "LAMP,QORDERDEF")],
)
def test_nir_rejects_an_image_type_without_the_instrument_lamp(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    inst: str,
    imageType: str,
) -> None:
    """The lamp test is a substring test on the one image type, per instrument."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=[imageType],
        imageTech=VALID_NIR_TECH,
        imageCat=["DISP_TAB_NIR"],
        arm="NIR",
        inst=inst,
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_TYPE_ERROR


def test_nir_rejects_an_unexpected_image_technique(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A technique outside the pinhole and image pair raises the NIR message."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["FLAT,LAMP"],
        imageTech=["ECHELLE,PINHOLE", "ECHELLE,SLIT"],
        imageCat=["DISP_TAB_NIR"],
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_ERROR


def test_nir_rejects_a_missing_dispersion_table(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The first-guess dispersion table must be catalogued for the NIR arm."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["FLAT,LAMP"],
        imageTech=VALID_NIR_TECH,
        imageCat=["DISP_TAB_VIS"],
        arm="NIR",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == NIR_ERROR


@pytest.mark.parametrize(
    ("inst", "imageType", "imageTech"),
    [
        ("SOXS", "FLAT,LAMP", ["ECHELLE,PINHOLE", "IMAGE"]),
        ("SOXS", "LAMP,DFLAT", ["ECHELLE,PINHOLE"]),
        ("XSH", "LAMP,ORDERDEF", ["ECHELLE,PINHOLE", "IMAGE"]),
    ],
)
def test_nir_accepts_a_lamp_set_with_its_dispersion_table_and_records_it(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    inst: str,
    imageType: str,
    imageTech: list[str],
) -> None:
    """A valid NIR set passes silently; the lamp-off frame is not required."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=[imageType],
        imageTech=imageTech,
        imageCat=["DISP_TAB_NIR"],
        arm="NIR",
        inst=inst,
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageType


# ---------------------------------------------------------------------------
# `verify_input_frames`, UVB AND VIS
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("inst", "imageTypes"),
    [
        ("SOXS", ["BIAS"]),
        ("SOXS", ["LAMP,FLAT", "LAMP,ORDERDEF"]),
        ("XSH", ["LAMP,FLAT"]),
    ],
)
def test_uvb_vis_rejects_any_image_type_outside_the_instrument_list(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    inst: str,
    imageTypes: list[str],
) -> None:
    """Every type is checked, and the message renders `{i}` literally."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
        arm="VIS",
        inst=inst,
    )

    # ACT
    message = _raised_message(recipe)

    # ASSERT
    assert message == UVB_VIS_TYPE_ERROR
    assert message.endswith("Found {i}")


@pytest.mark.parametrize(
    "imageCat",
    [["DISP_TAB_UVB"], ["MASTER_BIAS_UVB"], ["MASTER_BIAS_VIS", "DISP_TAB_VIS"]],
)
def test_uvb_vis_rejects_a_missing_bias_or_dispersion_table_for_its_own_arm(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageCat: list[str],
) -> None:
    """Both the master bias and the dispersion table must be catalogued for the arm."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=imageCat,
        arm="UVB",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == UVB_VIS_ERROR


@pytest.mark.parametrize(
    ("inst", "imageType"),
    [
        ("SOXS", "FLAT,LAMP"),
        ("SOXS", "LAMP,DFLAT"),
        ("SOXS", "LAMP,FLAT"),
        ("XSH", "LAMP,ORDERDEF"),
        ("XSH", "LAMP,DORDERDEF"),
        ("XSH", "LAMP,QORDERDEF"),
    ],
)
def test_uvb_vis_accepts_each_listed_lamp_type_and_records_it(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    inst: str,
    imageType: str,
) -> None:
    """A valid UVB/VIS set passes silently; the technique is not checked on this arm."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=[imageType],
        imageTech=["ECHELLE,SLIT", "IMAGE"],
        imageCat=["MASTER_BIAS_VIS", "DISP_TAB_VIS", "MASTER_DARK_VIS"],
        arm="VIS",
        inst=inst,
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageType


def test_the_image_type_check_runs_before_the_catalogue_check(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With both faults present, the image-type message is the one raised."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=[],
        arm="VIS",
    )

    # ACT / ASSERT
    assert _raised_message(recipe) == UVB_VIS_TYPE_ERROR


def test_a_rejection_prints_the_inventory_summary_and_leaves_the_image_type_unset(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The error path prints the summary, raises, and writes no `imageType`."""
    # ARRANGE
    recipe = _recipe_with_inventory(log)
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
        arm="VIS",
    )

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
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=["MASTER_BIAS_VIS", "DISP_TAB_VIS"],
        arm="VIS",
    )

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


def test_the_tuning_worker_writes_both_degrees_then_fails_on_an_undefined_self(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """DY-61: the worker reads `self.startNightDate` at module level and raises `NameError`.

    The two degree writes land first, and the detector is never built.
    """
    # ARRANGE
    recipeSettings: dict[str, Any] = {"detect-continuum": {}}
    built: list[object] = []
    monkeypatch.setattr(ORDER_MODULE, "detect_continuum", lambda **kwargs: built.append(kwargs))

    # ACT / ASSERT
    with pytest.raises(NameError, match="self"):
        ORDER_MODULE.parameterTuning(
            (3, 5),
            log=log,
            recipeSettings=recipeSettings,
            settings={},
            orderFrame=None,
            disp_map_table="DISP_TAB",
            orderPixelTable=None,
            qc=None,
            products=None,
            sofName="sof",
            binx=1,
            biny=1,
        )

    assert recipeSettings == {"detect-continuum": {"order-deg": 3, "disp-axis-deg": 5}}
    assert built == []
