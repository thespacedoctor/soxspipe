"""Characterization of `soxs_nod.__init__` and `soxs_nod.verify_input_frames`.

These tests pin what the recipe does today, before its constructor is split.
`soxs_nod` is the parent `soxs_offset` inherits from, and the tier 3 spec's
seam is the same one already used for the other tier 3 recipes: the
orchestration tests build a recipe with `__new__`, which skips `__init__`
entirely, so the constructor every tier 3 ticket edits is the code the
offline suite never reaches.

The constructor tests call the real `soxs_nod.__init__`. Everything the
recipe inherits from `base_recipe` is replaced at its own boundary - the
session lookup, the product-path prediction, the recipe logger, the basic
frame verification and the frame preparation - so what runs for real is the
recipe's own constructor and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call on a recipe built with
`__new__`. It branches on `self.recipeName` (nod vs. offset), and every
rejection message is pinned exactly, including the two currently-uncovered
`pass` branches that let `STD,FLUX` and `STD,TELLURIC` frames skip the
technique check in each branch, and the ordering rule that a later check in
the same call is skipped once an earlier one has already set the error.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_nod import soxs_nod

pytestmark = pytest.mark.unit

# THE CALIBRATIONS A NOD SET OF FRAMES MUST CARRY TO PASS VERIFICATION.
NOD_CATEGORIES = ["DISP_TAB_VIS", "ORDER_TAB_VIS", "DISP_IMAGE_VIS"]


def _nod_type_message(arm: str, found: str) -> str:
    """Return the rejection message the type loop and the nod technique loop render."""
    return (
        f"Found a {found} file. Input frames for soxspipe nod need to be an object/std nodding frames, "
        f"a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table "
        f"(ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm})."
    )


def _offset_type_message(arm: str, found: str) -> str:
    """Return the rejection message the offset technique loop renders."""
    return (
        f"Found a {found} file. Input frames for soxspipe offset need to be an object/std offset frames, "
        f"a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table "
        f"(ORDER_TAB_{arm}) and a master-flat (MASTER_FLAT_{arm})."
    )


def _missing_category_message(arm: str, missing: str) -> str:
    """Return the rejection message the calibration-category loop renders."""
    return (
        f"Input frames for soxspipe nod need to be an object/std nodding frames, a dispersion map image "
        f"(DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), an order-location table (ORDER_TAB_{arm}) and a "
        f"master-flat (MASTER_FLAT_{arm}). The sof file is missing a {missing} frame."
    )


# ============================================================================
# PART A: the real `soxs_nod.__init__` seam.
# ============================================================================


class FrameInventory:
    """Inventory exposing the `ImageFileCollection` seams the recipe uses."""

    def __init__(self, calls: list[str]) -> None:
        self.calls = calls
        self.sortedBy: list[str] | None = None
        self.summary = "SUMMARY TABLE"

    def sort(self, keywords: list[str]) -> None:
        self.sortedBy = list(keywords)
        self.calls.append("sort")


def _isolate(monkeypatch: pytest.MonkeyPatch, tmpPath: Path) -> None:
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
        lambda *_: (str(tmpPath / "reduced" / "nod.fits"), "2024-01-02"),
    )
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(
        toolkit,
        "utility_setup",
        lambda **_: (str(tmpPath / "qc"), str(tmpPath / "products")),
    )


def _stub_basics(monkeypatch: pytest.MonkeyPatch, calls: list[str]) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real basic verification is what records ``self.arm``, so the stub
    records it too. `verify_input_frames` reads that attribute immediately
    afterwards.
    """

    def fake_basics(self: soxs_nod) -> tuple[list[str], list[str], list[str]]:
        calls.append("verify")
        self.arm = "VIS"
        return ["OBJECT"], ["ECHELLE,SLIT,NODDING"], list(NOD_CATEGORIES)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


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
    settings: dict[str, Any] | None = None,
    recipeName: str = "soxs-nod",
) -> tuple[soxs_nod, dict[str, Any], FrameInventory, object]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, the keyword arguments `set_of_files` received,
      the inventory it returned, and the object frame preparation returned.
    """
    _isolate(monkeypatch, tmpPath)
    _stub_basics(monkeypatch, calls)

    # WRAP THE LOGGER'S PRINT SO THE VERIFICATION ANNOUNCEMENTS AND THE
    # SUMMARY BLOCK LAND IN THE SAME ORDERED LIST AS THE OTHER STEPS.
    originalPrint = log.print

    def recording_print(message: object) -> None:
        calls.append(f"print:{message}")
        originalPrint(message)

    monkeypatch.setattr(log, "print", recording_print)

    sofArguments: dict[str, Any] = {}
    inventory = FrameInventory(calls)
    preparedFrames = object()

    class StubSetOfFiles:
        """Set-of-files boundary returning the inventory under test."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.update(kwargs)
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            return inventory, {}

    def fake_prepare_frames(self: soxs_nod, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        return preparedFrames

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)

    recipe = soxs_nod(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_nod.sof"),
        verbose=verbose,
        turnOffMP=True,
        recipeName=recipeName,
    )
    return recipe, sofArguments, inventory, preparedFrames


def test_the_constructor_establishes_the_attributes_the_reduction_reads(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A real construction leaves the prepared frames, the supplement, the image type and the flags."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _, preparedFrames = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert recipe.inputFrames is preparedFrames
    assert recipe.supplementaryInput == {}
    assert recipe.imageType == "OBJECT"
    assert recipe.recipeName == "soxs-nod"
    assert recipe.verbose is False
    assert recipe.log is log


def test_the_constructor_runs_the_input_workflow_in_order(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The constructor collects, announces, verifies, sorts, then prepares, in that interleaved order."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _, _, inventory, _ = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert calls == [
        "set_of_files",
        "print:# VERIFYING INPUT FRAMES",
        "verify",
        "print:# VERIFYING INPUT FRAMES - ALL GOOD",
        "sort",
        "prepare_frames(save=False)",
    ]
    assert inventory.sortedBy == ["MJD-OBS"]


def test_set_of_files_receives_the_original_input_with_the_data_extension(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`set_of_files` receives the original input, the recipe name and the configured extension."""
    # ARRANGE
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_nod.sof")

    # ACT
    recipe, sofArguments, _, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        settings=_settings(tmp_path, **{"data-extension": 1}),
    )

    # ASSERT
    assert set(sofArguments) == {"log", "settings", "inputFrames", "recipeName", "ext"}
    assert sofArguments["inputFrames"] == sofPath
    assert sofArguments["recipeName"] == "soxs-nod"
    assert sofArguments["ext"] == 1
    assert sofArguments["log"] is recipe.log
    assert sofArguments["settings"] is recipe.settings


def test_the_recipe_name_override_reaches_set_of_files_and_the_recipe(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A `recipeName` override lands on the recipe and reaches `set_of_files`."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, sofArguments, _, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        recipeName="soxs-nod-custom",
    )

    # ASSERT
    assert recipe.recipeName == "soxs-nod-custom"
    assert sofArguments["recipeName"] == "soxs-nod-custom"


@pytest.mark.parametrize("saveIntermediateProducts", [True, False])
def test_saving_intermediate_products_reaches_prepare_frames(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    saveIntermediateProducts: bool,
) -> None:
    """The `save-intermediate-products` setting is the `save` argument of `prepare_frames`."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        settings=_settings(tmp_path, **{"save-intermediate-products": saveIntermediateProducts}),
    )

    # ASSERT
    assert calls[-1] == f"prepare_frames(save={saveIntermediateProducts})"


def test_a_verbose_construction_prints_the_raw_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Verbose mode prints the inventory summary after the verification announcements."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls, verbose=True)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == [
        "# VERIFYING INPUT FRAMES",
        "# VERIFYING INPUT FRAMES - ALL GOOD",
        "# RAW INPUT FRAMES - SUMMARY",
        "SUMMARY TABLE",
        "\n",
    ]


def test_a_quiet_construction_prints_only_the_verification_announcements(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Without verbose mode the summary block is not printed."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# VERIFYING INPUT FRAMES", "# VERIFYING INPUT FRAMES - ALL GOOD"]


# ============================================================================
# PART B: `verify_input_frames`, called directly on a recipe built with `__new__`.
# ============================================================================


class VerifiedFrames:
    """Minimal inventory carrying only the summary a rejection prints."""

    def __init__(self) -> None:
        self.summary = "SUMMARY TABLE"


def _unconstructed_recipe(log: Any, *, recipeName: str = "soxs-nod") -> soxs_nod:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.arm = None
    recipe.recipeName = recipeName
    recipe.kw = lambda keyword: keyword
    recipe.inputFrames = VerifiedFrames()
    return recipe


def _stub_basics_direct(
    monkeypatch: pytest.MonkeyPatch,
    *,
    imageTypes: list[str],
    imageTech: list[str],
    imageCat: list[str],
    arm: str = "VIS",
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real basic verification is what records ``self.arm``, so the stub
    records it too.
    """

    def fake_basics(self: soxs_nod) -> tuple[list[str], list[str], list[str]]:
        self.arm = arm
        return list(imageTypes), list(imageTech), list(imageCat)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


def test_a_disallowed_image_type_is_rejected_by_name(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An image type outside the allowed list is named in the first loop's message."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,SLIT,NODDING"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _nod_type_message("VIS", "BIAS")


def test_a_disallowed_technique_is_rejected_in_the_nod_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A technique outside the nod branch's list is named, and `ECHELLE,SLIT,OFFSET` is disallowed for nod."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-nod")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["OBJECT"],
        imageTech=["ECHELLE,SLIT,OFFSET"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _nod_type_message("VIS", "ECHELLE,SLIT,OFFSET")


def test_a_disallowed_technique_is_rejected_in_the_offset_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A recipe name containing "offset" switches to the offset branch's technique list and message."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-offset")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["OBJECT"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _offset_type_message("VIS", "ECHELLE,PINHOLE")


def test_the_offset_technique_passes_only_in_the_offset_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`ECHELLE,SLIT,OFFSET` is listed for the offset branch, unlike the nod branch (pinned above)."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-offset")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["OBJECT"],
        imageTech=["ECHELLE,SLIT,OFFSET"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "OBJECT"


def test_std_flux_and_std_telluric_skip_the_technique_check_in_the_offset_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Line 148: `STD,FLUX` and `STD,TELLURIC` frames pass the offset branch regardless of technique."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-offset")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["STD,FLUX", "STD,TELLURIC"],
        imageTech=["SOMETHING,UNLISTED", "SOMETHING,ELSE"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "STD,FLUX"


def test_std_flux_and_std_telluric_skip_the_technique_check_in_the_nod_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Line 161: `STD,FLUX` and `STD,TELLURIC` frames pass the nod branch regardless of technique."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-nod")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["STD,FLUX", "STD,TELLURIC"],
        imageTech=["SOMETHING,UNLISTED", "SOMETHING,ELSE"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "STD,FLUX"


@pytest.mark.parametrize(
    ("categories", "missing", "recipeName"),
    [
        (["ORDER_TAB_VIS", "DISP_IMAGE_VIS"], "DISP_TAB_VIS", "soxs-nod"),
        (["DISP_TAB_VIS", "DISP_IMAGE_VIS"], "ORDER_TAB_VIS", "soxs-nod"),
        (["DISP_TAB_VIS", "ORDER_TAB_VIS"], "DISP_IMAGE_VIS", "soxs-nod"),
        ([], "DISP_IMAGE_VIS", "soxs-nod"),
        ([], "DISP_IMAGE_VIS", "soxs-offset"),
    ],
)
def test_a_missing_calibration_category_is_named_and_the_last_missing_one_wins(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    categories: list[str],
    missing: str,
    recipeName: str,
) -> None:
    """The category loop does not stop at the first miss; when several are missing, the last one wins.

    The loop checks `DISP_TAB_`, `ORDER_TAB_`, then `DISP_IMAGE_` in that
    order and overwrites the error each time, so when all three are missing
    `DISP_IMAGE_` is the one named, not `DISP_TAB_`.

    The last case is the offset run of DY-125: unlike the image-type and
    technique checks above it, this loop never switches its wording, so an
    offset reduction is told about nodding input. Pinned as found.
    """
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName=recipeName)
    # THE TECHNIQUE EACH BRANCH ACCEPTS, SO THAT ONLY THE CATEGORY CHECK CAN FAIL.
    technique = "ECHELLE,SLIT,OFFSET" if "offset" in recipeName else "ECHELLE,SLIT,NODDING"
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["OBJECT"],
        imageTech=[technique],
        imageCat=categories,
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _missing_category_message("VIS", missing)


def test_a_disallowed_image_type_wins_over_a_missing_category(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Once the type loop sets an error, later loops are skipped, so the type message wins."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,SLIT,NODDING"],
        imageCat=[],
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _nod_type_message("VIS", "BIAS")


def test_a_disallowed_image_type_wins_over_a_missing_category_in_the_offset_branch(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The type loop's error survives into the offset branch, whose technique loop is therefore skipped.

    The message still names the nod recipe for an offset run, which is the same family of wording
    defect as DY-125.
    """
    # ARRANGE
    recipe = _unconstructed_recipe(log, recipeName="soxs-offset")
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=[],
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _nod_type_message("VIS", "BIAS")


def test_a_rejection_prints_the_error_banner_and_the_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Before raising, the recipe prints the error banner, the inventory summary and a blank line, in order."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=["ECHELLE,SLIT,NODDING"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# VERIFYING INPUT FRAMES - **ERROR**\n", "SUMMARY TABLE", ""]


def test_the_passing_path_sets_the_image_type_and_prints_nothing(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A frame set that clears every check sets `imageType` from the first entry and prints nothing."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics_direct(
        monkeypatch,
        imageTypes=["OBJECT", "LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT,NODDING", "IMAGE"],
        imageCat=list(NOD_CATEGORIES),
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "OBJECT"
    assert [message for level, message in log.messages if level == "print"] == []


# THE CONSTRUCTOR-RULE METHODS BOTH RECIPES OF THE NOD/OFFSET CHAIN NOW DEFINE.
SHARED_CONSTRUCTOR_STEPS = (
    "_collect_input_frames",
    "_verify_and_announce_input_frames",
    "_sort_and_report_input_frames",
)


@pytest.mark.parametrize("methodName", SHARED_CONSTRUCTOR_STEPS)
def test_the_offset_override_of_each_constructor_step_still_matches_the_nod_one(methodName: str) -> None:
    """`soxs_offset` overrides all three steps, so a nod construction of an offset recipe runs the copies.

    `soxs_nod.__init__` calls these three methods by name. On a `soxs_offset` instance each call
    dispatches to the offset override instead, so the two bodies must stay the same until DY-113
    deletes the copies. The comparison strips docstrings, so only the code has to match.
    """
    # ARRANGE
    from soxspipe.recipes.soxs_offset import soxs_offset

    # ACT
    nodBody = _method_body_ast(soxs_nod, methodName)
    offsetBody = _method_body_ast(soxs_offset, methodName)

    # ASSERT
    assert methodName in vars(soxs_nod)
    assert methodName in vars(soxs_offset)
    assert nodBody == offsetBody


def _method_body_ast(recipeClass: type, methodName: str) -> str:
    """Return the dumped AST of a method's body, with its docstring removed."""
    import ast
    import inspect
    import textwrap

    tree = ast.parse(textwrap.dedent(inspect.getsource(vars(recipeClass)[methodName])))
    body = tree.body[0].body
    if isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant):
        body = body[1:]
    return "\n".join(ast.dump(statement) for statement in body)
