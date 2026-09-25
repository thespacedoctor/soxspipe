"""Characterization of `soxs_disp_solution.__init__` and `soxs_disp_solution.verify_input_frames`.

These tests pin what the recipe does today, before its constructor and its
frame verification are split. They are the recipe-layer construction seam the
tier 3 spec asks for: the orchestration tests build a recipe with `__new__`,
which skips `__init__` entirely, so the constructor this ticket edits is code
the offline suite otherwise never reaches.

The constructor tests call the real `soxs_disp_solution.__init__`. Everything
the recipe inherits is replaced at its own boundary -- the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation -- so what runs for real is the recipe's own constructor
and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call, because every one of its eight
rejection branches is a line the suite does not reach. Each branch's exact
rendered message is pinned here, and those pins are the evidence that the
`F507`/`UP031` rewrite of the `% locals()` sites renders a byte-identical
string.
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
from soxspipe.recipes.soxs_disp_solution import soxs_disp_solution

pytestmark = pytest.mark.unit

# THE TWO GENERIC REJECTION MESSAGES, ONE PER ARM BRANCH. NEITHER NAMES THE
# OFFENDING FRAMES ON ITS OWN; THE NIR MIXED-TYPE BRANCH IS THE ONE EXCEPTION,
# APPENDING THE OFFENDING TYPES TO THIS BASE MESSAGE (SEE
# test_nir_rejects_mixed_input_image_types).
NIR_ERROR = "Input frames for soxspipe disp_solution need to be single pinhole lamp on and lamp off frames for NIR"
UVB_VIS_ERROR = (
    "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
    "and a master-bias and possibly a master dark for UVB/VIS"
)

# THE CLASSIFICATIONS A VALID NIR AND A VALID UVB/VIS SET CARRY.
VALID_NIR_TECH = ["ECHELLE,PINHOLE", "IMAGE"]
VALID_UVB_VIS_TECH = ["ECHELLE,PINHOLE"]

# THE PACKAGE `__init__` BINDS THIS NAME TO THE RECIPE CLASS, SO THE MODULE
# ITSELF HAS TO BE IMPORTED BY PATH.
DISPERSION_MODULE = import_module("soxspipe.recipes.soxs_disp_solution")


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


def _recipe_with_inventory(log: Any, *, arm: str) -> soxs_disp_solution:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    recipe.log = log
    recipe.inputFrames = FrameInventory()
    recipe.kw = lambda keyword: keyword
    recipe.arm = arm
    return recipe


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    *,
    imageTypes: list[str],
    imageTech: list[str],
    imageCat: list[str],
    arm: str,
    calls: list[str] | None = None,
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real `_verify_input_frames_basics` sets `self.arm` before it returns,
    and `verify_input_frames` reads that attribute, so the stub sets it too.
    """

    def fake_basics(self: soxs_disp_solution) -> tuple[list[str], list[str], list[str]]:
        if calls is not None:
            calls.append("verify")
        self.arm = arm
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
) -> tuple[soxs_disp_solution, dict[str, Any], FrameInventory]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, the keyword arguments `set_of_files` received,
      and the inventory the stubbed `set_of_files` handed back.
    """
    _isolate(monkeypatch, tmpPath, productPath=tmpPath / "reduced" / "disp_map.fits")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_VIS"],
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

    recipe = soxs_disp_solution(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_disp_solution.sof"),
        verbose=verbose,
        polyOrders=polyOrders,
        turnOffMP=True,
    )
    recipe.preparedFramesSentinel = preparedFrames
    return recipe, sofArguments, inventory


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
    assert recipe.recipeName == "soxs-disp-solution"
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
    sofPath = str(tmp_path / "2024-01-02_disp_solution.sof")

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
    assert "# RAW INPUT FRAMES - SUMMARY" in printed
    assert "SUMMARY TABLE" in printed
    assert "# VERIFYING INPUT FRAMES - ALL GOOD" in printed


# ---------------------------------------------------------------------------
# THE `polyOrders` ARGUMENT, WHICH ONLY THE CONSTRUCTOR VALIDATES
# ---------------------------------------------------------------------------


def test_a_four_digit_integer_poly_orders_is_kept_as_an_integer(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """An integer `polyOrders` survives the constructor unchanged."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=3454)

    # ASSERT
    assert recipe.polyOrders == 3454


def test_a_digit_string_poly_orders_is_coerced_to_an_integer(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A string of digits is converted, so the reduction later re-renders it."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders="3454")

    # ASSERT
    assert recipe.polyOrders == 3454


def test_a_non_numeric_poly_orders_is_rejected_with_a_type_error(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A value `int()` cannot take is logged and then rejected."""
    # ARRANGE
    calls: list[str] = []

    # ACT / ASSERT
    with pytest.raises(TypeError, match="THE poly VALUE NEEDS TO BE A 4 DIGIT INTEGER"):
        _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders="not-a-number")

    # THE FAILED COERCION IS SWALLOWED AND RECORDED BEFORE THE REJECTION.
    debugged = [message for level, message in log.messages if level == "debug"]
    assert any("`self.polyOrders = int(self.polyOrders)` failed" in message for message in debugged)


def test_a_list_poly_orders_is_rejected_with_a_type_error(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A `TypeError` from `int()` takes the same swallow-then-reject path."""
    # ARRANGE
    calls: list[str] = []

    # ACT / ASSERT
    with pytest.raises(TypeError, match="THE poly VALUE NEEDS TO BE A 4 DIGIT INTEGER"):
        _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=[3, 4, 5, 4])


def test_a_float_poly_orders_is_truncated_rather_than_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`int()` accepts a float, so a non-integer value passes the four-digit check."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=3454.9)

    # ASSERT
    assert recipe.polyOrders == 3454


def test_poly_orders_is_not_checked_for_four_digits(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The check is `isinstance(int)`, not a digit count, whatever the message says."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _ = _construct(log, monkeypatch, tmp_path, calls=calls, polyOrders=7)

    # ASSERT
    assert recipe.polyOrders == 7


# ---------------------------------------------------------------------------
# `verify_input_frames`, NIR
# ---------------------------------------------------------------------------


def test_nir_rejects_mixed_input_image_types(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """More than one image type raises, and the message names both types."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="NIR")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE", "WAVE,LAMP"],
        imageTech=VALID_NIR_TECH,
        imageCat=["CALIB"],
        arm="NIR",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == f"{NIR_ERROR}. Found LAMP,WAVE and WAVE,LAMP"


def test_nir_rejects_an_image_type_outside_the_pinhole_set(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A single unexpected image type raises the NIR message."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="NIR")
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=VALID_NIR_TECH,
        imageCat=["CALIB"],
        arm="NIR",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == NIR_ERROR


def test_nir_rejects_an_unexpected_image_technique(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A technique outside the pinhole and image pair raises the NIR message."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="NIR")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["ECHELLE,PINHOLE", "IMAGE", "ECHELLE,SLIT"],
        imageCat=["CALIB"],
        arm="NIR",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == NIR_ERROR


def test_nir_rejects_a_missing_lamp_off_frame(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """NIR needs both the pinhole and the lamp-off `IMAGE` technique present."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="NIR")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["ECHELLE,PINHOLE"],
        imageCat=["CALIB"],
        arm="NIR",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == NIR_ERROR


@pytest.mark.parametrize("imageType", ["LAMP,FMTCHK", "LAMP,WAVE", "WAVE,LAMP"])
def test_nir_accepts_each_pinhole_image_type_and_records_it(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageType: str,
) -> None:
    """A valid NIR set passes silently and records the image type."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="NIR")
    _stub_basics(
        monkeypatch,
        imageTypes=[imageType],
        imageTech=VALID_NIR_TECH,
        imageCat=["CALIB"],
        arm="NIR",
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageType


# ---------------------------------------------------------------------------
# `verify_input_frames`, UVB and VIS
# ---------------------------------------------------------------------------


def test_uvb_vis_rejects_an_image_type_outside_the_pinhole_set(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An unexpected image type raises the UVB/VIS message."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="VIS")
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_VIS"],
        arm="VIS",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == UVB_VIS_ERROR


def test_uvb_vis_rejects_a_missing_pinhole_technique(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Without `ECHELLE,PINHOLE` the UVB/VIS message is raised."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="UVB")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=["IMAGE"],
        imageCat=["MASTER_BIAS_UVB"],
        arm="UVB",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == UVB_VIS_ERROR


def test_uvb_vis_rejects_a_missing_master_bias_for_its_own_arm(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The master bias must be catalogued for the arm under reduction."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="VIS")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_UVB"],
        arm="VIS",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == UVB_VIS_ERROR


def test_uvb_vis_rejects_every_frame_when_one_image_type_is_wrong(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The UVB/VIS type check runs over every type, not only the first."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="VIS")
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE", "BIAS"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_VIS"],
        arm="VIS",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == UVB_VIS_ERROR


@pytest.mark.parametrize("arm", ["UVB", "VIS"])
def test_uvb_vis_accepts_a_pinhole_set_with_its_master_bias(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    arm: str,
) -> None:
    """A valid UVB/VIS set passes silently and records the image type."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm=arm)
    _stub_basics(
        monkeypatch,
        imageTypes=["LAMP,WAVE"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=[f"MASTER_BIAS_{arm}", f"MASTER_DARK_{arm}"],
        arm=arm,
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "LAMP,WAVE"


def test_a_rejection_prints_the_inventory_summary_before_it_raises(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The error path prints the summary so the user sees the offending frames."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="VIS")
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_VIS"],
        arm="VIS",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    printed = [message for level, message in log.messages if level == "print"]
    assert "# VERIFYING INPUT FRAMES - **ERROR**\n" in printed
    assert "SUMMARY TABLE" in printed


def test_the_tuning_worker_rewrites_the_degrees_and_returns_nothing(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One tuning permutation splits into the two degree pairs and records nothing."""
    # ARRANGE
    recipeSettings: dict[str, Any] = {}
    mapArguments: dict[str, Any] = {}

    class RecordingDispersionMap:
        def __init__(self, **kwargs: object) -> None:
            mapArguments.update(kwargs)

        def get(self) -> tuple[object, ...]:
            return None, None, (), None, None, None

    monkeypatch.setattr(commonutils, "create_dispersion_map", RecordingDispersionMap)

    # ACT
    returned = DISPERSION_MODULE.parameterTuning(
        (3, 4, 5, 6),
        log=log,
        recipeSettings=recipeSettings,
        settings={},
        pinholeFrame=None,
        qc=None,
        products=None,
        sofName="sof",
        lineDetectionTable="LINE-DETECTION-TABLE",
    )

    # ASSERT
    assert returned is None
    assert recipeSettings == {"order-deg": [3, 4], "wavelength-deg": [5, 6]}
    assert mapArguments["create2DMap"] is False
    assert mapArguments["startNightDate"] is False
    assert mapArguments["lineDetectionTable"] == "LINE-DETECTION-TABLE"


def test_a_rejection_leaves_the_image_type_attribute_unset(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`self.imageType` is written only on the passing path."""
    # ARRANGE
    recipe = _recipe_with_inventory(log, arm="VIS")
    _stub_basics(
        monkeypatch,
        imageTypes=["BIAS"],
        imageTech=VALID_UVB_VIS_TECH,
        imageCat=["MASTER_BIAS_VIS"],
        arm="VIS",
    )

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    assert not hasattr(recipe, "imageType")
