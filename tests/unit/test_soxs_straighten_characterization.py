"""Characterization of `soxs_straighten.__init__` and `soxs_straighten.verify_input_frames`.

These tests pin what the recipe does today, before its constructor is split.
They are the recipe-layer constructor seam the tier 3 spec asks for: the
orchestration tests build a recipe with `__new__`, which skips `__init__`
entirely, so the constructor every tier 3 ticket edits is the code the offline
suite never reaches.

The constructor tests call the real `soxs_straighten.__init__`. Everything the
recipe inherits is replaced at its own boundary — the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation — so what runs for real is the recipe's own constructor
and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call, because its one live rejection
branch is the line the suite does not reach: an arm with no full
dispersion/spatial solution among the supplementary input. That branch's
message is also the module's one `UP031` site, so the exact rendered string is
pinned here before the percent format is rewritten.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_straighten import soxs_straighten

pytestmark = pytest.mark.unit

# THE ARM EVERY SYNTHETIC FRAME IN THIS MODULE IS TAKEN WITH. THE RECIPE
# REJECTS A SET OF FRAMES WHOSE SUPPLEMENTARY INPUT CARRIES NO FULL
# DISPERSION/SPATIAL SOLUTION FOR IT.
UNIFORM_ARM = "VIS"

# THE SUPPLEMENTARY INPUT A SET OF FRAMES MUST CARRY TO PASS VERIFICATION.
FULL_SOLUTION = {UNIFORM_ARM: {"2D_MAP": "vis_2d_map.fits"}}


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


def _recipe_with_supplement(
    log: Any,
    supplementaryInput: dict[str, Any],
    *,
    arm: str = UNIFORM_ARM,
) -> soxs_straighten:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_straighten.__new__(soxs_straighten)
    recipe.log = log
    recipe.arm = arm
    recipe.inputFrames = FrameInventory()
    recipe.supplementaryInput = supplementaryInput
    recipe.kw = lambda keyword: keyword
    return recipe


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    calls: list[str] | None = None,
    *,
    imageTypes: list[str] | None = None,
    arm: str = UNIFORM_ARM,
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real basic verification is what records ``self.arm``, so the stub
    records it too. Straighten's own verification reads that attribute
    immediately afterwards.
    """
    resolvedTypes = ["LAMP,FMTCHK"] if imageTypes is None else imageTypes

    def fake_basics(self: object) -> tuple[list[str], list[str], list[str]]:
        if calls is not None:
            calls.append("verify")
        self.arm = arm
        return list(resolvedTypes), ["IMAGE"], ["CALIB"]

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
    inventory: FrameInventory,
    calls: list[str],
    verbose: bool = False,
    settings: dict[str, Any] | None = None,
    supplementaryInput: dict[str, Any] | None = None,
) -> tuple[soxs_straighten, dict[str, Any]]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, and the keyword arguments `set_of_files`
      received.
    """
    _isolate(monkeypatch, tmpPath, productPath=tmpPath / "reduced" / "straighten.fits")
    _stub_basics(monkeypatch, calls)
    inventory.calls = calls
    supplement = FULL_SOLUTION if supplementaryInput is None else supplementaryInput

    sofArguments: dict[str, Any] = {}
    preparedFrames = object()

    class StubSetOfFiles:
        """Set-of-files boundary returning the inventory under test."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.update(kwargs)
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            return inventory, supplement

    def fake_prepare_frames(self: object, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        return preparedFrames

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)

    recipe = soxs_straighten(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_straighten.sof"),
        verbose=verbose,
        turnOffMP=True,
    )
    recipe.preparedFramesSentinel = preparedFrames
    return recipe, sofArguments


def test_the_constructor_establishes_the_attributes_the_reduction_reads(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A real construction leaves the prepared frames, the supplement and the image type."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    recipe, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        inventory=inventory,
        calls=calls,
    )

    # ASSERT
    assert recipe.inputFrames is recipe.preparedFramesSentinel
    assert recipe.supplementaryInput == FULL_SOLUTION
    assert recipe.imageType == "LAMP,FMTCHK"
    assert recipe.recipeName == "soxs-straighten"
    assert recipe.verbose is False
    assert recipe.settings["data-extension"] == 0


def test_the_constructor_verifies_then_sorts_then_prepares(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The set of files is read, then verified, then sorted, then prepared."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

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
    inventory = FrameInventory()
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_straighten.sof")

    # ACT
    recipe, sofArguments = _construct(
        log,
        monkeypatch,
        tmp_path,
        inventory=inventory,
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
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(
        log,
        monkeypatch,
        tmp_path,
        inventory=inventory,
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
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(
        log,
        monkeypatch,
        tmp_path,
        inventory=inventory,
        calls=calls,
        verbose=True,
    )

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert "# RAW INPUT FRAMES - SUMMARY" in printed
    assert "SUMMARY TABLE" in printed


def test_an_arm_absent_from_the_supplementary_input_is_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An arm with no supplementary entry at all raises `TypeError` naming that arm."""
    # ARRANGE
    recipe = _recipe_with_supplement(log, {})
    _stub_basics(monkeypatch)

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == ("Need a full dispersion/spatial solution for VIS - none found with the input files")


def test_a_supplementary_input_without_a_two_dimensional_map_is_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """An arm present but carrying no `2D_MAP` raises the same `TypeError`."""
    # ARRANGE
    recipe = _recipe_with_supplement(log, {UNIFORM_ARM: {"DISP_TAB": "vis_disp_tab.fits"}})
    _stub_basics(monkeypatch)

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == ("Need a full dispersion/spatial solution for VIS - none found with the input files")


def test_the_rejection_message_renders_whichever_arm_the_frames_carry(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The arm reaches the message verbatim, whichever arm the frames were taken with."""
    # ARRANGE
    recipe = _recipe_with_supplement(log, {}, arm="NIR")
    _stub_basics(monkeypatch, arm="NIR")

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == ("Need a full dispersion/spatial solution for NIR - none found with the input files")


def test_a_full_solution_passes_and_records_the_image_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verification is silent when the arm carries a full dispersion/spatial solution."""
    # ARRANGE
    recipe = _recipe_with_supplement(log, FULL_SOLUTION)
    _stub_basics(monkeypatch)

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "LAMP,FMTCHK"


def test_a_rejection_prints_nothing_before_it_raises(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The missing-solution branch raises directly, without the error summary block.

    Straighten's `error` flag is never assigned a message, so the block that
    prints the offending frames is unreachable. This pins that the user sees
    only the exception, which is the behaviour the refactor must preserve.
    """
    # ARRANGE
    recipe = _recipe_with_supplement(log, {})
    _stub_basics(monkeypatch)

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    printed = [message for level, message in log.messages if level == "print"]
    assert "# VERIFYING INPUT FRAMES - **ERROR**\n" not in printed
    assert "SUMMARY TABLE" not in printed
