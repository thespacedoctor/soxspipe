"""Characterization of `soxs_mdark.__init__` and `soxs_mdark.verify_input_frames`.

These tests pin what the recipe does today, before its constructor is split.
They are the recipe-layer constructor seam the tier 3 spec asks for: the
orchestration tests build a recipe with `__new__`, which skips `__init__`
entirely, so the constructor every tier 3 ticket edits is the code the offline
suite never reaches.

The constructor tests call the real `soxs_mdark.__init__`. Everything the
recipe inherits is replaced at its own boundary — the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation — so what runs for real is the recipe's own constructor
and its own `verify_input_frames`.

`verify_input_frames` is tested by direct call, because its three rejection
branches are the lines the suite does not reach: mixed image types, a
non-`DARK` image type, and differing exposure times.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_mdark import soxs_mdark

pytestmark = pytest.mark.unit

RECIPE_MODULE = importlib.import_module("soxspipe.recipes.soxs_mdark")

# THE EXPOSURE TIME EVERY SYNTHETIC DARK FRAME CARRIES. THE RECIPE REJECTS A
# SET OF FRAMES THAT DOES NOT AGREE ON IT.
UNIFORM_EXPTIME = 60.0


class FrameInventory:
    """Inventory exposing the `ImageFileCollection` seams the recipe uses."""

    def __init__(
        self,
        *,
        imageTypes: list[str] | None = None,
        exptimes: list[float] | None = None,
    ) -> None:
        self.imageTypes = ["DARK"] if imageTypes is None else imageTypes
        self.exptimes = [UNIFORM_EXPTIME] if exptimes is None else exptimes
        self.sortedBy: list[str] | None = None
        self.summary = "SUMMARY TABLE"

    def values(self, keyword: str, *, unique: bool) -> list[Any]:
        assert unique is True
        assert keyword == "EXPTIME"
        return list(self.exptimes)

    def sort(self, keywords: list[str]) -> None:
        self.sortedBy = list(keywords)


def _recipe_with_inventory(log: Any, inventory: FrameInventory) -> soxs_mdark:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_mdark.__new__(soxs_mdark)
    recipe.log = log
    recipe.inputFrames = inventory
    recipe.kw = lambda keyword: keyword
    return recipe


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    inventory: FrameInventory,
) -> None:
    """Replace the inherited basic verification with the classifications it returns."""
    monkeypatch.setattr(
        base_recipe,
        "_verify_input_frames_basics",
        lambda self: (list(inventory.imageTypes), ["IMAGE"], ["CALIB"]),
    )


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
) -> tuple[soxs_mdark, dict[str, Any]]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, and the keyword arguments `set_of_files`
      received.
    """
    _isolate(monkeypatch, tmpPath, productPath=tmpPath / "reduced" / "mdark.fits")
    _stub_basics(monkeypatch, inventory)

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

    recipe = soxs_mdark(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_dark.sof"),
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
    assert recipe.supplementaryInput == "SUPPLEMENTARY"
    assert recipe.imageType == "DARK"
    assert recipe.recipeName == "soxs-mdark"
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
    assert calls == ["set_of_files", "prepare_frames(save=False)"]
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
    sofPath = str(tmp_path / "2024-01-02_dark.sof")

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
    assert calls == ["set_of_files", "prepare_frames(save=True)"]


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


def test_mixed_input_image_types_are_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two image types raise `TypeError` naming both, joined by "and"."""
    # ARRANGE
    inventory = FrameInventory(imageTypes=["DARK", "BIAS"])
    recipe = _recipe_with_inventory(log, inventory)
    _stub_basics(monkeypatch, inventory)

    # ACT / ASSERT
    with pytest.raises(TypeError, match="Input frames are a mix of DARK and BIAS"):
        recipe.verify_input_frames()


def test_a_non_dark_image_type_is_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A single non-`DARK` image type raises `TypeError`."""
    # ARRANGE
    inventory = FrameInventory(imageTypes=["BIAS"])
    recipe = _recipe_with_inventory(log, inventory)
    _stub_basics(monkeypatch, inventory)

    # ACT / ASSERT
    with pytest.raises(TypeError, match="Input frames not DARK frames"):
        recipe.verify_input_frames()


def test_differing_exposure_times_are_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Two exposure times raise `TypeError` naming both, rendered as strings."""
    # ARRANGE
    inventory = FrameInventory(exptimes=[60.0, 120.0])
    recipe = _recipe_with_inventory(log, inventory)
    _stub_basics(monkeypatch, inventory)

    # ACT / ASSERT
    with pytest.raises(
        TypeError,
        match="Input frames have differing exposure-times 60.0 and 120.0",
    ):
        recipe.verify_input_frames()


def test_a_uniform_set_of_darks_passes_and_records_the_image_type(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verification is silent for one image type and one exposure time."""
    # ARRANGE
    inventory = FrameInventory()
    recipe = _recipe_with_inventory(log, inventory)
    _stub_basics(monkeypatch, inventory)

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "DARK"


def test_a_rejection_prints_the_inventory_summary_before_it_raises(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The error path prints the summary so the user sees the offending frames."""
    # ARRANGE
    inventory = FrameInventory(imageTypes=["BIAS"])
    recipe = _recipe_with_inventory(log, inventory)
    _stub_basics(monkeypatch, inventory)

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    printed = [message for level, message in log.messages if level == "print"]
    assert "# VERIFYING INPUT FRAMES - **ERROR**\n" in printed
    assert "SUMMARY TABLE" in printed
