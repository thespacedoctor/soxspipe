"""Characterization of `soxs_offset.__init__`.

These tests pin what the recipe's constructor does today, before it is split.
They are the recipe-layer constructor seam the tier 3 spec asks for: the
orchestration tests build a recipe with `__new__`, and the public-contract
test replaces `base_recipe.__init__` wholesale, so the offset constructor has
never run against the real inherited constructors in the suite.

`soxs_offset` subclasses `soxs_nod`, not `base_recipe` directly, and its
constructor calls `soxs_nod.__init__` first. That parent constructor already
collects, verifies, sorts and prepares the input frames. The offset
constructor then does all four steps again, from the original input, so a
real construction prepares the frames twice (DY-113). The tests pin both
passes and that the second one wins.

Everything the recipe inherits from `base_recipe` is replaced at its own
boundary. `soxs_nod.__init__` and `soxs_nod.verify_input_frames` run for real,
because they are the offset recipe's own parent code.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_offset import soxs_offset

pytestmark = pytest.mark.unit

# THE CALIBRATIONS AN OFFSET SET OF FRAMES MUST CARRY TO PASS VERIFICATION.
OFFSET_CATEGORIES = ["DISP_TAB_VIS", "ORDER_TAB_VIS", "DISP_IMAGE_VIS"]

# ONE PASS OF THE INPUT WORKFLOW, AS EACH CONSTRUCTOR RUNS IT.
ONE_PASS = ["set_of_files", "verify", "sort", "prepare_frames(save=False)"]


class FrameInventory:
    """Inventory exposing the `ImageFileCollection` seams the recipe uses."""

    def __init__(self, name: str, calls: list[str]) -> None:
        self.name = name
        self.calls = calls
        self.sortedBy: list[str] | None = None
        self.summary = f"SUMMARY TABLE {name}"

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
        lambda *_: (str(tmpPath / "reduced" / "offset.fits"), "2024-01-02"),
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
    records it too. The nod verification the offset recipe inherits reads
    that attribute immediately afterwards.
    """

    def fake_basics(self: soxs_offset) -> tuple[list[str], list[str], list[str]]:
        calls.append("verify")
        self.arm = "VIS"
        return ["OBJECT"], ["ECHELLE,SLIT,OFFSET"], list(OFFSET_CATEGORIES)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


def _settings(tmpPath: Path, **overrides: Any) -> dict[str, Any]:
    """Return the minimum settings the recipe constructors read."""
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
) -> tuple[soxs_offset, list[dict[str, Any]], list[FrameInventory], list[object]]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, the keyword arguments each `set_of_files`
      received, the inventory each one returned, and the object each frame
      preparation returned, all in call order.
    """
    _isolate(monkeypatch, tmpPath)
    _stub_basics(monkeypatch, calls)

    sofArguments: list[dict[str, Any]] = []
    inventories: list[FrameInventory] = []
    preparedFrames: list[object] = []

    class StubSetOfFiles:
        """Set-of-files boundary returning a fresh inventory per construction pass."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.append(dict(kwargs))
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            inventory = FrameInventory(f"pass {len(sofArguments)}", calls)
            inventories.append(inventory)
            return inventory, {}

    def fake_prepare_frames(self: soxs_offset, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        prepared = object()
        preparedFrames.append(prepared)
        return prepared

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)

    recipe = soxs_offset(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_offset.sof"),
        verbose=verbose,
        turnOffMP=True,
    )
    return recipe, sofArguments, inventories, preparedFrames


def test_the_constructor_establishes_the_attributes_the_reduction_reads(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A real construction leaves the second pass's prepared frames, the supplement and the image type."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    recipe, _, _, preparedFrames = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert len(preparedFrames) == 2
    assert recipe.inputFrames is preparedFrames[1]
    assert recipe.supplementaryInput == {}
    assert recipe.imageType == "OBJECT"
    assert recipe.recipeName == "soxs-offset"
    assert recipe.sofName == "2024-01-02_offset"
    assert recipe.verbose is False


def test_the_parent_constructor_runs_the_whole_input_workflow_before_the_offset_one(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The nod constructor collects, verifies, sorts and prepares, then the offset constructor repeats it."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _, _, inventories, _ = _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    assert calls == ONE_PASS + ONE_PASS
    assert [inventory.sortedBy for inventory in inventories] == [["MJD-OBS"], ["MJD-OBS"]]


def test_both_passes_collect_from_the_original_input_with_the_data_extension(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Each `set_of_files` receives the original input, the offset recipe name and the configured extension."""
    # ARRANGE
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_offset.sof")

    # ACT
    recipe, sofArguments, _, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        calls=calls,
        settings=_settings(tmp_path, **{"data-extension": 1}),
    )

    # ASSERT
    assert len(sofArguments) == 2
    for arguments in sofArguments:
        assert set(arguments) == {"log", "settings", "inputFrames", "recipeName", "ext"}
        assert arguments["inputFrames"] == sofPath
        assert arguments["recipeName"] == "soxs-offset"
        assert arguments["ext"] == 1
        assert arguments["log"] is recipe.log
        assert arguments["settings"] is recipe.settings


def test_saving_intermediate_products_reaches_both_frame_preparations(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The `save-intermediate-products` setting is the `save` argument of both passes."""
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
    assert [call for call in calls if call.startswith("prepare_frames")] == [
        "prepare_frames(save=True)",
        "prepare_frames(save=True)",
    ]


def test_a_verbose_construction_prints_each_pass_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Verbose mode prints the verification announcements and a summary block once per pass."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls, verbose=True)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    onePass = [
        "# VERIFYING INPUT FRAMES",
        "# VERIFYING INPUT FRAMES - ALL GOOD",
        "# RAW INPUT FRAMES - SUMMARY",
    ]
    assert printed == [
        *onePass,
        "SUMMARY TABLE pass 1",
        "\n",
        *onePass,
        "SUMMARY TABLE pass 2",
        "\n",
    ]


def test_a_quiet_construction_prints_only_the_verification_announcements(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Without verbose mode each pass prints its two announcements and no summary."""
    # ARRANGE
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, calls=calls)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# VERIFYING INPUT FRAMES", "# VERIFYING INPUT FRAMES - ALL GOOD"] * 2
