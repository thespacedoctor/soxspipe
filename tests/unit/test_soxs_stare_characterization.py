"""Characterization of `soxs_stare.__init__` and `soxs_stare.verify_input_frames`.

These tests pin what the recipe does today, before its constructor is split.
They are the recipe-layer constructor seam the tier 3 spec asks for: the
orchestration tests build a recipe with `__new__`, and the public-contract
test replaces `base_recipe.__init__` wholesale, so the stare constructor has
never run against the real inherited constructor in the suite.

The constructor tests call the real `soxs_stare.__init__`. Everything the
recipe inherits is replaced at its own boundary — the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation — so what runs for real is the recipe's own constructor
and its own `verify_input_frames`.

Stare differs from the first two tier-3 recipes in three ways the tests pin:
it builds `set_of_files` with no `ext` argument, it reads its recipe settings
before collecting the frames, and it names its products from the set-of-files
name after the frames are prepared. Without a set-of-files name, that last
step calls a name the module never imports (DY-111), a branch that DY-90
masks today.

`verify_input_frames` is tested by direct call. Its NIR branch and its
error-reporting block are the lines the suite does not reach, and every
rejection message is pinned exactly.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from astropy.nddata import CCDData

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_stare import soxs_stare
from tests.factories import synthetic_ccd

pytestmark = pytest.mark.unit

# THE STARE RECIPE SETTINGS THE CONSTRUCTOR READS BEFORE IT COLLECTS FRAMES.
STARE_SETTINGS = {"sky-subtraction": {"subtract_sky": True}, "use_flat": True}

# THE CALIBRATIONS A UVB OR VIS SET OF FRAMES MUST CARRY TO PASS VERIFICATION.
UVB_VIS_CATEGORIES = ["MASTER_BIAS_VIS", "DISP_TAB_VIS"]


def _uvb_vis_message(arm: str, missing: str) -> str:
    """Return the rejection message the UVB/VIS branch renders for one frame."""
    return (
        f"Input frames for soxspipe stare need to be an object frame (OBJECT_{arm}), "
        f"a dispersion map image (DISP_IMAGE_{arm}), a dispersion map table (DISP_TAB_{arm}), "
        f"an order-location table (ORDER_TAB_{arm}), a master-bias (MASTER_BIAS_{arm}), "
        f"a master-flat (MASTER_FLAT_{arm}) and optionally a master dark (MASTER_DARK_{arm}) "
        f"for UVB/VIS. The sof file is missing a {missing} frame."
    )


# THE PART OF EVERY NIR REJECTION MESSAGE THAT NAMES THE EXPECTED INPUT.
NIR_EXPECTED_INPUT = (
    "Input frames for soxspipe stare need to be an object frame (OBJECT_NIR), "
    "a dispersion map image (DISP_IMAGE_NIR), a dispersion map table (DISP_TAB_NIR), "
    "an order-location table (ORDER_TAB_NIR), a master-flat (MASTER_FLAT_NIR) "
    "and master dark (MASTER_DARK_NIR) or off-frame for NIR."
)


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


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    calls: list[str] | None = None,
    *,
    imageTypes: list[str] | None = None,
    imageTechniques: list[str] | None = None,
    imageCategories: list[str] | None = None,
    arm: str = "VIS",
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real basic verification is what records ``self.arm``, so the stub
    records it too. Stare's own verification branches on that attribute
    immediately afterwards.
    """
    resolvedTypes = ["OBJECT"] if imageTypes is None else imageTypes
    resolvedTechniques = ["ECHELLE,SLIT,STARE"] if imageTechniques is None else imageTechniques
    resolvedCategories = UVB_VIS_CATEGORIES if imageCategories is None else imageCategories

    def fake_basics(self: soxs_stare) -> tuple[list[str], list[str], list[str]]:
        if calls is not None:
            calls.append("verify")
        self.arm = arm
        return list(resolvedTypes), list(resolvedTechniques), list(resolvedCategories)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


def _unconstructed_recipe(log: Any, *, settings: dict[str, Any] | None = None) -> soxs_stare:
    """Return an unconstructed recipe carrying only what verification reads."""
    recipe = soxs_stare.__new__(soxs_stare)
    recipe.log = log
    recipe.arm = None
    recipe.settings = {} if settings is None else settings
    recipe.inputFrames = FrameInventory()
    recipe.kw = lambda keyword: keyword
    return recipe


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
        lambda *_: (str(tmpPath / "reduced" / "stare.fits"), "2024-01-02"),
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
        "soxs-stare": dict(STARE_SETTINGS),
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
    inputFrames: str | list[str] | None = None,
    preparedFrames: object | None = None,
) -> tuple[soxs_stare, dict[str, Any]]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, and the keyword arguments `set_of_files`
      received.
    """
    _isolate(monkeypatch, tmpPath)
    _stub_basics(monkeypatch, calls)
    inventory.calls = calls

    sofArguments: dict[str, Any] = {}
    preparedFrames = object() if preparedFrames is None else preparedFrames
    realGetRecipeSettings = base_recipe.get_recipe_settings

    class StubSetOfFiles:
        """Set-of-files boundary returning the inventory under test."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.update(kwargs)
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            return inventory, {}

    def fake_prepare_frames(self: soxs_stare, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        return preparedFrames

    def recording_get_recipe_settings(self: soxs_stare) -> Any:
        calls.append("get_recipe_settings")
        return realGetRecipeSettings(self)

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)
    monkeypatch.setattr(base_recipe, "get_recipe_settings", recording_get_recipe_settings)

    recipe = soxs_stare(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_stare.sof") if inputFrames is None else inputFrames,
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
    """A real construction leaves the prepared frames, the settings, the name template and the flags."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    recipe, _ = _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

    # ASSERT
    assert recipe.inputFrames is recipe.preparedFramesSentinel
    assert recipe.supplementaryInput == {}
    assert recipe.imageType == "OBJECT"
    assert recipe.recipeName == "soxs-stare"
    assert recipe.recipeSettings is recipe.settings["soxs-stare"]
    assert recipe.recipeSettings == STARE_SETTINGS
    assert recipe.sofName == "2024-01-02_stare"
    assert recipe.filenameTemplate == "2024-01-02_stare.fits"
    assert recipe.generateReponseCurve is False
    assert recipe.verbose is False


def test_the_constructor_reads_settings_then_collects_verifies_sorts_and_prepares(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The recipe settings are read before the set of files, and the frames are prepared last."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

    # ASSERT
    assert calls == [
        "get_recipe_settings",
        "set_of_files",
        "verify",
        "sort",
        "prepare_frames(save=False)",
    ]
    assert inventory.sortedBy == ["MJD-OBS"]


def test_the_set_of_files_receives_no_data_extension(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Unlike the other recipes, stare builds `set_of_files` without an `ext` argument."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_stare.sof")

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
    assert set(sofArguments) == {"log", "settings", "inputFrames"}
    assert sofArguments["inputFrames"] == sofPath
    assert sofArguments["log"] is recipe.log
    assert sofArguments["settings"] is recipe.settings


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
    assert calls[-1] == "prepare_frames(save=True)"


def test_a_verbose_construction_prints_the_raw_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Verbose mode prints the inventory summary after the verification announcements."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls, verbose=True)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed[-4:] == [
        "# VERIFYING INPUT FRAMES",
        "# VERIFYING INPUT FRAMES - ALL GOOD",
        "# RAW INPUT FRAMES - SUMMARY",
        "SUMMARY TABLE",
    ]


def test_a_quiet_construction_prints_only_the_verification_announcements(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """Without verbose mode the summary is not printed."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert printed[-2:] == ["# VERIFYING INPUT FRAMES", "# VERIFYING INPUT FRAMES - ALL GOOD"]
    assert "SUMMARY TABLE" not in printed


def test_a_list_of_frames_fails_in_the_inherited_constructor(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A list of frames never reaches the stare constructor's own body.

    This pins a defect, DY-90: `base_recipe.__init__` reads
    `self.startNightDate`, which only a set-of-files path sets. Whoever fixes
    DY-90 moves this pin, and then meets DY-111.
    """
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT / ASSERT
    with pytest.raises(AttributeError, match="startNightDate"):
        _construct(
            log,
            monkeypatch,
            tmp_path,
            inventory=inventory,
            calls=calls,
            inputFrames=[str(tmp_path / "object.fits")],
        )

    assert calls == []


class FramesWithFiles:
    """Prepared-frame collection stub exposing the `files_filtered` seam `filenamer` reads through."""

    def __init__(self, *, paths: list[str]) -> None:
        self._paths = paths

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return list(self._paths)


def test_a_missing_set_of_files_name_names_products_from_the_first_prepared_frame(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """With no set-of-files name, the product name falls back to the first prepared frame.

    This moves the DY-111 pin: `filenamer` is now imported, and the branch
    reads the first frame of the prepared collection (`self.inputFrames`,
    the only frame-bearing attribute the constructor has set by this point)
    rather than the never-set `self.objectFrame`. The branch is still masked
    by DY-90 today, so the inherited product-path step is stubbed to leave
    `sofName` false while still setting `startNightDate`.
    """
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []
    preparedFrames = FramesWithFiles(paths=[str(tmp_path / "prepared_frame.fits")])
    firstFrame = synthetic_ccd()

    def unnamed_product_path(self: soxs_stare, log: Any, *_: object) -> None:
        self.sofName = False
        self.productPath = False
        self.startNightDate = "2024-01-02"
        self.log = log

    def fake_ccddata_read(path: str, **_: object) -> CCDData:
        assert path == str(tmp_path / "prepared_frame.fits")
        return firstFrame

    monkeypatch.setattr(base_recipe, "_resolve_product_path", unnamed_product_path)
    monkeypatch.setattr(CCDData, "read", staticmethod(fake_ccddata_read))

    # ACT
    recipe, _ = _construct(
        log,
        monkeypatch,
        tmp_path,
        inventory=inventory,
        calls=calls,
        preparedFrames=preparedFrames,
    )

    # ASSERT
    assert recipe.filenameTemplate == "2024.01.02T03.04.05.678_VIS_RO2_BIAS.fits"
    assert calls[-1] == "prepare_frames(save=False)"


def test_a_vis_object_frame_with_its_calibrations_passes(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Verification is silent for a UVB/VIS object frame with a master bias and dispersion table."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageTypes=["OBJECT", "BIAS"])

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == "OBJECT"
    assert [message for level, message in log.messages if level == "print"] == []


def test_an_unexpected_vis_image_type_is_rejected_by_name(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A UVB/VIS frame of an unlisted type raises the UVB/VIS message naming that type."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageTypes=["OBJECT", "LAMP,WAVE"])

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _uvb_vis_message("VIS", "LAMP,WAVE")


def test_the_last_unexpected_image_type_is_the_one_reported(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Each unexpected type overwrites the message, so the last one in the list is named."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageTypes=["LAMP,WAVE", "LAMP,FMTCHK"], arm="UVB")

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _uvb_vis_message("UVB", "LAMP,FMTCHK")


@pytest.mark.parametrize(
    ("categories", "missing"),
    [
        (["DISP_TAB_VIS"], "MASTER_BIAS_VIS"),
        (["MASTER_BIAS_VIS"], "DISP_TAB_VIS"),
        ([], "DISP_TAB_VIS"),
    ],
)
def test_a_vis_set_without_a_master_bias_or_dispersion_table_is_rejected(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    categories: list[str],
    missing: str,
) -> None:
    """A missing calibration is named; when both are missing, the dispersion table is."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageCategories=categories)

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == _uvb_vis_message("VIS", missing)


def test_a_rejection_prints_the_error_banner_and_the_frame_summary(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Before raising, the recipe prints the error banner, the inventory summary and a blank line."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageCategories=[])

    # ACT / ASSERT
    with pytest.raises(TypeError):
        recipe.verify_input_frames()

    printed = [message for level, message in log.messages if level == "print"]
    assert printed == ["# VERIFYING INPUT FRAMES - **ERROR**\n", "SUMMARY TABLE", ""]


@pytest.mark.parametrize(
    ("imageTypes", "imageTechniques"),
    [
        (["OBJECT", "DARK"], ["ECHELLE,SLIT,STARE", "IMAGE"]),
        (["STD,FLUX"], ["ECHELLE,SLIT,NODDING"]),
        (["OBJECT,ASYNC", "LAMP,FLAT"], ["ECHELLE,SLIT", "ECHELLE,MULTI-PINHOLE"]),
        (["STD,TELLURIC"], ["ECHELLE,SLIT,STARE"]),
    ],
)
def test_a_nir_set_of_listed_types_and_techniques_passes(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    imageTypes: list[str],
    imageTechniques: list[str],
) -> None:
    """The NIR branch accepts every listed type and technique and needs no calibration categories."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(
        monkeypatch,
        imageTypes=imageTypes,
        imageTechniques=imageTechniques,
        imageCategories=[],
        arm="NIR",
    )

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert recipe.imageType == imageTypes[0]


def test_an_unexpected_nir_image_type_is_rejected_by_name(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A NIR frame of an unlisted type raises the NIR message prefixed with the type it found."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageTypes=["OBJECT", "BIAS"], arm="NIR")

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == f"Found a BIAS file. {NIR_EXPECTED_INPUT}"


def test_an_unexpected_nir_technique_is_rejected_by_name(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A NIR frame of an unlisted technique raises the NIR message naming the technique."""
    # ARRANGE
    recipe = _unconstructed_recipe(log)
    _stub_basics(monkeypatch, imageTechniques=["ECHELLE,SLIT,OFFSET"], arm="NIR")

    # ACT / ASSERT
    with pytest.raises(TypeError) as raised:
        recipe.verify_input_frames()

    assert str(raised.value) == (f"{NIR_EXPECTED_INPUT} The sof file is missing a ECHELLE,SLIT,OFFSET frame.")


@pytest.mark.parametrize(
    ("settings", "shouldPass"),
    [({"PAE": True}, True), ({"PAE": False}, False), ({}, False)],
)
def test_the_pae_setting_admits_the_nir_flat_lamp_pinhole_frame(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    settings: dict[str, Any],
    shouldPass: bool,
) -> None:
    """Only a true `PAE` setting adds `FLAT,LAMP` and `ECHELLE,PINHOLE` to the NIR lists."""
    # ARRANGE
    recipe = _unconstructed_recipe(log, settings=settings)
    _stub_basics(
        monkeypatch,
        imageTypes=["FLAT,LAMP"],
        imageTechniques=["ECHELLE,PINHOLE"],
        imageCategories=[],
        arm="NIR",
    )

    # ACT / ASSERT
    if shouldPass:
        recipe.verify_input_frames()
        assert recipe.imageType == "FLAT,LAMP"
    else:
        with pytest.raises(TypeError) as raised:
            recipe.verify_input_frames()
        assert str(raised.value) == f"Found a FLAT,LAMP file. {NIR_EXPECTED_INPUT}"
