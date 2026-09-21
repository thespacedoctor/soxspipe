"""Characterization of `soxs_mflat.__init__` and `soxs_mflat.produce_product` row contracts.

These tests pin what the recipe does today, before commit 2 replaces its
inline timestamp and QC/product row construction with `toolkit.utcnow_string`,
`toolkit.append_qc` and `toolkit.append_product`. They are the recipe-layer
constructor seam the tier 3 spec asks for, and the product-row contracts that
the offline orchestration suite exercises through `_assert_product_row` but
never checks for column order.

The constructor tests call the real `soxs_mflat.__init__`. Everything the
recipe inherits is replaced at its own boundary -- the session lookup, the
product-path prediction, the recipe logger, the basic frame verification and
the frame preparation -- so what runs for real is the recipe's own
constructor and its own `verify_input_frames`, against a synthetic VIS frame
inventory (`LAMP,FLAT` plus `MASTER_BIAS_VIS` and `ORDER_TAB_VIS`, one
exposure time) built to pass verification.

The `produce_product` tests build the recipe with `__new__`, exactly as
`tests/integration/test_recipe_orchestration.py` does, and start `self.qc`
and `self.products` from the real declared tables
(`base_recipe._empty_qc_and_product_tables`), not the richer factory tables
those orchestration tests use. That choice matters: the factory tables
already carry a `product_desc` column in the middle of the row, so appending
a `product_desc` key never moves it. The real declared table has no such
column, so the first row that carries one pushes it to the end of the table
-- pinned here as the "real" column order the next commit must preserve.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.set_of_files as set_of_files_module
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_mflat import soxs_mflat
from tests.factories import pipeline_settings, prepared_fits, synthetic_ccd

# IMPORTED VIA `import_module`, NOT `import soxspipe.recipes.soxs_mflat`: THE
# PACKAGE `__init__` REBINDS THE ATTRIBUTE `soxspipe.recipes.soxs_mflat` TO
# THE CLASS ITSELF (`from .soxs_mflat import soxs_mflat`), SO A DOTTED
# `import ... as` WOULD RESOLVE TO THE CLASS, NOT THE MODULE, AND
# `monkeypatch.setattr` ON MODULE-LEVEL COLLABORATORS WOULD FAIL.
soxs_mflat_module = import_module("soxspipe.recipes.soxs_mflat")

pytestmark = pytest.mark.unit

# THE FORMAT EVERY `datetime.utcnow().strftime(...)` CALL RENDERS. FREEZING
# THE CONTRACT, NOT THE CLOCK: THIS HOLDS BEFORE AND AFTER COMMIT 2 SWAPS IN
# `toolkit.utcnow_string`.
TIMESTAMP_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}")

# THE EXPOSURE TIME EVERY SYNTHETIC LAMP-FLAT FRAME CARRIES. VERIFICATION
# REJECTS A SET OF LAMP FRAMES THAT DOES NOT AGREE ON IT.
UNIFORM_EXPTIME = 30.0


def _empty_products() -> pd.DataFrame:
    """Return the products table exactly as `base_recipe` declares it, empty."""
    return base_recipe._empty_qc_and_product_tables(None, pd)[1]


def _empty_qc() -> pd.DataFrame:
    """Return the QC table exactly as `base_recipe` declares it, empty."""
    return base_recipe._empty_qc_and_product_tables(None, pd)[0]


class FrameInventory:
    """Inventory exposing the `ImageFileCollection` seams the recipe uses.

    `filter` matches by the *value* a caller passes, not by the resolved FITS
    keyword name, so the double stands in for the real `ImageFileCollection`
    without needing to know what `keyword_lookup` resolves `DPR_TYPE` and
    `DPR_TECH` to for the "soxs" instrument under test.
    """

    def __init__(
        self,
        *,
        exptimesByType: dict[str, list[float]] | None = None,
        files: list[str] | None = None,
        exptimes: list[float] | None = None,
        calls: list[str] | None = None,
    ) -> None:
        self.calls = calls
        self.files = [] if files is None else files
        self._exptimes = [] if exptimes is None else exptimes
        self.exptimesByType = {"LAMP,FLAT": [UNIFORM_EXPTIME]} if exptimesByType is None else exptimesByType
        self.sortedBy: list[str] | None = None
        self.summary = "SUMMARY TABLE"

    def filter(self, **filters: str) -> FrameInventory:
        for lampType, exptimes in self.exptimesByType.items():
            if lampType in filters.values():
                return FrameInventory(files=[f"{lampType}.fits"], exptimes=exptimes)
        return FrameInventory(files=[])

    def values(self, keyword: str, unique: bool) -> list[float]:
        assert unique is True
        return list(dict.fromkeys(self._exptimes))

    def sort(self, keywords: list[str]) -> None:
        self.sortedBy = list(keywords)
        if self.calls is not None:
            self.calls.append("sort")


def _stub_basics(
    monkeypatch: pytest.MonkeyPatch,
    calls: list[str] | None = None,
    *,
    imageTypes: list[str] | None = None,
    imageCategories: list[str] | None = None,
    arm: str = "VIS",
    inst: str = "SOXS",
) -> None:
    """Replace the inherited basic verification with the classifications it returns.

    The real basic verification is what records `self.arm` and `self.inst`,
    so the stub records both -- mflat's own `verify_input_frames` branches on
    `self.arm` immediately and reads `self.inst` for the UVB/VIS D/Q-lamp
    check.
    """
    resolvedTypes = ["LAMP,FLAT"] if imageTypes is None else imageTypes
    resolvedCategories = ["MASTER_BIAS_VIS", "ORDER_TAB_VIS"] if imageCategories is None else imageCategories

    def fake_basics(self: soxs_mflat) -> tuple[list[str], list[str], list[str]]:
        if calls is not None:
            calls.append("verify")
        self.arm = arm
        self.inst = inst
        return list(resolvedTypes), ["ECHELLE,SLIT"], list(resolvedCategories)

    monkeypatch.setattr(base_recipe, "_verify_input_frames_basics", fake_basics)


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
        lambda *_: (str(tmpPath / "reduced" / "mflat.fits"), "2024-01-02"),
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
) -> tuple[soxs_mflat, dict[str, Any]]:
    """Construct a recipe for real against a stubbed set of files.

    **Return:**

    - the constructed recipe, and the keyword arguments `set_of_files`
      received.
    """
    _isolate(monkeypatch, tmpPath)
    _stub_basics(monkeypatch, calls)
    inventory.calls = calls

    sofArguments: dict[str, Any] = {}
    preparedFrames = object()

    class StubSetOfFiles:
        """Set-of-files boundary returning the inventory under test."""

        def __init__(self, **kwargs: object) -> None:
            sofArguments.update(kwargs)
            calls.append("set_of_files")

        def get(self) -> tuple[FrameInventory, dict[str, Any]]:
            return inventory, {}

    def fake_prepare_frames(self: soxs_mflat, save: bool = False) -> object:
        calls.append(f"prepare_frames(save={save})")
        return preparedFrames

    monkeypatch.setattr(set_of_files_module, "set_of_files", StubSetOfFiles)
    monkeypatch.setattr(base_recipe, "prepare_frames", fake_prepare_frames)

    recipe = soxs_mflat(
        log=log,
        settings=_settings(tmpPath) if settings is None else settings,
        inputFrames=str(tmpPath / "2024-01-02_mflat.sof"),
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
    """A real construction leaves the prepared frames, the supplement, and the image type."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    # ACT
    recipe, _ = _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

    # ASSERT
    assert recipe.inputFrames is recipe.preparedFramesSentinel
    assert recipe.supplementaryInput == {}
    assert recipe.imageType == "LAMP,FLAT"
    assert recipe.recipeName == "soxs-mflat"
    assert recipe.verbose is False
    assert recipe.settings["data-extension"] == 0


def test_the_constructor_collects_verifies_sorts_and_prepares_in_order(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """The frames are collected, verified, sorted, and only then prepared.

    The verification announcements are interleaved with those steps, so the
    recorded log printing itself into the same list pins the true order:
    collect, print, verify, print, sort, prepare.
    """
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []

    originalPrint = log.print

    def recording_print(message: object) -> None:
        calls.append(f"print:{message}")
        return originalPrint(message)

    monkeypatch.setattr(log, "print", recording_print)

    # ACT
    _construct(log, monkeypatch, tmp_path, inventory=inventory, calls=calls)

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


def test_the_set_of_files_receives_the_configured_data_extension(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`set_of_files` is built from the recipe's own log, settings, frames and extension."""
    # ARRANGE
    inventory = FrameInventory()
    calls: list[str] = []
    sofPath = str(tmp_path / "2024-01-02_mflat.sof")

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
    assert printed[-5:] == [
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


class _OrderTableOnlyCollection:
    """Image collection stub serving only the `ORDER_TAB_VIS` lookup `produce_product` needs."""

    def __init__(self, orderTablePath: str) -> None:
        self._orderTablePath = orderTablePath

    def filter(self, **_: object) -> _OrderTableOnlyCollection:
        return self

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return [self._orderTablePath]


class _RecordingSubtractBackground:
    """`subtract_background` collaborator double recording the arguments it received."""

    instances: list[_RecordingSubtractBackground] = []

    def __init__(self, **kwargs: object) -> None:
        self.kwargs = kwargs
        type(self).instances.append(self)

    def subtract(self) -> tuple[Any, Any, pd.DataFrame]:
        raise NotImplementedError


_PRODUCT_COLUMNS_WITH_DESC = [
    "soxspipe_recipe",
    "product_label",
    "file_name",
    "file_type",
    "obs_date_utc",
    "reduction_date_utc",
    "file_path",
    "label",
    "product_desc",
]


def _new_bare_recipe(
    log: Any,
    tmp_path: Path,
    orderPath: Path,
    *,
    subtractBackground: bool,
    calibratedFlatFiles: list[str],
    dFlatFiles: list[str],
    qFlatFiles: list[str],
    domeFlatFiles: list[str],
) -> soxs_mflat:
    """Return a `produce_product`-ready recipe carrying only what that method reads."""
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.kw = lambda keyword: "DATE-OBS" if keyword == "DATE_OBS" else keyword
    recipe.inputFrames = _OrderTableOnlyCollection(str(orderPath))
    recipe.recipeName = "soxs-mflat"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.recipeSettings = {"subtract_background": subtractBackground}
    recipe.sofName = "MASTER_FLAT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.binRatioX = 1
    recipe.binRatioY = 1
    recipe.calibratedFlatFiles = calibratedFlatFiles
    recipe.dFlatFiles = dFlatFiles
    recipe.qFlatFiles = qFlatFiles
    recipe.domeFlatFiles = domeFlatFiles
    recipe.qc = _empty_qc()
    recipe.products = _empty_products()
    return recipe


def _make_fake_edges(recipe: soxs_mflat, orderPath: Path) -> type:
    """Return a `detect_order_edges` double appending one tag-scoped ORDER_LOC row.

    Mirrors the real collaborator's own `pd.concat` call, `product_desc`
    included, so the column-order pin in the tests reflects what
    `detect_order_edges` really does to `self.products`, not a
    simplification of it.
    """

    class FakeEdges:
        def __init__(self, **kwargs: object) -> None:
            self._products = kwargs["productsTable"]
            self._qc = kwargs["qcTable"]
            self._tag = kwargs["tag"]

        def get(self) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
            row = {
                "soxspipe_recipe": recipe.recipeName,
                "product_label": f"ORDER_LOC{self._tag}",
                "product_desc": "table of coefficients from polynomial fits to order locations",
                "file_name": orderPath.name,
                "file_type": "FITS",
                "obs_date_utc": "2024-01-02T03:04:05.678",
                "reduction_date_utc": "2024-01-02T04:05:06",
                "file_path": str(orderPath),
                "label": "PROD",
            }
            products = pd.concat([self._products, pd.DataFrame([row])], ignore_index=True)
            return products, self._qc.copy(deep=True), pd.DataFrame({"order": [10], "count": [5]})

    return FakeEdges


def _stub_shared_collaborators(
    recipe: soxs_mflat,
    monkeypatch: pytest.MonkeyPatch,
    orderPath: Path,
) -> None:
    """Patch the module-level and bookkeeping collaborators every scenario below needs."""
    monkeypatch.setattr(recipe, "update_fits_keywords", lambda **k: None)
    monkeypatch.setattr(recipe, "report_output", lambda: recipe.qc)
    monkeypatch.setattr(recipe, "clean_up", lambda: None)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **k: None)
    monkeypatch.setattr(
        soxs_mflat_module,
        "unpack_order_table",
        lambda **k: (pd.DataFrame(), pd.DataFrame(), pd.DataFrame()),
    )
    monkeypatch.setattr(soxs_mflat_module, "generic_quality_checks", lambda **k: k["qcTable"])
    monkeypatch.setattr(soxs_mflat_module, "spectroscopic_image_quality_checks", lambda **k: k["qcTable"])
    monkeypatch.setattr(soxs_mflat_module, "detect_order_edges", _make_fake_edges(recipe, orderPath))


def _sequential_write_stub(
    paths: list[str],
) -> tuple[Any, list[tuple[tuple[object, ...], dict[str, object]]]]:
    """Return a `_write` double returning one path per call, and the calls it recorded."""
    calls: list[tuple[tuple[object, ...], dict[str, object]]] = []

    def fake_write(*args: object, **kwargs: object) -> str:
        calls.append((args, kwargs))
        return paths[len(calls) - 1]

    return fake_write, calls


def _assert_background_subtraction_products(
    recipe: soxs_mflat,
    *,
    backgroundProductPath: Path,
    tagProductPath: Path,
    finalProductPath: Path,
) -> None:
    """Pin the BKGROUND, tag-scoped MFLAT and final MFLAT rows and the table's column order."""
    assert list(recipe.products.columns) == _PRODUCT_COLUMNS_WITH_DESC
    assert list(recipe.products["product_label"]) == ["ORDER_LOC", "BKGROUND", "MFLAT", "MFLAT"]

    bkgRow = recipe.products.iloc[1]
    assert bkgRow["soxspipe_recipe"] == "soxs-mflat"
    assert bkgRow["file_name"] == "MASTER_FLAT_VIS_BKGROUND.fits"
    assert bkgRow["file_type"] == "FITS"
    assert bkgRow["obs_date_utc"] == recipe.dateObs
    assert TIMESTAMP_PATTERN.fullmatch(bkgRow["reduction_date_utc"])
    assert bkgRow["product_desc"] == "modelled scatter background light image (removed from master flat)"
    assert bkgRow["file_path"] == str(backgroundProductPath)
    assert bkgRow["label"] == "QC"

    tagRow = recipe.products.iloc[2]
    assert tagRow["file_name"] == "MASTER_FLAT_VIS.fits"
    assert tagRow["product_desc"] == "VIS master spectroscopic flat frame"
    assert tagRow["file_path"] == str(tagProductPath)
    assert tagRow["label"] == "PROD"

    finalRow = recipe.products.iloc[3]
    assert finalRow["file_name"] == "MASTER_FLAT_VIS_final.fits"
    assert finalRow["product_desc"] == "VIS master spectroscopic flat frame"
    assert finalRow["file_path"] == str(finalProductPath)
    assert finalRow["label"] == "PROD"


def test_single_lamp_background_subtraction_records_bkground_and_mflat_rows(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A background-subtracting single-lamp reduction records BKGROUND, then two MFLAT rows.

    `product_desc` is not a column `base_recipe` declares, so the first row
    that carries one -- here the BKGROUND row -- is what pushes it to the end
    of `recipe.products`, after `label`. The tag-scoped MFLAT row (written
    from `self.sofName + ".fits"`) and the final MFLAT row (written with no
    filename override) both carry the product label `"MFLAT"`, because a
    single-lamp reduction never enters the multi-lamp UV-stitch branch that
    would otherwise be the only consumer of the tag-scoped row: this is
    pinned as found, not as intended.
    """
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=131)
    flatFrame = synthetic_ccd(seed=132, prepared=True)
    firstStack = synthetic_ccd(seed=133, prepared=True)
    secondStack = synthetic_ccd(seed=134, prepared=True)
    backgroundSubtractedFlat = synthetic_ccd(seed=135, prepared=True)
    backgroundFrame = synthetic_ccd(seed=136, prepared=True)
    maskedFrame = synthetic_ccd(seed=137, prepared=True)

    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=True,
        calibratedFlatFiles=["flat-1.fits"],
        dFlatFiles=[],
        qFlatFiles=[],
        domeFlatFiles=[],
    )

    stackResults = iter([firstStack, secondStack])
    medianFlux = pd.DataFrame({"order": [10], "medianFlux": [100.0]})
    monkeypatch.setattr(recipe, "calibrate_frame_set", lambda: ([flatFrame], [], [], []))
    monkeypatch.setattr(recipe, "normalise_flats", lambda *a, **k: [flatFrame.copy()])
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **k: next(stackResults))
    monkeypatch.setattr(recipe, "mask_low_sens_pixels", lambda **k: (maskedFrame, medianFlux))
    _stub_shared_collaborators(recipe, monkeypatch, orderPath)

    _RecordingSubtractBackground.instances = []

    class FakeSubtractBackground(_RecordingSubtractBackground):
        def subtract(self) -> tuple[Any, Any, pd.DataFrame]:
            return backgroundFrame, backgroundSubtractedFlat, recipe.products

    monkeypatch.setattr(soxs_mflat_module, "subtract_background", FakeSubtractBackground)

    backgroundProductPath = tmp_path / "MASTER_FLAT_VIS_BKGROUND.fits"
    tagProductPath = tmp_path / "MASTER_FLAT_VIS.fits"
    finalProductPath = tmp_path / "MASTER_FLAT_VIS_final.fits"
    fake_write, writeCalls = _sequential_write_stub(
        [str(backgroundProductPath), str(tagProductPath), str(finalProductPath)]
    )
    monkeypatch.setattr(recipe, "_write", fake_write)

    # ACT
    returnedPath, returnedQc = recipe.produce_product()

    # ASSERT
    assert returnedPath == str(finalProductPath)
    assert returnedQc is recipe.qc
    _assert_background_subtraction_products(
        recipe,
        backgroundProductPath=backgroundProductPath,
        tagProductPath=tagProductPath,
        finalProductPath=finalProductPath,
    )

    # THE BACKGROUND FRAME'S HEADER IS A DEEPCOPY OF THE COMBINED FLAT'S
    # HEADER AT THE TIME OF WRITING, NOT THE SAME OBJECT.
    writtenBackgroundFrame = writeCalls[0][0][0]
    assert writtenBackgroundFrame is backgroundFrame
    assert writtenBackgroundFrame.header == backgroundSubtractedFlat.header
    assert writtenBackgroundFrame.header is not backgroundSubtractedFlat.header
    assert writeCalls[0][1]["filename"] == "MASTER_FLAT_VIS_BKGROUND.fits"
    assert writeCalls[1][1]["filename"] == "MASTER_FLAT_VIS.fits"
    assert "filename" not in writeCalls[2][1]

    receivedKwargs = _RecordingSubtractBackground.instances[0].kwargs
    assert receivedKwargs["frame"] is secondStack
    assert receivedKwargs["sofName"] == "MASTER_FLAT_VIS"
    assert receivedKwargs["recipeName"] == "soxs-mflat"
    assert receivedKwargs["lamp"] == ""


def _assert_multi_lamp_products(recipe: soxs_mflat) -> None:
    """Pin the tag-scoped product rows, their descriptions, and the table's column order."""
    assert list(recipe.products.columns) == _PRODUCT_COLUMNS_WITH_DESC
    assert list(recipe.products["product_label"]) == [
        "ORDER_LOC",
        "MFLAT",
        "ORDER_LOC_DLAMP",
        "MFLAT_DLAMP",
        "ORDER_LOC_QLAMP",
        "MFLAT_QLAMP",
        "MFLAT",
    ]

    dlampRow = recipe.products.loc[recipe.products["product_label"] == "MFLAT_DLAMP"].iloc[0]
    assert dlampRow["product_desc"] == "VIS master spectroscopic flat frame (DLAMP)"
    assert dlampRow["label"] == "PROD"

    qlampRow = recipe.products.loc[recipe.products["product_label"] == "MFLAT_QLAMP"].iloc[0]
    assert qlampRow["product_desc"] == "VIS master spectroscopic flat frame (QLAMP)"
    assert qlampRow["label"] == "PROD"

    unTaggedRows = recipe.products.loc[recipe.products["product_label"] == "MFLAT"]
    assert list(unTaggedRows["product_desc"]) == [
        "VIS master spectroscopic flat frame",
        "VIS master spectroscopic flat frame",
    ]

    for reductionDate in recipe.products["reduction_date_utc"]:
        assert TIMESTAMP_PATTERN.fullmatch(reductionDate) or reductionDate == "2024-01-02T04:05:06"


def test_multi_lamp_master_flat_records_tagged_dlamp_and_qlamp_product_rows(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Three lamp sets each write a tagged product row before the UV stitch step."""
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=231)
    productPath = tmp_path / "MASTER_FLAT_VIS.fits"
    frames = [synthetic_ccd(seed=seed, prepared=True) for seed in range(232, 235)]
    stitched = synthetic_ccd(seed=236, prepared=True)

    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=False,
        calibratedFlatFiles=["orderdef.fits"],
        dFlatFiles=["dorderdef.fits"],
        qFlatFiles=["qorderdef.fits"],
        domeFlatFiles=[],
    )

    def fake_normalise(*args: object, **kwargs: object) -> list[object]:
        lamp = str(kwargs["lamp"])
        frameIndex = {"": 0, "_DLAMP": 1, "_QLAMP": 2}[lamp]
        return [frames[frameIndex].copy()]

    medianFlux = pd.DataFrame({"order": [10], "medianFlux": [100.0]})
    monkeypatch.setattr(
        recipe,
        "calibrate_frame_set",
        lambda: ([frames[0]], [frames[1]], [frames[2]], []),
    )
    monkeypatch.setattr(recipe, "normalise_flats", fake_normalise)
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **k: k["frames"][0].copy())
    monkeypatch.setattr(recipe, "mask_low_sens_pixels", lambda **k: (k["frame"].copy(), medianFlux))
    monkeypatch.setattr(recipe, "_write", lambda *a, **k: str(productPath))
    monkeypatch.setattr(recipe, "stitch_uv_mflats", lambda *a, **k: stitched)
    _stub_shared_collaborators(recipe, monkeypatch, orderPath)

    # ACT
    returnedPath, returnedQc = recipe.produce_product()

    # ASSERT
    assert returnedPath == str(productPath)
    assert returnedQc is recipe.qc
    _assert_multi_lamp_products(recipe)
