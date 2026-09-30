"""Characterization of `soxs_nod.process_single_ab_nodding_cycle` and `stack_extractions`.

These two methods are due to split (`process_single_ab_nodding_cycle`, 216
lines) and to have their timestamp and two product rows converted to shared
helpers (`stack_extractions`). These tests pin call order, arguments, the
frames written to disk, and the exact numeric output of both methods before
that refactor, so a later commit can prove it changed nothing observable.

Where the pinned behaviour looks like a defect, the test name says so and it
is reported rather than "fixed" here — characterization tests describe what
the code does today, not what it should do.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits

from soxspipe.recipes import soxs_nod
from tests.factories import pipeline_settings, qc_row, qc_table, synthetic_ccd

pytestmark = pytest.mark.integration


def _cycle_recipe(
    log: Any,
    tmp_path: Path,
    *,
    recipeName: str = "soxs-nod",
    saveSingleFrameExtractions: object = False,
) -> soxs_nod:
    """Return an unconstructed recipe configured for `process_single_ab_nodding_cycle`."""
    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.recipeName = recipeName
    recipe.sofName = "SYNTHETIC"
    recipe.debug = False
    recipe.qc = pd.DataFrame()
    recipe.products = pd.DataFrame()
    recipe.settings = {"instrument": "soxs"}
    recipe.recipeSettings = {"save_single_frame_extractions": saveSingleFrameExtractions}
    recipe.twoDMap = "two-d-map.fits"
    recipe.dispMap = "dispersion-map.fits"
    recipe.startNightDate = "2024-01-02"
    recipe.turnOffMP = True
    recipe.arm = "VIS"
    productDir = tmp_path / "products"
    productDir.mkdir(exist_ok=True)
    recipe.productDir = str(productDir)
    return recipe


def _ab_frames(
    *,
    aOverrides: dict[str, object] | None = None,
    bOverrides: dict[str, object] | None = None,
) -> tuple[Any, Any]:
    """Return one deterministic A frame and one deterministic B frame."""
    aFrame = synthetic_ccd(shape=(3, 3), seed=11, prepared=True, headerOverrides=aOverrides or {})
    bFrame = synthetic_ccd(shape=(3, 3), seed=12, prepared=True, headerOverrides=bOverrides or {})
    return aFrame, bFrame


def _patch_extraction_collaborators(
    monkeypatch: pytest.MonkeyPatch,
    recipe: soxs_nod,
    *,
    joins: list[dict[int, int]] | None = None,
) -> tuple[list[dict[str, object]], list[tuple[str, dict[str, object]]]]:
    """Patch `horne_extraction` and the two QC functions, recording call order and kwargs."""
    commonutils = import_module("soxspipe.commonutils")
    nodModule = import_module("soxspipe.recipes.soxs_nod")
    extractCalls: list[dict[str, object]] = []
    qcCalls: list[tuple[str, dict[str, object]]] = []
    joinSequence = joins or [{10: 1}]

    class RecordingExtractor:
        def __init__(self, **kwargs: object) -> None:
            extractCalls.append(kwargs)

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[int, int], str]:
            index = len(extractCalls) - 1
            spectrum = pd.DataFrame({"WAVE": [500.0 + index]})
            join = joinSequence[min(index, len(joinSequence) - 1)]
            return recipe.qc, recipe.products, spectrum, join, f"extraction_{index}.fits"

    def record_generic_checks(**kwargs: object) -> Any:
        qcCalls.append(("generic", kwargs))
        return kwargs["qcTable"]

    def record_spectroscopic_checks(**kwargs: object) -> Any:
        qcCalls.append(("spectroscopic", kwargs))
        return kwargs["qcTable"]

    monkeypatch.setattr(commonutils, "horne_extraction", RecordingExtractor)
    monkeypatch.setattr(nodModule, "generic_quality_checks", record_generic_checks)
    monkeypatch.setattr(nodModule, "spectroscopic_image_quality_checks", record_spectroscopic_checks)
    return extractCalls, qcCalls


def _patch_detrend(monkeypatch: pytest.MonkeyPatch, recipe: soxs_nod) -> list[dict[str, object]]:
    """Replace `detrend` with a recorder that returns a copy of the frame it was given."""
    detrendCalls: list[dict[str, object]] = []

    def record_detrend(**kwargs: object) -> Any:
        detrendCalls.append(kwargs)
        return kwargs["inputFrame"].copy()

    monkeypatch.setattr(recipe, "detrend", record_detrend)
    return detrendCalls


def test_nod_cycle_writes_ab_and_ba_difference_frames_with_expected_data_and_headers(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The nod branch writes both directional differences, each stamped with its own source header."""
    recipe = _cycle_recipe(log, tmp_path)
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames(aOverrides={"OBJECT": "FRAME-A"}, bOverrides={"OBJECT": "FRAME-B"})

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    productDir = Path(recipe.productDir)
    assert sorted(path.name for path in productDir.iterdir()) == [
        "SYNTHETIC_AB_1.fits",
        "SYNTHETIC_BA_1.fits",
    ]
    with fits.open(productDir / "SYNTHETIC_AB_1.fits") as hdus:
        assert hdus[0].header["OBJECT"] == "FRAME-A"
        assert hdus[0].data == pytest.approx(aFrame.data - bFrame.data, rel=1e-12)
    with fits.open(productDir / "SYNTHETIC_BA_1.fits") as hdus:
        assert hdus[0].header["OBJECT"] == "FRAME-B"
        assert hdus[0].data == pytest.approx(bFrame.data - aFrame.data, rel=1e-12)


def test_nod_cycle_prints_the_flattened_cycle_banner(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The flattened-pass banner has a trailing space where the not-flattened text would go."""
    recipe = _cycle_recipe(log, tmp_path)
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    assert ("print", "\n# PROCESSING AB NODDING CYCLE 1 ") in log.messages


def test_nod_cycle_prints_not_flattened_extra_text_when_requested(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The not-flattened pass appends its efficiency-calculation explanation to the banner."""
    recipe = _cycle_recipe(log, tmp_path)
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=2,
        orderTablePath="orders.fits",
        notFlattened=True,
    )

    expected = "\n# PROCESSING AB NODDING CYCLE 2  (not flattened this time - needed to calculate efficiency)"
    assert ("print", expected) in log.messages


def test_offset_named_recipe_writes_only_the_onoff_frame_with_a_single_extraction(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A recipe name without `nod` takes the offset branch: one file, one extraction, no B spectrum."""
    recipe = _cycle_recipe(log, tmp_path, recipeName="soxs-offset")
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    _, mergedB, _ = recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    assert sorted(path.name for path in Path(recipe.productDir).iterdir()) == ["SYNTHETIC_ONOFF_1.fits"]
    assert ("print", "\n# PROCESSING ON-OFF OFFSET CYCLE 1 ") in log.messages
    assert mergedB is False
    assert len(extractCalls) == 1


def test_quality_checks_run_generic_then_spectroscopic_on_the_ab_difference_frame(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both QC functions see the A-minus-B frame, generic first, and only spectroscopic gets the order table."""
    recipe = _cycle_recipe(log, tmp_path)
    _, qcCalls = _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    assert [name for name, _ in qcCalls] == ["generic", "spectroscopic"]
    expectedData = aFrame.data - bFrame.data
    assert qcCalls[0][1]["frame"].data == pytest.approx(expectedData, rel=1e-12)
    assert qcCalls[1][1]["frame"].data == pytest.approx(expectedData, rel=1e-12)
    assert qcCalls[1][1]["orderTablePath"] == "orders.fits"


def test_nod_cycle_detrends_four_frames_in_order_when_flattening_with_a_master_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A real master flat with flattening requested detrends A-B, B-A, bFrame, then aFrame, in that order."""
    recipe = _cycle_recipe(log, tmp_path)
    detrendCalls = _patch_detrend(monkeypatch, recipe)
    extractCalls, _ = _patch_extraction_collaborators(monkeypatch, recipe, joins=[{10: 1}, {20: 2}])
    aFrame, bFrame = _ab_frames()
    masterFlat = synthetic_ccd(shape=(3, 3), seed=99, prepared=True)

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        masterFlat=masterFlat,
    )

    assert len(detrendCalls) == 4
    expectedInputs = [aFrame.data - bFrame.data, bFrame.data - aFrame.data, bFrame.data, aFrame.data]
    for call, expected in zip(detrendCalls, expectedInputs, strict=True):
        assert call["inputFrame"].data == pytest.approx(expected, rel=1e-12)
        assert call["master_flat"] is masterFlat
        assert call["order_table"] == "orders.fits"
    assert extractCalls[0]["skySubtractedFrame"].data == pytest.approx(detrendCalls[0]["inputFrame"].data, rel=1e-12)
    assert extractCalls[0]["subtractedFrame"].data == pytest.approx(detrendCalls[2]["inputFrame"].data, rel=1e-12)
    assert extractCalls[0]["unflattenedFrame"].data == pytest.approx(aFrame.data - bFrame.data, rel=1e-12)
    assert extractCalls[1]["subtractedFrame"].data == pytest.approx(detrendCalls[3]["inputFrame"].data, rel=1e-12)


def test_nod_cycle_skips_detrend_when_there_is_no_master_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A false master flat leaves the not-flattened differences as the extracted frames."""
    recipe = _cycle_recipe(log, tmp_path)
    detrendCalls: list[dict[str, object]] = []
    monkeypatch.setattr(recipe, "detrend", lambda **kwargs: detrendCalls.append(kwargs))
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits", masterFlat=False
    )

    assert detrendCalls == []


def test_nod_cycle_skips_detrend_when_not_flattened_even_with_a_master_flat(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`notFlattened=True` skips detrending regardless of a supplied master flat."""
    recipe = _cycle_recipe(log, tmp_path)
    detrendCalls: list[dict[str, object]] = []
    monkeypatch.setattr(recipe, "detrend", lambda **kwargs: detrendCalls.append(kwargs))
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()
    masterFlat = synthetic_ccd(shape=(3, 3), seed=99, prepared=True)

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        masterFlat=masterFlat,
        notFlattened=True,
    )

    assert detrendCalls == []


def test_offset_named_cycle_detrends_three_frames_without_a_ba_pass(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The offset branch detrends A-B, bFrame and aFrame; it has no B-A pass to detrend."""
    recipe = _cycle_recipe(log, tmp_path, recipeName="soxs-offset")
    detrendCalls = _patch_detrend(monkeypatch, recipe)
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()
    masterFlat = synthetic_ccd(shape=(3, 3), seed=99, prepared=True)

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        masterFlat=masterFlat,
    )

    assert len(detrendCalls) == 3


@pytest.mark.parametrize(
    ("saveSetting", "extractorSeesFalse", "productsReplaced"),
    [
        (True, False, True),
        (False, True, False),
        (None, False, False),
        (0, True, False),
        (1, False, True),
        ("yes", False, False),
    ],
)
def test_save_single_frame_extractions_e712_comparison_controls_the_products_table(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    saveSetting: object,
    extractorSeesFalse: bool,
    productsReplaced: bool,
) -> None:
    """`== False`/`== True` comparisons, not truthiness, decide the extractor's input and `self.products`.

    `save_single_frame_extractions` is compared with `==` rather than `is`, so `0` reads as `False` on the
    way in (the extractor gets `False`) but not on the way out (the result is dropped), and `1`/`None`/`"yes"`
    each land on a different side of that same asymmetry. This test pins the six-way split exactly as it
    behaves today; it looks like a footgun, not a specification.
    """
    recipe = _cycle_recipe(log, tmp_path, saveSingleFrameExtractions=saveSetting)
    originalProducts = pd.DataFrame({"marker": ["ORIGINAL"]})
    recipe.products = originalProducts
    productsSeen: list[object] = []
    commonutils = import_module("soxspipe.commonutils")
    nodModule = import_module("soxspipe.recipes.soxs_nod")
    calls = {"count": 0}

    class SequencedExtractor:
        def __init__(self, **kwargs: object) -> None:
            productsSeen.append(kwargs["productsTable"])
            calls["count"] += 1

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[int, int], str]:
            returned = pd.DataFrame({"marker": [f"RETURNED-{calls['count']}"]})
            return recipe.qc, returned, pd.DataFrame({"WAVE": [500.0]}), {10: 1}, "extraction.fits"

    monkeypatch.setattr(commonutils, "horne_extraction", SequencedExtractor)
    monkeypatch.setattr(nodModule, "generic_quality_checks", lambda **kwargs: kwargs["qcTable"])
    monkeypatch.setattr(nodModule, "spectroscopic_image_quality_checks", lambda **kwargs: kwargs["qcTable"])
    aFrame, bFrame = _ab_frames()

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    assert (productsSeen[0] is False) is extractorSeesFalse
    assert (recipe.products is originalProducts) is (not productsReplaced)
    if productsReplaced:
        assert recipe.products["marker"].tolist() == ["RETURNED-2"]


@pytest.mark.parametrize(("notFlattenedValue", "detrendRuns"), [(0, True), (None, False)])
def test_not_flattened_e712_comparison_controls_whether_detrend_runs(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    notFlattenedValue: object,
    detrendRuns: bool,
) -> None:
    """`notFlattened == False`, not `not notFlattened`, gates detrending: `0` runs it, `None` does not."""
    recipe = _cycle_recipe(log, tmp_path)
    detrendCalls = _patch_detrend(monkeypatch, recipe)
    _patch_extraction_collaborators(monkeypatch, recipe)
    aFrame, bFrame = _ab_frames()
    masterFlat = synthetic_ccd(shape=(3, 3), seed=99, prepared=True)

    recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame,
        bFrame=bFrame,
        locationSetIndex=1,
        orderTablePath="orders.fits",
        masterFlat=masterFlat,
        notFlattened=notFlattenedValue,
    )

    assert bool(detrendCalls) is detrendRuns


def test_process_cycle_returns_order_joins_from_the_last_extraction(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The returned `orderJoins` is the B-minus-A extraction's, overwriting the A-minus-B one."""
    recipe = _cycle_recipe(log, tmp_path)
    _patch_extraction_collaborators(monkeypatch, recipe, joins=[{10: 1}, {20: 2}])
    aFrame, bFrame = _ab_frames()

    _, _, orderJoins = recipe.process_single_ab_nodding_cycle(
        aFrame=aFrame, bFrame=bFrame, locationSetIndex=1, orderTablePath="orders.fits"
    )

    assert orderJoins == {20: 2}


@pytest.fixture
def passthrough_snr_qcs(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace `add_snr_efficiency_qcs` with a pass-through, which most stack tests do not exercise."""

    def pass_through(**kwargs: object) -> Any:
        return kwargs["qcTable"]

    monkeypatch.setattr(import_module("soxspipe.commonutils.toolkit"), "add_snr_efficiency_qcs", pass_through)


def _stack_recipe(log: Any, tmp_path: Path, *, recipeName: str = "soxs-nod") -> soxs_nod:
    """Return an unconstructed recipe configured for `stack_extractions`, with declared-empty tables."""
    recipe = soxs_nod.__new__(soxs_nod)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = recipeName
    recipe.sofName = "SYNTHETIC"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.productDir = str(tmp_path)
    # THE DECLARED EMPTY TABLES, NOT THE TEST-FACTORY ONES, SO THE PRODUCT-COLUMN
    # ORDER MATCHES WHAT A REAL RECIPE STARTS WITH.
    _, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.masterHeaderFrame = synthetic_ccd(shape=(3, 3), prepared=True)
    recipe.update_fits_keywords = lambda **_: None
    return recipe


def _spectrum_frames() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the two order-merged spectra whose medians and roundings this file pins."""
    first = pd.DataFrame(
        {
            "WAVE": [500.12344 * u.nm, 500.12346 * u.nm, 700.0 * u.nm],
            "FLUX_COUNTS": [9.0, 55.0, 42.0],
            "VARIANCE": [9.0, 4.0, 49.0],
            "FLUX_DENSITY_COUNTS": [1.111, 3.001, 9.5],
        }
    )
    second = pd.DataFrame(
        {
            "WAVE": [500.12345 * u.nm, 700.0 * u.nm],
            "FLUX_COUNTS": [15.6913, 58.0],
            "VARIANCE": [9.0, 49.0],
            "FLUX_DENSITY_COUNTS": [3.3334, 10.5],
        }
    )
    return first, second


def test_stack_extractions_computes_column_specific_rounded_medians(
    log: Any,
    tmp_path: Path,
    passthrough_snr_qcs: None,
) -> None:
    """4-dp WAVE rounding decides merge-vs-split, then each column keeps its own decimal precision.

    500.12344 and 500.12345 round to the same 4-dp key and merge (median); 500.12344 and 500.12346 do not
    and stay as two rows — even though both display as WAVE 500.12 at 2 dp, proving the two-row split, not
    just the displayed value. FLUX_COUNTS at 3 dp, SNR at 2 dp and FLUX_DENSITY_COUNTS at 3 dp each keep their
    own decimals rather than any other column's: the raw FLUX_COUNTS median 12.34565 would round to 12.35 at
    SNR's 2 dp, not the pinned 12.346.
    """
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    stacked, _ = recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert len(stacked) == 3
    assert stacked["WAVE"].tolist() == pytest.approx([500.12, 500.12, 700.0], rel=1e-12)
    assert stacked["FLUX_COUNTS"].tolist() == pytest.approx([12.346, 55.0, 50.0], rel=1e-12)
    assert stacked["FLUX_DENSITY_COUNTS"].tolist() == pytest.approx([2.222, 3.001, 10.0], rel=1e-12)
    assert stacked["SNR"].tolist() == pytest.approx([4.12, 27.5, 7.14], rel=1e-12)


def test_stack_extractions_writes_notflat_suffixed_files_when_requested(
    log: Any,
    tmp_path: Path,
    passthrough_snr_qcs: None,
) -> None:
    """`notFlattened=True` names both product files with the `_NOTFLAT` suffix."""
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    _, fitsPath = recipe.stack_extractions([first, second], notFlattened=True, orderJoins={10: 1})

    assert Path(fitsPath).name == "SYNTHETIC_EXTRACTED_MERGED_NOTFLAT.fits"
    assert (tmp_path / "SYNTHETIC_EXTRACTED_MERGED_NOTFLAT.txt").is_file()
    assert recipe.filenameTemplate == "SYNTHETIC.fits"


def test_stack_extractions_sets_the_reduced_header_keywords_and_date_obs(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    passthrough_snr_qcs: None,
) -> None:
    """The written header carries PRO_TYPE REDUCED and PRO_CATG SCI_SLIT_FLUX_<ARM>, and dateObs comes from it."""
    recipe = _stack_recipe(log, tmp_path)
    keywordCalls: list[dict[str, object]] = []
    monkeypatch.setattr(
        recipe, "update_fits_keywords", lambda **kwargs: keywordCalls.append(kwargs)
    )
    first, second = _spectrum_frames()

    _, fitsPath = recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert keywordCalls == [{"frame": recipe.masterHeaderFrame}]
    with fits.open(fitsPath) as hdus:
        assert hdus[0].header["ESO PRO TYPE"] == "REDUCED"
        assert hdus[0].header["ESO PRO CATG"] == "SCI_SLIT_FLUX_VIS"
    assert recipe.dateObs == recipe.masterHeaderFrame.header["DATE-OBS"]


def test_stack_extractions_writes_wave_in_angstrom_with_two_decimal_formatting(
    log: Any,
    tmp_path: Path,
    passthrough_snr_qcs: None,
) -> None:
    """The ASCII companion converts WAVE to Angstrom (x10) and formats it to 2 decimal places."""
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    _, fitsPath = recipe.stack_extractions([first, second], orderJoins={10: 1})

    asciiLines = Path(fitsPath.replace(".fits", ".txt")).read_text().splitlines()
    assert asciiLines[0] == "WAVE FLUX_COUNTS VARIANCE FLUX_DENSITY_COUNTS SNR"
    assert asciiLines[1].startswith("5001.20 ")
    assert asciiLines[3].startswith("7000.00 ")


def test_stack_extractions_appends_product_rows_that_push_desc_to_the_last_column(
    log: Any,
    tmp_path: Path,
    passthrough_snr_qcs: None,
) -> None:
    """Concatenating onto the declared empty table moves `product_desc` after `label`, not before `file_path`.

    The declared empty product table (`base_recipe._empty_qc_and_product_tables`) has no `product_desc`
    column, unlike the test factory's `product_table()`. `pd.concat` appends brand-new columns after the
    receiving frame's own columns, so the literal dict order in the source — `product_desc` before
    `file_path` — is not the column order a real recipe ends up with. This looks like a defect worth a
    human's attention before the split.
    """
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    _, fitsPath = recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert list(recipe.products.columns) == [
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
    tableRow, asciiRow = recipe.products.to_dict("records")
    assert tableRow["product_label"] == "EXTRACTED_MERGED_TABLE"
    assert tableRow["file_type"] == "FITS"
    assert tableRow["product_desc"] == "Table of the extracted source in each order. All nodding cycles combined."
    assert tableRow["file_path"] == fitsPath
    assert tableRow["label"] == "PROD"
    assert tableRow["obs_date_utc"] == recipe.dateObs
    assert asciiRow["obs_date_utc"] == recipe.dateObs
    assert asciiRow["product_label"] == "EXTRACTED_MERGED_ASCII"
    assert asciiRow["file_type"] == "TXT"
    assert asciiRow["product_desc"] == "Ascii version of extracted source spectrum"
    assert asciiRow["file_path"] == fitsPath.replace(".fits", ".txt")


def test_stack_extractions_stamps_both_product_rows_with_one_whole_second_timestamp(
    log: Any,
    tmp_path: Path,
    passthrough_snr_qcs: None,
) -> None:
    """By contract, not by patching the clock: one whole-second UTC string feeds both rows and `self.utcnow`.

    The next commit moves the timestamp source from a function-local `datetime.utcnow()` to
    `toolkit.utcnow_string()`; a test that patches `datetime` would only pass against one revision. Asserting
    the shape and the shared value, without touching the clock, survives that change.
    """
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    recipe.stack_extractions([first, second], orderJoins={10: 1})

    reductionDates = recipe.products["reduction_date_utc"].tolist()
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}", reductionDates[0])
    assert len(set(reductionDates)) == 1
    assert recipe.utcnow == reductionDates[0]


def test_stack_extractions_reads_the_clock_exactly_once(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Both product rows share one reading of `utcnow_string`, not one reading each.

    A format assertion cannot catch a second mint, so this stubs the clock at the name the module imported
    and counts the reads. It holds only after the conversion, which is why it is separate from the contract
    test written against the old source.
    """
    toolkit = import_module("soxspipe.commonutils.toolkit")
    nodModule = import_module("soxspipe.recipes.soxs_nod")
    monkeypatch.setattr(toolkit, "add_snr_efficiency_qcs", lambda **kwargs: kwargs["qcTable"])
    reads: list[str] = []

    def counting_clock(**kwargs: object) -> str:
        reads.append("read")
        return f"2024-01-02T03:04:0{len(reads)}"

    monkeypatch.setattr(nodModule, "utcnow_string", counting_clock)
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert reads == ["read"]
    assert recipe.products["reduction_date_utc"].tolist() == ["2024-01-02T03:04:01"] * 2


@pytest.mark.parametrize("recipeName", ["soxs-nod-std", "soxs-offset"])
def test_stack_extractions_stamps_product_rows_with_the_current_recipe_name(
    log: Any,
    tmp_path: Path,
    recipeName: str,
    passthrough_snr_qcs: None,
) -> None:
    """`stack_extractions` is inherited by `soxs_offset`; both recipe names flow into its product rows."""
    recipe = _stack_recipe(log, tmp_path, recipeName=recipeName)
    first, second = _spectrum_frames()

    recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert recipe.products["soxspipe_recipe"].tolist() == [recipeName, recipeName]


def test_stack_extractions_calls_add_snr_efficiency_qcs_and_writes_its_qc_to_the_header(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`add_snr_efficiency_qcs` sees the grouped frame and order joins, and its return reaches the FITS header."""
    toolkit = import_module("soxspipe.commonutils.toolkit")
    qcCalls: list[dict[str, object]] = []
    returnedQc = qc_row(
        recipeName="soxs-nod",
        name="SYNTHETIC QC",
        value=42.0,
        unit="electron",
        comment="Synthetic QC metric",
    )

    def fake_add_snr(**kwargs: object) -> pd.DataFrame:
        qcCalls.append(kwargs)
        return returnedQc

    monkeypatch.setattr(toolkit, "add_snr_efficiency_qcs", fake_add_snr)
    recipe = _stack_recipe(log, tmp_path)
    first, second = _spectrum_frames()

    _, fitsPath = recipe.stack_extractions([first, second], orderJoins={10: 1})

    assert len(qcCalls) == 1
    assert qcCalls[0]["orderJoins"] == {10: 1}
    assert qcCalls[0]["recipeName"] == "soxs-nod"
    assert qcCalls[0]["dateObs"] == recipe.dateObs
    assert list(qcCalls[0]["spectrumDF"]["SNR"]) == pytest.approx([4.12, 27.5, 7.14], rel=1e-12)
    assert recipe.qc is returnedQc
    with fits.open(fitsPath) as hdus:
        assert hdus[0].header["ESO QC SYNTHETIC QC"] == pytest.approx(42.0)
