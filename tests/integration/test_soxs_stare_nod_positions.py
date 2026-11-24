"""The two nod positions of a nodding SOF are extracted, combined and finished as one stare spectrum (DY-2290).

`_reduce_nod_positions` extracts the A stack, then the B stack, combines the two extractions the way the nod recipe
combines its cycles, and flux-calibrates, writes and plots the combined spectrum.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.table import Table

from soxspipe.commonutils.toolkit import append_qc
from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import position_stack, reduction_inputs, stack_extraction
from tests.factories import synthetic_ccd
from tests.integration.test_soxs_stare_nod_mode import MARKER, SOF_NAME, _tail_inputs, _tail_recipe

pytestmark = pytest.mark.integration

FLUX = {"A": 10.0, "B": 20.0}
# THE NOT-FLATTENED SPECTRUM OF A POSITION IS THIS MANY TIMES ITS FLAT-FIELDED SPECTRUM
NOT_FLAT_SCALE = 3.0
# A QC THE EXTRACTION OF EACH POSITION ADDS, AND ITS VALUE FOR EACH POSITION
POSITION_QC = "N ORDERS SLIT EDGE"
QC_VALUE = {"A": 1, "B": 2}
JOINS = {"A": {"10-11": 11}, "B": {"10-11": 12}}


def _spectrum(flux: float) -> pd.DataFrame:
    """Return an order-merged spectrum of one wavelength, shaped like the extraction's output."""
    return pd.DataFrame(
        {
            "WAVE": [500.0 * u.nm],
            "FLUX_COUNTS": [flux],
            "VARIANCE": [4.0],
            "SKY_COUNTS": [1.0],
            "SNR": [0.0],
            "FLUX_DENSITY_COUNTS": [flux],
        }
    )


def _stacks() -> list[position_stack]:
    stacks = []
    for position, seed in (("A", 10), ("B", 20)):
        notFlattened = synthetic_ccd(seed=seed, prepared=True)
        flattened = synthetic_ccd(seed=seed + 1, prepared=True)
        uncovered = np.zeros(notFlattened.data.shape, dtype=bool)
        stacks.append(
            position_stack(position, [synthetic_ccd(seed=seed + 2, prepared=True)], notFlattened, flattened, uncovered)
        )
    return stacks


def _recipe(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    *,
    notes: list[tuple[str, str]] | None = None,
    fallbacks: tuple[str, ...] = (),
) -> tuple[soxs_stare, list[Any]]:
    """Return a tail recipe whose extraction reports one spectrum for each position, and the extraction calls.

    A position in ``fallbacks`` finds no object trace, so its flat-fielded extraction is a boxcar and records so in the
    ``EXTRACT METHOD`` QC, as the Horne extraction does.
    """
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)
    recipe.productDir = str(tmp_path)
    extractions: list[Any] = []

    class Extractor:
        slitEdgeOrders: list = []

        def __init__(self, **kwargs: object) -> None:
            captured["extractor"].append(kwargs)
            self.position = str(kwargs["sofName"]).split("_")[-2]
            self.traceFallback = self.position in fallbacks
            self.isNotFlat = bool(kwargs.get("notFlattened"))
            self.extractionNotes = list(notes or [])
            self.flatPassHandoff = f"handoff-{self.position}"

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[str, int], str]:
            extractions.append(f"{self.position}-notflat" if self.isNotFlat else self.position)
            if not self.isNotFlat:
                recipe.qc = append_qc(
                    recipe.qc,
                    recipeName=recipe.recipeName,
                    qcName="EXTRACT METHOD",
                    qcValue="boxcar" if self.traceFallback else "optimal",
                    qcComment="Optimal, or boxcar if no trace",
                    obsDateUtc=recipe.dateObs,
                    reductionDateUtc="2024-01-02T03:04:05",
                    qcUnit=None,
                    toHeader=True,
                )
                recipe.qc = append_qc(
                    recipe.qc,
                    recipeName=recipe.recipeName,
                    qcName=POSITION_QC,
                    qcValue=QC_VALUE[self.position],
                    qcComment="A QC the extraction of a position adds",
                    obsDateUtc=recipe.dateObs,
                    reductionDateUtc="2024-01-02T03:04:05",
                    qcUnit=None,
                    toHeader=True,
                )
            flux = FLUX[self.position] * (NOT_FLAT_SCALE if self.isNotFlat else 1.0)
            name = f"notflat_{self.position}.fits" if self.isNotFlat else f"extracted_{self.position}.fits"
            return recipe.qc, recipe.products, _spectrum(flux), JOINS[self.position], name

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", Extractor)
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"), "add_snr_efficiency_qcs", lambda **kwargs: kwargs["qcTable"]
    )
    return recipe, extractions


def _spy_finish(recipe: soxs_stare, monkeypatch: pytest.MonkeyPatch) -> dict[str, list[tuple[Any, ...]]]:
    calls: dict[str, list[tuple[Any, ...]]] = {"fluxCal": [], "ascii": [], "plot": []}
    monkeypatch.setattr(
        recipe, "_flux_calibrate_spectrum", lambda *args: calls["fluxCal"].append(args) or "fluxcal.fits"
    )
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: calls["ascii"].append(args))
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: calls["plot"].append(args))
    return calls


def test_the_a_stack_is_extracted_before_the_b_stack_and_the_caller_keeps_neither(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, extractions = _recipe(log, tmp_path, monkeypatch)
    _spy_finish(recipe, monkeypatch)
    stacks = _stacks()

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT: THE LIST IS EMPTIED AS THE STACKS ARE EXTRACTED, SO A STACK IS FREED ONCE ITS POSITION IS DONE
    assert extractions == ["A", "B"]
    assert stacks == []


def test_the_two_extractions_are_combined_with_the_order_joins_of_the_last_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _spy_finish(recipe, monkeypatch)
    stacks = _stacks()
    stackA = stacks[0]
    combined: list[tuple[Any, ...]] = []

    def stack_extractions(dataFrameList: list[Any], notFlattened: bool = False, orderJoins: Any = None) -> Any:
        combined.append((dataFrameList, notFlattened, orderJoins, recipe.masterHeaderFrame))
        return Table({"WAVE": [500.0]}), "combined.fits"

    monkeypatch.setattr(recipe, "stack_extractions", stack_extractions)

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT: THE HEADER OF THE COMBINED PRODUCT COMES FROM THE FIRST POSITION
    [(spectra, notFlattened, orderJoins, headerFrame)] = combined
    assert [list(spectrum["FLUX_COUNTS"]) for spectrum in spectra] == [[FLUX["A"]], [FLUX["B"]]]
    assert notFlattened is False
    assert orderJoins == JOINS["B"]
    assert headerFrame is stackA.flattened


def test_the_combined_spectrum_is_what_flux_calibration_ascii_and_the_plots_get(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    calls = _spy_finish(recipe, monkeypatch)
    stacks = _stacks()
    stackedSpectrum = Table({"WAVE": [500.0]})
    monkeypatch.setattr(recipe, "stack_extractions", lambda *_, **__: (stackedSpectrum, "combined.fits"))
    allFrames = [synthetic_ccd(seed=1, prepared=True), synthetic_ccd(seed=2, prepared=True)]
    inputs = reduction_inputs(False, False, False, "order.fits", "disp.fits", "map.fits", "resp.fits")
    stackA = stacks[0]

    # ACT
    productPath, failure = recipe._reduce_nod_positions(stacks, allFrames, False, inputs)

    # ASSERT: THE RESPONSE USES THE COMBINED SPECTRUM, THE A HEADER AND THE FRAMES OF BOTH POSITIONS
    assert (productPath, failure) == ("combined.fits", False)
    assert calls["fluxCal"] == [("resp.fits", stackedSpectrum, stackA.flattened, allFrames)]
    assert calls["ascii"] == [("combined.fits", "fluxcal.fits")]
    assert calls["plot"] == [(stackedSpectrum, JOINS["B"], "fluxcal.fits", [], False)]


def test_the_extraction_notes_of_both_positions_are_tagged_and_plotted_together(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, notes=[("WARNING", "an order was skipped")])
    calls = _spy_finish(recipe, monkeypatch)
    monkeypatch.setattr(recipe, "stack_extractions", lambda *_, **__: (Table({"WAVE": [500.0]}), "combined.fits"))

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert calls["plot"][0][3] == [("WARNING", "[A] an order was skipped"), ("WARNING", "[B] an order was skipped")]


def test_the_final_product_is_the_mean_of_the_two_positions_under_the_unmarked_position_name(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _spy_finish(recipe, monkeypatch)

    # ACT
    productPath, _ = recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT: THE SOF NAME IS THE ONE THE RUN STARTED WITH, SO THE PRODUCT HAS NO POSITION IN ITS NAME
    assert Path(productPath).name == f"{SOF_NAME}{MARKER}_SPECTRUM.fits"
    assert list(Table.read(productPath)["FLUX_COUNTS"]) == pytest.approx([15.0])
    assert recipe.sofName == SOF_NAME
    labels = recipe.products["product_label"]
    assert labels[labels.str.startswith("SPECTRUM_TABLE")].tolist() == [
        "SPECTRUM_TABLE",
    ]


def test_a_slit_edge_spill_in_either_position_is_reported_for_the_combined_spectrum(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    spilled: list[bool] = []
    monkeypatch.setattr(recipe, "stack_extractions", lambda *_, **__: (Table({"WAVE": [500.0]}), "combined.fits"))
    monkeypatch.setattr(
        recipe, "_finish_extraction", lambda extraction, *_: spilled.append(extraction.slitEdgeSpill) or False
    )
    realExtract = recipe._extract_position

    def extract_position(stack: position_stack, inputs: Any, rawFrames: Any) -> stack_extraction:
        extraction = realExtract(stack, inputs, rawFrames)
        return extraction._replace(slitEdgeSpill=stack.position == "B")

    monkeypatch.setattr(recipe, "_extract_position", extract_position)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert spilled == [True]


def _standard(recipe: soxs_stare) -> None:
    recipe.generateReponseCurve = True


def test_a_standard_position_extracts_its_not_flattened_frame_with_its_own_flat_pass(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, extractions = _recipe(log, tmp_path, monkeypatch)
    _standard(recipe)
    stack = _stacks()[1]

    # ACT
    extraction = recipe._extract_position(stack, _tail_inputs(), False)

    # ASSERT: THE FLAT PASS RUNS FIRST, THEN THE NOT-FLATTENED PASS WITH THE HANDOFF OF THAT FLAT PASS
    assert extractions == ["B", "B-notflat"]
    assert extraction.responseInputs is None
    assert extraction.notFlatPath == "notflat_B.fits"
    assert list(extraction.notFlatSpectrum["FLUX_COUNTS"]) == [FLUX["B"] * NOT_FLAT_SCALE]
    assert recipe.sofName == SOF_NAME


def test_the_not_flattened_pass_of_a_position_is_named_under_that_position(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _standard(recipe)
    built: list[dict[str, Any]] = []
    horne = import_module("soxspipe.commonutils").horne_extraction

    def record(**kwargs: Any) -> Any:
        built.append(kwargs)
        return horne(**kwargs)

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", record)

    # ACT
    recipe._extract_position(_stacks()[0], _tail_inputs(), False)

    # ASSERT
    assert [call["sofName"] for call in built] == [f"{SOF_NAME}_A{MARKER}"] * 2
    assert [bool(call.get("notFlattened")) for call in built] == [False, True]
    assert built[1]["flatPass"] == "handoff-A"


def test_a_science_position_makes_no_not_flattened_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, extractions = _recipe(log, tmp_path, monkeypatch)

    # ACT
    extraction = recipe._extract_position(_stacks()[0], _tail_inputs(), False)

    # ASSERT
    assert extractions == ["A"]
    assert extraction.notFlatPath is None
    assert extraction.notFlatSpectrum is None


def test_a_standard_with_a_response_table_makes_no_not_flattened_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, extractions = _recipe(log, tmp_path, monkeypatch)
    _standard(recipe)
    inputs = reduction_inputs(False, False, False, "order.fits", "disp.fits", "map.fits", "resp.fits")

    # ACT
    recipe._extract_position(_stacks()[0], inputs, False)

    # ASSERT
    assert extractions == ["A"]


def test_the_not_flattened_extractions_are_combined_and_give_one_response_curve(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _standard(recipe)
    combined: list[tuple[Any, ...]] = []
    responses: list[tuple[Any, ...]] = []

    def stack_extractions(dataFrameList: list[Any], notFlattened: bool = False, orderJoins: Any = None) -> Any:
        combined.append(([list(spectrum["FLUX_COUNTS"]) for spectrum in dataFrameList], notFlattened, orderJoins))
        return Table({"WAVE": [500.0]}), "notflat.fits" if notFlattened else "flat.fits"

    monkeypatch.setattr(recipe, "stack_extractions", stack_extractions)
    monkeypatch.setattr(recipe, "_compute_response", lambda *args: responses.append(args) or "response failed")
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: None)
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: None)
    allFrames = [synthetic_ccd(seed=1, prepared=True)]

    # ACT
    productPath, failure = recipe._reduce_nod_positions(_stacks(), allFrames, False, _tail_inputs())

    # ASSERT: THE FLAT-FIELDED POSITIONS ARE COMBINED, THEN THE NOT-FLATTENED ONES, AND ONE RESPONSE IS BUILT FROM BOTH
    assert combined == [
        ([[FLUX["A"]], [FLUX["B"]]], False, JOINS["B"]),
        ([[FLUX["A"] * NOT_FLAT_SCALE], [FLUX["B"] * NOT_FLAT_SCALE]], True, JOINS["B"]),
    ]
    assert responses == [("flat.fits", "notflat.fits", JOINS["B"], allFrames)]
    assert (productPath, failure) == ("flat.fits", "response failed")


def test_a_standard_that_spills_over_the_slit_edge_in_one_position_fails_the_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _standard(recipe)
    monkeypatch.setattr(recipe, "stack_extractions", lambda *_, **__: (Table({"WAVE": [500.0]}), "combined.fits"))
    monkeypatch.setattr(recipe, "_compute_response", lambda *args: False)
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: None)
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: None)
    realExtract = recipe._extract_position

    def extract_position(stack: position_stack, inputs: Any, rawFrames: Any) -> stack_extraction:
        return realExtract(stack, inputs, rawFrames)._replace(slitEdgeSpill=stack.position == "A")

    monkeypatch.setattr(recipe, "_extract_position", extract_position)

    # ACT
    _, failure = recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert failure is True


def _extract_methods(recipe: soxs_stare) -> list[Any]:
    """Return the values of the ``EXTRACT METHOD`` QC rows in the order they were recorded."""
    return recipe.qc.loc[recipe.qc["qc_name"] == "EXTRACT METHOD", "qc_value"].tolist()


def _record_combination(recipe: soxs_stare, monkeypatch: pytest.MonkeyPatch) -> list[tuple[Any, ...]]:
    """Replace `stack_extractions` with a recorder of the spectra, order joins and header frame it gets."""
    calls: list[tuple[Any, ...]] = []

    def stack_extractions(dataFrameList: list[Any], notFlattened: bool = False, orderJoins: Any = None) -> Any:
        spectra = [None if spectrum is None else list(spectrum["FLUX_COUNTS"]) for spectrum in dataFrameList]
        calls.append((spectra, notFlattened, orderJoins, recipe.masterHeaderFrame))
        return Table({"WAVE": [500.0]}), "notflat.fits" if notFlattened else "flat.fits"

    monkeypatch.setattr(recipe, "stack_extractions", stack_extractions)
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: None)
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: None)
    return calls


def test_a_position_with_no_object_trace_is_left_out_and_the_other_position_is_used_alone(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=("B",))
    combined = _record_combination(recipe, monkeypatch)
    plots: list[tuple[Any, ...]] = []
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: plots.append(args))
    stackA = _stacks()[0]
    stacks = [stackA, _stacks()[1]]

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT: B CONTRIBUTES NOTHING, AND THE HEADER AND THE ORDER JOINS COME FROM THE POSITION THAT WAS KEPT
    [(spectra, _, orderJoins, headerFrame)] = combined
    assert spectra == [[FLUX["A"]], None]
    assert orderJoins == JOINS["A"]
    assert headerFrame is stackA.flattened
    message = (
        "No object trace was found in the B frames, so the stare-mode spectrum is the A position's extraction alone"
    )
    assert any(level == "warning" and message in text for level, text in log.messages)
    assert ("WARNING", message) in plots[0][3]


def test_the_header_and_order_joins_come_from_the_b_position_when_the_a_position_has_no_trace(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=("A",))
    combined = _record_combination(recipe, monkeypatch)
    stacks = _stacks()
    stackB = stacks[1]

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT
    [(spectra, _, orderJoins, headerFrame)] = combined
    assert spectra == [None, [FLUX["B"]]]
    assert orderJoins == JOINS["B"]
    assert headerFrame is stackB.flattened


def test_neither_position_having_a_trace_fails_as_nod_mode_does(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=("A", "B"))

    # ACT / ASSERT
    with pytest.raises(Exception, match="No object trace was found"):
        recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())
    assert any(level == "error" and "No object trace was found" in message for level, message in log.messages)


def test_the_extract_method_qc_says_optimal_when_one_position_fell_back_to_a_boxcar(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE B EXTRACTION RECORDS A BOXCAR ROW, BUT B IS LEFT OUT OF THE COMBINATION
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=("B",))
    _record_combination(recipe, monkeypatch)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert _extract_methods(recipe) == ["optimal"]


def test_the_extract_method_qc_is_left_as_the_extractions_recorded_it_when_both_positions_have_a_trace(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _record_combination(recipe, monkeypatch)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert _extract_methods(recipe) == ["optimal", "optimal"]


def test_a_standard_position_with_no_trace_makes_no_not_flattened_extraction_and_the_other_is_used(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, extractions = _recipe(log, tmp_path, monkeypatch, fallbacks=("B",))
    _standard(recipe)
    combined = _record_combination(recipe, monkeypatch)
    monkeypatch.setattr(recipe, "_compute_response", lambda *args: False)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT: B'S FLAT PASS RAN, ITS NOT-FLATTENED PASS DID NOT, AND THE NOT-FLATTENED COMBINATION HOLDS A ALONE
    assert extractions == ["A", "A-notflat", "B"]
    assert [call[:2] for call in combined] == [
        ([[FLUX["A"]], None], False),
        ([[FLUX["A"] * NOT_FLAT_SCALE], None], True),
    ]


def _qc_values(recipe: soxs_stare, name: str) -> list[Any]:
    return recipe.qc.loc[recipe.qc["qc_name"] == name, "qc_value"].tolist()


def _header_values(path: str, name: str) -> list[Any]:
    """Return the value of the ``ESO QC <name>`` card in each HDU of a FITS file that has one."""
    from astropy.io import fits

    key = f"ESO QC {name}".upper()
    with fits.open(path) as hdus:
        return [hdu.header[key] for hdu in hdus if key in hdu.header]


@pytest.mark.parametrize(("dropped", "kept"), [("B", "A"), ("A", "B")])
def test_the_qc_rows_of_a_position_with_no_trace_reach_neither_the_report_nor_the_product_header(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, dropped: str, kept: str
) -> None:
    # ARRANGE: THE REAL COMBINATION WRITES THE PRODUCT, WITH THE QC TABLE IN ITS HEADER
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=(dropped,))
    monkeypatch.setattr(recipe, "_write_ascii_spectrum", lambda *args: None)
    monkeypatch.setattr(recipe, "_plot_merged_spectrum_qcs", lambda *args: None)

    # ACT
    productPath, _ = recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT: ONLY THE KEPT POSITION'S ROWS SURVIVE, INCLUDING ITS OPTIMAL EXTRACT METHOD
    assert _qc_values(recipe, POSITION_QC) == [QC_VALUE[kept]]
    assert _extract_methods(recipe) == ["optimal"]
    assert _header_values(productPath, POSITION_QC) == [QC_VALUE[kept]]
    assert _header_values(productPath, "EXTRACT METHOD") == ["optimal"]


def test_the_qc_rows_of_both_positions_stay_when_both_have_a_trace(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    _record_combination(recipe, monkeypatch)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert _qc_values(recipe, POSITION_QC) == [QC_VALUE["A"], QC_VALUE["B"]]


def test_a_slit_edge_spill_in_a_position_left_out_does_not_count(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: B SPILLS OVER THE SLIT EDGE, BUT B HAS NO TRACE AND IS LEFT OUT
    recipe, _ = _recipe(log, tmp_path, monkeypatch)
    spilled: list[bool] = []
    monkeypatch.setattr(recipe, "stack_extractions", lambda *_, **__: (Table({"WAVE": [500.0]}), "combined.fits"))
    monkeypatch.setattr(
        recipe, "_finish_extraction", lambda extraction, *_: spilled.append(extraction.slitEdgeSpill) or False
    )
    realExtract = recipe._extract_position

    def extract_position(stack: position_stack, inputs: Any, rawFrames: Any) -> stack_extraction:
        extraction = realExtract(stack, inputs, rawFrames)
        return extraction._replace(slitEdgeSpill=stack.position == "B", traceFallback=stack.position == "B")

    monkeypatch.setattr(recipe, "_extract_position", extract_position)

    # ACT
    recipe._reduce_nod_positions(_stacks(), [], False, _tail_inputs())

    # ASSERT
    assert spilled == [False]
