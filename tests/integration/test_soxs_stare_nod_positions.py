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

from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import position_stack, reduction_inputs, stack_extraction
from tests.factories import synthetic_ccd
from tests.integration.test_soxs_stare_nod_mode import MARKER, SOF_NAME, _tail_inputs, _tail_recipe

pytestmark = pytest.mark.integration

FLUX = {"A": 10.0, "B": 20.0}
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
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, notes: list[tuple[str, str]] | None = None
) -> tuple[soxs_stare, list[Any]]:
    """Return a tail recipe whose extraction reports one spectrum for each position, and the extraction calls."""
    recipe, captured, _ = _tail_recipe(log, tmp_path, monkeypatch)
    recipe.productDir = str(tmp_path)
    extractions: list[Any] = []

    class Extractor:
        slitEdgeOrders: list = []
        flatPassHandoff = None
        traceFallback = False

        def __init__(self, **kwargs: object) -> None:
            captured["extractor"].append(kwargs)
            self.position = str(kwargs["sofName"]).split("_")[-2]
            self.extractionNotes = list(notes or [])

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[str, int], str]:
            extractions.append(self.position)
            return (
                recipe.qc,
                recipe.products,
                _spectrum(FLUX[self.position]),
                JOINS[self.position],
                f"extracted_{self.position}.fits",
            )

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
    assert Path(productPath).name == f"{SOF_NAME}{MARKER}_EXTRACTED_MERGED.fits"
    assert list(Table.read(productPath)["FLUX_COUNTS"]) == pytest.approx([15.0])
    assert recipe.sofName == SOF_NAME
    labels = recipe.products["product_label"]
    assert labels[labels.str.startswith("EXTRACTED_MERGED_TABLE")].tolist() == [
        "EXTRACTED_MERGED_TABLE",
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
