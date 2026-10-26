"""A nod, nod_std or offset reduction whose extractions found no object trace fails with a clear error (DY-1617).

`horne_extraction.extract()` returns a `None` spectrum when it finds no trace. Stacking used to hand that straight
to `pd.concat`, which raised pandas' `ValueError: All objects passed were None`.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

from soxspipe.recipes import soxs_nod, soxs_offset
from tests.integration.test_soxs_nod_cycle_characterization import _ab_frames, _spectrum_frames, _stack_recipe

pytestmark = pytest.mark.integration


@pytest.fixture
def passthrough_snr_qcs(monkeypatch: pytest.MonkeyPatch) -> None:
    """Replace `add_snr_efficiency_qcs` with a pass-through, which these tests do not exercise."""

    def pass_through(**kwargs: object) -> Any:
        return kwargs["qcTable"]

    monkeypatch.setattr(import_module("soxspipe.commonutils.toolkit"), "add_snr_efficiency_qcs", pass_through)


@pytest.fixture
def snr_qc_order_joins(monkeypatch: pytest.MonkeyPatch) -> list[Any]:
    """Record the `orderJoins` each call of `add_snr_efficiency_qcs` receives, and return the recording."""
    received: list[Any] = []

    def record(**kwargs: object) -> Any:
        received.append(kwargs["orderJoins"])
        return kwargs["qcTable"]

    monkeypatch.setattr(import_module("soxspipe.commonutils.toolkit"), "add_snr_efficiency_qcs", record)
    return received


def _cycle_loop_recipe(recipeClass: type, log: Any, tmp_path: Path, recipeName: str) -> Any:
    """Return a recipe whose per-cycle loop runs against stub frames and a cycle that finds no trace."""
    recipe = _stack_recipe(log, tmp_path, recipeName=recipeName)
    recipe.__class__ = recipeClass
    recipe.recipeSettings = {"use_flat": False}
    recipe.process_single_ab_nodding_cycle = lambda **_: (None, None, None)
    return recipe


def _stub_frames(count: int) -> tuple[list[Any], list[str]]:
    """Return `count` header-only stand-in frames, in observation order, with their names."""
    frames = [SimpleNamespace(header={"MJD-OBS": 60000.0 + index}) for index in range(count)]
    return frames, [f"frame_{index}.fits" for index in range(count)]


def test_nod_stack_raises_a_clear_error_when_neither_position_found_a_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")

    with pytest.raises(Exception, match="No object trace was found in either nod position"):
        recipe.stack_extractions([None, None], orderJoins=None)

    assert any("No object trace was found" in str(message) for _, message in log.messages)


def test_offset_stack_raises_a_clear_error_when_the_offset_frame_found_no_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-offset")

    with pytest.raises(Exception, match="No object trace was found in the offset frame"):
        recipe.stack_extractions([None], orderJoins=None)


def test_stack_writes_nothing_when_no_trace_was_found(log: Any, tmp_path: Path, passthrough_snr_qcs: None) -> None:
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")
    before = sorted(tmp_path.rglob("*"))

    with pytest.raises(Exception, match="No object trace"):
        recipe.stack_extractions([None, None], orderJoins=None)

    assert sorted(tmp_path.rglob("*")) == before


def test_stack_still_stacks_the_spectrum_of_the_one_position_that_found_a_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    """Only the stacking of `[A_df, None]` is pinned here; the SNR and efficiency QCs of that case are DY-1616's."""
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")
    first, _ = _spectrum_frames()

    stacked, fitsPath = recipe.stack_extractions([first, None], orderJoins={"10-11": 1})

    assert len(stacked) == 3
    assert Path(fitsPath).is_file()


def test_nod_cycle_loop_raises_the_clear_error_when_every_cycle_found_no_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _cycle_loop_recipe(soxs_nod, log, tmp_path, "soxs-nod")
    allFrameA, namesA = _stub_frames(2)
    allFrameB, namesB = _stub_frames(2)

    with pytest.raises(Exception, match="No object trace was found in either nod position"):
        recipe._reduce_each_nodding_cycle(allFrameA, allFrameB, namesA, namesB, "orders.fits", False)


def test_offset_cycle_loop_raises_the_clear_error_when_every_cycle_found_no_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _cycle_loop_recipe(soxs_offset, log, tmp_path, "soxs-offset")
    allFrameON, namesON = _stub_frames(2)
    allFrameOFF, namesOFF = _stub_frames(2)

    with pytest.raises(Exception, match="No object trace was found in the offset frame"):
        recipe._reduce_each_offset_cycle(allFrameON, allFrameOFF, namesON, namesOFF, "orders.fits", False)


def test_nod_cycle_loop_keeps_the_spectra_of_cycles_that_found_a_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _cycle_loop_recipe(soxs_nod, log, tmp_path, "soxs-nod")
    first, second = _spectrum_frames()
    cycleResults = iter([(None, None, None), (first, second, {"10-11": 1})])
    recipe.process_single_ab_nodding_cycle = lambda **_: next(cycleResults)
    allFrameA, namesA = _stub_frames(2)
    allFrameB, namesB = _stub_frames(2)

    stacked, orderJoins = recipe._reduce_each_nodding_cycle(allFrameA, allFrameB, namesA, namesB, "orders.fits", False)

    assert len(stacked) == 3
    assert orderJoins == {"10-11": 1}


def _run_nod_cycles(recipe: Any, cycleResults: list[tuple[Any, Any, Any]]) -> tuple[Any, Any]:
    """Run the per-cycle nod loop over one frame pair per entry of `cycleResults`."""
    results = iter(cycleResults)
    recipe.process_single_ab_nodding_cycle = lambda **_: next(results)
    allFrameA, namesA = _stub_frames(len(cycleResults))
    allFrameB, namesB = _stub_frames(len(cycleResults))
    return recipe._reduce_each_nodding_cycle(allFrameA, allFrameB, namesA, namesB, "orders.fits", False)


def _run_offset_cycles(recipe: Any, cycleResults: list[tuple[Any, Any, Any]]) -> tuple[Any, Any, Any]:
    """Run the per-cycle offset loop over one frame pair per entry of `cycleResults`."""
    results = iter(cycleResults)
    recipe.process_single_ab_nodding_cycle = lambda **_: next(results)
    allFrameON, namesON = _stub_frames(len(cycleResults))
    allFrameOFF, namesOFF = _stub_frames(len(cycleResults))
    return recipe._reduce_each_offset_cycle(allFrameON, allFrameOFF, namesON, namesOFF, "orders.fits", False)


def test_nod_cycle_loop_keeps_the_order_joins_when_the_last_cycle_found_no_trace(
    log: Any, tmp_path: Path, snr_qc_order_joins: list[Any]
) -> None:
    recipe = _cycle_loop_recipe(soxs_nod, log, tmp_path, "soxs-nod")
    first, second = _spectrum_frames()

    stacked, orderJoins = _run_nod_cycles(recipe, [(first, second, {"10-11": 1}), (None, None, None)])

    assert len(stacked) == 3
    assert orderJoins == {"10-11": 1}
    assert snr_qc_order_joins == [{"10-11": 1}]


def test_offset_cycle_loop_keeps_the_order_joins_when_the_last_cycle_found_no_trace(
    log: Any, tmp_path: Path, snr_qc_order_joins: list[Any]
) -> None:
    recipe = _cycle_loop_recipe(soxs_offset, log, tmp_path, "soxs-offset")
    first, _ = _spectrum_frames()

    stacked, _, orderJoins = _run_offset_cycles(recipe, [(first, None, {"10-11": 1}), (None, None, None)])

    assert len(stacked) == 3
    assert orderJoins == {"10-11": 1}
    assert snr_qc_order_joins == [{"10-11": 1}]


def test_nod_cycle_loop_stacks_when_each_position_found_a_trace_in_a_different_cycle(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _cycle_loop_recipe(soxs_nod, log, tmp_path, "soxs-nod")
    first, second = _spectrum_frames()

    stacked, _ = _run_nod_cycles(recipe, [(None, second, {"10-11": 1}), (first, None, {"10-11": 1})])
    expected, _ = recipe.stack_extractions([second, first], orderJoins={"10-11": 1})
    firstOnly, _ = recipe.stack_extractions([first], orderJoins={"10-11": 1})

    assert stacked["FLUX_COUNTS"].tolist() == pytest.approx(expected["FLUX_COUNTS"].tolist())
    assert stacked["FLUX_COUNTS"].tolist() != pytest.approx(firstOnly["FLUX_COUNTS"].tolist())


def test_offset_cycle_loop_stacks_when_only_a_later_cycle_found_a_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _cycle_loop_recipe(soxs_offset, log, tmp_path, "soxs-offset")
    first, _ = _spectrum_frames()

    stacked, _, _ = _run_offset_cycles(recipe, [(None, None, None), (first, None, {"10-11": 1})])
    expected, _ = recipe.stack_extractions([first], orderJoins={"10-11": 1})

    assert stacked["FLUX_COUNTS"].tolist() == pytest.approx(expected["FLUX_COUNTS"].tolist())


def _extracting_nod_recipe(
    monkeypatch: pytest.MonkeyPatch, log: Any, tmp_path: Path, extractions: list[tuple[Any, Any]]
) -> Any:
    """Return a nod recipe whose real cycle code extracts `extractions`, one `(spectrum, orderJoins)` per call.

    The extractor stands in for `horne_extraction`. Call 1 is the A-B extraction and call 2 the B-A extraction of
    the cycle, and so on. A `(None, None)` entry is an extraction that found no object trace.
    """
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")
    recipe.debug = False
    recipe.recipeSettings = {"save_single_frame_extractions": False, "use_flat": False}
    recipe.twoDMap = "two-d-map.fits"
    recipe.dispMap = "dispersion-map.fits"
    recipe.startNightDate = "2024-01-02"
    recipe.turnOffMP = True
    # produce_product SETS THIS FROM THE RAW FRAMES; FALSE LETS THE EXTRACTOR USE ITS OWN FRAME HEADER
    recipe.scienceHeaders = False
    callCount = {"count": 0}

    class StubExtractor:
        flatPassHandoff = None
        noddingSequence = ""
        extractionNotes: tuple = ()

        def __init__(self, **kwargs: object) -> None:
            pass

        def extract(self) -> tuple[Any, Any, Any, Any, Any]:
            spectrum, joins = extractions[callCount["count"]]
            callCount["count"] += 1
            return recipe.qc, recipe.products, spectrum, joins, None if spectrum is None else "extraction.fits"

    nodModule = import_module("soxspipe.recipes.soxs_nod")
    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", StubExtractor)
    monkeypatch.setattr(nodModule, "generic_quality_checks", lambda **kwargs: kwargs["qcTable"])
    monkeypatch.setattr(nodModule, "spectroscopic_image_quality_checks", lambda **kwargs: kwargs["qcTable"])
    return recipe


def _run_real_nod_cycle(recipe: Any) -> tuple[Any, Any]:
    """Run the per-cycle nod loop over one real A and B frame pair, so the joins travel the whole path."""
    aFrame, bFrame = _ab_frames(aOverrides={"MJD-OBS": 60000.0}, bOverrides={"MJD-OBS": 60000.1})
    return recipe._reduce_each_nodding_cycle([aFrame], [bFrame], ["a.fits"], ["b.fits"], "orders.fits", False)


def test_nod_cycle_keeps_the_a_spectrum_and_joins_when_only_the_a_position_found_a_trace(
    monkeypatch: pytest.MonkeyPatch, log: Any, tmp_path: Path, snr_qc_order_joins: list[Any]
) -> None:
    first, _ = _spectrum_frames()
    recipe = _extracting_nod_recipe(monkeypatch, log, tmp_path, [(first, {"10-11": 1}), (None, None)])

    stacked, orderJoins = _run_real_nod_cycle(recipe)

    assert orderJoins == {"10-11": 1}
    assert snr_qc_order_joins == [{"10-11": 1}]
    expected, _ = recipe.stack_extractions([first], orderJoins={"10-11": 1})
    assert stacked["FLUX_COUNTS"].tolist() == pytest.approx(expected["FLUX_COUNTS"].tolist())


def test_nod_cycle_keeps_the_b_spectrum_and_joins_when_only_the_b_position_found_a_trace(
    monkeypatch: pytest.MonkeyPatch, log: Any, tmp_path: Path, snr_qc_order_joins: list[Any]
) -> None:
    second, _ = _spectrum_frames()
    recipe = _extracting_nod_recipe(monkeypatch, log, tmp_path, [(None, None), (second, {"20-21": 2})])

    stacked, orderJoins = _run_real_nod_cycle(recipe)

    assert orderJoins == {"20-21": 2}
    assert snr_qc_order_joins == [{"20-21": 2}]
    expected, _ = recipe.stack_extractions([second], orderJoins={"20-21": 2})
    assert stacked["FLUX_COUNTS"].tolist() == pytest.approx(expected["FLUX_COUNTS"].tolist())


def test_nod_cycle_uses_the_b_joins_when_both_positions_found_a_trace(
    monkeypatch: pytest.MonkeyPatch, log: Any, tmp_path: Path, snr_qc_order_joins: list[Any]
) -> None:
    first, second = _spectrum_frames()
    recipe = _extracting_nod_recipe(monkeypatch, log, tmp_path, [(first, {"10-11": 1}), (second, {"20-21": 2})])

    _, orderJoins = _run_real_nod_cycle(recipe)

    assert orderJoins == {"20-21": 2}
    assert snr_qc_order_joins == [{"20-21": 2}]


def test_nod_cycle_raises_when_neither_position_found_a_trace(
    monkeypatch: pytest.MonkeyPatch, log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _extracting_nod_recipe(monkeypatch, log, tmp_path, [(None, None), (None, None)])

    with pytest.raises(Exception, match="No object trace was found in either nod position"):
        _run_real_nod_cycle(recipe)


def test_nod_stack_stacks_the_one_combined_extraction_of_a_shift_and_add_pass(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    """A shift-and-add pass hands the stack one spectrum, whose median is itself (DY-1674)."""
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")
    first, _ = _spectrum_frames()

    stacked, fitsPath = recipe.stack_extractions([first], orderJoins={"10-11": 1})

    assert list(stacked["FLUX_COUNTS"]) == [9.0, 55.0, 42.0]
    assert Path(fitsPath).is_file()


def test_nod_stack_raises_a_clear_error_when_the_combined_extraction_found_no_trace(
    log: Any, tmp_path: Path, passthrough_snr_qcs: None
) -> None:
    recipe = _stack_recipe(log, tmp_path, recipeName="soxs-nod")

    with pytest.raises(Exception, match="No object trace was found in either nod position"):
        recipe.stack_extractions([None], orderJoins=None)


@pytest.mark.parametrize(
    ("spectrumB", "expectedLength"),
    [(False, 1), (None, 2), ("spectrum", 2)],
)
def test_only_a_pass_that_made_one_extraction_stacks_one_spectrum(spectrumB: object, expectedLength: int) -> None:
    """A B spectrum of `False` means no B extraction was made; `None` means one was made and found no trace."""
    spectra = soxs_nod._spectra_to_stack("spectrum A", spectrumB)

    assert len(spectra) == expectedLength
    assert spectra[0] == "spectrum A"
