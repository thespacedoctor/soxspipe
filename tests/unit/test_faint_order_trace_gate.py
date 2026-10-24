"""The faint-order trace runs in `detect_continuum.get()` only where it is switched on (DY-1941)."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.detect_continuum import detect_continuum
from tests.unit._continuum_vis_scene import (
    TraceFrame,
    vis_dispersion_map,
    vis_transfer_image,
    vis_transfer_table,
)

pytestmark = pytest.mark.unit
continuumModule = importlib.import_module("soxspipe.commonutils.detect_continuum")
PACKAGE_DIR = Path(soxspipe.__file__).parent
FAINT_QC = "N ORDERS FAINT TRACE"
AMPLITUDES = {1.0: 100.0, 2.0: 0.5, 3.0: 100.0, 4.0: 100.0}
FLANKS = {1.0: 30.0, 2.0: 0.8, 3.0: 30.0, 4.0: 30.0}


def _shipped_block(recipe: str) -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return dict(settings[recipe]["vis"]["detect-continuum"])


def _block(**overrides: object) -> dict:
    """Return a `detect-continuum` block with the faint trace switched on, then apply `overrides`."""
    return {
        "faint-order-trace": True,
        "poly-fitting-residual-clipping-sigma": 4,
        "poly-clipping-iteration-limit": 5,
        "max-unsampled-trace-fraction": 0.2,
    } | overrides


class _Run:
    """What `get()` handed to the first global fit, and the detector it ran on."""

    def __init__(self, detector: detect_continuum, sampled: pd.DataFrame) -> None:
        self.detector = detector
        self.sampled = sampled
        self.fitted: list[pd.DataFrame] = []


def _run_get(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    *,
    block: dict,
    sampled: pd.DataFrame,
    image: np.ndarray,
    inst: str = "SOXS",
    arm: str = "VIS",
    recipe: str = "soxs-stare",
    suppliedTable: bool = False,
) -> _Run:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.inst, detector.arm, detector.recipeName = inst, arm, recipe
    detector.axisA, detector.axisB = "x", "y"
    detector.orderDeg = detector.axisBDeg = 1
    detector.recipeSettings = block
    detector.settings = {"tune-pipeline": False}
    detector.dateObs = "2024-01-02T00:00:00"
    detector.qc = pd.DataFrame({"qc_name": []})
    detector.products = False
    detector.noddingSequence = ""
    detector.notFlattened = ""
    detector.coeff_dict = {}
    detector.dispersion_map = vis_dispersion_map(tmp_path)
    detector.mapBinning = {"x": 1.0, "y": 1.0}
    detector.traceFrame = np.ma.masked_array(image, mask=np.zeros(image.shape, dtype=bool)).view(TraceFrame)
    detector.orderPixelTable = sampled if suppliedTable else False
    run = _Run(detector, sampled)

    def fit_trace(**kwargs: object) -> tuple[list[float], pd.DataFrame, pd.DataFrame]:
        run.fitted.append(kwargs["pixelList"])
        fitted = kwargs["pixelList"].copy()
        fitted[f"{kwargs['axisACol']}_fit_res"] = 0.0
        return [1.0, 2.0, 3.0, 4.0], fitted, pd.DataFrame()

    monkeypatch.setattr(detector, "sample_trace", lambda: (sampled, 100.0))
    monkeypatch.setattr(detector, "fit_global_polynomial", fit_trace)
    monkeypatch.setattr(detector, "plot_results", lambda **kwargs: "/qc/continuum.pdf")
    monkeypatch.setattr(detector, "build_order_metadata", lambda **kwargs: pd.DataFrame({"axis": ["x"]}))
    monkeypatch.setattr(detector, "write_order_table_to_file", lambda **kwargs: "/products/orders.fits")
    detector.get()
    return run


def _scene(faintOrders: tuple[float, ...]) -> tuple[pd.DataFrame, np.ndarray]:
    amplitudes = {order: (0.5 if order in faintOrders else 100.0) for order in AMPLITUDES}
    flanks = {order: (0.8 if order in faintOrders else 30.0) for order in AMPLITUDES}
    return vis_transfer_table(faintOrders).assign(mask=False), vis_transfer_image(amplitudes, flanks)


def _assert_untouched(run: _Run) -> None:
    firstFit = run.fitted[0]
    assert firstFit is run.sampled
    assert "faint_transfer" not in firstFit.columns
    assert FAINT_QC not in run.detector.qc["qc_name"].tolist()


def test_get_transfers_the_faint_order_when_the_trace_is_switched_on(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene((2.0,))
    original = sampled.copy()

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT: ONLY THE FAINT ORDER'S ROWS ARE MARKED AND CHANGED
    fitted = run.fitted[0]
    assert fitted["faint_transfer"].tolist() == (original["order"] == 2.0).tolist()
    brightRows = original["order"] != 2.0
    pd.testing.assert_frame_equal(fitted.loc[brightRows, original.columns], original.loc[brightRows])
    assert not fitted.loc[~brightRows, "cont_x"].equals(original.loc[~brightRows, "cont_x"])


def test_get_leaves_a_bright_only_table_as_the_same_object_with_no_new_column_or_qc(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene(())
    before = sampled.copy()

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    _assert_untouched(run)
    pd.testing.assert_frame_equal(run.fitted[0][before.columns], before)


GATE_OFF_CASES = [
    pytest.param(_block(), {"recipe": "soxs-order-centres"}, id="order-centres recipe"),
    pytest.param(_block(), {"arm": "NIR"}, id="NIR arm"),
    pytest.param(_block(), {"inst": "XSH"}, id="X-Shooter"),
    pytest.param({key: value for key, value in _block().items() if key != "faint-order-trace"}, {}, id="key missing"),
    pytest.param(_block(**{"faint-order-trace": False}), {}, id="switched off"),
    pytest.param(_shipped_block("soxs-stare-std"), {"recipe": "soxs-stare-std"}, id="stare-std block"),
    pytest.param(_shipped_block("soxs-nod-std"), {"recipe": "soxs-nod-std"}, id="nod-std block"),
    pytest.param(_shipped_block("soxs-offset-std"), {"recipe": "soxs-offset-std"}, id="offset-std block"),
    pytest.param(_block(), {"suppliedTable": True}, id="table supplied for tuning"),
]


@pytest.mark.parametrize(("block", "options"), GATE_OFF_CASES)
def test_get_changes_nothing_where_the_trace_is_off(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, block: dict, options: dict
) -> None:
    # ARRANGE: A FRAME WITH A FAINT u ORDER, WHICH WOULD BE TRANSFERRED IF THE TRACE WERE ON
    sampled, image = _scene((2.0,))
    before = sampled.copy()

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=block, sampled=sampled, image=image, **options)

    # ASSERT: THE SAME TABLE OBJECT REACHES THE FIT, WITH EVERY SAMPLED VALUE AS IT WAS
    _assert_untouched(run)
    pd.testing.assert_frame_equal(run.fitted[0][before.columns], before)


def test_the_trace_never_runs_for_a_gate_that_is_off(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene((2.0,))

    def forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("the faint-order transfer ran with the trace switched off")

    monkeypatch.setattr(continuumModule, "transfer_faint_orders", forbidden)

    # ACT
    block = _block(**{"faint-order-trace": False})
    run = _run_get(log, monkeypatch, tmp_path, block=block, sampled=sampled, image=image)

    # ASSERT
    _assert_untouched(run)


def test_sample_trace_keeps_the_frame_to_map_binning_for_the_faint_order_pixels(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE FRAME IS BINNED 2 IN x AND 3 IN y RELATIVE TO THE MAP
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.inst, detector.arm, detector.recipeName = "XSH", "VIS", "soxs-stare"
    detector.axisA, detector.axisB = "x", "y"
    detector.traceFrame = np.ma.masked_array(np.ones((32, 32)), mask=np.zeros((32, 32), dtype=bool)).view(TraceFrame)
    detector.traceFrame.header = {"BINX": 2, "BINY": 3}
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.recipeSettings = {"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5}
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    samples = pd.DataFrame(
        {"order": 10.0, "wavelength": np.linspace(500.0, 510.0, 20), "fit_x": 16.0, "fit_y": np.arange(4.0, 24.0)}
    )
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (samples.copy(), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    monkeypatch.setattr(
        detector,
        "fit_1d_gaussian_to_slices",
        lambda orderPixelTable, **_: orderPixelTable.assign(
            cont_x=orderPixelTable["fit_x"] - 2.0, cont_y=orderPixelTable["fit_y"], gauss_stddev=1.5
        ),
    )

    # ACT
    detector.sample_trace()

    # ASSERT
    assert detector.mapBinning == {"x": 2.0, "y": 3.0}


class _NoteLog:
    """A recording logger that also takes extraction notes, as the Horne extraction's logger does."""

    def __init__(self, log: object) -> None:
        self.log = log
        self.notes: list[str] = []

    def note(self, message: str) -> None:
        self.notes.append(message)

    def __getattr__(self, name: str) -> object:
        return getattr(self.log, name)


def test_get_records_the_count_of_transferred_orders_as_a_header_qc(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene((2.0,))

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    qc = run.detector.qc
    row = qc.loc[qc["qc_name"] == FAINT_QC].iloc[0]
    assert row["qc_value"] == 1
    assert bool(row["to_header"]) is True
    assert row["soxspipe_recipe"] == "soxs-stare"


def test_get_adds_the_qc_only_when_an_order_was_transferred(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE: U AND G ARE BOTH FAINT, SO NEITHER HAS A BRIGHT PARTNER
    sampled, image = _scene((2.0, 3.0))

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    _assert_untouched(run)


def test_get_prints_and_notes_one_line_for_each_transferred_order(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene((2.0, 4.0))
    noteLog = _NoteLog(log)

    # ACT
    _run_get(noteLog, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    assert len(noteLog.notes) == 2
    assert noteLog.notes[0].startswith("Order 2 (u) trace transferred from 3 (g): ")
    assert noteLog.notes[1].startswith("Order 4 (i) trace transferred from 1 (r): ")
    assert "/15 bins, deg " in noteLog.notes[0]
    printed = [message for level, message in log.messages if level == "print" and "trace transferred" in message]
    assert printed == noteLog.notes


def test_get_runs_without_a_note_method_on_the_logger(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    sampled, image = _scene((2.0,))

    # ACT
    run = _run_get(log, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    assert FAINT_QC in run.detector.qc["qc_name"].tolist()
    assert any("Order 2 (u) trace transferred from 3 (g)" in message for _, message in log.messages)


def test_get_says_so_when_a_faint_order_falls_back_to_the_anchor_track(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE: NOTHING AT ALL IS IN THE u ORDER
    sampled, _ = _scene((2.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0})
    noteLog = _NoteLog(log)

    # ACT
    _run_get(noteLog, monkeypatch, tmp_path, block=_block(), sampled=sampled, image=image)

    # ASSERT
    assert noteLog.notes == ["Order 2 (u) trace transferred from 3 (g): no bin kept, anchor track used"]
