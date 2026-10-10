"""Continuum-trace fitting never writes into a pandas slice of a caller's table (DY-1294).

Every test runs under ``mode.chained_assignment = "raise"`` and keeps the parent table alive, because pandas only
flags a write to a slice while the frame it was cut from still exists.
"""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_continuum import detect_continuum

pytestmark = pytest.mark.unit
continuumModule = importlib.import_module("soxspipe.commonutils.detect_continuum")

ROWS_PER_ORDER = 60
TRACE_STDDEV = 1.5


def _detector(log: object, *, orderDeg: int = 1, axisBDeg: int = 1) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = orderDeg
    detector.axisBDeg = axisBDeg
    detector.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 3.0,
        "poly-clipping-iteration-limit": 4,
    }
    return detector


def _two_order_trace() -> pd.DataFrame:
    """Two orders with x = 1 + 2y + 0.5 order, one NaN row and one outlier."""
    orders = np.repeat([10.0, 11.0], 20)
    yValues = np.tile(np.arange(20, dtype=float), 2)
    table = pd.DataFrame(
        {
            "source": np.arange(40),
            "order": orders,
            "cont_y": yValues,
            "cont_x": 1.0 + 2.0 * yValues + 0.5 * orders,
        }
    )
    table.loc[5, "cont_x"] = 500.0
    table.loc[25, "cont_y"] = np.nan
    return table


def test_global_fit_of_a_filtered_trace_table_never_writes_into_the_caller_slice(log: object) -> None:
    # ARRANGE: THE CALLER PASSES A FILTERED VIEW AND KEEPS THE PARENT ALIVE, AS THE ORDER-EDGE CODE DOES
    detector = _detector(log)
    parent = _two_order_trace()
    parent["pre-clipped"] = False
    caller = parent.loc[parent["source"] != 39]

    # ACT
    with pd.option_context("mode.chained_assignment", "raise"):
        coefficients, fitted, clipped = detector.fit_global_polynomial(caller)

    # ASSERT: x = 1 + 0.5 ORDER + 2 y, SO TERMS (ORDER^0 y^0, ORDER^0 y^1, ORDER^1 y^0, ORDER^1 y^1)
    np.testing.assert_allclose(coefficients, [1.0, 2.0, 0.5, 0.0], rtol=0, atol=1e-8)
    assert set(clipped["source"]) == {5, 25}
    assert {5, 25, 39}.isdisjoint(set(fitted["source"]))
    np.testing.assert_allclose(fitted["cont_x_fit_res"], 0.0, rtol=0, atol=1e-8)


def test_global_fit_of_a_dropna_edge_table_never_writes_into_the_caller_slice(log: object) -> None:
    # ARRANGE: THE ORDER-EDGE DETECTOR PASSES table.dropna(...) WITH NO cont_x OR pre-clipped COLUMN
    detector = _detector(log)
    parent = _two_order_trace().rename(columns={"cont_x": "xcoord_upper", "cont_y": "ycoord"})
    parent.loc[5, "xcoord_upper"] = 1.0 + 2.0 * parent.loc[5, "ycoord"] + 0.5 * parent.loc[5, "order"]
    edges = parent.dropna(axis="index", how="any", subset=["ycoord"])

    # ACT
    with pd.option_context("mode.chained_assignment", "raise"):
        coefficients, fitted, clipped = detector.fit_global_polynomial(
            edges, axisACol="xcoord_upper", axisBCol="ycoord"
        )

    # ASSERT
    np.testing.assert_allclose(coefficients, [1.0, 2.0, 0.5, 0.0], rtol=0, atol=1e-8)
    assert len(fitted.index) == 39
    np.testing.assert_allclose(fitted["xcoord_upper_fit"], fitted["xcoord_upper"], rtol=0, atol=1e-8)


def _shift_for(yValues: np.ndarray) -> np.ndarray:
    """The object sits 1.9 to 2.1 pixels below the dispersion-solution centre, varying smoothly along the order."""
    return 1.9 + 0.2 * (yValues % ROWS_PER_ORDER) / (ROWS_PER_ORDER - 1)


def _stddev_for(yValues: np.ndarray) -> np.ndarray:
    """The trace width varies smoothly between 1.4 and 1.6 pixels, so sigma-clipping sees a real spread."""
    return TRACE_STDDEV - 0.1 + 0.2 * ((7 * yValues) % ROWS_PER_ORDER) / (ROWS_PER_ORDER - 1)


def _gaussian_slice(**kwargs: object) -> tuple[np.ndarray, int, float]:
    """A cross-dispersion slice holding a Gaussian at the true trace position for this sample."""
    x = float(kwargs["x"])
    y = float(kwargs["y"])
    length = int(kwargs["length"])
    lengthOffset = int(x - length / 2)
    trueCentre = 16.0 - _shift_for(np.array([y]))[0]
    stddev = _stddev_for(np.array([y]))[0]
    pixels = np.arange(length, dtype=float)
    data = 100.0 * np.exp(-0.5 * ((pixels - (trueCentre - lengthOffset)) / stddev) ** 2)
    return data, lengthOffset, y


def _trace_sampler(log: object, monkeypatch: pytest.MonkeyPatch, *, inst: str, orders: list[float]) -> detect_continuum:
    detector = _detector(log)
    detector.inst = inst
    detector.recipeName = "soxs-stare"
    detector.traceFrame = CCDData(np.ones((32, 32)), unit=u.electron)
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.recipeSettings.update({"slice-length": 12, "peak-sigma-limit": 3.0, "slice-width": 5})
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    detector.sliceWidth = 5
    detector.detectorParams = {"dispersion-axis": "y"}
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    rowCount = ROWS_PER_ORDER * len(orders)
    samples = pd.DataFrame(
        {
            "order": np.repeat(orders, ROWS_PER_ORDER),
            "wavelength": np.linspace(500.0, 700.0, rowCount),
            "fit_x": np.full(rowCount, 16.0),
            "fit_y": np.arange(rowCount, dtype=float),
        }
    )
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (samples.copy(), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    monkeypatch.setattr(continuumModule, "cut_image_slice", _gaussian_slice)
    return detector


def test_soxs_vis_trace_sampling_never_writes_into_the_strided_probe_slices(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: THE VIS SAMPLER PROBES EVERY 10TH ROW OF EACH ORDER PAIR BEFORE THE FULL FIT
    detector = _trace_sampler(log, monkeypatch, inst="SOXS", orders=[1.0, 2.0, 3.0, 4.0])

    # ACT
    with pd.option_context("mode.chained_assignment", "raise"):
        result, detectionPercentage = detector.sample_trace()

    # ASSERT: EVERY SAMPLE IS FOUND AT ITS TRUE POSITION
    assert detectionPercentage == pytest.approx(100.0)
    assert result.groupby("order").size().to_dict() == {1.0: 60, 2.0: 60, 3.0: 60, 4.0: 60}
    expected = 16.0 - _shift_for(result["fit_y"].to_numpy())
    np.testing.assert_allclose(result["cont_x"], expected, rtol=0, atol=1e-3)
    np.testing.assert_allclose(result["gauss_stddev"], _stddev_for(result["fit_y"].to_numpy()), rtol=0, atol=1e-3)
    assert result["pre-clipped"].eq(False).all()


def test_trace_sampling_never_writes_into_the_strided_probe_slices_outside_soxs_vis(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    detector = _trace_sampler(log, monkeypatch, inst="XSH", orders=[10.0, 11.0])

    # ACT
    with pd.option_context("mode.chained_assignment", "raise"):
        result, detectionPercentage = detector.sample_trace()

    # ASSERT
    assert detectionPercentage == pytest.approx(100.0)
    expected = 16.0 - _shift_for(result["fit_y"].to_numpy())
    np.testing.assert_allclose(result["cont_x"], expected, rtol=0, atol=1e-3)
    assert result["pre-clipped"].eq(False).all()
