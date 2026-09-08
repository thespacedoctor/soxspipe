"""Analytic continuum-trace fitting contracts."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum

pytestmark = pytest.mark.unit
continuumModule = importlib.import_module("soxspipe.commonutils.detect_continuum")


def _detector(log: object, *, orderDeg: int = 0, axisBDeg: int = 1) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = orderDeg
    detector.axisBDeg = axisBDeg
    detector.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": 2.0,
        "poly-clipping-iteration-limit": 4,
    }
    return detector


def test_calculate_residuals_matches_exact_global_trace(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=1)
    pixels = pd.DataFrame(
        {
            "order": [10.0, 10.0, 11.0],
            "cont_y": [0.0, 1.0, 2.0],
        }
    )
    expectedFit = (
        1.0
        + 2.0 * pixels["cont_y"]
        + 3.0 * pixels["order"]
        + 4.0 * pixels["order"] * pixels["cont_y"]
    )
    pixels["cont_x"] = expectedFit - np.array([1.0, -2.0, 0.0])

    residuals, mean, std, median, fitted = detector.calculate_residuals(
        orderPixelTable=pixels,
        coeff=[1.0, 2.0, 3.0, 4.0],
        axisACol="cont_x",
        axisBCol="cont_y",
        orderCol="order",
    )

    np.testing.assert_allclose(fitted, expectedFit, rtol=1e-12, atol=1e-12)
    np.testing.assert_allclose(residuals, [1.0, -2.0, 0.0], rtol=1e-12, atol=1e-12)
    assert mean == pytest.approx(1.0)
    assert std == pytest.approx(np.std([1.0, -2.0, 0.0]))
    assert median == pytest.approx(1.0)


def test_global_fit_removes_nan_and_outlier_rows(log: object) -> None:
    detector = _detector(log)
    yValues = np.arange(12, dtype=float)
    pixels = pd.DataFrame(
        {
            "source": np.arange(12),
            "order": np.full(12, 10.0),
            "cont_y": yValues,
            "cont_x": 1.0 + 2.0 * yValues,
        }
    )
    pixels.loc[10, "cont_x"] = 100.0
    pixels.loc[11, "cont_y"] = np.nan

    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    np.testing.assert_allclose(coefficients, [1.0, 2.0], rtol=1e-6, atol=1e-6)
    assert set(fitted["source"]) == set(range(10))
    assert set(clipped["source"]) == {10, 11}
    np.testing.assert_allclose(fitted["cont_x_fit_res"], 0.0, rtol=0, atol=1e-6)


def test_global_fit_preserves_empty_input_failure(log: object) -> None:
    detector = _detector(log)
    empty = pd.DataFrame(columns=["order", "cont_y", "cont_x"])

    with pytest.raises(ValueError, match="cannot set a frame with no defined index"):
        detector.fit_global_polynomial(empty)


def test_gaussian_slice_fit_recovers_trace_and_leaves_low_signal_unfitted(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    detector.traceFrame = object()
    detector.sliceWidth = 3
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.debug = False
    detector.peakSigmaLimit = 3.0
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    samplePositions = np.arange(21, dtype=float)
    gaussian = 20.0 * np.exp(-0.5 * ((samplePositions - 10.0) / 2.0) ** 2)
    slices = iter([(gaussian, 5, 7.0), (np.zeros(21), 5, 8.0)])
    monkeypatch.setattr(
        continuumModule,
        "cut_image_slice",
        lambda **kwargs: next(slices),
    )
    pixels = pd.DataFrame({"fit_x": [15.0, 16.0], "fit_y": [7.0, 8.0]})

    result = detector.fit_1d_gaussian_to_slices(pixels.copy(), sliceLength=21)

    np.testing.assert_allclose(result.loc[0, "cont_x"], 15.0, rtol=0, atol=1e-6)
    np.testing.assert_allclose(result.loc[0, "cont_y"], 7.0, rtol=0, atol=1e-12)
    np.testing.assert_allclose(result.loc[0, "gauss_mean"], 10.0, rtol=0, atol=1e-6)
    np.testing.assert_allclose(result.loc[0, "gauss_stddev"], 2.0, rtol=0, atol=1e-6)
    assert np.isnan(result.loc[1, "cont_x"])
    assert np.isnan(result.loc[1, "gauss_mean"])
