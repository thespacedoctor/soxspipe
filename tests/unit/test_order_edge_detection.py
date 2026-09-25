"""Controlled cross-section contracts for order-edge detection."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_order_edges import detect_order_edges

pytestmark = pytest.mark.unit
edgeModule = importlib.import_module("soxspipe.commonutils.detect_order_edges")


def _detector(log: object) -> detect_order_edges:
    detector = object.__new__(detect_order_edges)
    detector.log = log
    detector.axisA = "x"
    detector.axisB = "y"
    detector.arm = "VIS"
    detector.flatFrame = object()
    detector.sliceWidth = 3
    detector.sliceLength = 41
    detector.minThresholdPercenage = 0.2
    detector.maxThresholdPercenage = 0.8
    detector.debug = False
    return detector


def _mock_slice(
    monkeypatch: pytest.MonkeyPatch,
    profile: np.ndarray,
) -> None:
    monkeypatch.setattr(
        edgeModule,
        "cut_image_slice",
        lambda **kwargs: (profile.copy(), 0, 0),
    )


def test_determine_order_flux_threshold_uses_smoothed_central_slice(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    profile = np.concatenate([np.zeros(10), np.ones(21), np.zeros(10)])
    _mock_slice(monkeypatch, profile)
    orderData = pd.Series({"order": 10.0})
    pixels = pd.DataFrame(
        {
            "order": [10.0, 10.0, 10.0],
            "xcoord_centre": [30.0, 30.0, 30.0],
            "ycoord": [10.0, 20.0, 30.0],
        }
    )

    result = detector.determine_order_flux_threshold(orderData.copy(), pixels)

    assert result["minThreshold"] == pytest.approx(0.2)
    assert result["maxThreshold"] == pytest.approx(0.8)
    assert result["maxvalue"] == pytest.approx(1.0)


def test_edge_positions_are_symmetric_for_symmetric_profile(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    profile = np.concatenate([np.zeros(10), np.ones(21), np.zeros(10)])
    _mock_slice(monkeypatch, profile)
    orderData = pd.Series({"order": 10.0, "xcoord_centre": 30.0, "ycoord": 20.0})

    result = detector.determine_lower_upper_edge_pixel_positions(orderData.copy())

    assert result["xcoord_lower"] == pytest.approx(19.2)
    assert result["xcoord_upper"] == pytest.approx(38.8)
    assert (result["xcoord_lower"] + result["xcoord_upper"]) / 2 == pytest.approx(29.0)


def test_edge_positions_remain_nan_for_invalid_flat_profile(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    _mock_slice(monkeypatch, np.ones(41))
    orderData = pd.Series({"order": 10.0, "xcoord_centre": 30.0, "ycoord": 20.0})

    result = detector.determine_lower_upper_edge_pixel_positions(orderData.copy())

    assert np.isnan(result["xcoord_lower"])
    assert np.isnan(result["xcoord_upper"])


def test_edge_positions_map_to_y_for_horizontal_dispersion(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    detector = _detector(log)
    detector.axisA = "y"
    detector.axisB = "x"
    profile = np.concatenate([np.zeros(10), np.ones(21), np.zeros(10)])
    _mock_slice(monkeypatch, profile)
    orderData = pd.Series({"order": 10.0, "xcoord": 20.0, "ycoord_centre": 30.0})

    result = detector.determine_lower_upper_edge_pixel_positions(orderData.copy())

    assert result["ycoord_lower"] == pytest.approx(19.2)
    assert result["ycoord_upper"] == pytest.approx(38.8)
