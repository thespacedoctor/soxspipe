"""Deterministic least-squares contracts for the global order x axis-B polynomial fit (DY-1222)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import LSTSQ_RCOND, detect_continuum

pytestmark = pytest.mark.unit

NO_CLIPPING_SIGMA = 1.0e6


def _detector(log: object, *, orderDeg: int, axisBDeg: int, clippingSigma: float = 2.0) -> detect_continuum:
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.arm = "VIS"
    detector.axisA = "x"
    detector.axisB = "y"
    detector.orderDeg = orderDeg
    detector.axisBDeg = axisBDeg
    detector.recipeSettings = {
        "poly-fitting-residual-clipping-sigma": clippingSigma,
        "poly-clipping-iteration-limit": 4,
    }
    return detector


def _design_matrix(pixels: pd.DataFrame, orderDeg: int, axisBDeg: int) -> np.ndarray:
    """Independent monomial design matrix, term order i (order) outer, j (axis B) inner."""
    orders = pixels["order"].to_numpy(dtype=float)
    yValues = pixels["cont_y"].to_numpy(dtype=float)
    return np.column_stack([orders**i * yValues**j for i in range(orderDeg + 1) for j in range(axisBDeg + 1)])


def _rank_deficient_pixels() -> pd.DataFrame:
    """Three distinct orders cannot support an order-degree-3 polynomial, so the 16-term design has rank 12."""
    rng = np.random.default_rng(1222)
    orders = np.repeat([10.0, 11.0, 12.0], 40)
    yValues = np.tile(np.linspace(0.0, 2000.0, 40), 3)
    xValues = 5.0 + 0.01 * yValues + 0.3 * orders + rng.normal(0.0, 0.05, orders.size)
    return pd.DataFrame({"order": orders, "cont_y": yValues, "cont_x": xValues})


def _minimum_norm_solution(pixels: pd.DataFrame, orderDeg: int, axisBDeg: int) -> np.ndarray:
    design = _design_matrix(pixels, orderDeg, axisBDeg)
    scale = np.linalg.norm(design, axis=0)
    solution = np.linalg.pinv(design / scale, rcond=LSTSQ_RCOND) @ pixels["cont_x"].to_numpy()
    return solution / scale


def test_rank_deficient_fit_returns_the_minimum_norm_solution(log: object) -> None:
    detector = _detector(log, orderDeg=3, axisBDeg=3, clippingSigma=NO_CLIPPING_SIGMA)
    pixels = _rank_deficient_pixels()
    design = _design_matrix(pixels, 3, 3)
    assert np.linalg.matrix_rank(design / np.linalg.norm(design, axis=0), tol=LSTSQ_RCOND) == 12

    coefficients, _, _ = detector.fit_global_polynomial(pixels.copy())

    expected = _minimum_norm_solution(pixels, 3, 3)
    np.testing.assert_allclose(design @ coefficients, design @ expected, rtol=0, atol=1e-6)
    np.testing.assert_allclose(coefficients, expected, rtol=1e-6, atol=1e-12)


def test_rank_deficient_fit_is_bit_identical_across_calls_and_allocations(log: object) -> None:
    detector = _detector(log, orderDeg=3, axisBDeg=3, clippingSigma=NO_CLIPPING_SIGMA)
    pixels = _rank_deficient_pixels()

    results = []
    ballast = []
    for padding in range(12):
        ballast.append(np.zeros(padding * 7919 + 1))
        coefficients, _, _ = detector.fit_global_polynomial(pixels.copy())
        results.append(np.asarray(coefficients).tobytes())

    assert len(set(results)) == 1


def test_full_rank_data_recovers_the_generating_coefficients(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=2, clippingSigma=NO_CLIPPING_SIGMA)
    orders = np.repeat([10.0, 11.0, 12.0, 13.0], 25)
    yValues = np.tile(np.linspace(0.0, 100.0, 25), 4)
    trueCoefficients = np.array([1.5, 0.25, -0.002, 0.3, -0.01, 0.0004])
    pixels = pd.DataFrame({"order": orders, "cont_y": yValues})
    pixels["cont_x"] = _design_matrix(pixels, 1, 2) @ trueCoefficients

    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    np.testing.assert_allclose(coefficients, trueCoefficients, rtol=1e-8, atol=1e-10)
    np.testing.assert_allclose(fitted["cont_x_fit"], pixels["cont_x"], rtol=1e-10, atol=1e-8)
    assert len(clipped) == 0


def test_sigma_clipping_removes_injected_outliers(log: object) -> None:
    detector = _detector(log, orderDeg=0, axisBDeg=1)
    yValues = np.arange(60, dtype=float)
    pixels = pd.DataFrame(
        {
            "source": np.arange(60),
            "order": np.full(60, 10.0),
            "cont_y": yValues,
            "cont_x": 1.0 + 2.0 * yValues,
        }
    )
    outliers = [7, 31, 52]
    pixels.loc[outliers, "cont_x"] += [60.0, -80.0, 70.0]

    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    np.testing.assert_allclose(coefficients, [1.0, 2.0], rtol=1e-8, atol=1e-8)
    assert set(outliers).issubset(set(clipped["source"]))
    assert set(outliers).isdisjoint(set(fitted["source"]))
    np.testing.assert_allclose(fitted["cont_x_fit_res"], 0.0, rtol=0, atol=1e-8)


def test_fewer_points_than_coefficients_returns_no_solution(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=1)
    pixels = pd.DataFrame(
        {
            "order": [10.0, 10.0, 11.0],
            "cont_y": [0.0, 1.0, 2.0],
            "cont_x": [1.0, 2.0, 3.0],
        }
    )

    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    assert coefficients is None
    pd.testing.assert_frame_equal(fitted[list(pixels.columns)], pixels)
    assert clipped is fitted


def test_exponents_included_path_matches_the_computed_path(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=2, clippingSigma=NO_CLIPPING_SIGMA)
    pixels = pd.DataFrame(
        {
            "order": np.repeat([10.0, 11.0, 12.0], 20),
            "cont_y": np.tile(np.linspace(0.0, 50.0, 20), 3),
        }
    )
    pixels["cont_x"] = 2.0 + 0.1 * pixels["cont_y"] + 0.5 * pixels["order"] + 0.01 * pixels["order"] * pixels["cont_y"]
    withPowers = pixels.copy()
    for j in range(3):
        withPowers[f"y_pow_{j}"] = withPowers["cont_y"].pow(j)
    for i in range(2):
        withPowers[f"order_pow_{i}"] = withPowers["order"].pow(i)

    computed, _, _ = detector.fit_global_polynomial(pixels.copy())
    precomputed, _, _ = detector.fit_global_polynomial(withPowers, exponentsIncluded=True)

    np.testing.assert_allclose(precomputed, computed, rtol=1e-9, atol=1e-12)


def test_non_finite_values_raise_instead_of_returning_garbage(log: object) -> None:
    detector = _detector(log, orderDeg=0, axisBDeg=1)
    pixels = pd.DataFrame(
        {
            "order": np.full(6, 10.0),
            "cont_y": np.arange(6, dtype=float),
            "stddev": [1.0, 1.1, np.nan, 1.2, 1.3, 1.4],
        }
    )

    with pytest.raises(ValueError, match="infs or NaNs"):
        detector.fit_global_polynomial(pixels, axisACol="stddev")
