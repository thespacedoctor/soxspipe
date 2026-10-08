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
    singularValues = np.linalg.svd(design / np.linalg.norm(design, axis=0), compute_uv=False)
    assert np.count_nonzero(singularValues > LSTSQ_RCOND * singularValues[0]) == 12

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


def test_exponents_included_flag_makes_the_fit_use_the_precomputed_power_columns(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=2, clippingSigma=NO_CLIPPING_SIGMA)
    pixels = pd.DataFrame(
        {
            "order": np.repeat([10.0, 11.0, 12.0], 20),
            "cont_y": np.tile(np.linspace(0.0, 50.0, 20), 3),
        }
    )
    pixels["cont_x"] = 2.0 + 0.1 * pixels["cont_y"] + 0.5 * pixels["order"] + 0.01 * pixels["order"] * pixels["cont_y"]
    # POWER COLUMNS BUILT FROM A SHIFTED AND SCALED BASIS, SO THEY DIFFER FROM A RECOMPUTATION FROM THE RAW COORDINATES
    shiftedY = pixels["cont_y"] / 10.0
    shiftedOrder = pixels["order"] - 10.0
    withPowers = pixels.copy()
    for j in range(3):
        withPowers[f"y_pow_{j}"] = shiftedY.pow(j)
    for i in range(2):
        withPowers[f"order_pow_{i}"] = shiftedOrder.pow(i)
    shiftedDesign = np.column_stack([shiftedOrder**i * shiftedY**j for i in range(2) for j in range(3)])
    expectedShifted = np.linalg.lstsq(shiftedDesign, pixels["cont_x"].to_numpy(), rcond=None)[0]

    computed, _, _ = detector.fit_global_polynomial(withPowers.copy())
    precomputed, _, _ = detector.fit_global_polynomial(withPowers.copy(), exponentsIncluded=True)

    np.testing.assert_allclose(precomputed, expectedShifted, rtol=1e-8, atol=1e-10)
    assert not np.allclose(precomputed, computed, rtol=1e-3, atol=1e-6)


def test_all_zero_design_columns_get_zero_coefficients_and_the_fit_still_matches(log: object) -> None:
    detector = _detector(log, orderDeg=1, axisBDeg=2, clippingSigma=NO_CLIPPING_SIGMA)
    pixels = pd.DataFrame(
        {
            "order": np.repeat([10.0, 11.0, 12.0], 4),
            "cont_y": np.zeros(12),
        }
    )
    pixels["cont_x"] = 2.0 + 0.5 * pixels["order"]

    coefficients, fitted, _ = detector.fit_global_polynomial(pixels.copy())

    # TERM ORDER IS ORDER-POWER OUTER, AXIS-B POWER INNER, SO EVERY AXIS-B POWER ABOVE ZERO IS AN ALL-ZERO COLUMN
    zeroColumns = [1, 2, 4, 5]
    assert np.isfinite(coefficients).all()
    np.testing.assert_array_equal(np.asarray(coefficients)[zeroColumns], 0.0)
    np.testing.assert_allclose(fitted["cont_x_fit"], pixels["cont_x"], rtol=0, atol=1e-9)


def test_nullable_dtype_missing_values_raise_instead_of_returning_garbage(log: object) -> None:
    detector = _detector(log, orderDeg=0, axisBDeg=1)
    pixels = pd.DataFrame(
        {
            "order": np.full(6, 10.0),
            "cont_y": np.arange(6, dtype=float),
            "stddev": pd.array([1.0, 1.1, pd.NA, 1.2, 1.3, 1.4], dtype="Float64"),
        }
    )

    with pytest.raises(ValueError, match="infs or NaNs"):
        detector.fit_global_polynomial(pixels, axisACol="stddev")


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


def _line_with_outliers() -> tuple[pd.DataFrame, list[int]]:
    """A noise-free straight line in one order with three gross outliers."""
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
    return pixels, outliers


def test_cap_exit_refits_the_rows_that_survive_the_last_clip_and_logs_it(log: object) -> None:
    # ARRANGE: ONE PASS ONLY, SO THE LOOP EXITS ON THE CAP WITH THE OUTLIERS JUST CLIPPED
    detector = _detector(log, orderDeg=0, axisBDeg=1)
    detector.recipeSettings["poly-clipping-iteration-limit"] = 1
    pixels, outliers = _line_with_outliers()

    # ACT
    coefficients, fitted, clipped = detector.fit_global_polynomial(pixels.copy())

    # ASSERT: THE COEFFICIENTS DESCRIBE THE SURVIVORS, THEIR RESIDUALS ARE CONSISTENT WITH THEM, AND THE CAP IS LOGGED
    assert set(outliers).issubset(set(clipped["source"]))
    np.testing.assert_allclose(coefficients, [1.0, 2.0], rtol=1e-8, atol=1e-8)
    np.testing.assert_allclose(fitted["cont_x_fit_res"], 0.0, rtol=0, atol=1e-8)
    capMessages = [message for level, message in log.messages if level == "info" and "iteration limit" in message]
    assert len(capMessages) == 1
    assert f"{len(clipped)} rows still being clipped" in capMessages[0]


def test_converged_fit_does_not_report_a_cap_exit(log: object) -> None:
    # ARRANGE: PLENTY OF PASSES, SO A PASS THAT CLIPS NOTHING ENDS THE LOOP
    detector = _detector(log, orderDeg=0, axisBDeg=1)
    detector.recipeSettings["poly-clipping-iteration-limit"] = 10
    pixels, _ = _line_with_outliers()

    # ACT
    detector.fit_global_polynomial(pixels.copy())

    # ASSERT
    assert not [message for level, message in log.messages if "iteration limit" in message]


def test_cap_exit_returns_no_solution_when_the_survivors_cannot_constrain_the_fit(log: object) -> None:
    # ARRANGE: FOUR POINTS ON A LINE WITH ONE GROSS OUTLIER, ONE PASS, AND THREE COEFFICIENTS
    detector = _detector(log, orderDeg=0, axisBDeg=2, clippingSigma=1.0)
    detector.recipeSettings["poly-clipping-iteration-limit"] = 1
    pixels = pd.DataFrame(
        {
            "order": np.full(4, 10.0),
            "cont_y": [0.0, 1.0, 2.0, 3.0],
            "cont_x": [1.0, 3.0, 5.0, 500.0],
        }
    )

    # ACT
    coefficients, kept, clipped = detector.fit_global_polynomial(pixels.copy())

    # ASSERT: FEWER SURVIVORS THAN COEFFICIENTS MEANS NO SOLUTION, AS ON ANY OTHER PASS
    assert len(clipped) > 0
    assert coefficients is None
    assert len(kept) < 3
