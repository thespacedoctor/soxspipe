"""Analytic tests for the public polynomial evaluators."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_allclose, assert_array_equal

from soxspipe.commonutils.polynomials import (
    chebyshev_order_wavelength_polynomials,
    chebyshev_order_xy_polynomials,
    chebyshev_xy_polynomial,
)

pytestmark = pytest.mark.unit


def test_order_wavelength_polynomial_matches_worked_example(log: object) -> None:
    points = pd.DataFrame(
        {
            "order": [1.0, 2.0],
            "wavelength": [10.0, 20.0],
            "slit_position": [0.0, 1.0],
        }
    )
    polynomial = chebyshev_order_wavelength_polynomials(
        log=log,
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=0,
    ).poly

    actual = polynomial(points, 1.0, 2.0, 3.0, 4.0)

    assert_allclose(actual, [64.0, 207.0])


def test_order_wavelength_polynomial_accepts_precomputed_axis_powers(
    log: object,
) -> None:
    points = pd.DataFrame(
        {
            "order_pow_x_0": [1.0, 1.0],
            "order_pow_x_1": [1.0, 2.0],
            "wavelength_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_1": [10.0, 20.0],
            "slit_position_pow_x_0": [1.0, 1.0],
        }
    )
    polynomial = chebyshev_order_wavelength_polynomials(
        log=log,
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=0,
        exponentsIncluded=True,
        axis="x",
    ).poly

    assert_allclose(polynomial(points, 1.0, 2.0, 3.0, 4.0), [64.0, 207.0])


def test_order_wavelength_polynomial_preserves_empty_shape(log: object) -> None:
    points = pd.DataFrame(columns=["order", "wavelength", "slit_position"])
    polynomial = chebyshev_order_wavelength_polynomials(
        log=log,
        orderDeg=1,
        wavelengthDeg=1,
        slitDeg=0,
    ).poly

    result = polynomial(points, 1.0, 2.0, 3.0, 4.0)

    assert_array_equal(result, np.zeros(0))


def test_xy_polynomial_accepts_an_array(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=2).poly

    assert_allclose(polynomial([0.0, 1.0, 2.0], 1.0, 2.0, 3.0), [1.0, 6.0, 17.0])


def test_order_xy_polynomial_matches_worked_example(log: object) -> None:
    points = pd.DataFrame({"order": [1.0, 2.0], "y": [10.0, 20.0]})
    polynomial = chebyshev_order_xy_polynomials(
        log=log,
        orderDeg=1,
        axisBDeg=1,
        axisB="y",
        axisBCol="y",
        orderCol="order",
    ).poly

    assert_allclose(polynomial(points, 1.0, 2.0, 3.0, 4.0), [64.0, 207.0])
