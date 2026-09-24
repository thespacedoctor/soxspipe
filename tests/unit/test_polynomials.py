"""Analytic tests for the public polynomial evaluators."""

from __future__ import annotations

from collections.abc import Callable

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


# CHARACTERIZATION FIXTURES: NON-TRIVIAL DEGREES AND NON-INTEGER INPUTS, SO
# THE PINNED VALUES EXERCISE EVERY TERM AND THE COEFFICIENT CONSUMPTION ORDER
ORDERS = [11.0, 14.0, 17.0, 20.0]
WAVELENGTHS = [0.52, 0.61, 0.73, 0.88]
SLIT_POSITIONS = [-2.5, -0.75, 1.25, 3.0]
Y_VALUES = [100.5, 700.25, 1300.0, 1900.75]

ORDER_WAVELENGTH_COEFFS = [
    0.37, -0.172, 0.149333, -0.0535, 0.126, 0.003333, 0.130857, 0.04475,
    0.145111, 0.08, 0.163636, 0.112167, 0.184462, 0.142571, 0.206667, 0.171875,
    0.229765, 0.200444, 0.253474, 0.2285, 0.277619, 0.256182, 0.302087, 0.283583,
]  # fmt: skip
ORDER_WAVELENGTH_EXPECTED = [
    -72.767460242592,
    38.800235616650994,
    432.802390821393,
    1337.228486976384,
]

XY_COEFFS = [3.25, -0.0125, 4.5e-6, -7.5e-10]
XY_INPUT = [0.0, 512.5, 1024.0, 2047.75]
XY_EXPECTED = [3.25, -2.0752553710937502, -5.636714368000001, -9.917206654738283]

ORDER_XY_COEFFS = [
    0.85, 0.566667, -0.425, -0.34, 0.283333, 0.242857,
    -0.2125, -0.188889, 0.17, 0.154545, -0.141667, -0.130769,
]  # fmt: skip
ORDER_XY_EXPECTED = [
    -18714618.91937275,
    -9840781723.206762,
    -90907299471.31004,
    -387701668378.98254,
]


def _order_wavelength_raw_table() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "order": ORDERS,
            "wavelength": WAVELENGTHS,
            "slit_position": SLIT_POSITIONS,
        }
    )


def _order_wavelength_power_table(axisPrefix: str) -> pd.DataFrame:
    columns = {}
    for name, values, degree in (
        ("order", ORDERS, 2),
        ("wavelength", WAVELENGTHS, 3),
        ("slit_position", SLIT_POSITIONS, 1),
    ):
        for power in range(degree + 1):
            columns[f"{name}_pow_{axisPrefix}{power}"] = np.array(values) ** power
    return pd.DataFrame(columns)


def _order_wavelength_poly(log: object, **kwargs: object) -> Callable[..., np.ndarray]:
    return chebyshev_order_wavelength_polynomials(log=log, orderDeg=2, wavelengthDeg=3, slitDeg=1, **kwargs).poly


def _order_xy_power_table() -> pd.DataFrame:
    columns = {}
    for power in range(3):
        columns[f"order_pow_{power}"] = np.array(ORDERS) ** power
    for power in range(4):
        columns[f"y_pow_{power}"] = np.array(Y_VALUES) ** power
    return pd.DataFrame(columns)


def test_order_wavelength_polynomial_pins_a_full_degree_fit(log: object) -> None:
    polynomial = _order_wavelength_poly(log)

    result = polynomial(_order_wavelength_raw_table(), *ORDER_WAVELENGTH_COEFFS)

    assert result == pytest.approx(ORDER_WAVELENGTH_EXPECTED, rel=1e-12, abs=0)


@pytest.mark.parametrize("axis", ["x", "Y"])
def test_order_wavelength_polynomial_pins_precomputed_axis_powers(log: object, axis: str) -> None:
    polynomial = _order_wavelength_poly(log, exponentsIncluded=True, axis=axis)
    table = _order_wavelength_power_table(f"{axis.lower()}_")

    result = polynomial(table, *ORDER_WAVELENGTH_COEFFS)

    assert result == pytest.approx(ORDER_WAVELENGTH_EXPECTED, rel=1e-12, abs=0)


def test_order_wavelength_polynomial_without_axis_reads_unprefixed_powers(
    log: object,
) -> None:
    polynomial = _order_wavelength_poly(log, exponentsIncluded=True)

    result = polynomial(_order_wavelength_power_table(""), *ORDER_WAVELENGTH_COEFFS)

    assert result == pytest.approx(ORDER_WAVELENGTH_EXPECTED, rel=1e-12, abs=0)


def test_order_wavelength_polynomial_treats_none_exponents_flag_as_precomputed(
    log: object,
) -> None:
    # `== False` IS FALSE FOR NONE, SO NONE TAKES THE PRECOMPUTED-POWERS BRANCH
    polynomial = _order_wavelength_poly(log, exponentsIncluded=None)

    result = polynomial(_order_wavelength_power_table(""), *ORDER_WAVELENGTH_COEFFS)

    assert result == pytest.approx(ORDER_WAVELENGTH_EXPECTED, rel=1e-12, abs=0)
    with pytest.raises(KeyError):
        polynomial(_order_wavelength_raw_table(), *ORDER_WAVELENGTH_COEFFS)


@pytest.mark.parametrize("flag", [0, np.False_])
def test_order_wavelength_polynomial_treats_false_equal_flags_as_raw(log: object, flag: object) -> None:
    polynomial = _order_wavelength_poly(log, exponentsIncluded=flag)

    result = polynomial(_order_wavelength_raw_table(), *ORDER_WAVELENGTH_COEFFS)

    assert result == pytest.approx(ORDER_WAVELENGTH_EXPECTED, rel=1e-12, abs=0)


def test_order_wavelength_polynomial_rejects_a_wrong_coefficient_count(
    log: object,
) -> None:
    polynomial = _order_wavelength_poly(log)

    with pytest.raises(ValueError):
        polynomial(_order_wavelength_raw_table(), *ORDER_WAVELENGTH_COEFFS[:-1])


def test_order_wavelength_polynomial_logs_entry_and_exit_at_debug(
    log: object,
) -> None:
    polynomial = _order_wavelength_poly(log)

    polynomial(_order_wavelength_raw_table(), *ORDER_WAVELENGTH_COEFFS)

    assert log.messages == [
        ("debug", "starting the ``poly`` method"),
        ("debug", "completed the ``poly`` method"),
    ]


def test_xy_polynomial_pins_a_cubic_on_an_array(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=3).poly

    result = polynomial(XY_INPUT, *XY_COEFFS)

    assert result == pytest.approx(XY_EXPECTED, rel=1e-12, abs=0)


def test_xy_polynomial_reads_the_named_column_of_a_dataframe(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=3, yCol="ypix").poly
    table = pd.DataFrame({"ypix": XY_INPUT, "other": [9.0, 9.0, 9.0, 9.0]})

    result = polynomial(table, *XY_COEFFS)

    assert result == pytest.approx(XY_EXPECTED, rel=1e-12, abs=0)


def test_xy_polynomial_reads_precomputed_y_powers(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=3, yCol="ypix", exponentsIncluded=True).poly
    table = pd.DataFrame({f"y_pow_{power}": np.array(XY_INPUT) ** power for power in range(4)} | {"ypix": XY_INPUT})

    result = polynomial(table, *XY_COEFFS)

    assert result == pytest.approx(XY_EXPECTED, rel=1e-12, abs=0)


def test_xy_polynomial_ignores_surplus_coefficients(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=3).poly

    result = polynomial(XY_INPUT, *XY_COEFFS, 123.0)

    assert result == pytest.approx(XY_EXPECTED, rel=1e-12, abs=0)


def test_xy_polynomial_logs_entry_and_exit_at_info(log: object) -> None:
    polynomial = chebyshev_xy_polynomial(log=log, y_deg=3).poly

    polynomial(XY_INPUT, *XY_COEFFS)

    assert log.messages == [
        ("info", "starting the ``poly`` method"),
        ("info", "completed the ``poly`` method"),
    ]


def test_order_xy_polynomial_pins_a_full_degree_fit(log: object) -> None:
    polynomial = chebyshev_order_xy_polynomials(
        log=log, orderDeg=2, axisBDeg=3, axisB="y", axisBCol="y", orderCol="order"
    ).poly
    table = pd.DataFrame({"order": ORDERS, "y": Y_VALUES})

    result = polynomial(table, *ORDER_XY_COEFFS)

    assert result == pytest.approx(ORDER_XY_EXPECTED, rel=1e-12, abs=0)


def test_order_xy_polynomial_reads_precomputed_powers_for_axis_b(
    log: object,
) -> None:
    polynomial = chebyshev_order_xy_polynomials(log=log, orderDeg=2, axisBDeg=3, axisB="y", exponentsIncluded=True).poly

    result = polynomial(_order_xy_power_table(), *ORDER_XY_COEFFS)

    assert result == pytest.approx(ORDER_XY_EXPECTED, rel=1e-12, abs=0)


def test_order_xy_polynomial_treats_none_exponents_flag_as_precomputed(
    log: object,
) -> None:
    # `== False` IS FALSE FOR NONE, SO NONE TAKES THE PRECOMPUTED-POWERS BRANCH
    polynomial = chebyshev_order_xy_polynomials(log=log, orderDeg=2, axisBDeg=3, axisB="y", exponentsIncluded=None).poly

    result = polynomial(_order_xy_power_table(), *ORDER_XY_COEFFS)

    assert result == pytest.approx(ORDER_XY_EXPECTED, rel=1e-12, abs=0)


@pytest.mark.parametrize("flag", [0, np.False_])
def test_order_xy_polynomial_treats_false_equal_flags_as_raw(log: object, flag: object) -> None:
    polynomial = chebyshev_order_xy_polynomials(
        log=log,
        orderDeg=2,
        axisBDeg=3,
        axisBCol="y",
        orderCol="order",
        exponentsIncluded=flag,
    ).poly
    table = pd.DataFrame({"order": ORDERS, "y": Y_VALUES})

    result = polynomial(table, *ORDER_XY_COEFFS)

    assert result == pytest.approx(ORDER_XY_EXPECTED, rel=1e-12, abs=0)


def test_order_xy_polynomial_truncates_float_degrees(log: object) -> None:
    polynomial = chebyshev_order_xy_polynomials(
        log=log,
        orderDeg=2.9,
        axisBDeg=3.2,
        axisBCol="y",
        orderCol="order",
    ).poly
    table = pd.DataFrame({"order": ORDERS, "y": Y_VALUES})

    result = polynomial(table, *ORDER_XY_COEFFS)

    assert result == pytest.approx(ORDER_XY_EXPECTED, rel=1e-12, abs=0)
