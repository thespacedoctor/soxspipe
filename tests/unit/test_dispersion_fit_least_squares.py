"""Exact least-squares contracts for the dispersion and spatial solution x/y fits (DY-1525)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials

pytestmark = pytest.mark.unit

VIS_ORDER_CENTRES = {1.0: 900.0, 2.0: 760.0, 3.0: 660.0, 4.0: 580.0}


def _mapper(log: object, recipeName: str = "soxs-disp-solution") -> create_dispersion_map:
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.recipeName = recipeName
    mapper.settings = {"workspace-root-dir": "/nonexistent"}
    return mapper


def _add_exponents(table: pd.DataFrame, orderDeg: int, wavelengthDeg: int, slitDeg: int) -> pd.DataFrame:
    """Add the precomputed power columns ``fit_polynomials`` adds before calling the x/y fit."""
    powers = {}
    for axis in "xy":
        powers.update({f"order_pow_{axis}_{i}": table["order"] ** i for i in range(orderDeg + 1)})
        powers.update({f"wavelength_pow_{axis}_{j}": table["wavelength"] ** j for j in range(wavelengthDeg + 1)})
        powers.update({f"slit_position_pow_{axis}_{k}": table["slit_position"] ** k for k in range(slitDeg + 1)})
    return table.assign(**powers)


def _vis_like_lines(seed: int, slitPositions: tuple[float, ...] = (0.0,)) -> pd.DataFrame:
    """About 130 arc lines on four VIS-like orders, with smooth x/y traces and 0.2 px centroid noise."""
    rng = np.random.default_rng(seed)
    lineCount = 130
    orders = rng.choice(list(VIS_ORDER_CENTRES), lineCount)
    centres = np.array([VIS_ORDER_CENTRES[order] for order in orders])
    wavelengths = centres * (1 + rng.uniform(-0.08, 0.08, lineCount))
    slit = np.tile(np.asarray(slitPositions, dtype=float), lineCount)
    orders = np.repeat(orders, len(slitPositions))
    centres = np.repeat(centres, len(slitPositions))
    wavelengths = np.repeat(wavelengths, len(slitPositions))
    noiseCount = orders.size
    observedX = 300 * orders + 0.002 * (wavelengths - 700) ** 2 + 4 * slit + rng.normal(0, 0.2, noiseCount)
    observedY = 48000 * (wavelengths / centres - 1) + 2000 + 3 * orders**2 + 0.5 * slit + rng.normal(0, 0.2, noiseCount)
    return pd.DataFrame(
        {
            "order": orders,
            "wavelength": wavelengths,
            "slit_position": slit,
            "observed_x": observedX,
            "observed_y": observedY,
        }
    )


def _optimum_ssr(table: pd.DataFrame, column: str, orderDeg: int, wavelengthDeg: int, slitDeg: int) -> float:
    """Sum of squared residuals at the least-squares minimum, from an independent monomial design matrix."""
    orders = table["order"].to_numpy(dtype=float)
    wavelengths = table["wavelength"].to_numpy(dtype=float)
    slit = table["slit_position"].to_numpy(dtype=float)
    design = np.column_stack(
        [
            orders**i * wavelengths**j * slit**k
            for i in range(orderDeg + 1)
            for j in range(wavelengthDeg + 1)
            for k in range(slitDeg + 1)
        ]
    )
    observed = table[column].to_numpy(dtype=float)
    scale = np.linalg.norm(design, axis=0)
    solution = np.linalg.lstsq(design / scale, observed, rcond=None)[0] / scale
    return float(np.sum((observed - design @ solution) ** 2))


def _fit(log: object, table: pd.DataFrame, orderDeg: int, wavelengthDeg: int, slitDeg: int, **mapperKwargs: str):
    polys = {
        axis: chebyshev_order_wavelength_polynomials(
            log=log,
            orderDeg=orderDeg,
            wavelengthDeg=wavelengthDeg,
            slitDeg=slitDeg,
            exponentsIncluded=True,
            axis=axis,
        ).poly
        for axis in "xy"
    }
    fitTable = _add_exponents(table, orderDeg, wavelengthDeg, slitDeg)
    xcoeff, ycoeff, failure = _mapper(log, **mapperKwargs)._fit_xy_coefficients(
        fitTable, polys["x"], polys["y"], orderDeg, wavelengthDeg, slitDeg
    )
    return polys, fitTable, xcoeff, ycoeff, failure


@pytest.mark.parametrize("seed", [1, 2])
def test_single_pinhole_fit_reaches_the_least_squares_minimum_on_both_axes(log: object, seed: int) -> None:
    # ARRANGE: THE SOXS VIS DISPERSION-SOLUTION DEGREES (ORDER 3, WAVELENGTH 5), WHERE CURVE_FIT STOPPED 8% SHORT
    table = _vis_like_lines(seed)

    # ACT
    polys, fitTable, xcoeff, ycoeff, failure = _fit(log, table, orderDeg=3, wavelengthDeg=5, slitDeg=0)

    # ASSERT
    assert failure is None
    for axis, coeff in (("x", xcoeff), ("y", ycoeff)):
        ssr = float(np.sum((table[f"observed_{axis}"].to_numpy() - polys[axis](fitTable, *coeff)) ** 2))
        assert ssr == pytest.approx(_optimum_ssr(table, f"observed_{axis}", 3, 5, 0), rel=1e-6)


def test_multi_pinhole_fit_reaches_the_least_squares_minimum_with_a_slit_term(log: object) -> None:
    # ARRANGE: THE SOXS VIS SPATIAL-SOLUTION DEGREES (ORDER 3, WAVELENGTH 5, SLIT 1) OVER NINE PINHOLES
    table = _vis_like_lines(3, slitPositions=tuple(np.linspace(-5.0, 5.0, 9)))

    # ACT
    polys, fitTable, xcoeff, ycoeff, failure = _fit(
        log, table, orderDeg=3, wavelengthDeg=5, slitDeg=1, recipeName="soxs-spatial-solution"
    )

    # ASSERT
    assert failure is None
    for axis, coeff in (("x", xcoeff), ("y", ycoeff)):
        ssr = float(np.sum((table[f"observed_{axis}"].to_numpy() - polys[axis](fitTable, *coeff)) ** 2))
        assert ssr == pytest.approx(_optimum_ssr(table, f"observed_{axis}", 3, 5, 1), rel=1e-6)


def test_a_tiny_centroid_perturbation_barely_moves_the_fitted_positions(log: object) -> None:
    # ARRANGE: NUDGE EVERY OBSERVED POSITION BY AT MOST 1E-6 PX
    table = _vis_like_lines(1)
    nudged = table.assign(
        observed_x=table["observed_x"] + np.random.default_rng(7).uniform(-1e-6, 1e-6, len(table)),
        observed_y=table["observed_y"] + np.random.default_rng(8).uniform(-1e-6, 1e-6, len(table)),
    )

    # ACT
    polys, fitTable, xcoeff, ycoeff, _ = _fit(log, table, orderDeg=3, wavelengthDeg=5, slitDeg=0)
    _, nudgedTable, nudgedX, nudgedY, _ = _fit(log, nudged, orderDeg=3, wavelengthDeg=5, slitDeg=0)

    # ASSERT: AN EXACT SOLVE IS LINEAR IN THE DATA, SO THE MODEL MOVES NO MORE THAN THE DATA DID
    np.testing.assert_allclose(polys["x"](nudgedTable, *nudgedX), polys["x"](fitTable, *xcoeff), rtol=0, atol=1e-5)
    np.testing.assert_allclose(polys["y"](nudgedTable, *nudgedY), polys["y"](fitTable, *ycoeff), rtol=0, atol=1e-5)


@pytest.mark.parametrize(
    ("column", "expected"), [("observed_x", ("xerror", None, None, None)), ("observed_y", (None, "yerror", None, None))]
)
def test_a_non_finite_observed_position_returns_the_axis_failure_tuple(
    log: object, column: str, expected: tuple[str | None, ...]
) -> None:
    # ARRANGE
    table = _vis_like_lines(1)
    table.loc[5, column] = np.nan

    # ACT
    _, _, xcoeff, ycoeff, failure = _fit(log, table, orderDeg=3, wavelengthDeg=5, slitDeg=0)

    # ASSERT
    assert (xcoeff, ycoeff) == (None, None)
    assert failure == expected
