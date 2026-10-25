"""The trace-failure gate of `detect_continuum.get` judges the measured rows, not the synthetic ones (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum

pytestmark = pytest.mark.unit

MEASURED_ROWS = 40
SYNTHETIC_ROWS = 60
MEASURED_RESIDUAL_PX = 3.0
GATE_LIMIT_PX = 2.0


def _residual_table(withTransferColumn: bool) -> pd.DataFrame:
    """Return fitted rows: measured ones 3 px off the trace, and more synthetic ones exactly on it."""
    count = MEASURED_ROWS + SYNTHETIC_ROWS
    isSynthetic = np.arange(count) >= MEASURED_ROWS
    table = pd.DataFrame(
        {
            "order": np.where(isSynthetic, 2.0, 3.0),
            "cont_x": np.linspace(10.0, 200.0, count),
            "cont_y": np.linspace(10.0, 2000.0, count),
            "gauss_stddev": 1.5,
            "pre-clipped": False,
        }
    )
    return table.assign(faint_transfer=isSynthetic) if withTransferColumn else table


def _residuals(withTransferColumn: bool) -> np.ndarray:
    """Return the fit residual of each row: the measured rows 3 px off, the synthetic rows (if marked) on the trace."""
    isSynthetic = np.arange(MEASURED_ROWS + SYNTHETIC_ROWS) >= MEASURED_ROWS
    return np.where(isSynthetic & withTransferColumn, 0.0, MEASURED_RESIDUAL_PX)


def _run_get(log: object, table: pd.DataFrame, residuals: np.ndarray) -> tuple:
    """Run `get()` on `table` with a stubbed fit that reports `residuals` for the trace fit."""
    detector = object.__new__(detect_continuum)
    detector.log, detector.inst, detector.arm = log, "SOXS", "VIS"
    detector.recipeName = "soxs-order-centres"
    detector.axisA, detector.axisB = "x", "y"
    detector.orderDeg, detector.axisBDeg = 1, 1
    detector.coeff_dict = {}
    detector.recipeSettings = {}
    detector.settings = {"tune-pipeline": False}
    detector.qc = pd.DataFrame()
    detector.products = False
    detector.orderPixelTable = table
    fitCalls = []

    def fit_global_polynomial(pixelList: pd.DataFrame, **kwargs: object) -> tuple:
        fitCalls.append(kwargs["axisACol"])
        if kwargs["axisACol"] != "cont_x":
            raise RuntimeError("the gate let the trace through to the width fit")
        return np.zeros(4), pixelList.assign(cont_x_fit_res=residuals), pixelList.iloc[0:0]

    detector.fit_global_polynomial = fit_global_polynomial
    return detector.get(), fitCalls


def test_synthetic_rows_do_not_hide_a_failed_trace_from_the_mean_residual_gate(log: object) -> None:
    # ARRANGE: THE MEASURED ROWS AVERAGE 3 PX, ABOVE THE 2 PX LIMIT, BUT THE SYNTHETIC ROWS PULL THE MEAN OF ALL ROWS
    # DOWN TO 1.2 PX
    table = _residual_table(withTransferColumn=True)
    residuals = _residuals(withTransferColumn=True)
    assert np.abs(residuals).mean() < GATE_LIMIT_PX

    # ACT
    result, fitCalls = _run_get(log, table, residuals)

    # ASSERT: THE GATE FIRES ON EVERY ATTEMPT AND THE RECIPE GIVES UP, AS FOR A TABLE OF MEASURED ROWS ONLY
    assert result[0] is None
    assert fitCalls == ["cont_x"] * 5


def test_a_table_without_the_transfer_column_is_gated_on_all_its_rows(log: object) -> None:
    # ARRANGE: THE SAME RESIDUALS WITH NO FAINT_TRANSFER COLUMN, SO THE MEAN OVER ALL ROWS IS THE DEVELOP BEHAVIOUR
    table = _residual_table(withTransferColumn=False)

    # ACT
    result, fitCalls = _run_get(log, table, _residuals(withTransferColumn=False))

    # ASSERT
    assert result[0] is None
    assert fitCalls == ["cont_x"] * 5
