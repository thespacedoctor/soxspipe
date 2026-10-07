"""Contracts for the slit-position drift QC series and plot."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_drift_qc import (
    SLIT_DRIFT_FIT_SAMPLES,
    plot_slit_drift_qc,
    slit_drift_series,
)
from tests.unit._plot_spies import spy_figures

pytestmark = pytest.mark.unit

MEAN_SLIT_ARCSEC = 0.25


def _trace_table() -> pd.DataFrame:
    """Two orders: order 10 (red) drifts linearly with noise, order 11 (blue) is flat."""
    return pd.DataFrame(
        {
            "order": [10, 10, 10, 11, 11],
            "wavelength": [600.0, 601.0, 602.0, 500.0, 501.0],
            "slit_position": [1.1, 1.9, 3.2, 0.5, 0.5],
        }
    )


def _two_order_series(**overrides: object) -> list[dict]:
    arguments = {
        "orderPixelTable": _trace_table(),
        "orders": [10, 11],
        "centreCoeffs": [np.array([1.0, -599.0]), np.array([0.5])],
        "wlMinMax": [(600.0, 602.0), (500.0, 501.0)],
        "fallbackFlags": [False, False],
    }
    arguments.update(overrides)
    return slit_drift_series(**arguments)


def test_residual_is_trace_point_minus_fitted_polynomial() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order10 = next(entry for entry in series if entry["order"] == 10)
    np.testing.assert_allclose(order10["traceSlit"], [1.1, 1.9, 3.2])
    np.testing.assert_allclose(order10["residual"], [0.1, -0.1, 0.2], atol=1e-12)
    assert order10["isFallback"] is False


def test_non_finite_trace_points_are_dropped() -> None:
    # ARRANGE
    table = pd.DataFrame(
        {
            "order": [10, 10, 10, 10],
            "wavelength": [600.0, np.nan, 602.0, 603.0],
            "slit_position": [1.0, 2.0, np.inf, 4.0],
        }
    )

    # ACT
    series = slit_drift_series(table, [10], [np.array([1.0, -599.0])], [(600.0, 603.0)], [False])

    # ASSERT
    np.testing.assert_allclose(series[0]["traceWavelength"], [600.0, 603.0])
    np.testing.assert_allclose(series[0]["traceSlit"], [1.0, 4.0])
    assert len(series[0]["residual"]) == 2


def test_only_points_for_the_order_are_used() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order11 = next(entry for entry in series if entry["order"] == 11)
    np.testing.assert_allclose(order11["traceWavelength"], [500.0, 501.0])


def test_fallback_order_has_no_trace_and_constant_fit_over_its_wavelength_range() -> None:
    # ARRANGE
    table = _trace_table()
    table = table[table["order"] == 10]

    # ACT
    series = slit_drift_series(
        table,
        [10, 11],
        [np.array([1.0, -599.0]), np.array([MEAN_SLIT_ARCSEC])],
        [(600.0, 602.0), (500.0, 501.0)],
        [False, True],
    )

    # ASSERT
    fallback = next(entry for entry in series if entry["order"] == 11)
    assert fallback["isFallback"] is True
    assert len(fallback["traceWavelength"]) == 0
    assert len(fallback["residual"]) == 0
    assert len(fallback["fitWavelength"]) == SLIT_DRIFT_FIT_SAMPLES
    assert fallback["fitWavelength"][0] == pytest.approx(500.0)
    assert fallback["fitWavelength"][-1] == pytest.approx(501.0)
    np.testing.assert_allclose(fallback["fitSlit"], MEAN_SLIT_ARCSEC)


def test_fit_curve_is_the_polynomial_sampled_across_the_order_range() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order10 = next(entry for entry in series if entry["order"] == 10)
    expected = np.polyval([1.0, -599.0], order10["fitWavelength"])
    np.testing.assert_allclose(order10["fitSlit"], expected)


def test_series_is_sorted_blue_to_red_by_minimum_wavelength() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    assert [entry["order"] for entry in series] == [11, 10]


def test_plot_writes_a_pdf_and_leaves_no_figure_open(tmp_path: Path) -> None:
    # ARRANGE
    series = _two_order_series()
    filePath = tmp_path / "SLIT_DRIFT.pdf"
    plt.close("all")

    # ACT
    result = plot_slit_drift_qc(series, MEAN_SLIT_ARCSEC, "title", str(filePath))

    # ASSERT
    assert result == str(filePath)
    assert filePath.read_bytes().startswith(b"%PDF")
    assert plt.get_fignums() == []


def test_plot_has_two_labelled_panels_and_a_mean_line(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    closed: list[tuple] = []
    realClose = plt.close

    def recording_close(*args: object) -> None:
        closed.append(args)
        realClose(*args)

    monkeypatch.setattr(plt, "close", recording_close)

    # ACT
    plot_slit_drift_qc(_two_order_series(), MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plot.pdf"))

    # ASSERT
    axes = figures[0].axes
    assert len(axes) == 2
    mainAxis, residualAxis = axes
    assert residualAxis.get_xlabel() == "wavelength (nm)"
    assert mainAxis.get_ylabel() == "slit position (arcsec)"
    assert residualAxis.get_ylabel() == "residual (arcsec)"
    dashed = [line for line in mainAxis.get_lines() if line.get_linestyle() == "--"]
    assert len(dashed) == 1
    assert list(dashed[0].get_ydata()) == [MEAN_SLIT_ARCSEC, MEAN_SLIT_ARCSEC]
    assert "former single centre" in dashed[0].get_label()
    assert "#555" not in dashed[0].get_label()
    assert f"{MEAN_SLIT_ARCSEC:.3f}" in dashed[0].get_label()
    assert closed == [(figures[0],)]


def test_plot_marks_fallback_order_with_a_span_and_text(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    series = _two_order_series(fallbackFlags=[False, True], orderPixelTable=_trace_table()[lambda t: t["order"] == 10])

    # ACT
    plot_slit_drift_qc(series, MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plot.pdf"))

    # ASSERT
    mainAxis = figures[0].axes[0]
    texts = [text.get_text() for text in mainAxis.texts]
    assert "order 11: no usable trace fit, fell back to mean" in texts
    assert len(mainAxis.patches) == 1


def test_plot_closes_figure_when_saving_fails(tmp_path: Path) -> None:
    # ARRANGE
    plt.close("all")
    missingDirectory = tmp_path / "missing" / "plot.pdf"

    # ACT / ASSERT
    with pytest.raises(FileNotFoundError):
        plot_slit_drift_qc(_two_order_series(), MEAN_SLIT_ARCSEC, "title", str(missingDirectory))
    assert plt.get_fignums() == []
