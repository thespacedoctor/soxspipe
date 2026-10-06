"""Characterization of the sky-model fitting quick-look plot in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit._plot_spies import quiet_show, spy_figures

pytestmark = pytest.mark.unit

PIXEL_COUNT = 200


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.settings = {}
    subtractor.arm = "VIS"
    return subtractor


def _skylines(*, empty: bool = False) -> pd.DataFrame:
    if empty:
        return pd.DataFrame({"WAVELENGTH": [], "FLUX": [], "ISOLATED": []})
    return pd.DataFrame(
        {
            "WAVELENGTH": [502.0, 505.0, 508.0],
            "FLUX": [10.0, 20.0, 30.0],
            "ISOLATED": [True, False, True],
        }
    )


def _order_pixels() -> pd.DataFrame:
    index = np.arange(PIXEL_COUNT)
    flux = 100.0 + 20.0 * np.sin(index / 7.0)
    flux[::37] = -80.0
    skyLine = np.where(index % 11 == 0, "line", np.where(index % 13 == 0, "node", "")).astype(object)
    skyLine[skyLine == ""] = False
    return pd.DataFrame(
        {
            "wavelength": np.linspace(500.0, 510.0, PIXEL_COUNT),
            "flux": flux,
            "flagged_object_clipped": index % 41 == 0,
            "flagged_all_clipped": index % 29 == 0,
            "flagged_bad_pixel_clipped": index % 53 == 0,
            "flagged_edge_clipped": index < 3,
            "flux_percentile_smoothed": 100.0 + np.cos(index / 5.0),
            "flux_minus_smoothed_residual": 4.0 * np.sin(index / 3.0),
            "flux_minus_smoothed_residual_upper_limit": np.full(PIXEL_COUNT, 12.0),
            "flux_upper_limit": np.full(PIXEL_COUNT, 112.0),
            "flagged_noisy_region": index % 17 == 0,
            "flagged_sky_line": skyLine,
            "flagged_bspline_clipped": index % 19 == 0,
            "sky_residuals": 3.0 * np.abs(np.cos(index / 4.0)),
            "sky_residual_floor": np.full(PIXEL_COUNT, 1.5),
            "sky_residual_rolling_average": 2.0 * np.abs(np.sin(index / 9.0)),
        }
    )


def _spline() -> tuple:
    index = np.arange(PIXEL_COUNT)
    wavelength = np.linspace(500.0, 510.0, PIXEL_COUNT)
    return splrep(wavelength, 100.0 + 20.0 * np.sin(index / 7.0), k=3, s=PIXEL_COUNT * 400)


@pytest.fixture
def plot_recorders(monkeypatch: pytest.MonkeyPatch) -> dict[str, list]:
    # THE PLOT PAUSES FOR A SECOND AND SHOWS WINDOWS; RECORD BOTH INSTEAD
    pauses: list[tuple] = []
    monkeypatch.setattr(plt, "pause", lambda *args, **kwargs: pauses.append(args))
    return {
        "figures": spy_figures(monkeypatch),
        "shows": quiet_show(monkeypatch),
        "pauses": pauses,
    }


@pytest.fixture
def skylines(monkeypatch: pytest.MonkeyPatch) -> list[tuple]:
    # THE STATIC SKYLINE TABLE LIVES IN UNDISTRIBUTED CALIBRATION DATA
    calls: list[tuple] = []

    def fake_skylines(log: Any, settings: Any, arm: str) -> pd.DataFrame:
        calls.append((settings, arm))
        return _skylines()

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", fake_skylines)
    return calls


def test_clipping_view_without_a_spline_pauses_and_returns_before_the_model(
    log: Any, plot_recorders: dict[str, list], skylines: list[tuple]
) -> None:
    """Without a spline the plot shows the clipping stage only, pauses, and closes."""
    subtractor = _subtractor(log)

    result = subtractor.plot_order_skymodel_fitting_quicklook(_order_pixels(), None, title="clipping")

    assert result is None
    assert skylines == [({}, "VIS")]
    assert plot_recorders["pauses"] == [(1,)]
    assert plot_recorders["shows"] == []
    [figure] = plot_recorders["figures"]
    fluxAxis, residualAxis = figure.axes
    assert fluxAxis.get_title() == "clipping"
    assert fluxAxis.get_ylim() == (
        pytest.approx(88.00143995810838, rel=1e-12, abs=0),
        pytest.approx(254.50353330762323, rel=1e-12, abs=0),
    )
    assert fluxAxis.get_xlim() == (
        pytest.approx(500.0502512562814, rel=1e-12, abs=0),
        pytest.approx(510.0, rel=1e-12, abs=0),
    )
    assert residualAxis.get_ylim() == (
        pytest.approx(-13.959204742696885, rel=1e-12, abs=0),
        pytest.approx(14.151629608944601, rel=1e-12, abs=0),
    )
    assert [collection.get_label() for collection in fluxAxis.collections] == [
        "unclipped",
        "flux percentile smoothed (windowed)",
        "bad pixels clipped",
        "object",
    ]
    assert [collection.get_offsets().shape[0] for collection in fluxAxis.collections] == [189, 189, 4, 5]
    assert [collection.get_label() for collection in residualAxis.collections] == ["flux residuals (windowed)"]
    # THREE SKYLINE axvline CALLS PLUS TWO EMPTY LEGEND PROXIES PER AXIS
    assert len(fluxAxis.lines) == 5
    assert len(residualAxis.lines) == 5
    assert fluxAxis.get_legend() is None


def test_fitting_view_with_a_spline_draws_the_model_and_shows_the_figure(
    log: Any, plot_recorders: dict[str, list], skylines: list[tuple]
) -> None:
    """With a spline the plot adds sky-line, noise, and model layers and calls show."""
    subtractor = _subtractor(log)

    result = subtractor.plot_order_skymodel_fitting_quicklook(
        _order_pixels(), _spline(), knots=np.array([503.0, 507.0])
    )

    assert result is None
    assert plot_recorders["pauses"] == []
    assert len(plot_recorders["shows"]) == 1
    [figure] = plot_recorders["figures"]
    fluxAxis, residualAxis = figure.axes
    assert fluxAxis.get_title() == ""
    assert residualAxis.get_ylim() == (
        pytest.approx(0.0, rel=1e-12, abs=0),
        pytest.approx(8.345956191113675, rel=1e-12, abs=0),
    )
    assert [collection.get_label() for collection in fluxAxis.collections] == [
        "unclipped",
        "predicted skyline",
        "predicted skyline wings/peak",
        "bspline clipped",
        "noisy region",
    ]
    assert [collection.get_offsets().shape[0] for collection in fluxAxis.collections] == [189, 18, 14, 11, 12]
    assert [collection.get_label() for collection in residualAxis.collections] == [
        "flux residuals (windowed)",
        "bspline knots",
    ]
    modelLine = fluxAxis.lines[-1]
    assert modelLine.get_label() == "sky model bspline"
    # THE MODEL IS SAMPLED ON A ONE-MILLION POINT GRID
    assert modelLine.get_xdata().shape == (1000000,)
    assert [text.get_text() for text in fluxAxis.get_legend().get_texts()] == [
        "calibration skyline",
        "skyline",
        "unclipped",
        "predicted skyline",
        "predicted skyline wings/peak",
        "bspline clipped",
        "noisy region",
        "sky model bspline",
    ]
    assert [text.get_text() for text in residualAxis.get_legend().get_texts()] == [
        "calibration skyline",
        "skyline",
        "flux residuals (windowed)",
        "sky residual floor",
        "sky residual rolling average",
        "bspline knots",
    ]


def test_fitting_view_without_knots_omits_the_knot_markers(
    log: Any, plot_recorders: dict[str, list], skylines: list[tuple]
) -> None:
    """The boolean default for `knots` skips the knot scatter on the residual panel."""
    subtractor = _subtractor(log)

    subtractor.plot_order_skymodel_fitting_quicklook(_order_pixels(), _spline())

    [figure] = plot_recorders["figures"]
    _, residualAxis = figure.axes
    assert [collection.get_label() for collection in residualAxis.collections] == ["flux residuals (windowed)"]


def test_an_empty_skyline_table_adds_no_skyline_legend_entries(
    log: Any, plot_recorders: dict[str, list], monkeypatch: pytest.MonkeyPatch
) -> None:
    """Without listed skylines there are no vertical lines and no skyline proxies."""
    monkeypatch.setattr(toolkit, "get_skylines_dataframe", lambda log, settings, arm: _skylines(empty=True))
    subtractor = _subtractor(log)

    subtractor.plot_order_skymodel_fitting_quicklook(_order_pixels(), _spline())

    [figure] = plot_recorders["figures"]
    fluxAxis, residualAxis = figure.axes
    assert [line.get_label() for line in fluxAxis.lines] == ["sky model bspline"]
    assert [text.get_text() for text in residualAxis.get_legend().get_texts()] == [
        "flux residuals (windowed)",
        "sky residual floor",
        "sky residual rolling average",
    ]


def test_limit_statistics_receive_only_unclipped_values_above_floor(
    log: Any,
    plot_recorders: dict[str, list],
    skylines: list[tuple],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Limit statistics exclude both clipped pixels and values at or below -50."""
    import astropy.stats

    realStats = astropy.stats.sigma_clipped_stats
    sampleSizes: list[int] = []

    def record_stats(data: Any, **kwargs: Any) -> tuple:
        sampleSizes.append(len(data))
        return realStats(data, **kwargs)

    monkeypatch.setattr(astropy.stats, "sigma_clipped_stats", record_stats)
    subtractor = _subtractor(log)
    pixels = _order_pixels()
    # Include an explicitly unclipped value below the floor in each panel.
    assert not pixels.loc[1, ["flagged_all_clipped", "flagged_object_clipped"]].any()
    pixels.loc[1, ["flux", "flux_minus_smoothed_residual", "sky_residuals"]] = -50.0
    # THE PLOT'S OWN CLIPPED MASK IS ALL-CLIPPED OR OBJECT-CLIPPED
    assert int((~(pixels["flagged_all_clipped"] | pixels["flagged_object_clipped"])).sum()) == 189

    subtractor.plot_order_skymodel_fitting_quicklook(pixels, None)
    subtractor.plot_order_skymodel_fitting_quicklook(pixels, _spline())

    # Five other unclipped flux rows are below the floor in the synthetic data.
    assert sampleSizes == [183, 188, 183, 188]
