"""Characterization of the slit-illumination profile fit in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.unit._plot_spies import quiet_show, spy_figures

pytestmark = pytest.mark.unit

PIXEL_COUNT = 2000


def _subtractor(log: Any) -> subtract_sky:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.debug = True
    subtractor.binx = 2
    subtractor.biny = 2
    subtractor.recipeSettings = {"sky-subtraction": {"slit_illumination_order": 1}}
    return subtractor


def _slit_profile() -> pd.DataFrame:
    """A linear slit profile with seeded noise; every fiftieth pixel is clipped."""
    slitPosition = np.linspace(-1.0, 1.0, PIXEL_COUNT)
    return pd.DataFrame(
        {
            "order": 11,
            "slit_position": slitPosition,
            "sky_subtracted_flux": 20.0 + 3.0 * slitPosition + np.random.default_rng(8).normal(0.0, 0.5, PIXEL_COUNT),
            "residual_windowed_std": np.full(PIXEL_COUNT, 0.5),
            "flagged_all_clipped": np.arange(PIXEL_COUNT) % 50 == 0,
        }
    )


def test_a_debug_binned_profile_plots_the_polynomial_and_spline_but_applies_no_correction(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """2x2 binning gives 125 points per starter knot; the fits are drawn and the ratio stays 1."""
    figures = spy_figures(monkeypatch)
    shows = quiet_show(monkeypatch)
    pixels = _slit_profile()

    corrected = _subtractor(log).cross_dispersion_flux_normaliser(pixels.copy())

    assert (corrected["slit_normalisation_ratio"] == 1).all()
    pd.testing.assert_series_equal(corrected["sky_subtracted_flux"], pixels["sky_subtracted_flux"])
    assert len(shows) == 1
    [figure] = figures
    [axis] = figure.axes
    assert axis.get_title() == "Slit Illumination Profile"
    dataLine, polynomialLine, splineLine = axis.lines
    # 40 OF THE 2000 PIXELS ARE CLIPPED BEFORE PLOTTING
    assert dataLine.get_xdata().shape == (1960,)
    assert polynomialLine.get_xdata()[[0, -1]].tolist() == [-1.0, 1.0]
    assert polynomialLine.get_xdata().shape == (100,)
    assert polynomialLine.get_ydata()[[0, 50, -1]].tolist() == [
        pytest.approx(16.99133001094896, rel=1e-12, abs=0),
        pytest.approx(20.017965227092574, rel=1e-12, abs=0),
        pytest.approx(22.984067738913314, rel=1e-12, abs=0),
    ]
    assert splineLine.get_label() == "sky model"
    assert splineLine.get_xdata().shape == (10000,)
    assert splineLine.get_ydata()[[0, 5000, -1]].tolist() == [
        pytest.approx(16.973806180278967, rel=1e-12, abs=0),
        pytest.approx(19.939925492462898, rel=1e-12, abs=0),
        pytest.approx(22.9564173596589, rel=1e-12, abs=0),
    ]
    assert float(splineLine.get_ydata().sum()) == pytest.approx(199877.58182724216, rel=1e-12, abs=0)
