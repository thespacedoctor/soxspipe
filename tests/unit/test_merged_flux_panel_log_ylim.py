"""The log-scale flux panel of the merged-spectrum QC plot gets a positive lower limit, so matplotlib never warns."""

from __future__ import annotations

import warnings

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from soxspipe.commonutils.toolkit import _draw_merged_flux_panels  # noqa: E402

pytestmark = pytest.mark.unit


def _draw(flux: np.ndarray) -> plt.Axes:
    """Draw the two flux panels for a spectrum and return the log-scale panel."""
    merged = pd.DataFrame({"WAVE": np.linspace(500.0, 600.0, flux.size), "FLUX_COUNTS": flux})
    fig = plt.figure()
    gs = fig.add_gridspec(2, 1)
    try:
        _, middlePanel = _draw_merged_flux_panels(fig, gs, merged, "vis", "spectrum.fits", fluxCalibrated=False)
        return middlePanel
    finally:
        plt.close(fig)


def test_a_spectrum_with_zero_and_negative_flux_gets_a_positive_log_lower_limit_without_a_warning() -> None:
    # ARRANGE
    # A SMOOTH RAMP FROM -20 TO 100, SO THE 3-SIGMA LOWER CLIP KEEPS THE ZERO AND NEGATIVE SAMPLES
    flux = np.linspace(-20.0, 100.0, 241)

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        middlePanel = _draw(flux)

    # ASSERT
    lower, upper = middlePanel.get_ylim()
    assert lower == pytest.approx(0.5 * flux[flux > 0].min())
    assert upper == pytest.approx(200.0)


def test_a_spectrum_with_only_positive_flux_keeps_half_its_minimum_as_the_log_lower_limit() -> None:
    # ARRANGE
    flux = np.linspace(40.0, 100.0, 200)

    # ACT
    middlePanel = _draw(flux)

    # ASSERT
    assert middlePanel.get_ylim() == pytest.approx((20.0, 200.0))


def test_a_spectrum_with_no_positive_flux_sets_no_log_limit() -> None:
    # ARRANGE
    flux = np.linspace(-100.0, -10.0, 200)

    # ACT
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        _draw(flux)

    # ASSERT
    assert [warning for warning in raised if "non-positive ylim" in str(warning.message)] == []
