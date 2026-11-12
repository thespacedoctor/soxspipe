"""The sky fit scales the master flat's gradient by the per-order factor k of SOXS NIR (DY-2150).

``_slit_normalisation_ratio`` takes the factor of the order from ``gradientScales``, which ``subtract_sky`` loads from
the detector parameters of the object frame's arm, and passes it to ``flat_slit_ratio``. These tests read the slope of
the ratio it returns along the slit.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.flat_slit_gradient import load_gradient_scales
from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

COLUMNS = 120
START_WAVELENGTH = 1700.0
DISPERSION = 0.37
KEPT_RANGE = (-5.0, 5.0)
# THE MEASURED NIR GRADIENT OF THE FLAT (DY-1990), IN FRACTION PER ARCSEC
FLAT_SLOPE = -0.0035
SETTING_KEY = "flat-slit-gradient"


def _subtractor(log: Any, instrument: str, arm: str, order: int, isOn: bool = True) -> subtract_sky:
    """A subtract_sky with only the attributes ``_slit_normalisation_ratio`` reads."""
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    detectorParams = detector_lookup(log=log, settings={"instrument": instrument}).get(arm)
    subtractor.gradientScales = load_gradient_scales(detectorParams)
    subtractor.recipeSettings = {"sky-subtraction": {SETTING_KEY: isOn}}
    subtractor.flatFitRanges = {order: KEPT_RANGE}
    return subtractor


def _order_map(order: int) -> pd.DataFrame:
    """One order whose master flat has a straight slit profile, with every pixel usable by the sky fit."""
    slit = np.tile(np.round(np.arange(-120, 121) * 0.05, 2), COLUMNS)
    column = np.repeat(np.arange(COLUMNS), 241)
    return pd.DataFrame(
        {
            "order": float(order),
            "slit_position": slit,
            "wavelength": START_WAVELENGTH + DISPERSION * column,
            "flagged_all_clipped": np.abs(slit) >= 4.0,
            "flat": (500.0 + 10.0 * column) * (1.0 + FLAT_SLOPE * slit),
        }
    )


def _ratio_slope(log: Any, instrument: str, arm: str, order: int) -> float:
    """The slope along the slit of the ratio ``_slit_normalisation_ratio`` gives, per arcsec."""
    imageMapOrder = _order_map(order)
    ratio = _subtractor(log, instrument, arm, order)._slit_normalisation_ratio(imageMapOrder)
    isFitted = imageMapOrder["slit_position"].between(*KEPT_RANGE).to_numpy()
    return np.polyfit(imageMapOrder["slit_position"].to_numpy()[isFitted], ratio[isFitted], 1)[0]


def test_ratio_gradient_is_k_times_the_flat_gradient_for_a_soxs_nir_order_in_the_table(log: Any) -> None:
    # ARRANGE
    unscaledSlope = _ratio_slope(log, "xsh", "NIR", 17)

    # ACT
    slope = _ratio_slope(log, "soxs", "NIR", 17)

    # ASSERT: ORDER 17 HAS K = 0.6
    assert slope == pytest.approx(0.6 * unscaledSlope, rel=0.01)


def test_ratio_gradient_is_the_flat_gradient_for_a_soxs_nir_order_with_k_of_one(log: Any) -> None:
    # ARRANGE
    unscaledSlope = _ratio_slope(log, "xsh", "NIR", 12)

    # ACT
    slope = _ratio_slope(log, "soxs", "NIR", 12)

    # ASSERT
    assert slope == pytest.approx(unscaledSlope, rel=1e-9)
    assert abs(unscaledSlope) > 0.003


@pytest.mark.parametrize(("instrument", "arm"), [("soxs", "VIS"), ("xsh", "NIR"), ("xsh", "VIS")])
def test_ratio_gradient_is_the_flat_gradient_for_vis_and_xshooter_in_a_table_order(
    log: Any, instrument: str, arm: str
) -> None:
    # ACT
    slope = _ratio_slope(log, instrument, arm, 17)

    # ASSERT: THE RATIO IS 1 / P(s), WHOSE SLOPE IS MINUS THE SLOPE OF THE FLAT PROFILE
    assert slope == pytest.approx(-FLAT_SLOPE, rel=0.02)


def test_the_log_gives_the_k_used_in_each_order(log: Any) -> None:
    # ACT
    _ratio_slope(log, "soxs", "NIR", 17)
    _ratio_slope(log, "soxs", "NIR", 12)

    # ASSERT
    messages = [message for _, message in log.messages if "the master flat gives a gradient" in message]
    assert messages[0].startswith("order 17:") and messages[0].endswith("scaled by k = 0.60")
    assert messages[1].startswith("order 12:") and messages[1].endswith("scaled by k = 1.00")
