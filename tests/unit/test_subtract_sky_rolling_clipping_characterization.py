"""Characterization of rolling-window object clipping in subtract_sky (DY-40)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.subtract_sky import subtract_sky

pytestmark = pytest.mark.unit

PIXEL_COUNT = 2000
# PINNING EVERY ROW AS A LITERAL IS IMPRACTICAL; THESE ROWS SAMPLE THE FIRST
# OBJECT FLANK, A CLIPPED ROW, THE THIRD OBJECT AND THE LAST ROW
SAMPLED_ROWS = [5, 400, 1000, 1999]


def _subtractor(log: Any, *, arm: str) -> tuple[subtract_sky, list[str]]:
    subtractor = subtract_sky.__new__(subtract_sky)
    subtractor.log = log
    subtractor.arm = arm
    subtractor.debug = True
    subtractor.stopSubtraction = False
    iterations: list[str] = []
    # DEBUG MODE PLOTS EVERY ITERATION; RECORD THE ITERATION LINE OF EACH TITLE INSTEAD
    subtractor.plot_order_skymodel_fitting_quicklook = lambda frame, spline, title=None: iterations.append(
        title.split("\n")[1]
    )
    return subtractor, iterations


def _object_order(level: float) -> pd.DataFrame:
    """Seeded shot noise at `level` plus four gaussian objects 4 to 16 pixels wide."""
    rng = np.random.default_rng(5)
    index = np.arange(PIXEL_COUNT)
    flux = level + rng.normal(0.0, 1.0, PIXEL_COUNT) * np.sqrt(level)
    for position, width in enumerate((4, 8, 12, 16), start=1):
        flux += 30.0 * np.sqrt(level) * np.exp(-0.5 * ((index - position * PIXEL_COUNT / 5) / width) ** 2)
    return pd.DataFrame(
        {
            "order": 10,
            "flux": flux,
            "flagged_all_clipped": False,
            "flagged_object_clipped": False,
        }
    )


def _clip(subtractor: subtract_sky, pixels: pd.DataFrame, *, sigma: float = 3.0, iterations: int = 10) -> Any:
    return subtractor.rolling_window_clipping(pixels, windowSize=11, sigma_clip_limit=sigma, max_iterations=iterations)


def _messages(log: Any, *levels: str) -> list[tuple[str, str]]:
    return [(level, message) for level, message in log.messages if level in levels]


def _approx_or_nan(values: list[float]) -> list[Any]:
    return [pytest.approx(value, rel=1e-12, abs=0, nan_ok=True) for value in values]


def test_a_faint_vis_order_retries_clipping_after_five_low_yield_iterations(log: Any) -> None:
    """Under 5% clipped at iteration 5 restarts clipping four times before it converges.

    Each restart clears the object flags but not `flagged_all_clipped`, so 430
    pixels end up excluded while only 212 are flagged as object.
    """
    # SUSPICIOUS (WRONG SCIENCE): RESET PIXELS STAY EXCLUDED, FILED AS DY-595
    subtractor, iterations = _subtractor(log, arm="VIS")

    clipped = _clip(subtractor, _object_order(100.0))

    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
        "iteration 1 - 0.7% clipped",
        "iteration 2 - 1.1% clipped",
        "iteration 3 - 1.6% clipped",
        "iteration 4 - 1.8% clipped",
        "iteration 5 - 2.0% clipped",
        "iteration 1 - 1.3% clipped",
        "iteration 2 - 1.7% clipped",
        "iteration 3 - 1.8% clipped",
        "iteration 4 - 1.9% clipped",
        "iteration 5 - 1.9% clipped",
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.6% clipped",
        "iteration 3 - 3.2% clipped",
        "iteration 4 - 3.5% clipped",
        "iteration 5 - 3.8% clipped",
        "iteration 1 - 4.7% clipped",
        "iteration 2 - 7.8% clipped",
        "iteration 3 - 9.2% clipped",
        "iteration 4 - 10.1% clipped",
        "iteration 5 - 10.4% clipped",
        "iteration 6 - 10.6% clipped",
    ]
    assert subtractor.stopSubtraction is False
    assert int(clipped["flagged_object_clipped"].sum()) == 212
    assert int(clipped["flagged_all_clipped"].sum()) == 430
    assert _messages(log, "print", "warning") == [("print", "\tORDER 10: 212 pixels clipped in total = 10.6%)")]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(1603.8692183469568, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 439
    assert globalSigma.iloc[SAMPLED_ROWS].tolist() == _approx_or_nan(
        [1.0136609843031115, np.nan, 0.6894122443363749, np.nan]
    )
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(150988.5760495704, rel=1e-12, abs=0)
    assert int(clipped["flux_percentile_smoothed"].isna().sum()) == 321


def test_a_bright_vis_order_stops_sky_subtraction_at_iteration_five(log: Any) -> None:
    """A smoothed sky above 2500 with under 5% clipped returns None and raises the stop flag."""
    subtractor, iterations = _subtractor(log, arm="VIS")

    result = _clip(subtractor, _object_order(3000.0))

    assert result is None
    assert subtractor.stopSubtraction is True
    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
    ]
    assert _messages(log, "print", "warning") == [
        ("warning", "OBJECT IS LIKELY VERY BRIGHT - STOPPING SKY-SUBTRACTION TO AVOID CLIPPING TOO MANY PIXELS")
    ]


def test_the_same_bright_order_in_uvb_iterates_to_convergence(log: Any) -> None:
    """The retry and bright-object stop are VIS-only; UVB keeps clipping until nothing changes."""
    subtractor, iterations = _subtractor(log, arm="UVB")

    clipped = _clip(subtractor, _object_order(3000.0))

    assert subtractor.stopSubtraction is False
    assert iterations == [
        "iteration 1 - 1.8% clipped",
        "iteration 2 - 2.2% clipped",
        "iteration 3 - 2.7% clipped",
        "iteration 4 - 3.0% clipped",
        "iteration 5 - 3.2% clipped",
        "iteration 6 - 3.4% clipped",
        "iteration 7 - 3.5% clipped",
        "iteration 8 - 3.7% clipped",
    ]
    assert int(clipped["flagged_object_clipped"].sum()) == 74
    assert int(clipped["flagged_all_clipped"].sum()) == 74
    assert _messages(log, "print", "warning") == [("print", "\tORDER 10: 74 pixels clipped in total = 3.7%)")]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(820.1670639502325, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 83
    assert globalSigma.iloc[SAMPLED_ROWS].tolist() == _approx_or_nan(
        [0.5072931128660357, np.nan, 0.22461213396868954, np.nan]
    )
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(5924705.139971411, rel=1e-12, abs=0)


def test_more_than_85_percent_clipped_logs_zero_percent_but_keeps_pixels_excluded(log: Any) -> None:
    """A negative sigma clips 99.1% of pixels; the object flags reset but the exclusions remain."""
    # SUSPICIOUS (WRONG SCIENCE): "CLIPPING 0% INSTEAD" LEAVES 1982 PIXELS EXCLUDED, FILED AS DY-595
    subtractor, iterations = _subtractor(log, arm="UVB")

    clipped = _clip(subtractor, _object_order(100.0), sigma=-3.0, iterations=3)

    assert iterations == ["iteration 1 - 99.1% clipped"]
    assert int(clipped["flagged_object_clipped"].sum()) == 0
    assert int(clipped["flagged_all_clipped"].sum()) == 1982
    assert _messages(log, "print", "warning") == [
        ("print", "\tORDER 10: 1982 pixels clipped in total = 99.1%)"),
        ("warning", "ORDER 10: More than 85% of pixels flagged to be clipped (99.1%). Clipping 0% instead."),
    ]
    globalSigma = clipped["residual_global_sigma"]
    assert float(globalSigma.sum()) == pytest.approx(7.526984597702706, rel=1e-12, abs=0)
    assert int(globalSigma.isna().sum()) == 1991
    assert float(clipped["flux_percentile_smoothed"].sum()) == pytest.approx(860.4901487716307, rel=1e-12, abs=0)
