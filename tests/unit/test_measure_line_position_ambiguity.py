"""Single-source arc-line selection measures from the real first guess and rejects ambiguous stamps (DY-2207)."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData
from astropy.table import Table

from soxspipe.commonutils.create_dispersion_map import (
    AMBIGUITY_DISTANCE_RATIO,
    create_dispersion_map,
    measure_line_position,
)

pytestmark = pytest.mark.unit

FRAME_SHAPE = (40, 60)
LINE_Y = 20.0
LINE_SIGMA_PX = 0.9
NOISE_ADU = 1.0


def _arc_frame(lineXs: list[float], peaks: list[float] | None = None) -> np.ndarray:
    """Return a noisy frame with one Gaussian arc line per x position, all at `LINE_Y`."""
    peaks = peaks or [1000.0] * len(lineXs)
    y, x = np.mgrid[: FRAME_SHAPE[0], : FRAME_SHAPE[1]]
    frame = np.random.default_rng(7).normal(100.0, NOISE_ADU, FRAME_SHAPE)
    for lineX, peak in zip(lineXs, peaks, strict=True):
        frame += peak * np.exp(-((x - lineX) ** 2 + (y - LINE_Y) ** 2) / (2 * LINE_SIGMA_PX**2))
    return frame


def _stamp_at(frame: np.ndarray, guessX: float, guessY: float, windowHalf: int) -> tuple:
    """Cut a stamp around a first guess the way `detect_pinhole_arc_lines` does."""
    xlow = max(int(np.round(guessX - windowHalf)), 0)
    xup = min(int(np.round(guessX + windowHalf)), frame.shape[1])
    ylow = max(int(np.round(guessY - windowHalf)), 0)
    yup = min(int(np.round(guessY + windowHalf)), frame.shape[0])
    data = frame[ylow:yup, xlow:xup]
    stamp = CCDData(data, unit=u.adu, mask=np.zeros_like(data, dtype=bool))
    return (stamp, xlow, xup, ylow, yup, guessX, guessY)


def _measure(stampInfo: tuple, log: object, windowHalf: int) -> list[dict]:
    """Run the final single-source pass on one stamp."""
    return measure_line_position(
        stampInfo,
        log,
        windowHalf=windowHalf,
        iraf=False,
        sigmaLimit=3,
        iteration=2,
        exclude_border=False,
        returnAll=False,
    )


@pytest.mark.parametrize("guessX", [28.40, 28.60, 28.80])
def test_first_guesses_either_side_of_a_rounding_boundary_get_the_same_answer(log: object, guessX: float) -> None:
    """Two lines 3.9 px apart with the guess between them are ambiguous whatever the stamp origin rounds to."""
    # ARRANGE
    frame = _arc_frame([26.55, 30.45])

    # ACT
    result = _measure(_stamp_at(frame, guessX, LINE_Y, windowHalf=5), log, windowHalf=5)

    # ASSERT
    assert len(result) == 1
    assert np.isnan(result[0]["observed_x"])
    assert np.isnan(result[0]["observed_y"])


def test_a_stamp_clipped_at_the_frame_edge_measures_from_the_real_first_guess(log: object) -> None:
    """At the frame edge the stamp centre is far from the guess, so the line at the guess must still win."""
    # ARRANGE
    frame = _arc_frame([3.0, 10.5])
    stampInfo = _stamp_at(frame, 3.5, LINE_Y, windowHalf=8)

    # ACT
    result = _measure(stampInfo, log, windowHalf=8)

    # ASSERT
    assert stampInfo[1] == 0
    assert len(result) == 1
    assert result[0]["observed_x"] == pytest.approx(3.0, abs=0.2)
    assert result[0]["observed_y"] == pytest.approx(LINE_Y, abs=0.2)


@pytest.mark.parametrize(
    ("guessX", "isAmbiguous"),
    [
        # D1 = 1.2, D2 = 2.7: RATIO 2.25 IS ABOVE THE CONSTANT
        (28.18, False),
        # D1 = 1.5, D2 = 2.4: RATIO 1.6 IS ABOVE THE CONSTANT, SO THE NEAREST LINE IS KEPT
        (28.48, False),
        # D1 = 1.64, D2 = 2.26: RATIO 1.38 IS BELOW THE CONSTANT
        (28.62, True),
    ],
)
def test_two_sources_are_ambiguous_only_below_the_distance_ratio(
    log: object,
    guessX: float,
    isAmbiguous: bool,
) -> None:
    """The nearest source is kept unless the second-nearest is within the ratio of its distance."""
    assert AMBIGUITY_DISTANCE_RATIO == 1.5
    # ARRANGE
    frame = _arc_frame([26.98, 30.88])
    nearest = abs(guessX - 26.98)
    secondNearest = abs(guessX - 30.88)
    assert (secondNearest / nearest < AMBIGUITY_DISTANCE_RATIO) is isAmbiguous

    # ACT
    result = _measure(_stamp_at(frame, guessX, LINE_Y, windowHalf=5), log, windowHalf=5)

    # ASSERT
    assert len(result) == 1
    assert result[0]["ambiguous"] is isAmbiguous
    if isAmbiguous:
        assert np.isnan(result[0]["observed_x"])
        assert np.isnan(result[0]["observed_y"])
    else:
        assert result[0]["observed_x"] == pytest.approx(26.98, abs=0.2)


def test_one_source_is_kept_wherever_the_first_guess_lies_in_the_stamp(log: object) -> None:
    """A lone line is kept as before, even when the guess is off-centre."""
    # ARRANGE
    frame = _arc_frame([30.0])

    # ACT
    result = _measure(_stamp_at(frame, 27.6, LINE_Y + 1.0, windowHalf=5), log, windowHalf=5)

    # ASSERT
    assert len(result) == 1
    assert result[0]["observed_x"] == pytest.approx(30.0, abs=0.2)
    assert result[0]["observed_y"] == pytest.approx(LINE_Y, abs=0.2)
    assert result[0]["flux"] > 0


def test_a_source_exactly_at_the_first_guess_is_never_ambiguous(log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    """With d1 = 0 any second source is farther by an infinite ratio, so the source at the guess is kept."""
    # ARRANGE
    photutilsModule = importlib.import_module("photutils")
    sources = Table(
        {
            "xcentroid": [2.5, 2.0],
            "ycentroid": [2.0, 2.0],
            "sharpness": [0.1, 0.2],
            "roundness1": [0.1, 0.2],
            "roundness2": [0.1, 0.2],
            "npix": [4, 5],
            "sky": [0.0, 0.0],
            "peak": [10.0, 100.0],
            "flux": [10.0, 100.0],
        }
    )
    monkeypatch.setattr(photutilsModule, "DAOStarFinder", lambda **_: lambda *args, **kwargs: sources)
    stamp = CCDData(np.arange(36.0).reshape(6, 6), unit=u.adu)

    # ACT
    result = _measure((stamp, 4, 10, 8, 14, 6.0, 10.0), log, windowHalf=2)

    # ASSERT
    assert len(result) == 1
    assert result[0]["observed_x"] == pytest.approx(6.0)
    assert result[0]["observed_y"] == pytest.approx(10.0)


def test_detect_pinhole_arc_lines_rejects_ambiguous_stamps_as_non_detections_and_counts_them(log: object) -> None:
    """The line table gets NaN for an ambiguous stamp, a position for a clear one, and the pass logs the count."""
    # ARRANGE
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.pinholeFrameMasked = np.ma.array(_arc_frame([26.98, 30.88]), mask=np.zeros(FRAME_SHAPE, dtype=bool))
    mapper.windowHalf = 5
    mapper.debug = False
    inputLines = pd.DataFrame(
        {
            "detector_x": [0.0, 0.0],
            "detector_y": [0.0, 0.0],
            "detector_x_shifted": [28.60, 27.20],
            "detector_y_shifted": [LINE_Y, LINE_Y],
        }
    )

    # ACT
    result = mapper.detect_pinhole_arc_lines(
        inputLines,
        iraf=False,
        sigmaLimit=3,
        iteration=2,
        brightest=False,
        exclude_border=False,
        returnAll=False,
    )

    # ASSERT
    assert np.isnan(result["observed_x"][0][0])
    assert result["observed_x"][1][0] == pytest.approx(26.98, abs=0.2)
    assert ("debug", "1 of 2 arc-line stamps rejected as ambiguous (iteration 2)") in log.messages


def test_brightest_mode_gives_no_detection_when_no_source_has_a_positive_peak(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The brightest path reports a plain non-detection, never an ambiguous stamp."""
    # ARRANGE
    photutilsModule = importlib.import_module("photutils")
    sources = Table(
        {
            "xcentroid": [1.0, 3.0],
            "ycentroid": [1.0, 3.0],
            "sharpness": [0.1, 0.2],
            "roundness1": [0.1, 0.2],
            "roundness2": [0.1, 0.2],
            "npix": [4, 5],
            "sky": [0.0, 0.0],
            "peak": [-1.0, -2.0],
            "flux": [10.0, 100.0],
        }
    )
    monkeypatch.setattr(photutilsModule, "DAOStarFinder", lambda **_: lambda *args, **kwargs: sources)
    stamp = CCDData(np.arange(36.0).reshape(6, 6), unit=u.adu)

    # ACT
    result = measure_line_position(
        (stamp, 4, 10, 8, 14, 6.0, 10.0),
        log,
        windowHalf=2,
        iraf=False,
        sigmaLimit=3,
        iteration=2,
        brightest=True,
        returnAll=False,
    )

    # ASSERT
    assert result == []
