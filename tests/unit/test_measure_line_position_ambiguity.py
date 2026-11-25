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


def _measure(stampInfo: tuple, log: object, windowHalf: int, **kwargs: object) -> list[dict]:
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
        **kwargs,
    )


def _bare_mapper(log: object, **attributes: object) -> create_dispersion_map:
    """Return a `create_dispersion_map` with no constructor run and only the given attributes set."""
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    for name, value in attributes.items():
        setattr(mapper, name, value)
    return mapper


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


def test_ambiguity_distance_ratio_is_1_5() -> None:
    """Dave ruled 1.5 on DY-2207 after the real-data A/B (2.0 dropped too many good lines)."""
    # ASSERT
    assert AMBIGUITY_DISTANCE_RATIO == 1.5


@pytest.mark.parametrize(
    "guessX",
    [
        # D1 = 1.2, D2 = 2.7: RATIO 2.25
        28.18,
        # D1 = 1.5, D2 = 2.4: RATIO 1.6, ABOVE 1.5 BUT BELOW THE EARLIER 2.0
        28.48,
    ],
)
def test_the_nearest_source_is_kept_when_the_second_is_at_least_the_ratio_farther(log: object, guessX: float) -> None:
    """The nearest source is kept when the second-nearest lies at least the ratio times farther from the guess."""
    # ARRANGE
    frame = _arc_frame([26.98, 30.88])

    # ACT
    result = _measure(_stamp_at(frame, guessX, LINE_Y, windowHalf=5), log, windowHalf=5)

    # ASSERT
    assert len(result) == 1
    assert result[0]["ambiguous"] is False
    assert result[0]["observed_x"] == pytest.approx(26.98, abs=0.2)


def test_two_sources_below_the_distance_ratio_give_no_detection(log: object) -> None:
    """With D1 = 1.64 and D2 = 2.26 (ratio 1.38, below 1.5) the stamp is ambiguous and gives NaN."""
    # ARRANGE
    frame = _arc_frame([26.98, 30.88])

    # ACT
    result = _measure(_stamp_at(frame, 28.62, LINE_Y, windowHalf=5), log, windowHalf=5)

    # ASSERT
    assert len(result) == 1
    assert result[0]["ambiguous"] is True
    assert np.isnan(result[0]["observed_x"])
    assert np.isnan(result[0]["observed_y"])


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
    mapper = _bare_mapper(
        log,
        pinholeFrameMasked=np.ma.array(_arc_frame([26.98, 30.88]), mask=np.zeros(FRAME_SHAPE, dtype=bool)),
        windowHalf=5,
        debug=False,
    )
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


MULTIPIN_LINE_XS = [21.0, 26.0, 31.0, 36.0]
# THE GUESS IS 2.4 PX FROM THE PINHOLE AT 26 AND 2.6 PX FROM THE ONE AT 31: RATIO 1.08, INSIDE THE 1.5 LIMIT
MULTIPIN_GUESS_X = 28.4


def test_multipinhole_stamp_is_ambiguous_by_default(log: object) -> None:
    """Four pinholes 5 px apart with the guess near a midpoint are rejected when the ambiguity test is on."""
    # ARRANGE
    frame = _arc_frame(MULTIPIN_LINE_XS)

    # ACT
    result = _measure(_stamp_at(frame, MULTIPIN_GUESS_X, LINE_Y, windowHalf=10), log, windowHalf=10)

    # ASSERT
    assert len(result) == 1
    assert result[0]["ambiguous"] is True
    assert np.isnan(result[0]["observed_x"])


def test_multipinhole_stamp_keeps_the_nearest_source_when_the_ambiguity_test_is_off(log: object) -> None:
    """With a first-guess map the nearest pinhole to the true first guess is kept, however close the next one is."""
    # ARRANGE
    frame = _arc_frame(MULTIPIN_LINE_XS)

    # ACT
    result = _measure(
        _stamp_at(frame, MULTIPIN_GUESS_X, LINE_Y, windowHalf=10),
        log,
        windowHalf=10,
        shouldTestAmbiguity=False,
    )

    # ASSERT
    assert len(result) == 1
    assert result[0]["ambiguous"] is False
    assert result[0]["observed_x"] == pytest.approx(MULTIPIN_LINE_XS[1], abs=0.2)
    assert result[0]["observed_y"] == pytest.approx(LINE_Y, abs=0.2)


def test_switching_the_ambiguity_test_off_still_measures_from_the_true_first_guess(log: object) -> None:
    """A stamp clipped at the frame edge keeps the source at the guess, not the one at the stamp centre."""
    # ARRANGE
    frame = _arc_frame([3.0, 10.5])
    stampInfo = _stamp_at(frame, 3.5, LINE_Y, windowHalf=8)

    # ACT
    result = _measure(stampInfo, log, windowHalf=8, shouldTestAmbiguity=False)

    # ASSERT
    assert stampInfo[1] == 0
    assert result[0]["observed_x"] == pytest.approx(3.0, abs=0.2)


def test_detect_pinhole_arc_lines_counts_no_ambiguous_stamps_when_the_test_is_off(log: object) -> None:
    """The spatial-solution pass keeps the nearest pinhole of an otherwise ambiguous stamp and logs zero rejections."""
    # ARRANGE
    mapper = _bare_mapper(
        log,
        pinholeFrameMasked=np.ma.array(_arc_frame(MULTIPIN_LINE_XS), mask=np.zeros(FRAME_SHAPE, dtype=bool)),
        windowHalf=10,
        debug=False,
    )
    inputLines = pd.DataFrame(
        {
            "detector_x": [0.0],
            "detector_y": [0.0],
            "detector_x_shifted": [MULTIPIN_GUESS_X],
            "detector_y_shifted": [LINE_Y],
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
        shouldTestAmbiguity=False,
    )

    # ASSERT
    assert result["observed_x"][0][0] == pytest.approx(MULTIPIN_LINE_XS[1], abs=0.2)
    assert ("debug", "0 of 1 arc-line stamps rejected as ambiguous (iteration 2)") in log.messages


class _StopAfterDetectionError(Exception):
    """Raised by the fake detector to end `get()` once every detection pass has been seen."""


DETECTION_PASSES = 3


@pytest.mark.parametrize("tightFit", [False, True])
@pytest.mark.parametrize(
    ("firstGuessMap", "expectedShouldTest"),
    [("first-guess.fits", False), (False, True)],
)
def test_get_tests_ambiguity_on_every_pass_only_without_a_first_guess_map(
    log: object,
    firstGuessMap: object,
    expectedShouldTest: bool,
    tightFit: bool,
) -> None:
    """A first-guess map (spatial solution) turns the ambiguity test off on every pass; without one it stays on."""
    # ARRANGE
    lines = pd.DataFrame(
        {
            "wavelength": [500.0],
            "order": [10],
            "slit_position": [0.0],
            "detector_x": [100.0],
            "detector_y": [20.0],
            "observed_x": [100.2],
            "observed_y": [20.1],
        }
    )
    mapper = _bare_mapper(
        log,
        lineDetectionTable=False,
        firstGuessMap=firstGuessMap,
        recipeSettings={"pinhole-detection-thres-sigma": 3},
        windowSize=10,
        arm="VIS",
        inst="SOXS",
        debug=False,
    )
    seenKwargs: list[dict[str, object]] = []

    def detect(*args: object, **kwargs: object) -> pd.DataFrame:
        seenKwargs.append(kwargs)
        if len(seenKwargs) == DETECTION_PASSES:
            raise _StopAfterDetectionError
        return lines.copy()

    mapper._initialize_recipe_settings = lambda: (False, tightFit, 1, 1, 0)
    mapper.get_predicted_line_list = lambda: lines.copy()
    mapper._prepare_pinhole_frame = lambda: "prepared-frame"
    mapper._generate_quicklook_and_debug_plots = lambda *args: None
    mapper.detect_pinhole_arc_lines = detect
    mapper._explode_multiple_detections = lambda table: table
    mapper._calculate_position_differences = lambda table: table
    mapper._handle_multipin_hole_big_shift = lambda table, *args: table
    mapper._find_and_apply_cluster_shift = lambda table, mask: (table, 0.2, 0.1)
    mapper._calculate_order_shift_statistics = lambda table, mask: (0.2, 0.1, 0.01, 0.01, 0.3, 0.02)

    # ACT
    with pytest.raises(_StopAfterDetectionError):
        mapper.get()

    # ASSERT
    assert [kwargs["iteration"] for kwargs in seenKwargs] == [0, 1, 2]
    assert [kwargs["shouldTestAmbiguity"] for kwargs in seenKwargs] == [expectedShouldTest] * DETECTION_PASSES
