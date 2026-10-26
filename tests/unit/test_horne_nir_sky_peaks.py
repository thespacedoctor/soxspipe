"""Contracts for the NIR sky-line peak detector, which works on the raw sky flux (DY-1968)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.horne_extraction import NIR_SKY_PEAK_MIN_SEPARATION_PX, horne_extraction
from soxspipe.commonutils.nir_sky_refit import _seeded_match, nir_order_sky_lines, order_geometry
from tests.unit._horne_nir_scene import DISPERSION, ORDERS, build_extractor, build_orders

pytestmark = pytest.mark.unit

FIXTURE_PATH = Path(__file__).parent / "data" / "nod024342_o15_sky.npz"
PATTERN_AMPLITUDE = 1.0
LINE_SIGMA_PX = 1.5
SAMPLES = 400
# ALTERNATING +/- AMPLITUDE IS DETERMINISTIC, AND UNLIKE WHITE NOISE IT NEVER MAKES A PEAK PASS THE THRESHOLD.
# ITS FIRST DIFFERENCES ARE +/- 2 * AMPLITUDE, SO THE DETECTOR ESTIMATES THE NOISE AS SQRT(2) * AMPLITUDE, NOT AMPLITUDE
NOISE_PATTERN = np.where(np.arange(SAMPLES) % 2 == 0, PATTERN_AMPLITUDE, -PATTERN_AMPLITUDE)
LINE_CENTRES = (60, 140, 220, 330)
# A MASKED GAP BETWEEN THE LINES AT 140 AND 220
GAP_START = 160
GAP_END = 180
# WHITE NOISE OF UNIT SIGMA ON A SLOPE; LINE HEIGHTS ARE IN UNITS OF THAT SIGMA
GAUSSIAN_SEED = 7
GAUSSIAN_SLOPE = 0.5
GAUSSIAN_LINES = {60: 12.0, 140: 6.0, 220: 2.0, 330: 12.0}
# THE ORDER 15 SKY OF NOD FRAME 024342 HAS 14 REFERENCE LINES, MEASURED AT HALF MAXIMUM.
# THIS DETECTOR FINDS 21 PEAKS ON IT; THE 7 EXTRA PEAKS WERE NOT IDENTIFIED.
# THE FLOOR OF 15 IS A SANITY BOUND ON THE DETECTOR'S OWN COUNT: ABOVE THE 14 REFERENCE LINES, BELOW ITS 21.
# IT DOES NOT PROVE THAT EACH REFERENCE LINE IS FOUND; THE CATALOGUE MATCH TEST BELOW COVERS THAT.
# THE OLD DETECTOR KEPT 1 TO 5 ARTEFACTS
MIN_O15_PEAKS = 15
MAX_OLD_O15_PEAKS = 5
MIN_O15_GLOBAL_SEED_MATCHES = 10


def _detector() -> horne_extraction:
    return object.__new__(horne_extraction)


def _spectrum(
    lineHeights: dict[int, float], slope: float = 0.0, offset: float = 100.0
) -> tuple[pd.Series, pd.Series, pd.Series]:
    pixels = np.arange(SAMPLES, dtype=float)
    sky = offset + slope * pixels + NOISE_PATTERN
    for centre, height in lineHeights.items():
        sky = sky + height * np.exp(-0.5 * ((pixels - centre) / LINE_SIGMA_PX) ** 2)
    wave = pd.Series(1000.0 + 0.05 * pixels)
    return wave, pd.Series(sky), pd.Series(np.ones(SAMPLES, dtype=bool))


def _detect(wave: pd.Series, sky: pd.Series, valid: pd.Series):
    return _detector()._detect_nir_sky_peaks(wave, sky, valid)


def test_finds_lines_standing_above_the_prominence_threshold() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum(dict.fromkeys(LINE_CENTRES, 12.0 * PATTERN_AMPLITUDE))

    # ACT
    skyVals, peaks, waveVals = _detect(wave, sky, valid)

    # ASSERT
    assert sorted(peaks.tolist()) == list(LINE_CENTRES)
    np.testing.assert_allclose(skyVals, sky.to_numpy())
    np.testing.assert_allclose(waveVals, wave.to_numpy())


def test_ignores_bumps_below_the_prominence_threshold() -> None:
    # ARRANGE
    faint = 1.0 * PATTERN_AMPLITUDE
    wave, sky, valid = _spectrum({60: 12.0 * PATTERN_AMPLITUDE, 140: faint, 220: faint})

    # ACT
    _, peaks, _ = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.tolist() == [60]


def test_ignores_a_smooth_continuum_slope() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum({}, slope=5.0)

    # ACT
    _, peaks, _ = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.size == 0


def test_returns_no_peaks_for_a_constant_sky() -> None:
    # ARRANGE
    wave = pd.Series(1000.0 + 0.05 * np.arange(SAMPLES))
    sky = pd.Series(np.full(SAMPLES, 50.0))

    # ACT
    _, peaks, _ = _detect(wave, sky, pd.Series(np.ones(SAMPLES, dtype=bool)))

    # ASSERT
    assert peaks.size == 0


def test_returns_no_peaks_when_the_valid_samples_are_all_nan() -> None:
    # ARRANGE
    wave = pd.Series(1000.0 + 0.05 * np.arange(SAMPLES))
    sky = pd.Series(np.full(SAMPLES, np.nan))

    # ACT
    skyVals, peaks, waveVals = _detect(wave, sky, pd.Series(np.ones(SAMPLES, dtype=bool)))

    # ASSERT
    assert peaks.size == 0
    assert skyVals.size == waveVals.size == SAMPLES


def test_returns_no_peaks_for_a_single_valid_sample() -> None:
    # ARRANGE
    wave, sky, _ = _spectrum({})
    valid = pd.Series(np.zeros(SAMPLES, dtype=bool))
    valid.iloc[10] = True

    # ACT
    _, peaks, waveVals = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.size == 0
    assert waveVals.size == 1


def test_uses_only_the_valid_samples() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum(dict.fromkeys(LINE_CENTRES, 12.0 * PATTERN_AMPLITUDE))
    valid.iloc[:50] = False
    sky.iloc[:50] = np.nan

    # ACT
    _, peaks, waveVals = _detect(wave, sky, valid)

    # ASSERT
    assert waveVals.size == SAMPLES - 50
    assert sorted(peaks.tolist()) == [centre - 50 for centre in LINE_CENTRES]


def test_threshold_follows_the_difference_noise_of_gaussian_sky() -> None:
    # ARRANGE
    trueSigma = 1.0
    noise = np.random.default_rng(GAUSSIAN_SEED).normal(0.0, trueSigma, SAMPLES)
    sky = 100.0 + GAUSSIAN_SLOPE * np.arange(SAMPLES) + noise
    for centre, height in GAUSSIAN_LINES.items():
        sky[centre] += height * trueSigma
    wave = pd.Series(1000.0 + 0.05 * np.arange(SAMPLES))
    valid = pd.Series(np.ones(SAMPLES, dtype=bool))

    # ACT
    _, peaks, _ = _detect(wave, pd.Series(sky), valid)

    # ASSERT
    # THE THRESHOLD IS 5 * TRUE SIGMA. NO DIVISION BY SQRT(2) WOULD RAISE IT TO 7.1 AND LOSE THE 6 SIGMA LINE.
    # THE STD OF THE SLOPED SKY ITSELF WOULD RAISE IT FAR ABOVE EVERY LINE.
    assert peaks.tolist() == [60, 140, 330]


def test_keeps_only_the_higher_of_two_lines_closer_than_the_minimum_separation() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum({})
    sky.iloc[100] += 20.0
    sky.iloc[100 + NIR_SKY_PEAK_MIN_SEPARATION_PX - 2] += 15.0

    # ACT
    _, peaks, _ = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.tolist() == [100]


def test_keeps_two_lines_separated_by_exactly_the_minimum_separation() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum({})
    sky.iloc[100] += 20.0
    sky.iloc[100 + NIR_SKY_PEAK_MIN_SEPARATION_PX] += 15.0

    # ACT
    _, peaks, _ = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.tolist() == [100, 100 + NIR_SKY_PEAK_MIN_SEPARATION_PX]


def test_maps_peaks_to_the_right_wavelengths_across_an_interior_masked_gap() -> None:
    # ARRANGE
    wave, sky, valid = _spectrum(dict.fromkeys(LINE_CENTRES, 12.0 * PATTERN_AMPLITUDE))
    valid.iloc[GAP_START:GAP_END] = False
    sky.iloc[GAP_START:GAP_END] = np.nan

    # ACT
    skyVals, peaks, waveVals = _detect(wave, sky, valid)

    # ASSERT
    assert skyVals.size == waveVals.size == SAMPLES - (GAP_END - GAP_START)
    np.testing.assert_allclose(waveVals[peaks], wave.iloc[list(LINE_CENTRES)].to_numpy())


def test_nir_order_lines_are_detected_without_the_shared_continuum_detector(log, tmp_path, monkeypatch) -> None:
    # ARRANGE
    extractor = build_extractor(log, tmp_path)

    def fail(*args, **kwargs):
        raise AssertionError("the NIR global fit must not use the shared peak detector")

    monkeypatch.setattr(extractor, "_detect_sky_peaks", fail)

    # ACT
    geometry, orderLines = extractor._nir_order_lines(build_orders(), ORDERS, dict.fromkeys(ORDERS, DISPERSION))

    # ASSERT
    assert sorted(orderLines) == ORDERS
    assert sorted(geometry) == ORDERS


# REAL ORDER 15 SKY OF NOD FRAME 024342


@pytest.fixture(scope="module")
def o15() -> dict:
    with np.load(FIXTURE_PATH, allow_pickle=False) as fixture:
        return {key: fixture[key] for key in fixture.files}


def _o15_series(o15: dict) -> tuple[pd.Series, pd.Series, pd.Series]:
    return pd.Series(o15["wave"]), pd.Series(o15["sky"]), pd.Series(o15["valid"])


def test_finds_the_real_o15_sky_lines(o15: dict) -> None:
    # ARRANGE
    wave, sky, valid = _o15_series(o15)

    # ACT
    _, peaks, _ = _detect(wave, sky, valid)

    # ASSERT
    assert peaks.size >= MIN_O15_PEAKS


def test_the_shared_detector_misses_the_real_o15_sky_lines(o15: dict) -> None:
    # ARRANGE
    wave, sky, valid = _o15_series(o15)
    objectFlux = pd.Series(np.zeros(wave.size))

    # ACT
    _, _, peaks, _, _ = _detector()._detect_sky_peaks(wave, sky, objectFlux, valid)

    # ASSERT
    assert peaks.size <= MAX_OLD_O15_PEAKS


def test_the_real_o15_lines_match_the_catalogue_about_the_global_seed(o15: dict) -> None:
    # ARRANGE
    wave, sky, valid = _o15_series(o15)
    detector = _detector()
    skyVals, peaks, waveVals = _detect(wave, sky, valid)
    peakNm = np.asarray(detector._measure_peak_centres(waveVals, skyVals, peaks), dtype=float)
    lines = nir_order_sky_lines(
        order=15,
        peakNm=peakNm,
        catalogueNm=o15["catalogueNm"],
        pixelScaleNm=np.interp(peakNm, waveVals, o15["pixelScaleNm"]),
        toleranceNm=float(o15["toleranceNm"]),
    )
    geometry = order_geometry(*(float(x) for x in o15["geometry"]))

    # ACT
    match = _seeded_match(lines, geometry, float(o15["globalSeedPx"]), float(o15["ownSeedPx"]))

    # ASSERT
    assert match.globalLines >= MIN_O15_GLOBAL_SEED_MATCHES
    assert not match.keepsOwnSeed
