"""The order merge drops a pixel whose variance is far above its neighbours' before it resamples (DY-2078).

The merge interpolates flux linearly, so one pixel with a variance millions of times the local median still carries its
full flux into the merged spectrum and leaks into the grid points either side of it. The pixel is dropped from the flux,
variance and sky columns alike, and the number dropped is logged. The factor is the recipe setting
``merge-variance-clip-factor``.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.horne_extraction import horne_extraction

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
PIXEL_COUNT = 400
SPIKE_INDEX = 200
SPIKE_FLUX = -2.26e6
SPIKE_VARIANCE = 4.9e12
SPIKE_SKY = 5.0e5
SETTING = "merge-variance-clip-factor"


def _extractor(log: object, arm: str = "NIR", recipeSettings: dict | None = None) -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.arm = arm
    extractor.kw = lambda key: key
    extractor.notFlattened = ""
    extractor.recipeSettings = {} if recipeSettings is None else recipeSettings
    return extractor


def _orders(varianceScale: float = 1.0) -> pd.DataFrame:
    """One NIR order of straight-line flux and sky, so a dropped pixel is recovered exactly by interpolation."""
    index = np.arange(PIXEL_COUNT)
    return pd.DataFrame(
        {
            "order": 10,
            "wavelengthMean": 1000.0 + 0.03 * index,
            "pixelScaleNm": 0.03,
            "extractedFluxOptimal": 1000.0 + 2.0 * index,
            "extractedFluxBoxcar": 1000.0 + 2.0 * index,
            "extractedFluxBoxcarRobust": 1000.0 + 2.0 * index,
            "varianceSpectrum": 1000.0 * varianceScale,
            "skyFlux": 50.0 + 0.1 * index,
        }
    )


def _with_spike(orders: pd.DataFrame, index: int = SPIKE_INDEX, variance: float = SPIKE_VARIANCE) -> pd.DataFrame:
    spiked = orders.copy()
    spiked.loc[index, ["extractedFluxOptimal", "extractedFluxBoxcar", "extractedFluxBoxcarRobust"]] = SPIKE_FLUX
    spiked.loc[index, "varianceSpectrum"] = variance
    spiked.loc[index, "skyFlux"] = SPIKE_SKY
    return spiked


def _columns(merged: pd.DataFrame) -> dict[str, np.ndarray]:
    return {name: merged[name].values.value for name in ("WAVE", "FLUX_COUNTS", "VARIANCE", "SKY_COUNTS")}


def test_one_huge_variance_pixel_does_not_change_the_merged_spectrum(log: object) -> None:
    # ARRANGE: THE SAME ORDER WITH AND WITHOUT THE PIXEL. THE FLUX IS A STRAIGHT LINE, SO THE MERGE ACROSS THE GAP MUST
    # GIVE THE CLEAN VALUES
    clean = _orders()
    spiked = _with_spike(clean)

    # ACT
    mergedClean, _ = _extractor(log).merge_extracted_orders(clean.copy())
    mergedSpiked, _ = _extractor(log).merge_extracted_orders(spiked)

    # ASSERT
    cleanColumns, spikedColumns = _columns(mergedClean), _columns(mergedSpiked)
    nearPixel = np.abs(cleanColumns["WAVE"] - (1000.0 + 0.03 * SPIKE_INDEX)) < 0.2
    assert nearPixel.sum() >= 3
    np.testing.assert_allclose(
        spikedColumns["FLUX_COUNTS"][nearPixel], cleanColumns["FLUX_COUNTS"][nearPixel], rtol=1e-5
    )
    for name, values in cleanColumns.items():
        np.testing.assert_allclose(spikedColumns[name], values, rtol=1e-5, err_msg=name)


def test_the_dropped_pixel_is_removed_from_flux_variance_and_sky_alike(log: object) -> None:
    # ARRANGE: THE SPIKE CARRIES A HUGE SKY VALUE TOO. THE REFERENCE IS THE TABLE WITH THAT ROW LEFT OUT
    spiked = _with_spike(_orders())
    withoutRow = spiked.drop(index=SPIKE_INDEX).reset_index(drop=True)

    # ACT
    merged, _ = _extractor(log).merge_extracted_orders(spiked)
    reference, _ = _extractor(log).merge_extracted_orders(withoutRow)

    # ASSERT
    for name, values in _columns(reference).items():
        np.testing.assert_array_equal(_columns(merged)[name], values, err_msg=name)


def test_an_ordinary_variance_range_gives_a_bit_identical_merged_spectrum(log: object) -> None:
    # ARRANGE: VARIANCE RISES BY A FACTOR OF 500 ALONG THE ORDER AND HAS A 50-FOLD BUMP, AS A SKYLINE DOES
    orders = _orders()
    orders["varianceSpectrum"] = 1000.0 * (1.0 + 499.0 * np.linspace(0.0, 1.0, PIXEL_COUNT) ** 3)
    orders.loc[100:103, "varianceSpectrum"] *= 50.0
    wave = orders["wavelengthMean"].values
    expected = {
        "FLUX_COUNTS": orders["extractedFluxOptimal"].values,
        "VARIANCE": orders["varianceSpectrum"].values,
        "SKY_COUNTS": orders["skyFlux"].values,
    }

    # ACT
    merged, _ = _extractor(log).merge_extracted_orders(orders.copy())

    # ASSERT: EVERY PIXEL IS KEPT, SO EACH COLUMN IS THE PLAIN LINEAR RESAMPLE OF THE INPUT
    columns = _columns(merged)
    for name, values in expected.items():
        resampled = np.interp(columns["WAVE"], wave, values, left=0.0, right=0.0).astype(np.float32)
        np.testing.assert_array_equal(columns[name], resampled, err_msg=name)


def test_the_merge_logs_how_many_pixels_it_dropped(log: object) -> None:
    # ARRANGE
    spiked = _with_spike(_orders())

    # ACT
    _extractor(log).merge_extracted_orders(spiked)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert any(f"1 of {PIXEL_COUNT} pixels" in message for message in printed)


def test_the_merge_logs_zero_when_no_pixel_is_dropped(log: object) -> None:
    # ACT
    _extractor(log).merge_extracted_orders(_orders())

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert any(f"0 of {PIXEL_COUNT} pixels" in message for message in printed)


def test_a_smaller_factor_in_the_recipe_settings_drops_a_milder_pixel(log: object) -> None:
    # ARRANGE: A PIXEL WITH 5000 TIMES THE LOCAL MEDIAN VARIANCE IS BELOW THE DEFAULT FACTOR
    milder = _with_spike(_orders(), variance=5.0e6)

    # ACT
    keptByDefault, _ = _extractor(log).merge_extracted_orders(milder.copy())
    droppedByFactor, _ = _extractor(log, recipeSettings={SETTING: 1000.0}).merge_extracted_orders(milder.copy())

    # ASSERT
    nearPixel = np.abs(_columns(keptByDefault)["WAVE"] - (1000.0 + 0.03 * SPIKE_INDEX)) < 0.1
    assert _columns(keptByDefault)["FLUX_COUNTS"][nearPixel].min() < -1.0e5
    assert _columns(droppedByFactor)["FLUX_COUNTS"][nearPixel].min() > 0.0


def test_an_infinite_factor_switches_the_mask_off(log: object) -> None:
    # ACT
    merged, _ = _extractor(log, recipeSettings={SETTING: float("inf")}).merge_extracted_orders(_with_spike(_orders()))

    # ASSERT
    assert _columns(merged)["FLUX_COUNTS"].min() < -1.0e5


@pytest.mark.parametrize("badFactor", [1.0, 0.0, -5.0, "1e6", True, None, float("nan")])
def test_a_factor_that_is_not_a_number_above_one_stops_the_merge(log: object, badFactor: object) -> None:
    # ACT AND ASSERT
    with pytest.raises(ValueError, match=SETTING):
        _extractor(log, recipeSettings={SETTING: badFactor}).merge_extracted_orders(_orders())


def test_a_dropped_end_pixel_does_not_leave_a_zero_bin_at_the_end_of_the_spectrum(log: object) -> None:
    # ARRANGE: THE LAST PIXEL OF THE ORDER IS THE SPIKE
    spiked = _with_spike(_orders(), index=PIXEL_COUNT - 1)

    # ACT
    merged, _ = _extractor(log).merge_extracted_orders(spiked)

    # ASSERT: NO GRID POINT LIES BEYOND THE LAST KEPT PIXEL, AND NO BIN IS ZERO
    columns = _columns(merged)
    assert columns["WAVE"].max() <= 1000.0 + 0.03 * (PIXEL_COUNT - 2)
    assert (columns["VARIANCE"] > 0.0).all()


def test_the_flux_density_of_the_not_flattened_merge_skips_the_dropped_pixel(log: object) -> None:
    # ARRANGE
    spiked = _with_spike(_orders())
    withoutRow = spiked.drop(index=SPIKE_INDEX).reset_index(drop=True)
    notFlat, reference = _extractor(log), _extractor(log)
    notFlat.notFlattened = reference.notFlattened = "_NOTFLAT"

    # ACT
    merged, _ = notFlat.merge_extracted_orders(spiked)
    expected, _ = reference.merge_extracted_orders(withoutRow)

    # ASSERT
    np.testing.assert_array_equal(
        merged["FLUX_DENSITY_COUNTS"].values.value, expected["FLUX_DENSITY_COUNTS"].values.value
    )


def test_the_dropped_pixel_stays_in_the_extracted_orders_table(log: object) -> None:
    # ARRANGE
    spiked = _with_spike(_orders())

    # ACT
    _extractor(log).merge_extracted_orders(spiked)

    # ASSERT: THE EXTRACTION IS NOT CHANGED, ONLY WHAT THE MERGE RESAMPLES
    assert len(spiked) == PIXEL_COUNT
    assert (spiked["varianceSpectrum"] == SPIKE_VARIANCE).sum() == 1


@pytest.mark.parametrize(
    ("fileName", "recipes"),
    [
        (
            "soxs_default_settings.yaml",
            ["soxs-stare", "soxs-stare-std", "soxs-nod", "soxs-nod-std", "soxs-offset", "soxs-offset-std"],
        ),
        ("xsh_default_settings.yaml", ["soxs-stare", "soxs-nod"]),
    ],
)
def test_every_extracting_recipe_block_ships_the_merge_variance_factor(fileName: str, recipes: list[str]) -> None:
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ASSERT
    for recipe in recipes:
        assert settings[recipe][SETTING] == 1.0e6, recipe
