"""The SOXS order merge puts the spectrum on the fixed log-wavelength grid of its arm (DY-2348).

The grid step comes from ``log-wavelength-step-kms`` in ``soxs_detector_parameters.yaml``, per unbinned pixel, times
the binning along the dispersion axis. X-Shooter keeps its linear grid and ``np.interp`` (see
``test_horne_extraction_helpers``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml
from astropy.io import fits

from soxspipe.commonutils.horne_extraction import horne_extraction
from soxspipe.commonutils.log_wavelength_grid import LOG_GRID_ORIGIN_NM, SPEED_OF_LIGHT_KMS
from tests.factories import product_table, soxs_merge_extractor
from tests.unit.test_horne_extraction_helpers import _configure_synthetic_orchestration

pytestmark = pytest.mark.unit

VIS_STEP_KMS = 17.83
NIR_STEP_KMS = 13.16


def _extractor(
    log: object, arm: str = "NIR", cards: dict | None = None, notFlattened: str = "", carrySkyCounts: bool = True
) -> horne_extraction:
    return soxs_merge_extractor(log, arm, cards, notFlattened, carrySkyCounts)


def _order(
    order: int, start: float, count: int, flux: float, scale: float = 0.05, variance: float = 4.0
) -> pd.DataFrame:
    wavelength = start + scale * np.arange(count)
    return pd.DataFrame(
        {
            "order": order,
            "wavelengthMean": wavelength,
            "pixelScaleNm": scale,
            "extractedFluxOptimal": flux,
            "extractedFluxBoxcar": flux,
            "extractedFluxBoxcarRobust": flux,
            "varianceSpectrum": variance,
            "skyFlux": 1.5,
        }
    )


def _wave(merged: pd.DataFrame) -> np.ndarray:
    return merged["WAVE"].values.value


def _column(merged: pd.DataFrame, name: str) -> np.ndarray:
    return merged[name].values.value


def test_nir_merged_wavelengths_are_points_of_the_grid_of_the_measured_nir_step(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))

    # ASSERT: EVERY POINT IS 1000 * exp(k * dv / c) FOR AN INTEGER k, AND NEIGHBOURS DIFFER BY dv
    wave = _wave(merged)
    index = np.log(wave / LOG_GRID_ORIGIN_NM) * SPEED_OF_LIGHT_KMS / NIR_STEP_KMS
    np.testing.assert_allclose(index, np.round(index), atol=1e-8)
    np.testing.assert_allclose(wave[1:] / wave[:-1], np.exp(NIR_STEP_KMS / SPEED_OF_LIGHT_KMS), rtol=1e-12)
    np.testing.assert_allclose(_column(merged, "FLUX_COUNTS"), 50.0, rtol=1e-6)


def test_the_merge_reads_the_step_from_the_detector_parameters_yaml(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")
    extractor.detectorParams = {**extractor.detectorParams, "log-wavelength-step-kms": 30.0}

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))

    # ASSERT
    wave = _wave(merged)
    np.testing.assert_allclose(wave[1:] / wave[:-1], np.exp(30.0 / SPEED_OF_LIGHT_KMS), rtol=1e-12)


def test_the_shipped_yaml_holds_the_measured_steps() -> None:
    # ARRANGE
    import soxspipe

    path = f"{soxspipe.__path__[0]}/resources/soxs_detector_parameters.yaml"

    # ACT
    with open(path) as stream:
        parameters = yaml.safe_load(stream)

    # ASSERT
    assert parameters["VIS"]["log-wavelength-step-kms"] == VIS_STEP_KMS
    assert parameters["NIR"]["log-wavelength-step-kms"] == NIR_STEP_KMS


@pytest.mark.parametrize(
    ("cards", "binning"),
    [
        ({"ESO DET BINX": 1, "ESO DET BINY": 1}, 1),
        ({"ESO DET BINX": 2, "ESO DET BINY": 1}, 1),
        ({"ESO DET BINX": 1, "ESO DET BINY": 2}, 2),
        ({"ESO DET BINX": 2, "ESO DET BINY": 2}, 2),
    ],
    ids=["1x1", "2x1-is-unbinned-along-the-dispersion", "1x2", "2x2"],
)
def test_the_vis_step_is_multiplied_by_biny_and_not_binx(log: object, cards: dict, binning: int) -> None:
    # ARRANGE
    extractor = _extractor(log, "VIS", cards)

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(1, 600.0, 400, 50.0, scale=0.04))

    # ASSERT
    wave = _wave(merged)
    np.testing.assert_allclose(wave[1:] / wave[:-1], np.exp(VIS_STEP_KMS * binning / SPEED_OF_LIGHT_KMS), rtol=1e-12)


def test_a_vis_frame_without_a_biny_card_is_an_error(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "VIS", {"ESO DET BINX": 1})

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="ESO DET BINY"):
        extractor.merge_extracted_orders(_order(1, 600.0, 400, 50.0, scale=0.04))


def test_a_nir_frame_with_no_binning_cards_is_unbinned(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR", {})

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))

    # ASSERT
    wave = _wave(merged)
    np.testing.assert_allclose(wave[1:] / wave[:-1], np.exp(NIR_STEP_KMS / SPEED_OF_LIGHT_KMS), rtol=1e-12)


def test_the_merged_wave_is_not_rounded_to_two_decimals(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))

    # ASSERT: A 13 KM/S BIN IS 0.044 NM WIDE AT 1000 NM, SO MOST POINTS HAVE MORE THAN 2 DECIMALS
    wave = _wave(merged)
    assert (np.abs(wave - np.round(wave, 2)) > 1e-6).any()


def test_the_overlap_of_two_orders_is_cut_at_the_join_with_no_zigzag(log: object) -> None:
    # ARRANGE: THE RED ORDER (11) RUNS 1010 TO 1030 NM AT 100 E-, THE BLUE ORDER (12) 1000 TO 1020 NM AT 200 E-
    extractor = _extractor(log, "NIR")
    orders = pd.concat([_order(11, 1010.0, 400, 100.0), _order(12, 1000.0, 400, 200.0)], ignore_index=True)

    # ACT
    merged, joins = extractor.merge_extracted_orders(orders)

    # ASSERT: ONE HARD STEP FROM 200 TO 100, AT THE JOIN, AND NO BIN HOLDS A MIX
    flux = _column(merged, "FLUX_COUNTS").astype(float).round(3)
    assert set(np.unique(flux)) == {100.0, 200.0}
    assert np.count_nonzero(np.diff(flux)) == 1
    assert list(joins) == ["11-12"]
    assert joins["11-12"] == pytest.approx(1015.0, abs=0.05)
    stepAt = np.flatnonzero(np.diff(flux))[0]
    assert _wave(merged)[stepAt] < joins["11-12"] < _wave(merged)[stepAt + 1]


def test_the_gap_between_orders_is_left_empty_and_not_interpolated(log: object) -> None:
    # ARRANGE: THE BLUE ORDER ENDS AT 1020 NM AND THE RED ORDER STARTS AT 1031.5 NM, AS ORDERS 11 AND 10 DO IN THE NIR
    extractor = _extractor(log, "NIR")
    orders = pd.concat([_order(10, 1031.5, 400, 150.0), _order(11, 1000.0, 400, 0.0)], ignore_index=True)

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT: THE GAP BINS HOLD 0, WHERE LINEAR INTERPOLATION BETWEEN THE ORDER ENDS WOULD GIVE A RAMP FROM 0 TO 150
    wave = _wave(merged)
    flux = _column(merged, "FLUX_COUNTS")
    variance = _column(merged, "VARIANCE")
    inGap = (wave > 1020.1) & (wave < 1031.4)
    assert inGap.sum() > 200
    assert (flux[inGap] == 0.0).all()
    assert (variance[inGap] == 0.0).all()
    assert (_column(merged, "SKY_COUNTS")[inGap] == 0.0).all()
    assert np.isnan(_column(merged, "SNR")[inGap]).all()
    np.testing.assert_allclose(flux[wave > 1031.6], 150.0, rtol=1e-6)


def test_the_merge_keeps_the_extracted_orders_table_untouched(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR", notFlattened="_NOTFLAT")
    orders = pd.concat([_order(11, 1010.0, 400, 100.0), _order(12, 1000.0, 400, 200.0)], ignore_index=True)
    before = orders.copy()

    # ACT
    extractor.merge_extracted_orders(orders)

    # ASSERT
    pd.testing.assert_frame_equal(orders, before)


def test_the_merged_variance_is_the_sum_of_the_squared_weights_times_the_variances(log: object) -> None:
    # ARRANGE: PIXELS 0.05 NM WIDE ON A GRID OF ABOUT 0.044 NM BINS: A BIN TAKES 1 OR 2 PIXELS
    extractor = _extractor(log, "NIR")

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0, variance=9.0))

    # ASSERT: ONE PIXEL GIVES 9, TWO PIXELS WITH WEIGHTS a AND 1 - a GIVE 9 * (a^2 + (1 - a)^2), BETWEEN 4.5 AND 9
    variance = _column(merged, "VARIANCE")
    assert variance.max() <= 9.0 * (1 + 1e-6)
    assert variance.min() >= 4.5 * (1 - 1e-6)
    assert variance.min() < 9.0


def test_snr_is_the_flux_over_the_square_root_of_the_variance(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0, variance=9.0))

    # ASSERT
    np.testing.assert_allclose(
        _column(merged, "SNR"), _column(merged, "FLUX_COUNTS") / np.sqrt(_column(merged, "VARIANCE")), rtol=1e-6
    )


def test_a_pixel_with_a_huge_variance_is_dropped_and_leaves_no_flux_behind(log: object) -> None:
    # ARRANGE: A SPIKE OF -2E6 E- WITH A VARIANCE OF 5E12 IN THE MIDDLE OF A FLAT ORDER (DY-2078)
    extractor = _extractor(log, "NIR")
    orders = _order(10, 1000.0, 400, 1000.0)
    orders.loc[200, ["extractedFluxOptimal", "varianceSpectrum", "skyFlux"]] = [-2.0e6, 5.0e12, 5.0e5]

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT: NO BIN CARRIES THE SPIKE. THE BINS THAT ONLY IT OVERLAPPED HAVE NO DATA AND HOLD 0
    flux = _column(merged, "FLUX_COUNTS")
    assert flux.min() >= 0.0
    assert flux.max() == pytest.approx(1000.0, rel=1e-6)
    assert _column(merged, "VARIANCE").max() < 1.0e4
    assert _column(merged, "SKY_COUNTS").max() < 10.0


def test_sky_counts_stay_nan_when_the_extraction_had_no_sky_frame(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")
    orders = pd.concat([_order(10, 1031.5, 400, 150.0), _order(11, 1000.0, 400, 0.0)], ignore_index=True)
    orders["skyFlux"] = np.nan

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT
    assert np.isnan(_column(merged, "SKY_COUNTS")).all()


def test_no_sky_column_is_made_when_the_extractor_carries_no_sky_counts(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR", carrySkyCounts=False)

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))

    # ASSERT
    assert list(merged.columns) == ["WAVE", "FLUX_COUNTS", "VARIANCE", "SNR"]


@pytest.mark.parametrize(("notFlattened", "hasDensity"), [("", False), ("_NOTFLAT", True)])
def test_flux_density_is_made_only_for_the_not_flattened_spectrum(
    log: object, notFlattened: str, hasDensity: bool
) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR", notFlattened=notFlattened)

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 12.0, scale=0.06))

    # ASSERT
    assert ("FLUX_DENSITY_COUNTS" in merged.columns) is hasDensity


def test_the_flux_density_is_resampled_with_the_same_normalized_weights(log: object) -> None:
    # ARRANGE: 12 E- IN A 0.06 NM PIXEL IS 200 E-/NM EVERYWHERE, SO ANY WEIGHTED MEAN IS 200
    extractor = _extractor(log, "NIR", notFlattened="_NOTFLAT")

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 12.0, scale=0.06))

    # ASSERT
    np.testing.assert_allclose(_column(merged, "FLUX_DENSITY_COUNTS"), 200.0, rtol=1e-6)


def test_the_integrated_flux_of_a_ramp_is_kept_through_the_merge(log: object) -> None:
    # ARRANGE: A RAMP FROM 100 TO 500 E- ACROSS ONE ORDER
    extractor = _extractor(log, "NIR")
    orders = _order(10, 1000.0, 400, 0.0)
    orders["extractedFluxOptimal"] = np.linspace(100.0, 500.0, 400)

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT: THE DENSITY INTEGRATED OVER WAVELENGTH IS THE SAME BEFORE AND AFTER, FOR THE BINS INSIDE THE ORDER
    wave = _wave(merged)
    binWidth = wave * NIR_STEP_KMS / SPEED_OF_LIGHT_KMS
    integrated = (_column(merged, "FLUX_COUNTS") * binWidth).sum()
    lower = wave[0] * np.exp(-0.5 * NIR_STEP_KMS / SPEED_OF_LIGHT_KMS)
    upper = wave[-1] * np.exp(0.5 * NIR_STEP_KMS / SPEED_OF_LIGHT_KMS)
    pixelLower = orders["wavelengthMean"].values - 0.025
    pixelUpper = orders["wavelengthMean"].values + 0.025
    overlap = np.clip(np.minimum(pixelUpper, upper) - np.maximum(pixelLower, lower), 0, None)
    assert integrated == pytest.approx((orders["extractedFluxOptimal"].values * overlap).sum(), rel=1e-4)


def test_an_instrument_other_than_soxs_keeps_the_linear_grid(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")
    extractor.inst = "XSHOOTER"

    # ACT
    merged, _ = extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0, scale=0.06))

    # ASSERT: THE 0.06 NM STEP OF THE X-SHOOTER NIR MERGE
    np.testing.assert_allclose(np.diff(_wave(merged)), 0.06, rtol=1e-9)


def test_the_written_merged_spectrum_keeps_the_exact_grid_and_the_empty_gap_bins(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """Smoke test: the SOXS merge and the real FITS writer, with no writer mocked (DY-2348)."""
    # ARRANGE: THE REAL MERGE OF TWO NIR ORDERS WITH A GAP, AND THE REAL FITS WRITER
    extractor, _, _ = _configure_synthetic_orchestration(log, monkeypatch, products=product_table().iloc[0:0])
    soxsExtractor = _extractor(log, "NIR")
    orders = pd.concat([_order(10, 1031.5, 400, 150.0), _order(11, 1000.0, 400, 80.0)], ignore_index=True)
    extractor.arm = "NIR"
    extractor.detectorParams = soxsExtractor.detectorParams
    extractor.merge_extracted_orders = lambda _: soxsExtractor.merge_extracted_orders(orders)
    extractor.filenameTemplate = "SYNTHETIC.fits"
    extractor.productDir = str(tmp_path)
    extractor.noddingSequence = ""
    extractor.settings = {"instrument": "soxs"}

    # ACT
    _, _, merged, _, _ = extractor.extract()

    # ASSERT
    with fits.open(tmp_path / "SYNTHETIC_SPECTRUM.fits") as hdus:
        table = hdus["SPECTRUM"].data
        tableHeader = hdus["SPECTRUM"].header
    wave = table["WAVE"]
    np.testing.assert_array_equal(wave, merged["WAVE"].values.value)
    index = np.log(wave / LOG_GRID_ORIGIN_NM) * SPEED_OF_LIGHT_KMS / NIR_STEP_KMS
    np.testing.assert_allclose(index, np.round(index), atol=1e-8)
    inGap = (wave > 1020.1) & (wave < 1031.4)
    assert (table["FLUX_COUNTS"][inGap] == 0).all()
    assert np.isnan(table["SNR"][inGap]).all()
    assert tableHeader["EXTNAME"] == "SPECTRUM"


def test_the_soxs_merge_logs_how_many_pixels_the_variance_mask_dropped(log: object) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")
    orders = _order(10, 1000.0, 400, 1000.0)
    orders.loc[200, "varianceSpectrum"] = 5.0e12

    # ACT
    extractor.merge_extracted_orders(orders)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert any("1 of 400 pixels" in message for message in printed)


def test_the_variance_clip_setting_applies_to_the_soxs_merge(log: object) -> None:
    # ARRANGE: A PIXEL WITH 5000 TIMES THE LOCAL MEDIAN VARIANCE IS BELOW THE DEFAULT FACTOR OF 1E6
    orders = _order(10, 1000.0, 400, 1000.0, variance=4000.0)
    orders.loc[200, ["extractedFluxOptimal", "varianceSpectrum"]] = [-2.0e6, 2.0e7]
    strict = _extractor(log, "NIR")
    strict.recipeSettings = {"merge-variance-clip-factor": 1000.0}

    # ACT
    kept, _ = _extractor(log, "NIR").merge_extracted_orders(orders)
    dropped, _ = strict.merge_extracted_orders(orders)

    # ASSERT
    assert _column(kept, "FLUX_COUNTS").min() < -1.0e5
    assert _column(dropped, "FLUX_COUNTS").min() >= 0.0


@pytest.mark.parametrize("badFactor", [1.0, True, None, float("nan")])
def test_a_variance_clip_factor_that_is_not_a_number_above_one_stops_the_soxs_merge(
    log: object, badFactor: object
) -> None:
    # ARRANGE
    extractor = _extractor(log, "NIR")
    extractor.recipeSettings = {"merge-variance-clip-factor": badFactor}

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="merge-variance-clip-factor"):
        extractor.merge_extracted_orders(_order(10, 1000.0, 400, 50.0))


def test_an_empty_extraction_table_stops_the_soxs_merge_on_the_wavelength_column(log: object) -> None:
    # ARRANGE
    empty = _order(10, 1000.0, 3, 50.0).iloc[0:0]

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="wavelengthMean"):
        _extractor(log, "NIR").merge_extracted_orders(empty)


def test_an_order_with_no_pixel_left_changes_nothing_in_the_orders_beside_it(log: object) -> None:
    # ARRANGE: THE MIDDLE ORDER HAS EVERY PIXEL MASKED, SO IT IS ABSENT FROM THE TABLE (DY-2290)
    extractor = _extractor(log, "NIR")
    low = _order(10, 1000.0, 400, 10.0)
    high = _order(12, 1015.0, 400, 30.0)
    empty = _order(11, 1010.0, 400, 20.0).iloc[0:0]

    # ACT
    withEmpty, _ = extractor.merge_extracted_orders(pd.concat([low, empty, high], ignore_index=True))
    without, _ = extractor.merge_extracted_orders(pd.concat([low, high], ignore_index=True))

    # ASSERT: EVERY COLUMN THE FLUX CALIBRATION READS IS FINITE
    for column in ("WAVE", "FLUX_COUNTS", "VARIANCE"):
        assert np.isfinite(_column(withEmpty, column)).all()
        np.testing.assert_array_equal(_column(withEmpty, column), _column(without, column))


def test_the_merged_bins_lie_inside_the_extracted_wavelength_range(log: object) -> None:
    # ARRANGE: THE DATA START 0.004 NM INTO A BIN, SO THE FIRST BIN IS ONLY PART COVERED
    extractor = _extractor(log, "NIR")
    orders = _order(10, 1000.004, 400, 50.0)

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT: ALL THE MERGED BINS (CENTRES) LIE INSIDE THE PIXEL EDGES, AND EVERY BIN HAS DATA AND VARIANCE
    wave = _wave(merged)
    assert wave.min() >= 1000.004 - 0.025
    assert wave.max() <= 1000.004 + 0.05 * 399 + 0.025
    assert (_column(merged, "VARIANCE") > 0).all()


def test_a_constant_flux_density_is_flat_across_a_join_of_orders_with_different_pixel_scales(log: object) -> None:
    # ARRANGE: THE RED ORDER HAS 0.06 NM PIXELS AND THE BLUE ORDER 0.04 NM PIXELS. 200 E-/NM IS 12 E- IN A 0.06 NM
    # PIXEL AND 8 E- IN A 0.04 NM PIXEL, SO THE DENSITY IS THE SAME ON BOTH SIDES OF THE JOIN (DY-2339)
    extractor = _extractor(log, "NIR", notFlattened="_NOTFLAT")
    red = _order(11, 1010.0, 300, 12.0, scale=0.06)
    blue = _order(12, 1000.0, 400, 8.0, scale=0.04)

    # ACT
    merged, joins = extractor.merge_extracted_orders(pd.concat([red, blue], ignore_index=True))

    # ASSERT: THE FLUX COUNTS STEP AT THE JOIN, THE DENSITY DOES NOT
    wave = _wave(merged)
    assert wave.min() < joins["11-12"] < wave.max()
    np.testing.assert_allclose(_column(merged, "FLUX_DENSITY_COUNTS"), 200.0, rtol=1e-6)
    assert _column(merged, "FLUX_COUNTS").max() / _column(merged, "FLUX_COUNTS").min() > 1.4


def test_the_no_data_bins_stay_reachable_after_the_merge_without_a_new_column(log: object) -> None:
    # ARRANGE: THE SAME GAP AS test_the_gap_between_orders_is_left_empty_and_not_interpolated
    extractor = _extractor(log, "NIR")
    orders = pd.concat([_order(10, 1031.5, 400, 150.0), _order(11, 1000.0, 400, 0.0)], ignore_index=True)

    # ACT
    merged, _ = extractor.merge_extracted_orders(orders)

    # ASSERT: THE MASK IS THE SAME LENGTH AS THE SPECTRUM, TRUE EXACTLY WHERE NO PIXEL OVERLAPS THE BIN
    wave = _wave(merged)
    assert extractor.noDataBins.shape == wave.shape
    inGap = (wave > 1020.1) & (wave < 1031.4)
    assert extractor.noDataBins[inGap].all()
    assert (_column(merged, "VARIANCE")[extractor.noDataBins] == 0.0).all()
    assert (_column(merged, "VARIANCE")[~extractor.noDataBins] > 0.0).all()
    assert list(merged.columns) == ["WAVE", "FLUX_COUNTS", "VARIANCE", "SKY_COUNTS", "SNR"]
