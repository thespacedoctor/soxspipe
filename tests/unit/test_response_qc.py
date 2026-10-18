"""Behaviour contracts for the response-fit QC measures."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.response_fit import fit_response_curve
from soxspipe.commonutils.response_qc import measure_response_qcs, response_qc_rows, response_qc_warnings
from tests.factories import qc_table

pytestmark = pytest.mark.unit

FIT_SETTINGS = {
    "bin_width_nm": 2,
    "knot_spacing_nm": 40,
    "clip_sigma_low": 4,
    "clip_sigma_high": 3,
    "max_iteration": 10,
    "min_bin_samples": 5,
    "telluric_regions": [[550.0, 700.0]],
    "stellar_lines": [],
    "line_half_width_nm": 2.5,
    "selfcal_bands": [[360, 380], [450, 500], [780, 840], [2000, 2100]],
    "clip_frac_max": 0.15,
    "max_gap_nm_max": 40,
    "selfcal_res_max": 0.05,
}


def _truth(wave: np.ndarray) -> np.ndarray:
    return np.exp(1.0 + 1.5 * np.exp(-(wave - 350.0) / 200.0))


def _fit(noise: float = 0.01, outlierBlock: bool = False):
    wave = np.linspace(350.0, 850.0, 10000)
    ratio = _truth(wave) * (1.0 + noise * np.random.default_rng(3).standard_normal(wave.size))
    if outlierBlock:
        ratio[(wave > 400.0) & (wave < 420.0)] *= 1.6
        ratio[(wave > 440.0) & (wave < 460.0)] *= 0.5
    return fit_response_curve(wave, ratio, FIT_SETTINGS)


def test_the_fit_object_models_the_truth_and_flags_the_telluric_gap() -> None:
    # ARRANGE AND ACT
    fit = _fit()

    # ASSERT
    probe = np.array([360.0, 500.0, 800.0])
    np.testing.assert_allclose(fit.evaluate(probe), _truth(probe), rtol=0.01)
    assert not fit.mask[(fit.wave > 560.0) & (fit.wave < 690.0)].any()
    assert fit.keptMask.all()


def test_the_widest_gap_is_the_masked_telluric_stretch() -> None:
    # ARRANGE
    fit = _fit()

    # ACT
    measures = measure_response_qcs(fit, FIT_SETTINGS["selfcal_bands"])

    # ASSERT
    assert measures["maxGap"] == pytest.approx(150.0, abs=5.0)


def test_clipped_bins_raise_the_clip_fraction() -> None:
    # ARRANGE
    fit = _fit(outlierBlock=True)

    # ACT
    measures = measure_response_qcs(fit, FIT_SETTINGS["selfcal_bands"])

    # ASSERT
    expected = 1.0 - fit.keptMask.sum() / fit.keptMask.size
    assert measures["clipFrac"] == pytest.approx(expected)
    assert measures["clipFrac"] > 0.0


def test_selfcal_residual_is_the_median_calibrated_over_catalogue_flux_minus_one_and_skips_empty_bands() -> None:
    # ARRANGE
    fit = _fit(noise=0.01)
    inBand = fit.mask & (fit.wave > 780.0) & (fit.wave < 840.0)

    # ACT
    measures = measure_response_qcs(fit, FIT_SETTINGS["selfcal_bands"])

    # ASSERT
    bands = [band for band, _ in measures["selfCal"]]
    assert [780, 840] in bands
    assert [2000, 2100] not in bands
    value = {tuple(band): value for band, value in measures["selfCal"]}[(780, 840)]
    expected = np.median(fit.evaluate(fit.wave[inBand]) / fit.ratio[inBand]) - 1.0
    assert value == pytest.approx(expected)


def test_selfcal_residual_is_unbiased_when_low_signal_counts_go_negative() -> None:
    # ARRANGE: SIGNAL-TO-NOISE OF 1.5 PER SAMPLE IN THE RED, SO MANY EXTRACTED COUNTS ARE NEGATIVE
    wave = np.linspace(350.0, 850.0, 20000)
    noise = np.where(wave > 780.0, 1.0 / 1.5, 0.01)
    counts = (1.0 + noise * np.random.default_rng(5).standard_normal(wave.size)) / _truth(wave)
    fit = fit_response_curve(wave, 1.0 / counts, FIT_SETTINGS)

    # ACT
    measures = measure_response_qcs(fit, [[780, 840]])

    # ASSERT: THE RESIDUAL REFLECTS THE MODEL ERROR (SMALL), NOT THE NOISE
    model = fit.evaluate(np.array([810.0]))[0]
    modelError = model / _truth(np.array([810.0]))[0] - 1.0
    assert measures["selfCal"][0][1] == pytest.approx(modelError, abs=0.05)


def test_qc_rows_follow_the_production_schema() -> None:
    # ARRANGE
    measures = {"clipFrac": 0.04, "maxGap": 152.0, "selfCal": [([360, 380], 0.012)]}

    # ACT
    rows = response_qc_rows(qc_table().iloc[0:0], measures, "soxs-stare", "2024-01-02T03:04:05")

    # ASSERT
    assert rows["qc_name"].tolist() == ["RESP BIN CLIP FRAC", "RESP MAX GAP", "RESP SELFCAL RES 360-380"]
    assert rows["qc_value"].tolist() == [0.04, 152.0, 0.012]
    assert rows["to_header"].all()
    assert set(qc_table().columns) <= set(rows.columns)
    assert rows["soxspipe_recipe"].unique().tolist() == ["soxs-stare"]
    assert isinstance(rows, pd.DataFrame)


def test_out_of_tolerance_measures_produce_one_warning_each() -> None:
    # ARRANGE
    measures = {"clipFrac": 0.4, "maxGap": 250.0, "selfCal": [([360, 380], -0.2), ([450, 500], 0.01)]}

    # ACT
    warnings = response_qc_warnings(measures, FIT_SETTINGS)

    # ASSERT
    assert len(warnings) == 3
    assert any("clipped" in message for message in warnings)
    assert any("gap" in message for message in warnings)
    assert any("360-380" in message for message in warnings)


def test_in_tolerance_measures_produce_no_warnings() -> None:
    # ARRANGE
    measures = {"clipFrac": 0.01, "maxGap": 30.0, "selfCal": [([360, 380], 0.01)]}

    # ACT AND ASSERT
    assert response_qc_warnings(measures, FIT_SETTINGS) == []


def test_qc_rows_into_an_empty_schema_table_keep_the_schema_without_a_pandas_warning() -> None:
    # ARRANGE
    import warnings

    from soxspipe.commonutils.response_function import QC_COLUMNS

    emptyTable = pd.DataFrame(columns=QC_COLUMNS)
    measures = {"clipFrac": 0.04, "maxGap": 152.0, "selfCal": [([360, 380], 0.012)]}

    # ACT
    with warnings.catch_warnings():
        warnings.simplefilter("error", FutureWarning)
        rows = response_qc_rows(emptyTable, measures, "soxs-nod", "2024-01-02T03:04:05")

    # ASSERT
    assert list(rows.columns) == list(QC_COLUMNS)
    assert len(rows) == 3
