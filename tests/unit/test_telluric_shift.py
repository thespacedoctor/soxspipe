"""Contracts for the telluric-line wavelength shift of a star spectrum (DY-1969)."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.telluric_shift import (
    MIN_BAND_LINES,
    fit_line_centroid,
    load_telluric_lines,
    measure_telluric_shift,
    measure_line_shifts,
    plot_telluric_shift_qc,
    record_telluric_qcs,
    summarise_bands,
    summarise_telluric_shift,
)

pytestmark = pytest.mark.unit

PIXEL_NM = 0.02
LINE_SIGMA_NM = 0.1
HALF_WINDOW_NM = 0.8


def _spectrum(centreNm: float, depth: float, noise: float, seed: int, slope: float = 0.0) -> tuple:
    wavelength = np.arange(1498.0, 1502.0, PIXEL_NM)
    continuum = 1.0 + slope * (wavelength - 1500.0)
    flux = continuum * (1.0 - depth * np.exp(-0.5 * ((wavelength - centreNm) / LINE_SIGMA_NM) ** 2))
    flux = flux + np.random.default_rng(seed).normal(0.0, noise, wavelength.size)
    return wavelength, flux


def test_fit_recovers_the_known_shift_of_a_line_within_its_error() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.05, depth=0.4, noise=0.005, seed=1)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is not None
    assert fit.centroidNm - 1500.0 == pytest.approx(0.05, abs=3 * fit.errorNm)
    assert 0.0 < fit.errorNm < 0.005


def test_fit_recovers_the_shift_on_a_sloped_continuum() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1499.96, depth=0.3, noise=0.004, seed=2, slope=0.4)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is not None
    assert fit.centroidNm - 1500.0 == pytest.approx(-0.04, abs=3 * fit.errorNm)


def test_a_fit_that_runs_into_the_shift_limit_is_rejected() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.35, depth=0.4, noise=0.005, seed=8)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM, maxShiftNm=0.2)

    # ASSERT
    assert fit is None


def test_a_fit_whose_width_runs_into_its_limit_is_rejected() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.0, depth=0.4, noise=0.005, seed=9)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM, sigmaGuessNm=0.02)

    # ASSERT
    assert fit is None


def test_fit_returns_none_when_the_window_holds_too_few_samples() -> None:
    # ARRANGE
    wavelength = np.array([1499.9, 1500.0, 1500.1])

    # ACT
    fit = fit_line_centroid(wavelength, np.ones(3), 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is None


def test_fit_returns_none_when_the_window_has_no_positive_continuum() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.0, depth=0.4, noise=0.005, seed=10)

    # ACT
    fit = fit_line_centroid(wavelength, -flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is None


def test_fit_returns_none_when_the_optimiser_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.0, depth=0.4, noise=0.005, seed=11)

    def failing(*args: object, **kwargs: object) -> None:
        raise RuntimeError("Optimal parameters not found")

    monkeypatch.setattr("scipy.optimize.curve_fit", failing)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is None


def test_fit_returns_none_when_the_centre_has_no_finite_error(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.0, depth=0.4, noise=0.005, seed=13)
    params = np.array([1500.0, 0.1, 0.4, 1.0, 0.0])

    def unconstrained(*args: object, **kwargs: object) -> tuple:
        return params, np.full((5, 5), np.inf)

    monkeypatch.setattr("scipy.optimize.curve_fit", unconstrained)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is None


def test_a_noise_only_window_gives_no_usable_centroid() -> None:
    # ARRANGE
    wavelength, flux = _spectrum(1500.0, depth=0.0, noise=0.05, seed=3)

    # ACT
    fit = fit_line_centroid(wavelength, flux, 1500.0, HALF_WINDOW_NM)

    # ASSERT
    assert fit is None or fit.errorNm > 0.02 or fit.depth < 0.05


def _line_rows(band: str, shiftsNm: list[float], errorNm: float = 0.005) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "band": band,
            "lineNm": 1500.0 + np.arange(len(shiftsNm)),
            "shiftNm": shiftsNm,
            "errorNm": errorNm,
            "depth": 0.3,
            "order": 10,
        }
    )


def test_a_band_with_the_minimum_number_of_good_lines_is_measured() -> None:
    # ARRANGE
    lines = _line_rows("H2O_1130", [0.09, 0.10, 0.11, 0.10, 0.09][:MIN_BAND_LINES])

    # ACT
    bands = summarise_bands(lines)

    # ASSERT
    row = bands.set_index("band").loc["H2O_1130"]
    assert bool(row["isMeasured"]) is True
    assert row["nLines"] == MIN_BAND_LINES
    assert row["shiftNm"] == pytest.approx(0.10)


def test_a_band_one_line_short_of_the_minimum_is_not_measured() -> None:
    # ARRANGE
    lines = _line_rows("H2O_1130", [0.10] * (MIN_BAND_LINES - 1))

    # ACT
    bands = summarise_bands(lines)

    # ASSERT
    row = bands.set_index("band").loc["H2O_1130"]
    assert bool(row["isMeasured"]) is False
    assert row["nLines"] == MIN_BAND_LINES - 1
    assert np.isnan(row["shiftNm"])


def test_lines_with_a_centroid_error_above_the_limit_do_not_count_towards_the_gate() -> None:
    # ARRANGE
    good = _line_rows("CO2_2000", [0.0] * (MIN_BAND_LINES - 1))
    noisy = _line_rows("CO2_2000", [0.5], errorNm=0.05)

    # ACT
    bands = summarise_bands(pd.concat([good, noisy], ignore_index=True))

    # ASSERT
    row = bands.set_index("band").loc["CO2_2000"]
    assert bool(row["isMeasured"]) is False
    assert row["nLines"] == MIN_BAND_LINES - 1


def test_overall_shift_is_the_median_of_the_lines_in_measured_bands_only() -> None:
    # ARRANGE
    measured = _line_rows("H2O_1130", [0.08, 0.09, 0.10, 0.11, 0.12])
    unmeasured = _line_rows("H2O_1450", [0.90, 0.95])

    # ACT
    result = summarise_telluric_shift(pd.concat([measured, unmeasured], ignore_index=True))

    # ASSERT
    assert result.shiftNm == pytest.approx(0.10)
    assert result.nLines == 5
    expectedError = 1.4826 * 0.01 / np.sqrt(5)
    assert result.errorNm == pytest.approx(expectedError)
    assert result.bands.set_index("band")["isMeasured"].to_dict() == {"H2O_1130": True, "H2O_1450": False}


def test_overall_shift_is_missing_when_no_band_is_measured() -> None:
    # ARRANGE
    lines = _line_rows("O2_A", [0.1, 0.1])

    # ACT
    result = summarise_telluric_shift(lines)

    # ASSERT
    assert np.isnan(result.shiftNm)
    assert np.isnan(result.errorNm)
    assert result.nLines == 0


LINE_CENTRES_NM = np.array([1115.0, 1120.0, 1125.0, 1130.0, 1135.0, 1140.0])


def _order_frame(order: int, startNm: float, stopNm: float, shiftNm: float, noise: float, seed: int) -> pd.DataFrame:
    wavelength = np.arange(startNm, stopNm, 0.02)
    flux = np.ones_like(wavelength)
    for centre in LINE_CENTRES_NM:
        flux -= 0.4 * np.exp(-0.5 * ((wavelength - centre - shiftNm) / 0.1) ** 2)
    flux += np.random.default_rng(seed).normal(0.0, noise, wavelength.size)
    return pd.DataFrame(
        {
            "order": order,
            "wavelengthMean": wavelength,
            "extractedFluxOptimal": flux * 1000.0,
            "varianceSpectrum": (noise * 1000.0) ** 2,
        }
    )


def _line_table() -> pd.DataFrame:
    return pd.DataFrame({"BAND": "H2O_1130", "WAVELENGTH": LINE_CENTRES_NM, "STRENGTH": 1e-24, "SPECIES": "H2O"})


def test_measure_line_shifts_recovers_the_shift_of_every_line_in_the_spectrum() -> None:
    # ARRANGE
    spectrum = _order_frame(10, 1110.0, 1145.0, shiftNm=-0.08, noise=0.004, seed=4)

    # ACT
    shifts = measure_line_shifts(spectrum, _line_table())

    # ASSERT
    assert len(shifts) == len(LINE_CENTRES_NM)
    assert shifts["band"].eq("H2O_1130").all()
    assert shifts["shiftNm"].to_numpy() == pytest.approx(-0.08, abs=0.01)
    assert (shifts["errorNm"] < 0.01).all()


def test_a_line_too_close_to_the_edge_of_an_order_for_its_window_is_not_fitted() -> None:
    # ARRANGE
    spectrum = _order_frame(10, 1110.0, 1145.0, shiftNm=0.0, noise=0.004, seed=14)
    lines = _line_table().iloc[[0]].assign(WAVELENGTH=1110.3)

    # ACT
    shifts = measure_line_shifts(spectrum, lines)

    # ASSERT
    assert shifts.empty


def test_measure_telluric_shift_reports_every_band_of_the_line_list() -> None:
    # ARRANGE
    spectrum = _order_frame(10, 1110.0, 1145.0, shiftNm=-0.08, noise=0.004, seed=15)
    lines = pd.concat([_line_table(), _line_table().assign(BAND="H2O_1450", WAVELENGTH=1450.0)], ignore_index=True)

    # ACT
    result = measure_telluric_shift(spectrum, lines)

    # ASSERT
    assert result.bands.set_index("band")["isMeasured"].to_dict() == {"H2O_1130": True, "H2O_1450": False}
    assert result.shiftNm == pytest.approx(-0.08, abs=0.01)


def test_an_order_with_too_few_samples_is_skipped() -> None:
    # ARRANGE
    tiny = _order_frame(10, 1110.0, 1145.0, shiftNm=0.0, noise=0.004, seed=12).iloc[:4]

    # ACT
    shifts = measure_line_shifts(tiny, _line_table())

    # ASSERT
    assert shifts.empty


def test_a_line_in_two_overlapping_orders_is_reported_once_from_the_better_fit() -> None:
    # ARRANGE
    noisy = _order_frame(10, 1110.0, 1145.0, shiftNm=0.30, noise=0.08, seed=5)
    clean = _order_frame(11, 1110.0, 1145.0, shiftNm=0.05, noise=0.003, seed=6)

    # ACT
    shifts = measure_line_shifts(pd.concat([noisy, clean], ignore_index=True), _line_table())

    # ASSERT
    assert len(shifts) == len(LINE_CENTRES_NM)
    assert shifts["order"].eq(11).all()
    assert shifts["shiftNm"].to_numpy() == pytest.approx(0.05, abs=0.01)


def _recorded(result) -> pd.DataFrame:
    return record_telluric_qcs(
        pd.DataFrame(),
        result,
        recipeName="soxs-nod-std",
        dateObs="2026-01-11T02:43:42",
        reductionDateUtc="2026-11-26T10:00:00",
    ).set_index("qc_name")


def test_the_three_header_qcs_carry_the_shift_its_error_and_the_line_count() -> None:
    # ARRANGE
    lines = _line_rows("H2O_1130", [0.08, 0.09, 0.10, 0.11, 0.12])
    result = summarise_telluric_shift(lines)

    # ACT
    qc = _recorded(result)

    # ASSERT
    assert qc.loc["TELL SHIFT", "qc_value"] == pytest.approx(0.10)
    assert qc.loc["TELL SHIFT", "qc_unit"] == "nm"
    assert qc.loc["TELL SHIFT ERR", "qc_value"] == pytest.approx(result.errorNm, abs=1e-5)
    assert qc.loc["TELL NLINES", "qc_value"] == 5
    assert qc.loc[["TELL SHIFT", "TELL SHIFT ERR", "TELL NLINES"], "to_header"].all()


def test_the_per_band_values_go_to_the_qc_table_and_not_to_the_header() -> None:
    # ARRANGE
    lines = pd.concat(
        [_line_rows("H2O_1130", [0.08, 0.09, 0.10, 0.11, 0.12]), _line_rows("H2O_1450", [0.9, 0.95])],
        ignore_index=True,
    )

    # ACT
    qc = _recorded(summarise_telluric_shift(lines))

    # ASSERT
    assert qc.loc["TELL SHIFT H2O_1130", "qc_value"] == pytest.approx(0.10)
    assert qc.loc["TELL NLINES H2O_1130", "qc_value"] == 5
    assert qc.loc["TELL SHIFT H2O_1450", "qc_value"] == "not measured"
    assert qc.loc["TELL SHIFT ERR H2O_1450", "qc_value"] == "not measured"
    assert qc.loc["TELL NLINES H2O_1450", "qc_value"] == 2
    perBand = [name for name in qc.index if name.endswith(("H2O_1130", "H2O_1450"))]
    assert not qc.loc[perBand, "to_header"].any()


def test_a_run_with_no_measured_band_writes_a_zero_line_count_and_no_header_shift() -> None:
    # ARRANGE
    result = summarise_telluric_shift(_line_rows("O2_A", [0.1, 0.1]))

    # ACT
    qc = _recorded(result)

    # ASSERT
    assert qc.loc["TELL NLINES", "qc_value"] == 0
    assert qc.loc["TELL NLINES", "to_header"]
    for name in ("TELL SHIFT", "TELL SHIFT ERR"):
        assert qc.loc[name, "qc_value"] == "not measured"
        assert not qc.loc[name, "to_header"]


NIR_BANDS = {"H2O_1130", "H2O_1450", "H2O_1800", "CO2_2000"}
VIS_BANDS = {"O2_B", "H2O_720", "O2_A", "H2O_940"}


@pytest.mark.parametrize(("arm", "bands"), [("NIR", NIR_BANDS), ("VIS", VIS_BANDS)])
def test_the_shipped_line_list_holds_only_the_bands_of_the_requested_arm(
    log: object, arm: str, bands: set[str]
) -> None:
    # ARRANGE
    settings = {"instrument": "soxs"}

    # ACT
    lines = load_telluric_lines(log, settings, arm)

    # ASSERT
    assert lines is not None
    assert set(lines["BAND"]) <= bands
    assert lines["WAVELENGTH"].is_monotonic_increasing
    assert lines["WAVELENGTH"].dtype == float


def test_an_instrument_without_a_line_list_has_no_telluric_lines(log: object) -> None:
    # ARRANGE
    settings = {"instrument": "xsh"}

    # ACT
    lines = load_telluric_lines(log, settings, "NIR")

    # ASSERT
    assert lines is None


def test_the_qc_plot_is_written_and_leaves_no_figure_open(tmp_path: Path) -> None:
    # ARRANGE
    lines = pd.concat(
        [_line_rows("H2O_1130", [0.08, 0.09, 0.10, 0.11, 0.12]), _line_rows("H2O_1450", [0.9, 0.95])], ignore_index=True
    )
    result = summarise_telluric_shift(lines)
    filePath = tmp_path / "tell.pdf"
    figuresBefore = plt.get_fignums()

    # ACT
    written = plot_telluric_shift_qc(result, str(filePath), "NIR telluric-line shift")

    # ASSERT
    assert written == str(filePath)
    assert filePath.read_bytes().startswith(b"%PDF")
    assert plt.get_fignums() == figuresBefore


def test_the_qc_plot_is_written_when_no_line_was_fitted(tmp_path: Path) -> None:
    # ARRANGE
    result = summarise_telluric_shift(
        measure_line_shifts(_order_frame(10, 1110.0, 1145.0, 0.0, 0.004, 7).iloc[:5], _line_table())
    )
    filePath = tmp_path / "empty.pdf"

    # ACT
    plot_telluric_shift_qc(result, str(filePath), "NIR telluric-line shift")

    # ASSERT
    assert filePath.read_bytes().startswith(b"%PDF")
