"""Contracts for the slit-position drift QC series and plot."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest
from astropy.io.fits import Header

from soxspipe.commonutils.slit_drift_qc import (
    CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES,
    SLIT_DRIFT_EDGE_TRIM_FRACTION,
    SLIT_DRIFT_FIT_SAMPLES,
    format_time_offset,
    pinhole_offset_labels,
    plot_slit_drift_qc,
    rotator_angle_offset,
    slit_drift_arcsec,
    slit_drift_series,
    slit_zero_crossings,
)
from tests.unit._plot_spies import spy_figures

pytestmark = pytest.mark.unit

MEAN_SLIT_ARCSEC = 0.25


def _trace_table() -> pd.DataFrame:
    """Two orders: order 10 (red) drifts linearly with noise, order 11 (blue) is flat."""
    return pd.DataFrame(
        {
            "order": [10, 10, 10, 11, 11],
            "wavelength": [600.0, 601.0, 602.0, 500.0, 501.0],
            "slit_position": [1.1, 1.9, 3.2, 0.5, 0.5],
        }
    )


def _two_order_series(**overrides: object) -> list[dict]:
    arguments = {
        "orderPixelTable": _trace_table(),
        "orders": [10, 11],
        "centreCoeffs": [np.array([1.0, -599.0]), np.array([0.5])],
        "wlMinMax": [(600.0, 602.0), (500.0, 501.0)],
        "fallbackFlags": [False, False],
    }
    arguments.update(overrides)
    return slit_drift_series(**arguments)


def test_residual_is_trace_point_minus_fitted_polynomial() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order10 = next(entry for entry in series if entry["order"] == 10)
    np.testing.assert_allclose(order10["traceSlit"], [1.1, 1.9, 3.2])
    np.testing.assert_allclose(order10["residual"], [0.1, -0.1, 0.2], atol=1e-12)
    assert order10["isFallback"] is False


def test_non_finite_trace_points_are_dropped() -> None:
    # ARRANGE
    table = pd.DataFrame(
        {
            "order": [10, 10, 10, 10],
            "wavelength": [600.0, np.nan, 602.0, 603.0],
            "slit_position": [1.0, 2.0, np.inf, 4.0],
        }
    )

    # ACT
    series = slit_drift_series(table, [10], [np.array([1.0, -599.0])], [(600.0, 603.0)], [False])

    # ASSERT
    np.testing.assert_allclose(series[0]["traceWavelength"], [600.0, 603.0])
    np.testing.assert_allclose(series[0]["traceSlit"], [1.0, 4.0])
    assert len(series[0]["residual"]) == 2


def test_only_points_for_the_order_are_used() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order11 = next(entry for entry in series if entry["order"] == 11)
    np.testing.assert_allclose(order11["traceWavelength"], [500.0, 501.0])


def test_fallback_order_has_no_trace_and_constant_fit_over_its_wavelength_range() -> None:
    # ARRANGE
    table = _trace_table()
    table = table[table["order"] == 10]

    # ACT
    series = slit_drift_series(
        table,
        [10, 11],
        [np.array([1.0, -599.0]), np.array([MEAN_SLIT_ARCSEC])],
        [(600.0, 602.0), (500.0, 501.0)],
        [False, True],
    )

    # ASSERT
    fallback = next(entry for entry in series if entry["order"] == 11)
    assert fallback["isFallback"] is True
    assert len(fallback["traceWavelength"]) == 0
    assert len(fallback["residual"]) == 0
    assert len(fallback["fitWavelength"]) == SLIT_DRIFT_FIT_SAMPLES
    assert fallback["fitWavelength"][0] == pytest.approx(500.0)
    assert fallback["fitWavelength"][-1] == pytest.approx(501.0)
    np.testing.assert_allclose(fallback["fitSlit"], MEAN_SLIT_ARCSEC)


def test_fit_curve_is_the_polynomial_sampled_across_the_order_range() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    order10 = next(entry for entry in series if entry["order"] == 10)
    expected = np.polyval([1.0, -599.0], order10["fitWavelength"])
    np.testing.assert_allclose(order10["fitSlit"], expected)


def test_series_is_sorted_blue_to_red_by_minimum_wavelength() -> None:
    # ARRANGE / ACT
    series = _two_order_series()

    # ASSERT
    assert [entry["order"] for entry in series] == [11, 10]


def test_plot_writes_a_pdf_and_leaves_no_figure_open(tmp_path: Path) -> None:
    # ARRANGE
    series = _two_order_series()
    filePath = tmp_path / "SLIT_DRIFT.pdf"
    plt.close("all")

    # ACT
    result = plot_slit_drift_qc(series, MEAN_SLIT_ARCSEC, "title", str(filePath))

    # ASSERT
    assert result == str(filePath)
    assert filePath.read_bytes().startswith(b"%PDF")
    assert plt.get_fignums() == []


def test_plot_has_two_labelled_panels_and_a_mean_line(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    closed: list[tuple] = []
    realClose = plt.close

    def recording_close(*args: object) -> None:
        closed.append(args)
        realClose(*args)

    monkeypatch.setattr(plt, "close", recording_close)

    # ACT
    plot_slit_drift_qc(_two_order_series(), MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plot.pdf"))

    # ASSERT
    axes = figures[0].axes
    assert len(axes) == 2
    mainAxis, residualAxis = axes
    assert residualAxis.get_xlabel() == "wavelength (nm)"
    assert mainAxis.get_ylabel() == "slit position (arcsec)"
    assert residualAxis.get_ylabel() == "residual (arcsec)"
    dashed = [line for line in mainAxis.get_lines() if line.get_linestyle() == "--"]
    assert len(dashed) == 1
    assert list(dashed[0].get_ydata()) == [MEAN_SLIT_ARCSEC, MEAN_SLIT_ARCSEC]
    assert "former single centre" in dashed[0].get_label()
    assert "#555" not in dashed[0].get_label()
    assert f"{MEAN_SLIT_ARCSEC:.3f}" in dashed[0].get_label()
    assert closed == [(figures[0],)]


def test_plot_marks_fallback_order_with_a_span_and_text(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    series = _two_order_series(fallbackFlags=[False, True], orderPixelTable=_trace_table()[lambda t: t["order"] == 10])

    # ACT
    plot_slit_drift_qc(series, MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plot.pdf"))

    # ASSERT
    mainAxis = figures[0].axes[0]
    texts = [text.get_text() for text in mainAxis.texts]
    assert "order 11: no usable trace fit, fell back to mean" in texts
    assert len(mainAxis.patches) == 1


def test_plot_closes_figure_when_saving_fails(tmp_path: Path) -> None:
    # ARRANGE
    plt.close("all")
    missingDirectory = tmp_path / "missing" / "plot.pdf"

    # ACT / ASSERT
    with pytest.raises(FileNotFoundError):
        plot_slit_drift_qc(_two_order_series(), MEAN_SLIT_ARCSEC, "title", str(missingDirectory))
    assert plt.get_fignums() == []


def _entry(wavelength: list[float], slit: list[float], isFallback: bool = False) -> dict:
    """A minimal series entry carrying only the keys the metric functions read."""
    return {
        "fitWavelength": np.array(wavelength, dtype=float),
        "fitSlit": np.array(slit, dtype=float),
        "isFallback": isFallback,
    }


def test_drift_is_trimmed_reddest_fitted_slit_position_minus_trimmed_bluest() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [0.5, 1.0]), _entry([600.0, 610.0], [2.0, 3.5])]
    trim = SLIT_DRIFT_EDGE_TRIM_FRACTION
    blueSlit = 0.5 + trim * 0.5
    redSlit = 3.5 - trim * 1.5

    # ACT
    drift = slit_drift_arcsec(series)

    # ASSERT
    assert drift == pytest.approx(redSlit - blueSlit)


def test_drift_is_negative_when_trace_moves_towards_lower_slit_positions() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [2.0, 1.5]), _entry([600.0, 610.0], [1.0, 0.5])]
    trim = SLIT_DRIFT_EDGE_TRIM_FRACTION
    expected = (0.5 + trim * 0.5) - (2.0 - trim * 0.5)

    # ACT
    drift = slit_drift_arcsec(series)

    # ASSERT
    assert drift == pytest.approx(expected)
    assert drift < 0


def test_drift_skips_fallback_orders_at_either_end() -> None:
    # ARRANGE
    series = [
        _entry([400.0, 410.0], [9.0, 9.0], isFallback=True),
        _entry([500.0, 510.0], [0.5, 1.0]),
        _entry([600.0, 610.0], [2.0, 3.5]),
        _entry([700.0, 710.0], [-9.0, -9.0], isFallback=True),
    ]
    trim = SLIT_DRIFT_EDGE_TRIM_FRACTION

    # ACT
    drift = slit_drift_arcsec(series)

    # ASSERT
    assert drift == pytest.approx((3.5 - trim * 1.5) - (0.5 + trim * 0.5))


def test_drift_of_a_single_order_is_its_trimmed_end_to_end_change() -> None:
    # ARRANGE
    wavelength = np.linspace(500.0, 600.0, 101)
    series = [_entry(list(wavelength), list((wavelength - 500.0) / 100.0))]

    # ACT
    drift = slit_drift_arcsec(series)

    # ASSERT
    assert drift == pytest.approx(1.0 - 2 * SLIT_DRIFT_EDGE_TRIM_FRACTION)


def test_drift_ignores_a_spike_within_the_trimmed_edge_of_the_range() -> None:
    # ARRANGE
    wavelength = np.linspace(500.0, 600.0, 101)
    slit = (wavelength - 500.0) / 100.0
    spiked = slit.copy()
    spiked[:4] = 50.0
    spiked[-4:] = -50.0

    # ACT
    drift = slit_drift_arcsec([_entry(list(wavelength), list(spiked))])

    # ASSERT
    assert drift == pytest.approx(slit_drift_arcsec([_entry(list(wavelength), list(slit))]))


def test_drift_uses_the_bluest_and_reddest_ends_even_when_one_order_is_nested_inside_another() -> None:
    # ARRANGE
    wavelength = np.linspace(500.0, 600.0, 101)
    wide = _entry(list(wavelength), list((wavelength - 500.0) / 100.0))
    nested = _entry([540.0, 560.0], [9.0, 9.0])

    # ACT
    drift = slit_drift_arcsec([wide, nested])

    # ASSERT
    assert drift == pytest.approx(1.0 - 2 * SLIT_DRIFT_EDGE_TRIM_FRACTION)


def test_drift_is_finite_when_an_end_wavelength_is_nan() -> None:
    # ARRANGE
    wavelength = np.linspace(500.0, 600.0, 101)
    wavelength[0] = np.nan
    series = [_entry(list(wavelength), list((np.linspace(500.0, 600.0, 101) - 500.0) / 100.0))]

    # ACT
    drift = slit_drift_arcsec(series)

    # ASSERT
    assert drift == pytest.approx(0.891)


def test_drift_is_none_when_the_fitted_values_are_not_finite() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [np.nan, np.nan])]

    # ACT / ASSERT
    assert slit_drift_arcsec(series) is None


def test_drift_is_none_when_every_order_fell_back() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [0.25, 0.25], isFallback=True)]

    # ACT / ASSERT
    assert slit_drift_arcsec(series) is None


def test_drift_is_none_for_an_empty_series() -> None:
    # ACT / ASSERT
    assert slit_drift_arcsec([]) is None


def test_zero_crossing_inside_one_order_is_linearly_interpolated() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0, 520.0], [-1.0, -0.5, 1.5])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([512.5])


def test_a_sign_change_between_orders_across_a_gap_gives_no_crossing() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [-2.0, -1.0]), _entry([530.0, 540.0], [3.0, 4.0])]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_overlapping_orders_with_opposite_sign_edges_do_not_create_false_crossings() -> None:
    # ARRANGE
    series = [_entry([500.0, 520.0], [-2.0, -1.0]), _entry([510.0, 530.0], [1.0, 2.0])]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_a_fallback_order_between_two_valid_orders_is_not_bridged() -> None:
    # ARRANGE
    series = [
        _entry([500.0, 510.0], [-1.0, -1.0]),
        _entry([520.0, 530.0], [0.25, 0.25], isFallback=True),
        _entry([540.0, 550.0], [1.0, 1.0]),
    ]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_a_run_of_exact_zero_samples_counts_as_one_crossing_at_its_bluest_sample() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0, 520.0, 530.0], [-1.0, 0.0, 0.0, 1.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([510.0])


def test_no_zero_crossing_returns_an_empty_list() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [0.5, 1.0]), _entry([600.0, 610.0], [2.0, 3.5])]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_two_zero_crossings_are_returned_blue_first() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0], [-1.0, 1.0]), _entry([600.0, 610.0], [1.0, -1.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([505.0, 605.0])


def test_crossing_exactly_on_a_sample_counts_once() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0, 520.0], [-1.0, 0.0, 1.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([510.0])


def test_zero_crossings_ignore_fallback_orders() -> None:
    # ARRANGE
    series = [
        _entry([500.0, 510.0], [-1.0, 1.0], isFallback=True),
        _entry([600.0, 610.0], [1.0, 2.0]),
    ]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_an_exact_zero_sample_at_a_non_finite_wavelength_is_not_a_crossing() -> None:
    # ARRANGE
    series = [_entry([500.0, np.nan, 520.0, 530.0], [-1.0, 0.0, -1.0, 1.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([525.0])


def test_a_sign_change_across_a_non_finite_wavelength_is_not_a_crossing() -> None:
    # ARRANGE
    series = [_entry([500.0, np.nan, 520.0], [-1.0, 1.0, 2.0])]

    # ACT / ASSERT
    assert slit_zero_crossings(series) == []


def test_a_zero_at_a_non_finite_wavelength_does_not_hide_the_next_exact_zero() -> None:
    # ARRANGE
    series = [_entry([np.nan, 510.0], [0.0, 0.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == pytest.approx([510.0])


def test_opposite_signed_infinite_slit_values_give_no_crossing() -> None:
    # ARRANGE
    series = [_entry([500.0, 510.0, 520.0], [-np.inf, np.inf, 1.0])]

    # ACT
    crossings = slit_zero_crossings(series)

    # ASSERT
    assert crossings == []


def test_zero_crossings_of_an_empty_series_is_empty() -> None:
    # ACT / ASSERT
    assert slit_zero_crossings([]) == []


def test_plot_with_drift_and_zero_crossing_writes_pdf_and_annotates_both_axes(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    filePath = tmp_path / "annotated.pdf"

    # ACT
    plot_slit_drift_qc(
        _two_order_series(),
        MEAN_SLIT_ARCSEC,
        "title",
        str(filePath),
        driftArcsec=2.25,
        zeroWavelength=550.0,
    )

    # ASSERT
    assert filePath.read_bytes().startswith(b"%PDF")
    mainAxis, residualAxis = figures[0].axes
    legendLabels = [text.get_text() for text in mainAxis.get_legend().get_texts()]
    assert "slit drift (red - blue) = +2.250 arcsec" in legendLabels
    assert "trace crosses slit 0 at 550.0 nm" in [text.get_text() for text in mainAxis.texts]
    for axis in (mainAxis, residualAxis):
        verticals = [line for line in axis.get_lines() if list(line.get_xdata()) == [550.0, 550.0]]
        assert len(verticals) == 1
        assert verticals[0].get_linestyle() == ":"
    horizontals = [line for line in mainAxis.get_lines() if list(line.get_ydata()) == [0, 0]]
    assert len(horizontals) == 1


def test_plot_without_new_keywords_has_no_drift_legend_or_zero_crossing_label(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)

    # ACT
    plot_slit_drift_qc(_two_order_series(), MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plain.pdf"))

    # ASSERT
    mainAxis = figures[0].axes[0]
    legendLabels = [text.get_text() for text in mainAxis.get_legend().get_texts()]
    assert not any("slit drift" in label for label in legendLabels)
    assert list(mainAxis.texts) == []
    assert [line for line in mainAxis.get_lines() if list(line.get_ydata()) == [0, 0]] == []


@pytest.mark.parametrize(
    ("seconds", "expected"),
    [
        (42 * 60.0, "+42.0 min"),
        (59 * 60.0, "+59.0 min"),
        (61 * 60.0, "+1.0 h"),
        (23 * 3600.0, "+23.0 h"),
        (25 * 3600.0, "+1.0 d"),
        (-42 * 60.0, "-42.0 min"),
        (3599.99, "+1.0 h"),
        (-3599.99, "-1.0 h"),
        (3600.0, "+1.0 h"),
        (-3600.0, "-1.0 h"),
        (86399.0, "+1.0 d"),
        (86400.0, "+1.0 d"),
        (-86400.0, "-1.0 d"),
        (-0.4, "+0.0 min"),
        (-5.3 * 3600.0, "-5.3 h"),
        (2.1 * 86400.0, "+2.1 d"),
        (0.0, "+0.0 min"),
    ],
)
def test_time_offset_uses_the_unit_that_fits_its_size_and_keeps_its_sign(seconds: float, expected: str) -> None:
    # ARRANGE / ACT
    text = format_time_offset(seconds)

    # ASSERT
    assert text == expected


@pytest.mark.parametrize(
    ("scienceAngle", "pinholeAngle", "expected"),
    [
        (179.0, -179.0, -2.0),
        (-179.0, 179.0, 2.0),
        (10.0, 5.0, 5.0),
        (5.0, 10.0, -5.0),
        (180.0, 0.0, 180.0),
        (0.0, 180.0, 180.0),
    ],
)
def test_rotator_offset_is_wrapped_into_minus_180_to_180(
    scienceAngle: float, pinholeAngle: float, expected: float
) -> None:
    # ARRANGE / ACT
    offset = rotator_angle_offset(scienceAngle, pinholeAngle)

    # ASSERT
    assert offset == pytest.approx(expected)


KW = {"MJDOBS": "MJD-OBS", "DATE_OBS": "DATE-OBS", "ABSROT": "ESO ADA ABSROT END"}.get


def _offset_header(**cards: object) -> Header:
    return Header(cards)


def test_labels_give_time_and_rotator_offsets_when_both_headers_carry_them() -> None:
    # ARRANGE
    science = _offset_header(**{"MJD-OBS": 60000.5, "ESO ADA ABSROT END": 10.0})
    pinhole = _offset_header(**{"MJD-OBS": 60000.5 - 42.0 / 1440.0, "ESO ADA ABSROT END": 13.2})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == [
        "science - multi-pinhole obs time = +42.0 min",
        "science - multi-pinhole rotator angle (ABSROT) = -3.2 deg",
    ]


def test_labels_fall_back_to_the_obs_date_when_mjd_is_missing() -> None:
    # ARRANGE
    science = _offset_header(**{"DATE-OBS": "2024-01-02T05:00:00"})
    pinhole = _offset_header(**{"DATE-OBS": "2024-01-02T03:00:00"})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole obs time = +2.0 h"]


def test_labels_drop_only_the_time_entry_when_one_side_has_no_usable_time() -> None:
    # ARRANGE
    science = _offset_header(**{"DATE-OBS": "not a date", "ESO ADA ABSROT END": 1.0})
    pinhole = _offset_header(**{"MJD-OBS": 60000.0, "ESO ADA ABSROT END": 0.0})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole rotator angle (ABSROT) = +1.0 deg"]


def test_labels_drop_only_the_rotator_entry_when_the_science_absrot_is_missing() -> None:
    # ARRANGE
    science = _offset_header(**{"MJD-OBS": 60001.0})
    pinhole = _offset_header(**{"MJD-OBS": 60000.0, "ESO ADA ABSROT END": 4.0})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole obs time = +1.0 d"]


def test_labels_assume_the_parked_rotator_angle_when_the_pinhole_header_has_no_absrot() -> None:
    # ARRANGE
    science = _offset_header(**{"MJD-OBS": 60001.0, "ESO ADA ABSROT END": 84.9})
    pinhole = _offset_header(**{"MJD-OBS": 60000.0})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == [
        "science - multi-pinhole obs time = +1.0 d",
        "science - multi-pinhole rotator angle (ABSROT) = +84.9 deg "
        f"(multi-pinhole ABSROT assumed {CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES:.1f} deg)",
    ]


def test_labels_use_a_real_pinhole_absrot_over_the_parked_angle() -> None:
    # ARRANGE
    science = _offset_header(**{"ESO ADA ABSROT END": 84.9})
    pinhole = _offset_header(**{"ESO ADA ABSROT END": 80.0})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole rotator angle (ABSROT) = +4.9 deg"]


def test_labels_treat_a_boolean_pinhole_absrot_as_missing() -> None:
    # ARRANGE
    science = _offset_header(**{"ESO ADA ABSROT END": 10.0})
    pinhole = _offset_header(**{"ESO ADA ABSROT END": True})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == [
        "science - multi-pinhole rotator angle (ABSROT) = +10.0 deg "
        f"(multi-pinhole ABSROT assumed {CALIBRATION_ROTATOR_PARK_ANGLE_DEGREES:.1f} deg)"
    ]


def test_labels_fall_back_to_the_obs_date_when_mjd_is_a_boolean() -> None:
    # ARRANGE
    science = _offset_header(**{"MJD-OBS": True, "DATE-OBS": "2024-01-02T05:00:00"})
    pinhole = _offset_header(**{"MJD-OBS": np.bool_(True), "DATE-OBS": "2024-01-02T03:00:00"})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole obs time = +2.0 h"]


def test_labels_skip_header_values_too_large_for_a_float() -> None:
    # ARRANGE
    science = {"MJD-OBS": 10**400, "DATE-OBS": "2024-01-02T05:00:00", "ESO ADA ABSROT END": 10**400}
    pinhole = {"MJD-OBS": 10**400, "DATE-OBS": "2024-01-02T03:00:00", "ESO ADA ABSROT END": 10.0}

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == ["science - multi-pinhole obs time = +2.0 h"]


def test_labels_drop_the_time_entry_when_the_obs_date_has_no_time_of_day() -> None:
    # ARRANGE
    science = _offset_header(**{"DATE-OBS": "2024-01-02"})
    pinhole = _offset_header(**{"DATE-OBS": "2024-01-01T03:00:00"})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == []


def test_labels_drop_offsets_that_are_not_finite() -> None:
    # ARRANGE
    science = _offset_header(**{"MJD-OBS": 1e308, "ESO ADA ABSROT END": 1e308})
    pinhole = _offset_header(**{"MJD-OBS": -1e308, "ESO ADA ABSROT END": -1e308})

    # ACT
    labels = pinhole_offset_labels(science, pinhole, KW)

    # ASSERT
    assert labels == []


def test_labels_are_empty_when_nothing_usable_is_in_either_header() -> None:
    # ARRANGE / ACT
    labels = pinhole_offset_labels(_offset_header(), _offset_header(), KW)

    # ASSERT
    assert labels == []


def test_labels_are_empty_when_the_headers_are_not_headers() -> None:
    # ARRANGE / ACT
    labels = pinhole_offset_labels(None, None, KW)

    # ASSERT
    assert labels == []


def test_plot_legend_lists_each_extra_label(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    figures = spy_figures(monkeypatch)
    extraLabels = ["science - multi-pinhole obs time = +42.0 min", "science - multi-pinhole rotator angle = -3.2 deg"]

    # ACT
    plot_slit_drift_qc(
        _two_order_series(), MEAN_SLIT_ARCSEC, "title", str(tmp_path / "plot.pdf"), extraLabels=extraLabels
    )

    # ASSERT
    legendTexts = [text.get_text() for text in figures[0].axes[0].get_legend().get_texts()]
    for label in extraLabels:
        assert label in legendTexts
