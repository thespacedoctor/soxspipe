"""Contracts for the slit-illumination QC plot and the rotator spread helper (DY-1280)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_illumination import (
    FLAT_EDGE_COLUMNS,
    OFFSET_COLUMNS,
    PROFILE_COLUMNS,
    STATUS_BORROWED,
    STATUS_FAILED,
    STATUS_FITTED,
    slit_illumination_result,
    top_hat,
)
from soxspipe.commonutils.slit_illumination_qc import (
    plot_slit_illumination_qc,
    read_rotator_angle,
    rotator_angle_spread,
)

pytestmark = pytest.mark.unit

SLIT = np.round(np.arange(-60, 61) * 0.1, 1)


def _profile_rows(order, offset):
    profile = top_hat(SLIT, 1.0, -5.0 + offset, 5.0 + offset, 0.15)
    return pd.DataFrame(
        {
            "order": order,
            "slit": SLIT,
            "profile": profile,
            "profileErr": np.full(SLIT.size, 0.05),
            "nLines": 10,
        },
        columns=PROFILE_COLUMNS,
    )


def _offset_row(order, offset, status, reason=""):
    return {
        "order": order,
        "offset": offset,
        "offsetErr": 0.02,
        "scienceWidth": 0.15,
        "nBins": 100,
        "status": status,
        "reason": reason,
        "borrowedFrom": np.nan,
    }


def _result(offsetRows, edgeRows, profiles):
    return slit_illumination_result(
        offsets=pd.DataFrame(offsetRows, columns=OFFSET_COLUMNS),
        flatEdges=pd.DataFrame(edgeRows, columns=FLAT_EDGE_COLUMNS),
        profiles=profiles,
        isApplied=True,
        reason="",
        overlapMarginArcsec=0.3,
    )


def _edge_row(order, lower=-5.0, upper=5.0):
    return {"order": order, "lower": lower, "upper": upper, "width": 0.15, "reason": ""}


def test_writes_a_non_empty_pdf_for_a_mix_of_fitted_and_failed_orders(tmp_path):
    result = _result(
        [
            _offset_row(11, 0.8, STATUS_FITTED),
            _offset_row(12, np.nan, STATUS_FAILED, "too few skyline bins"),
            _offset_row(13, -0.4, STATUS_BORROWED),
        ],
        [_edge_row(11), _edge_row(12), _edge_row(13)],
        pd.concat([_profile_rows(11, 0.8), _profile_rows(13, -0.4)], ignore_index=True),
    )
    filePath = tmp_path / "qc.pdf"

    returned = plot_slit_illumination_qc(result, "title", str(filePath))

    assert returned == str(filePath)
    assert filePath.read_bytes().startswith(b"%PDF")
    assert filePath.stat().st_size > 1000


def test_draws_one_panel_per_order_present_in_the_profiles_or_the_offsets(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    figures = []
    realSubplots = plt.subplots

    def recording_subplots(*args, **kwargs):
        fig, axes = realSubplots(*args, **kwargs)
        figures.append((fig, np.atleast_1d(axes)))
        return fig, axes

    monkeypatch.setattr(plt, "subplots", recording_subplots)
    result = _result(
        [_offset_row(11, 0.8, STATUS_FITTED), _offset_row(12, np.nan, STATUS_FAILED, "no lines")],
        [_edge_row(11), _edge_row(12)],
        _profile_rows(11, 0.8),
    )

    plot_slit_illumination_qc(result, "title", str(tmp_path / "qc.pdf"))

    [(_, axes)] = figures
    assert axes.size == 2


def test_handles_an_order_with_no_flat_edges(tmp_path):
    result = _result(
        [_offset_row(11, 0.8, STATUS_FITTED), _offset_row(12, 0.5, STATUS_FITTED)],
        [_edge_row(11)],
        pd.concat([_profile_rows(11, 0.8), _profile_rows(12, 0.5)], ignore_index=True),
    )
    filePath = tmp_path / "qc.pdf"

    plot_slit_illumination_qc(result, "title", str(filePath))

    assert filePath.stat().st_size > 1000


def test_handles_an_order_with_an_offset_but_no_profile_rows(tmp_path):
    result = _result(
        [_offset_row(11, 0.8, STATUS_BORROWED)],
        [_edge_row(11)],
        pd.DataFrame(columns=PROFILE_COLUMNS),
    )
    filePath = tmp_path / "qc.pdf"

    plot_slit_illumination_qc(result, "title", str(filePath))

    assert filePath.stat().st_size > 1000


def test_handles_an_order_with_profile_rows_but_no_offset_row_and_no_flat_edges(tmp_path):
    result = _result([], [], _profile_rows(11, 0.8))
    filePath = tmp_path / "qc.pdf"

    plot_slit_illumination_qc(result, "title", str(filePath))

    assert filePath.stat().st_size > 1000


def _panel_texts(tmp_path, monkeypatch, result):
    """The title and legend labels of every panel drawn for ``result``."""
    import matplotlib.pyplot as plt

    figures = []
    realSubplots = plt.subplots

    def recording_subplots(*args, **kwargs):
        fig, axes = realSubplots(*args, **kwargs)
        figures.append(axes)
        return fig, axes

    monkeypatch.setattr(plt, "subplots", recording_subplots)
    plot_slit_illumination_qc(result, "title", str(tmp_path / "qc.pdf"))
    [axes] = figures
    return [(axis.get_title(), [text.get_text() for text in axis.get_legend().get_texts()]) for axis in axes.ravel()]


def test_each_panel_names_the_method_that_measured_its_order(tmp_path, monkeypatch):
    skylineRow = {**_offset_row(11, 0.8, STATUS_FITTED), "method": "skyline"}
    continuumRow = {**_offset_row(12, 0.5, STATUS_FITTED), "method": "continuum"}
    result = _result(
        [skylineRow, continuumRow],
        [_edge_row(11), _edge_row(12)],
        pd.concat([_profile_rows(11, 0.8), _profile_rows(12, 0.5)], ignore_index=True),
    )

    [(skylineTitle, skylineLabels), (continuumTitle, continuumLabels)] = _panel_texts(tmp_path, monkeypatch, result)

    assert skylineTitle == "order 11: fitted (skyline)"
    assert "skyline profile" in skylineLabels
    assert continuumTitle == "order 12: fitted (continuum)"
    assert "continuum profile" in continuumLabels


def test_a_panel_without_a_recorded_method_keeps_its_plain_title(tmp_path, monkeypatch):
    result = _result([_offset_row(11, 0.8, STATUS_FITTED)], [_edge_row(11)], _profile_rows(11, 0.8))

    [(title, labels)] = _panel_texts(tmp_path, monkeypatch, result)

    assert title == "order 11: fitted"
    assert "skyline profile" in labels


def test_hides_the_unused_panels_of_the_last_row(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    figures = []
    realSubplots = plt.subplots

    def recording_subplots(*args, **kwargs):
        fig, axes = realSubplots(*args, **kwargs)
        figures.append(axes)
        return fig, axes

    monkeypatch.setattr(plt, "subplots", recording_subplots)
    orders = (11, 12, 13, 14)
    result = _result(
        [_offset_row(order, 0.1, STATUS_FITTED) for order in orders],
        [_edge_row(order) for order in orders],
        pd.concat([_profile_rows(order, 0.1) for order in orders], ignore_index=True),
    )

    plot_slit_illumination_qc(result, "title", str(tmp_path / "qc.pdf"))

    [axes] = figures
    assert [axis.get_visible() for axis in axes.ravel()] == [True] * 4 + [False] * 2


def test_closes_the_figure_when_drawing_fails(tmp_path, monkeypatch):
    import matplotlib.pyplot as plt

    result = _result([_offset_row(11, 0.8, STATUS_FITTED)], [_edge_row(11)], _profile_rows(11, 0.8))
    before = set(plt.get_fignums())

    with pytest.raises(OSError):
        plot_slit_illumination_qc(result, "title", str(tmp_path / "missing-dir" / "qc.pdf"))

    assert set(plt.get_fignums()) == before


@pytest.mark.parametrize(
    ("angles", "expected"),
    [
        ([10.0, 10.0], 0.0),
        ([10.0, 25.0, 17.0], 15.0),
        ([179.0, -179.0], 2.0),
        ([-170.0, 170.0, 180.0], 20.0),
        ([90.0], 0.0),
        ([], 0.0),
        ([10.0, float("nan"), 30.0], 20.0),
    ],
    ids=["same", "ordinary", "wrap", "wrap-three", "single", "none", "nan-ignored"],
)
def test_rotator_angle_spread_is_the_smallest_arc_holding_every_angle(angles, expected):
    assert rotator_angle_spread(angles) == pytest.approx(expected)


@pytest.mark.parametrize(
    ("header", "expected"),
    [
        ({"ABSROT": 12.5}, 12.5),
        ({"ABSROT": "-3.5"}, -3.5),
        ({}, None),
        ({"ABSROT": "junk"}, None),
        ({"ABSROT": float("nan")}, None),
        ({"ABSROT": float("inf")}, None),
        ({"ABSROT": True}, None),
        ({"ABSROT": np.bool_(False)}, None),
        ({"ABSROT": 10**400}, None),
    ],
    ids=["number", "numeric-string", "missing", "unparseable", "nan", "inf", "bool", "numpy-bool", "overflow"],
)
def test_read_rotator_angle_gives_a_finite_float_or_none(header, expected):
    assert read_rotator_angle(header, "ABSROT") == expected
