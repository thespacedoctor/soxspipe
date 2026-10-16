"""The QC plot of slit-alignment shifts (DY-1612)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.slit_alignment import SHIFT_COLUMNS
from soxspipe.commonutils.slit_alignment_qc import plot_slit_alignment_qc
from tests.unit._plot_spies import probe_on_savefig

pytestmark = pytest.mark.unit


def _shift_table(rows: list[tuple[str, float, int, float, bool]]) -> pd.DataFrame:
    """Return a shift table from ``(frame, mjd, order, shift, is_fallback)`` rows."""
    return pd.DataFrame(
        [
            {
                "frame": frame,
                "mjd": mjd,
                "order": order,
                "position_arcsec": np.nan,
                "shift_arcsec": shift,
                "is_fallback": isFallback,
                "is_reference": mjd == min(row[1] for row in rows),
            }
            for frame, mjd, order, shift, isFallback in rows
        ],
        columns=SHIFT_COLUMNS,
    )


def test_the_plot_is_a_pdf_with_one_line_per_order_and_a_marker_for_each_fallback(tmp_path, monkeypatch) -> None:
    # ARRANGE: TWO ORDERS IN THREE FRAMES. THE SECOND FRAME'S ORDER 2 SHIFT WAS NOT MEASURED
    shifts = _shift_table(
        [
            ("a.fits", 60000.0, 1, 0.0, False),
            ("a.fits", 60000.0, 2, 0.0, False),
            ("b.fits", 60000.1, 1, 0.4, False),
            ("b.fits", 60000.1, 2, 0.4, True),
            ("c.fits", 60000.2, 1, 0.5, False),
            ("c.fits", 60000.2, 2, 0.1, False),
        ]
    )
    probes = probe_on_savefig(
        monkeypatch,
        lambda figure: [(line.get_label(), list(line.get_ydata()), line.get_marker()) for line in figure.axes[0].lines],
    )
    filePath = tmp_path / "shifts.pdf"

    # ACT
    plot_slit_alignment_qc(shifts, "Slit alignment", str(filePath))

    # ASSERT
    assert filePath.read_bytes().startswith(b"%PDF")
    lines = probes[0]
    assert [label for label, _, _ in lines if label.startswith("order")] == ["order 1", "order 2"]
    fallbackLine = next(line for line in lines if line[0] == "fallback")
    assert fallbackLine[1] == [0.4]


def test_a_whole_frame_shift_is_drawn_as_one_line_with_the_frame_names_on_the_axis(tmp_path, monkeypatch) -> None:
    # ARRANGE
    shifts = _shift_table(
        [
            ("a.fits", 60000.0, -1, 0.0, False),
            ("b.fits", 60000.1, -1, 0.4, False),
        ]
    )
    probes = probe_on_savefig(
        monkeypatch,
        lambda figure: (
            [line.get_label() for line in figure.axes[0].lines],
            [label.get_text() for label in figure.axes[0].get_xticklabels()],
        ),
    )

    # ACT
    plot_slit_alignment_qc(shifts, "Slit alignment", str(tmp_path / "shifts.pdf"))

    # ASSERT
    labels, tickLabels = probes[0]
    assert [label for label in labels if not label.startswith("_")] == ["frame"]
    assert tickLabels == ["a.fits", "b.fits"]
