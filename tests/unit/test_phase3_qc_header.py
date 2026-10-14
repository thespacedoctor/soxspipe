"""Unit tests for `phase3.add_qc_to_header`, the one writer of QC rows into FITS headers."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.phase3 import add_qc_to_header

pytestmark = pytest.mark.unit


def _qc(rows: list[tuple[str, object, bool]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"qc_name": n, "qc_value": v, "qc_comment": f"comment {n}", "to_header": h} for n, v, h in rows],
    ).astype({"qc_value": object})


def test_add_qc_to_header_writes_flagged_rows_with_their_comments() -> None:
    # ARRANGE
    header = fits.Header()
    qc = _qc([("ron", 3.2, True), ("lines", np.int64(42), True), ("hidden", 1.0, False)])

    # ACT
    add_qc_to_header(header, qc)

    # ASSERT
    assert header["ESO QC RON"] == pytest.approx(3.2)
    assert header.comments["ESO QC RON"] == "comment ron"
    assert header["ESO QC LINES"] == 42
    assert "ESO QC HIDDEN" not in header


@pytest.mark.parametrize("value", [float("nan"), np.float64("nan"), float("inf"), np.float32("-inf")])
def test_add_qc_to_header_skips_non_finite_values_that_fits_headers_cannot_hold(value: object) -> None:
    # ARRANGE
    header = fits.Header()
    qc = _qc([("ron", 3.2, True), ("inner order pix mean", value, True)])

    # ACT
    add_qc_to_header(header, qc)

    # ASSERT
    assert header["ESO QC RON"] == pytest.approx(3.2)
    assert "ESO QC INNER ORDER PIX MEAN" not in header
