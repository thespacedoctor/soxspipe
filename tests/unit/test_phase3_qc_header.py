"""Unit tests for `phase3.add_qc_to_header`, the one writer of QC rows into FITS headers."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.horne_extraction import EXTRACT_METHOD_QC_COMMENT
from soxspipe.commonutils.phase3 import add_qc_to_header
from soxspipe.commonutils.toolkit import append_qc

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


def test_add_qc_to_header_writes_the_extract_method_text_qc_and_it_survives_a_file_round_trip(tmp_path: Path) -> None:
    # ARRANGE
    qc = append_qc(
        pd.DataFrame(),
        recipeName="soxs-stare",
        qcName="EXTRACT METHOD",
        qcValue="boxcar",
        qcComment=EXTRACT_METHOD_QC_COMMENT,
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-03T00:00:00",
        qcUnit=None,
        toHeader=True,
    )
    header = fits.Header()
    filePath = tmp_path / "extract_method.fits"

    # ACT
    with warnings.catch_warnings():
        # A TRUNCATED COMMENT WOULD RAISE
        warnings.filterwarnings("error", message="Card is too long")
        add_qc_to_header(header, qc)
        fits.PrimaryHDU(header=header).writeto(filePath)

    # ASSERT
    cards = fits.getheader(filePath)
    assert cards["ESO QC EXTRACT METHOD"] == "boxcar"
    assert cards.comments["ESO QC EXTRACT METHOD"] == EXTRACT_METHOD_QC_COMMENT
