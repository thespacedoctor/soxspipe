"""Unit tests for `phase3.add_qc_to_header`, the one writer of QC rows into FITS headers."""

from __future__ import annotations

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.table import Table

from soxspipe.commonutils.horne_extraction import EXTRACT_METHOD_QC_COMMENT
from soxspipe.commonutils.phase3 import (
    TRAILING_SEPARATORS,
    add_qc_to_header,
    fit_card_comment,
    write_fits_table_to_disk,
)
from soxspipe.commonutils.toolkit import append_qc
from tests.factories import pipeline_settings

pytestmark = pytest.mark.unit

CARD_LENGTH = 80
CARD_TOO_LONG = "Card is too long"
SHORT_COMMENT = "read noise"
# "HIERARCH ESO QC LINES = 42 / " LEAVES EXACTLY 51 CHARACTERS
FITTING_COMMENT = "number of arc lines detected and fitted in this arm"
# ASTROPY CUTS THIS ONE MID-WORD, AFTER "ALONG THE DIS"
CUT_COMMENT = "Shift of the sky lines along the dispersion axis in pixels"
CUT_AT_WORD = "Shift of the sky lines along the"
LONG_NAME = "INNER ORDER PIXEL MEAN FLUX LEVEL OF THE WHOLE FRAME"
LONG_NAME_COMMENT = "description of the frame"
FITTING_ROWS = (
    ("RON", 3.2, SHORT_COMMENT),
    ("LINES", 42, FITTING_COMMENT),
    ("SKY SHIFT O15 PX", -0.034, CUT_COMMENT),
    (LONG_NAME, 1, LONG_NAME_COMMENT),
)


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


def _commented_qc(rows: tuple[tuple[str, object, str], ...]) -> pd.DataFrame:
    # AN OBJECT COLUMN FROM THE START, SO AN INTEGER VALUE STAYS AN INTEGER AND IS NOT WRITTEN AS "42.0"
    return pd.DataFrame(
        {
            "qc_name": [n for n, _, _ in rows],
            "qc_value": pd.Series([v for _, v, _ in rows], dtype=object),
            "qc_comment": [c for _, _, c in rows],
            "to_header": [True] * len(rows),
        },
    )


def _is_mid_word_cut(written: str, full: str) -> bool:
    """True when ``written`` is a shortened ``full`` that does not end at a word boundary."""
    if written in ("", full):
        return False
    # THE FITTER DROPS A SEPARATOR LEFT AT THE CUT, SO "AXIS" IS A WORD-BOUNDARY CUT OF "AXIS, IN"
    rest = full[len(written) :]
    return not (full.startswith(written) and rest[:1] in tuple(TRAILING_SEPARATORS))


def test_add_qc_to_header_fits_each_comment_into_the_80_character_card() -> None:
    # ARRANGE
    header = fits.Header()
    qc = _commented_qc(FITTING_ROWS)
    assert len(f"HIERARCH ESO QC LINES = 42 / {FITTING_COMMENT}") == CARD_LENGTH

    # ACT
    add_qc_to_header(header, qc)

    # ASSERT
    assert header.comments["ESO QC RON"] == SHORT_COMMENT
    assert header.comments["ESO QC LINES"] == FITTING_COMMENT
    assert header.comments["ESO QC SKY SHIFT O15 PX"] == CUT_AT_WORD
    assert header.comments[f"ESO QC {LONG_NAME}"] == ""


def test_qc_cards_written_to_a_file_raise_no_card_too_long_warning_and_keep_their_values(tmp_path: Path) -> None:
    # ARRANGE
    header = fits.Header()
    qc = _commented_qc(FITTING_ROWS)
    filePath = tmp_path / "qc_cards.fits"

    # ACT
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        add_qc_to_header(header, qc)
        fits.PrimaryHDU(header=header).writeto(filePath)

    # ASSERT
    tooLong = [w for w in raised if issubclass(w.category, VerifyWarning) and CARD_TOO_LONG in str(w.message)]
    assert tooLong == []
    cards = fits.getheader(filePath)
    for name, value, _ in FITTING_ROWS:
        assert cards[f"ESO QC {name}"] == pytest.approx(value)


def test_add_qc_to_header_leaves_the_qc_table_comments_unchanged() -> None:
    # ARRANGE
    qc = _commented_qc(FITTING_ROWS)

    # ACT
    add_qc_to_header(fits.Header(), qc)

    # ASSERT
    assert list(qc["qc_comment"]) == [c for _, _, c in FITTING_ROWS]


def test_a_written_product_holds_no_qc_comment_cut_mid_word(tmp_path: Path, log: object) -> None:
    # ARRANGE
    qc = _commented_qc(FITTING_ROWS)
    filePath = tmp_path / "SOXS_product.fits"

    # ACT
    write_fits_table_to_disk(
        log=log,
        settings=pipeline_settings(tmp_path),
        header=fits.Header(),
        tables=[Table({"WAVE": [500.0, 500.1]})],
        filePath=filePath,
        qc=qc,
    )

    # ASSERT
    cards = fits.getheader(filePath)
    for name, value, comment in FITTING_ROWS:
        assert cards[f"ESO QC {name}"] == pytest.approx(value)
        assert not _is_mid_word_cut(cards.comments[f"ESO QC {name}"], comment)


@pytest.mark.parametrize("comment", [None, ""])
def test_fit_card_comment_returns_a_missing_comment_unchanged(comment: str | None) -> None:
    # ACT AND ASSERT
    assert fit_card_comment("ESO QC RON", 3.2, comment) == comment


def test_fit_card_comment_drops_a_separator_left_at_the_cut() -> None:
    # ARRANGE
    # "HIERARCH ESO QC SKY SHIFT O15 PX = -0.034 / " LEAVES 36 CHARACTERS; THE CUT FALLS AFTER "AXIS,"
    comment = "Sky shift measured along the axis, in pixels of the detector"

    # ACT
    fitted = fit_card_comment("ESO QC SKY SHIFT O15 PX", -0.034, comment)

    # ASSERT
    assert fitted == "Sky shift measured along the axis"


@pytest.mark.parametrize(
    ("name", "value"),
    [
        # A STRING TOO LONG FOR ONE CARD IS WRITTEN AS CONTINUE CARDS
        ("EXTRACT METHOD", "x" * 60),
        # "HIERARCH ESO QC <NAME> = 12345678 / X" IS 81 CHARACTERS, SO EVEN THE PLACEHOLDER IS CUT
        ("A LONG NAME OF A QC THAT IS LONG ENOUGH TO FILL IT", 12345678),
    ],
)
def test_fit_card_comment_returns_none_when_the_value_fills_the_card(name: str, value: object) -> None:
    # ACT
    fitted = fit_card_comment(f"ESO QC {name}", value, "method")

    # ASSERT
    assert fitted is None


def test_a_short_string_qc_card_raises_no_card_too_long_warning(tmp_path: Path) -> None:
    """Astropy pads a short string value under a short keyword, so the comment room is smaller than it looks."""
    # ARRANGE
    header = fits.Header()
    comment = "alpha beta gamma delta epsilon zeta eta theta iota kappa lambda"
    qc = _commented_qc((("S", "boxcar", comment),))
    filePath = tmp_path / "string_qc.fits"

    # ACT
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        add_qc_to_header(header, qc)
        fits.PrimaryHDU(header=header).writeto(filePath)

    # ASSERT
    tooLong = [w for w in raised if issubclass(w.category, VerifyWarning) and CARD_TOO_LONG in str(w.message)]
    assert tooLong == []
    cards = fits.getheader(filePath)
    assert cards["ESO QC S"] == "boxcar"
    assert not _is_mid_word_cut(cards.comments["ESO QC S"], comment)
