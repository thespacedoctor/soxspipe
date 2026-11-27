"""The ``SPEC_SYE`` card of a merged spectrum: the largest ``SKY SHIFT O<n> ERR`` over its orders (DY-2359).

The QC rows come from the extraction (one per order, already the largest over the frames and cycles of a stack). The
card is left out, with one warning, when there is no row to take it from, which is the case for every NIR spectrum
until the NIR error formulas exist.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.toolkit import append_qc
from soxspipe.commonutils.wavelength_keywords import (
    spec_sye_cards,
    with_merged_spectrum_cards,
    with_wavelength_keywords,
)
from tests.factories.disp_tab import disp_tab_with_qcs

pytestmark = pytest.mark.unit

SOXS_SETTINGS = {"instrument": "soxs"}
XSH_SETTINGS = {"instrument": "xsh"}
DATE = "2026-01-02T03:04:05"


def _qc(errors: dict[int, Any] | None = None, **others: Any) -> pd.DataFrame:
    """A QC table with a ``SKY SHIFT O<n> ERR`` row for each order in ``errors``, and any other QC given by name."""
    table = pd.DataFrame()
    for order, value in (errors or {}).items():
        table = append_qc(
            table,
            "soxs-stare",
            f"SKY SHIFT O{order} ERR",
            value,
            "c",
            DATE,
            DATE,
            qcUnit="nm",
            toHeader=True,
            qcOrder=order,
        )
    for name, value in others.items():
        table = append_qc(table, "soxs-stare", name.replace("_", " "), value, "c", DATE, DATE, toHeader=True)
    return table


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_spec_sye_is_the_largest_sky_shift_error_over_the_orders(log: Any) -> None:
    # ARRANGE
    qc = _qc({1: 0.0182, 2: 0.024, 3: 0.0285, 4: 0.0031})

    # ACT
    cards = spec_sye_cards(qc, "VIS", log)

    # ASSERT
    assert [(keyword, value) for keyword, value, _ in cards] == [("SPEC_SYE", 0.0285)]
    assert "nm" in cards[0][2]
    assert _warnings(log) == []


def test_spec_sye_ignores_every_other_qc_including_the_shift_and_source_cards(log: Any) -> None:
    # ARRANGE
    qc = _qc({4: 0.003}, **{"SKY_SHIFT_O4": 0.9, "SKY_SHIFT_O4_SRC": "measured", "SKY_SHIFT_O4_PX": 7.0, "RON": 50.0})

    # ACT
    cards = spec_sye_cards(qc, "VIS", log)

    # ASSERT
    assert cards[0][1] == 0.003


def test_spec_sye_ignores_an_error_that_is_not_a_finite_number(log: Any) -> None:
    # ARRANGE
    qc = _qc({1: float("nan"), 2: "n/a", 3: 0.02})

    # ACT
    cards = spec_sye_cards(qc, "VIS", log)

    # ASSERT
    assert cards[0][1] == 0.02


@pytest.mark.parametrize(
    "qc",
    [None, False, pd.DataFrame(), _qc(RON=1.0), _qc({1: float("nan")})],
    ids=["none", "false", "empty", "no-error-rows", "no-finite-error"],
)
def test_no_error_row_gives_no_card_and_one_warning(log: Any, qc: Any) -> None:
    # ACT
    cards = spec_sye_cards(qc, "NIR", log)

    # ASSERT
    assert cards == []
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "SPEC_SYE" in warnings[0] and "NIR" in warnings[0] and "SKY SHIFT O<n> ERR" in warnings[0]


def test_the_stamp_gives_a_copy_with_spec_sye_and_leaves_the_header_alone(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0"})
    dispMap = disp_tab_with_qcs(tmp_path, **{"WAVE RMS": 0.0036, "WAVE NLINES": 121})

    # ACT
    stamped = with_wavelength_keywords(header, dispMap, "VIS", kw, SOXS_SETTINGS, log, qc=_qc({1: 0.01, 4: 0.02}))

    # ASSERT
    assert "SPEC_SYE" not in header
    assert stamped["SPEC_SYE"] == 0.02
    assert stamped.comments["SPEC_SYE"]


def test_the_stamp_removes_a_stale_spec_sye_that_it_cannot_recompute(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0", "SPEC_SYE": 9.0})

    # ACT
    stamped = with_wavelength_keywords(header, None, "NIR", kw, SOXS_SETTINGS, log, qc=_qc())

    # ASSERT
    assert "SPEC_SYE" not in stamped


def test_x_shooter_gets_no_spec_sye_and_no_warning(log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=XSH_SETTINGS).get
    header = fits.Header({"OBJECT": "star"})

    # ACT
    stamped = with_wavelength_keywords(header, None, "VIS", kw, XSH_SETTINGS, log, qc=_qc({1: 0.01}))

    # ASSERT
    assert stamped == header
    assert _warnings(log) == []


def test_the_merged_spectrum_stamp_passes_the_qc_on(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0"})

    # ACT
    stamped = with_merged_spectrum_cards(header, None, None, "VIS", kw, SOXS_SETTINGS, log, qc=_qc({2: 0.024}))

    # ASSERT
    assert stamped["SPEC_SYE"] == 0.024


def test_the_card_survives_a_strict_phase_3_write_with_its_comment_whole(tmp_path: Path, log: Any) -> None:
    from astropy.table import Table

    from soxspipe.commonutils.phase3 import PER_ARM_MODE, write_fits_table_to_disk

    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    header = with_wavelength_keywords(
        fits.Header({kw("SLIT_VIS"): "SLIT1.0"}), None, "VIS", kw, SOXS_SETTINGS, log, qc=_qc({3: 0.0285})
    )
    filePath = tmp_path / "SOXS_20260102T030405_STARE_ph3.fits"

    # ACT
    write_fits_table_to_disk(
        log=log,
        settings=SOXS_SETTINGS,
        header=header,
        tables=[Table({"WAVE": [500.0]})],
        filePath=filePath,
        phase3=PER_ARM_MODE,
        extname="SPECTRUM",
        strict=True,
    )

    # ASSERT
    written = fits.getheader(filePath)
    assert written.comments["SPEC_SYE"] == header.comments["SPEC_SYE"]
    assert math.isclose(written["SPEC_SYE"], 0.0285, rel_tol=1e-12)
