"""The wavelength-calibration and resolution cards of a merged spectrum (DY-2358).

`LAMRMS`, `LAMNLIN` and `SPEC_ERR` are read from the `WAVE RMS` and `WAVE NLINES` QCs of the `DISP_TAB` the science
frames were reduced with. `SPEC_RES` comes from the nominal `spec-res` table of the detector parameters, by arm and
slit. A mandatory card that cannot be computed is left out with one warning, which is how an arm fails (DY-2336).
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.wavelength_keywords import (
    spec_res_cards,
    wavelength_accuracy_cards,
    with_wavelength_keywords,
)

pytestmark = pytest.mark.unit

SOXS_SETTINGS = {"instrument": "soxs"}
XSH_SETTINGS = {"instrument": "xsh"}
WAVE_RMS_NM = 0.0036
WAVE_NLINES = 121


def _values(cards: list[tuple[str, Any, str]]) -> dict[str, Any]:
    return {keyword: value for keyword, value, _ in cards}


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _disp_tab(tmp_path: Path, **qcCards: Any) -> str:
    """Write a `DISP_TAB` with the given `ESO QC` cards, and return its path."""
    header = fits.Header({f"HIERARCH ESO QC {name}": value for name, value in qcCards.items()})
    path = tmp_path / "SOXS_DISP_TAB_VIS.fits"
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU.from_columns([])]).writeto(path)
    return str(path)


def _nominal(arm: str, slit: str) -> int:
    log = type("Quiet", (), {"debug": staticmethod(lambda *a, **k: None)})()
    return detector_lookup(log=log, settings=SOXS_SETTINGS).get(arm)["spec-res"][slit]


def test_the_accuracy_cards_come_from_the_wave_qcs_of_the_dispersion_map(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES})

    # ACT
    cards = wavelength_accuracy_cards(path, log)

    # ASSERT
    values = _values(cards)
    assert values["LAMRMS"] == WAVE_RMS_NM
    assert values["LAMNLIN"] == WAVE_NLINES
    assert isinstance(values["LAMNLIN"], int)
    assert values["SPEC_ERR"] == pytest.approx(WAVE_RMS_NM / math.sqrt(WAVE_NLINES), rel=1e-12)
    assert _warnings(log) == []


def test_a_float_line_count_is_written_as_an_integer(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE QC TABLE HOLDS ITS NUMBERS AS FLOATS, SO THE CARD MAY READ 121.0
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": 121.0})

    # ACT
    values = _values(wavelength_accuracy_cards(path, log))

    # ASSERT
    assert values["LAMNLIN"] == 121
    assert isinstance(values["LAMNLIN"], int)


def test_a_dispersion_map_without_the_wave_qcs_gives_no_cards_and_one_warning(tmp_path: Path, log: Any) -> None:
    # ARRANGE: A MAP MADE BEFORE THE QCs EXISTED HAS OTHER QCs BUT NOT THESE
    path = _disp_tab(tmp_path, **{"X RES MEDIAN": 0.1})

    # ACT
    cards = wavelength_accuracy_cards(path, log)

    # ASSERT
    assert cards == []
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "WAVE RMS" in warnings[0] and "LAMRMS" in warnings[0] and "SPEC_ERR" in warnings[0]


@pytest.mark.parametrize("path", [None, False, ""])
def test_no_dispersion_map_gives_no_cards_and_one_warning(log: Any, path: Any) -> None:
    # ACT
    cards = wavelength_accuracy_cards(path, log)

    # ASSERT
    assert cards == []
    assert len(_warnings(log)) == 1


def test_an_unreadable_dispersion_map_gives_no_cards_and_one_warning(tmp_path: Path, log: Any) -> None:
    # ACT
    cards = wavelength_accuracy_cards(str(tmp_path / "missing.fits"), log)

    # ASSERT
    assert cards == []
    assert len(_warnings(log)) == 1


@pytest.mark.parametrize("nlines", [0, -3])
def test_a_line_count_below_one_gives_no_cards_and_one_warning(tmp_path: Path, log: Any, nlines: int) -> None:
    # ARRANGE
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": nlines})

    # ACT
    cards = wavelength_accuracy_cards(path, log)

    # ASSERT
    assert cards == []
    assert len(_warnings(log)) == 1


def test_a_wave_qc_that_is_not_a_number_gives_no_cards_and_one_warning(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    path = _disp_tab(tmp_path, **{"WAVE RMS": "n/a", "WAVE NLINES": WAVE_NLINES})

    # ACT
    cards = wavelength_accuracy_cards(path, log)

    # ASSERT
    assert cards == []
    assert len(_warnings(log)) == 1


@pytest.mark.parametrize(
    ("arm", "slit"), [("VIS", "SLIT1.0"), ("VIS", "SLIT5.0"), ("NIR", "SLIT0.5"), ("NIR", "SLIT1.5")]
)
def test_spec_res_is_the_nominal_value_of_the_arm_and_slit(log: Any, arm: str, slit: str) -> None:
    # ACT
    cards = spec_res_cards(arm, slit, SOXS_SETTINGS, log)

    # ASSERT
    assert _values(cards) == {"SPEC_RES": _nominal(arm, slit)}
    assert _warnings(log) == []


def test_the_nir_slit_name_may_carry_the_jh_filter_suffix(log: Any) -> None:
    # ACT
    cards = spec_res_cards("NIR", "SLIT1.5JH", SOXS_SETTINGS, log)

    # ASSERT
    assert _values(cards) == {"SPEC_RES": _nominal("NIR", "SLIT1.5")}


@pytest.mark.parametrize("slit", ["BLANK", "UNDEFINED", "SLIT9.9", "", None])
def test_a_slit_with_no_table_row_gives_no_spec_res_and_one_warning_that_fails_the_arm(log: Any, slit: Any) -> None:
    # ACT
    cards = spec_res_cards("NIR", slit, SOXS_SETTINGS, log)

    # ASSERT
    assert cards == []
    warnings = _warnings(log)
    assert len(warnings) == 1
    assert "SPEC_RES" in warnings[0] and "arm fails" in warnings[0]


def test_spec_res_is_never_the_pinhole_resolving_power(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE MAP CARRIES A PINHOLE RESOLVING POWER, WHICH IS FAR ABOVE THE SLIT VALUE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES, "R PIN MEDIAN": 9400.0})
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0"})

    # ACT
    stamped = with_wavelength_keywords(header, path, "VIS", kw, SOXS_SETTINGS, log)

    # ASSERT
    assert stamped["SPEC_RES"] == _nominal("VIS", "SLIT1.0")
    assert stamped["SPEC_RES"] != 9400


def test_the_stamp_gives_a_copy_with_every_card_and_leaves_the_header_alone(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES})
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.5"})

    # ACT
    stamped = with_wavelength_keywords(header, path, "VIS", kw, SOXS_SETTINGS, log)

    # ASSERT
    assert [key for key in ("SPEC_RES", "LAMRMS", "LAMNLIN", "SPEC_ERR") if key in header] == []
    assert stamped["SPEC_RES"] == _nominal("VIS", "SLIT1.5")
    assert stamped["LAMRMS"] == WAVE_RMS_NM
    assert stamped["LAMNLIN"] == WAVE_NLINES
    assert all(stamped.comments[key] for key in ("SPEC_RES", "LAMRMS", "LAMNLIN", "SPEC_ERR"))


def test_the_stamp_removes_a_stale_card_that_it_cannot_recompute(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE HEADER ALREADY HOLDS CARDS FROM AN EARLIER STAMP, AND THIS MAP HAS NO WAVE QCs
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    path = _disp_tab(tmp_path)
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0", "LAMRMS": 1.0, "LAMNLIN": 5, "SPEC_ERR": 0.4})

    # ACT
    stamped = with_wavelength_keywords(header, path, "VIS", kw, SOXS_SETTINGS, log)

    # ASSERT
    assert [key for key in ("LAMRMS", "LAMNLIN", "SPEC_ERR") if key in stamped] == []
    assert "SPEC_RES" in stamped


def test_the_stamp_reads_the_slit_card_of_the_arm(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    path = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES})
    header = fits.Header({kw("SLIT_VIS"): "SLIT1.0", kw("SLIT_NIR"): "SLIT0.5"})

    # ACT
    nir = with_wavelength_keywords(header, path, "NIR", kw, SOXS_SETTINGS, log)

    # ASSERT
    assert nir["SPEC_RES"] == _nominal("NIR", "SLIT0.5")


def test_x_shooter_keeps_its_header_and_gets_no_warning(tmp_path: Path, log: Any) -> None:
    # ARRANGE: X-SHOOTER HAS NO spec-res TABLE, AND ITS CODE PATHS DO NOT CHANGE
    kw = keyword_lookup(log=log, settings=XSH_SETTINGS).get
    header = fits.Header({"OBJECT": "star"})

    # ACT
    stamped = with_wavelength_keywords(header, None, "VIS", kw, XSH_SETTINGS, log)

    # ASSERT
    assert stamped == header
    assert stamped is not header
    assert _warnings(log) == []


def test_the_cards_survive_a_strict_phase_3_write_with_their_comments_whole(tmp_path: Path, log: Any) -> None:
    from astropy.table import Table

    from soxspipe.commonutils.phase3 import PER_ARM_MODE, write_fits_table_to_disk

    # ARRANGE: A REAL TABLE FILE WRITTEN THE WAY A PHASE 3 PRODUCT IS, WHERE A CARD THAT DOES NOT FIT RAISES
    kw = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    dispMap = _disp_tab(tmp_path, **{"WAVE RMS": WAVE_RMS_NM, "WAVE NLINES": WAVE_NLINES})
    header = with_wavelength_keywords(fits.Header({kw("SLIT_VIS"): "SLIT1.0"}), dispMap, "VIS", kw, SOXS_SETTINGS, log)
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
    for card in ("SPEC_RES", "LAMRMS", "LAMNLIN", "SPEC_ERR"):
        assert written.comments[card] == header.comments[card]
    assert written["LAMNLIN"] == WAVE_NLINES
    assert written["SPEC_ERR"] == pytest.approx(WAVE_RMS_NM / math.sqrt(WAVE_NLINES), rel=1e-12)
