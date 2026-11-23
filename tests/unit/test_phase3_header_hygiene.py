"""Unit tests for the Phase 3 header rules of `basic_header_scrubbing` and `write_fits_table_to_disk` (DY-2344)."""

from __future__ import annotations

from pathlib import Path

import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.__version__ import __version__
from soxspipe.commonutils.phase3 import (
    JOINED_MODE,
    PER_ARM_MODE,
    basic_header_scrubbing,
    write_fits_table_to_disk,
)
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.unit

CCDPROC_CARDS = (
    "SUBBIAS",
    "SUBTRACT_BIAS",
    "SUBDARK",
    "FLATCOR",
    "FLAT_CORRECT",
    "GAINCOR",
    "TRIMIM",
    "CREATVAR",
    "COSMICRAY_LACOSMIC",
)
LONG_VALUE_LENGTH = 69
MAX_VALUE_LENGTH = 68


def _science_header() -> fits.Header:
    """A raw-frame-like header with one card of every group the writer treats differently."""
    header = instrument_header()
    header["TELESCOP"] = "ESO-NTT"
    header["ORIGFILE"] = "SOXS.2026-01-02T03:04:05.678.fits"
    header["BUNIT"] = "electron"
    header["SXSPRE"] = "2026-01-02T03:04:05.678901"
    for card in CCDPROC_CARDS:
        header[f"HIERARCH {card}"] = "ccdproc step"
    header["HIERARCH ESO OBS ID"] = 42
    header["HIERARCH ESO TPL START"] = "2026-01-02T03:04:05"
    header["HIERARCH ESO TEL TARG ALPHA"] = 123456.789
    header["HIERARCH ESO DET NDIT"] = 2
    header["HIERARCH ESO INS SLIT1 NAME"] = "SLIT1.0"
    header["HIERARCH ESO SEQ ARM"] = "VIS"
    header["HIERARCH ESO QC RON"] = 3.2
    header["HIERARCH ESO PRO CATG"] = "SCI_SLIT_FLUX_VIS"
    return header


def _scrub(tmp_path: Path, log: object, header: fits.Header, **kwargs: object) -> fits.Header:
    return basic_header_scrubbing(log=log, settings=pipeline_settings(tmp_path), header=header, **kwargs)


def _write(tmp_path: Path, log: object, header: fits.Header, **kwargs: object) -> Path:
    outputPath = tmp_path / "SOXS_product.fits"
    write_fits_table_to_disk(
        log=log,
        settings=pipeline_settings(tmp_path),
        header=header,
        tables=kwargs.pop("tables", [Table({"WAVE": [500.0, 500.1]})]),
        filePath=outputPath,
        **kwargs,
    )
    return outputPath


def test_per_arm_mode_keeps_eso_and_qc_cards_and_drops_dpr(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header, mode=PER_ARM_MODE)

    # ASSERT
    keptCards = [
        "ESO OBS ID",
        "ESO TPL START",
        "ESO TEL TARG ALPHA",
        "ESO DET NDIT",
        "ESO INS SLIT1 NAME",
        "ESO SEQ ARM",
        "ESO QC RON",
        "ESO PRO CATG",
    ]
    assert all(card in scrubbed for card in keptCards)
    assert not [card for card in scrubbed if card.startswith("ESO DPR")]


def test_per_arm_is_the_default_mode(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header)

    # ASSERT
    assert "ESO DET NDIT" in scrubbed
    assert "ESO QC RON" in scrubbed


def test_joined_mode_drops_det_ins_seq_and_qc_cards(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header, mode=JOINED_MODE)

    # ASSERT
    droppedPrefixes = ("ESO DET ", "ESO INS ", "ESO SEQ ", "ESO QC ")
    assert not [card for card in scrubbed if card.startswith(droppedPrefixes)]


def test_joined_mode_keeps_obs_tpl_and_tel_cards_of_the_input_header(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header, mode=JOINED_MODE)

    # ASSERT
    assert scrubbed["ESO OBS ID"] == 42
    assert scrubbed["ESO TPL START"] == "2026-01-02T03:04:05"
    assert scrubbed["ESO TEL TARG ALPHA"] == pytest.approx(123456.789)
    assert scrubbed["TELESCOP"] == "ESO-NTT"


def test_unknown_mode_raises_value_error(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="mode"):
        _scrub(tmp_path, log, header, mode="stacked")


def test_image_header_keeps_bunit_and_sxspre_by_default(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header)

    # ASSERT
    assert scrubbed["BUNIT"] == "electron"
    assert "SXSPRE" in scrubbed


def test_table_header_drops_bunit_sxspre_and_ccdproc_cards(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    scrubbed = _scrub(tmp_path, log, header, isTable=True)

    # ASSERT
    assert "BUNIT" not in scrubbed
    assert "SXSPRE" not in scrubbed
    assert not [card for card in CCDPROC_CARDS if card in scrubbed]


def test_written_table_drops_bunit_sxspre_and_ccdproc_cards_from_every_hdu(tmp_path: Path, log: object) -> None:
    # ARRANGE
    table = Table({"WAVE": [500.0, 500.1]})
    table.meta.update({"BUNIT": "nm", "SXSPRE": "2026-01-02T03:04:05.678901", "SUBBIAS": "subtracted"})

    # ACT
    outputPath = _write(tmp_path, log, _science_header(), tables=[table])

    # ASSERT
    with fits.open(outputPath) as hdus:
        for hdu in hdus:
            assert "BUNIT" not in hdu.header
            assert "SXSPRE" not in hdu.header
            assert not [card for card in CCDPROC_CARDS if card in hdu.header]


def test_written_table_keeps_eso_and_qc_cards_in_per_arm_mode(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    outputPath = _write(tmp_path, log, header, mode=PER_ARM_MODE)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert hdus[0].header["ESO QC RON"] == pytest.approx(3.2)
        assert hdus[0].header["ESO DET NDIT"] == 2
        assert "ESO DPR TYPE" not in hdus[0].header


def test_written_table_drops_det_ins_seq_and_qc_cards_in_joined_mode(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    outputPath = _write(tmp_path, log, header, mode=JOINED_MODE)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert hdus[0].header["ESO OBS ID"] == 42
        assert hdus[0].header["ESO TEL TARG ALPHA"] == pytest.approx(123456.789)
        assert not [card for card in hdus[0].header if card.startswith(("ESO DET ", "ESO INS ", "ESO SEQ "))]


def test_joined_mode_drops_the_qc_rows_the_writer_is_asked_to_add(tmp_path: Path, log: object) -> None:
    # ARRANGE
    import pandas as pd

    qc = pd.DataFrame([{"qc_name": "RON", "qc_value": 3.2, "qc_comment": "Read noise", "to_header": True}])

    # ACT
    outputPath = _write(tmp_path, log, instrument_header(), mode=JOINED_MODE, qc=qc)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert "ESO QC RON" not in hdus[0].header


def test_origfile_is_the_name_of_the_file_written(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert hdus[0].header["ORIGFILE"] == "SOXS_product.fits"


def test_named_table_has_extname_and_the_primary_has_none(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["EXTNAME"] = "RAW_EXTNAME"

    # ACT
    outputPath = _write(tmp_path, log, header, extname="SPECTRUM")

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert "EXTNAME" not in hdus[0].header
        assert hdus[1].header["EXTNAME"] == "SPECTRUM"


def test_table_is_left_unnamed_when_no_extname_is_given(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert "EXTNAME" not in hdus[0].header
        assert "EXTNAME" not in hdus[1].header


def test_identification_keywords_are_set_on_the_primary_header(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["ORIGIN"] = "SOMEWHERE"
    header["REFERENC"] = "2020A&A...000A...0X"

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    with fits.open(outputPath) as hdus:
        primary = hdus[0].header
        assert primary["ORIGIN"] == "ESO"
        assert primary["PROCSOFT"] == f"soxspipe/v{__version__}"
        assert primary["REFERENC"] == ""


def test_telescop_comes_from_the_raw_frame_header(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["TELESCOP"] = "ESO-VLT-U1"

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert hdus[0].header["TELESCOP"] == "ESO-VLT-U1"


def test_value_of_68_characters_is_written(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["OBJECT"] = "A" * MAX_VALUE_LENGTH

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    with fits.open(outputPath) as hdus:
        assert hdus[0].header["OBJECT"] == "A" * MAX_VALUE_LENGTH


def test_value_over_68_characters_raises_and_writes_no_file(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["OBJECT"] = "A" * LONG_VALUE_LENGTH

    # ACT
    with pytest.raises(ValueError, match="OBJECT"):
        _write(tmp_path, log, header)

    # ASSERT
    assert list(tmp_path.glob("*.fits")) == []


def test_long_value_on_a_table_extension_raises(tmp_path: Path, log: object) -> None:
    # ARRANGE
    table = Table({"WAVE": [500.0]})
    table.meta["HIERARCH LONG"] = "B" * LONG_VALUE_LENGTH

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="LONG"):
        _write(tmp_path, log, _science_header(), tables=[table])


def test_header_with_a_continue_card_raises(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header.append(fits.Card("CONTINUE", "'more text&'"))

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="CONTINUE"):
        _write(tmp_path, log, header)


def test_hierarch_card_that_astropy_would_write_as_continue_raises(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["HIERARCH ESO PRO REC1 RAW1 NAME"] = "C" * MAX_VALUE_LENGTH

    # ACT AND ASSERT
    with pytest.raises(ValueError, match="CONTINUE"):
        _write(tmp_path, log, header)


def test_long_comment_card_does_not_raise(tmp_path: Path, log: object) -> None:
    # ARRANGE
    header = _science_header()
    header["COMMENT"] = "D" * 70

    # ACT
    outputPath = _write(tmp_path, log, header)

    # ASSERT
    assert outputPath.exists()
