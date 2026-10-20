"""Contracts for skipping the Horne-extraction skyline tuning when a SOXS arc is attached (DY-1665, DY-1666)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.extraction_notes import extraction_note_logger
from soxspipe.commonutils.ug_wavelength_transfer import (
    MAX_ATTACHED_ARC_MINUTES,
    MAX_ATTACHED_ARC_ROTATOR_DEGREES,
    attached_arc_offsets,
    is_attached_arc,
)
from tests.unit._horne_vis_scene import KEYWORDS, build_extractor, qc_by_name

pytestmark = pytest.mark.unit

MINUTES_PER_DAY = 1440.0
# SCIENCE FRAME OF THE REAL-DATA CHECK (20260626T044348, EXPTIME 180 S) AND ITS ARC (EXPTIME 30 S)
SCIENCE_MJD = 61217.19708617
SCIENCE_ABSROT = -0.57699
SCIENCE_EXPTIME = 180.0
ARC_EXPTIME = 30.0
SKIPPED_QC = "skipped: attached arc"
TUNING_QC = "SKY TUNING"


def _header(
    minutesFromScience: float = 0.0,
    absrot: float | None = SCIENCE_ABSROT,
    mjd: bool = True,
    exptime: float | None = SCIENCE_EXPTIME,
) -> fits.Header:
    header = fits.Header()
    if mjd:
        header[KEYWORDS["MJDOBS"]] = SCIENCE_MJD + minutesFromScience / MINUTES_PER_DAY
    if exptime is not None:
        header[KEYWORDS["EXPTIME"]] = exptime
    if absrot is not None:
        header[KEYWORDS["ABSROT"]] = absrot
        header[KEYWORDS["PAR_ANG_START"]] = -85.867
    return header


def _arc(gapMinutes: float = 3.0, degrees: float = 1.0, after: bool = True, **kwargs: object) -> fits.Header:
    """An arc ``gapMinutes`` after the end of the science frame, or ending ``gapMinutes`` before its start."""
    start = SCIENCE_EXPTIME / 60.0 + gapMinutes if after else -gapMinutes - ARC_EXPTIME / 60.0
    return _header(minutesFromScience=start, absrot=SCIENCE_ABSROT + degrees, exptime=ARC_EXPTIME, **kwargs)


def _offsets(arc: fits.Header, science: fits.Header | list[fits.Header] | None = None) -> tuple[float, float] | None:
    scienceHeaders = science if isinstance(science, list) else [_header() if science is None else science]
    return attached_arc_offsets(scienceHeaders, arc, KEYWORDS.__getitem__)


def _attached(arc: fits.Header, science: fits.Header | list[fits.Header] | None = None) -> bool:
    scienceHeaders = science if isinstance(science, list) else [_header() if science is None else science]
    return is_attached_arc(scienceHeaders, arc, KEYWORDS.__getitem__)


def _nod_block(frameCount: int = 4, cadenceMinutes: float = 4.0, degreesPerFrame: float = 2.0) -> list[fits.Header]:
    """A nod sequence of 180 s frames, one every ``cadenceMinutes``, with the rotator turning between frames."""
    return [
        _header(minutesFromScience=index * cadenceMinutes, absrot=SCIENCE_ABSROT + index * degreesPerFrame)
        for index in range(frameCount)
    ]


def test_an_arc_3_minutes_and_1_degree_away_is_attached() -> None:
    assert _attached(_arc(3.0, 1.0)) is True


def test_the_real_data_arc_after_its_science_frame_is_1_19_minutes_away() -> None:
    # ARRANGE: ARC 20260626T044759 STARTS 4.19 MIN AFTER THE 180 S SCIENCE FRAME 20260626T044348 STARTS
    arc = _header(minutesFromScience=4.19, absrot=SCIENCE_ABSROT + 1.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, _ = _offsets(arc)

    # ASSERT
    assert minutes == pytest.approx(1.19, abs=1e-4)


def test_an_arc_before_the_science_frame_is_measured_from_the_arc_end() -> None:
    # ARRANGE: THE 30 S ARC STARTS 4 MIN BEFORE THE SCIENCE FRAME, SO IT ENDS 3.5 MIN BEFORE IT
    arc = _header(minutesFromScience=-4.0, absrot=SCIENCE_ABSROT + 1.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, _ = _offsets(arc)

    # ASSERT
    assert minutes == pytest.approx(3.5, abs=1e-4)


def test_the_offsets_are_the_gap_minutes_and_wrap_aware_degrees() -> None:
    # ARRANGE
    science = _header(absrot=179.0)
    arc = _header(minutesFromScience=-4.5, absrot=-179.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, degrees = _offsets(arc, science)

    # ASSERT
    assert minutes == pytest.approx(4.0, abs=1e-3)
    assert degrees == pytest.approx(2.0)


@pytest.mark.parametrize("after", [True, False])
def test_an_arc_exactly_at_the_minutes_limit_is_attached(after: bool) -> None:
    assert _attached(_arc(MAX_ATTACHED_ARC_MINUTES, 1.0, after=after)) is True


@pytest.mark.parametrize("after", [True, False])
def test_an_arc_just_beyond_the_minutes_limit_is_not_attached(after: bool) -> None:
    assert _attached(_arc(MAX_ATTACHED_ARC_MINUTES + 0.01, 1.0, after=after)) is False


def test_an_arc_7_5_minutes_after_the_science_start_is_attached_when_the_frame_is_long() -> None:
    # ARRANGE: START TO START IS BEYOND THE LIMIT, END TO START (4.5 MIN) IS NOT
    arc = _header(minutesFromScience=7.5, absrot=SCIENCE_ABSROT + 1.0, exptime=ARC_EXPTIME)

    # ACT / ASSERT
    assert _attached(arc) is True


def test_an_arc_that_overlaps_the_science_frame_is_0_minutes_away() -> None:
    # ARRANGE: THE ARC STARTS 1 MIN INTO THE 3 MIN SCIENCE EXPOSURE
    arc = _header(minutesFromScience=1.0, absrot=SCIENCE_ABSROT + 1.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, _ = _offsets(arc)

    # ASSERT
    assert minutes == 0.0


def test_an_arc_after_a_long_nod_sequence_is_measured_from_the_last_frame_end() -> None:
    # ARRANGE: FOUR FRAMES OVER 15 MIN; THE ARC STARTS 2 MIN AFTER THE LAST ENDS, 17 MIN AFTER THE FIRST STARTS
    block = _nod_block()
    arc = _header(minutesFromScience=12.0 + 3.0 + 2.0, absrot=SCIENCE_ABSROT + 6.0 + 1.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, degrees = _offsets(arc, block)

    # ASSERT
    assert minutes == pytest.approx(2.0, abs=1e-4)
    assert degrees == pytest.approx(1.0)
    assert _attached(arc, block) is True


def test_an_arc_before_a_nod_sequence_is_compared_with_the_first_frame() -> None:
    # ARRANGE: THE ARC ENDS 2 MIN BEFORE THE FIRST FRAME, 1 DEGREE FROM IT AND 7 DEGREES FROM THE LAST
    block = _nod_block()
    arc = _arc(2.0, 1.0, after=False)

    # ACT
    minutes, degrees = _offsets(arc, block)

    # ASSERT
    assert minutes == pytest.approx(2.0, abs=1e-4)
    assert degrees == pytest.approx(1.0)


def test_an_arc_inside_a_nod_sequence_is_0_minutes_away_and_compared_with_the_nearest_frame() -> None:
    # ARRANGE: THE ARC (11.5 TO 12 MIN) FALLS BETWEEN FRAME 3 (ENDS AT 11 MIN) AND FRAME 4 (STARTS AT 12 MIN)
    block = _nod_block()
    arc = _header(minutesFromScience=11.5, absrot=SCIENCE_ABSROT + 6.0, exptime=ARC_EXPTIME)

    # ACT
    minutes, degrees = _offsets(arc, block)

    # ASSERT
    assert minutes == 0.0
    assert degrees == pytest.approx(0.0)


def test_the_science_block_frames_may_come_in_any_order() -> None:
    # ARRANGE
    block = _nod_block()
    arc = _header(minutesFromScience=17.0, absrot=SCIENCE_ABSROT + 7.0, exptime=ARC_EXPTIME)

    # ACT / ASSERT
    assert _offsets(arc, list(reversed(block))) == pytest.approx(_offsets(arc, block))


def test_a_nir_frame_ends_after_its_whole_sequence_not_after_one_dit() -> None:
    # ARRANGE: A REAL NIR FRAME (2026-01-08T02:07:54) HAS EXPTIME 600 S (ONE DIT) AND A 2413 S SEQUENCE (NDIT 4)
    science = _header(exptime=600.0)
    science[KEYWORDS["SEQ_EXPTIME"]] = 2413.0882
    arc = _header(minutesFromScience=2413.0882 / 60.0 + 2.0, absrot=SCIENCE_ABSROT + 1.0, exptime=15.0)

    # ACT
    minutes, _ = _offsets(arc, science)

    # ASSERT
    assert minutes == pytest.approx(2.0, abs=1e-4)


def test_a_nir_arc_ends_after_its_whole_sequence() -> None:
    # ARRANGE: THE ARC'S SEQUENCE LASTS 120 S, ITS EXPTIME (ONE DIT) 15 S; IT STARTS 6 MIN BEFORE THE SCIENCE
    arc = _header(minutesFromScience=-6.0, absrot=SCIENCE_ABSROT + 1.0, exptime=15.0)
    arc[KEYWORDS["SEQ_EXPTIME"]] = 120.0

    # ACT
    minutes, _ = _offsets(arc)

    # ASSERT
    assert minutes == pytest.approx(4.0, abs=1e-4)


def test_an_arc_at_the_rotator_limit_is_not_attached() -> None:
    assert _attached(_arc(3.0, MAX_ATTACHED_ARC_ROTATOR_DEGREES)) is False


def test_an_arc_beyond_the_rotator_limit_is_not_attached() -> None:
    assert _attached(_arc(3.0, 6.0)) is False


def test_rotator_angles_either_side_of_the_wrap_are_attached() -> None:
    # ARRANGE: THE ARC STARTS 3 MIN AFTER THE 3 MIN SCIENCE FRAME ENDS
    arc = _header(minutesFromScience=6.0, absrot=-179.0, exptime=ARC_EXPTIME)

    # ACT / ASSERT
    assert _attached(arc, science=_header(absrot=179.0)) is True


@pytest.mark.parametrize("after", [True, False])
def test_an_arc_before_or_after_the_science_frame_is_attached(after: bool) -> None:
    assert _attached(_arc(3.0, 1.0, after=after)) is True


def test_an_arc_without_a_rotator_angle_is_not_attached() -> None:
    assert _attached(_header(minutesFromScience=6.0, absrot=None, exptime=ARC_EXPTIME)) is False


def test_a_science_frame_without_a_rotator_angle_is_not_attached() -> None:
    assert _attached(_arc(3.0, 1.0), science=_header(absrot=None)) is False


def test_a_science_block_with_one_frame_without_a_rotator_angle_is_not_attached() -> None:
    # ARRANGE
    block = _nod_block()
    block[1] = _header(minutesFromScience=4.0, absrot=None)

    # ACT / ASSERT
    assert _offsets(_arc(2.0, 1.0, after=False), block) is None


def test_an_empty_science_block_is_not_attached() -> None:
    assert _offsets(_arc(3.0, 1.0), []) is None


def test_an_arc_with_the_parked_placeholder_rotator_angle_is_not_attached() -> None:
    # ARRANGE
    arc = _header(minutesFromScience=6.0, absrot=0.0, exptime=ARC_EXPTIME)
    arc[KEYWORDS["PAR_ANG_START"]] = 0.0

    # ACT / ASSERT
    assert _attached(arc, science=_header(absrot=1.0)) is False


@pytest.mark.parametrize("missingOn", ["arc", "science"])
@pytest.mark.parametrize("key", ["MJDOBS", "EXPTIME"])
def test_a_frame_without_an_observation_date_or_exposure_time_is_not_attached(missingOn: str, key: str) -> None:
    # ARRANGE
    arc = _arc(3.0, 1.0)
    science = _header()
    del (arc if missingOn == "arc" else science)[KEYWORDS[key]]

    # ACT / ASSERT
    assert _offsets(arc, science) is None
    assert _attached(arc, science=science) is False


# A FITS HEADER CANNOT HOLD NaN OR inf, SO A BOOLEAN AND A STRING ARE THE BAD VALUES REAL HEADERS CAN CARRY
@pytest.mark.parametrize("value", [True, "UNKNOWN"])
@pytest.mark.parametrize("key", ["MJDOBS", "ABSROT", "EXPTIME"])
@pytest.mark.parametrize("badOn", ["arc", "science"])
def test_a_boolean_or_text_date_angle_or_exposure_time_is_not_attached(badOn: str, key: str, value: object) -> None:
    # ARRANGE
    arc = _arc(3.0, 1.0)
    science = _header()
    (arc if badOn == "arc" else science)[KEYWORDS[key]] = value

    # ACT / ASSERT
    assert _offsets(arc, science) is None
    assert _attached(arc, science=science) is False


def test_a_negative_exposure_time_is_not_attached() -> None:
    assert _offsets(_arc(3.0, 1.0), _header(exptime=-1.0)) is None


def _gate_extractor(
    log: object,
    tmp_path: Path,
    arc: fits.Header,
    inst: str = "SOXS",
    science: fits.Header | None = None,
    block: list[fits.Header] | None = None,
) -> object:
    tmp_path.mkdir(parents=True, exist_ok=True)
    scienceHeader = _header() if science is None else science
    extractor = build_extractor(log, tmp_path, dict(scienceHeader), inst=inst)
    extractor.skySubtractedFrame.header = scienceHeader
    extractor.scienceHeaders = [scienceHeader] if block is None else block
    with fits.open(extractor.dispersionMap, mode="update") as hdus:
        hdus[0].header.update(arc)
    extractor.calls = []

    def spy(extractedOrdersDF: pd.DataFrame, arm: str, byOrder: bool = True) -> pd.DataFrame:
        extractor.calls.append(arm)
        return extractedOrdersDF.assign(tuned=True)

    extractor.tune_wavelength_calibration_to_skylines = spy
    return extractor


def _spectrum() -> pd.DataFrame:
    return pd.DataFrame({"order": [1, 1], "wavelengthMean": [500.0, 500.05]})


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_an_attached_soxs_arc_skips_tuning_and_records_it(log: object, tmp_path: Path, arm: str) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(4.0, 0.5))
    spectrum = _spectrum()

    # ACT
    result = extractor._tune_or_skip(spectrum, arm)

    # ASSERT
    assert extractor.calls == []
    assert "tuned" not in result.columns
    assert qc_by_name(extractor, TUNING_QC)["qc_value"] == SKIPPED_QC
    assert not [name for name in extractor.qc["qc_name"] if name.startswith("SKY SHIFT")]
    printed = [message for level, message in log.messages if level == "print"]
    assert any("Attached arc (Δt 4.0 min, Δrot 0.5°): skyline wavelength tuning skipped" in m for m in printed)


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_a_soxs_arc_that_is_not_attached_is_tuned(log: object, tmp_path: Path, arm: str) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(30.0, 0.5))

    # ACT
    result = extractor._tune_or_skip(_spectrum(), arm)

    # ASSERT
    assert extractor.calls == [arm]
    assert result["tuned"].all()
    assert qc_by_name(extractor, TUNING_QC)["qc_value"] == "skylines"


def test_the_gate_reads_the_arc_date_from_the_dispersion_map(log: object, tmp_path: Path) -> None:
    # ARRANGE
    attachedExtractor = _gate_extractor(log, tmp_path / "near", _arc(6.0, 0.5))
    movedExtractor = _gate_extractor(log, tmp_path / "far", _arc(16.0, 0.5))

    # ACT
    attachedExtractor._tune_or_skip(_spectrum(), "VIS")
    movedExtractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert attachedExtractor.calls == []
    assert movedExtractor.calls == ["VIS"]


def test_the_gate_reads_the_science_rotator_angle_from_the_science_frame(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(4.0, 0.5), science=_header(absrot=SCIENCE_ABSROT + 20.0))

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == ["VIS"]


def test_the_gate_measures_an_arc_after_a_nod_sequence_from_the_last_frame(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE EXTRACTED FRAME CARRIES THE FIRST FRAME'S HEADER; THE ARC STARTS 17 MIN AFTER IT, 2 MIN AFTER THE
    # LAST FRAME ENDS
    block = _nod_block()
    arc = _header(minutesFromScience=17.0, absrot=SCIENCE_ABSROT + 7.0, exptime=ARC_EXPTIME)
    extractor = _gate_extractor(log, tmp_path, arc, science=block[0], block=block)

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == []
    printed = [message for level, message in log.messages if level == "print"]
    assert any("Attached arc (Δt 2.0 min, Δrot 1.0°)" in m for m in printed)


def test_the_gate_tunes_when_the_arc_is_far_from_every_frame_of_the_block(log: object, tmp_path: Path) -> None:
    # ARRANGE: THE ARC STARTS 8 MIN AFTER THE LAST FRAME ENDS
    block = _nod_block()
    arc = _header(minutesFromScience=23.0, absrot=SCIENCE_ABSROT + 7.0, exptime=ARC_EXPTIME)
    extractor = _gate_extractor(log, tmp_path, arc, science=block[0], block=block)

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == ["VIS"]


def test_a_dispersion_map_without_a_readable_header_is_tuned(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(4.0, 0.5))
    extractor.dispersionMap = str(tmp_path / "missing.fits")

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == ["VIS"]


def test_an_x_shooter_extraction_is_tuned_even_with_an_attached_looking_arc(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(4.0, 0.5), inst="XSHOOTER")

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == ["VIS"]
    assert qc_by_name(extractor, TUNING_QC)["qc_value"] == "skylines"


def test_the_skip_is_added_to_the_extraction_notes(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _gate_extractor(extraction_note_logger(log), tmp_path, _arc(4.0, 0.5))

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.extractionNotes == [
        ("NOTE", "Attached arc (Δt 4.0 min, Δrot 0.5°): skyline wavelength tuning skipped")
    ]


def test_no_qc_is_recorded_when_the_extractor_collects_none(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _gate_extractor(log, tmp_path, _arc(4.0, 0.5))
    extractor.qc = False

    # ACT
    extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.qc is False
    assert extractor.calls == []
