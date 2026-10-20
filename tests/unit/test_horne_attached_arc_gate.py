"""Contracts for skipping the Horne-extraction skyline tuning when a SOXS arc is attached (DY-1665)."""

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
# SCIENCE FRAME OF THE REAL-DATA CHECK (20260626T044348)
SCIENCE_MJD = 61217.19708617
SCIENCE_ABSROT = -0.57699
SKIPPED_QC = "skipped: attached arc"
TUNING_QC = "SKY TUNING"


def _header(minutesFromScience: float = 0.0, absrot: float | None = SCIENCE_ABSROT, mjd: bool = True) -> fits.Header:
    header = fits.Header()
    if mjd:
        header[KEYWORDS["MJDOBS"]] = SCIENCE_MJD + minutesFromScience / MINUTES_PER_DAY
    if absrot is not None:
        header[KEYWORDS["ABSROT"]] = absrot
        header[KEYWORDS["PAR_ANG_START"]] = -85.867
    return header


def _arc(minutes: float = 3.0, degrees: float = 1.0, **kwargs: object) -> fits.Header:
    return _header(minutesFromScience=minutes, absrot=SCIENCE_ABSROT + degrees, **kwargs)


def _attached(arc: fits.Header, science: fits.Header | None = None) -> bool:
    return is_attached_arc(_header() if science is None else science, arc, KEYWORDS.__getitem__)


def test_an_arc_3_minutes_and_1_degree_away_is_attached() -> None:
    assert _attached(_arc(3.0, 1.0)) is True


def test_the_offsets_are_the_absolute_minutes_and_wrap_aware_degrees() -> None:
    # ARRANGE
    science = _header(absrot=179.0)
    arc = _header(minutesFromScience=-4.0, absrot=-179.0)

    # ACT
    minutes, degrees = attached_arc_offsets(science, arc, KEYWORDS.__getitem__)

    # ASSERT
    assert minutes == pytest.approx(4.0, abs=1e-3)
    assert degrees == pytest.approx(2.0)


def test_an_arc_exactly_at_the_minutes_limit_is_attached() -> None:
    assert _attached(_arc(MAX_ATTACHED_ARC_MINUTES, 1.0)) is True


def test_an_arc_beyond_the_minutes_limit_is_not_attached() -> None:
    assert _attached(_arc(MAX_ATTACHED_ARC_MINUTES + 1.0, 1.0)) is False


def test_an_arc_at_the_rotator_limit_is_not_attached() -> None:
    assert _attached(_arc(3.0, MAX_ATTACHED_ARC_ROTATOR_DEGREES)) is False


def test_an_arc_beyond_the_rotator_limit_is_not_attached() -> None:
    assert _attached(_arc(3.0, 6.0)) is False


def test_rotator_angles_either_side_of_the_wrap_are_attached() -> None:
    assert _attached(_header(absrot=-179.0), science=_header(absrot=179.0)) is True


@pytest.mark.parametrize("minutes", [-3.0, 3.0])
def test_an_arc_before_or_after_the_science_frame_is_attached(minutes: float) -> None:
    assert _attached(_arc(minutes, 1.0)) is True


def test_an_arc_without_a_rotator_angle_is_not_attached() -> None:
    assert _attached(_header(minutesFromScience=3.0, absrot=None)) is False


def test_a_science_frame_without_a_rotator_angle_is_not_attached() -> None:
    assert _attached(_arc(3.0, 1.0), science=_header(absrot=None)) is False


def test_an_arc_with_the_parked_placeholder_rotator_angle_is_not_attached() -> None:
    # ARRANGE
    arc = _header(minutesFromScience=3.0, absrot=0.0)
    arc[KEYWORDS["PAR_ANG_START"]] = 0.0

    # ACT / ASSERT
    assert _attached(arc, science=_header(absrot=1.0)) is False


@pytest.mark.parametrize("missingOn", ["arc", "science"])
def test_a_frame_without_an_observation_date_is_not_attached(missingOn: str) -> None:
    # ARRANGE
    arc = _arc(3.0, 1.0, mjd=missingOn != "arc")
    science = _header(mjd=missingOn != "science")

    # ACT / ASSERT
    assert attached_arc_offsets(science, arc, KEYWORDS.__getitem__) is None
    assert _attached(arc, science=science) is False


# A FITS HEADER CANNOT HOLD NaN OR inf, SO A BOOLEAN AND A STRING ARE THE BAD VALUES REAL HEADERS CAN CARRY
@pytest.mark.parametrize("value", [True, "UNKNOWN"])
@pytest.mark.parametrize("key", ["MJDOBS", "ABSROT"])
@pytest.mark.parametrize("badOn", ["arc", "science"])
def test_a_boolean_or_text_date_or_angle_is_not_attached(badOn: str, key: str, value: object) -> None:
    # ARRANGE
    arc = _arc(3.0, 1.0)
    science = _header()
    (arc if badOn == "arc" else science)[KEYWORDS[key]] = value

    # ACT / ASSERT
    assert attached_arc_offsets(science, arc, KEYWORDS.__getitem__) is None
    assert _attached(arc, science=science) is False


def _gate_extractor(
    log: object, tmp_path: Path, arc: fits.Header, inst: str = "SOXS", science: fits.Header | None = None
) -> object:
    tmp_path.mkdir(parents=True, exist_ok=True)
    scienceHeader = _header() if science is None else science
    extractor = build_extractor(log, tmp_path, dict(scienceHeader), inst=inst)
    extractor.skySubtractedFrame.header = scienceHeader
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
