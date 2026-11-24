"""*tests for the developer script that measures the per-arm sky-shift uncertainty of orders that are not tuned*

The script reads the `SKY SHIFT O<n>` QCs of reduced merged spectra. A shift is the correction that skyline tuning
applied to an order, so it is also the error that remains when the order is not tuned. The shifts of frames whose arc
was far from the science frame give the `untuned` value, and those of frames whose arc was attached give the
`attached-arc` value. Each value is the largest over orders of the RMS about zero, in nm.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest
import yaml
from astropy.io import fits

from tests.unit._tools_loader import load_tool

shift_tool = load_tool("measure_sky_shift_uncertainty")

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]

MINUTE = 1.0 / 1440.0
SCIENCE_MJD = 60000.5
EXPOSURE_S = 600.0


def _science_header(
    arm: str = "VIS",
    shifts: dict[int, float] | None = None,
    mjd: float = SCIENCE_MJD,
    rotator: float = 40.0,
    arcName: str = "arc_SSOL.fits",
    tuning: str | None = "skylines",
) -> fits.Header:
    """*the header of a reduced merged spectrum with sky-shift QCs and the name of its arc*"""
    header = fits.Header()
    header["ESO SEQ ARM"] = arm
    header["MJD-OBS"] = mjd
    header["EXPTIME"] = EXPOSURE_S
    header["ESO ADA ABSROT END"] = rotator
    header["ESO TEL PARANG START"] = 12.0
    header["ESO PRO REC1 CAL1 CATG"] = "DISP_TAB_VIS"
    header["ESO PRO REC1 CAL1 NAME"] = arcName
    if tuning is not None:
        header["ESO QC SKY TUNING"] = tuning
    for order, shift in (shifts or {}).items():
        header[f"ESO QC SKY SHIFT O{order}"] = shift
        header[f"ESO QC SKY SHIFT O{order} PX"] = shift / 0.03
    return header


def _arc_header(minutesAfterScience: float = 60.0, rotator: float | None = 40.0) -> fits.Header:
    """*an arc that starts the given minutes after the science exposure ends*"""
    header = fits.Header()
    header["MJD-OBS"] = SCIENCE_MJD + EXPOSURE_S / 86400.0 + minutesAfterScience * MINUTE
    header["EXPTIME"] = 60.0
    if rotator is not None:
        header["ESO ADA ABSROT END"] = rotator
        header["ESO TEL PARANG START"] = 12.0
    return header


def _sample(kind: str, shifts: dict[int, float], arm: str = "VIS", mjd: float = SCIENCE_MJD) -> object:
    """*a tuned frame's shifts, as the script builds them*"""
    return shift_tool.ShiftSample(arm=arm, kind=kind, mjd=mjd, shifts=shifts)


def test_read_sky_shifts_returns_the_nm_shift_of_each_order_and_skips_the_other_cards() -> None:
    # ARRANGE
    header = _science_header(shifts={1: 0.1, 4: -0.2})
    header["ESO QC SKY SHIFT O2 ERR"] = 0.04

    # ACT
    shifts = shift_tool.read_sky_shifts(header)

    # ASSERT
    assert shifts == {1: 0.1, 4: -0.2}


def test_read_sky_shifts_is_empty_for_a_frame_that_did_not_tune() -> None:
    # ARRANGE
    header = _science_header(tuning="skipped: attached arc")

    # ACT
    shifts = shift_tool.read_sky_shifts(header)

    # ASSERT
    assert shifts == {}


def test_rms_about_zero_is_the_root_mean_square_so_a_common_offset_counts() -> None:
    # ARRANGE / ACT
    spread, offset = shift_tool.rms_about_zero([3.0, 4.0]), shift_tool.rms_about_zero([0.2, 0.2, 0.2])

    # ASSERT
    assert spread == pytest.approx(math.sqrt(12.5))
    assert offset == pytest.approx(0.2)


def test_rms_about_zero_raises_with_no_shift() -> None:
    # ARRANGE / ACT / ASSERT
    with pytest.raises(shift_tool.ToolError, match="no shift"):
        shift_tool.rms_about_zero([])


def test_find_arc_matches_a_header_name_that_the_card_truncated() -> None:
    # ARRANGE
    candidates = [Path("/n1/20260101T0100_VIS_1X1_1_SSOL.fits"), Path("/n1/20260101T0200_VIS_1X1_1_SSOL.fits")]
    header = _science_header(arcName="20260101T0200_VIS_1")

    # ACT
    arc = shift_tool.find_arc(header, candidates)

    # ASSERT
    assert arc == candidates[1]


def test_find_arc_returns_none_when_no_candidate_matches() -> None:
    # ARRANGE
    header = _science_header(arcName="other.fits")

    # ACT
    arc = shift_tool.find_arc(header, [Path("/n1/arc_SSOL.fits")])

    # ASSERT
    assert arc is None


def test_find_arc_raises_when_the_name_matches_two_files() -> None:
    # ARRANGE
    candidates = [Path("/n1/arc_SSOL.fits"), Path("/n2/arc_SSOL.fits")]

    # ACT / ASSERT
    with pytest.raises(shift_tool.ToolError, match="more than one"):
        shift_tool.find_arc(_science_header(), candidates)


@pytest.mark.parametrize(
    ("minutes", "rotator", "kind"),
    [
        (3.0, 42.0, "attached-arc"),
        (7.0, 44.9, "attached-arc"),
        (8.0, 42.0, "untuned"),
        (3.0, 50.0, "untuned"),
        (600.0, 42.0, "untuned"),
    ],
)
def test_classify_frame_splits_on_the_attached_arc_limits(minutes: float, rotator: float, kind: str) -> None:
    # ARRANGE
    science = _science_header(rotator=40.0)
    arc = _arc_header(minutesAfterScience=minutes, rotator=rotator)

    # ACT
    result = shift_tool.classify_frame(science, arc)

    # ASSERT
    assert result == kind


def test_classify_frame_calls_an_arc_with_no_usable_rotator_angle_untuned() -> None:
    # ARRANGE: A DAYTIME ARC CAN NEVER BE ATTACHED (DY-1665)
    arc = _arc_header(minutesAfterScience=3.0, rotator=None)

    # ACT
    kind = shift_tool.classify_frame(_science_header(), arc)

    # ASSERT
    assert kind == "untuned"


def test_uncertainty_by_arm_takes_the_largest_per_order_rms_as_the_value_and_keeps_the_pooled_rms() -> None:
    # ARRANGE
    samples = [
        _sample("untuned", {1: 0.3, 2: 0.4}),
        _sample("untuned", {1: -0.3, 2: -0.4}),
        _sample("attached-arc", {1: 0.01, 4: -0.01}),
    ]

    # ACT
    result = shift_tool.uncertainty_by_arm(samples)

    # ASSERT
    untuned = result["VIS"]["untuned"]
    assert untuned.valueNm == pytest.approx(0.4)
    assert untuned.pooledNm == pytest.approx(math.sqrt((0.09 + 0.16) / 2))
    assert (untuned.shiftCount, untuned.frameCount) == (4, 2)
    assert untuned.perOrderNm == {1: pytest.approx(0.3), 2: pytest.approx(0.4)}
    assert result["VIS"]["attached-arc"].valueNm == pytest.approx(0.01)


def test_uncertainty_by_arm_keeps_the_arms_apart() -> None:
    # ARRANGE
    samples = [_sample("untuned", {1: 0.3}), _sample("untuned", {1: 0.5}, arm="NIR")]

    # ACT
    result = shift_tool.uncertainty_by_arm(samples)

    # ASSERT
    assert result["VIS"]["untuned"].valueNm == pytest.approx(0.3)
    assert result["NIR"]["untuned"].valueNm == pytest.approx(0.5)


def test_format_yaml_prints_the_worst_order_value_and_the_pooled_rms_in_a_comment() -> None:
    # ARRANGE
    samples = [
        _sample("untuned", {1: 0.3, 2: 0.4}, mjd=60000.5),
        _sample("untuned", {1: -0.3, 2: -0.4}, mjd=60001.5),
        _sample("attached-arc", {1: 0.01, 4: -0.01}, mjd=60001.6),
    ]
    result = shift_tool.uncertainty_by_arm(samples)

    # ACT
    text = shift_tool.format_yaml(result, version="0.19.1", today="2026-11-22")

    # ASSERT
    assert text == (
        "VIS:\n"
        "  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 3 FRAMES ON 2 NIGHTS (MJD 60000 TO 60001)\n"
        "  # VALUE: LARGEST PER-ORDER RMS ABOUT ZERO OF THE SKY SHIFTS (SPEC_SYE IS A MAXIMUM OVER ORDERS)\n"
        "  # UNTUNED: 4 SHIFTS FROM 2 FRAMES; POOLED RMS 0.3536; PER ORDER O1 0.3000, O2 0.4000 nm\n"
        "  # ATTACHED-ARC: 2 SHIFTS FROM 1 FRAMES; POOLED RMS 0.0100; PER ORDER O1 0.0100, O4 0.0100 nm\n"
        "  sky-shift-uncertainty-nm:\n"
        "    untuned: 0.4\n"
        "    attached-arc: 0.01\n"
    )


def test_format_yaml_says_when_a_kind_has_no_frame() -> None:
    # ARRANGE
    result = shift_tool.uncertainty_by_arm([_sample("untuned", {1: 0.3})])

    # ACT
    text = shift_tool.format_yaml(result, version="0.19.1", today="2026-11-22")

    # ASSERT
    assert "# ATTACHED-ARC: NO FRAME; NOT MEASURED" in text
    assert yaml.safe_load(text) == {"VIS": {"sky-shift-uncertainty-nm": {"untuned": 0.3}}}


def _write_product(path: Path, header: fits.Header) -> Path:
    """*write a header-only FITS file*"""
    fits.PrimaryHDU(header=header).writeto(path)
    return path


@pytest.fixture
def workspace(tmp_path: Path) -> tuple[Path, Path]:
    """*three tuned frames (one with an attached arc), one frame that skipped tuning, and the two arcs*"""
    products, arcs = tmp_path / "products", tmp_path / "arcs"
    products.mkdir()
    arcs.mkdir()
    _write_product(arcs / "near_SSOL.fits", _arc_header(minutesAfterScience=3.0, rotator=41.0))
    _write_product(arcs / "day_SSOL.fits", _arc_header(minutesAfterScience=600.0, rotator=None))
    day = {"arcName": "day_SSOL.fits"}
    _write_product(products / "a_SPECTRUM.fits", _science_header(shifts={1: 0.3, 4: 0.4}, **day))
    _write_product(products / "b_SPECTRUM.fits", _science_header(shifts={1: -0.3, 4: -0.4}, **day))
    near = {"arcName": "near_SSOL.fits"}
    _write_product(products / "c_SPECTRUM.fits", _science_header(shifts={1: 0.01, 4: -0.01}, **near))
    _write_product(products / "d_SPECTRUM.fits", _science_header(tuning="skipped: attached arc", **near))
    return products, arcs


def test_main_prints_the_untuned_and_attached_arc_values(
    workspace: tuple[Path, Path], capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    products, arcs = workspace
    monkeypatch.setattr(shift_tool, "today_iso", lambda: "2026-11-22")

    # ACT
    status = shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    # ASSERT
    captured = capsys.readouterr()
    values = yaml.safe_load(captured.out)["VIS"]["sky-shift-uncertainty-nm"]
    assert status == 0
    assert values == {"untuned": pytest.approx(0.4), "attached-arc": 0.01}
    assert "d_SPECTRUM.fits" in captured.err


def test_main_raises_when_no_product_holds_a_sky_shift(tmp_path: Path) -> None:
    # ARRANGE
    _write_product(tmp_path / "p.fits", _science_header(tuning="skipped: switched off"))
    _write_product(tmp_path / "arc_SSOL.fits", _arc_header())
    arguments = ["--products", str(tmp_path / "p.fits"), "--arcs", str(tmp_path / "arc_SSOL.fits")]

    # ACT / ASSERT
    with pytest.raises(shift_tool.ToolError, match="no sky shift"):
        shift_tool.main(arguments)


def test_main_skips_a_product_whose_arc_is_not_among_the_given_arcs(
    workspace: tuple[Path, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    products, arcs = workspace
    (arcs / "day_SSOL.fits").unlink()

    # ACT
    shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    # ASSERT
    assert "a_SPECTRUM.fits" in capsys.readouterr().err


def test_main_skips_and_names_a_product_that_cannot_be_read(
    workspace: tuple[Path, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    products, arcs = workspace
    (products / "e_SPECTRUM.fits").write_text("not a FITS file")

    # ACT
    status = shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    # ASSERT
    captured = capsys.readouterr()
    assert status == 0
    assert "skipped" in captured.err
    assert "e_SPECTRUM.fits" in captured.err


def test_main_skips_and_names_a_product_whose_arc_cannot_be_read(
    workspace: tuple[Path, Path], capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    products, arcs = workspace
    (arcs / "day_SSOL.fits").write_text("not a FITS file")

    # ACT
    shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    # ASSERT
    assert "day_SSOL.fits" in capsys.readouterr().err


def test_main_names_the_file_and_card_when_a_tuned_product_has_no_mjd(workspace: tuple[Path, Path]) -> None:
    # ARRANGE
    products, arcs = workspace
    path = products / "a_SPECTRUM.fits"
    with fits.open(path, mode="update") as hdus:
        del hdus[0].header["MJD-OBS"]

    # ACT / ASSERT
    with pytest.raises(shift_tool.ToolError, match=re.escape(f"{path} has no MJD-OBS card")):
        shift_tool.main(["--products", str(path), "--arcs", str(arcs / "day_SSOL.fits")])


def test_the_records_hold_their_shifts_and_mjds_read_only() -> None:
    # ARRANGE
    sample = _sample("untuned", {1: 0.3})
    uncertainty = shift_tool.uncertainty_by_arm([sample])["VIS"]["untuned"]

    # ACT / ASSERT
    with pytest.raises(TypeError):
        sample.shifts[2] = 0.1
    with pytest.raises(TypeError):
        uncertainty.perOrderNm[2] = 0.1
    assert isinstance(uncertainty.mjds, tuple)
