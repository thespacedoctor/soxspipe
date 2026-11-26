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
common = shift_tool.common

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
        header[shift_tool.shift_card(order)] = shift
        header[shift_tool.shift_px_card(order)] = shift / 0.03
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


def _vis_header(
    shifts: dict[int, float],
    rSource: str | None = "measured",
    ugSource: str | None = "model",
    ugComment: str = "predicted from r/i shifts and ABSROT",
    iSnr: float | None = 12.0,
    **kwargs: object,
) -> fits.Header:
    """*a VIS product with the cards that say where each order's shift came from*"""
    header = _science_header(shifts=shifts, **kwargs)
    if rSource is not None:
        header["ESO QC R SHIFT SOURCE"] = rSource
    if ugSource is not None:
        header["ESO QC UG SHIFT SOURCE"] = (ugSource, ugComment)
    if iSnr is not None:
        header[shift_tool.line_snr_card(4)] = iSnr
    return header


VIS_SHIFTS = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}


def test_unmeasured_orders_names_the_r_order_when_its_source_is_none() -> None:
    # ARRANGE
    header = _vis_header(VIS_SHIFTS, rSource="none", ugSource="none", ugComment="order 1: 0 sky lines (min 2)")

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert 1 in leftOut
    assert "R SHIFT SOURCE" in leftOut[1]


def test_unmeasured_orders_names_the_u_and_g_orders_when_their_source_is_none() -> None:
    # ARRANGE
    header = _vis_header({1: 0.02, 2: 0.0, 3: 0.0, 4: 0.03}, ugSource="none", ugComment="order 4: 6 sky lines (min 10)")

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [2, 3]
    assert "UG SHIFT SOURCE" in leftOut[2]


def test_unmeasured_orders_names_the_i_order_when_its_sky_line_snr_is_below_the_minimum() -> None:
    # ARRANGE
    header = _vis_header(
        {1: 0.01, 2: 0.0, 3: 0.0, 4: 0.0}, ugSource="none", ugComment="order 4: sky-line SNR 1.5 (min 5.4)", iSnr=1.5
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [2, 3, 4]
    assert "SNR" in leftOut[4]


def test_unmeasured_orders_names_the_i_order_when_it_has_no_sky_line_snr_card() -> None:
    # ARRANGE
    header = _vis_header(
        {1: 0.01, 2: 0.0, 3: 0.0, 4: 0.0}, ugSource="none", ugComment="order 4: no sky shift", iSnr=None
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert 4 in leftOut


def test_unmeasured_orders_names_the_i_order_when_the_pipeline_kept_too_few_sky_lines_for_it() -> None:
    # ARRANGE: THE SNR IS HIGH, BUT 0 LINES WERE KEPT, SO THE PIPELINE APPLIED NO SHIFT (THE EG274 FRAME OF DY-2378)
    header = _vis_header(
        {1: -0.002, 2: 0.0, 3: 0.0, 4: 0.0}, ugSource="none", ugComment="order 4: 0 sky lines (min 10)", iSnr=8.6
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [2, 3, 4]
    assert "sky lines" in leftOut[4]


def test_unmeasured_orders_names_an_i_order_whose_rounded_snr_card_sits_at_the_limit_and_whose_shift_is_zero() -> None:
    # ARRANGE: THE CARD HOLDS ROUND(SNR, 1), SO 5.4 MAY COME FROM A TRUE 5.36 THAT THE PIPELINE REFUSED (IT TESTS < 5.4)
    header = _vis_header(
        {1: 0.01, 2: 0.0, 3: 0.0, 4: 0.0}, ugSource="none", ugComment="order 4: 12 sky lines (min 10)", iSnr=5.4
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert 4 in leftOut


def test_unmeasured_orders_counts_an_i_order_at_the_snr_limit_when_its_shift_is_not_zero() -> None:
    # ARRANGE
    header = _vis_header(
        {1: 0.01, 2: 0.0, 3: 0.0, 4: 0.012}, ugSource="none", ugComment="order 4: 12 sky lines (min 10)", iSnr=5.4
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [2, 3]


def test_unmeasured_orders_counts_a_non_zero_shift_whatever_the_source_cards_say() -> None:
    # ARRANGE: A NON-ZERO SHIFT IS ALWAYS A MEASUREMENT, EVEN IF A SOURCE CARD SAYS none
    header = _vis_header(
        {1: 0.03, 2: 0.02, 3: -0.02, 4: 0.04},
        rSource="none",
        ugSource="none",
        ugComment="order 1: 0 sky lines",
        iSnr=0.1,
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert leftOut == {}


def test_the_line_count_pattern_parses_the_comment_the_pipeline_builds_for_a_failed_i_order() -> None:
    # ARRANGE: THE SCRIPT READS THE WORDING OF _order_guard_reason, SO A CHANGE THERE MUST FAIL HERE
    common.ensure_repo_on_path()
    from soxspipe.commonutils.ug_wavelength_transfer import VIS_I_ORDER, _order_guard_reason

    # ACT
    comment = _order_guard_reason(VIS_I_ORDER, {"nLines": 1, "lineSnr": 20.0})
    match = shift_tool.SKY_LINES_REASON_PATTERN.match(comment)

    # ASSERT
    assert match is not None
    assert (int(match.group(1)), int(match.group(2))) == (VIS_I_ORDER, 1)


def test_unmeasured_orders_keeps_an_i_order_with_a_few_lines_that_the_pipeline_applied() -> None:
    # ARRANGE: 7 LINES ARE BELOW THE 10 THE U/G MODEL NEEDS, BUT ABOVE THE 2 THAT APPLY THE i SHIFT ITSELF
    header = _vis_header(
        {1: -0.049, 2: 0.0, 3: 0.0, 4: -0.043},
        ugSource="none",
        ugComment="order 4: 7 sky lines (min 10)",
        iSnr=7.7,
    )

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [2, 3]


def test_unmeasured_orders_counts_a_measured_shift_of_exactly_zero() -> None:
    # ARRANGE: r MEASURED, u/g MODELLED AND i MEASURED, EACH AT EXACTLY 0.000 nm
    header = _vis_header(VIS_SHIFTS, rSource="measured", ugSource="model", iSnr=40.0)

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert leftOut == {}


def test_unmeasured_orders_counts_a_predicted_and_a_combined_r_shift() -> None:
    # ARRANGE
    predicted = _vis_header({1: 0.05, 2: 0.0, 3: 0.0, 4: 0.04}, rSource="predicted")
    combined = _vis_header({1: 0.05, 2: 0.0, 3: 0.0, 4: 0.04}, rSource="combined")

    # ACT
    leftOutPredicted = shift_tool.unmeasured_orders(predicted, "VIS", shift_tool.read_sky_shifts(predicted))
    leftOutCombined = shift_tool.unmeasured_orders(combined, "VIS", shift_tool.read_sky_shifts(combined))

    # ASSERT
    assert leftOutPredicted == {}
    assert leftOutCombined == {}


def test_unmeasured_orders_leaves_out_nothing_for_a_vis_product_that_records_no_source() -> None:
    # ARRANGE: A BUILD BEFORE THE SOURCE CARDS GIVES NO WAY TO TELL, SO EVERY SHIFT COUNTS
    header = _vis_header(VIS_SHIFTS, rSource=None, ugSource=None, iSnr=None)

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "VIS", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert leftOut == {}


def test_unmeasured_orders_names_a_nir_order_with_a_zero_shift_and_no_px_card() -> None:
    # ARRANGE: THE NIR GLOBAL FIT RECORDS 0 AND NO PX CARD FOR AN ORDER WITHOUT A CORRECTION
    header = _science_header(arm="NIR", shifts={11: 0.02})
    header[shift_tool.shift_card(12)] = 0.0

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "NIR", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert sorted(leftOut) == [12]
    assert "PX" in leftOut[12]


def test_unmeasured_orders_counts_a_nir_zero_shift_that_has_a_px_card() -> None:
    # ARRANGE: A SHIFT UNDER 0.0005 nm ROUNDS TO 0.0, WITH THE PX CARD STILL HOLDING THE MEASUREMENT
    header = _science_header(arm="NIR", shifts={11: 0.02, 12: 0.0})
    header[shift_tool.shift_px_card(12)] = 0.008

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "NIR", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert leftOut == {}


def test_unmeasured_orders_counts_a_nir_shift_that_has_no_px_card_but_is_not_zero() -> None:
    # ARRANGE
    header = _science_header(arm="NIR", shifts={11: 0.02})
    del header[shift_tool.shift_px_card(11)]

    # ACT
    leftOut = shift_tool.unmeasured_orders(header, "NIR", shift_tool.read_sky_shifts(header))

    # ASSERT
    assert leftOut == {}


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


def _attached_workspace(tmp_path: Path, frames: dict[str, fits.Header]) -> list[str]:
    """*the given products, each with the near (attached) arc, and the command line that reads them*"""
    products, arcs = tmp_path / "products", tmp_path / "arcs"
    products.mkdir()
    arcs.mkdir()
    _write_product(arcs / "near_SSOL.fits", _arc_header(minutesAfterScience=3.0, rotator=41.0))
    for name, header in frames.items():
        _write_product(products / name, header)
    return ["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")]


ATTACHED_ARC_KWARGS = {"arcName": "near_SSOL.fits"}


def test_main_leaves_unmeasured_orders_out_of_the_rms_and_names_each_one(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE: THE FIRST FRAME HAS u, g AND i UNMEASURED (EXACTLY 0), THE SECOND MEASURED 0.01 nm IN EVERY ORDER
    unmeasured = _vis_header(
        {1: -0.002, 2: 0.0, 3: 0.0, 4: 0.0},
        ugSource="none",
        ugComment="order 4: 0 sky lines (min 10)",
        iSnr=8.6,
        **ATTACHED_ARC_KWARGS,
    )
    measured = _vis_header({1: 0.01, 2: 0.01, 3: 0.01, 4: 0.01}, **ATTACHED_ARC_KWARGS)
    arguments = _attached_workspace(tmp_path, {"a_SPECTRUM.fits": unmeasured, "b_SPECTRUM.fits": measured})

    # ACT
    shift_tool.main(arguments)

    # ASSERT
    captured = capsys.readouterr()
    assert yaml.safe_load(captured.out)["VIS"]["sky-shift-uncertainty-nm"]["attached-arc"] == 0.01
    assert "# ATTACHED-ARC: 5 SHIFTS FROM 2 FRAMES" in captured.out
    for order in (2, 3, 4):
        assert re.search(rf"left out order {order} .*a_SPECTRUM\.fits", captured.err)
    assert "b_SPECTRUM.fits" not in captured.err


def test_main_counts_a_measured_shift_of_exactly_zero(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    # ARRANGE: THE SECOND FRAME MEASURED OR MODELLED EVERY ORDER AT EXACTLY 0.000 nm
    shifted = _vis_header({1: 0.01, 2: 0.01, 3: 0.01, 4: 0.01}, **ATTACHED_ARC_KWARGS)
    zero = _vis_header(VIS_SHIFTS, iSnr=40.0, **ATTACHED_ARC_KWARGS)
    arguments = _attached_workspace(tmp_path, {"a_SPECTRUM.fits": shifted, "b_SPECTRUM.fits": zero})

    # ACT
    shift_tool.main(arguments)

    # ASSERT
    captured = capsys.readouterr()
    assert yaml.safe_load(captured.out)["VIS"]["sky-shift-uncertainty-nm"]["attached-arc"] == round(0.01 / 2**0.5, 4)
    assert "left out" not in captured.err


def test_main_skips_a_frame_with_no_measured_order_and_names_its_orders(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE: THE SECOND FRAME IS A FAINT-SKY FRAME WHERE r, u, g AND i ARE ALL UNMEASURED (EXACT ZEROS)
    measured = _vis_header({1: 0.01, 2: 0.01, 3: 0.01, 4: 0.01}, **ATTACHED_ARC_KWARGS)
    faint = _vis_header(
        VIS_SHIFTS,
        rSource="none",
        ugSource="none",
        ugComment="order 1: 0 sky lines (min 2)",
        iSnr=0.1,
        **ATTACHED_ARC_KWARGS,
    )
    arguments = _attached_workspace(tmp_path, {"a_SPECTRUM.fits": measured, "b_SPECTRUM.fits": faint})

    # ACT
    shift_tool.main(arguments)

    # ASSERT
    captured = capsys.readouterr()
    assert "# ATTACHED-ARC: 4 SHIFTS FROM 1 FRAMES" in captured.out
    assert re.search(r"skipped .*b_SPECTRUM\.fits: no sky shift was measured", captured.err)
    assert len(re.findall(r"left out order \d .*b_SPECTRUM\.fits", captured.err)) == 4


def test_main_raises_when_every_order_of_every_product_is_unmeasured(tmp_path: Path) -> None:
    # ARRANGE
    faint = _vis_header(
        VIS_SHIFTS,
        rSource="none",
        ugSource="none",
        ugComment="order 1: 0 sky lines (min 2)",
        iSnr=0.1,
        **ATTACHED_ARC_KWARGS,
    )
    arguments = _attached_workspace(tmp_path, {"a_SPECTRUM.fits": faint})

    # ACT / ASSERT
    with pytest.raises(shift_tool.ToolError, match="no sky shift"):
        shift_tool.main(arguments)


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
