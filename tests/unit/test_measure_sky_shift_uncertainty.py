"""*tests for the developer script that measures the per-arm sky-shift uncertainty of orders that are not tuned*

The script reads the `SKY SHIFT O<n>` QCs of reduced spectra. A shift is the correction that skyline tuning applied to
an order, so it is also the error that remains when the order is not tuned. The shifts of frames whose arc was far from
the science frame give the `untuned` value, and those of frames whose arc was attached give the `attached-arc` value.
Both are the RMS about zero of the shifts, in nm.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import pytest
import yaml
from astropy.io import fits

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import measure_sky_shift_uncertainty as shift_tool  # noqa: E402
import measurement_common as common  # noqa: E402

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]

MINUTE = 1.0 / 1440.0
SCIENCE_MJD = 60000.5
EXPOSURE_S = 600.0


def _science_header(arm="VIS", shifts=None, mjd=SCIENCE_MJD, rotator=40.0, arcName="arc_SSOL.fits", tuning="skylines"):
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


def _arc_header(minutesAfterScience=60.0, rotator=40.0):
    """*an arc that starts the given minutes after the science exposure ends*"""
    header = fits.Header()
    header["MJD-OBS"] = SCIENCE_MJD + EXPOSURE_S / 86400.0 + minutesAfterScience * MINUTE
    header["EXPTIME"] = 60.0
    if rotator is not None:
        header["ESO ADA ABSROT END"] = rotator
        header["ESO TEL PARANG START"] = 12.0
    return header


def test_read_sky_shifts_returns_the_nm_shift_of_each_order_and_skips_the_other_cards():
    header = _science_header(shifts={1: 0.1, 4: -0.2})
    header["ESO QC SKY SHIFT O2 ERR"] = 0.04

    assert shift_tool.read_sky_shifts(header) == {1: 0.1, 4: -0.2}


def test_read_sky_shifts_is_empty_for_a_frame_that_did_not_tune():
    assert shift_tool.read_sky_shifts(_science_header(tuning="skipped: attached arc")) == {}


def test_rms_about_zero_is_the_root_mean_square_so_a_common_offset_counts():
    assert shift_tool.rms_about_zero([3.0, 4.0]) == pytest.approx(math.sqrt(12.5))
    assert shift_tool.rms_about_zero([0.2, 0.2, 0.2]) == pytest.approx(0.2)


def test_rms_about_zero_raises_with_no_shift():
    with pytest.raises(shift_tool.ToolError, match="no shift"):
        shift_tool.rms_about_zero([])


def test_find_arc_matches_a_header_name_that_the_card_truncated():
    candidates = [Path("/n1/20260101T0100_VIS_1X1_1_SSOL.fits"), Path("/n1/20260101T0200_VIS_1X1_1_SSOL.fits")]
    header = _science_header(arcName="20260101T0200_VIS_1")

    assert shift_tool.find_arc(header, candidates) == candidates[1]


def test_find_arc_returns_none_when_no_candidate_matches():
    assert shift_tool.find_arc(_science_header(arcName="other.fits"), [Path("/n1/arc_SSOL.fits")]) is None


def test_find_arc_raises_when_the_name_matches_two_files():
    candidates = [Path("/n1/arc_SSOL.fits"), Path("/n2/arc_SSOL.fits")]

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
def test_classify_frame_splits_on_the_attached_arc_limits(minutes, rotator, kind):
    science = _science_header(rotator=40.0)
    arc = _arc_header(minutesAfterScience=minutes, rotator=rotator)

    assert shift_tool.classify_frame(science, arc) == kind


def test_classify_frame_calls_an_arc_with_no_usable_rotator_angle_untuned():
    # A DAYTIME ARC CAN NEVER BE ATTACHED (DY-1665)
    assert shift_tool.classify_frame(_science_header(), _arc_header(minutesAfterScience=3.0, rotator=None)) == "untuned"


def _sample(kind, shifts, arm="VIS", mjd=SCIENCE_MJD):
    return shift_tool.ShiftSample(arm=arm, kind=kind, mjd=mjd, shifts=shifts)


def test_uncertainty_by_arm_pools_the_shifts_of_all_orders_and_frames_of_a_kind():
    samples = [
        _sample("untuned", {1: 0.3, 2: 0.4}),
        _sample("untuned", {1: -0.3, 2: -0.4}),
        _sample("attached-arc", {1: 0.01, 4: -0.01}),
    ]

    result = shift_tool.uncertainty_by_arm(samples, aggregate="pooled")

    untuned = result["VIS"]["untuned"]
    assert untuned.valueNm == pytest.approx(math.sqrt((0.09 + 0.16) / 2))
    assert (untuned.shiftCount, untuned.frameCount) == (4, 2)
    assert untuned.perOrderNm == {1: pytest.approx(0.3), 2: pytest.approx(0.4)}
    assert result["VIS"]["attached-arc"].valueNm == pytest.approx(0.01)


def test_uncertainty_by_arm_worst_order_takes_the_largest_per_order_rms():
    samples = [_sample("untuned", {1: 0.3, 2: 0.4}), _sample("untuned", {1: -0.3, 2: -0.4})]

    result = shift_tool.uncertainty_by_arm(samples, aggregate="worst-order")

    assert result["VIS"]["untuned"].valueNm == pytest.approx(0.4)


def test_uncertainty_by_arm_keeps_the_arms_apart():
    samples = [_sample("untuned", {1: 0.3}), _sample("untuned", {1: 0.5}, arm="NIR")]

    result = shift_tool.uncertainty_by_arm(samples, aggregate="pooled")

    assert result["VIS"]["untuned"].valueNm == pytest.approx(0.3)
    assert result["NIR"]["untuned"].valueNm == pytest.approx(0.5)


def test_uncertainty_by_arm_rejects_an_unknown_aggregate():
    with pytest.raises(shift_tool.ToolError, match="aggregate"):
        shift_tool.uncertainty_by_arm([_sample("untuned", {1: 0.3})], aggregate="mean")


def test_format_yaml_prints_both_values_with_the_provenance():
    samples = [
        _sample("untuned", {1: 0.3, 2: 0.4}, mjd=60000.5),
        _sample("untuned", {1: -0.3, 2: -0.4}, mjd=60001.5),
        _sample("attached-arc", {1: 0.01, 4: -0.01}, mjd=60001.6),
    ]
    result = shift_tool.uncertainty_by_arm(samples, aggregate="pooled")

    text = shift_tool.format_yaml(result, "pooled", version="0.19.1", today="2026-11-22")

    assert text == (
        "VIS:\n"
        "  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 3 FRAMES ON 2 NIGHTS (MJD 60000 TO 60001)\n"
        "  # RMS ABOUT ZERO OF THE SKY SHIFTS, POOLED OVER ORDERS AND FRAMES\n"
        "  # UNTUNED: 4 SHIFTS FROM 2 FRAMES; PER ORDER O1 0.3000, O2 0.4000 nm\n"
        "  # ATTACHED-ARC: 2 SHIFTS FROM 1 FRAMES; PER ORDER O1 0.0100, O4 0.0100 nm\n"
        "  sky-shift-uncertainty-nm:\n"
        "    untuned: 0.3536\n"
        "    attached-arc: 0.01\n"
    )


def test_format_yaml_says_when_a_kind_has_no_frame():
    result = shift_tool.uncertainty_by_arm([_sample("untuned", {1: 0.3})], aggregate="pooled")

    text = shift_tool.format_yaml(result, "pooled", version="0.19.1", today="2026-11-22")

    assert "# ATTACHED-ARC: NO FRAME; NOT MEASURED" in text
    assert yaml.safe_load(text) == {"VIS": {"sky-shift-uncertainty-nm": {"untuned": 0.3}}}


def _write_product(path: Path, header: fits.Header) -> Path:
    fits.PrimaryHDU(header=header).writeto(path)
    return path


def _write_arc(path: Path, **kwargs) -> Path:
    return _write_product(path, _arc_header(**kwargs))


@pytest.fixture
def workspace(tmp_path):
    """*three tuned frames (one with an attached arc), one frame that skipped tuning and one that is not VIS-tuned*"""
    products, arcs = tmp_path / "products", tmp_path / "arcs"
    products.mkdir()
    arcs.mkdir()
    _write_arc(arcs / "near_SSOL.fits", minutesAfterScience=3.0, rotator=41.0)
    _write_arc(arcs / "day_SSOL.fits", minutesAfterScience=600.0, rotator=None)
    _write_product(
        products / "a_EXTRACTED_MERGED.fits", _science_header(shifts={1: 0.3, 4: 0.4}, arcName="day_SSOL.fits")
    )
    _write_product(
        products / "b_EXTRACTED_MERGED.fits", _science_header(shifts={1: -0.3, 4: -0.4}, arcName="day_SSOL.fits")
    )
    _write_product(
        products / "c_EXTRACTED_MERGED.fits", _science_header(shifts={1: 0.01, 4: -0.01}, arcName="near_SSOL.fits")
    )
    _write_product(
        products / "d_EXTRACTED_MERGED.fits",
        _science_header(arcName="near_SSOL.fits", tuning="skipped: attached arc"),
    )
    return products, arcs


def test_main_prints_the_untuned_and_attached_arc_values(workspace, capsys, monkeypatch):
    products, arcs = workspace
    monkeypatch.setattr(shift_tool, "today_iso", lambda: "2026-11-22")

    status = shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    captured = capsys.readouterr()
    values = yaml.safe_load(captured.out)["VIS"]["sky-shift-uncertainty-nm"]
    assert status == 0
    assert values == {"untuned": pytest.approx(math.sqrt((0.09 + 0.16) / 2), abs=1e-4), "attached-arc": 0.01}
    assert "d_EXTRACTED_MERGED.fits" in captured.err


def test_main_raises_when_no_product_holds_a_sky_shift(tmp_path):
    _write_product(tmp_path / "p.fits", _science_header(tuning="skipped: switched off"))
    _write_arc(tmp_path / "arc_SSOL.fits")

    with pytest.raises(shift_tool.ToolError, match="no sky shift"):
        shift_tool.main(["--products", str(tmp_path / "p.fits"), "--arcs", str(tmp_path / "arc_SSOL.fits")])


def test_main_skips_a_product_whose_arc_is_not_among_the_given_arcs(workspace, capsys):
    products, arcs = workspace
    (arcs / "day_SSOL.fits").unlink()

    shift_tool.main(["--products", str(products / "*.fits"), "--arcs", str(arcs / "*.fits")])

    assert "a_EXTRACTED_MERGED.fits" in capsys.readouterr().err


def test_the_keyword_aliases_used_for_the_arc_test_exist_in_the_soxs_keywords():
    aliases = common.load_yaml(common.KEYWORDS_PATH)

    assert {"MJDOBS", "EXPTIME", "SEQ_EXPTIME", "ABSROT", "PAR_ANG_START"} <= set(aliases)
