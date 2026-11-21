"""*tests for the developer script that measures the nominal SPEC_RES of each arm and slit*

The line measurement itself is tested in `test_line_fitting.py`. These tests cover what the script adds around it: the
choice of frames, the pairing with a DISP_IMAGE, the trim, the table and its yaml, and the sky-line cross-check.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.io import fits

from tests.unit._arc_scene import (
    Scene,
    raw_header,
    sky_header,
    write_disp_image,
    write_raw_arc,
    write_sky_pair,
)
from tests.unit._tools_loader import load_tool

spec_tool = load_tool("measure_spec_res")
frames_tool = load_tool("spec_res_frames")
line_fitting = load_tool("line_fitting")

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]


@pytest.fixture
def full_frame_trim(monkeypatch: pytest.MonkeyPatch) -> None:
    """*the synthetic frames are not the size of a SOXS detector, so keep every pixel*"""
    monkeypatch.setattr(frames_tool, "science_slices", lambda arm, binX, binY: (slice(None), slice(None)))


def _lines(*powers: float) -> list:
    """*fitted lines with the given resolving powers*"""
    return [line_fitting.LineFit(1, 500.0, 0.1, 3.0, power) for power in powers]


def _measurements() -> list:
    """*two SLIT1.0 frames on two nights and one SLIT0.5 frame*"""
    return [
        spec_tool.FrameMeasurement("VIS", "SLIT1.0", 60000.4, _lines(4800, 5000, 5200)),
        spec_tool.FrameMeasurement("VIS", "SLIT1.0", 60001.4, _lines(4900, 5100, 5300, 5500)),
        spec_tool.FrameMeasurement("VIS", "SLIT0.5", 60000.5, _lines(9000, 9400)),
    ]


def test_derive_wide_slit_scales_the_one_arcsec_value_by_the_inverse_width() -> None:
    # ARRANGE
    values = {"SLIT0.5": 8965, "SLIT1.0": 4906, "SLIT1.5": 3271}

    # ACT
    table = spec_tool.derive_wide_slit(values)

    # ASSERT
    assert table["SLIT5.0"] == 981
    assert table["SLIT1.0"] == 4906


def test_derive_wide_slit_adds_nothing_when_the_one_arcsec_slit_was_not_measured() -> None:
    # ARRANGE
    values = {"SLIT0.5": 8965}

    # ACT
    table = spec_tool.derive_wide_slit(values)

    # ASSERT
    assert "SLIT5.0" not in table


def test_slit_width_raises_for_a_name_that_is_not_a_slit() -> None:
    # ARRANGE / ACT / ASSERT
    with pytest.raises(spec_tool.ToolError, match="SLIT<width>"):
        spec_tool.slit_width("PINHOLE")


def test_summarise_takes_the_median_over_all_lines_of_a_slit_and_counts_them() -> None:
    # ARRANGE
    measurements = _measurements()

    # ACT
    slit = spec_tool.summarise(measurements)["VIS"]["SLIT1.0"]

    # ASSERT
    assert slit.resolvingPower == 5100
    assert slit.lineCount == 7
    assert slit.frameCount == 2
    assert slit.nightMedians == {60000: 5000, 60001: 5200}


def test_summarise_leaves_out_a_frame_that_found_no_lines() -> None:
    # ARRANGE
    measurements = [*_measurements(), spec_tool.FrameMeasurement("VIS", "SLIT1.5", 60002.4, [])]

    # ACT
    summary = spec_tool.summarise(measurements)

    # ASSERT
    assert "SLIT1.5" not in summary["VIS"]


def test_summarise_is_empty_when_no_frame_found_a_line() -> None:
    # ARRANGE
    measurements = [spec_tool.FrameMeasurement("VIS", "SLIT1.0", 60000.4, [])]

    # ACT
    summary = spec_tool.summarise(measurements)

    # ASSERT
    assert summary == {}


def test_format_yaml_prints_the_table_in_slit_order_with_the_derived_wide_slit() -> None:
    # ARRANGE
    summary = spec_tool.summarise(_measurements())

    # ACT
    text = spec_tool.format_yaml(summary, version="0.19.1", today="2026-11-22")

    # ASSERT
    assert text == (
        "VIS:\n"
        "  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 3 SLIT ARCS ON 2 NIGHTS (MJD 60000 TO 60001)\n"
        "  # SLIT0.5: 2 LINES FROM 1 FRAMES; MEDIAN R PER NIGHT 9200 (MJD 60000)\n"
        "  # SLIT1.0: 7 LINES FROM 2 FRAMES; MEDIAN R PER NIGHT 5000 (MJD 60000), 5200 (MJD 60001)\n"
        "  # SLIT5.0 IS DERIVED, NOT MEASURED: THE SLIT1.0 VALUE TIMES 1.0 / 5.0\n"
        "  spec-res:\n"
        "    SLIT0.5: 9200\n"
        "    SLIT1.0: 5100\n"
        "    SLIT5.0: 1020\n"
    )


def test_format_yaml_block_loads_back_as_integer_resolving_powers() -> None:
    # ARRANGE
    summary = spec_tool.summarise(_measurements())

    # ACT
    block = yaml.safe_load(spec_tool.format_yaml(summary, version="0.19.1", today="2026-11-22"))

    # ASSERT
    assert block == {"VIS": {"spec-res": {"SLIT0.5": 9200, "SLIT1.0": 5100, "SLIT5.0": 1020}}}


def test_format_yaml_says_when_the_wide_slit_cannot_be_derived() -> None:
    # ARRANGE
    measurements = [spec_tool.FrameMeasurement("VIS", "SLIT0.5", 60000.5, _lines(9000, 9400))]
    summary = spec_tool.summarise(measurements)

    # ACT
    text = spec_tool.format_yaml(summary, version="0.19.1", today="2026-11-22")

    # ASSERT
    assert "# SLIT5.0 NOT DERIVED: NO SLIT1.0 MEASUREMENT" in text
    assert "SLIT5.0:" not in text.replace("# SLIT5.0", "")


def test_format_yaml_prints_the_sky_cross_check_as_comments_and_leaves_the_table_alone() -> None:
    # ARRANGE
    summary = spec_tool.summarise(_measurements())
    sky = spec_tool.summarise([spec_tool.FrameMeasurement("VIS", "SLIT1.0", 60000.45, _lines(5355, 5355))])

    # ACT
    withSky = spec_tool.format_yaml(summary, "0.19.1", "2026-11-22", skySummary=sky)
    without = spec_tool.format_yaml(summary, "0.19.1", "2026-11-22")

    # ASSERT
    assert (
        "  # SKY-LINE CROSS-CHECK, NOT IN THE TABLE: SLIT1.0 R 5355 (2 LINES FROM 1 FRAMES, SKY/ARC 1.05)\n" in withSky
    )
    assert yaml.safe_load(withSky) == yaml.safe_load(without)


def test_format_yaml_prints_a_sky_slit_with_no_arc_value_without_a_ratio() -> None:
    # ARRANGE
    summary = spec_tool.summarise(_measurements())
    sky = spec_tool.summarise([spec_tool.FrameMeasurement("VIS", "SLIT1.5", 60000.45, _lines(3500))])

    # ACT
    text = spec_tool.format_yaml(summary, "0.19.1", "2026-11-22", skySummary=sky)

    # ASSERT
    assert "NOT IN THE TABLE: SLIT1.5 R 3500 (1 LINES FROM 1 FRAMES)\n" in text


def _write_night(tmp_path: Path) -> tuple[Path, Path]:
    """*one night: SLIT1.0 (3 px) and SLIT1.5 (4.5 px) arcs, a SLIT5.0 arc, a pinhole arc and one DISP_IMAGE*"""
    rawDir, mapDir = tmp_path / "raw", tmp_path / "maps"
    rawDir.mkdir()
    mapDir.mkdir()
    write_raw_arc(rawDir / "slit10.fits", Scene(fwhmPx=3.0, seed=1), raw_header(slit="SLIT1.0", mjd=60000.40))
    write_raw_arc(rawDir / "slit15.fits", Scene(fwhmPx=4.5, seed=2), raw_header(slit="SLIT1.5", mjd=60000.41))
    write_raw_arc(rawDir / "slit50.fits", Scene(fwhmPx=6.0, seed=3), raw_header(slit="SLIT5.0", mjd=60000.42))
    write_raw_arc(rawDir / "pinhole.fits", Scene(fwhmPx=2.0, seed=4), raw_header(tech="ECHELLE,PINHOLE", mjd=60000.43))
    write_disp_image(mapDir / "night_DISP_IMAGE.fits", Scene(), mjd=60000.1)
    return rawDir, mapDir


def test_main_prints_the_measured_table_with_the_derived_wide_slit(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    scene = Scene()
    expected10 = scene.expected_resolving_power(np.median(scene.centres))

    # ACT
    status = spec_tool.main(["--arcs", str(rawDir / "*.fits"), "--disp-images", str(mapDir / "*_DISP_IMAGE.fits")])

    # ASSERT
    captured = capsys.readouterr()
    table = yaml.safe_load(captured.out)["VIS"]["spec-res"]
    assert status == 0
    assert list(table) == ["SLIT1.0", "SLIT1.5", "SLIT5.0"]
    assert table["SLIT1.0"] == pytest.approx(expected10, rel=0.03)
    assert table["SLIT1.5"] == pytest.approx(expected10 / 1.5, rel=0.04)
    assert table["SLIT5.0"] == round(table["SLIT1.0"] / 5)
    assert "slit50.fits" in captured.err
    assert "pinhole.fits" in captured.err


def test_main_adds_the_sky_cross_check_as_comments_and_the_table_does_not_change(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE: THE SKY LINES ARE 20% WIDER THAN THE ARC LINES
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    stareDir = tmp_path / "stare"
    stareDir.mkdir()
    skyScene = Scene(fwhmPx=3.6, minSlit=1.9, seed=5)
    write_sky_pair(stareDir / "obs1_SKYSUB.fits", skyScene, sky_header(slit="SLIT1.0"))
    arcArguments = ["--arcs", str(rawDir / "*.fits"), "--disp-images", str(mapDir / "*_DISP_IMAGE.fits")]
    spec_tool.main(arcArguments)
    tableWithoutSky = yaml.safe_load(capsys.readouterr().out)

    # ACT
    spec_tool.main([*arcArguments, "--sky-products", str(stareDir / "*_SKYSUB.fits")])

    # ASSERT
    output = capsys.readouterr().out
    assert yaml.safe_load(output) == tableWithoutSky
    match = re.search(
        r"SKY-LINE CROSS-CHECK, NOT IN THE TABLE: SLIT1.0 R (\d+) \(\d+ LINES FROM 1 FRAMES, SKY/ARC ([\d.]+)\)", output
    )
    assert match is not None
    assert float(match.group(2)) == pytest.approx(1 / 1.2, rel=0.06)


def test_main_names_a_sky_frame_whose_model_is_missing(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    stareDir = tmp_path / "stare"
    stareDir.mkdir()
    write_sky_pair(stareDir / "obs1_SKYSUB.fits", Scene(minSlit=1.9), sky_header())
    (stareDir / "obs1_SKYMODEL.fits").unlink()

    # ACT
    spec_tool.main(
        [
            "--arcs",
            str(rawDir / "*.fits"),
            "--disp-images",
            str(mapDir / "*_DISP_IMAGE.fits"),
            "--sky-products",
            str(stareDir / "*_SKYSUB.fits"),
        ]
    )

    # ASSERT
    captured = capsys.readouterr()
    assert "obs1_SKYSUB.fits" in captured.err
    assert "SKY-LINE CROSS-CHECK" not in captured.out


def test_main_names_the_arcs_that_have_no_disp_image_in_range_and_stops_when_none_is_paired(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], full_frame_trim: None
) -> None:
    # ARRANGE
    rawDir, mapDir = _write_night(tmp_path)
    arguments = [
        "--arcs",
        str(rawDir / "slit10.fits"),
        "--disp-images",
        str(mapDir / "*.fits"),
        "--max-gap-days",
        "0.1",
    ]

    # ACT / ASSERT
    with pytest.raises(spec_tool.ToolError, match="paired"):
        spec_tool.main(arguments)
    assert "slit10.fits" in capsys.readouterr().err


def test_main_raises_when_no_arc_line_is_measured(tmp_path: Path, full_frame_trim: None) -> None:
    # ARRANGE
    rawDir, mapDir = _write_night(tmp_path)
    arguments = ["--arcs", str(rawDir / "slit10.fits"), "--disp-images", str(mapDir / "*.fits"), "--saturation", "3000"]

    # ACT / ASSERT
    with pytest.raises(spec_tool.ToolError, match="no arc line"):
        spec_tool.main(arguments)


def test_the_records_hold_their_lines_mjds_and_night_medians_read_only() -> None:
    # ARRANGE
    measurement = spec_tool.FrameMeasurement("VIS", "SLIT1.0", 60000.4, _lines(5000))
    slit = spec_tool.summarise([measurement])["VIS"]["SLIT1.0"]

    # ACT / ASSERT
    assert isinstance(measurement.lines, tuple)
    assert isinstance(slit.mjds, tuple)
    with pytest.raises(TypeError):
        slit.nightMedians[60001] = 1


def test_main_prints_the_arc_table_and_exits_zero_when_the_sky_path_fails(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE: THE SKY PATTERN MATCHES A FILE THAT IS NOT A _SKYSUB.fits PRODUCT
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    stareDir = tmp_path / "stare"
    stareDir.mkdir()
    (stareDir / "obs1_FLUXCAL.fits").write_bytes(b"")
    arcArguments = ["--arcs", str(rawDir / "*.fits"), "--disp-images", str(mapDir / "*_DISP_IMAGE.fits")]
    spec_tool.main(arcArguments)
    tableWithoutSky = yaml.safe_load(capsys.readouterr().out)

    # ACT
    status = spec_tool.main([*arcArguments, "--sky-products", str(stareDir / "*.fits")])

    # ASSERT
    captured = capsys.readouterr()
    assert status == 0
    assert yaml.safe_load(captured.out) == tableWithoutSky
    assert "skipped the sky-line cross-check" in captured.err
    assert "obs1_FLUXCAL.fits" in captured.err


def test_main_prints_the_arc_table_when_a_sky_frame_does_not_match_its_maps(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE: A SKY FRAME ONE ROW SHORT OF THE MAPS
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    stareDir = tmp_path / "stare"
    stareDir.mkdir()
    skyPath = write_sky_pair(stareDir / "obs1_SKYSUB.fits", Scene(minSlit=1.9), sky_header())
    for name in ("obs1_SKYSUB.fits", "obs1_SKYMODEL.fits"):
        with fits.open(stareDir / name, mode="update") as hdus:
            hdus[0].data = hdus[0].data[:-1]
            if len(hdus) > 1:
                hdus[1].data = hdus[1].data[:-1]

    # ACT
    status = spec_tool.main(
        [
            "--arcs",
            str(rawDir / "*.fits"),
            "--disp-images",
            str(mapDir / "*_DISP_IMAGE.fits"),
            "--sky-products",
            str(skyPath),
        ]
    )

    # ASSERT
    captured = capsys.readouterr()
    assert status == 0
    assert "spec-res:" in captured.out
    assert "SKY-LINE CROSS-CHECK" not in captured.out
    assert "obs1_SKYSUB.fits" in captured.err


def test_main_names_a_sky_arm_that_has_no_arc_arm_and_leaves_it_out_of_the_output(
    tmp_path: Path, capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, full_frame_trim: None
) -> None:
    # ARRANGE: VIS ARCS, BUT A NIR SKY FRAME WITH ITS OWN NIR DISP_IMAGE
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    nirScene = Scene(fwhmPx=3.0, minSlit=1.9).transposed()
    write_disp_image(mapDir / "nir_DISP_IMAGE.fits", nirScene, mjd=60000.1, arm="NIR")
    stareDir = tmp_path / "stare"
    stareDir.mkdir()
    write_sky_pair(stareDir / "obs1_SKYSUB.fits", nirScene, sky_header(arm="NIR"))

    # ACT
    status = spec_tool.main(
        [
            "--arcs",
            str(rawDir / "*.fits"),
            "--disp-images",
            str(mapDir / "*_DISP_IMAGE.fits"),
            "--sky-products",
            str(stareDir / "*_SKYSUB.fits"),
        ]
    )

    # ASSERT
    captured = capsys.readouterr()
    assert status == 0
    assert "NIR arm not recorded" in captured.err
    assert "NIR" not in yaml.safe_load(captured.out)
    assert "SKY-LINE CROSS-CHECK" not in captured.out
