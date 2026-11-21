"""*tests for the frame reading, pairing and measurement behind the SPEC_RES script of tools/*

The line measurement itself is tested in `test_line_fitting.py`. These tests cover the frames that the script reads:
the slit arcs, the reduced stare frames with their sky models, the `DISP_IMAGE` files, their pairing and the trim.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from tests.unit._arc_scene import (
    LINE_COLUMNS,
    Scene,
    raw_header,
    sky_header,
    write_disp_image,
    write_raw_arc,
    write_sky_pair,
)
from tests.unit._tools_loader import load_tool

frames_tool = load_tool("spec_res_frames")

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]


@pytest.fixture
def full_frame_trim(monkeypatch: pytest.MonkeyPatch) -> None:
    """*the synthetic frames are not the size of a SOXS detector, so keep every pixel*"""
    monkeypatch.setattr(frames_tool, "science_slices", lambda arm, binX, binY: (slice(None), slice(None)))


def test_is_slit_arc_accepts_a_slit_arc_in_either_dpr_type_order() -> None:
    # ARRANGE
    lampWave, waveLamp = raw_header(), raw_header(dprType="WAVE,LAMP")

    # ACT / ASSERT
    assert frames_tool.is_slit_arc(lampWave)
    assert frames_tool.is_slit_arc(waveLamp)


def test_is_slit_arc_rejects_a_pinhole_arc_and_a_science_frame() -> None:
    # ARRANGE
    pinhole, science = raw_header(tech="ECHELLE,PINHOLE"), raw_header(dprType="OBJECT")

    # ACT / ASSERT
    assert not frames_tool.is_slit_arc(pinhole)
    assert not frames_tool.is_slit_arc(science)


@pytest.mark.parametrize(
    ("arm", "binX", "binY", "rows", "columns"),
    [
        ("VIS", 1, 1, (0, 4096), (50, 858)),
        ("VIS", 2, 1, (0, 4096), (25, 429)),
        ("VIS", 1, 2, (0, 2048), (50, 858)),
        ("NIR", 1, 1, (4, 2044), (4, 2044)),
    ],
)
def test_science_slices_follow_the_detector_parameters_and_the_binning(
    arm: str, binX: int, binY: int, rows: tuple[int, int], columns: tuple[int, int]
) -> None:
    # ARRANGE / ACT
    rowSlice, columnSlice = frames_tool.science_slices(arm, binX, binY)

    # ASSERT
    assert (rowSlice.start, rowSlice.stop) == rows
    assert (columnSlice.start, columnSlice.stop) == columns


def test_read_arc_ref_returns_the_arm_slit_mjd_and_binning(tmp_path: Path) -> None:
    # ARRANGE
    path = write_raw_arc(tmp_path / "arc.fits", Scene(), raw_header(slit="SLIT1.5", mjd=60003.25))

    # ACT
    arc = frames_tool.read_arc_ref(path)

    # ASSERT
    assert (arc.arm, arc.slit, arc.mjd, arc.binX, arc.binY) == ("VIS", "SLIT1.5", 60003.25, 1, 1)


def test_read_arc_ref_returns_none_for_a_frame_that_is_not_a_slit_arc(tmp_path: Path) -> None:
    # ARRANGE
    path = write_raw_arc(tmp_path / "pin.fits", Scene(), raw_header(tech="ECHELLE,PINHOLE"))

    # ACT
    arc = frames_tool.read_arc_ref(path)

    # ASSERT
    assert arc is None


@pytest.mark.parametrize("missingCard", ["ESO DET BINX", "ESO DET BINY", "MJD-OBS", "ESO INS VISE NAME"])
def test_read_arc_ref_names_the_file_and_the_card_when_a_header_card_is_missing(
    tmp_path: Path, missingCard: str
) -> None:
    # ARRANGE
    header = raw_header()
    del header[missingCard]
    path = write_raw_arc(tmp_path / "arc.fits", Scene(), header)

    # ACT / ASSERT
    with pytest.raises(frames_tool.ToolError, match=re.escape(f"{path} has no {missingCard} card")):
        frames_tool.read_arc_ref(path)


def test_read_map_ref_names_the_file_when_the_mjd_card_is_missing(tmp_path: Path) -> None:
    # ARRANGE
    path = write_disp_image(tmp_path / "map_DISP_IMAGE.fits", Scene())
    with fits.open(path, mode="update") as hdus:
        del hdus[0].header["MJD-OBS"]

    # ACT / ASSERT
    with pytest.raises(frames_tool.ToolError, match=re.escape(f"{path} has no MJD-OBS card")):
        frames_tool.read_map_ref(path)


def test_pair_frames_with_maps_takes_the_nearest_map_of_the_same_arm() -> None:
    # ARRANGE
    arcs = [frames_tool.ArcRef(Path("a"), "VIS", "SLIT1.0", 60000.4, 1, 1)]
    maps = [
        frames_tool.MapRef(Path("far"), "VIS", 59999.9),
        frames_tool.MapRef(Path("near"), "VIS", 60000.1),
        frames_tool.MapRef(Path("nir"), "NIR", 60000.4),
    ]

    # ACT
    pairs, unpaired = frames_tool.pair_frames_with_maps(arcs, maps, maxGapDays=1.0)

    # ASSERT
    assert [(arc.path, mapRef.path) for arc, mapRef in pairs] == [(Path("a"), Path("near"))]
    assert unpaired == []


def test_pair_frames_with_maps_leaves_a_frame_with_no_map_in_range_unpaired() -> None:
    # ARRANGE
    arcs = [frames_tool.ArcRef(Path("a"), "VIS", "SLIT1.0", 60010.0, 1, 1)]
    maps = [frames_tool.MapRef(Path("m"), "VIS", 60000.1)]

    # ACT
    pairs, unpaired = frames_tool.pair_frames_with_maps(arcs, maps, maxGapDays=1.0)

    # ASSERT
    assert pairs == []
    assert [arc.path for arc in unpaired] == [Path("a")]


def test_measure_arc_reads_the_frame_and_the_maps_and_returns_the_lines(tmp_path: Path, full_frame_trim: None) -> None:
    # ARRANGE
    scene = Scene(fwhmPx=3.0)
    arc = frames_tool.read_arc_ref(write_raw_arc(tmp_path / "arc.fits", scene, raw_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_arc(arc, mapRef, saturation=40000.0)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS)


def test_measure_arc_leaves_out_lines_that_reach_the_saturation_level(tmp_path: Path, full_frame_trim: None) -> None:
    # ARRANGE
    scene = Scene()
    arc = frames_tool.read_arc_ref(write_raw_arc(tmp_path / "arc.fits", scene, raw_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_arc(arc, mapRef, saturation=3000.0)

    # ASSERT
    assert lines == []


def test_measure_arc_raises_when_the_trimmed_frame_does_not_match_the_maps(tmp_path: Path) -> None:
    # ARRANGE: NO `full_frame_trim`, SO THE REAL VIS TRIM IS LARGER THAN THE SYNTHETIC FRAME AND THE SHAPES DISAGREE
    scene = Scene()
    arc = frames_tool.read_arc_ref(write_raw_arc(tmp_path / "arc.fits", scene, raw_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT / ASSERT
    with pytest.raises(frames_tool.ToolError, match="shape"):
        frames_tool.measure_arc(arc, mapRef, saturation=40000.0)


def test_read_sky_ref_finds_the_sky_model_beside_the_sky_subtracted_frame(tmp_path: Path) -> None:
    # ARRANGE
    path = write_sky_pair(tmp_path / "stare_SKYSUB.fits", Scene(), sky_header(slit="SLIT1.5"))

    # ACT
    sky = frames_tool.read_sky_ref(path)

    # ASSERT
    assert (sky.arm, sky.slit, sky.mjd) == ("VIS", "SLIT1.5", 60000.45)
    assert sky.modelPath == tmp_path / "stare_SKYMODEL.fits"


def test_read_sky_ref_returns_none_when_the_sky_model_is_missing(tmp_path: Path) -> None:
    # ARRANGE
    path = write_sky_pair(tmp_path / "stare_SKYSUB.fits", Scene(), sky_header())
    (tmp_path / "stare_SKYMODEL.fits").unlink()

    # ACT
    sky = frames_tool.read_sky_ref(path)

    # ASSERT
    assert sky is None


def test_read_sky_ref_raises_for_a_file_that_is_not_a_skysub_product(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "stare_FLUXCAL.fits"

    # ACT / ASSERT
    with pytest.raises(frames_tool.ToolError, match="_SKYSUB.fits"):
        frames_tool.read_sky_ref(path)


def test_measure_sky_adds_the_sky_model_back_and_measures_lines_2_to_5_arcsec_from_the_centre(tmp_path: Path) -> None:
    # ARRANGE: SKY LINES ONLY MORE THAN 1.9 ARCSEC FROM THE CENTRE, ON A FRAME THE SIZE OF THE MAPS
    scene = Scene(fwhmPx=3.0, minSlit=1.9)
    sky = frames_tool.read_sky_ref(write_sky_pair(tmp_path / "stare_SKYSUB.fits", scene, sky_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_sky(sky, mapRef)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS)
    assert np.median([line.resolvingPower for line in lines]) == pytest.approx(
        scene.expected_resolving_power(np.median(scene.centres)), rel=0.03
    )


def test_measure_sky_leaves_out_a_line_that_the_quality_extension_flags(tmp_path: Path) -> None:
    # ARRANGE: ONE FLAGGED PIXEL AT +3 ARCSEC, INSIDE THE FIT WINDOW OF THE FOURTH LINE
    scene = Scene(fwhmPx=3.0, minSlit=1.9)
    row = int(np.argmin(np.abs(scene.slit[:, 0] - 3.0)))
    badPixel = scene.pixel_at(scene.centres[3] + 0.02, row)
    sky = frames_tool.read_sky_ref(write_sky_pair(tmp_path / "stare_SKYSUB.fits", scene, sky_header(), badPixel))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_sky(sky, mapRef)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS) - 1


def test_measure_sky_ignores_lines_that_sit_inside_2_arcsec_of_the_slit_centre(tmp_path: Path) -> None:
    # ARRANGE: THE TARGET REGION CARRIES A LINE, THE SKY REGION DOES NOT
    scene = Scene(fwhmPx=3.0, maxSlit=1.95)
    sky = frames_tool.read_sky_ref(write_sky_pair(tmp_path / "stare_SKYSUB.fits", scene, sky_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_sky(sky, mapRef)

    # ASSERT
    assert lines == []


def test_read_arc_ref_raises_naming_the_file_when_the_binning_card_is_not_a_number(tmp_path: Path) -> None:
    # ARRANGE
    header = raw_header()
    header["ESO DET BINX"] = "two"
    path = write_raw_arc(tmp_path / "arc.fits", Scene(), header)

    # ACT / ASSERT
    with pytest.raises(
        frames_tool.ToolError, match=re.escape(f"{path} has a ESO DET BINX card that is not a whole number")
    ):
        frames_tool.read_arc_ref(path)


def test_measure_sky_uses_the_positive_side_of_the_slit_only(tmp_path: Path) -> None:
    # ARRANGE: THE PROTOTYPE USED 2 < SLIT < 5 ONLY, SO LINES ON THE NEGATIVE SIDE ONLY ARE NOT FOUND
    scene = Scene(fwhmPx=3.0, minSlit=1.9, side=-1)
    sky = frames_tool.read_sky_ref(write_sky_pair(tmp_path / "stare_SKYSUB.fits", scene, sky_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_sky(sky, mapRef)

    # ASSERT
    assert lines == []


def test_measure_sky_finds_the_lines_of_the_positive_side_of_the_slit(tmp_path: Path) -> None:
    # ARRANGE
    scene = Scene(fwhmPx=3.0, minSlit=1.9, side=1)
    sky = frames_tool.read_sky_ref(write_sky_pair(tmp_path / "stare_SKYSUB.fits", scene, sky_header()))
    mapRef = frames_tool.read_map_ref(write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    # ACT
    lines = frames_tool.measure_sky(sky, mapRef)

    # ASSERT
    assert len(lines) == len(LINE_COLUMNS)
