"""*tests for the developer script that measures the nominal SPEC_RES of each arm and slit*

The line measurement runs on synthetic arc frames: a wavelength map that is linear along the dispersion axis, a slit
map in arcsec across it, and Gaussian emission lines of a known FWHM. The resolving power the script returns must match
`centre / FWHM` of the lines that were put in.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.io import fits

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import measure_spec_res as spec_tool  # noqa: E402

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]

ROWS = 60
COLUMNS = 2000
NM_PER_PIXEL = 0.03
START_NM = 500.0
BACKGROUND = 100.0
AMPLITUDE = 5000.0
FWHM_PER_SIGMA = 2.3548
LINE_COLUMNS = list(range(150, 1900, 150))
ROW_TILT = 0.37


class Scene:
    """*a synthetic arc frame with the maps that a DISP_IMAGE carries; dispersion along x*

    The wavelength map shifts by `ROW_TILT` pixels per row, so the pixels of one wavelength bin come from many columns,
    as they do on a real frame. The lines are Gaussian in wavelength, so they follow that map.
    """

    def __init__(self, fwhmPx: float = 3.0, lineColumns=LINE_COLUMNS, minSlit: float = 0.0, seed: int = 0):
        columns = np.arange(COLUMNS)
        rows = np.arange(ROWS)
        self.wavelength = START_NM + NM_PER_PIXEL * (columns[None, :] + ROW_TILT * rows[:, None])
        slitAxis = np.linspace(-6.0, 6.0, ROWS)
        self.slit = np.tile(slitAxis[:, None], (1, COLUMNS))
        self.order = np.ones((ROWS, COLUMNS), dtype=int)
        self.fwhmPx = fwhmPx
        self.centres = [START_NM + NM_PER_PIXEL * (column + ROW_TILT * ROWS / 2) for column in lineColumns]
        rng = np.random.default_rng(seed)
        data = BACKGROUND + rng.normal(0.0, 1.0, (ROWS, COLUMNS))
        sigmaNm = fwhmPx * NM_PER_PIXEL / FWHM_PER_SIGMA
        for centre in self.centres:
            line = AMPLITUDE * np.exp(-0.5 * ((self.wavelength - centre) / sigmaNm) ** 2)
            data += line * (np.abs(self.slit) >= minSlit)
        self.data = data
        self.good = np.ones((ROWS, COLUMNS), dtype=bool)

    def pixel_at(self, centreNm: float, row: int) -> tuple[int, int]:
        """*the (row, column) of the pixel in a row whose wavelength is closest to a wavelength*"""
        return row, int(round((centreNm - START_NM) / NM_PER_PIXEL - ROW_TILT * row))

    def measure(self, **overrides):
        arguments = {
            "data": self.data,
            "good": self.good,
            "wavelength": self.wavelength,
            "slit": self.slit,
            "order": self.order,
            "dispersionAxis": "x",
        }
        arguments.update(overrides)
        return spec_tool.measure_lines(**arguments)

    def expected_resolving_power(self, centreNm: float) -> float:
        return centreNm / (self.fwhmPx * NM_PER_PIXEL)


def test_measure_lines_recovers_the_resolving_power_of_isolated_lines():
    scene = Scene(fwhmPx=3.0)

    lines = scene.measure()

    assert len(lines) == len(LINE_COLUMNS)
    for line in lines:
        assert line.order == 1
        assert line.fwhmPx == pytest.approx(3.0, rel=0.03)
        assert line.resolvingPower == pytest.approx(scene.expected_resolving_power(line.centreNm), rel=0.03)


def test_measure_lines_reports_a_wider_line_as_a_lower_resolving_power():
    narrow = Scene(fwhmPx=2.5).measure()
    wide = Scene(fwhmPx=5.0).measure()

    assert np.median([line.resolvingPower for line in narrow]) == pytest.approx(
        2 * np.median([line.resolvingPower for line in wide]), rel=0.06
    )


def test_measure_lines_works_when_the_dispersion_axis_is_y():
    scene = Scene(fwhmPx=3.0)

    lines = scene.measure(
        data=scene.data.T,
        good=scene.good.T,
        wavelength=scene.wavelength.T,
        slit=scene.slit.T,
        order=scene.order.T,
        dispersionAxis="y",
    )

    assert len(lines) == len(LINE_COLUMNS)
    assert lines[0].fwhmPx == pytest.approx(3.0, rel=0.03)


def test_measure_lines_leaves_out_a_line_with_a_bad_pixel_in_its_fit_window():
    scene = Scene()
    scene.good[scene.pixel_at(scene.centres[3] + 0.02, row=ROWS // 2)] = False

    lines = scene.measure()

    assert len(lines) == len(LINE_COLUMNS) - 1
    assert all(abs(line.centreNm - scene.centres[3]) > 0.05 for line in lines)


def test_measure_lines_leaves_out_two_lines_that_blend():
    scene = Scene(lineColumns=[300, 306, 900, 1500])

    lines = scene.measure()

    assert sorted(round(line.centreNm) for line in lines) == sorted(
        round(c) for c in (scene.centres[2], scene.centres[3])
    )


def test_measure_lines_ignores_light_outside_the_slit_range():
    scene = Scene(minSlit=4.0)

    assert scene.measure() == []


def test_measure_lines_returns_nothing_when_no_pixel_belongs_to_an_order():
    scene = Scene()

    assert scene.measure(order=np.zeros_like(scene.order)) == []


def test_median_resolving_power_is_the_median_over_lines():
    lines = [spec_tool.LineFit(1, 500.0, 0.1, 3.0, power) for power in (4000.0, 5000.0, 9000.0)]

    assert spec_tool.median_resolving_power(lines) == 5000.0


def test_derive_wide_slit_scales_the_one_arcsec_value_by_the_inverse_width():
    values = {"SLIT0.5": 8965, "SLIT1.0": 4906, "SLIT1.5": 3271}

    table = spec_tool.derive_wide_slit(values)

    assert table["SLIT5.0"] == 981
    assert table["SLIT1.0"] == 4906


def test_derive_wide_slit_adds_nothing_when_the_one_arcsec_slit_was_not_measured():
    assert "SLIT5.0" not in spec_tool.derive_wide_slit({"SLIT0.5": 8965})


def _raw_header(arm="VIS", slit="SLIT1.0", mjd=60000.4, tech="ECHELLE,SLIT", dprType="LAMP,WAVE", binning=(1, 1)):
    header = fits.Header()
    header["ESO SEQ ARM"] = arm
    header["ESO INS VISE NAME" if arm == "VIS" else "ESO INS NISE NAME"] = slit
    header["ESO DPR TECH"] = tech
    header["ESO DPR TYPE"] = dprType
    header["ESO DET BINX"], header["ESO DET BINY"] = binning
    header["MJD-OBS"] = mjd
    return header


def _write_raw_arc(path: Path, scene: Scene, **cards) -> Path:
    fits.PrimaryHDU(scene.data.astype(np.float32), header=_raw_header(**cards)).writeto(path)
    return path


def _write_disp_image(path: Path, scene: Scene, mjd: float = 60000.1, arm: str = "VIS") -> Path:
    header = fits.Header()
    header["ESO SEQ ARM"] = arm
    header["MJD-OBS"] = mjd
    primary = fits.PrimaryHDU(scene.wavelength.astype(np.float32), header=header)
    primary.header["EXTNAME"] = "WAVELENGTH"
    slitHdu = fits.ImageHDU(scene.slit.astype(np.float32), name="SLIT")
    orderHdu = fits.ImageHDU(scene.order.astype(np.int32), name="ORDER")
    fits.HDUList([primary, slitHdu, orderHdu]).writeto(path)
    return path


@pytest.fixture
def full_frame_trim(monkeypatch):
    """*the synthetic frames are not the size of a SOXS detector, so keep every pixel*"""
    monkeypatch.setattr(spec_tool, "science_slices", lambda arm, binX, binY: (slice(None), slice(None)))


def test_is_slit_arc_accepts_a_slit_arc_in_either_dpr_type_order():
    assert spec_tool.is_slit_arc(_raw_header())
    assert spec_tool.is_slit_arc(_raw_header(dprType="WAVE,LAMP"))


def test_is_slit_arc_rejects_a_pinhole_arc_and_a_science_frame():
    assert not spec_tool.is_slit_arc(_raw_header(tech="ECHELLE,PINHOLE"))
    assert not spec_tool.is_slit_arc(_raw_header(dprType="OBJECT"))


@pytest.mark.parametrize(
    ("arm", "binX", "binY", "rows", "columns"),
    [
        ("VIS", 1, 1, (0, 4096), (50, 858)),
        ("VIS", 2, 1, (0, 4096), (25, 429)),
        ("VIS", 1, 2, (0, 2048), (50, 858)),
        ("NIR", 1, 1, (4, 2044), (4, 2044)),
    ],
)
def test_science_slices_follow_the_detector_parameters_and_the_binning(arm, binX, binY, rows, columns):
    rowSlice, columnSlice = spec_tool.science_slices(arm, binX, binY)

    assert (rowSlice.start, rowSlice.stop) == rows
    assert (columnSlice.start, columnSlice.stop) == columns


def test_read_arc_ref_returns_the_arm_slit_mjd_and_binning(tmp_path):
    path = _write_raw_arc(tmp_path / "arc.fits", Scene(), slit="SLIT1.5", mjd=60003.25, binning=(1, 1))

    arc = spec_tool.read_arc_ref(path)

    assert (arc.arm, arc.slit, arc.mjd, arc.binX, arc.binY) == ("VIS", "SLIT1.5", 60003.25, 1, 1)


def test_read_arc_ref_returns_none_for_a_frame_that_is_not_a_slit_arc(tmp_path):
    path = _write_raw_arc(tmp_path / "pin.fits", Scene(), tech="ECHELLE,PINHOLE")

    assert spec_tool.read_arc_ref(path) is None


def test_pair_arcs_with_maps_takes_the_nearest_map_of_the_same_arm():
    arcs = [spec_tool.ArcRef(Path("a"), "VIS", "SLIT1.0", 60000.4, 1, 1)]
    maps = [
        spec_tool.MapRef(Path("far"), "VIS", 59999.9),
        spec_tool.MapRef(Path("near"), "VIS", 60000.1),
        spec_tool.MapRef(Path("nir"), "NIR", 60000.4),
    ]

    pairs, unpaired = spec_tool.pair_arcs_with_maps(arcs, maps, maxGapDays=1.0)

    assert [(arc.path, mapRef.path) for arc, mapRef in pairs] == [(Path("a"), Path("near"))]
    assert unpaired == []


def test_pair_arcs_with_maps_leaves_an_arc_with_no_map_in_range_unpaired():
    arcs = [spec_tool.ArcRef(Path("a"), "VIS", "SLIT1.0", 60010.0, 1, 1)]
    maps = [spec_tool.MapRef(Path("m"), "VIS", 60000.1)]

    pairs, unpaired = spec_tool.pair_arcs_with_maps(arcs, maps, maxGapDays=1.0)

    assert pairs == []
    assert [arc.path for arc in unpaired] == [Path("a")]


def test_measure_arc_reads_the_frame_and_the_maps_and_returns_the_lines(tmp_path, full_frame_trim):
    scene = Scene(fwhmPx=3.0)
    arc = spec_tool.read_arc_ref(_write_raw_arc(tmp_path / "arc.fits", scene))
    mapRef = spec_tool.read_map_ref(_write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    lines = spec_tool.measure_arc(arc, mapRef, saturation=40000.0)

    assert len(lines) == len(LINE_COLUMNS)


def test_measure_arc_leaves_out_lines_that_reach_the_saturation_level(tmp_path, full_frame_trim):
    scene = Scene()
    arc = spec_tool.read_arc_ref(_write_raw_arc(tmp_path / "arc.fits", scene))
    mapRef = spec_tool.read_map_ref(_write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    assert spec_tool.measure_arc(arc, mapRef, saturation=3000.0) == []


def test_measure_arc_raises_when_the_trimmed_frame_does_not_match_the_maps(tmp_path):
    # NO `full_frame_trim`: THE REAL VIS TRIM IS LARGER THAN THE SYNTHETIC FRAME, SO THE SHAPES DISAGREE
    scene = Scene()
    arc = spec_tool.read_arc_ref(_write_raw_arc(tmp_path / "arc.fits", scene))
    mapRef = spec_tool.read_map_ref(_write_disp_image(tmp_path / "map_DISP_IMAGE.fits", scene))

    with pytest.raises(spec_tool.ToolError, match="shape"):
        spec_tool.measure_arc(arc, mapRef, saturation=40000.0)


def _lines(*powers: float) -> list:
    return [spec_tool.LineFit(1, 500.0, 0.1, 3.0, power) for power in powers]


def _measurements():
    return [
        spec_tool.ArcMeasurement("VIS", "SLIT1.0", 60000.4, _lines(4800, 5000, 5200)),
        spec_tool.ArcMeasurement("VIS", "SLIT1.0", 60001.4, _lines(4900, 5100, 5300, 5500)),
        spec_tool.ArcMeasurement("VIS", "SLIT0.5", 60000.5, _lines(9000, 9400)),
    ]


def test_summarise_takes_the_median_over_all_lines_of_a_slit_and_counts_them():
    summary = spec_tool.summarise(_measurements())

    slit = summary["VIS"]["SLIT1.0"]
    assert slit.resolvingPower == 5100
    assert slit.lineCount == 7
    assert slit.frameCount == 2
    assert slit.nightMedians == {60000: 5000, 60001: 5200}


def test_summarise_leaves_out_a_frame_that_found_no_lines():
    measurements = [*_measurements(), spec_tool.ArcMeasurement("VIS", "SLIT1.5", 60002.4, [])]

    assert "SLIT1.5" not in spec_tool.summarise(measurements)["VIS"]


def test_summarise_raises_when_no_frame_found_a_line():
    with pytest.raises(spec_tool.ToolError, match="no arc line"):
        spec_tool.summarise([spec_tool.ArcMeasurement("VIS", "SLIT1.0", 60000.4, [])])


def test_format_yaml_prints_the_table_in_slit_order_with_the_derived_wide_slit():
    summary = spec_tool.summarise(_measurements())

    text = spec_tool.format_yaml(summary, version="0.19.1", today="2026-11-22")

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


def test_format_yaml_block_loads_back_as_integer_resolving_powers():
    summary = spec_tool.summarise(_measurements())

    block = yaml.safe_load(spec_tool.format_yaml(summary, version="0.19.1", today="2026-11-22"))

    assert block == {"VIS": {"spec-res": {"SLIT0.5": 9200, "SLIT1.0": 5100, "SLIT5.0": 1020}}}


def test_format_yaml_says_when_the_wide_slit_cannot_be_derived():
    measurements = [spec_tool.ArcMeasurement("VIS", "SLIT0.5", 60000.5, _lines(9000, 9400))]

    text = spec_tool.format_yaml(spec_tool.summarise(measurements), version="0.19.1", today="2026-11-22")

    assert "# SLIT5.0 NOT DERIVED: NO SLIT1.0 MEASUREMENT" in text
    assert "SLIT5.0:" not in text.replace("# SLIT5.0", "")


def _write_night(tmp_path: Path) -> tuple[Path, Path]:
    """*one night: SLIT1.0 (3 px) and SLIT1.5 (4.5 px) arcs, a SLIT5.0 arc, a pinhole arc and one DISP_IMAGE*"""
    rawDir, mapDir = tmp_path / "raw", tmp_path / "maps"
    rawDir.mkdir()
    mapDir.mkdir()
    _write_raw_arc(rawDir / "slit10.fits", Scene(fwhmPx=3.0, seed=1), slit="SLIT1.0", mjd=60000.40)
    _write_raw_arc(rawDir / "slit15.fits", Scene(fwhmPx=4.5, seed=2), slit="SLIT1.5", mjd=60000.41)
    _write_raw_arc(rawDir / "slit50.fits", Scene(fwhmPx=6.0, seed=3), slit="SLIT5.0", mjd=60000.42)
    _write_raw_arc(rawDir / "pinhole.fits", Scene(fwhmPx=2.0, seed=4), tech="ECHELLE,PINHOLE", mjd=60000.43)
    _write_disp_image(mapDir / "night_DISP_IMAGE.fits", Scene(), mjd=60000.1)
    return rawDir, mapDir


def test_main_prints_the_measured_table_with_the_derived_wide_slit(tmp_path, capsys, monkeypatch, full_frame_trim):
    rawDir, mapDir = _write_night(tmp_path)
    monkeypatch.setattr(spec_tool, "today_iso", lambda: "2026-11-22")
    scene = Scene()
    expected10 = scene.expected_resolving_power(np.median(scene.centres))

    status = spec_tool.main(["--arcs", str(rawDir / "*.fits"), "--disp-images", str(mapDir / "*_DISP_IMAGE.fits")])

    captured = capsys.readouterr()
    table = yaml.safe_load(captured.out)["VIS"]["spec-res"]
    assert status == 0
    assert list(table) == ["SLIT1.0", "SLIT1.5", "SLIT5.0"]
    assert table["SLIT1.0"] == pytest.approx(expected10, rel=0.03)
    assert table["SLIT1.5"] == pytest.approx(expected10 / 1.5, rel=0.04)
    assert table["SLIT5.0"] == round(table["SLIT1.0"] / 5)
    assert "slit50.fits" in captured.err
    assert "pinhole.fits" in captured.err


def test_main_names_the_arcs_that_have_no_disp_image_in_range_and_stops_when_none_is_paired(
    tmp_path, capsys, full_frame_trim
):
    rawDir, mapDir = _write_night(tmp_path)

    with pytest.raises(spec_tool.ToolError, match="paired"):
        spec_tool.main(
            ["--arcs", str(rawDir / "slit10.fits"), "--disp-images", str(mapDir / "*.fits"), "--max-gap-days", "0.1"]
        )

    assert "slit10.fits" in capsys.readouterr().err
