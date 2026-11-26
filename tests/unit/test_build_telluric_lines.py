"""*tests for the developer script that builds the telluric line list*

The HITRAN fetch is replaced by a synthetic span with one isolated water line, so the build runs offline in seconds.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from tests.unit._tools_loader import load_tool

build_tool = load_tool("build_telluric_lines")

pytestmark = pytest.mark.unit

LINE_NM = 1130.0
# A WATER LINE OF ABOUT 15 PERCENT DEPTH AFTER SMOOTHING TO R = 3500, SO IT PASSES THE DEPTH LIMITS AT EVERY COLUMN
LINE_STRENGTH = 3e-23
RESOLVING_POWERS = {"VIS": (3500.0, 5000.0), "NIR": (3500.0, 5000.0)}


def _selection(minBandLines: int = 1, resolvingPowers: dict | None = None) -> object:
    return build_tool.Selection(
        resolvingPowers=RESOLVING_POWERS if resolvingPowers is None else resolvingPowers,
        minDepth=0.04,
        maxDepth=0.85,
        maxSpreadNm=0.015,
        maxOffsetNm=0.05,
        mergeNm=0.03,
        minBandLines=minBandLines,
    )


def _span(lineNms: list[float], strengths: list[float], molecule: int = 1) -> dict:
    """A fetched span that holds one line per entry, at the given air wavelengths."""
    airNm = np.array(lineNms, dtype=float)
    nu = np.array([1e7 / build_tool.vac_to_air_inverse(wavelength) for wavelength in airNm])
    count = nu.size
    return {
        "nu": nu,
        "sw": np.array(strengths, dtype=float),
        "gamma_air": np.full(count, 0.08),
        "n_air": np.full(count, 0.7),
        "elower": np.full(count, 100.0),
        "delta_air": np.zeros(count),
        "molecule": molecule,
    }


def _tables(lineNms: list[float] | None = None, strengths: list[float] | None = None) -> dict:
    lineNms = [LINE_NM] if lineNms is None else lineNms
    strengths = [LINE_STRENGTH] * len(lineNms) if strengths is None else strengths
    return {"H2O_VIS": _span([], []), "O2_VIS": _span([], [], molecule=7), "H2O_NIR": _span(lineNms, strengths)}


def test_the_resolving_powers_run_from_the_lowest_to_the_highest_in_equal_steps_of_the_logarithm() -> None:
    # ACT
    powers = build_tool.resolving_power_grid(954.0, 5000.0, 3)

    # ASSERT
    assert powers == (954.0, 2184.0, 5000.0)


@pytest.mark.parametrize(("lowest", "highest", "count"), [(5000.0, 954.0, 3), (0.0, 5000.0, 3), (954.0, 5000.0, 1)])
def test_resolving_powers_that_cannot_span_a_range_are_refused(lowest: float, highest: float, count: int) -> None:
    # ACT AND ASSERT
    with pytest.raises(ValueError, match="cannot span resolving powers"):
        build_tool.resolving_power_grid(lowest, highest, count)


def test_vacuum_wavelengths_become_air_wavelengths_with_the_standard_refractivity() -> None:
    # ACT
    airNm = build_tool.vac_to_air_nm(1000.0)

    # ASSERT
    assert airNm == pytest.approx(999.726, abs=1e-3)


def test_the_inverse_of_the_air_conversion_recovers_the_vacuum_wavelength() -> None:
    # ACT
    vacNm = build_tool.vac_to_air_inverse(1000.0)

    # ASSERT
    assert float(build_tool.vac_to_air_nm(vacNm)) == pytest.approx(1000.0, abs=1e-9)


def test_a_voigt_profile_has_unit_area() -> None:
    # ARRANGE
    offsets = np.arange(-200.0, 200.0, 0.001)

    # ACT
    profile = build_tool.voigt_profile(offsets, doppler=0.012, lorentz=0.06)

    # ASSERT
    assert np.trapezoid(profile, offsets) == pytest.approx(1.0, abs=5e-3)


def test_smoothing_keeps_the_continuum_and_fills_in_a_narrow_dip() -> None:
    # ARRANGE
    grid = np.arange(8000.0, 8100.0, build_tool.GRID_STEP)
    signal = np.ones_like(grid)
    signal[np.abs(grid - 8050.0) < 0.02] = 0.0

    # ACT
    smoothed = build_tool.convolve_to_resolution(grid, signal, 2000.0)

    # ASSERT
    assert smoothed[np.argmin(np.abs(grid - 8025.0))] == pytest.approx(1.0, abs=0.01)
    assert smoothed[np.argmin(np.abs(grid - 8050.0))] > 0.5


def test_only_the_strong_lines_inside_the_band_are_candidates() -> None:
    # ARRANGE
    span = _span([1100.0, 1120.0, 1130.0, 1200.0], [3e-23, 3e-23, 1e-30, 3e-23])

    # ACT
    candidates = build_tool.band_candidates([span], 1110.0, 1160.0)

    # ASSERT
    assert [round(wavelength) for wavelength, _, _ in candidates] == [1120]
    assert candidates[0][2] == "H2O"


def test_of_two_features_closer_than_the_merge_distance_only_the_stronger_stays() -> None:
    # ARRANGE
    features = [(1130.00, 1130.0, 1.0, "H2O"), (1130.01, 1130.01, 3.0, "H2O"), (1131.0, 1131.0, 1.0, "H2O")]

    # ACT
    kept = build_tool.merge_close_features(features, 0.03)

    # ASSERT
    assert [feature[1] for feature in kept] == [1130.01, 1131.0]


@pytest.mark.parametrize(
    ("centres", "depths", "expected"),
    [
        ([1130.0, 1130.005], [0.1, 0.2], True),
        ([1130.0, 1130.005], [0.02, 0.2], False),
        ([1130.0, 1130.005], [0.1, 0.9], False),
        ([1130.0, 1130.03], [0.1, 0.2], False),
        ([1130.08, 1130.085], [0.1, 0.2], False),
    ],
)
def test_a_line_passes_when_its_depth_spread_and_offset_are_inside_the_limits(
    centres: list[float], depths: list[float], expected: bool
) -> None:
    # ACT
    passes = build_tool.passes_selection(np.array(centres), np.array(depths), 1130.0, _selection())

    # ASSERT
    assert passes is expected


def test_an_isolated_line_is_kept_at_its_own_wavelength() -> None:
    # ACT
    kept = build_tool.select_band(_tables(), ("H2O_1130", 1110.0, 1160.0, "NIR", ("H2O_NIR",)), _selection())

    # ASSERT
    assert len(kept) == 1
    effective, hitranNm, strength, species = kept[0]
    assert effective == pytest.approx(LINE_NM, abs=0.01)
    assert hitranNm == pytest.approx(LINE_NM, abs=1e-6)
    assert species == "H2O"


def test_a_band_with_no_hitran_line_keeps_nothing() -> None:
    # ACT
    kept = build_tool.select_band(_tables([]), ("H2O_1130", 1110.0, 1160.0, "NIR", ("H2O_NIR",)), _selection())

    # ASSERT
    assert kept == []


def test_two_lines_that_blend_move_the_fitted_centre_and_are_not_kept_as_isolated_lines() -> None:
    # ARRANGE
    tables = _tables([LINE_NM, LINE_NM + 0.35], [LINE_STRENGTH, LINE_STRENGTH * 0.6])

    # ACT
    kept = build_tool.select_band(tables, ("H2O_1130", 1110.0, 1160.0, "NIR", ("H2O_NIR",)), _selection())

    # ASSERT
    assert kept == []


def test_a_band_with_fewer_lines_than_the_minimum_is_dropped_and_named() -> None:
    # ARRANGE
    rows = {"A": [(1,)] * 4, "B": [(1,)] * 5, "C": []}

    # ACT
    kept, dropped = build_tool.drop_sparse_bands(rows, 5)

    # ASSERT
    assert list(kept) == ["B"]
    assert dropped == ["A", "C"]


def test_the_build_drops_the_bands_with_too_few_lines_and_reports_them() -> None:
    # ACT
    table, dropped = build_tool.build_line_list(_tables(), _selection(minBandLines=1))

    # ASSERT
    assert list(table["BAND"]) == ["H2O_1130"]
    assert list(table["ARM"]) == ["NIR"]
    assert dropped == ["O2_B", "H2O_720", "H2O_940", "H2O_1450", "H2O_1800"]


def test_the_header_cards_record_the_model_settings_and_the_selection_limits() -> None:
    # ACT
    cards = {keyword: value for keyword, value, _ in build_tool.model_header_cards(_selection(5), ["O2_B"])}

    # ASSERT
    assert cards["HIERARCH MODEL TEMPERATURE"] == 270.0
    assert cards["HIERARCH MODEL PRESSURE"] == 0.75
    assert cards["HIERARCH MODEL AIRMASS"] == 1.3
    assert cards["HIERARCH MODEL H2O COLUMN"] == 1.3e22
    assert cards["HIERARCH MODEL H2O COLUMN FACTORS"] == "0.6,1,1.6"
    assert cards["HIERARCH MODEL RESOLVING POWERS VIS"] == "3500,5000"
    assert cards["HIERARCH MODEL RESOLVING POWERS NIR"] == "3500,5000"
    assert cards["HIERARCH SELECT MIN DEPTH"] == 0.04
    assert cards["HIERARCH SELECT MAX DEPTH"] == 0.85
    assert cards["HIERARCH SELECT MAX SPREAD"] == 0.015
    assert cards["HIERARCH SELECT MAX OFFSET"] == 0.05
    assert cards["HIERARCH SELECT MERGE DISTANCE"] == 0.03
    assert cards["HIERARCH SELECT MIN BAND LINES"] == 5
    assert cards["HIERARCH SELECT DROPPED BANDS"] == "O2_B"


def test_the_written_file_holds_the_lines_and_the_cards_and_stays_small(tmp_path: Path) -> None:
    # ARRANGE
    table, dropped = build_tool.build_line_list(_tables(), _selection())
    path = tmp_path / "telluric_lines.fits"

    # ACT
    build_tool.write_line_list(table, _selection(), dropped, str(path))

    # ASSERT
    with fits.open(path) as hdus:
        header = hdus["TELLURIC_LINES"].header
        assert header["HIERARCH MODEL TEMPERATURE"] == 270.0
        assert header["HIERARCH SELECT DROPPED BANDS"] == ",".join(dropped)
        assert list(hdus["TELLURIC_LINES"].data["BAND"]) == ["H2O_1130"]
    assert path.stat().st_size < 20000


def test_hapi_fetches_only_the_spans_that_the_cache_lacks(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    # ARRANGE
    fetched: list[str] = []
    cached = {"H2O_VIS"}
    columns = {"nu": [10000.0], "sw": [1e-24], "gamma_air": [0.08], "n_air": [0.7], "elower": [1.0], "delta_air": [0.0]}
    fakeHapi = types.SimpleNamespace(
        db_begin=lambda cache: None,
        tableList=lambda: list(cached),
        fetch=lambda name, *args: fetched.append(name),
        getColumn=lambda name, column: columns[column],
    )
    monkeypatch.setitem(sys.modules, "hapi", fakeHapi)

    # ACT
    tables = build_tool.fetch_tables(str(tmp_path / "cache"))

    # ASSERT
    assert fetched == ["O2_VIS", "H2O_NIR"]
    assert tables["O2_VIS"]["molecule"] == 7
    assert tables["H2O_VIS"]["nu"].tolist() == [10000.0]


def test_the_command_line_builds_the_file_and_prints_the_bands_and_the_resolving_powers(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # ARRANGE
    monkeypatch.setattr(build_tool, "fetch_tables", lambda cache: _tables())
    out = tmp_path / "lines.fits"
    argv = ["--cache", str(tmp_path), "--out", str(out), "--lowest-r-vis", "3500", "--lowest-r-nir", "3500"]
    argv += ["--highest-r", "5000", "--min-band-lines", "1"]

    # ACT
    status = build_tool.main(argv)

    # ASSERT
    printed = capsys.readouterr().out
    assert status == 0
    assert "NIR resolving powers: 3500, 4183, 5000" in printed
    assert "H2O_1130" in printed
    assert "dropped bands: O2_B, H2O_720, H2O_940, H2O_1450, H2O_1800" in printed
    assert len(Table.read(out, format="fits")) == 1
