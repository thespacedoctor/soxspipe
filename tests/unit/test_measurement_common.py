"""*tests for the helpers shared by the developer measurement scripts in tools/*"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from astropy.io import fits

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import measurement_common as common  # noqa: E402

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]


def _touch(directory: Path, *names: str) -> None:
    for name in names:
        (directory / name).write_bytes(b"")


def test_expand_paths_returns_sorted_unique_matches_of_every_pattern(tmp_path):
    _touch(tmp_path, "b.fits", "a.fits", "c.txt")

    paths = common.expand_paths([str(tmp_path / "*.fits"), str(tmp_path / "a.fits")])

    assert paths == [tmp_path / "a.fits", tmp_path / "b.fits"]


def test_expand_paths_expands_a_leading_tilde(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    _touch(tmp_path, "x.fits")

    assert common.expand_paths(["~/x.fits"]) == [tmp_path / "x.fits"]


def test_expand_paths_raises_when_a_pattern_matches_nothing(tmp_path):
    with pytest.raises(common.ToolError, match="no file matches"):
        common.expand_paths([str(tmp_path / "*.fits")])


def _header(**cards) -> fits.Header:
    header = fits.Header()
    for key, value in cards.items():
        header[key.replace("_", " ")] = value
    return header


def test_arm_of_reads_the_sequencer_arm_card():
    assert common.arm_of(_header(ESO_SEQ_ARM="NIR")) == "NIR"


def test_arm_of_raises_for_a_missing_or_unknown_arm():
    with pytest.raises(common.ToolError, match="ESO SEQ ARM"):
        common.arm_of(_header())
    with pytest.raises(common.ToolError, match="UVB"):
        common.arm_of(_header(ESO_SEQ_ARM="UVB"))


@pytest.mark.parametrize(
    ("arm", "cards", "binning"),
    [
        ("VIS", {"ESO_DET_BINX": 2, "ESO_DET_BINY": 1}, 2),
        ("NIR", {"ESO_DET_BINX": 1, "ESO_DET_BINY": 1}, 1),
    ],
)
def test_spectral_binning_reads_the_card_of_the_dispersion_axis(arm, cards, binning):
    assert common.spectral_binning(_header(**cards), arm) == binning


def test_spectral_binning_uses_the_y_card_when_the_dispersion_axis_is_y():
    # THE NIR ARM DISPERSES ALONG y, SO ITS SPECTRAL BINNING IS BINY
    assert common.spectral_binning(_header(ESO_DET_BINX=4, ESO_DET_BINY=2), "NIR") == 2


def test_spectral_binning_raises_when_the_card_is_missing_or_below_one():
    with pytest.raises(common.ToolError, match="ESO DET BINX"):
        common.spectral_binning(_header(), "VIS")
    with pytest.raises(common.ToolError, match="ESO DET BINX"):
        common.spectral_binning(_header(ESO_DET_BINX=0), "VIS")


def test_night_count_counts_distinct_utc_dates():
    assert common.night_count([60000.1, 60000.9, 60001.5, 60003.0]) == 3


def test_provenance_comments_name_the_subject_the_nights_and_the_version():
    lines = common.provenance_comments("FRAMES", [60000.1, 60001.2, 60001.4], "0.20.0", "2026-11-22")

    assert lines == [
        "# MEASURED 2026-11-22 WITH soxspipe 0.20.0 FROM 3 FRAMES ON 2 NIGHTS (MJD 60000 TO 60001)",
    ]


def test_run_main_returns_the_function_status():
    assert common.run_main(lambda argv: 0, []) == 0


def test_run_main_reports_a_tool_error_on_stderr_and_returns_two(capsys):
    def failing(argv):
        raise common.ToolError("no file matches x")

    status = common.run_main(failing, [])

    assert status == 2
    assert "no file matches x" in capsys.readouterr().err


def test_soxspipe_version_reads_the_generated_version_file(tmp_path):
    versionFile = tmp_path / "__version__.py"
    versionFile.write_text('__version__ = "1.2.3"\n')

    assert common.soxspipe_version(versionFile) == "1.2.3"


def test_soxspipe_version_is_unknown_when_the_file_is_missing_or_has_no_version(tmp_path):
    emptyFile = tmp_path / "__version__.py"
    emptyFile.write_text("# NOTHING\n")

    assert common.soxspipe_version(tmp_path / "missing.py") == "unknown"
    assert common.soxspipe_version(emptyFile) == "unknown"


def test_today_iso_is_an_iso_date():
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", common.today_iso())


def test_ensure_repo_on_path_puts_the_checkout_on_the_import_path(monkeypatch):
    monkeypatch.setattr(sys, "path", [entry for entry in sys.path if entry != str(common.REPO_ROOT)])

    common.ensure_repo_on_path()

    assert sys.path[0] == str(common.REPO_ROOT)


def test_read_header_returns_the_primary_header(tmp_path):
    path = tmp_path / "a.fits"
    fits.PrimaryHDU(header=_header(ESO_SEQ_ARM="VIS")).writeto(path)

    assert common.read_header(path)["ESO SEQ ARM"] == "VIS"


def test_read_header_raises_a_tool_error_for_a_file_that_is_not_fits(tmp_path):
    path = tmp_path / "not_fits.fits"
    path.write_text("plain text")

    with pytest.raises(common.ToolError, match="cannot read"):
        common.read_header(path)
