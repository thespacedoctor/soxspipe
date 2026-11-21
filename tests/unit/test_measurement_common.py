"""*tests for the helpers shared by the developer measurement scripts in tools/*"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest
from astropy.io import fits

from tests.unit._tools_loader import load_tool

common = load_tool("measurement_common")

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]


def _touch(directory: Path, *names: str) -> None:
    """*create empty files*"""
    for name in names:
        (directory / name).write_bytes(b"")


def _header(**cards: object) -> fits.Header:
    """*a header from keyword arguments, with underscores in the names read as spaces*"""
    header = fits.Header()
    for key, value in cards.items():
        header[key.replace("_", " ")] = value
    return header


def test_expand_paths_returns_sorted_unique_matches_of_every_pattern(tmp_path: Path) -> None:
    # ARRANGE
    _touch(tmp_path, "b.fits", "a.fits", "c.txt")

    # ACT
    paths = common.expand_paths([str(tmp_path / "*.fits"), str(tmp_path / "a.fits")])

    # ASSERT
    assert paths == [tmp_path / "a.fits", tmp_path / "b.fits"]


def test_expand_paths_expands_a_leading_tilde(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    monkeypatch.setenv("HOME", str(tmp_path))
    _touch(tmp_path, "x.fits")

    # ACT
    paths = common.expand_paths(["~/x.fits"])

    # ASSERT
    assert paths == [tmp_path / "x.fits"]


def test_expand_paths_raises_when_a_pattern_matches_nothing(tmp_path: Path) -> None:
    # ARRANGE
    pattern = str(tmp_path / "*.fits")

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match="no file matches"):
        common.expand_paths([pattern])


def test_arm_of_reads_the_sequencer_arm_card() -> None:
    # ARRANGE
    header = _header(ESO_SEQ_ARM="NIR")

    # ACT
    arm = common.arm_of(header)

    # ASSERT
    assert arm == "NIR"


def test_arm_of_raises_for_a_missing_or_unknown_arm() -> None:
    # ARRANGE
    missing, unknown = _header(), _header(ESO_SEQ_ARM="UVB")

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match="ESO SEQ ARM"):
        common.arm_of(missing)
    with pytest.raises(common.ToolError, match="UVB"):
        common.arm_of(unknown)


def test_arm_of_names_the_file_when_the_arm_card_is_missing(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "a.fits"

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match=re.escape(f"{path} has no ESO SEQ ARM card")):
        common.arm_of(_header(), path)


@pytest.mark.parametrize(
    ("arm", "cards", "binning"),
    [
        ("VIS", {"ESO_DET_BINX": 2, "ESO_DET_BINY": 1}, 2),
        ("NIR", {"ESO_DET_BINX": 1, "ESO_DET_BINY": 1}, 1),
        # THE NIR ARM DISPERSES ALONG y, SO ITS SPECTRAL BINNING IS BINY
        ("NIR", {"ESO_DET_BINX": 4, "ESO_DET_BINY": 2}, 2),
    ],
)
def test_spectral_binning_reads_the_card_of_the_dispersion_axis(arm: str, cards: dict, binning: int) -> None:
    # ARRANGE
    header = _header(**cards)

    # ACT
    result = common.spectral_binning(header, arm)

    # ASSERT
    assert result == binning


def test_spectral_binning_raises_naming_the_file_and_card_when_the_card_is_missing(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "a.fits"

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match=re.escape(f"{path} has no ESO DET BINX card")):
        common.spectral_binning(_header(), "VIS", path)


def test_spectral_binning_raises_when_the_card_is_below_one() -> None:
    # ARRANGE
    header = _header(ESO_DET_BINX=0)

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match="ESO DET BINX"):
        common.spectral_binning(header, "VIS")


def test_required_card_returns_the_value_of_a_present_card() -> None:
    # ARRANGE
    header = _header(ESO_INS_VISE_NAME="SLIT1.0")

    # ACT
    value = common.required_card(header, "ESO INS VISE NAME")

    # ASSERT
    assert value == "SLIT1.0"


def test_read_mjd_obs_returns_the_mjd_as_a_float() -> None:
    # ARRANGE
    header = fits.Header()
    header["MJD-OBS"] = 60001.25

    # ACT
    mjd = common.read_mjd_obs(header)

    # ASSERT
    assert mjd == 60001.25


def test_read_mjd_obs_raises_naming_the_file_and_card_when_the_card_is_missing(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "frame.fits"

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match=re.escape(f"{path} has no MJD-OBS card")):
        common.read_mjd_obs(fits.Header(), path)


def test_read_mjd_obs_raises_when_the_card_is_not_a_number() -> None:
    # ARRANGE
    header = fits.Header()
    header["MJD-OBS"] = "tonight"

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match="not a number"):
        common.read_mjd_obs(header)


def test_night_count_counts_distinct_utc_dates() -> None:
    # ARRANGE
    mjds = [60000.1, 60000.9, 60001.5, 60003.0]

    # ACT
    nights = common.night_count(mjds)

    # ASSERT
    assert nights == 3


def test_provenance_comments_name_the_subject_the_nights_and_the_version() -> None:
    # ARRANGE
    mjds = [60000.1, 60001.2, 60001.4]

    # ACT
    lines = common.provenance_comments("FRAMES", mjds, "0.20.0", "2026-11-22")

    # ASSERT
    assert lines == ["# MEASURED 2026-11-22 WITH soxspipe 0.20.0 FROM 3 FRAMES ON 2 NIGHTS (MJD 60000 TO 60001)"]


def test_arm_block_indents_the_comments_and_the_body_under_the_arm() -> None:
    # ARRANGE
    comments, body = ["# NOTE"], ["key:", "  inner: 1"]

    # ACT
    lines = common.arm_block("VIS", comments, body)

    # ASSERT
    assert lines == ["VIS:", "  # NOTE", "  key:", "    inner: 1"]


def test_join_blocks_ends_the_text_with_a_newline() -> None:
    # ARRANGE
    lines = ["VIS:", "  key: 1"]

    # ACT
    text = common.join_blocks(lines)

    # ASSERT
    assert text == "VIS:\n  key: 1\n"


def test_soxspipe_version_reads_the_generated_version_file(tmp_path: Path) -> None:
    # ARRANGE
    versionFile = tmp_path / "__version__.py"
    versionFile.write_text('__version__ = "1.2.3"\n')

    # ACT
    version = common.soxspipe_version(versionFile)

    # ASSERT
    assert version == "1.2.3"


def test_soxspipe_version_is_unknown_when_the_file_is_missing_or_has_no_version(tmp_path: Path) -> None:
    # ARRANGE
    emptyFile = tmp_path / "__version__.py"
    emptyFile.write_text("# NOTHING\n")

    # ACT
    missing, empty = common.soxspipe_version(tmp_path / "missing.py"), common.soxspipe_version(emptyFile)

    # ASSERT
    assert (missing, empty) == ("unknown", "unknown")


def test_detector_parameters_and_keyword_aliases_are_read_once() -> None:
    # ARRANGE / ACT
    first, second = common.detector_parameters(), common.detector_parameters()
    aliases = common.keyword_aliases()

    # ASSERT
    assert first is second
    assert aliases is common.keyword_aliases()
    assert {"MJDOBS", "EXPTIME", "SEQ_EXPTIME", "ABSROT", "PAR_ANG_START"} <= set(aliases)


def test_today_iso_is_an_iso_date() -> None:
    # ARRANGE / ACT
    today = common.today_iso()

    # ASSERT
    assert re.fullmatch(r"\d{4}-\d{2}-\d{2}", today)


def test_ensure_repo_on_path_puts_the_checkout_on_the_import_path(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    monkeypatch.setattr(sys, "path", [entry for entry in sys.path if entry != str(common.REPO_ROOT)])

    # ACT
    common.ensure_repo_on_path()

    # ASSERT
    assert sys.path[0] == str(common.REPO_ROOT)


def test_read_header_returns_the_primary_header(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "a.fits"
    fits.PrimaryHDU(header=_header(ESO_SEQ_ARM="VIS")).writeto(path)

    # ACT
    header = common.read_header(path)

    # ASSERT
    assert header["ESO SEQ ARM"] == "VIS"


def test_read_header_raises_a_tool_error_for_a_file_that_is_not_fits(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "not_fits.fits"
    path.write_text("plain text")

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match="cannot read"):
        common.read_header(path)


def test_run_main_returns_the_function_status() -> None:
    # ARRANGE / ACT
    status = common.run_main(lambda argv: 0, [])

    # ASSERT
    assert status == 0


def test_run_main_reports_a_tool_error_on_stderr_and_returns_two(capsys: pytest.CaptureFixture[str]) -> None:
    # ARRANGE
    def failing(argv: list[str]) -> int:
        raise common.ToolError("no file matches x")

    # ACT
    status = common.run_main(failing, [])

    # ASSERT
    assert status == 2
    assert "no file matches x" in capsys.readouterr().err


def test_spectral_binning_raises_when_the_card_is_not_a_whole_number(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "a.fits"
    header = _header(ESO_DET_BINX="two")

    # ACT / ASSERT
    with pytest.raises(common.ToolError, match=re.escape(f"{path} has a ESO DET BINX card that is not a whole number")):
        common.spectral_binning(header, "VIS", path)


@pytest.mark.parametrize(("axis", "card"), [("x", "ESO_DET_BINX"), ("y", "ESO_DET_BINY")])
def test_binning_along_reads_the_card_of_the_axis(axis: str, card: str) -> None:
    # ARRANGE
    header = _header(ESO_DET_BINX=3, ESO_DET_BINY=5)

    # ACT
    binning = common.binning_along(header, axis)

    # ASSERT
    assert binning == (3 if card.endswith("X") else 5)


def test_the_cached_yaml_reads_are_read_only_all_the_way_down() -> None:
    # ARRANGE
    parameters, aliases = common.detector_parameters(), common.keyword_aliases()

    # ACT / ASSERT
    with pytest.raises(TypeError):
        parameters["VIS"] = {}
    with pytest.raises(TypeError):
        parameters["VIS"]["binning"] = [2, 2]
    with pytest.raises(TypeError):
        aliases["MJDOBS"] = "X"
    assert parameters["VIS"]["binning"] == (1, 1)


def test_frozen_mapping_is_a_read_only_copy() -> None:
    # ARRANGE
    original = {1: 0.3}

    # ACT
    copy = common.frozen_mapping(original)
    original[2] = 0.4

    # ASSERT
    assert dict(copy) == {1: 0.3}
    with pytest.raises(TypeError):
        copy[3] = 1.0
