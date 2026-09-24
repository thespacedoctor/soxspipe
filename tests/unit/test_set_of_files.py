"""Characterization tests pinning the current behaviour of `set_of_files.py`.

These tests describe what the module does today, including behaviour that
looks like a defect (DY-73, DY-74, DY-75, and others noted inline). They must
all pass unchanged against the current implementation -- they are not
RED/GREEN tests, they exist to freeze behaviour ahead of a refactor.
"""

from __future__ import annotations

from pathlib import Path

import pytest
from astropy.io import fits

from soxspipe.commonutils.set_of_files import (
    set_of_files,
)
from tests.factories import prepared_fits, raw_fits, sof_file

pytestmark = pytest.mark.unit


def _settings(
    tmp_path: Path,
    *,
    default: list[str] | None = None,
    verbose: list[str] | None = None,
    nodding_extras: list[str] | None = None,
) -> dict[str, object]:
    """Return a minimal settings mapping with overridable summary-keys lists."""
    from tests.factories import pipeline_settings

    return pipeline_settings(
        tmp_path,
        overrides={
            "summary-keys": {
                "default": [] if default is None else default,
                "verbose": [] if verbose is None else verbose,
                "nodding_extras": [] if nodding_extras is None else nodding_extras,
            }
        },
    )


# ---------------------------------------------------------------------------
# set_of_files.__init__
# ---------------------------------------------------------------------------


def test_init_defaults_when_only_log_and_settings_are_given(
    tmp_path: Path,
    log: object,
) -> None:
    settings = _settings(tmp_path, default=["DPR_TYPE"], verbose=["SEQ_ARM"])

    sof = set_of_files(log=log, settings=settings)

    assert sof.inputFrames == []
    assert sof.verbose is True
    assert sof.recipeName is False
    assert sof.ext == 0
    # VERBOSE DEFAULTS TO TRUE, SO THE VERBOSE KEYS ARE SELECTED
    assert sof.keys == ["ESO SEQ ARM", "file"]


def test_init_default_settings_of_false_raises_type_error(log: object) -> None:
    with pytest.raises(TypeError, match="'bool' object is not subscriptable"):
        set_of_files(log=log)


def test_init_verbose_true_selects_the_verbose_summary_keys(
    tmp_path: Path, log: object
) -> None:
    settings = _settings(tmp_path, default=["DPR_TYPE"], verbose=["SEQ_ARM"])

    sof = set_of_files(log=log, settings=settings, verbose=True)

    assert sof.keys == ["ESO SEQ ARM", "file"]


def test_init_verbose_false_selects_the_default_summary_keys(
    tmp_path: Path, log: object
) -> None:
    settings = _settings(tmp_path, default=["DPR_TYPE"], verbose=["SEQ_ARM"])

    sof = set_of_files(log=log, settings=settings, verbose=False)

    assert sof.keys == ["ESO DPR TYPE", "file"]


def test_init_nodding_recipe_does_not_mutate_the_shared_settings_list(
    tmp_path: Path, log: object
) -> None:
    """DY-73: `keys += nodding_extras` used to mutate the settings dict's own list.

    Every `set_of_files` instantiation with `recipeName="soxs-nod"` used to
    append the nodding-extras keys onto `settings["summary-keys"]["verbose"]`
    again, in place, so the settings object grew unboundedly across repeated
    use. `set_of_files.__init__` now builds its per-instance key list from a
    new list, so the settings dict's own list is left untouched no matter how
    many times `set_of_files` is instantiated.
    """
    settings = _settings(tmp_path, verbose=["DPR_TYPE"], nodding_extras=["SEQ_ARM"])
    verboseKeysList = settings["summary-keys"]["verbose"]

    for _ in range(3):
        set_of_files(
            log=log, settings=settings, verbose=True, recipeName="soxs-nod"
        )

    assert settings["summary-keys"]["verbose"] is verboseKeysList
    assert verboseKeysList == ["DPR_TYPE"]


def test_init_keys_is_a_new_list_not_the_settings_list(
    tmp_path: Path, log: object
) -> None:
    settings = _settings(tmp_path, verbose=["DPR_TYPE"])

    sof = set_of_files(log=log, settings=settings, verbose=True)

    assert sof.keys is not settings["summary-keys"]["verbose"]


def test_init_expands_a_leading_tilde_in_input_frames_string(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    fakeHome = str(tmp_path / "fakehome")
    monkeypatch.setenv("HOME", fakeHome)
    settings = _settings(tmp_path)

    sof = set_of_files(
        log=log, settings=settings, inputFrames="~/x/y.fits", verbose=False
    )

    # NOTE THE DOUBLE SLASH: `home + "/" + inputFrames[1:]` KEEPS THE LEADING
    # "/" FROM THE ORIGINAL "~/x/y.fits" STRING
    assert sof.inputFrames == f"{fakeHome}//x/y.fits"


def test_init_current_session_is_none_without_a_workspace_session(
    tmp_path: Path, log: object
) -> None:
    settings = _settings(tmp_path)

    sof = set_of_files(log=log, settings=settings, verbose=False)

    assert sof.currentSession is None


# ---------------------------------------------------------------------------
# set_of_files.get()
# ---------------------------------------------------------------------------


def test_get_expands_a_leading_tilde_set_after_construction(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    fakeHome = tmp_path / "fakehome"
    fakeHome.mkdir()
    monkeypatch.setenv("HOME", str(fakeHome))
    framesDir = fakeHome / "frames_dir"
    framesDir.mkdir()
    raw_fits(framesDir / "frame.fits")
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    sof = set_of_files(log=log, settings=settings, verbose=False)
    # __init__'s OWN EXPANSION ALREADY RAN AGAINST THE DEFAULT `[]`, SO SET
    # THE TILDE PATH AFTERWARDS TO EXERCISE `get()`'s SEPARATE EXPANSION
    sof.inputFrames = "~/frames_dir"

    collection, _supplementary = sof.get()

    assert list(collection.summary["file"]) == ["frame.fits"]


def test_get_directory_ext_missing_keys_fall_back_to_the_primary_header(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    firstPath = prepared_fits(framesDir / "a.fits", seed=1)
    secondPath = prepared_fits(framesDir / "b.fits", seed=2)
    fits.setval(firstPath, "ESO DPR TYPE", value="BIAS")
    fits.setval(secondPath, "ESO DPR TYPE", value="OBJECT")
    settings = _settings(tmp_path, default=["DPR_TYPE", "SEQ_ARM"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(framesDir), verbose=False, ext=1
    ).get()

    assert list(collection.summary["file"]) == ["a.fits", "b.fits"]
    assert list(collection.summary["ESO DPR TYPE"]) == ["BIAS", "OBJECT"]
    assert list(collection.summary["ESO SEQ ARM"]) == ["VIS", "VIS"]


def test_get_directory_ext_all_keys_present_skips_the_primary_fallback(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    prepared_fits(framesDir / "a.fits")
    # NAXIS1/NAXIS2 ARE PRESENT ON THE EXTENSION HEADER ITSELF, SO NO
    # `missingKeys` REMAIN AND THE PRIMARY-HEADER JOIN IS SKIPPED
    settings = _settings(tmp_path, default=["NAXIS1", "NAXIS2"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(framesDir), verbose=False, ext=1
    ).get()

    assert collection.summary.colnames == ["file", "NAXIS1", "NAXIS2"]


def test_get_directory_excludes_hidden_files_from_supplementary(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    raw_fits(framesDir / "a.fits")
    (framesDir / ".VIS_HIDDEN.csv").write_text("secret", encoding="utf-8")
    (framesDir / "VIS_2D_MAP.csv").write_text("x\n1\n", encoding="utf-8")
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    _collection, supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(framesDir), verbose=False
    ).get()

    assert supplementary == {
        "VIS": {"2D_MAP": str(framesDir / "VIS_2D_MAP.csv")}
    }


def test_get_sof_session_rewrites_reduced_paths_to_the_session_directory(
    tmp_path: Path, log: object
) -> None:
    reducedDir = tmp_path / "sessions" / "mysession" / "reduced"
    reducedDir.mkdir(parents=True)
    raw_fits(reducedDir / "a.fits")
    (reducedDir / "VIS_ORDER_LOCATIONS.csv").write_text("x\n1\n", encoding="utf-8")
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        "./reduced/a.fits BIAS_VIS\n./reduced/VIS_ORDER_LOCATIONS.csv\n",
        encoding="utf-8",
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    sof = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    )
    sof.currentSession = "mysession"

    collection, supplementary = sof.get()

    assert list(collection.summary["file"]) == ["a.fits"]
    assert supplementary == {
        "VIS": {
            "ORDER_LOCATIONS": "./sessions/mysession/reduced/VIS_ORDER_LOCATIONS.csv"
        }
    }


def test_get_sof_ignores_comment_and_blank_lines(tmp_path: Path, log: object) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        f"# a comment line\n\n{framePath} BIAS_VIS\n", encoding="utf-8"
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert list(collection.summary["file"]) == ["frame.fits"]


def test_get_sof_ignores_supplementary_lines_of_three_characters_or_fewer(
    tmp_path: Path, log: object
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    inputPath = tmp_path / "input.sof"
    # "abc" IS A THREE-CHARACTER NON-FITS LINE: `len(l) > 3` EXCLUDES IT
    inputPath.write_text(f"{framePath} BIAS_VIS\nabc\n", encoding="utf-8")
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    _collection, supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert supplementary == {}


def test_get_sof_crlf_supplementary_line_keeps_its_trailing_carriage_return(
    tmp_path: Path, log: object
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    (tmp_path / "VIS_ORDER_LOCATIONS.csv").write_text("x\n1\n", encoding="utf-8")
    inputPath = tmp_path / "input.sof"
    # WRITTEN IN BINARY MODE SO THE CRLF SURVIVES: `codecs.open(..., mode="r")`
    # DOES NOT TRANSLATE NEWLINES
    inputPath.write_bytes(
        f"{framePath} BIAS_VIS\r\nVIS_ORDER_LOCATIONS.csv\r\n".encode()
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    with pytest.raises(FileNotFoundError) as raisedError:
        set_of_files(
            log=log, settings=settings, inputFrames=str(inputPath), verbose=False
        ).get()

    assert str(raisedError.value) == (
        "the input file `VIS_ORDER_LOCATIONS.csv\r` does not appear to exist"
    )


def test_get_sof_crlf_fits_only_line_still_resolves(
    tmp_path: Path, log: object
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    inputPath = tmp_path / "input.sof"
    inputPath.write_bytes(f"{framePath} BIAS_VIS\r\n".encode())
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert list(collection.summary["file"]) == ["frame.fits"]


def test_get_sof_supplementary_files_pollute_the_collection_but_not_the_summary(
    tmp_path: Path, log: object
) -> None:
    """DY-74: `fitsFiles.extend(supplementaryFilepaths)` mutates `fitsFiles`.

    The SOF branch feeds the (now supplementary-polluted) `fitsFiles` list
    straight into `ImageFileCollection(filenames=fitsFiles, ...)`, so a
    supplementary file that ccdproc cannot parse as FITS ends up in
    `collection.files` even though it never appears in `collection.summary`.
    This is a known defect (DY-74); it is pinned here, not fixed.
    """
    framePath = raw_fits(tmp_path / "a.fits")
    supplementaryPath = tmp_path / "VIS_ORDER_LOCATIONS.csv"
    supplementaryPath.write_text("x\n1\n", encoding="utf-8")
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        f"{framePath} BIAS_VIS\n{supplementaryPath}\n", encoding="utf-8"
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert collection.files == ["a.fits", "VIS_ORDER_LOCATIONS.csv"]
    assert list(collection.summary["file"]) == ["a.fits"]


def test_get_sof_supplementary_in_a_different_directory_forces_no_location(
    tmp_path: Path, log: object
) -> None:
    """DY-74, continued: a supplementary file elsewhere ruins `location`.

    Because the supplementary path is folded into `fitsFiles` before the
    common-location check, a supplementary file living outside the frame
    directory makes `len(set(locations)) != 1`, so `location` becomes `None`
    and even the genuine FITS frame is reported with its full path.
    """
    framesDir = tmp_path / "frames"
    suppDir = tmp_path / "supp"
    framesDir.mkdir()
    suppDir.mkdir()
    framePath = raw_fits(framesDir / "a.fits")
    supplementaryPath = suppDir / "VIS_ORDER_LOCATIONS.csv"
    supplementaryPath.write_text("x\n1\n", encoding="utf-8")
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        f"{framePath.resolve()} BIAS_VIS\n{supplementaryPath.resolve()}\n",
        encoding="utf-8",
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert collection.location == ""
    assert list(collection.summary["file"]) == [str(framePath.resolve())]


def test_get_sof_ext_all_keys_present_skips_the_primary_fallback(
    tmp_path: Path, log: object
) -> None:
    preparedPath = prepared_fits(tmp_path / "a.fits")
    inputPath = sof_file(tmp_path / "input.sof", [(preparedPath, "BIAS_VIS")])
    settings = _settings(tmp_path, default=["NAXIS1", "NAXIS2"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False, ext=1
    ).get()

    assert collection.summary.colnames == ["file", "NAXIS1", "NAXIS2"]


def test_get_sof_no_common_location_uses_full_paths(
    tmp_path: Path, log: object
) -> None:
    firstDir = tmp_path / "dirA"
    secondDir = tmp_path / "dirB"
    firstDir.mkdir()
    secondDir.mkdir()
    firstPath = raw_fits(firstDir / "first.fits", seed=1)
    secondPath = raw_fits(secondDir / "second.fits", seed=2)
    inputPath = sof_file(
        tmp_path / "input.sof",
        [(secondPath.resolve(), "BIAS_VIS"), (firstPath.resolve(), "BIAS_VIS")],
    )
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(inputPath), verbose=False
    ).get()

    assert list(collection.summary["file"]) == [
        str(secondPath.resolve()),
        str(firstPath.resolve()),
    ]


def test_get_list_drops_dot_prefixed_supplementary_paths(
    tmp_path: Path, log: object
) -> None:
    """DY-75: a relative supplementary path starting with "." is dropped.

    The list-input branch filters supplementary candidates with
    `f[0] != "."`, so a path like "./VIS_DISP_MAP.csv" is excluded from
    `supplementaryFilepaths` entirely, rather than being resolved relative to
    the current directory. This is a known defect (DY-75); it is pinned
    here, not fixed.
    """
    framePath = raw_fits(tmp_path / "frame.fits")
    (tmp_path / "VIS_DISP_MAP.csv").write_text("x\n1\n", encoding="utf-8")
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    sof = set_of_files(log=log, settings=settings, verbose=False)
    sof.inputFrames = [str(framePath), "./VIS_DISP_MAP.csv"]

    _collection, supplementary = sof.get()

    assert supplementary == {}


def test_get_list_no_common_location_uses_full_paths(
    tmp_path: Path, log: object
) -> None:
    firstDir = tmp_path / "dirA"
    secondDir = tmp_path / "dirB"
    firstDir.mkdir()
    secondDir.mkdir()
    firstPath = raw_fits(firstDir / "first.fits", seed=1)
    secondPath = raw_fits(secondDir / "second.fits", seed=2)
    settings = _settings(tmp_path, default=["DPR_TYPE"])

    collection, _supplementary = set_of_files(
        log=log,
        settings=settings,
        inputFrames=[str(secondPath), str(firstPath)],
        verbose=False,
    ).get()

    assert list(collection.summary["file"]) == [str(secondPath), str(firstPath)]


def test_get_list_ext_all_keys_present_skips_the_primary_fallback(
    tmp_path: Path, log: object
) -> None:
    preparedPath = prepared_fits(tmp_path / "a.fits")
    settings = _settings(tmp_path, default=["NAXIS1", "NAXIS2"])

    collection, _supplementary = set_of_files(
        log=log,
        settings=settings,
        inputFrames=[str(preparedPath)],
        verbose=False,
        ext=1,
    ).get()

    assert "filename" in collection.summary.colnames
    assert "NAXIS1" in collection.summary.colnames
    assert "NAXIS2" in collection.summary.colnames


def test_get_list_adds_filename_column_and_prepends_key_on_every_call(
    tmp_path: Path, log: object
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(
        log=log, settings=settings, inputFrames=[str(framePath)], verbose=False
    )

    sof.get()
    assert sof.keys == ["filename", "ESO DPR TYPE", "file"]

    # A SECOND `get()` CALL PREPENDS "filename" AGAIN, GIVING A DUPLICATE ENTRY
    sof.get()
    assert sof.keys == ["filename", "filename", "ESO DPR TYPE", "file"]


# ---------------------------------------------------------------------------
# _generate_sof_file_from_directory
# ---------------------------------------------------------------------------


def test_generate_sof_expands_a_leading_tilde_in_directory_and_sof_path(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    fakeHome = tmp_path / "fakehome"
    fakeHome.mkdir()
    monkeypatch.setenv("HOME", str(fakeHome))
    framesDir = fakeHome / "frames_dir"
    framesDir.mkdir()
    framePath = raw_fits(framesDir / "a.fits")
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(log=log, settings=settings, verbose=False)

    result = sof._generate_sof_file_from_directory("~/frames_dir", "~/out/frames.sof")

    assert result == str(fakeHome / "out" / "frames.sof")
    assert Path(result).read_text(encoding="utf-8") == f"{framePath.resolve()} BIAS_VIS\n"


def test_generate_sof_skips_non_fits_and_subdirectory_entries(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    framePath = raw_fits(framesDir / "a.fits")
    (framesDir / "notes.txt").write_text("ignore me", encoding="utf-8")
    (framesDir / "subdir").mkdir()
    outputPath = tmp_path / "out.sof"
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(log=log, settings=settings, verbose=False)

    sof._generate_sof_file_from_directory(str(framesDir), str(outputPath))

    assert outputPath.read_text(encoding="utf-8") == f"{framePath.resolve()} BIAS_VIS\n"


def test_generate_sof_appends_an_int_truncated_binning_suffix(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    framePath = raw_fits(framesDir / "a.fits")
    fits.setval(framePath, "ESO DET BINX", value=2.7)
    fits.setval(framePath, "ESO DET BINY", value=3.9)
    outputPath = tmp_path / "out.sof"
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(log=log, settings=settings, verbose=False)

    sof._generate_sof_file_from_directory(str(framesDir), str(outputPath))

    # int(2.7) == 2 AND int(3.9) == 3: TRUNCATION, NOT ROUNDING
    assert outputPath.read_text(encoding="utf-8") == (
        f"{framePath.resolve()} BIAS_VIS_2x3\n"
    )


def test_generate_sof_reuses_an_already_existing_output_directory(
    tmp_path: Path, log: object
) -> None:
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    raw_fits(framesDir / "a.fits")
    outputDir = tmp_path / "out"
    outputDir.mkdir()
    outputPath = outputDir / "frames.sof"
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(log=log, settings=settings, verbose=False)

    result = sof._generate_sof_file_from_directory(str(framesDir), str(outputPath))

    assert result == str(outputPath)
    assert outputPath.exists()


# ---------------------------------------------------------------------------
# create_supplementary_file_dictionary
# ---------------------------------------------------------------------------


def test_create_supplementary_dictionary_matches_several_arms_and_kinds(
    tmp_path: Path, log: object
) -> None:
    settings = _settings(tmp_path, default=["DPR_TYPE"])
    sof = set_of_files(log=log, settings=settings, verbose=False)

    result = sof.create_supplementary_file_dictionary(
        [
            # MIXED CASE, AND A FILENAME MATCHING TWO ARMS AND TWO KINDS AT ONCE
            "nir_VIS_disp_map_order_locations.csv",
            "vis_2D_MAP.CSV",
        ]
    )

    assert result == {
        "NIR": {
            "DISP_MAP": "nir_VIS_disp_map_order_locations.csv",
            "ORDER_LOCATIONS": "nir_VIS_disp_map_order_locations.csv",
        },
        "VIS": {
            "DISP_MAP": "nir_VIS_disp_map_order_locations.csv",
            "ORDER_LOCATIONS": "nir_VIS_disp_map_order_locations.csv",
            "2D_MAP": "vis_2D_MAP.CSV",
        },
    }
