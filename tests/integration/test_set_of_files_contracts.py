"""Synthetic integration tests for accepted set-of-files inputs."""

from __future__ import annotations

from pathlib import Path

import pytest
from astropy.io import fits

from soxspipe.commonutils.set_of_files import set_of_files
from tests.factories import pipeline_settings, prepared_fits, raw_fits, sof_file

pytestmark = pytest.mark.integration


def _settings(tmp_path: Path) -> dict[str, object]:
    return pipeline_settings(
        tmp_path,
        overrides={
            "summary-keys": {
                "default": ["DPR_TYPE", "SEQ_ARM"],
                "verbose": ["DPR_TYPE", "SEQ_ARM"],
                "nodding_extras": [],
            }
        },
    )


def test_directory_input_is_ordered_and_maps_supplementary_files(
    tmp_path: Path,
    log: object,
) -> None:
    framesPath = tmp_path / "frames"
    framesPath.mkdir()
    raw_fits(framesPath / "b.fits", seed=2)
    raw_fits(framesPath / "a.fits", seed=1)
    dispersionPath = framesPath / "VIS_DISP_MAP.csv"
    dispersionPath.write_text("synthetic", encoding="utf-8")

    collection, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(framesPath),
        verbose=False,
    ).get()

    assert list(collection.summary["file"]) == ["a.fits", "b.fits"]
    assert supplementary == {"VIS": {"DISP_MAP": str(dispersionPath)}}


def test_list_input_preserves_caller_order(tmp_path: Path, log: object) -> None:
    firstPath = raw_fits(tmp_path / "first.fits", seed=1)
    secondPath = raw_fits(tmp_path / "second.fits", seed=2)

    collection, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=[str(secondPath), str(firstPath)],
        verbose=False,
    ).get()

    assert list(collection.summary["filename"]) == ["second.fits", "first.fits"]
    assert supplementary == {}


def test_list_input_classifies_each_supplementary_file_type(
    tmp_path: Path,
    log: object,
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    dispersionPath = tmp_path / "VIS_DISP_MAP.csv"
    ordersPath = tmp_path / "VIS_ORDER_CENTRES.csv"
    twoDPath = tmp_path / "VIS_2D_MAP.csv"
    for path in (dispersionPath, ordersPath, twoDPath):
        path.touch()

    _, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=[
            str(framePath),
            str(dispersionPath),
            str(ordersPath),
            str(twoDPath),
        ],
        verbose=False,
    ).get()

    assert supplementary == {
        "VIS": {
            "DISP_MAP": str(dispersionPath),
            "ORDER_LOCATIONS": str(ordersPath),
            "2D_MAP": str(twoDPath),
        }
    }


def test_sof_input_preserves_member_order(tmp_path: Path, log: object) -> None:
    firstPath = raw_fits(tmp_path / "first.fits", seed=1)
    secondPath = raw_fits(tmp_path / "second.fits", seed=2)
    inputPath = sof_file(
        tmp_path / "input.sof",
        [(secondPath, "BIAS_VIS"), (firstPath, "BIAS_VIS")],
    )

    collection, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(inputPath),
        verbose=False,
    ).get()

    assert list(collection.summary["file"]) == ["second.fits", "first.fits"]
    assert supplementary == {}


def test_sof_input_rejects_a_missing_member(tmp_path: Path, log: object) -> None:
    missingPath = tmp_path / "missing.fits"
    inputPath = sof_file(tmp_path / "input.sof", [(missingPath, "BIAS_VIS")])

    with pytest.raises(
        FileNotFoundError,
        match=r"input file `.*missing\.fits` does not appear to exist",
    ):
        set_of_files(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=str(inputPath),
            verbose=False,
        ).get()


def test_sof_input_maps_tagged_supplementary_members(
    tmp_path: Path,
    log: object,
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    ordersPath = tmp_path / "VIS_ORDER_LOCATIONS.csv"
    ordersPath.touch()
    inputPath = sof_file(
        tmp_path / "input.sof",
        [
            (framePath, "BIAS_VIS"),
            (ordersPath, "ORDER_TAB_VIS"),
        ],
    )

    _, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(inputPath),
        verbose=False,
    ).get()

    assert supplementary == {
        "VIS": {"ORDER_LOCATIONS": str(ordersPath)}
    }


def test_sof_input_accepts_an_untagged_supplementary_path(
    tmp_path: Path,
    log: object,
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    dispersionPath = tmp_path / "VIS DISP_MAP.csv"
    dispersionPath.touch()
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        f"{framePath} BIAS_VIS\n{dispersionPath}\n",
        encoding="utf-8",
    )

    _, supplementary = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(inputPath),
        verbose=False,
    ).get()

    assert supplementary == {
        "VIS": {"DISP_MAP": str(dispersionPath)}
    }


def test_sof_input_does_not_redirect_a_missing_untagged_path(
    tmp_path: Path,
    log: object,
) -> None:
    framePath = raw_fits(tmp_path / "frame.fits")
    existingPrefix = tmp_path / "VIS"
    existingPrefix.touch()
    missingPath = tmp_path / "VIS VIS_DISP_MAP"
    inputPath = tmp_path / "input.sof"
    inputPath.write_text(
        f"{framePath} BIAS_VIS\n{missingPath}\n",
        encoding="utf-8",
    )

    with pytest.raises(FileNotFoundError, match=str(missingPath)):
        set_of_files(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=str(inputPath),
            verbose=False,
        ).get()


def test_prepared_extension_uses_primary_header_fallback(
    tmp_path: Path,
    log: object,
) -> None:
    preparedPath = prepared_fits(tmp_path / "prepared.fits")

    collection, _ = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=[str(preparedPath)],
        verbose=False,
        ext=1,
    ).get()

    assert list(collection.summary["ESO DPR TYPE"]) == ["BIAS"]
    assert list(collection.summary["ESO SEQ ARM"]) == ["VIS"]


def test_prepared_list_fallback_preserves_header_and_filename_order(
    tmp_path: Path,
    log: object,
) -> None:
    firstPath = prepared_fits(tmp_path / "a.fits", seed=1)
    secondPath = prepared_fits(tmp_path / "b.fits", seed=2)
    fits.setval(firstPath, "ESO DPR TYPE", value="BIAS")
    fits.setval(secondPath, "ESO DPR TYPE", value="OBJECT")

    collection, _ = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=[str(secondPath), str(firstPath)],
        verbose=False,
        ext=1,
    ).get()

    assert list(collection.summary["filename"]) == ["b.fits", "a.fits"]
    assert list(collection.summary["ESO DPR TYPE"]) == ["OBJECT", "BIAS"]


def test_prepared_sof_fallback_preserves_member_and_header_order(
    tmp_path: Path,
    log: object,
) -> None:
    firstPath = prepared_fits(tmp_path / "a.fits", seed=1)
    secondPath = prepared_fits(tmp_path / "b.fits", seed=2)
    fits.setval(firstPath, "ESO DPR TYPE", value="BIAS")
    fits.setval(secondPath, "ESO DPR TYPE", value="OBJECT")
    inputPath = sof_file(
        tmp_path / "input.sof",
        [(secondPath, "OBJECT_VIS"), (firstPath, "BIAS_VIS")],
    )

    collection, _ = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        inputFrames=str(inputPath),
        verbose=False,
        ext=1,
    ).get()

    assert list(collection.summary["file"]) == ["b.fits", "a.fits"]
    assert list(collection.summary["ESO DPR TYPE"]) == ["OBJECT", "BIAS"]


def test_generate_sof_from_directory_is_sorted_and_uses_header_categories(
    tmp_path: Path,
    log: object,
) -> None:
    framesPath = tmp_path / "frames"
    framesPath.mkdir()
    laterPath = raw_fits(framesPath / "b.fits", seed=2)
    earlierPath = raw_fits(framesPath / "a.fits", seed=1)
    outputPath = tmp_path / "inventory" / "frames.sof"

    result = set_of_files(
        log=log,
        settings=_settings(tmp_path),
        verbose=False,
    )._generate_sof_file_from_directory(str(framesPath), str(outputPath))

    assert result == str(outputPath)
    assert outputPath.read_text(encoding="utf-8").splitlines() == [
        f"{earlierPath.resolve()} BIAS_VIS",
        f"{laterPath.resolve()} BIAS_VIS",
    ]


@pytest.mark.parametrize("invalidInput", [None, 17, ""])
def test_invalid_input_has_a_precise_type_error(
    tmp_path: Path,
    log: object,
    invalidInput: object,
) -> None:
    with pytest.raises(TypeError, match="'inputFrames' should be"):
        set_of_files(
            log=log,
            settings=_settings(tmp_path),
            inputFrames=invalidInput,
            verbose=False,
        ).get()
