"""The `filepath` of a predicted merged spectrum is the path the recipe writes it to (DY-2343)."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils import toolkit
from soxspipe.commonutils.stare_mode import predicted_nod_product_path
from tests.factories import raw_group_table, workspace_organiser

pytestmark = pytest.mark.unit

SOF_MAP_PATH = Path(soxspipe.__file__).parent / "resources" / "soxs_sof_map.yaml"
NIGHT = "2024-01-01"
NOD_STEM = "20240102T030405_VIS_1_NOD_OBJ_SLIT1_0_5_5_0S_SOXS_SN2024abc"
MERGED_SUFFIX = "_EXTRACTED_MERGED.fits"


def _products(filterName):
    """Return the predicted-product definitions of one filter of the real SOF map."""
    with open(SOF_MAP_PATH) as stream:
        return yaml.safe_load(stream)[filterName]["products"]


def _predicted_rows(tmp_path, log, filterName, recipe, *, sofName=None, settings=None):
    """Run `predict_product_frames` on one raw group and return the `product_frames` rows."""
    organiser = workspace_organiser(tmp_path, log=log)
    organiser.conn = sqlite3.connect(":memory:")
    if settings is not None:
        organiser.settings = settings
    rawGroups = raw_group_table().assign(recipe=recipe, sof=sofName or f"{NOD_STEM}.sof")

    organiser.predict_product_frames(_products(filterName), rawGroups, recipe)

    return pd.read_sql("SELECT * FROM product_frames", organiser.conn)


@pytest.mark.parametrize(
    ("filterName", "recipe", "recipeName", "stem"),
    [
        ("science_stare", "stare_obj", "soxs-stare", "20240102T030405_VIS_1_STARE_OBJ_SLIT1_0_SOXS_SN2024abc"),
        ("science_nodding", "nod_obj", "soxs-nod", NOD_STEM),
        ("science_offset", "offset_obj", "soxs-offset", "20240102T030405_VIS_1_OFFSET_OBJ_SLIT1_0_SOXS_SN2024abc"),
    ],
)
def test_merged_spectrum_filepath_is_in_the_directory_the_recipe_writes_to(
    tmp_path, log, filterName, recipe, recipeName, stem
) -> None:
    # ARRANGE
    _, productDir = toolkit.utility_setup(log, {"workspace-root-dir": str(tmp_path / "out")}, recipeName, NIGHT)
    writtenDirName = Path(productDir).name

    # ACT
    rows = _predicted_rows(tmp_path, log, filterName, recipe, sofName=f"{stem}.sof")

    # ASSERT
    assert rows["filepath"].tolist() == [f"./reduced/{NIGHT}/{writtenDirName}/{stem}{MERGED_SUFFIX}"]
    assert rows["file"].tolist() == [f"{stem}{MERGED_SUFFIX}"]


def test_the_stare_mode_path_of_a_nod_sof_is_the_one_the_skip_check_predicts(tmp_path) -> None:
    # GUARD: THE ORGANISER TEST BELOW TAKES THE WRITTEN NAME FROM `predicted_nod_product_path`, SO PIN THAT IT IS THE
    # STARE RECIPE'S OWN PATH
    settings = {"workspace-root-dir": str(tmp_path / "out"), "soxs-nod": {"reduce-as-stare": True}}

    # ACT
    predictedPath = predicted_nod_product_path(settings, f"{NOD_STEM}.sof")

    # ASSERT
    stareRecipePath, _ = toolkit.predict_product_path(
        f"{NOD_STEM}.sof", "soxs-stare", "_STAREMODE", outputDir=settings["workspace-root-dir"]
    )
    assert predictedPath == stareRecipePath


def test_nod_sof_reduced_as_stare_records_the_stare_mode_file_that_is_written(tmp_path, log) -> None:
    # ARRANGE
    settings = {"workspace-root-dir": str(tmp_path / "out"), "soxs-nod": {"reduce-as-stare": True}}
    writtenPath = predicted_nod_product_path(settings, f"{NOD_STEM}.sof")

    # ACT
    rows = _predicted_rows(tmp_path, log, "science_nodding", "nod_obj", settings=settings)

    # ASSERT
    writtenFile = Path(writtenPath).name
    assert writtenFile == f"{NOD_STEM}_STAREMODE{MERGED_SUFFIX}"
    assert rows["file"].tolist() == [writtenFile]
    assert rows["filepath"].tolist() == [f"./reduced/{NIGHT}/soxs-stare/{writtenFile}"]


def test_nod_sof_warns_and_keeps_nod_names_when_the_session_settings_are_not_loaded(tmp_path, log) -> None:
    # ACT
    rows = _predicted_rows(tmp_path, log, "science_nodding", "nod_obj")

    # ASSERT
    assert rows["filepath"].tolist() == [f"./reduced/{NIGHT}/soxs-nod/{NOD_STEM}{MERGED_SUFFIX}"]
    assert [level for level, message in log.messages if "settings are not loaded" in message] == ["warning"]


def test_stare_sof_does_not_warn_when_the_session_settings_are_not_loaded(tmp_path, log) -> None:
    # ACT
    _predicted_rows(tmp_path, log, "science_stare", "stare_obj")

    # ASSERT
    assert not [level for level, message in log.messages if level == "warning"]


def test_nod_sof_reduced_as_stare_keeps_its_other_product_columns(tmp_path, log) -> None:
    # ARRANGE
    asNod = _predicted_rows(tmp_path, log, "science_nodding", "nod_obj")

    # ACT
    asStare = _predicted_rows(
        tmp_path, log, "science_nodding", "nod_obj", settings={"soxs-nod": {"reduce-as-stare": True}}
    )

    # ASSERT
    unchanged = [column for column in asNod.columns if column not in ("file", "filepath")]
    pd.testing.assert_frame_equal(asStare[unchanged], asNod[unchanged])


def test_nod_sof_stays_in_nod_directory_when_the_settings_do_not_ask_for_stare_mode(tmp_path, log) -> None:
    # ACT
    rows = _predicted_rows(
        tmp_path, log, "science_nodding", "nod_obj", settings={"soxs-nod": {"reduce-as-stare": False}}
    )

    # ASSERT
    assert rows["filepath"].tolist() == [f"./reduced/{NIGHT}/soxs-nod/{NOD_STEM}{MERGED_SUFFIX}"]


@pytest.mark.parametrize(
    ("filterName", "recipe", "directory"),
    [
        ("flux_standard_stare", "stare_std", "soxs-stare-std"),
        ("flux_standard_nodding", "nod_std", "soxs-nod-std"),
        ("flux_standard_offset", "offset_std", "soxs-offset-std"),
    ],
)
def test_standard_star_product_rows_are_unchanged(tmp_path, log, filterName, recipe, directory) -> None:
    # ACT
    rows = _predicted_rows(tmp_path, log, filterName, recipe, sofName=f"{NOD_STEM}.sof")

    # ASSERT
    assert set(rows["filepath"]) == {
        f"./reduced/{NIGHT}/{directory}/{NOD_STEM}{MERGED_SUFFIX}",
        f"./reduced/{NIGHT}/{directory}/{NOD_STEM}_RESP.fits",
    }
