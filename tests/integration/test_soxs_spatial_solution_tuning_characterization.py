"""Characterization of `soxs_spatial_solution.produce_product`: the tuning branch and the failure points.

Split from `test_soxs_spatial_solution_reduction_characterization.py` to keep
both files under 800 lines; the builders live there. The failure-point tests
pin where a missing multi-pinhole frame, a missing order table and a missing
dispersion table raise, because a split could move any of them.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_spatial_solution
from tests.integration.test_mapping_recipe_orchestration import _route
from tests.integration.test_soxs_spatial_solution_reduction_characterization import (
    _patch_reduction,
    _vis_recipe,
)

pytestmark = pytest.mark.integration

# THE TUNING GRID: FOUR ORDER AND WAVELENGTH DEGREES EACH, THREE SLIT DEGREES.
TUNING_PERMUTATIONS = 4 * 4 * 4 * 4 * 3 * 3

SPATIAL_MODULE = import_module("soxspipe.recipes.soxs_spatial_solution")


def _tuning_recipe(
    log: Any,
    tmpPath: Path,
    monkeypatch: pytest.MonkeyPatch,
    **options: Any,
) -> tuple[soxs_spatial_solution, list[str], dict[str, Any], dict[str, Any]]:
    """Return a recipe configured for the `tune-pipeline` branch, with tuning stubbed."""
    recipe = _vis_recipe(log, tmpPath, settingOverrides={"tune-pipeline": True}, **options)
    calls, captured = _patch_reduction(recipe, monkeypatch, mapPath=tmpPath / "MAP.fits", mapImagePath=None)
    tuningArguments: dict[str, Any] = {}

    def fake_fmultiprocess(**kwargs: object) -> list[object]:
        tuningArguments.update(kwargs)
        calls.append("fmultiprocess")
        return []

    monkeypatch.setattr(import_module("fundamentals"), "fmultiprocess", fake_fmultiprocess)
    return recipe, calls, captured, tuningArguments


def _in_fresh_directory(tmpPath: Path, monkeypatch: pytest.MonkeyPatch, name: str) -> Path:
    """Run the test from its own directory, because tuning removes `residuals.txt` relative to it."""
    directory = tmpPath / name
    directory.mkdir()
    monkeypatch.chdir(directory)
    return directory


# ---------------------------------------------------------------------------
# THE PARAMETER TUNING BRANCH
# ---------------------------------------------------------------------------


def test_tuning_fits_the_line_list_once_then_sweeps_the_polynomial_grid(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """One map without a 2D image gives the line list, then 2304 permutations are tuned."""
    # ARRANGE
    directory = _in_fresh_directory(tmp_path, monkeypatch, "tuning-run")
    (directory / "residuals.txt").write_text("stale residuals\n")
    recipe, calls, captured, tuningArguments = _tuning_recipe(log, tmp_path, monkeypatch, debug=True)

    # ACT
    returned = recipe.produce_product()

    # ASSERT
    assert returned == (None, None, None)
    assert calls == ["detrend", "keywords", "create_dispersion_map", "map_get", "fmultiprocess"]
    assert not (directory / "residuals.txt").exists()

    (mapArguments,) = captured["maps"]
    assert set(mapArguments) == {
        "log",
        "settings",
        "recipeSettings",
        "pinholeFrame",
        "firstGuessMap",
        "orderTable",
        "qcTable",
        "productsTable",
        "sofName",
        "create2DMap",
        "startNightDate",
        "debug",
    }
    assert mapArguments["create2DMap"] is False
    assert mapArguments["debug"] is True
    assert mapArguments["startNightDate"] == "2024-01-02"

    assert tuningArguments["function"] is SPATIAL_MODULE.parameterTuning
    assert len(tuningArguments["inputArray"]) == TUNING_PERMUTATIONS
    assert tuningArguments["inputArray"][0] == (2, 2, 2, 2, 1, 1)
    assert tuningArguments["inputArray"][-1] == (5, 5, 5, 5, 3, 3)
    assert tuningArguments["poolSize"] is False
    assert tuningArguments["timeout"] == 360000
    assert tuningArguments["turnOffMP"] is True
    assert tuningArguments["mute"] is True
    assert tuningArguments["progressBar"] is True
    assert tuningArguments["log"] is recipe.log
    assert tuningArguments["recipeSettings"] is recipe.recipeSettings
    assert tuningArguments["settings"] is recipe.settings
    assert tuningArguments["multiPinholeFrame"] is recipe.multiPinholeFrame
    assert tuningArguments["disp_map_table"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert tuningArguments["order_table"] == str(tmp_path / "ORDER_TAB_VIS.fits")
    assert tuningArguments["qc"] is recipe.qc
    assert tuningArguments["products"] is recipe.products
    assert tuningArguments["sofName"] == "synthetic-spatial-solution"
    assert tuningArguments["lineDetectionTable"] == "LINE-DETECTION-TABLE"


def test_tuning_skips_the_overrides_the_quicklook_and_the_tables(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`polyOrders` and the debug switch come after the tuning return, so neither applies."""
    # ARRANGE
    _in_fresh_directory(tmp_path, monkeypatch, "tuning-untouched")
    recipe, calls, _, _ = _tuning_recipe(log, tmp_path, monkeypatch, polyOrders=345435, debug=True)
    originalProducts = recipe.products.copy(deep=True)
    originalQc = recipe.qc.copy(deep=True)

    # ACT
    recipe.produce_product()

    # ASSERT
    pd.testing.assert_frame_equal(recipe.products, originalProducts)
    pd.testing.assert_frame_equal(recipe.qc, originalQc)
    assert not {"quicklook", "report", "clean_up"} & set(calls)
    assert recipe.polyOrders == 345435
    assert recipe.recipeSettings == {"use_flat": True}
    assert recipe.create2DMap is True
    assert recipe.slit_arc is None


def test_tuning_without_a_stale_residuals_file_records_the_failed_removal(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A missing `residuals.txt` is swallowed and logged rather than raised."""
    # ARRANGE
    _in_fresh_directory(tmp_path, monkeypatch, "tuning-clean")
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert debugged[0] == "starting the ``produce_product`` method"
    assert debugged[1].startswith("produce_product: `os.remove('residuals.txt')` failed, continuing: ")
    assert len(debugged) == 2


def test_tuning_announces_itself_on_standard_output(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The tuning banner is a bare `print`, not a logger call."""
    # ARRANGE
    _in_fresh_directory(tmp_path, monkeypatch, "tuning-banner")
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert "\n\nTUNING SOXSPIPE\n\n" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# FAILURE POINTS A SPLIT COULD MOVE
# ---------------------------------------------------------------------------


def test_a_missing_multi_pinhole_frame_fails_on_its_date_before_any_detrend(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The frame placeholder is `False`, so reading its header raises `AttributeError`."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    del recipe.inputFrames.filePathsByFilters[_route(DPR_TYPE="WAVE,LAMP", DPR_TECH="ECHELLE,MULTI-PINHOLE")]
    calls, _ = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT / ASSERT
    with pytest.raises(AttributeError, match="header"):
        recipe.produce_product()

    assert calls == []
    assert not hasattr(recipe, "dateObs")


def test_a_missing_order_table_fails_on_its_index_before_any_detrend(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The order table is indexed `[0]`, after the slit arc is set and before the detrend."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    del recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="ORDER_TAB_VIS")]
    calls, _ = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT / ASSERT
    with pytest.raises(IndexError):
        recipe.produce_product()

    assert calls == []
    assert recipe.slit_arc is None
    assert recipe.dateObs is not None


def test_a_missing_dispersion_table_fails_at_the_map_after_every_override(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The table name is bound only inside its loop, so it fails unbound, at the map call."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=345435, debug=True)
    del recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")]
    calls, _ = _patch_reduction(recipe, monkeypatch, mapPath=tmp_path / "MAP.fits", mapImagePath=None)

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError, match="disp_map_table"):
        recipe.produce_product()

    assert calls == ["detrend", "keywords"]
    assert recipe.recipeSettings["slit-deg"] == [3, 5]
    assert recipe.create2DMap is False


def test_a_missing_dispersion_table_fails_after_tuning_clears_the_residuals(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """On the tuning branch the unbound table is read only after `residuals.txt` is removed."""
    # ARRANGE
    directory = _in_fresh_directory(tmp_path, monkeypatch, "tuning-no-table")
    (directory / "residuals.txt").write_text("stale residuals\n")
    recipe, calls, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)
    del recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")]

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError, match="disp_map_table"):
        recipe.produce_product()

    assert calls == ["detrend", "keywords"]
    assert not (directory / "residuals.txt").exists()
