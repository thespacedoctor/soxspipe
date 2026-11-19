"""The recipe constructor predicts its product where the recipe writes it, under `workspace-root-dir` (DY-2241)."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy.io import fits

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.recipes import base_recipe
from soxspipe.recipes.base_recipe import RecipeProductExistsError

pytestmark = pytest.mark.unit

SOF_STEM = "20251211T123456_NOD_OBJ"
PRODUCT_NAME = f"{SOF_STEM}_EXTRACTED_MERGED.fits"
NIGHT = "2025-12-10"
SESSION = "S1"


class WorkspaceOrganiser:
    """Report the session of a prepared workspace in the current working directory."""

    def __init__(self, **_: object) -> None:
        return None

    def session_list(self, *, silent: bool) -> tuple[str, list[str]]:
        return SESSION, [SESSION]

    def close(self) -> None:
        return None


@pytest.fixture(name="loggerPaths")
def record_logger_paths(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> list[str]:
    """Run the real path prediction and `utility_setup`, and record where the recipe log would be written."""
    paths: list[str] = []
    monkeypatch.setattr(commonutils, "data_organiser", WorkspaceOrganiser)

    def record_logger_path(receivedLog: Any, productPath: str) -> Any:
        paths.append(productPath)
        return receivedLog

    monkeypatch.setattr(toolkit, "add_recipe_logger", record_logger_path)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    workspace = tmp_path / "workspace"
    workspace.mkdir()
    monkeypatch.chdir(workspace)
    return paths


def _write_fits(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([fits.PrimaryHDU(np.zeros((8, 8), dtype="float32"))]).writeto(path)
    return path


def _product_under(root: Path, layout: str) -> Path:
    return root / layout / NIGHT / "soxs-nod" / PRODUCT_NAME


def _construct(log: Any, settings: dict[str, Any], overwrite: bool = False) -> base_recipe:
    return base_recipe(
        log=log,
        settings={"instrument": "soxs", **settings},
        inputFrames=f"{SOF_STEM}.sof",
        recipeName="soxs-nod",
        overwrite=overwrite,
        turnOffMP=True,
    )


def test_a_complete_product_under_the_output_directory_refuses_the_run(
    log: Any, loggerPaths: list[str], tmp_path: Path
) -> None:
    # ARRANGE
    outputDir = tmp_path / "output"
    _write_fits(_product_under(outputDir, "reduced"))

    # ACT
    with pytest.raises(RecipeProductExistsError):
        _construct(log, {"workspace-root-dir": str(outputDir)})

    # ASSERT
    assert loggerPaths == []


def test_an_error_log_under_the_output_directory_refuses_the_run(
    log: Any, loggerPaths: list[str], tmp_path: Path
) -> None:
    # ARRANGE
    outputDir = tmp_path / "output"
    errorLog = _product_under(outputDir, "reduced").with_name(f"{SOF_STEM}_EXTRACTED_MERGED_ERROR.log")
    errorLog.parent.mkdir(parents=True)
    errorLog.write_text("previous failure\n")

    # ACT / ASSERT
    with pytest.raises(RecipeProductExistsError):
        _construct(log, {"workspace-root-dir": str(outputDir)})


def test_a_product_only_in_the_working_directory_does_not_refuse_a_run_with_an_output_directory(
    log: Any, loggerPaths: list[str], tmp_path: Path
) -> None:
    # ARRANGE
    outputDir = tmp_path / "output"
    _write_fits(_product_under(Path("."), f"sessions/{SESSION}/reduced"))

    # ACT
    recipe = _construct(log, {"workspace-root-dir": str(outputDir)})

    # ASSERT
    assert recipe.productPath == str(_product_under(outputDir, "reduced"))
    assert loggerPaths == [recipe.productPath]


@pytest.mark.parametrize(
    "settings",
    [
        pytest.param({}, id="settings-without-a-root-as-the-reducer-and-the-session-file-give"),
        pytest.param({"workspace-root-dir": "./"}, id="the-default-root"),
    ],
)
def test_a_run_without_an_output_directory_still_looks_in_the_working_directory_session(
    log: Any, loggerPaths: list[str], settings: dict[str, Any]
) -> None:
    # ARRANGE
    _write_fits(_product_under(Path("."), f"sessions/{SESSION}/reduced"))

    # ACT
    with pytest.raises(RecipeProductExistsError):
        _construct(log, settings)

    # ASSERT
    assert loggerPaths == []


def test_a_run_without_an_output_directory_predicts_the_session_product_path_and_log_location(
    log: Any, loggerPaths: list[str]
) -> None:
    # ACT
    recipe = _construct(log, {})

    # ASSERT
    assert recipe.productPath == f"./sessions/{SESSION}/reduced/{NIGHT}/soxs-nod/{PRODUCT_NAME}"
    assert loggerPaths == [recipe.productPath]


def test_overwrite_runs_over_a_product_under_the_output_directory(
    log: Any, loggerPaths: list[str], tmp_path: Path
) -> None:
    # ARRANGE
    outputDir = tmp_path / "output"
    _write_fits(_product_under(outputDir, "reduced"))

    # ACT
    recipe = _construct(log, {"workspace-root-dir": str(outputDir)}, overwrite=True)

    # ASSERT
    assert recipe.productPath == str(_product_under(outputDir, "reduced"))
