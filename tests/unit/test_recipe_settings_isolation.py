"""A recipe's changes to its own settings must not reach the settings dictionary it was given (DY-1724).

The serial reducer passes one settings dictionary to every recipe it runs in a process.
`base_recipe.get_recipe_settings` returns each recipe's own copy of its block, with the arm values merged in.
"""

from __future__ import annotations

import copy
import importlib
from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_stare import soxs_stare

pytestmark = pytest.mark.unit

stareModule = importlib.import_module("soxspipe.recipes.soxs_stare")


def _shared_settings() -> dict[str, Any]:
    """Return the settings a serial `reduce all` run hands to every stare recipe."""
    return {
        "instrument": "soxs",
        "soxs-stare": {
            "horne-extraction-slit-length": 15,
            "sky-subtraction": {"subtract_sky": True},
            "vis": {
                "horne-extraction-slit-length": 15,
                "sky-subtraction": {"subtract_sky": True},
            },
            "nir": {"sky-subtraction": {"subtract_sky": True}},
        },
        "soxs-disp-solution": {
            "order-deg": [4, 4],
            "wavelength-deg": [4, 4],
            "vis": {"pixel-window-size": 10},
        },
    }


def _recipe(cls: type, log: Any, settings: dict[str, Any], recipeName: str, arm: str) -> Any:
    """Return a recipe holding its settings the way its constructor leaves them."""
    recipe = cls.__new__(cls)
    recipe.log = log
    recipe.settings = settings
    recipe.recipeName = recipeName
    recipe.arm = arm
    recipe.recipeSettings = recipe.get_recipe_settings()
    return recipe


def _frame() -> CCDData:
    return CCDData(np.ones((3, 3)), unit=u.electron, uncertainty=StdDevUncertainty(np.ones((3, 3))))


def test_the_stare_sky_model_fallback_leaves_the_shared_settings_unchanged(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    settings = _shared_settings()
    snapshot = copy.deepcopy(settings)

    class SubtractSkyWithNoModel:
        def __init__(self, **_: Any) -> None:
            return None

        def subtract(self) -> tuple[None, None, None, str, str]:
            return None, None, None, "qc", "products"

    monkeypatch.setattr(stareModule, "subtract_sky", SubtractSkyWithNoModel)
    recipe = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")
    recipe.subtractSky = True
    recipe.qc, recipe.products = "qc", "products"
    recipe.sofName, recipe.startNightDate, recipe.debug = "stare", "2024-01-02", False
    recipe._flag_excluded_pixels = lambda frame, mask: frame

    # ACT
    recipe._subtract_sky(_frame(), "map.fits", "disp.fits", "orders.fits")

    # ASSERT: THE FALLBACK RAN, AND ONLY THE RECIPE'S OWN COPY SAW IT
    assert recipe.recipeSettings["sky-subtraction"]["subtract_sky"] is False
    assert recipe.recipeSettings["horne-extraction-slit-length"] == 30.0
    assert settings == snapshot
    assert settings["soxs-stare"]["vis"]["sky-subtraction"]["subtract_sky"] is True


def test_a_later_recipe_reads_the_settings_the_first_recipe_read(log: Any) -> None:
    # ARRANGE
    settings = _shared_settings()
    first = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")

    # ACT
    first.recipeSettings["sky-subtraction"]["subtract_sky"] = False
    first.recipeSettings["horne-extraction-slit-length"] = 30.0
    second = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")

    # ASSERT
    assert second.recipeSettings["sky-subtraction"]["subtract_sky"] is True
    assert second.recipeSettings["horne-extraction-slit-length"] == 15


def test_two_recipes_of_one_arm_share_no_nested_settings_objects(log: Any) -> None:
    # ARRANGE
    settings = _shared_settings()

    # ACT
    first = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")
    second = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")

    # ASSERT
    assert first.recipeSettings == second.recipeSettings
    assert first.recipeSettings["sky-subtraction"] is not second.recipeSettings["sky-subtraction"]
    assert first.recipeSettings["sky-subtraction"] is not settings["soxs-stare"]["sky-subtraction"]
    assert first.recipeSettings["sky-subtraction"] is not settings["soxs-stare"]["vis"]["sky-subtraction"]


def test_a_command_line_polynomial_override_leaves_the_shared_settings_unchanged(log: Any) -> None:
    # ARRANGE
    settings = _shared_settings()
    snapshot = copy.deepcopy(settings)
    recipe = _recipe(base_recipe, log, settings, "soxs-disp-solution", "VIS")

    # ACT: WHAT soxs_disp_solution DOES WITH --poly
    recipe.recipeSettings["order-deg"] = [3, 3]
    recipe.recipeSettings["wavelength-deg"] = [5, 5]

    # ASSERT
    assert settings == snapshot


def test_the_arm_values_are_merged_into_the_recipe_copy_only(log: Any) -> None:
    # ARRANGE
    settings = _shared_settings()
    settings["soxs-stare"]["vis"]["horne-extraction-slit-length"] = 12

    # ACT
    recipe = _recipe(soxs_stare, log, settings, "soxs-stare", "VIS")

    # ASSERT
    assert recipe.recipeSettings["horne-extraction-slit-length"] == 12
    assert settings["soxs-stare"]["horne-extraction-slit-length"] == 15


def test_a_recipe_with_no_settings_block_gets_false(log: Any) -> None:
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.log, recipe.settings, recipe.recipeName, recipe.arm = log, {}, False, "VIS"

    # ACT / ASSERT
    assert recipe.get_recipe_settings() is False


def test_the_header_records_the_parameters_the_recipe_ran_with(log: Any) -> None:
    # ARRANGE
    settings = _shared_settings()
    recipe = _recipe(base_recipe, log, settings, "soxs-disp-solution", "VIS")
    recipe.recipeSettings["order-deg"] = [3, 3]
    frame = CCDData(np.ones((2, 2)), unit=u.electron, header=fits.Header())

    # ACT
    recipe._stamp_recipe_parameter_records(frame)

    # ASSERT
    stamped = {
        frame.header[f"ESO PRO REC1 PARAM{n} NAME"]: frame.header[f"ESO PRO REC1 PARAM{n} VALUE"]
        for n in range(1, 10)
        if f"ESO PRO REC1 PARAM{n} NAME" in frame.header
    }
    assert stamped["order-deg"] == "3, 3"
