"""Failure-propagation contracts for the recipe parameter-tuning workers."""

from __future__ import annotations

import importlib
import sys
from typing import Any

import pytest

pytestmark = pytest.mark.unit


class _FailingMapper:
    """Stands in for a pipeline object whose ``get`` call fails."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    def get(self) -> Any:
        raise RuntimeError("tuning iteration failed")


def _tuning_kwargs(**overrides: Any) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "recipeSettings": {},
        "settings": {},
        "qc": None,
        "products": None,
        "sofName": "sof",
        "lineDetectionTable": None,
    }
    kwargs.update(overrides)
    return kwargs


def test_disp_solution_parameter_tuning_propagates_iteration_failures(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # THE PACKAGE __init__ BINDS THE NAME TO THE RECIPE CLASS, NOT THE MODULE
    recipeModule = importlib.import_module("soxspipe.recipes.soxs_disp_solution")

    # THE WORKER RE-IMPORTS THE CLASS FROM THE PACKAGE, SO PATCH IT THERE
    commonutils = sys.modules["soxspipe.commonutils"]
    monkeypatch.setattr(commonutils, "create_dispersion_map", _FailingMapper)

    with pytest.raises(RuntimeError, match="tuning iteration failed"):
        recipeModule.parameterTuning(
            (1, 1, 1, 1),
            log=log,
            pinholeFrame=None,
            **_tuning_kwargs(),
        )


def test_order_centres_parameter_tuning_propagates_iteration_failures(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # THE PACKAGE __init__ BINDS THE NAME TO THE RECIPE CLASS, NOT THE MODULE
    recipeModule = importlib.import_module("soxspipe.recipes.soxs_order_centres")

    # THE MODULE IMPORTS THE DETECTOR BY NAME AT IMPORT TIME, SO PATCH IT ON THE MODULE
    monkeypatch.setattr(recipeModule, "detect_continuum", _FailingMapper)

    with pytest.raises(RuntimeError, match="tuning iteration failed"):
        recipeModule.parameterTuning(
            (3, 5),
            log=log,
            recipeSettings={"detect-continuum": {}},
            settings={},
            orderFrame=None,
            disp_map_table=None,
            orderPixelTable=None,
            qc=None,
            products=None,
            sofName="sof",
            binx=1,
            biny=1,
            startNightDate="2024-01-02",
        )


def test_spatial_solution_parameter_tuning_propagates_iteration_failures(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # THE PACKAGE __init__ BINDS THE NAME TO THE RECIPE CLASS, NOT THE MODULE
    recipeModule = importlib.import_module("soxspipe.recipes.soxs_spatial_solution")

    # THE WORKER RE-IMPORTS THE CLASS FROM THE PACKAGE, SO PATCH IT THERE
    commonutils = sys.modules["soxspipe.commonutils"]
    monkeypatch.setattr(commonutils, "create_dispersion_map", _FailingMapper)

    with pytest.raises(RuntimeError, match="tuning iteration failed"):
        recipeModule.parameterTuning(
            (3, 4, 5, 4, 3, 5),
            log=log,
            multiPinholeFrame=None,
            disp_map_table=None,
            order_table=None,
            debug=False,
            **_tuning_kwargs(),
        )


# ---------------------------------------------------------------------------
# DY-1535: THE TUNING POOL HONOURS turnOffMP INSIDE BULK WORKERS
# ---------------------------------------------------------------------------


class _StubMapper:
    """Stands in for ``create_dispersion_map`` and returns an empty line-detection result."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    def get(self) -> tuple[None, None, None, None, None, None]:
        return None, None, None, None, None, None


class _StubDetector:
    """Stands in for ``detect_continuum`` and returns an empty trace sample."""

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        pass

    def sample_trace(self) -> tuple[None, float]:
        return None, 100.0


def _bare_recipe(recipeModule: Any, className: str, log: Any, turnOffMP: bool, debug: bool) -> Any:
    recipe = object.__new__(getattr(recipeModule, className))
    recipe.log = log
    recipe.turnOffMP = turnOffMP
    recipe.debug = debug
    recipe.settings = {}
    recipe.recipeSettings = {}
    recipe.qc = None
    recipe.products = None
    recipe.sofName = "sof"
    recipe.startNightDate = "2024-01-02"
    recipe.pinholeFrame = None
    recipe.multiPinholeFrame = None
    recipe.orderFrame = None
    return recipe


def _run_tuning(recipeName: str, recipe: Any) -> None:
    if recipeName == "soxs_disp_solution":
        recipe._tune_dispersion_parameters()
    elif recipeName == "soxs_order_centres":
        recipe._tune_order_centre_parameters(iter([(3, 5)]), None, 1, 1)
    else:
        recipe._tune_spatial_parameters(iter([(3, 4, 5, 4, 3, 5)]), None, None)


@pytest.mark.parametrize("recipeName", ["soxs_disp_solution", "soxs_order_centres", "soxs_spatial_solution"])
@pytest.mark.parametrize(
    ("turnOffMP", "debug", "expectedTurnOffMP"),
    [(True, False, True), (False, True, True), (False, False, False)],
)
def test_parameter_tuning_pool_runs_serially_when_multiprocessing_is_turned_off(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Any,
    recipeName: str,
    turnOffMP: bool,
    debug: bool,
    expectedTurnOffMP: bool,
) -> None:
    # ARRANGE: STUB THE LINE DETECTION AND CAPTURE THE POOL CALL INSTEAD OF STARTING IT
    import fundamentals

    monkeypatch.chdir(tmp_path)
    recipeModule = importlib.import_module(f"soxspipe.recipes.{recipeName}")
    commonutils = sys.modules["soxspipe.commonutils"]
    monkeypatch.setattr(commonutils, "create_dispersion_map", _StubMapper)
    if hasattr(recipeModule, "create_dispersion_map"):
        monkeypatch.setattr(recipeModule, "create_dispersion_map", _StubMapper)
    if hasattr(recipeModule, "detect_continuum"):
        monkeypatch.setattr(recipeModule, "detect_continuum", _StubDetector)
    poolCalls: list[dict[str, Any]] = []
    # THE RECIPES IMPORT fmultiprocess INSIDE THE TUNING METHOD, SO PATCHING THE PACKAGE ATTRIBUTE REACHES THEM
    monkeypatch.setattr(fundamentals, "fmultiprocess", lambda **kwargs: poolCalls.append(kwargs))
    recipe = _bare_recipe(recipeModule, recipeName, log, turnOffMP=turnOffMP, debug=debug)

    # ACT
    _run_tuning(recipeName, recipe)

    # ASSERT: A BULK WORKER (turnOffMP=True) MUST NOT START A DAEMONIC CHILD POOL
    assert len(poolCalls) == 1
    assert poolCalls[0]["turnOffMP"] is expectedTurnOffMP
