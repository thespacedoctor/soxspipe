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
