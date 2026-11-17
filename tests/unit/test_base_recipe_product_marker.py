"""The product-name marker of `base_recipe` names products and leaves database keys alone (DY-2194)."""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories import instrument_header, pipeline_settings, qc_table

pytestmark = pytest.mark.unit

MARKER = "_STAREMODE"


def _bare_recipe(log: Any, tmpPath: Path, marker: str | None = None) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmpPath)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.sofName = "20240102T030405_NOD_OBJ"
    recipe.recipeName = "soxs-stare"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = qc_table().iloc[0:0].copy()
    if marker is not None:
        recipe.productNameMarker = marker
    return recipe


def test_the_product_stem_is_the_sof_name_when_no_marker_is_set(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path)

    # ACT / ASSERT
    assert recipe.productNameMarker == ""
    assert recipe.productStem == "20240102T030405_NOD_OBJ"


def test_the_product_stem_adds_the_marker_to_the_sof_name(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path, MARKER)

    # ACT / ASSERT
    assert recipe.productStem == f"20240102T030405_NOD_OBJ{MARKER}"
    assert recipe.sofName == "20240102T030405_NOD_OBJ"


def _frame() -> CCDData:
    return CCDData(
        np.full((4, 4), 7.0, dtype=np.float32),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros((4, 4), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((4, 4)), unit=u.electron),
    )


def test_write_names_a_default_product_with_the_marker(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path, MARKER)

    # ACT
    filepath = recipe._write(_frame(), str(tmp_path))

    # ASSERT
    assert Path(filepath).name == f"20240102T030405_NOD_OBJ{MARKER}.fits"
    assert Path(filepath).exists()


def test_write_names_a_default_product_without_a_marker_as_before(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path)

    # ACT
    filepath = recipe._write(_frame(), str(tmp_path))

    # ASSERT
    assert Path(filepath).name == "20240102T030405_NOD_OBJ.fits"


def test_the_cosmic_ray_mask_file_carries_the_marker_and_the_database_name_does_not(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path, MARKER)
    recipe.qcDir = str(tmp_path)
    recipe.cosmicRayMasks = [
        {"fileName": "raw_1.fits", "mjdObs": 60000.5, "dateObs": "2024-01-02T03:04:05", "mask": np.ones((2, 2), "u1")}
    ]

    # ACT
    recipe._write_cosmic_ray_mask_product()

    # ASSERT
    assert recipe.cosmicRayProduct["fileName"] == f"20240102T030405_NOD_OBJ{MARKER}_LACOSMIC_MASK.fits"
    assert (tmp_path / recipe.cosmicRayProduct["fileName"]).exists()
    assert recipe.cosmicRayProduct["sofName"] == "20240102T030405_NOD_OBJ"


def test_the_scattered_light_products_are_named_with_the_marker(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _bare_recipe(log, tmp_path, MARKER)
    recipe.products = pd.DataFrame()
    calls: list[dict[str, object]] = []
    frame = _frame()

    class RecordingBackground:
        def __init__(self, **kwargs: object) -> None:
            calls.append(kwargs)

        def subtract(self) -> tuple[CCDData, CCDData, pd.DataFrame]:
            return frame, frame, pd.DataFrame()

    monkeypatch.setattr(import_module("soxspipe.recipes.base_recipe"), "subtract_background", RecordingBackground)

    # ACT
    recipe._subtract_scattered_light(frame, "order_table")

    # ASSERT
    assert len(calls) == 1
    assert calls[0]["sofName"] == f"20240102T030405_NOD_OBJ{MARKER}"


def _isolate(monkeypatch: pytest.MonkeyPatch, tmpPath: Path, predictCalls: list[tuple]) -> None:
    class IsolatedOrganiser:
        def __init__(self, **_: object) -> None:
            return None

        def session_list(self, *, silent: bool) -> tuple[bool, list[str]]:
            return False, []

        def close(self) -> None:
            return None

    def predict(*args: object) -> tuple[str, str]:
        predictCalls.append(args)
        return str(tmpPath / "reduced" / "p.fits"), "2024-01-02"

    monkeypatch.setattr(commonutils, "data_organiser", IsolatedOrganiser)
    monkeypatch.setattr(toolkit, "predict_product_path", predict)
    monkeypatch.setattr(toolkit, "add_recipe_logger", lambda receivedLog, _: receivedLog)
    monkeypatch.setattr(toolkit, "get_calibrations_path", lambda **_: "calibrations")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmpPath / "qc"), str(tmpPath / "products")))


def _construct(log: Any, tmpPath: Path, **kwargs: object) -> base_recipe:
    return base_recipe(
        log=log,
        settings={"workspace-root-dir": str(tmpPath), "instrument": "soxs"},
        inputFrames=str(tmpPath / "20240102T030405_NOD_OBJ.sof"),
        recipeName="soxs-stare",
        turnOffMP=True,
        **kwargs,
    )


def test_the_constructor_predicts_the_product_path_with_the_marker(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    predictCalls: list[tuple] = []
    _isolate(monkeypatch, tmp_path, predictCalls)

    # ACT
    recipe = _construct(log, tmp_path, productNameMarker=MARKER)

    # ASSERT
    assert predictCalls == [(str(tmp_path / "20240102T030405_NOD_OBJ.sof"), "soxs-stare", MARKER)]
    assert recipe.productNameMarker == MARKER
    assert recipe.sofName == "20240102T030405_NOD_OBJ"


def test_the_constructor_predicts_an_unmarked_product_path_by_default(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    predictCalls: list[tuple] = []
    _isolate(monkeypatch, tmp_path, predictCalls)

    # ACT
    recipe = _construct(log, tmp_path)

    # ASSERT
    assert predictCalls == [(str(tmp_path / "20240102T030405_NOD_OBJ.sof"), "soxs-stare", "")]
    assert recipe.productStem == "20240102T030405_NOD_OBJ"
