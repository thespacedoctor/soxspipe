"""The recipes merge a detector's QC and product tables without the pandas concat `FutureWarning` (DY-2409).

`pd.concat` warns when an empty or all-NA table joins a table with data, because pandas 3 changes the result dtypes.
`toolkit.concat_tables` keeps the pandas 2 dtypes and raises no warning. Each test promotes that warning to an error.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
from soxspipe.recipes.soxs_order_centres import soxs_order_centres
from soxspipe.recipes.soxs_spatial_solution import soxs_spatial_solution
from tests.factories import prepared_fits, synthetic_ccd
from tests.unit.test_soxs_mflat_characterization import (
    _empty_products,
    _empty_qc,
    _new_bare_recipe,
    _stub_shared_collaborators,
)

pytestmark = [
    pytest.mark.unit,
    pytest.mark.filterwarnings("error:The behavior of DataFrame concatenation:FutureWarning"),
]

orderCentresModule = import_module("soxspipe.recipes.soxs_order_centres")


def _qc_row_with_a_missing_order() -> pd.DataFrame:
    """Return a one-row QC table whose ``qc_order`` column is all missing, as a recipe-level QC row is."""
    return pd.DataFrame(
        {
            "soxspipe_recipe": ["soxs-test"],
            "qc_name": ["RON"],
            "qc_value": [1.5],
            "qc_unit": ["e-"],
            "qc_order": [None],
            "qc_comment": ["read noise"],
            "obs_date_utc": ["2024-01-02T03:04:05"],
            "reduction_date_utc": ["2024-01-02T04:05:06"],
        }
    )


def _product_row() -> pd.DataFrame:
    """Return a one-row product table."""
    return pd.DataFrame(
        {
            "soxspipe_recipe": ["soxs-test"],
            "product_label": ["MAP"],
            "file_name": ["map.fits"],
            "file_type": ["FITS"],
            "obs_date_utc": ["2024-01-02T03:04:05"],
            "reduction_date_utc": ["2024-01-02T04:05:06"],
            "file_path": ["/products/map.fits"],
            "label": ["PROD"],
        }
    )


def _assert_the_rows_were_merged(recipe: Any, qcTable: pd.DataFrame, productsTable: pd.DataFrame) -> None:
    """Assert the detector's rows are in the recipe's tables, with the values and the columns unchanged."""
    assert recipe.qc["qc_name"].tolist() == qcTable["qc_name"].tolist()
    assert recipe.qc["qc_value"].tolist() == qcTable["qc_value"].tolist()
    assert recipe.qc["qc_value"].dtype == qcTable["qc_value"].dtype
    assert list(recipe.qc.columns) == list(_empty_qc().columns)
    assert recipe.products["product_label"].tolist()[:1] == productsTable["product_label"].tolist()
    assert list(recipe.products.columns)[: len(_empty_products().columns)] == list(_empty_products().columns)


def test_the_order_centres_recipe_merges_the_detector_tables_without_a_concat_warning(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    qcTable = _qc_row_with_a_missing_order()
    productsTable = _product_row()

    class Detector:
        def __init__(self, **_: object) -> None:
            return None

        def get(self) -> tuple[Any, ...]:
            return (None, qcTable, productsTable, pd.DataFrame(), pd.DataFrame(), pd.DataFrame())

    monkeypatch.setattr(orderCentresModule, "detect_continuum", Detector)
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.orderFrame = None
    recipe.settings = {}
    recipe.recipeSettings = {}
    recipe.sofName = "SYNTHETIC"
    recipe.startNightDate = "2024-01-02"
    recipe.qc = _empty_qc()
    recipe.products = _empty_products()

    # ACT
    with pytest.raises(ArithmeticError, match="Could not converge"):
        recipe._fit_order_centres(lambda keyword: keyword, "disp_map.fits", 1, 1)

    # ASSERT
    _assert_the_rows_were_merged(recipe, qcTable, productsTable)


def test_the_spatial_solution_recipe_merges_the_mapper_tables_without_a_concat_warning(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    qcTable = _qc_row_with_a_missing_order()
    productsTable = _product_row()
    mapPath = str(tmp_path / "map.fits")

    def mapper(**_: object) -> tuple[Any, ...]:
        return mapPath, None, None, qcTable, productsTable, pd.DataFrame()

    monkeypatch.setattr(commonutils, "create_dispersion_map", lambda **kwargs: SimpleNamespace(get=lambda: mapper()))
    monkeypatch.setattr(import_module("soxspipe.commonutils.toolkit"), "quicklook_image", lambda **_: None)
    recipe = soxs_spatial_solution.__new__(soxs_spatial_solution)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = "soxs-spatial-solution"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.settings = {}
    recipe.recipeSettings = {}
    recipe.multiPinholeFrame = None
    recipe.sofName = "SYNTHETIC"
    recipe.startNightDate = "2024-01-02"
    recipe.create2DMap = False
    recipe.slit_arc = None
    recipe.debug = False
    recipe.turnOffMP = True
    recipe.qc = _empty_qc()
    recipe.products = _empty_products()

    # ACT
    recipe._fit_spatial_solution("disp_map.fits", "order_table.fits")

    # ASSERT
    _assert_the_rows_were_merged(recipe, qcTable, productsTable)


def test_the_master_flat_recipe_merges_the_edge_detection_qc_without_a_concat_warning(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    orderPath = prepared_fits(tmp_path / "ORDER_TAB_VIS.fits", seed=331)
    flatFrame = synthetic_ccd(seed=332, prepared=True)
    qcTable = _qc_row_with_a_missing_order()
    recipe = _new_bare_recipe(
        log,
        tmp_path,
        orderPath,
        subtractBackground=False,
        calibratedFlatFiles=["flat-1.fits"],
        dFlatFiles=[],
        qFlatFiles=[],
        domeFlatFiles=[],
    )
    medianFlux = pd.DataFrame({"order": [10], "medianFlux": [100.0]})
    monkeypatch.setattr(recipe, "calibrate_frame_set", lambda: ([flatFrame], [], [], []))
    monkeypatch.setattr(recipe, "normalise_flats", lambda *a, **k: [flatFrame.copy()])
    monkeypatch.setattr(recipe, "clip_and_stack", lambda **k: k["frames"][0].copy())
    monkeypatch.setattr(recipe, "mask_low_sens_pixels", lambda **k: (k["frame"].copy(), medianFlux))
    monkeypatch.setattr(recipe, "_write", lambda *a, **k: str(tmp_path / "MASTER_FLAT_VIS.fits"))
    monkeypatch.setattr(recipe, "_detect_lamp_order_edges", lambda *a, **k: (qcTable, str(orderPath)))
    _stub_shared_collaborators(recipe, monkeypatch, orderPath)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert recipe.qc["qc_name"].tolist() == ["RON"]
    assert recipe.qc["qc_value"].dtype == qcTable["qc_value"].dtype
    assert list(recipe.qc.columns) == list(_empty_qc().columns)


FLUX_CALIBRATION_SITES = [
    ("soxs_nod", "soxs_nod", "_flux_calibrate_nod_stack", ("resp.fits", None, [], [])),
    ("soxs_offset", "soxs_offset", "_flux_calibrate_stack", ("resp.fits", None, [])),
    ("soxs_stare", "soxs_stare", "_flux_calibrate_spectrum", ("resp.fits", None, None, [])),
]


def _flux_calibrating_recipe(
    moduleName: str,
    className: str,
    arguments: tuple[Any, ...],
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[Any, tuple[Any, ...], list[dict[str, Any]], pd.DataFrame]:
    """Build a recipe whose flux calibrator is a double.

    Return the recipe, the arguments for its flux-calibration method, the keyword arguments the double was built
    with, and the product table the double returns.
    """
    module = import_module(f"soxspipe.recipes.{moduleName}")
    sequenceModule = import_module("soxspipe.commonutils.sequence_airmass")
    productsTable = _product_row()
    constructorArguments: list[dict[str, Any]] = []
    # A FLOAT COLUMN JOINING AN EMPTY OBJECT COLUMN IS WHAT MAKES pd.concat WARN
    productsTable["product_rank"] = [1.0]

    class FluxCalibrator:
        def __init__(self, **kwargs: object) -> None:
            constructorArguments.append(kwargs)

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            return "fluxcal.fits", productsTable

    monkeypatch.setattr(commonutils, "flux_calibration", FluxCalibrator)
    monkeypatch.setattr(import_module("soxspipe.commonutils.flux_calibration"), "flux_calibration", FluxCalibrator)
    monkeypatch.setattr(sequenceModule, "mean_airmass", lambda *_: 1.2)
    monkeypatch.setattr(sequenceModule, "require_equal_exposure_times", lambda *_: None)
    monkeypatch.setattr(module, "get_calibrations_path", lambda **_: "/cal")
    monkeypatch.setattr(
        module,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"extinction": "extinction.dat"}),
        raising=False,
    )
    recipe = getattr(module, className).__new__(getattr(module, className))
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = lambda keyword: keyword
    recipe.settings = {}
    recipe.detectorParams = {"extinction": "extinction.dat"}
    recipe.recipeName = "soxs-test"
    recipe.startNightDate = "2024-01-02"
    recipe.sofName = "SYNTHETIC"
    recipe.debug = False
    recipe.qc = _empty_qc()
    recipe.products = _empty_products().assign(product_rank=pd.Series([], dtype="object"))
    frame = SimpleNamespace(header={"EXPTIME": 30.0})
    callArguments = tuple([frame] if argument == [] else argument for argument in arguments)
    combined = SimpleNamespace(header={"EXPTIME": 30.0})
    callArguments = tuple(
        combined if argument is None and index == 2 else argument for index, argument in enumerate(callArguments)
    )
    return recipe, callArguments, constructorArguments, productsTable


@pytest.mark.parametrize(("moduleName", "className", "method", "arguments"), FLUX_CALIBRATION_SITES)
def test_the_flux_calibrating_recipes_merge_the_flux_calibration_products_without_a_concat_warning(
    moduleName: str,
    className: str,
    method: str,
    arguments: tuple[Any, ...],
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    recipe, callArguments, _, _ = _flux_calibrating_recipe(moduleName, className, arguments, log, monkeypatch)

    # ACT
    getattr(recipe, method)(*callArguments)

    # ASSERT
    assert recipe.products["product_label"].tolist() == ["MAP"]
    assert recipe.products.index.tolist() == [0]
    assert recipe.products["product_rank"].tolist() == [1.0]
    assert recipe.products["product_rank"].dtype == "float64"


@pytest.mark.parametrize(("moduleName", "className", "method", "arguments"), FLUX_CALIBRATION_SITES)
def test_the_flux_calibrating_recipes_hand_their_qc_table_to_the_calibrator_for_the_spectrum_header(
    moduleName: str,
    className: str,
    method: str,
    arguments: tuple[Any, ...],
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    recipe, callArguments, constructorArguments, _ = _flux_calibrating_recipe(
        moduleName, className, arguments, log, monkeypatch
    )

    # ACT
    getattr(recipe, method)(*callArguments)

    # ASSERT
    assert constructorArguments[0]["qc"] is recipe.qc
