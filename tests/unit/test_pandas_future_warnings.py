"""Sites that raised a pandas ``FutureWarning`` and would break in pandas 3 (DY-2280).

Every test runs with the two pandas ``FutureWarning`` messages promoted to
errors: the incompatible-dtype assignment that pandas 3 turns into an error,
and the ``pd.concat`` onto an empty or all-missing table that changes the
result dtypes.
"""

from __future__ import annotations

from importlib import import_module
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from scipy.interpolate import splrep

from soxspipe.commonutils import toolkit
from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_disp_solution import soxs_disp_solution
from tests.unit.test_create_dispersion_map_helpers import _mapper
from tests.unit.test_sky_subtraction import _subtractor

DISPERSION_MODULE = import_module("soxspipe.recipes.soxs_disp_solution")

# THE TWO PANDAS WARNINGS THESE SITES RAISED. pandas REPORTS THEM AT THE CALLING LINE,
# SO A MODULE FILTER CANNOT SELECT THEM
pytestmark = [
    pytest.mark.unit,
    pytest.mark.filterwarnings("error:Setting an item of incompatible dtype:FutureWarning"),
    pytest.mark.filterwarnings("error:The behavior of DataFrame concatenation:FutureWarning"),
]

NOT_FLAGGED = "none"


def _empty_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return the empty QC and product tables every recipe starts with."""
    return base_recipe._empty_qc_and_product_tables(None, pd)


def _qc_row_kwargs(**overrides: Any) -> dict[str, Any]:
    return {
        "recipeName": "soxs-disp-solution",
        "qcName": "RON",
        "qcValue": 1.5,
        "qcComment": "A comment",
        "obsDateUtc": "2024-01-02T03:04:05",
        "reductionDateUtc": "2024-01-02T04:05:06",
        **overrides,
    }


# ---------------------------------------------------------------------------
# subtract_sky: flagged_sky_line
# ---------------------------------------------------------------------------


def test_sky_line_flag_column_holds_strings_with_an_explicit_not_flagged_sentinel(log: object) -> None:
    # ARRANGE: THE FIXTURE OF test_residual_floor_normalises_errors_and_flags_local_excess
    subtractor = _subtractor(log)
    wavelength = np.linspace(500.0, 501.0, 32)
    pixels = pd.DataFrame(
        {
            "wavelength": wavelength,
            "flux": np.full(32, 10.0),
            "error": np.ones(32),
            "flagged_all_clipped": False,
            "residual_windowed_long_median": np.zeros(32),
            "flux_windowed_long_median": np.full(32, 10.0),
            "sky_residual_rolling_average": np.full(32, 2.0),
        }
    )
    pixels.loc[10, "sky_residual_rolling_average"] = 3.0
    spline = splrep(wavelength, np.full(32, 8.0), k=1)

    # ACT
    result, _ = subtractor.determine_residual_floor(pixels, spline, iteration=1)

    # ASSERT
    flags = result["flagged_sky_line"]
    assert flags.map(type).eq(str).all()
    assert flags[10] == "line"
    assert (flags.drop(index=10) == NOT_FLAGGED).all()


# ---------------------------------------------------------------------------
# toolkit: append_qc
# ---------------------------------------------------------------------------


def test_append_qc_onto_the_empty_qc_table_keeps_the_declared_columns_and_dtypes() -> None:
    qcTable, _ = _empty_tables()

    result = toolkit.append_qc(qcTable, **_qc_row_kwargs(qcOrder=3, toHeader=True))

    assert list(result.columns) == list(qcTable.columns)
    assert result["qc_value"].dtype == np.float64
    assert result["qc_order"].dtype == object
    assert result["to_header"].dtype == object
    assert result["qc_unit"].dtype == object
    assert result.loc[0, "qc_order"] == 3


def test_append_qc_onto_a_float_column_with_a_missing_order_does_not_warn() -> None:
    qcTable, _ = _empty_tables()
    first = toolkit.append_qc(qcTable, **_qc_row_kwargs(qcOrder=3))

    result = toolkit.append_qc(first, **_qc_row_kwargs(qcName="FLUX", qcOrder=None, qcUnit=None))

    assert len(result) == 2
    assert result["qc_value"].dtype == np.float64


# ---------------------------------------------------------------------------
# create_dispersion_map.calculate_residuals
# ---------------------------------------------------------------------------


def test_dispersion_residual_qcs_append_onto_the_empty_qc_table_without_a_future_warning(log: object) -> None:
    mapper = _mapper(log)
    mapper.arcFrame = False
    mapper.qc, _ = _empty_tables()
    table = pd.DataFrame(
        {
            "order": [1.0, 1.0],
            "wavelength": [10.0, 20.0],
            "slit_position": [0.0, 0.0],
            "order_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_1": [10.0, 20.0],
            "slit_position_pow_x_0": [1.0, 1.0],
            "order_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_1": [10.0, 20.0],
            "slit_position_pow_y_0": [1.0, 1.0],
            "observed_x": [17.0, 40.0],
            "observed_y": [1.0, 5.0],
            "detector_x": [20.0, 40.0],
            "detector_y": [5.0, 5.0],
            "fwhm_pin_px": [2.0, 2.0],
        }
    )

    mapper.calculate_residuals(
        table,
        xcoeff=[0.0, 2.0],
        ycoeff=[5.0, 0.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
        pixelRange=True,
        writeQCs=True,
    )

    assert len(mapper.qc) > 0
    assert mapper.qc["qc_value"].dtype == np.float64
    assert mapper.qc["qc_order"].dtype == object


# ---------------------------------------------------------------------------
# soxs_disp_solution._fit_dispersion_solution
# ---------------------------------------------------------------------------


def test_dispersion_recipe_merges_the_map_tables_into_empty_recipe_tables_without_a_future_warning(
    log: object, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE DISPERSION MAP RETURNS ONE QC ROW AND ONE PRODUCT ROW
    qcTable, productTable = _empty_tables()
    mapQc = toolkit.append_qc(qcTable, **_qc_row_kwargs(qcOrder=3, toHeader=True))
    mapProducts = toolkit.append_product(
        productTable,
        recipeName="soxs-disp-solution",
        productLabel="DISP_MAP",
        fileName="map.fits",
        filePath="map.fits",
        productDesc="A map",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    class FakeMapper:
        def __init__(self, **_: object) -> None:
            return None

        def get(self) -> tuple:
            return "map.fits", "image.pdf", [], mapQc, mapProducts, pd.DataFrame()

    monkeypatch.setattr(DISPERSION_MODULE, "create_dispersion_map", FakeMapper)
    recipe = soxs_disp_solution.__new__(soxs_disp_solution)
    recipe.log = log
    recipe.polyOrders = False
    recipe.recipeSettings = {}
    recipe.settings = {}
    recipe.pinholeFrame = SimpleNamespace(header={"DATE-OBS": "2024-01-02T03:04:05"})
    recipe.sofName = "sof"
    recipe.startNightDate = "2024-01-02"
    recipe.turnOffMP = True
    recipe.recipeName = "soxs-disp-solution"
    recipe.arm = "VIS"
    recipe.qc, recipe.products = _empty_tables()
    recipe.report_output = lambda: recipe.qc
    recipe.clean_up = lambda: None

    # ACT
    productPath, reported = recipe._fit_dispersion_solution(lambda keyword: "DATE-OBS")

    # ASSERT
    assert productPath == "map.fits"
    assert len(reported) == 1
    assert reported["qc_value"].dtype == np.float64
    assert reported["qc_order"].dtype == object
    assert list(recipe.products["product_label"]) == ["DISP_MAP", "DISP_MAP"]


# ---------------------------------------------------------------------------
# toolkit.concat_tables
# ---------------------------------------------------------------------------


def test_concat_tables_returns_an_empty_table_when_every_table_is_empty() -> None:
    qcTable, _ = _empty_tables()

    result = toolkit.concat_tables([qcTable, qcTable])

    assert result.empty
    assert list(result.columns) == list(qcTable.columns)


def test_concat_tables_keeps_the_index_of_the_non_empty_tables_when_not_ignored() -> None:
    rows = pd.DataFrame({"a": [1.5, 2.5]}, index=[4, 7])

    result = toolkit.concat_tables([pd.DataFrame(), rows, rows])

    assert list(result.index) == [4, 7, 4, 7]


def test_concat_tables_gives_an_all_missing_column_the_dtype_of_the_other_table() -> None:
    held = pd.DataFrame({"value": [1.5], "text": ["a"]})
    missing = pd.DataFrame({"value": [None], "text": [None]})

    result = toolkit.concat_tables([held, missing], ignore_index=True)

    assert result["value"].dtype == np.float64
    assert result["text"].dtype == object
    assert result["value"].isna().tolist() == [False, True]


def test_concat_tables_makes_an_integer_column_object_when_the_other_table_has_it_all_missing() -> None:
    held = pd.DataFrame({"order": [3]})
    missing = pd.DataFrame({"order": [None]})

    result = toolkit.concat_tables([held, missing], ignore_index=True)

    assert result["order"].dtype == object
    assert result.loc[0, "order"] == 3


def test_concat_tables_adds_the_columns_of_an_empty_table_after_its_own_order() -> None:
    empty = pd.DataFrame({"b": pd.Series([], dtype=object), "a": pd.Series([], dtype=object)})
    rows = pd.DataFrame({"a": [1.5], "c": ["x"]})

    result = toolkit.concat_tables([empty, rows], ignore_index=True)

    assert list(result.columns) == ["b", "a", "c"]
    assert result["b"].dtype == object
    assert result["a"].dtype == np.float64


@pytest.mark.parametrize("floatFirst", [True, False])
def test_concat_tables_gives_an_all_missing_column_float_when_any_other_table_holds_floats(floatFirst: bool) -> None:
    floats = pd.DataFrame({"value": [1.5]})
    integers = pd.DataFrame({"value": [3]})
    donors = [floats, integers] if floatFirst else [integers, floats]
    missing = pd.DataFrame({"value": [None]})

    result = toolkit.concat_tables([*donors, missing], ignore_index=True)

    assert result["value"].dtype == np.float64
    assert result["value"].isna().tolist() == [False, False, True]


@pytest.mark.parametrize("objectFirst", [True, False])
def test_concat_tables_gives_an_all_missing_column_object_when_a_donor_is_object_and_none_is_float(
    objectFirst: bool,
) -> None:
    objects = pd.DataFrame({"value": ["a"]})
    integers = pd.DataFrame({"value": [3]})
    donors = [objects, integers] if objectFirst else [integers, objects]
    missing = pd.DataFrame({"value": [None]})

    result = toolkit.concat_tables([*donors, missing], ignore_index=True)

    assert result["value"].dtype == object


def test_concat_tables_rejects_duplicate_column_labels_naming_them() -> None:
    duplicated = pd.DataFrame([[1.5, 2.5, 3.5]], columns=["a", "b", "a"])

    with pytest.raises(ValueError, match="duplicate column labels.*'a'"):
        toolkit.concat_tables([pd.DataFrame({"a": [1.0]}), duplicated])


def test_concat_tables_rebuilds_an_empty_table_column_on_the_index_of_the_rows() -> None:
    empty = pd.DataFrame({"a": pd.Series([], dtype=object), "b": pd.Series([], dtype=object)})
    rows = pd.DataFrame({"a": [1.5, 2.5]}, index=[4, 7])

    result = toolkit.concat_tables([empty, rows])

    assert list(result.index) == [4, 7]
    assert list(result.columns) == ["a", "b"]
    assert result["b"].isna().all()
    assert result["b"].dtype == object
