"""Behavioural contracts for the shared QC/product/plot helper functions
and the ``base_recipe`` delegators that forward to them."""

from __future__ import annotations

import re
from datetime import datetime
from typing import Any

import pandas as pd
import pytest

from soxspipe.commonutils import toolkit
from soxspipe.recipes import base_recipe
from tests.factories import product_table, qc_table

pytestmark = pytest.mark.unit


def _recipe(log: Any) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-mbias"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = pd.DataFrame()
    recipe.products = pd.DataFrame()
    return recipe


# ---------------------------------------------------------------------------
# utcnow_string
# ---------------------------------------------------------------------------


def test_utcnow_string_returns_second_precision_timestamp() -> None:
    result = toolkit.utcnow_string()

    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}$", result)


def test_utcnow_string_adds_fractional_seconds_when_microseconds_requested() -> None:
    result = toolkit.utcnow_string(microseconds=True)

    assert re.match(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{6}$", result)


def test_utcnow_string_matches_naive_utcnow_strftime_shape() -> None:
    reference = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

    result = toolkit.utcnow_string()

    assert len(result) == len(reference)
    assert result[10] == "T" == reference[10]
    assert "+" not in result
    assert "Z" not in result


# ---------------------------------------------------------------------------
# append_qc
# ---------------------------------------------------------------------------


def test_omitted_sentinel_has_a_readable_repr() -> None:
    assert repr(toolkit.OMITTED) == "<omitted>"


def test_append_qc_returns_new_frame_and_leaves_input_unchanged() -> None:
    original = qc_table().iloc[0:0].copy()

    result = toolkit.append_qc(
        original,
        recipeName="soxs-mbias",
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert len(original) == 0
    assert len(result) == 1


def test_append_qc_omits_optional_columns_when_not_passed() -> None:
    result = toolkit.append_qc(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert "qc_unit" not in result.columns
    assert "to_header" not in result.columns
    assert "qc_order" not in result.columns


def test_append_qc_includes_optional_columns_when_passed_explicitly_none_or_false() -> None:
    result = toolkit.append_qc(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
        qcUnit=None,
        qcOrder=None,
        toHeader=False,
    )

    assert "qc_unit" in result.columns
    assert result.loc[0, "qc_unit"] is None or pd.isna(result.loc[0, "qc_unit"])
    assert "qc_order" in result.columns
    assert bool(result.loc[0, "to_header"]) is False


def test_append_qc_row_column_order_matches_predeclared_table_for_fresh_frame() -> None:
    result = toolkit.append_qc(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
        qcUnit="electron",
        qcOrder=3,
        toHeader=True,
    )

    assert list(result.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_unit",
        "qc_order",
        "qc_comment",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]


def test_append_qc_preserves_numeric_dtype_of_qc_value_column() -> None:
    seedTable = qc_table()

    result = toolkit.append_qc(
        seedTable,
        recipeName="soxs-mbias",
        qcName="STRUCTX",
        qcValue=0.0042,
        qcComment="Slope of bias in X",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
        qcUnit=None,
        qcOrder=None,
        toHeader=True,
    )

    assert pd.api.types.is_float_dtype(result["qc_value"])
    assert len(result) == 2


# ---------------------------------------------------------------------------
# append_product
# ---------------------------------------------------------------------------


def test_append_product_returns_new_frame_and_leaves_input_unchanged() -> None:
    original = product_table().iloc[0:0].copy()

    result = toolkit.append_product(
        original,
        recipeName="soxs-mbias",
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert len(original) == 0
    assert len(result) == 1


def test_append_product_omits_optional_columns_when_not_passed() -> None:
    result = toolkit.append_product(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert "file_type" not in result.columns
    assert "label" not in result.columns


def test_append_product_includes_optional_columns_when_passed_explicitly_none() -> None:
    result = toolkit.append_product(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
        fileType=None,
        label=None,
    )

    assert "file_type" in result.columns
    assert result.loc[0, "file_type"] is None or pd.isna(result.loc[0, "file_type"])
    assert "label" in result.columns
    assert result.loc[0, "label"] is None or pd.isna(result.loc[0, "label"])


def test_append_product_row_column_order_matches_call_site_shape() -> None:
    result = toolkit.append_product(
        pd.DataFrame(),
        recipeName="soxs-mbias",
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        obsDateUtc="2024-01-02T03:04:05",
        reductionDateUtc="2024-01-02T04:05:06",
        fileType="FITS",
        label="PROD",
    )

    assert list(result.columns) == [
        "soxspipe_recipe",
        "product_label",
        "file_name",
        "file_type",
        "obs_date_utc",
        "reduction_date_utc",
        "product_desc",
        "file_path",
        "label",
    ]


# ---------------------------------------------------------------------------
# save_qc_plot
# ---------------------------------------------------------------------------


def test_save_qc_plot_passes_bbox_inches_when_tight(monkeypatch: pytest.MonkeyPatch) -> None:
    import matplotlib.pyplot as plt

    calls: list[tuple[tuple, dict]] = []
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: calls.append((args, kwargs)))

    result = toolkit.save_qc_plot("plot.pdf", bboxInches="tight")

    assert result == "plot.pdf"
    args, kwargs = calls[0]
    assert args == ("plot.pdf",)
    assert kwargs == {"dpi": 120, "format": "pdf", "bbox_inches": "tight"}


def test_save_qc_plot_omits_bbox_inches_when_none(monkeypatch: pytest.MonkeyPatch) -> None:
    import matplotlib.pyplot as plt

    calls: list[tuple[tuple, dict]] = []
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: calls.append((args, kwargs)))

    toolkit.save_qc_plot("plot.pdf", bboxInches=None)

    _, kwargs = calls[0]
    assert "bbox_inches" not in kwargs
    assert kwargs == {"dpi": 120, "format": "pdf"}


def test_save_qc_plot_honours_non_default_dpi(monkeypatch: pytest.MonkeyPatch) -> None:
    import matplotlib.pyplot as plt

    calls: list[tuple[tuple, dict]] = []
    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: calls.append((args, kwargs)))

    toolkit.save_qc_plot("plot.pdf", dpi=360, bboxInches=None)

    _, kwargs = calls[0]
    assert kwargs["dpi"] == 360


def test_save_qc_plot_returns_the_file_path(monkeypatch: pytest.MonkeyPatch) -> None:
    import matplotlib.pyplot as plt

    monkeypatch.setattr(plt, "savefig", lambda *args, **kwargs: None)

    result = toolkit.save_qc_plot("plot.pdf")

    assert result == "plot.pdf"


# ---------------------------------------------------------------------------
# base_recipe.add_qc delegator
# ---------------------------------------------------------------------------


def test_add_qc_appends_row_using_recipe_name_and_dateobs_by_default(log: Any) -> None:
    recipe = _recipe(log)

    result = recipe.add_qc(
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert result is None
    assert len(recipe.qc) == 1
    assert recipe.qc.loc[0, "soxspipe_recipe"] == "soxs-mbias"
    assert recipe.qc.loc[0, "obs_date_utc"] == "2024-01-02T03:04:05"
    assert "qc_unit" not in recipe.qc.columns


def test_add_qc_honours_explicit_recipe_name_and_obs_date_overrides(log: Any) -> None:
    recipe = _recipe(log)

    recipe.add_qc(
        qcName="RON",
        qcValue=1.2,
        qcComment="Read noise",
        reductionDateUtc="2024-01-02T04:05:06",
        recipeName="soxs-stare",
        obsDateUtc="2099-01-01T00:00:00",
    )

    assert recipe.qc.loc[0, "soxspipe_recipe"] == "soxs-stare"
    assert recipe.qc.loc[0, "obs_date_utc"] == "2099-01-01T00:00:00"


# ---------------------------------------------------------------------------
# base_recipe.add_product delegator
# ---------------------------------------------------------------------------


def test_add_product_appends_row_using_recipe_name_and_dateobs_by_default(log: Any) -> None:
    recipe = _recipe(log)

    result = recipe.add_product(
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        reductionDateUtc="2024-01-02T04:05:06",
    )

    assert result is None
    assert len(recipe.products) == 1
    assert recipe.products.loc[0, "soxspipe_recipe"] == "soxs-mbias"
    assert recipe.products.loc[0, "obs_date_utc"] == "2024-01-02T03:04:05"
    assert "file_type" not in recipe.products.columns


def test_add_product_honours_explicit_recipe_name_and_obs_date_overrides(log: Any) -> None:
    recipe = _recipe(log)

    recipe.add_product(
        productLabel="MBIAS",
        fileName="MASTER_BIAS.fits",
        filePath="products/MASTER_BIAS.fits",
        productDesc="Master bias frame",
        reductionDateUtc="2024-01-02T04:05:06",
        recipeName="soxs-stare",
        obsDateUtc="2099-01-01T00:00:00",
    )

    assert recipe.products.loc[0, "soxspipe_recipe"] == "soxs-stare"
    assert recipe.products.loc[0, "obs_date_utc"] == "2099-01-01T00:00:00"
