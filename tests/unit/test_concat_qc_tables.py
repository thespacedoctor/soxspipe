"""Contracts for ``base_recipe._concat_qc_tables``, which joins the QC tables of several reduction passes (DY-1723)."""

import pandas as pd
import pytest

from soxspipe.recipes.base_recipe import base_recipe

pytestmark = pytest.mark.unit


def test_no_tables_give_an_empty_table():
    # ACT
    qcTable = base_recipe._concat_qc_tables([])

    # ASSERT
    assert isinstance(qcTable, pd.DataFrame)
    assert qcTable.empty


def test_empty_tables_are_skipped_so_they_do_not_decide_the_dtypes():
    # ARRANGE
    filled = pd.DataFrame({"qc_name": ["A", "B"], "qc_value": [1.5, 2.5]})
    empty = pd.DataFrame({"qc_name": pd.Series(dtype=object), "qc_value": pd.Series(dtype=object)})

    # ACT
    qcTable = base_recipe._concat_qc_tables([empty, filled, empty])

    # ASSERT
    assert qcTable["qc_name"].tolist() == ["A", "B"]
    assert qcTable["qc_value"].dtype == float
    assert qcTable.index.tolist() == [0, 1]


def test_only_empty_tables_give_an_empty_table_with_their_columns():
    # ARRANGE
    empty = pd.DataFrame({"qc_name": pd.Series(dtype=object)})

    # ACT
    qcTable = base_recipe._concat_qc_tables([empty, empty])

    # ASSERT
    assert qcTable.empty
    assert qcTable.columns.tolist() == ["qc_name"]
