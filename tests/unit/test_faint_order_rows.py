"""The one strict mask of synthetic rows that the fit, the statistics and the plot share (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_rows import transferred_row_mask

pytestmark = pytest.mark.unit


def test_transferred_row_mask_is_none_for_a_table_without_the_column() -> None:
    # ARRANGE
    table = pd.DataFrame({"order": [1.0, 2.0]})

    # ACT / ASSERT
    assert transferred_row_mask(table) is None


def test_transferred_row_mask_follows_a_bool_column() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": [False, True, True]})

    # ACT
    mask = transferred_row_mask(table)

    # ASSERT
    assert mask.dtype == bool
    assert mask.tolist() == [False, True, True]


@pytest.mark.parametrize(
    "column",
    [
        pd.Series([np.nan, True, None, False], dtype="object"),
        pd.Series([np.nan, 1.0, np.nan, 0.0], dtype="float64"),
        pd.Series([pd.NA, True, pd.NA, False], dtype="boolean"),
    ],
    ids=["object", "float", "nullable_boolean"],
)
def test_transferred_row_mask_reads_a_missing_value_as_not_transferred(column: pd.Series) -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": column})

    # ACT
    mask = transferred_row_mask(table)

    # ASSERT: ONLY THE ROW MARKED TRUE IS SYNTHETIC. A NAN OR NONE IS NOT, WHICH `astype(bool)` WOULD GET WRONG
    assert mask.dtype == bool
    assert mask.tolist() == [False, True, False, False]


def test_transferred_row_mask_does_not_read_a_string_as_true() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": pd.Series(["False", True, "x"], dtype="object")})

    # ACT / ASSERT
    assert transferred_row_mask(table).tolist() == [False, True, False]


def test_transferred_row_mask_does_not_read_the_string_true_as_true() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": pd.Series(["True", "False", "true"], dtype="object")})

    # ACT / ASSERT
    assert transferred_row_mask(table).tolist() == [False, False, False]


def test_transferred_row_mask_reads_the_integer_one_as_set_and_every_other_integer_as_not_set() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": [0, 1, 2, -1]})

    # ACT
    mask = transferred_row_mask(table)

    # ASSERT: THE FLAG IS SET WHEN IT EQUALS TRUE, WHICH INCLUDES 1
    assert mask.dtype == bool
    assert mask.tolist() == [False, True, False, False]


def test_transferred_row_mask_gives_a_bool_mask_for_a_string_dtype_column() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": pd.Series(["True", "False", pd.NA], dtype="string")})

    # ACT
    mask = transferred_row_mask(table)

    # ASSERT
    assert mask.dtype == bool
    assert mask.tolist() == [False, False, False]


def test_transferred_row_mask_is_an_empty_bool_array_for_an_empty_table_with_the_column() -> None:
    # ARRANGE
    table = pd.DataFrame({"faint_transfer": pd.Series([], dtype=bool)})

    # ACT
    mask = transferred_row_mask(table)

    # ASSERT
    assert mask.dtype == bool
    assert mask.shape == (0,)
