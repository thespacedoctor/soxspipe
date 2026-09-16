"""Contracts keeping reduction table ordering independent of the sort's tie-breaking.

`pandas.DataFrame.sort_values` defaults to an unstable quicksort, and NumPy
dispatches float sorts to SIMD kernels whose permutation of equal keys depends on
the CPU. Sorting a reduction table on a non-unique key therefore yields different
row orders on different machines for identical input, which changes summation
order in the least-squares fits downstream and flips lines across the
sigma-clipping threshold.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

pytestmark = pytest.mark.unit


def _line_table_with_duplicate_wavelengths() -> pd.DataFrame:
    """Return a line table whose wavelengths repeat across orders and slit positions."""
    # ARRANGE: NINE ROWS SHARING THREE WAVELENGTHS, AS A MULTI-PINHOLE FRAME DOES
    return pd.DataFrame(
        {
            "wavelength": [1100.0, 1100.0, 1100.0, 1050.0, 1050.0, 1050.0, 1200.0, 1200.0, 1200.0],
            "order": [15, 15, 15, 16, 16, 16, 14, 14, 14],
            "slit_index": [0, 1, 2, 0, 1, 2, 0, 1, 2],
            "observed_x": [10.5, 11.5, 12.5, 20.5, 21.5, 22.5, 30.5, 31.5, 32.5],
        }
    )


def test_sorting_on_wavelength_alone_leaves_tied_rows_in_an_unpinned_order():
    """A non-unique sort key does not determine the row order of tied rows."""
    # ARRANGE
    table = _line_table_with_duplicate_wavelengths()
    shuffled = table.iloc[[2, 0, 1, 5, 3, 4, 8, 6, 7]].reset_index(drop=True)

    # ACT: THE SAME ROWS AND VALUES, PRESENTED IN A DIFFERENT INCOMING ORDER
    sortedTable = table.sort_values(["wavelength"], kind="stable").reset_index(drop=True)
    sortedShuffled = shuffled.sort_values(["wavelength"], kind="stable").reset_index(drop=True)

    # ASSERT: THE KEY ALONE CANNOT PIN THE RESULT, SO THE TWO ORDERS STILL DIFFER
    assert not sortedTable["slit_index"].equals(sortedShuffled["slit_index"])


def test_sorting_on_the_full_key_pins_the_row_order_whatever_the_incoming_order():
    """A total sort key yields one row order regardless of how the rows arrive."""
    # ARRANGE
    table = _line_table_with_duplicate_wavelengths()
    shuffled = table.iloc[[2, 0, 1, 5, 3, 4, 8, 6, 7]].reset_index(drop=True)
    sortKey = ["wavelength", "order", "slit_index"]

    # ACT
    sortedTable = table.sort_values(sortKey, kind="stable").reset_index(drop=True)
    sortedShuffled = shuffled.sort_values(sortKey, kind="stable").reset_index(drop=True)

    # ASSERT
    pd.testing.assert_frame_equal(sortedTable, sortedShuffled)


def test_the_dispersion_map_sorts_its_line_table_on_a_total_key():
    """The wavelength sort in create_dispersion_map names a key that is unique per row."""
    # ARRANGE
    import inspect
    from pathlib import Path

    from soxspipe.commonutils.create_dispersion_map import create_dispersion_map

    # THE PACKAGE RE-EXPORTS THE CLASS UNDER THE MODULE'S OWN NAME, SO THE SOURCE
    # FILE IS RESOLVED FROM THE CLASS RATHER THAN FROM A MODULE ATTRIBUTE
    source = Path(inspect.getsourcefile(create_dispersion_map)).read_text(encoding="utf-8")

    # ACT: THE SORT THAT FEEDS THE POLYNOMIAL FIT
    offending = 'orderPixelTable.sort_values(["wavelength"], inplace=True)'

    # ASSERT
    assert offending not in source, (
        "sorting the line table on wavelength alone leaves tied rows unpinned; "
        "sort on wavelength, order and slit_index with kind='stable'"
    )


def test_a_stable_sort_preserves_the_incoming_order_of_tied_rows():
    """`kind='stable'` is what makes tied rows keep their incoming order."""
    # ARRANGE
    values = np.array([1100.0, 1050.0, 1100.0, 1050.0])
    frame = pd.DataFrame({"wavelength": values, "tag": ["a", "b", "c", "d"]})

    # ACT
    result = frame.sort_values(["wavelength"], kind="stable")["tag"].tolist()

    # ASSERT
    assert result == ["b", "d", "a", "c"]
