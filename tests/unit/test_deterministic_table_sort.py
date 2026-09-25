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


def _static_calibration(relativePath: str):
    """Return the path to a shipped static calibration, or None when it is absent."""
    from pathlib import Path

    import soxspipe

    candidate = Path(soxspipe.__file__).parent / "resources" / "static_calibrations" / relativePath
    return candidate if candidate.exists() else None


def _dispersion_map_source() -> str:
    """Return the on-disk source of the create_dispersion_map module."""
    import inspect
    from pathlib import Path

    from soxspipe.commonutils.create_dispersion_map import create_dispersion_map

    return Path(inspect.getsourcefile(create_dispersion_map)).read_text(encoding="utf-8")


def _calibration_shaped_line_table() -> pd.DataFrame:
    """Return the duplicate-key pair that the shipped NIR arc line list actually contains.

    `ArHgNeXe_clean_within_2.0pixel.fits` holds two rows sharing
    `(wavelength=1588.31, order=12, slit_index=0)` that differ only in their
    detector position, so those three columns are not a total key on real data.
    """
    # ARRANGE: THE TWO REAL ROWS, PLUS ONE ORDINARY ROW EITHER SIDE OF THEM
    return pd.DataFrame(
        {
            "wavelength": [1588.31, 1588.31, 1500.0, 1600.0],
            "order": [12, 12, 12, 12],
            "slit_index": [0, 0, 0, 0],
            "slit_position": [-5.74545, -5.74545, -5.74545, -5.74545],
            "detector_x": [531.281742, 534.379648, 400.0, 600.0],
            "detector_y": [465.973724, 465.659446, 300.0, 500.0],
        }
    )


def test_the_shipped_arc_line_list_repeats_wavelength_order_and_slit_index():
    """The three-column key is not total on the data the pipeline actually reduces."""
    # ARRANGE
    from astropy.table import Table

    lineList = _static_calibration("soxs/ArHgNeXe_clean_within_2.0pixel.fits")
    if lineList is None:
        pytest.skip("the shipped NIR arc line list is not installed")

    # ACT
    lines = Table.read(lineList).to_pandas()
    duplicated = lines.duplicated(subset=["wavelength", "order", "slit_index"], keep=False)

    # ASSERT: IF THIS EVER BECOMES EMPTY THE THREE-COLUMN KEY WOULD SUFFICE, BUT
    # THE SORT SHOULD STILL NAME THE DETECTOR POSITION SO IT CANNOT REGRESS
    assert duplicated.any(), "expected the shipped line list to contain a repeated (wavelength, order, slit_index)"


def test_sorting_a_calibration_shaped_tie_is_independent_of_the_incoming_order():
    """Two lines sharing wavelength, order and slit index still sort to one fixed order."""
    # ARRANGE
    table = _calibration_shaped_line_table()
    shuffled = table.iloc[[1, 0, 3, 2]].reset_index(drop=True)
    sortKey = ["wavelength", "order", "slit_index", "detector_x", "detector_y"]

    # ACT
    sortedTable = table.sort_values(sortKey, kind="stable").reset_index(drop=True)
    sortedShuffled = shuffled.sort_values(sortKey, kind="stable").reset_index(drop=True)

    # ASSERT
    pd.testing.assert_frame_equal(sortedTable, sortedShuffled)


def _sorter():
    """Return the dispersion map's sort helper bound to a stub carrying only a logger."""
    import logging
    from types import SimpleNamespace

    from soxspipe.commonutils.create_dispersion_map import create_dispersion_map

    stub = SimpleNamespace(log=logging.getLogger("test_deterministic_table_sort"))
    return lambda frame: create_dispersion_map._sort_line_table_on_a_total_key(stub, frame)


def test_the_dispersion_map_pins_the_order_of_a_calibration_shaped_tie():
    """Two lines agreeing on wavelength, order and slit index still sort to one order."""
    # ARRANGE
    sort = _sorter()
    table = _calibration_shaped_line_table()
    shuffled = table.iloc[[1, 0, 3, 2]].reset_index(drop=True)

    # ACT
    sortedTable = sort(table).reset_index(drop=True)
    sortedShuffled = sort(shuffled).reset_index(drop=True)

    # ASSERT: THE DETECTOR POSITION BREAKS THE TIE, SO BOTH ARRIVE AT THE SAME ORDER
    pd.testing.assert_frame_equal(sortedTable, sortedShuffled)


def test_the_dispersion_map_pins_the_order_of_rows_alike_on_every_physical_column():
    """Rows identical on every physical column are ordered by their incoming position."""
    # ARRANGE: TWO ROWS THE PHYSICAL KEY CANNOT SEPARATE AT ALL
    sort = _sorter()
    table = pd.DataFrame(
        {
            "wavelength": [1588.31, 1588.31, 1500.0],
            "order": [12, 12, 12],
            "slit_index": [0, 0, 0],
            "detector_x": [531.281742, 531.281742, 400.0],
            "detector_y": [465.973724, 465.973724, 300.0],
            "tag": ["first", "second", "other"],
        }
    )

    # ACT
    result = sort(table)

    # ASSERT: THE INPUT ORDER DECIDES, AND IT IS THE SAME EVERY RUN
    assert result["tag"].tolist() == ["other", "first", "second"]


def test_the_sort_helper_does_not_leave_its_tie_breaker_behind():
    """The sort adds no column to the table it returns."""
    # ARRANGE
    sort = _sorter()
    table = _calibration_shaped_line_table()

    # ACT
    result = sort(table)

    # ASSERT
    assert list(result.columns) == list(table.columns)


def test_the_dispersion_map_sort_key_names_the_detector_position():
    """The line-table sort key includes the columns that break a calibration tie."""
    # ARRANGE
    source = _dispersion_map_source()

    # ACT
    # THE KEY IS DECLARED AS A TUPLE AND THEN FILTERED AGAINST THE TABLE'S COLUMNS, SO THE
    # NAMES AND THE COMPREHENSION THAT CONSUMES THEM CAN SIT ON SEPARATE LINES
    candidateKeys = [line for line in source.splitlines() if "candidateKey = (" in line]
    filteredKeys = [line for line in source.splitlines() if "physicalKey = [" in line]

    # ASSERT
    assert candidateKeys, "expected create_dispersion_map to declare an explicit sort key"
    assert filteredKeys, "expected the sort key to be filtered against the table's columns"
    assert all("detector_x" in line and "detector_y" in line for line in candidateKeys), (
        "wavelength, order and slit_index repeat in the shipped arc line list, so the "
        "sort key must also name the detector position to be total"
    )
