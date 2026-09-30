"""Synthetic FITS contracts for unpacking order-location tables."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.toolkit import unpack_order_table
from tests.factories import order_table_fits

pytestmark = pytest.mark.integration


def _order_table(tmp_path: Path) -> Path:
    polynomials = pd.DataFrame(
        [
            {
                "degorder_cent": 1,
                "degy_cent": 1,
                "cent_00": 1.0,
                "cent_01": 2.0,
                "cent_10": 3.0,
                "cent_11": 4.0,
                "std_00": 2.0,
                "std_01": 0.0,
                "std_10": 0.0,
                "std_11": 0.0,
            }
        ]
    )
    metadata = pd.DataFrame(
        {
            "order": [10, 11],
            "ymin": [0.0, 2.0],
            "ymax": [8.0, 10.0],
        }
    )
    return order_table_fits(
        tmp_path / "orders.fits",
        polynomials=polynomials,
        metadata=metadata,
    )


def test_unpack_order_table_evaluates_polynomial_and_rounds_even_delta(
    tmp_path: Path,
    log: object,
) -> None:
    orderPath = _order_table(tmp_path)

    polynomialTable, pixelTable, metadataTable = unpack_order_table(
        log=log,
        orderTablePath=str(orderPath),
        pixelDelta=2,
    )

    assert list(pixelTable["ycoord"]) == [0, 3, 6, 2, 5, 8]
    expectedCentres = 1 + 2 * pixelTable["ycoord"] + 3 * pixelTable["order"]
    expectedCentres += 4 * pixelTable["order"] * pixelTable["ycoord"]
    np.testing.assert_allclose(
        pixelTable["xcoord_centre"], expectedCentres, rtol=1e-12, atol=1e-12
    )
    assert polynomialTable.iloc[0]["cent_11"] == 4.0
    assert list(metadataTable["order"]) == [10, 11]


def test_unpack_order_table_filters_order_and_applies_binning(
    tmp_path: Path,
    log: object,
) -> None:
    orderPath = _order_table(tmp_path)

    _, pixelTable, metadataTable = unpack_order_table(
        log=log,
        orderTablePath=str(orderPath),
        pixelDelta=1,
        binx=2,
        biny=2,
        order=11,
    )

    assert set(pixelTable["order"]) == {11}
    assert list(pixelTable["ycoord"]) == [1, 2, 3, 4]
    unbinnedY = pixelTable["ycoord"] * 2
    expectedCentres = 1 + 2 * unbinnedY + 3 * 11 + 4 * 11 * unbinnedY
    np.testing.assert_allclose(
        pixelTable["xcoord_centre"], expectedCentres / 2, rtol=1e-12, atol=1e-12
    )
    np.testing.assert_allclose(pixelTable["std"], 1.0, rtol=1e-12, atol=1e-12)
    assert metadataTable.iloc[0]["ymin"] == 1.0
    assert metadataTable.iloc[0]["ymax"] == 5.0
