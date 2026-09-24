from types import SimpleNamespace

import numpy as np
import pandas as pd

from soxspipe.commonutils.image_transformer import image_transformer


class _Log(object):

    def debug(self, *args, **kwargs):
        pass


def test_setup_order_table_discards_out_of_bounds_slices():
    transformer = object.__new__(image_transformer)
    transformer.log = _Log()
    transformer.orderNums = [1]
    transformer.amins = [0]
    transformer.amaxs = [10]
    transformer.waveLengthMin = [100]
    transformer.waveLengthMax = [200]
    transformer.uniqueOrders = np.array([1])
    transformer.sliceOrders = []
    transformer.axisA = "x"
    transformer.axisB = "y"
    transformer.dispersionAxis = "x"
    transformer.slitHalfLength = 2
    transformer.zoomFactor = 2
    transformer.twoDMap = {"WAVELENGTH": SimpleNamespace(data=np.zeros((5, 5)))}
    transformer.orderPixelTable = pd.DataFrame(
        {
            "order": [1, 1, 1],
            "xcoord_centre": [1.0, 3.0, 3.0],
            "ycoord": [2.0, 3.0, 5.0],
        }
    )

    transformer._setup_order_table_dataframes()

    assert transformer.sliceOrders == [1]
    assert len(transformer.orderSlices) == 1
    assert len(transformer.orderSlices[0]) == 1
    assert transformer.orderSlices[0]["xcoord_centre"].tolist() == [3.0]
    assert transformer.orderSlices[0]["ycoord"].tolist() == [3.0]
    assert transformer.wlMinMax == [(100, 200)]

    row_indexes, column_indexes = transformer.subPixelIndexes[0]
    assert len(row_indexes) == 1
    assert len(column_indexes) == 1
    assert len(row_indexes[0]) == 8
    assert len(column_indexes[0]) == 8
