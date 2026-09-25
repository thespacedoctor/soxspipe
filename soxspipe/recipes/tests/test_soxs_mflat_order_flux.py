import importlib

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes.soxs_mflat import soxs_mflat

mflat_module = importlib.import_module("soxspipe.recipes.soxs_mflat")


class _Log:
    def debug(self, *args, **kwargs):
        pass

    def print(self, *args, **kwargs):
        pass


def _recipe(axis):
    recipe = object.__new__(soxs_mflat)
    recipe.log = _Log()
    recipe.axisA = axis
    recipe.axisB = "y" if axis == "x" else "x"
    recipe.binx = 1
    recipe.biny = 1
    recipe.recipeSettings = {"low-sensitivity-clipping-sigma": 5}
    return recipe


def _order_table(axis):
    if axis == "x":
        return pd.DataFrame(
            {
                "order": [10, 10],
                "xcoord_edgeup": [4, 4],
                "xcoord_edgelow": [1, 1],
                "ycoord": [2, 3],
            }
        )
    return pd.DataFrame(
        {
            "order": [10, 10],
            "ycoord_edgeup": [4, 4],
            "ycoord_edgelow": [1, 1],
            "xcoord": [2, 3],
        }
    )


def _frame(axis, mask_extreme=False, mask_all=False):
    data = np.zeros((6, 6), dtype=float)
    if axis == "x":
        data[2, 1:4] = [1, 2, 999]
        data[3, 1:4] = [3, 4, 5]
        sampled = (slice(2, 4), slice(1, 4))
        extreme = (2, 3)
    else:
        data[1:4, 2] = [1, 2, 999]
        data[1:4, 3] = [3, 4, 5]
        sampled = (slice(1, 4), slice(2, 4))
        extreme = (3, 2)

    mask = np.zeros_like(data, dtype=bool)
    if mask_extreme:
        mask[extreme] = True
    if mask_all:
        mask[sampled] = True
    return CCDData(data, mask=mask, unit=u.electron)


def _measure_median(monkeypatch, axis, frame):
    order_table = _order_table(axis)
    monkeypatch.setattr(
        mflat_module,
        "unpack_order_table",
        lambda **kwargs: (None, order_table, None),
    )
    monkeypatch.setattr(mflat_module, "quicklook_image", lambda **kwargs: None)

    _, medians = _recipe(axis).mask_low_sens_pixels(
        frame,
        "unused.fits",
        returnMedianOrderFlux=True,
        writeQC=False,
    )
    return medians.loc[0, "medianFlux"]


@pytest.mark.parametrize("axis", ["x", "y"])
def test_order_median_excludes_pre_masked_extreme_pixel(monkeypatch, axis):
    assert _measure_median(monkeypatch, axis, _frame(axis, mask_extreme=True)) == 3


@pytest.mark.parametrize("axis", ["x", "y"])
def test_order_median_is_unchanged_without_masked_samples(monkeypatch, axis):
    assert _measure_median(monkeypatch, axis, _frame(axis)) == 3.5


@pytest.mark.parametrize("axis", ["x", "y"])
def test_fully_masked_order_raises_with_order_number(monkeypatch, axis):
    with pytest.raises(ValueError, match=r"order 10: no valid sampled pixels"):
        _measure_median(monkeypatch, axis, _frame(axis, mask_all=True))
