"""Contracts for ``slit_illumination_mask.exclude_mask_from_result``, the sky mask of a measurement (DY-1280)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.slit_illumination import (
    FLAT_EDGE_COLUMNS,
    OFFSET_COLUMNS,
    PROFILE_COLUMNS,
    STATUS_FITTED,
    _read_disp_map,
    slit_illumination_result,
)
from soxspipe.commonutils.slit_illumination_mask import exclude_mask_from_result, overlap_exclude_mask
from tests.factories.slit_illumination import synthetic_illumination_scene

SOXS_KEYWORDS = {"SEQ_ARM": "ESO SEQ ARM", "WIN_BINX": "ESO DET BINX", "WIN_BINY": "ESO DET BINY"}
kw = SOXS_KEYWORDS.__getitem__

pytestmark = pytest.mark.unit

MARGIN_ARCSEC = 0.3
OFFSET_ARCSEC = 1.0


def _result(orders, isApplied=True):
    offsets = pd.DataFrame(
        [
            {
                "order": order,
                "offset": OFFSET_ARCSEC,
                "offsetErr": 0.02,
                "scienceWidth": 0.15,
                "nBins": 100,
                "status": STATUS_FITTED,
                "reason": "",
                "borrowedFrom": np.nan,
            }
            for order in orders
        ],
        columns=OFFSET_COLUMNS,
    )
    flatEdges = pd.DataFrame(
        [{"order": order, "lower": -5.0, "upper": 5.0, "width": 0.15, "reason": ""} for order in orders],
        columns=FLAT_EDGE_COLUMNS,
    )
    return slit_illumination_result(
        offsets=offsets,
        flatEdges=flatEdges,
        profiles=pd.DataFrame(columns=PROFILE_COLUMNS),
        isApplied=isApplied,
        reason="" if isApplied else "no order has a usable offset",
        overlapMarginArcsec=MARGIN_ARCSEC,
    )


def test_returns_none_when_the_measurement_was_not_applied(log, tmp_path):
    scene = synthetic_illumination_scene(tmp_path)

    mask = exclude_mask_from_result(log, _result(scene.orders, isApplied=False), scene.twoDMapPath, 1, 1, kw)

    assert mask is None


def test_returns_the_overlap_mask_of_the_map_images_at_the_frame_binning(log, tmp_path):
    scene = synthetic_illumination_scene(tmp_path)
    result = _result(scene.orders)
    _, slit, orderMap = _read_disp_map(scene.twoDMapPath, 2, 2, kw)
    expected = overlap_exclude_mask(log, slit, orderMap, result.flatEdges, result.offsets, MARGIN_ARCSEC)

    mask = exclude_mask_from_result(log, result, scene.twoDMapPath, 2, 2, kw)

    assert mask.shape == (scene.slit.shape[0] // 2, scene.slit.shape[1] // 2)
    assert mask.any()
    np.testing.assert_array_equal(mask, expected)


def test_keeps_a_map_already_at_the_frame_binning_at_its_own_size(log, tmp_path):
    # ARRANGE: THE SCENE'S MAP RECORDED AT 2x2 BINNING, USED WITH A 2x2-BINNED FRAME (DY-1260)
    scene = synthetic_illumination_scene(tmp_path)
    fits.setval(scene.twoDMapPath, "ESO DET BINX", value=2, ext=0)
    fits.setval(scene.twoDMapPath, "ESO DET BINY", value=2, ext=0)

    # ACT
    mask = exclude_mask_from_result(log, _result(scene.orders), scene.twoDMapPath, 2, 2, kw)

    # ASSERT
    assert mask.shape == scene.slit.shape


def test_excludes_only_the_slit_end_the_offset_leaves_unlit(log, tmp_path):
    scene = synthetic_illumination_scene(tmp_path)

    mask = exclude_mask_from_result(log, _result(scene.orders), scene.twoDMapPath, 1, 1, kw)

    # AN OFFSET OF +1 ARCSEC KEEPS [-4 + 0.3, 5 - 0.3] OF A FLAT LIT FROM -5 TO 5 ARCSEC
    keptSlit = scene.slit[~mask & np.isfinite(scene.slit)]
    assert keptSlit.min() >= -3.7 - 1e-3
    assert keptSlit.max() <= 4.7 + 1e-3
    assert not mask[~np.isfinite(scene.orderMap)].any()
