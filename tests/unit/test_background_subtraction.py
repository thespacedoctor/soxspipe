"""Small-array contracts for scattered-background subtraction."""

from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.subtract_background import subtract_background
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.unit
backgroundModule = importlib.import_module("soxspipe.commonutils.subtract_background")


def _frame(shape: tuple[int, int] = (32, 32), value: float = 10.0) -> CCDData:
    return CCDData(
        np.full(shape, value, dtype=np.float32),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros(shape, dtype=bool),
        uncertainty=StdDevUncertainty(np.ones(shape), unit=u.electron),
    )


def _subtractor(tmp_path: Path, log: object) -> subtract_background:
    worker = object.__new__(subtract_background)
    worker.log = log
    worker.frame = _frame()
    worker.axisA = "x"
    worker.axisB = "y"
    worker.settings = pipeline_settings(tmp_path)
    worker.settings["background-subtraction"] = {
        "bspline-deg": 3,
        "gaussian-blur-sigma": 1,
    }
    worker.kw = keyword_lookup(log=log, settings=worker.settings).get
    worker.arm = "VIS"
    worker.orderTable = str(tmp_path / "orders.fits")
    worker.sofName = False
    worker.products = False
    worker.qcDir = str(tmp_path)
    return worker


def test_mask_order_locations_masks_expanded_orders(
    tmp_path: Path, log: object
) -> None:
    worker = _subtractor(tmp_path, log)
    orderPixels = pd.DataFrame(
        {
            "order": [10, 11],
            "ycoord": [5, 5],
            "xcoord_edgeup": [24.0, 12.0],
            "xcoord_edgelow": [20.0, 8.0],
        }
    )

    worker.mask_order_locations(orderPixels)

    assert worker.frame.mask[5, 6:14].all()
    assert worker.frame.mask[5, 18:26].all()
    assert not worker.frame.mask[5, 5]
    assert not worker.frame.mask[5, 14:18].any()
    assert not worker.frame.mask[4].any()


def test_subtract_restores_input_mask_and_subtracts_background(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    worker = _subtractor(tmp_path, log)
    worker.frame = _frame(shape=(8, 8))
    worker.frame.mask[0, 0] = True
    originalMask = worker.frame.mask.copy()
    background = _frame(shape=(8, 8), value=2.0)
    orderPixels = pd.DataFrame(
        {
            "order": [10],
            "ycoord": [4],
            "xcoord_edgeup": [6.0],
            "xcoord_edgelow": [2.0],
        }
    )
    monkeypatch.setattr(
        backgroundModule,
        "unpack_order_table",
        lambda **kwargs: (pd.DataFrame(), orderPixels, pd.DataFrame()),
    )
    monkeypatch.setattr(backgroundModule, "quicklook_image", lambda **kwargs: None)
    monkeypatch.setattr(worker, "mask_order_locations", lambda table: None)
    monkeypatch.setattr(worker, "create_background_image", lambda **kwargs: background)

    backgroundFrame, subtracted, products = worker.subtract()

    assert backgroundFrame is background
    np.testing.assert_allclose(subtracted.data, 8.0, rtol=0, atol=1e-6)
    np.testing.assert_array_equal(worker.frame.mask, originalMask)
    np.testing.assert_array_equal(subtracted.mask, originalMask)
    assert products is False
