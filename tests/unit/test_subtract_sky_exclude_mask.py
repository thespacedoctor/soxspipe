"""Pixels passed to ``subtract_sky`` in ``excludeMask`` stay out of the sky model and come out flagged (DY-1280).

Pixels where the master flat is low (DY-1517) are excluded the same way.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")

FRAME_SHAPE = (3, 3)
# ONE MAP ROW PER PIXEL OF THE TOP TWO DETECTOR ROWS; ORDER 12 IS DROPPED BY THE NIR JH BLOCKING FILTER
MAP_TABLE = pd.DataFrame(
    {
        "order": [12, 12, 12, 13, 13, 13],
        "x": [0, 1, 2, 0, 1, 2],
        "y": [0, 0, 0, 1, 1, 1],
        "wavelength": [1500.0, 1501.0, 1502.0, 1600.0, 1601.0, 1602.0],
        "slit_position": [-5.0, 0.0, 5.0, -5.0, 0.0, 5.0],
    }
)
# DETECTOR (ROW, COLUMN) PIXELS TO EXCLUDE: ONE IN EACH ORDER
EXCLUDED_PIXELS = [(0, 2), (1, 0)]


@pytest.fixture(autouse=True)
def collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Replace the file-reading and plotting collaborators; the keyword and detector lookups stay real."""
    monkeypatch.setattr(
        subtractSkyModule,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (MAP_TABLE.copy(), np.zeros(FRAME_SHAPE, dtype=bool)),
    )
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))


def _frame(slitName: str = "SLIT1.0") -> CCDData:
    header = instrument_header(arm="VIS", overrides={"ESO INS VISE NAME": slitName})
    return CCDData(
        np.full(FRAME_SHAPE, 10.0),
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(FRAME_SHAPE, 2.0)),
        meta=header,
    )


def _exclude_mask() -> np.ndarray:
    excludeMask = np.zeros(FRAME_SHAPE, dtype=bool)
    for row, column in EXCLUDED_PIXELS:
        excludeMask[row, column] = True
    return excludeMask


def _construct(log: Any, frame: CCDData, **kwargs: Any) -> subtract_sky:
    return subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={
            "sky-subtraction": {
                "bspline_order": 3,
                "clip-slit-edge-fraction": 0.1,
                "aggressive_object_masking": False,
                "sky_model_qc_plot": False,
            }
        },
        objectFrame=frame,
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        **kwargs,
    )


def _map_pixels(subtractor: subtract_sky) -> set[tuple[int, int]]:
    return set(zip(subtractor.mapDF["y"].tolist(), subtractor.mapDF["x"].tolist(), strict=True))


def test_without_an_exclude_mask_every_map_pixel_is_kept(log: Any) -> None:
    subtractor = _construct(log, _frame())

    assert len(subtractor.mapDF) == len(MAP_TABLE)


def test_the_excluded_pixels_are_dropped_from_the_map_and_no_others(log: Any) -> None:
    subtractor = _construct(log, _frame(), excludeMask=_exclude_mask())

    assert len(subtractor.mapDF) == len(MAP_TABLE) - 2
    assert _map_pixels(subtractor) == {(0, 0), (0, 1), (1, 1), (1, 2)}


def test_the_exclude_mask_is_applied_after_the_jh_order_filter(log: Any) -> None:
    header = instrument_header(arm="NIR", overrides={"ESO INS NISE NAME": "JH0.9"})
    frame = CCDData(np.full(FRAME_SHAPE, 10.0), unit=u.electron, meta=header)

    subtractor = _construct(log, frame, excludeMask=_exclude_mask())

    assert _map_pixels(subtractor) == {(1, 1), (1, 2)}


def test_the_caller_exclude_mask_is_not_modified(log: Any) -> None:
    excludeMask = _exclude_mask()

    _construct(log, _frame(), excludeMask=excludeMask)

    np.testing.assert_array_equal(excludeMask, _exclude_mask())


def test_an_exclude_mask_of_the_wrong_shape_raises(log: Any) -> None:
    with pytest.raises(ValueError, match="shape"):
        _construct(log, _frame(), excludeMask=np.zeros((2, 3), dtype=bool))


def test_an_exclude_mask_that_excludes_every_map_pixel_raises(log: Any) -> None:
    with pytest.raises(ValueError, match="every pixel"):
        _construct(log, _frame(), excludeMask=np.ones(FRAME_SHAPE, dtype=bool))


def test_an_excluded_pixel_comes_out_flagged_and_zero_with_its_uncertainty_kept(log: Any) -> None:
    subtractor = _construct(log, _frame(), excludeMask=_exclude_mask())
    # THE FIT STAGES RETURN THE PIXELS THEY ARE GIVEN WITH A FIXED SKY, SO ONLY THE MAP DECIDES WHAT IS MODELLED
    subtractor.get_over_sampled_sky_from_order = lambda frame, **_: frame.assign(sky_model=4.0, sky_subtracted_flux=6.0)
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = lambda frame: (
        frame.assign(error=2.0),
        (None, None, 3),
        np.array([1500.0]),
        np.array([1.0, 2.0]),
        0.5,
    )

    model, subtracted, _, _, _ = subtractor.subtract()

    for pixel in [(0, 2), (1, 0)]:
        assert model.data[pixel] == 0
        assert subtracted.data[pixel] == 0
        assert model.mask[pixel]
        assert subtracted.mask[pixel]
        assert subtracted.uncertainty.array[pixel] == 2.0
    for pixel in [(0, 0), (0, 1), (1, 1), (1, 2)]:
        assert model.data[pixel] == 4.0
        assert not subtracted.mask[pixel]


def _master_flat(rows: list[list[float]]) -> CCDData:
    return CCDData(np.array(rows, dtype=float), unit=u.dimensionless_unscaled)


GOOD_FLAT_ROWS = [[1.0, 1.0, 1.0], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]]


def test_without_a_master_flat_the_map_has_no_flat_column_and_nothing_is_excluded(log: Any) -> None:
    subtractor = _construct(log, _frame(), masterFlat=None)

    assert "flat" not in subtractor.mapDF.columns
    assert len(subtractor.mapDF) == len(MAP_TABLE)
    assert not [message for level, message in log.messages if "low-flat" in message]


def test_the_map_carries_the_master_flat_value_of_each_pixel(log: Any) -> None:
    flat = _master_flat([[1.0, 0.9, 0.8], [0.7, 0.6, 0.5], [9.0, 9.0, 9.0]])

    subtractor = _construct(log, _frame(), masterFlat=flat)

    expected = {(0, 0): 1.0, (0, 1): 0.9, (0, 2): 0.8, (1, 0): 0.7, (1, 1): 0.6, (1, 2): 0.5}
    actual = {
        (y, x): value
        for y, x, value in zip(subtractor.mapDF["y"], subtractor.mapDF["x"], subtractor.mapDF["flat"], strict=True)
    }
    assert actual == expected


def test_a_pixel_below_the_low_flat_fraction_of_its_order_is_excluded(log: Any) -> None:
    # ORDER 12 IS ROW 0, ORDER 13 IS ROW 1; 0.01 IS BELOW 2% OF THE 95TH PERCENTILE OF ORDER 12 (ABOUT 1.0)
    flat = _master_flat([[1.0, 1.0, 0.01], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

    subtractor = _construct(log, _frame(), masterFlat=flat)

    assert _map_pixels(subtractor) == {(0, 0), (0, 1), (1, 0), (1, 1), (1, 2)}


def test_the_low_flat_test_is_per_order_so_a_uniformly_low_order_is_kept(log: Any) -> None:
    # ORDER 13 IS A THOUSAND TIMES FAINTER THAN ORDER 12 BUT IS NOT LOW RELATIVE TO ITSELF
    flat = _master_flat([[1.0, 1.0, 1.0], [0.001, 0.001, 0.001], [1.0, 1.0, 1.0]])

    subtractor = _construct(log, _frame(), masterFlat=flat)

    assert len(subtractor.mapDF) == len(MAP_TABLE)


def test_a_flat_that_is_not_finite_and_positive_is_excluded(log: Any) -> None:
    flat = _master_flat([[1.0, np.nan, np.inf], [0.0, -2.0, 1.0], [1.0, 1.0, 1.0]])

    subtractor = _construct(log, _frame(), masterFlat=flat)

    assert _map_pixels(subtractor) == {(0, 0), (1, 2)}


def test_low_flat_pixels_are_excluded_as_well_as_those_of_the_exclude_mask(log: Any) -> None:
    flat = _master_flat([[1.0, 1.0, 1.0], [0.001, 1.0, 1.0], [1.0, 1.0, 1.0]])
    excludeMask = np.zeros(FRAME_SHAPE, dtype=bool)
    excludeMask[0, 2] = True

    subtractor = _construct(log, _frame(), masterFlat=flat, excludeMask=excludeMask)

    assert _map_pixels(subtractor) == {(0, 0), (0, 1), (1, 1), (1, 2)}
    assert excludeMask.sum() == 1


def test_a_master_flat_of_the_wrong_shape_raises(log: Any) -> None:
    with pytest.raises(ValueError, match="master flat"):
        _construct(log, _frame(), masterFlat=CCDData(np.ones((2, 3)), unit=u.dimensionless_unscaled))


def test_an_exclude_mask_of_the_wrong_shape_still_raises_with_a_master_flat(log: Any) -> None:
    with pytest.raises(ValueError, match="exclude mask shape"):
        _construct(
            log,
            _frame(),
            masterFlat=_master_flat(GOOD_FLAT_ROWS),
            excludeMask=np.zeros((2, 3), dtype=bool),
        )


def test_the_low_flat_pixel_count_is_logged_once_and_the_exclusion_message_is_generic(log: Any) -> None:
    flat = _master_flat([[1.0, 1.0, 0.01], [1.0, 1.0, 1.0], [1.0, 1.0, 1.0]])

    _construct(log, _frame(), masterFlat=flat)

    messages = [message for _, message in log.messages]
    assert len([message for message in messages if "low-flat" in message]) == 1
    assert any("Excluding 1 low-flat pixels" in message for message in messages)
    assert any("Excluding 1 pixels from the sky model" in message for message in messages)
    assert not any("slit-end" in message for message in messages)
