"""Partly illuminated order-edge pixels are masked by `soxs_mflat.mask_low_sens_pixels`."""

from __future__ import annotations

import importlib
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData

from soxspipe.recipes.soxs_mflat import soxs_mflat

pytestmark = pytest.mark.unit

mflatModule = importlib.import_module("soxspipe.recipes.soxs_mflat")

# ROW LEVELS VARY SO THE FRAME-WIDE SIGMA CLIP DOES NOT FLAG THE DIMMED EDGES ITSELF
ROW_LEVELS = [6.0, 8.0, 10.0, 12.0, 14.0]
EDGE_DIM_FRACTION = 0.17
SEGMENT_LOW = 2
SEGMENT_UP = 10
SHAPE = (5, 12)
SETTING_KEY = "order-edge-min-flat-fraction"
# A DIMMED PIXEL AWAY FROM THE SEGMENT ENDS: THE RULE IS RELATIVE-LOW-FLAT, NOT POSITIONAL
INTERIOR_DIMMED = (3, 5)
INTERIOR_DIM_FRACTION = 0.3


def _recipe(log: Any, axis: str, settings: dict[str, float | None] | None = None) -> soxs_mflat:
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.axisA = axis
    recipe.axisB = "y" if axis == "x" else "x"
    recipe.binx = 1
    recipe.biny = 1
    recipe.recipeName = "soxs-mflat"
    recipe.recipeSettings = {"low-sensitivity-clipping-sigma": 5, **(settings or {})}
    recipe.qc = pd.DataFrame()
    recipe.dateObs = "2024-01-01T00:00:00"
    return recipe


def _order_table(axis: str) -> pd.DataFrame:
    rows = len(ROW_LEVELS)
    return pd.DataFrame(
        {
            "order": [10] * rows,
            f"{axis}coord_edgeup": [SEGMENT_UP] * rows,
            f"{axis}coord_edgelow": [SEGMENT_LOW] * rows,
            f"{'y' if axis == 'x' else 'x'}coord": list(range(rows)),
        }
    )


def _frame(axis: str) -> CCDData:
    """Build a plateau per order row segment with a dimmed pixel at each end and one in the interior."""
    data = np.ones(SHAPE, dtype=float)
    for b, level in enumerate(ROW_LEVELS):
        data[b, SEGMENT_LOW:SEGMENT_UP] = level
        data[b, SEGMENT_LOW] = level * EDGE_DIM_FRACTION
        data[b, SEGMENT_UP - 1] = level * EDGE_DIM_FRACTION
    row, col = INTERIOR_DIMMED
    data[row, col] = ROW_LEVELS[row] * INTERIOR_DIM_FRACTION
    if axis == "y":
        data = data.T.copy()
    return CCDData(data, mask=np.zeros(data.shape, dtype=bool), unit=u.electron)


def _run(
    monkeypatch: pytest.MonkeyPatch,
    log: Any,
    axis: str,
    frame: CCDData,
    settings: dict[str, float | None] | None = None,
) -> tuple[CCDData, soxs_mflat]:
    orderTable = _order_table(axis)
    monkeypatch.setattr(mflatModule, "unpack_order_table", lambda **kwargs: (None, orderTable, None))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    recipe = _recipe(log, axis, settings)
    return recipe.mask_low_sens_pixels(frame, "orders.fits", writeQC=True), recipe


def _expected_dimmed_mask(axis: str) -> np.ndarray:
    mask = np.zeros(SHAPE, dtype=bool)
    mask[:, SEGMENT_LOW] = True
    mask[:, SEGMENT_UP - 1] = True
    mask[INTERIOR_DIMMED] = True
    return mask if axis == "x" else mask.T


@pytest.mark.parametrize("axis", ["x", "y"])
def test_dimmed_edge_and_interior_pixels_are_masked_and_plateau_pixels_are_not(
    monkeypatch: pytest.MonkeyPatch, log: Any, axis: str
) -> None:
    # ARRANGE
    frame = _frame(axis)

    # ACT
    masked, _ = _run(monkeypatch, log, axis, frame)

    # ASSERT
    np.testing.assert_array_equal(masked.mask, _expected_dimmed_mask(axis))


@pytest.mark.parametrize("axis", ["x", "y"])
def test_pixel_above_the_edge_fraction_is_not_masked(monkeypatch: pytest.MonkeyPatch, log: Any, axis: str) -> None:
    # ARRANGE
    frame = _frame(axis)
    row, col = 2, SEGMENT_LOW
    position = (row, col) if axis == "x" else (col, row)
    frame.data[position] = ROW_LEVELS[row] * 0.6

    # ACT
    masked, _ = _run(monkeypatch, log, axis, frame)

    # ASSERT
    assert not masked.mask[position]
    assert masked.mask.sum() == _expected_dimmed_mask(axis).sum() - 1


def test_order_edge_setting_override_changes_the_threshold(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")
    frame.data[2, SEGMENT_LOW] = ROW_LEVELS[2] * 0.6

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame, {SETTING_KEY: 0.7})

    # ASSERT
    assert masked.mask[2, SEGMENT_LOW]


def test_order_edge_threshold_of_zero_disables_edge_masking(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame, {SETTING_KEY: 0.0})

    # ASSERT
    assert not masked.mask.any()


def test_edge_flagging_does_not_change_the_low_sensitivity_qc_value(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE: A SIGMA OF 0.5 MAKES THE FRAME-WIDE CLIP FLAG SOME DIMMED EDGES ON ITS OWN
    settings = {"low-sensitivity-clipping-sigma": 0.5}
    _, withoutEdges = _run(monkeypatch, log, "x", _frame("x"), {**settings, SETTING_KEY: 0.0})

    # ACT
    _, withEdges = _run(monkeypatch, log, "x", _frame("x"), settings)

    # ASSERT
    lowSensWithout = withoutEdges.qc.loc[withoutEdges.qc["qc_name"] == "N LOW SENS", "qc_value"].iloc[0]
    lowSensWith = withEdges.qc.loc[withEdges.qc["qc_name"] == "N LOW SENS", "qc_value"].iloc[0]
    assert lowSensWithout > 0
    assert lowSensWith == lowSensWithout


def test_edge_pixel_count_is_logged_beside_the_low_sensitivity_line(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")

    # ACT
    _run(monkeypatch, log, "x", frame)

    # ASSERT
    printed = [message for level, message in log.messages if level == "print"]
    assert "        11 order-edge pixels added to bad-pixel mask" in printed


def test_flat_data_values_are_unchanged_by_edge_masking(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")
    expected = frame.data.copy()

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame)

    # ASSERT
    np.testing.assert_array_equal(masked.data, expected)


def test_pixels_already_masked_are_excluded_from_the_segment_median(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE: THREE USABLE PIXELS REMAIN IN ROW 2, TWO DIMMED AND ONE PLATEAU, SO THEIR MEDIAN IS THE DIMMED VALUE
    # (INCLUDING THE MASKED PLATEAU PIXELS THE MEDIAN WOULD BE THE PLATEAU LEVEL AND THE DIMMED PIXELS WOULD BE FLAGGED)
    frame = _frame("x")
    frame.mask[2, SEGMENT_LOW + 2 : SEGMENT_UP - 1] = True

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame)

    # ASSERT
    assert not masked.mask[2, SEGMENT_LOW]
    assert not masked.mask[2, SEGMENT_UP - 1]
    assert not masked.mask[2, SEGMENT_LOW + 1]


def test_segment_with_non_positive_median_is_not_edge_masked(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")
    frame.data[2, SEGMENT_LOW:SEGMENT_UP] = -1.0

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame)

    # ASSERT
    assert not masked.mask[2].any()


def test_segment_with_too_few_usable_pixels_is_not_edge_masked(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE: LEAVE TWO USABLE PIXELS IN ROW 2, THE DIMMED EDGE AND ONE PLATEAU PIXEL
    frame = _frame("x")
    frame.mask[2, SEGMENT_LOW + 2 : SEGMENT_UP] = True

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame)

    # ASSERT
    assert not masked.mask[2, SEGMENT_LOW]


def test_zero_fraction_does_not_mask_negative_pixels(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE: A NEGATIVE PIXEL IN A SEGMENT WITH A POSITIVE MEDIAN, AS BACKGROUND SUBTRACTION CAN PRODUCE
    frame = _frame("x")
    frame.data[2, 5] = -1.0
    frame.data[INTERIOR_DIMMED] = ROW_LEVELS[INTERIOR_DIMMED[0]]

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame, {SETTING_KEY: 0})

    # ASSERT
    assert not masked.mask.any()


def test_blank_fraction_setting_falls_back_to_the_default(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    frame = _frame("x")

    # ACT
    masked, _ = _run(monkeypatch, log, "x", frame, {SETTING_KEY: None})

    # ASSERT
    np.testing.assert_array_equal(masked.mask, _expected_dimmed_mask("x"))


def test_logged_edge_count_excludes_pixels_already_counted_as_low_sensitivity(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE: A SIGMA OF 0.5 MAKES THE FRAME-WIDE CLIP FLAG SOME OF THE SAME PIXELS THE EDGE RULE FLAGS
    settings = {"low-sensitivity-clipping-sigma": 0.5}

    # ACT
    masked, recipe = _run(monkeypatch, log, "x", _frame("x"), settings)

    # ASSERT
    lowSens = int(recipe.qc.loc[recipe.qc["qc_name"] == "N LOW SENS", "qc_value"].iloc[0])
    newlyAdded = int(masked.mask.sum()) - lowSens
    printed = [message for level, message in log.messages if level == "print"]
    assert lowSens > 0
    assert newlyAdded < _expected_dimmed_mask("x").sum()
    assert f"        {newlyAdded} order-edge pixels added to bad-pixel mask" in printed


def _qc_row(recipe: soxs_mflat, qcName: str) -> pd.Series:
    return recipe.qc.loc[recipe.qc["qc_name"] == qcName].iloc[0]


@pytest.mark.parametrize("axis", ["x", "y"])
def test_edge_only_mask_is_stored_and_counted_as_n_order_edge_qc(
    monkeypatch: pytest.MonkeyPatch, log: Any, axis: str
) -> None:
    # ARRANGE
    frame = _frame(axis)

    # ACT
    _, recipe = _run(monkeypatch, log, axis, frame)

    # ASSERT
    np.testing.assert_array_equal(recipe.orderEdgeOnlyMask, _expected_dimmed_mask(axis))
    row = _qc_row(recipe, "N ORDER EDGE")
    assert row["qc_value"] == float(_expected_dimmed_mask(axis).sum())
    assert row["qc_unit"] == "pixels"
    assert row["qc_comment"] == "Number of partly illuminated order-edge pixels masked in master flat"
    assert row["obs_date_utc"] == recipe.dateObs


def test_edge_only_mask_excludes_low_sensitivity_and_already_masked_pixels(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE: A SIGMA OF 0.5 MAKES THE CLIP FLAG SOME DIMMED EDGES, AND ONE MORE IS MASKED BEFORE THE CALL
    frame = _frame("x")
    frame.mask[0, SEGMENT_LOW] = True

    # ACT
    masked, recipe = _run(monkeypatch, log, "x", frame, {"low-sensitivity-clipping-sigma": 0.5})

    # ASSERT
    lowSens = int(_qc_row(recipe, "N LOW SENS")["qc_value"])
    edgeOnly = recipe.orderEdgeOnlyMask
    assert lowSens > 0
    assert not edgeOnly[0, SEGMENT_LOW]
    assert edgeOnly.sum() < _expected_dimmed_mask("x").sum() - 1
    assert not (edgeOnly & ~masked.mask).any()
    assert int(_qc_row(recipe, "N ORDER EDGE")["qc_value"]) == edgeOnly.sum()


def test_edge_qc_is_not_recorded_when_write_qc_is_off(monkeypatch: pytest.MonkeyPatch, log: Any) -> None:
    # ARRANGE
    monkeypatch.setattr(mflatModule, "unpack_order_table", lambda **kwargs: (None, _order_table("x"), None))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    recipe = _recipe(log, "x")

    # ACT
    recipe.mask_low_sens_pixels(_frame("x"), "orders.fits", writeQC=False)

    # ASSERT
    assert recipe.qc.empty
    np.testing.assert_array_equal(recipe.orderEdgeOnlyMask, _expected_dimmed_mask("x"))


def test_coldpix_qc_excludes_the_edge_flags_but_keeps_detector_defects(
    monkeypatch: pytest.MonkeyPatch, log: Any
) -> None:
    # ARRANGE: ONE DETECTOR DEFECT OUTSIDE THE ORDER SEGMENTS, PLUS THE DIMMED EDGE PIXELS
    from soxspipe.commonutils import toolkit

    monkeypatch.setattr(toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key))
    frame = _frame("x")
    frame.mask[0, 0] = True
    frame.header["SEQ_ARM"] = "VIS"
    frame.header["DATE_OBS"] = "2024-01-01T00:00:00"
    masked, recipe = _run(monkeypatch, log, "x", frame)

    # ACT
    counted = toolkit.generic_quality_checks(
        log, masked, {}, "soxs-mflat", pd.DataFrame(), excludeMask=recipe.orderEdgeOnlyMask
    )

    # ASSERT
    assert masked.mask.sum() == _expected_dimmed_mask("x").sum() + 1
    assert counted["qc_value"].tolist() == [1, round(1 / masked.mask.size, 6)]
