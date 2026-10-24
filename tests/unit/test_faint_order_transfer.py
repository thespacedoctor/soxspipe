"""Transferring the trace of a bright VIS order to its faint neighbour on a synthetic three-source scene (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.faint_order_trace import faint_trace_settings, make_pixels_at, transfer_faint_orders
from tests.unit._continuum_vis_scene import (
    VIS_FLANK_SEPARATION_ARCSEC,
    VIS_SEEING_SIGMA_PX,
    vis_dispersion_map,
    vis_slit_scale,
    vis_target_slit,
    vis_transfer_image,
    vis_transfer_table,
    vis_true_centre,
)

pytestmark = pytest.mark.unit

SETTINGS = faint_trace_settings.from_block({"faint-order-trace": True})
BRIGHT = 100.0
FAINT_PEAK = 0.5
FAINT_FLANK = 0.8
AMPLITUDES = {1.0: BRIGHT, 2.0: FAINT_PEAK, 3.0: BRIGHT, 4.0: BRIGHT}
FLANKS = {1.0: 30.0, 2.0: FAINT_FLANK, 3.0: 30.0, 4.0: 30.0}


class _SilentLog:
    def debug(self, *args: object, **kwargs: object) -> None:
        return None


@pytest.fixture(scope="module", name="pixelsAt")
def vis_pixels_at(tmp_path_factory: pytest.TempPathFactory) -> object:
    mapPath = vis_dispersion_map(tmp_path_factory.mktemp("vis_map"))
    return make_pixels_at(_SilentLog(), mapPath, {"x": 1.0, "y": 1.0}, "x", "y")


@pytest.fixture(scope="module", name="faintUScene")
def faint_u_scene() -> tuple:
    """The u order at a per-pixel SNR of 0.5 with flanking sources brighter than the target, and bright g, r, i."""
    return vis_transfer_image(AMPLITUDES, FLANKS), vis_transfer_table((2.0,))


def _transfer(scene: tuple, pixelsAt: object, offsets: dict | None = None) -> tuple:
    image, table = scene
    return transfer_faint_orders(table, image, pixelsAt, offsets, SETTINGS, "x", "y")


def test_transfer_puts_the_faint_trace_within_a_third_of_a_pixel_of_the_truth(
    faintUScene: tuple, pixelsAt: object
) -> None:
    # ARRANGE: THE TARGET'S PEAK IS HALF THE NOISE, AND A SOURCE 3.7 ARCSEC AWAY IS BRIGHTER THAN IT
    _, table = faintUScene

    # ACT
    result, reports = _transfer(faintUScene, pixelsAt)

    # ASSERT
    uRows = result[result["order"] == 2.0]
    truth = vis_true_centre(2.0, uRows["wavelength"].to_numpy())
    error = uRows["cont_x"].to_numpy() - truth
    flankDistancePx = abs(vis_slit_scale(2.0)) * VIS_FLANK_SEPARATION_ARCSEC
    assert np.abs(error).max() < 0.3
    assert np.abs(error).max() < flankDistancePx / 10
    assert len(reports) == 1
    assert (reports[0].order, reports[0].anchorOrder) == (2.0, 3.0)
    assert reports[0].binsKept >= 8
    assert reports[0].isAnchorTrack is False
    assert len(result) == len(table)


def test_transfer_replaces_the_faint_rows_one_for_one_and_marks_them(faintUScene: tuple, pixelsAt: object) -> None:
    # ARRANGE
    _, table = faintUScene
    wavelength = table.loc[table["order"] == 2.0, "wavelength"].to_numpy()

    # ACT
    result, _ = _transfer(faintUScene, pixelsAt)

    # ASSERT: SAME ROWS, SAME WAVELENGTHS AND MAP POSITIONS, AND THE SYNTHETIC ROWS ARE EXACT AND MARKED
    uRows = result[result["order"] == 2.0]
    np.testing.assert_array_equal(uRows["wavelength"], wavelength)
    pd.testing.assert_series_equal(uRows["fit_x"], table.loc[table["order"] == 2.0, "fit_x"])
    assert result["faint_transfer"].tolist() == (result["order"] == 2.0).tolist()
    assert not uRows["pre-clipped"].any()
    assert uRows["gauss_amplitude"].isna().all()
    np.testing.assert_allclose(uRows["centre_shift"], uRows["fit_x"] - uRows["cont_x"])
    assert uRows[["cont_x", "cont_y", "gauss_stddev"]].notna().all().all()
    assert (uRows["gauss_stddev"] > 0.5 * VIS_SEEING_SIGMA_PX).all()


def test_transfer_leaves_the_bright_orders_and_the_input_table_untouched(faintUScene: tuple, pixelsAt: object) -> None:
    # ARRANGE
    _, table = faintUScene
    before = table.copy()

    # ACT
    result, _ = _transfer(faintUScene, pixelsAt)

    # ASSERT
    pd.testing.assert_frame_equal(table, before)
    brightRows = result["order"] != 2.0
    pd.testing.assert_frame_equal(result.loc[brightRows].drop(columns="faint_transfer"), before.loc[brightRows])


def test_transfer_returns_the_same_table_when_no_order_is_faint(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table(())
    image = vis_transfer_image(AMPLITUDES | {2.0: BRIGHT}, FLANKS)

    # ACT
    result, reports = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ASSERT
    assert result is table
    assert reports == []
    assert "faint_transfer" not in result.columns


def test_transfer_returns_the_same_table_when_both_orders_of_a_group_are_faint(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table((2.0, 3.0))
    image = vis_transfer_image(AMPLITUDES, FLANKS)

    # ACT
    result, reports = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ASSERT: U AND G ARE BOTH FAINT, SO THEY HAVE NO ANCHOR, AND NOTHING IS CARRIED ACROSS FROM R
    assert result is table
    assert reports == []


def test_transfer_anchors_i_on_r(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table((4.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: BRIGHT, 4.0: FAINT_PEAK}, FLANKS | {4.0: FAINT_FLANK})

    # ACT
    result, reports = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ASSERT
    iRows = result[result["order"] == 4.0]
    error = iRows["cont_x"].to_numpy() - vis_true_centre(4.0, iRows["wavelength"].to_numpy())
    assert (reports[0].order, reports[0].anchorOrder) == (4.0, 1.0)
    assert np.abs(error).max() < 0.5


def test_transfer_falls_back_to_the_anchor_track_when_no_bin_has_a_signal(pixelsAt: object) -> None:
    # ARRANGE: THE FAINT ORDER HOLDS NOTHING BUT NOISE
    table = vis_transfer_table((2.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0})

    # ACT
    result, reports = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ASSERT: EVERY ROW SITS ON THE SAME SLIT POSITION AS THE BRIGHT NEIGHBOUR AT THE SHARED EDGE
    uRows = result[result["order"] == 2.0]
    slit = (uRows["cont_x"] - uRows["fit_x"]) / vis_slit_scale(2.0)
    assert reports[0].binsKept == 0
    assert reports[0].isAnchorTrack is True
    assert slit.std() < 1e-3
    assert slit.iloc[0] == pytest.approx(vis_target_slit(600.0), abs=0.05)


def test_transfer_moves_the_anchor_by_the_difference_of_the_slit_offsets(pixelsAt: object) -> None:
    # ARRANGE: THE SKY SAYS THE U SLIT IMAGE IS 0.5 ARCSEC ABOVE THE G ONE, WITH NOTHING IN U TO FIT
    table = vis_transfer_table((2.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0})
    plain, _ = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ACT
    shifted, _ = transfer_faint_orders(table, image, pixelsAt, {2.0: 0.5, 3.0: 0.0}, SETTINGS, "x", "y")

    # ASSERT: ANCHOR = S_BRIGHT - DELTA_BRIGHT + DELTA_FAINT, SO THE FAINT TRACK MOVES 0.5 ARCSEC
    uPlain = plain.loc[plain["order"] == 2.0, "cont_x"].to_numpy()
    uShifted = shifted.loc[shifted["order"] == 2.0, "cont_x"].to_numpy()
    np.testing.assert_allclose(uShifted - uPlain, 0.5 * vis_slit_scale(2.0), atol=0.01)


def test_transfer_ignores_offsets_that_miss_one_of_the_two_orders(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table((2.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0})
    plain, _ = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ACT
    partial, _ = transfer_faint_orders(table, image, pixelsAt, {2.0: 0.5}, SETTINGS, "x", "y")

    # ASSERT
    pd.testing.assert_frame_equal(partial, plain)


def test_transfer_holds_the_trace_flat_beyond_the_last_kept_bin(pixelsAt: object) -> None:
    # ARRANGE: THE TARGET VANISHES IN THE SECOND HALF OF THE U ORDER, SO NO BIN THERE CAN BE KEPT
    image, table = vis_transfer_image(AMPLITUDES, FLANKS), vis_transfer_table((2.0,))
    noise = np.random.default_rng(3).normal(0.0, 1.0, (1000, 120))
    image[1010:2010, 90:210] = noise[:1000]

    # ACT
    result, reports = transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")

    # ASSERT: THE POLYNOMIAL IS NOT EXTRAPOLATED, SO THE LAST HALF OF THE ORDER FOLLOWS ONE CONSTANT OFFSET
    uRows = result[result["order"] == 2.0]
    lastHalf = uRows.loc[uRows["fit_y"] > 1300, "cont_x"]
    assert 4 <= reports[0].binsKept < reports[0].binCount
    assert np.ptp(lastHalf.to_numpy()) < 1e-6


def test_transfer_rows_are_refitted_exactly_by_the_global_polynomial(
    faintUScene: tuple, pixelsAt: object, log: object
) -> None:
    # ARRANGE: A DEGREE-3 FIT TO THE u ROWS ALONE, WITH THE SAME CLIPPING SETTINGS AS THE SHIPPED STARE BLOCK
    result, _ = _transfer(faintUScene, pixelsAt)
    detector = object.__new__(detect_continuum)
    detector.log, detector.arm, detector.axisA, detector.axisB = log, "VIS", "x", "y"
    detector.orderDeg, detector.axisBDeg = 0, 3
    detector.recipeSettings = {"poly-fitting-residual-clipping-sigma": 4, "poly-clipping-iteration-limit": 5}
    uRows = result[result["order"] == 2.0]

    # ACT
    _, fitted, clipped = detector.fit_global_polynomial(uRows)

    # ASSERT: NOTHING IS CLIPPED, AND THE RESIDUALS ARE ROUNDING ERROR
    assert len(clipped) == 0
    assert len(fitted) == len(uRows)
    assert np.abs(fitted["cont_x_fit_res"]).max() < 1e-3
