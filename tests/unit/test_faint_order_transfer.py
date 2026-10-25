"""Transferring the trace of a bright VIS order to its faint neighbour on a synthetic three-source scene (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.faint_order_map import make_pixels_at
from soxspipe.commonutils.faint_order_trace import faint_trace_settings, transfer_faint_orders
from tests.unit._continuum_vis_scene import (
    VIS_FLANK_SEPARATION_ARCSEC,
    VIS_SEEING_SIGMA_PX,
    vis_dispersion_map,
    vis_slit_scale,
    vis_target_slit,
    vis_transfer_image,
    vis_transfer_table,
    vis_true_centre,
    vis_wavelength,
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


def _transfer_rows(table: pd.DataFrame, image: np.ndarray, pixelsAt: object, offsets: dict | None = None) -> tuple:
    return transfer_faint_orders(table, image, pixelsAt, slitOffsets=offsets, settings=SETTINGS, axisA="x", axisB="y")


def _transfer(scene: tuple, pixelsAt: object, offsets: dict | None = None) -> tuple:
    image, table = scene
    return _transfer_rows(table, image, pixelsAt, offsets)


def test_transfer_takes_the_offsets_settings_and_axes_by_keyword_only(faintUScene: tuple, pixelsAt: object) -> None:
    # ARRANGE
    image, table = faintUScene

    # ACT / ASSERT: SEVEN POSITIONAL ARGUMENTS INVITE A SWAPPED PAIR OF AXES
    with pytest.raises(TypeError):
        transfer_faint_orders(table, image, pixelsAt, None, SETTINGS, "x", "y")


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
    assert reports[0].is_transferred is True
    assert len(result) == len(table)


@pytest.mark.parametrize("fitOffsetArcsec", [0.5, 1.5], ids=["half_arcsec", "one_and_half_arcsec"])
def test_transfer_anchors_on_the_true_source_when_fit_x_sits_off_the_map_slit_zero(
    fitOffsetArcsec: float, pixelsAt: object
) -> None:
    # ARRANGE: THE SECOND SAMPLING PASS LEFT FIT_X OF EVERY ORDER OFF THE MAP'S SLIT ZERO, SO CENTRE_SHIFT IS NOT
    # MEASURED FROM IT. THE G TRACE STARTS AT 600 NM, THE EDGE IT SHARES WITH U
    image = vis_transfer_image(AMPLITUDES, FLANKS)
    table = vis_transfer_table((2.0,), fitOffsetArcsec=fitOffsetArcsec)

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: THE ANCHOR IS THE TARGET'S SLIT POSITION, AND THE U TRACE LANDS ON THE TARGET
    assert reports[0].slitArcsec == pytest.approx(float(vis_target_slit(600.0)), abs=0.05)
    uRows = result[result["order"] == 2.0]
    error = uRows["cont_x"].to_numpy() - vis_true_centre(2.0, uRows["wavelength"].to_numpy())
    assert np.abs(error).max() < 0.3


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
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT
    assert result is table
    assert reports == []
    assert "faint_transfer" not in result.columns


def test_transfer_returns_the_same_table_when_both_orders_of_a_group_are_faint(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table((2.0, 3.0))
    image = vis_transfer_image(AMPLITUDES, FLANKS)

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: U AND G ARE BOTH FAINT, SO THEY HAVE NO ANCHOR, AND NOTHING IS CARRIED ACROSS FROM R
    assert result is table
    assert reports == []


def _noisy_g_table(scatterPx: float) -> pd.DataFrame:
    """Return the table with u faint and g's centres scattered by `scatterPx`, as in VV4580Sgr."""
    table = vis_transfer_table((2.0,))
    rng = np.random.default_rng(4580)
    isG = (table["order"] == 3.0).to_numpy()
    noise = rng.normal(0.0, scatterPx, int(isG.sum()))
    table.loc[isG, "cont_x"] = table.loc[isG, "cont_x"] + noise
    table.loc[isG, "centre_shift"] = table.loc[isG, "centre_shift"] - noise
    return table


def test_transfer_anchors_a_faint_u_on_a_faint_but_usable_g_and_leaves_g_untouched(pixelsAt: object) -> None:
    # ARRANGE: G IS DETECTED EVERYWHERE BUT SCATTERS BY 2.5 PX, WHICH IS FAINT (OVER 1.5) AND STILL AN ANCHOR (UNDER 3)
    table = _noisy_g_table(2.5)
    image = vis_transfer_image(AMPLITUDES, FLANKS)

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: U IS TRANSFERRED FROM G, AND G, R AND I KEEP THEIR OWN SAMPLED ROWS
    assert [(report.order, report.anchorOrder) for report in reports] == [(2.0, 3.0)]
    kept = result["order"] != 2.0
    pd.testing.assert_frame_equal(result.loc[kept].drop(columns="faint_transfer"), table.loc[kept])
    assert result["faint_transfer"].tolist() == (result["order"] == 2.0).tolist()


def test_transfer_gives_the_same_table_when_both_orders_are_over_the_anchor_limit(pixelsAt: object) -> None:
    # ARRANGE: G SCATTERS BY 3.5 PX, OVER THE 3 PX ANCHOR LIMIT, AND U IS FAINT TOO
    table = _noisy_g_table(3.5)
    image = vis_transfer_image(AMPLITUDES, FLANKS)

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT
    assert result is table
    assert reports == []


def test_transfer_anchors_i_on_r(pixelsAt: object) -> None:
    # ARRANGE
    table = vis_transfer_table((4.0,))
    image = vis_transfer_image(AMPLITUDES | {2.0: BRIGHT, 4.0: FAINT_PEAK}, FLANKS | {4.0: FAINT_FLANK})

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT
    iRows = result[result["order"] == 4.0]
    error = iRows["cont_x"].to_numpy() - vis_true_centre(4.0, iRows["wavelength"].to_numpy())
    assert (reports[0].order, reports[0].anchorOrder) == (4.0, 1.0)
    assert np.abs(error).max() < 0.5


def test_transfer_leaves_the_faint_rows_as_measured_when_no_bin_has_a_signal(pixelsAt: object) -> None:
    # ARRANGE: THE FAINT ORDER HOLDS NOTHING BUT NOISE
    table = vis_transfer_table((2.0,))
    before = table.copy()
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0})

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: NO BIN WAS MEASURED, SO THERE IS NOTHING TO TRANSFER. THE SAME TABLE COMES BACK, WITH NO MARKER COLUMN
    assert result is table
    pd.testing.assert_frame_equal(table, before)
    assert "faint_transfer" not in result.columns
    assert [(report.order, report.anchorOrder, report.binsKept) for report in reports] == [(2.0, 3.0, 0)]
    assert reports[0].is_transferred is False


def test_transfer_leaves_a_blank_order_alone_and_still_transfers_the_other_faint_order(pixelsAt: object) -> None:
    # ARRANGE: U HOLDS NOTHING, AND I IS FAINT BUT HAS A SOURCE
    table = vis_transfer_table((2.0, 4.0))
    before = table.loc[table["order"] == 2.0]
    image = vis_transfer_image(AMPLITUDES | {2.0: 0.0, 4.0: FAINT_PEAK}, FLANKS | {2.0: 0.0, 4.0: FAINT_FLANK})

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: ONLY I IS MARKED, AND THE U ROWS ARE THE MEASURED ONES
    assert [(report.order, report.is_transferred) for report in reports] == [(2.0, False), (4.0, True)]
    assert result["faint_transfer"].tolist() == (result["order"] == 4.0).tolist()
    pd.testing.assert_frame_equal(result.loc[result["order"] == 2.0].drop(columns="faint_transfer"), before)


def test_transfer_still_transfers_an_order_that_keeps_one_bin_as_a_constant_offset(pixelsAt: object) -> None:
    # ARRANGE: THE TARGET SURVIVES IN THE FIRST BIN OF THE U ORDER ONLY
    image, table = vis_transfer_image(AMPLITUDES | {2.0: 0.0}, FLANKS | {2.0: 0.0}), vis_transfer_table((2.0,))
    rows = np.arange(10, 130)
    centre = vis_true_centre(2.0, vis_wavelength(2.0, rows))[:, np.newaxis]
    columns = np.arange(image.shape[1], dtype=float)
    image[rows] += 30.0 * np.exp(-0.5 * ((columns[np.newaxis, :] - centre) / VIS_SEEING_SIGMA_PX) ** 2)

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: AT LEAST ONE BIN WAS KEPT, THE CURVE IS A CONSTANT, AND THE ORDER IS MARKED
    assert 1 <= reports[0].binsKept <= 3
    assert reports[0].degree == 0
    assert reports[0].is_transferred is True
    assert result.loc[result["order"] == 2.0, "faint_transfer"].all()


def test_transfer_moves_the_anchor_by_the_difference_of_the_slit_offsets(pixelsAt: object) -> None:
    # ARRANGE: THE SKY SAYS THE U SLIT IMAGE IS 0.5 ARCSEC ABOVE THE G ONE
    image, table = vis_transfer_image(AMPLITUDES, FLANKS), vis_transfer_table((2.0,))
    _, plain = _transfer_rows(table, image, pixelsAt)

    # ACT
    _, shifted = _transfer_rows(table, image, pixelsAt, {2.0: 0.5, 3.0: 0.0})

    # ASSERT: ANCHOR = S_BRIGHT - DELTA_BRIGHT + DELTA_FAINT
    assert shifted[0].slitArcsec - plain[0].slitArcsec == pytest.approx(0.5, abs=1e-6)


def test_transfer_ignores_offsets_that_miss_one_of_the_two_orders(pixelsAt: object) -> None:
    # ARRANGE
    image, table = vis_transfer_image(AMPLITUDES, FLANKS), vis_transfer_table((2.0,))
    _, plain = _transfer_rows(table, image, pixelsAt)

    # ACT
    _, partial = _transfer_rows(table, image, pixelsAt, {2.0: 0.5})

    # ASSERT
    assert partial[0].slitArcsec == plain[0].slitArcsec


def test_transfer_holds_the_trace_flat_beyond_the_last_kept_bin(pixelsAt: object) -> None:
    # ARRANGE: THE TARGET VANISHES IN THE SECOND HALF OF THE U ORDER, SO NO BIN THERE CAN BE KEPT
    image, table = vis_transfer_image(AMPLITUDES, FLANKS), vis_transfer_table((2.0,))
    noise = np.random.default_rng(3).normal(0.0, 1.0, (1000, 120))
    image[1010:2010, 90:210] = noise[:1000]

    # ACT
    result, reports = _transfer_rows(table, image, pixelsAt)

    # ASSERT: THE POLYNOMIAL IS NOT EXTRAPOLATED, SO THE LAST HALF OF THE ORDER FOLLOWS ONE CONSTANT OFFSET
    uRows = result[result["order"] == 2.0]
    lastHalf = uRows.loc[uRows["fit_y"] > 1300, "cont_x"]
    assert 4 <= reports[0].binsKept < reports[0].binCount
    assert np.ptp(lastHalf.to_numpy()) < 1e-6


def test_transfer_rows_are_refitted_exactly_by_the_global_polynomial(
    faintUScene: tuple, pixelsAt: object, log: object
) -> None:
    # ARRANGE: A DEGREE-3 FIT TO THE U ROWS ALONE, WITH THE SAME CLIPPING SETTINGS AS THE SHIPPED STARE BLOCK
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
