"""The anchor measurement refuses degenerate input and the caller leaves the faint order as measured (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_anchor import bin_geometry, bright_edge_anchor, profile_edge_slit
from soxspipe.commonutils.faint_order_map import make_pixels_at
from soxspipe.commonutils.faint_order_profile import fit_fixed_width_gaussian
from soxspipe.commonutils.faint_order_trace import faint_trace_settings, transfer_faint_orders
from tests.unit._continuum_vis_scene import (
    VIS_FRAME_SHAPE,
    VIS_ORDER_ROWS,
    VIS_SEEING_SIGMA_PX,
    vis_dispersion_map,
    vis_slit_scale,
    vis_target_slit,
    vis_transfer_image,
    vis_transfer_table,
    vis_wavelength,
)

pytestmark = pytest.mark.unit

SETTINGS = faint_trace_settings.from_block({"faint-order-trace": True})
BRIGHT = 100.0
FAINT_PEAK = 0.5
FLANK_SEPARATION_ARCSEC = 2.4
FLANK_PEAK = 30.0
TARGET_PEAK = 20.0
ALL_ORDERS = (1.0, 2.0, 3.0, 4.0)


class _SilentLog:
    def debug(self, *args: object, **kwargs: object) -> None:
        return None


@pytest.fixture(scope="module", name="pixelsAt")
def vis_pixels_at(tmp_path_factory: pytest.TempPathFactory) -> object:
    mapPath = vis_dispersion_map(tmp_path_factory.mktemp("vis_map"))
    return make_pixels_at(_SilentLog(), mapPath, {"x": 1.0, "y": 1.0}, "x", "y")


def _stub_pixels_at(pixelsPerArcsec: float):
    """Return a `pixelsAt` for a map with order `n` at 100 n px and a slit of `pixelsPerArcsec` px per arcsec."""

    def pixels_at(order: float, wavelengths: np.ndarray, slitArcsec: float | np.ndarray) -> tuple:
        wavelengths = np.asarray(wavelengths, dtype=float)
        slit = np.broadcast_to(np.asarray(slitArcsec, dtype=float), wavelengths.shape)
        return 100.0 * order + pixelsPerArcsec * slit, wavelengths - 390.0

    return pixels_at


def _anchor_rows(order: float = 3.0, count: int = 20, sigmaPx: float = 1.5) -> pd.DataFrame:
    """Return `count` clean, accepted rows of one anchor order, 2 px from the map's slit-zero track."""
    wavelength = np.linspace(400.0, 500.0, count)
    return pd.DataFrame(
        {
            "order": order,
            "wavelength": wavelength,
            "cont_x": 100.0 * order + 2.0 + 0.001 * (wavelength - 450.0),
            "gauss_stddev": sigmaPx,
            "pre-clipped": False,
        }
    )


def test_bright_edge_anchor_raises_a_clear_error_for_a_table_with_no_rows() -> None:
    # ARRANGE
    rows = _anchor_rows().iloc[0:0]

    # ACT / ASSERT: THE ORDER CANNOT BE READ FROM AN EMPTY TABLE, BUT THE ERROR IS THE DOCUMENTED ONE
    with pytest.raises(ValueError, match="at least 4 accepted samples, got 0"):
        bright_edge_anchor(rows, 400.0, _stub_pixels_at(-4.0), axisA="x")


@pytest.mark.parametrize("pixelsPerArcsec", [0.0, np.nan, np.inf], ids=["zero", "nan", "inf"])
def test_bright_edge_anchor_raises_an_error_naming_the_order_and_edge_for_an_unusable_slit_scale(
    pixelsPerArcsec: float,
) -> None:
    # ARRANGE
    rows = _anchor_rows(order=3.0)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=r"order 3\.0 has no usable slit scale at 400\.0 nm"):
        bright_edge_anchor(rows, 400.0, _stub_pixels_at(pixelsPerArcsec), axisA="x")


@pytest.mark.parametrize("sigmaPx", [np.nan, 0.0, -1.5], ids=["nan", "zero", "negative"])
def test_bright_edge_anchor_raises_an_error_naming_the_order_and_edge_for_an_unusable_width(sigmaPx: float) -> None:
    # ARRANGE
    rows = _anchor_rows(order=3.0, sigmaPx=sigmaPx)

    # ACT / ASSERT
    with pytest.raises(ValueError, match=r"order 3\.0 has no usable width at 500\.0 nm"):
        bright_edge_anchor(rows, 500.0, _stub_pixels_at(-4.0), axisA="x")


def test_a_usable_anchor_is_unchanged_by_the_guards() -> None:
    # ARRANGE
    rows = _anchor_rows()

    # ACT
    slitArcsec, sigmaArcsec = bright_edge_anchor(rows, 400.0, _stub_pixels_at(-4.0), axisA="x")

    # ASSERT: 2 PX FROM THE ZERO TRACK AT 450 NM, DRIFTING 0.001 PX/NM; WIDTH 1.5 PX AT 4 PX PER ARCSEC
    assert slitArcsec == pytest.approx(-(2.0 - 0.05) / 4.0, abs=1e-6)
    assert sigmaArcsec == pytest.approx(1.5 / 4.0)


def test_transfer_leaves_a_faint_order_as_measured_when_its_anchor_has_no_usable_width(pixelsAt: object) -> None:
    # ARRANGE: THE G ORDER IS CLEAN BUT ITS GAUSSIAN WIDTHS ARE ALL NaN, SO ITS ANCHOR GEOMETRY CANNOT BE BUILT
    table = vis_transfer_table((2.0,))
    table.loc[table["order"] == 3.0, "gauss_stddev"] = np.nan
    image = vis_transfer_image({1.0: BRIGHT, 2.0: FAINT_PEAK, 3.0: BRIGHT, 4.0: BRIGHT}, dict.fromkeys(ALL_ORDERS, 0.0))
    before = table.copy()

    # ACT
    result, reports = transfer_faint_orders(
        table, image, pixelsAt, slitOffsets=None, settings=SETTINGS, axisA="x", axisB="y"
    )

    # ASSERT: THE ORDER IS REPORTED AS NOT TRANSFERRED AND THE TABLE COMES BACK AS THE SAME OBJECT, UNTOUCHED
    assert result is table
    pd.testing.assert_frame_equal(table, before)
    assert [(report.order, report.anchorOrder, report.is_transferred) for report in reports] == [(2.0, 3.0, False)]
    assert "anchor unusable" in reports[0].describe("2 (u)", "3 (g)")


def test_transfer_still_transfers_the_other_band_group_when_one_anchor_is_unusable(pixelsAt: object) -> None:
    # ARRANGE: BOTH THE U AND I ORDERS ARE FAINT, AND THE G ANCHOR IS UNUSABLE BUT THE R ANCHOR IS FINE
    table = vis_transfer_table((2.0, 4.0))
    table.loc[table["order"] == 3.0, "gauss_stddev"] = np.nan
    image = vis_transfer_image(
        {1.0: BRIGHT, 2.0: FAINT_PEAK, 3.0: BRIGHT, 4.0: FAINT_PEAK},
        dict.fromkeys(ALL_ORDERS, 0.0) | {2.0: 0.8, 4.0: 0.8},
    )

    # ACT
    result, reports = transfer_faint_orders(
        table, image, pixelsAt, slitOffsets=None, settings=SETTINGS, axisA="x", axisB="y"
    )

    # ASSERT
    assert [(report.order, report.is_transferred) for report in reports] == [(2.0, False), (4.0, True)]
    assert result["faint_transfer"].tolist() == (result["order"] == 4.0).tolist()


def _edge_geometry(order: float, edgeNm: float, curveSlit: float) -> bin_geometry:
    return bin_geometry(order, curveSlit, 0.4, 1.5, edgeNm, 3.0)


def _clean_image() -> np.ndarray:
    return vis_transfer_image(dict.fromkeys(ALL_ORDERS, 20.0), dict.fromkeys(ALL_ORDERS, 0.0))


@pytest.mark.parametrize(
    "wavelengths",
    [np.array([]), np.array([np.nan, np.nan, np.nan])],
    ids=["empty", "all_nan"],
)
def test_profile_edge_slit_returns_none_when_the_anchor_has_no_wavelengths(
    wavelengths: np.ndarray, pixelsAt: object
) -> None:
    # ARRANGE / ACT
    slit = profile_edge_slit(
        _clean_image(),
        pixelsAt,
        _edge_geometry(3.0, 600.0, float(vis_target_slit(600.0))),
        anchorWavelengths=wavelengths,
        widthNm=10.0,
        axisA="x",
    )

    # ASSERT: THE CALLER FALLS BACK TO THE CURVE ESTIMATE
    assert slit is None


def test_fit_fixed_width_gaussian_returns_none_for_an_empty_profile() -> None:
    # ARRANGE / ACT / ASSERT
    assert fit_fixed_width_gaussian(np.array([]), np.array([]), 1.5, 3.0, 3.0) is None


def _image_with_flank(flankSeparationArcsec: float) -> np.ndarray:
    """Return the clean scene with a source of `FLANK_PEAK` `flankSeparationArcsec` above the target in order 3."""
    image = vis_transfer_image(dict.fromkeys(ALL_ORDERS, TARGET_PEAK), dict.fromkeys(ALL_ORDERS, 0.0))
    rows = np.arange(10, 10 + VIS_ORDER_ROWS)
    columns = np.arange(VIS_FRAME_SHAPE[1], dtype=float)
    slit = vis_target_slit(vis_wavelength(3.0, rows)) + flankSeparationArcsec
    centre = (-30.0 + 90.0 * 3.0 + vis_slit_scale(3.0) * slit)[:, np.newaxis]
    image[rows] += FLANK_PEAK * np.exp(-0.5 * ((columns[np.newaxis, :] - centre) / VIS_SEEING_SIGMA_PX) ** 2)
    return image


def test_profile_edge_slit_does_not_land_on_a_2_4_arcsec_flank_when_the_curve_estimate_is_biased_towards_it(
    pixelsAt: object,
) -> None:
    # ARRANGE: A SOURCE BRIGHTER THAN THE TARGET SITS 2.4 ARCSEC AWAY, AND THE CURVE ESTIMATE IS 1 ARCSEC TOWARDS IT
    truth = float(vis_target_slit(600.0))
    geometry = _edge_geometry(3.0, 600.0, truth + 1.0)

    # ACT
    slit = profile_edge_slit(
        _image_with_flank(FLANK_SEPARATION_ARCSEC),
        pixelsAt,
        geometry,
        anchorWavelengths=np.linspace(600.0, 800.0, 150),
        widthNm=10.0,
        axisA="x",
    )

    # ASSERT: THE FIT EITHER FAILS, AND THE CALLER FALLS BACK, OR IT FINDS THE TARGET RATHER THAN THE FLANK
    assert slit is None or abs(slit - truth) < 0.3


def test_profile_edge_slit_tracks_the_edge_position_when_the_trace_slopes_in_wavelength(pixelsAt: object) -> None:
    # ARRANGE: THE TARGET DRIFTS 0.01 ARCSEC PER NM, TWENTY TIMES THE SCENE'S DRIFT, ACROSS A 18.6 NM STRETCH
    slopeArcsecPerNm = 0.01
    rows = np.arange(10, 10 + VIS_ORDER_ROWS)
    columns = np.arange(VIS_FRAME_SHAPE[1], dtype=float)
    image = vis_transfer_image(dict.fromkeys(ALL_ORDERS, 0.0), dict.fromkeys(ALL_ORDERS, 0.0))
    truth = float(vis_target_slit(600.0))
    slit = truth + slopeArcsecPerNm * (vis_wavelength(3.0, rows) - 600.0)
    centre = (-30.0 + 90.0 * 3.0 + vis_slit_scale(3.0) * slit)[:, np.newaxis]
    image[rows] += TARGET_PEAK * np.exp(-0.5 * ((columns[np.newaxis, :] - centre) / VIS_SEEING_SIGMA_PX) ** 2)

    # ACT
    measured = profile_edge_slit(
        image,
        pixelsAt,
        _edge_geometry(3.0, 600.0, truth),
        anchorWavelengths=np.linspace(600.0, 800.0, 150),
        widthNm=18.6,
        axisA="x",
    )

    # ASSERT: THE STRETCH IS CENTRED 9.3 NM FROM THE EDGE, SO THE TRACE HAS DRIFTED BY 0.093 ARCSEC AT ITS MIDDLE
    assert measured is not None
    assert measured == pytest.approx(truth, abs=0.15)
