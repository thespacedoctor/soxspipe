"""A faint anchor is measured at the shared edge from one collapsed profile, not from its noisy curve (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_anchor import bin_geometry, bright_edge_anchor, profile_edge_slit
from soxspipe.commonutils.faint_order_map import make_pixels_at
from soxspipe.commonutils.faint_order_trace import faint_trace_settings, transfer_faint_orders, transfer_report
from tests.unit._continuum_vis_scene import (
    vis_dispersion_map,
    vis_target_slit,
    vis_transfer_image,
    vis_transfer_table,
)

pytestmark = pytest.mark.unit

SETTINGS = faint_trace_settings.from_block({"faint-order-trace": True})
BRIGHT = 100.0
ANCHOR_PEAK = 3.0
FLANK_PEAK = 30.0
FAINT_PEAK = 0.5
# THE ANCHOR'S SAMPLES NEAREST THE SHARED EDGE ARE PUSHED THIS FAR (PIXELS) AND THE REST SCATTER BY 2.5 PX
EDGE_BIAS_PX = 4.0
EDGE_BIAS_SAMPLES = 25
SCATTER_PX = 2.5
# (FAINT ORDER, ANCHOR ORDER, SHARED EDGE IN NM, WHICH END OF THE ANCHOR ORDER FACES THE FAINT ORDER)
PAIRS = {"u_from_g": (2.0, 3.0, 600.0, "blue"), "i_from_r": (4.0, 1.0, 1000.0, "red")}


class _SilentLog:
    def debug(self, *args: object, **kwargs: object) -> None:
        return None


@pytest.fixture(scope="module", name="pixelsAt")
def vis_pixels_at(tmp_path_factory: pytest.TempPathFactory) -> object:
    mapPath = vis_dispersion_map(tmp_path_factory.mktemp("vis_map"))
    return make_pixels_at(_SilentLog(), mapPath, {"x": 1.0, "y": 1.0}, "x", "y")


def _scene(pair: str, *, anchorPeak: float = ANCHOR_PEAK) -> tuple[pd.DataFrame, np.ndarray]:
    """Return the table and frame in which `pair`'s anchor order is faint, biased at the shared edge, and faint-sourced.

    The anchor's sampled centres scatter by 2.5 px, and the ones nearest the shared edge are also pushed 4 px, which
    drags a degree-2 curve through them about an arcsecond from the truth at that edge. The anchor order also holds
    two flanking sources as bright as the target, 3.7 arcsec either side, as VV4580Sgr does.
    """
    faint, anchor, _, edge = PAIRS[pair]
    table = vis_transfer_table((faint,))
    isAnchor = (table["order"] == anchor).to_numpy()
    noise = np.random.default_rng(4580).normal(0.0, SCATTER_PX, int(isAnchor.sum()))
    biased = slice(0, EDGE_BIAS_SAMPLES) if edge == "blue" else slice(-EDGE_BIAS_SAMPLES, None)
    noise[biased] += EDGE_BIAS_PX
    table.loc[isAnchor, "cont_x"] = table.loc[isAnchor, "cont_x"] + noise
    table.loc[isAnchor, "centre_shift"] = table.loc[isAnchor, "centre_shift"] - noise
    amplitudes = {1.0: BRIGHT, 2.0: BRIGHT, 3.0: BRIGHT, 4.0: BRIGHT} | {faint: FAINT_PEAK, anchor: anchorPeak}
    flanks = dict.fromkeys(amplitudes, 0.0) | {anchor: FLANK_PEAK}
    return table, vis_transfer_image(amplitudes, flanks)


def _transfer(table: pd.DataFrame, image: np.ndarray, pixelsAt: object, offsets: dict | None = None) -> tuple:
    return transfer_faint_orders(table, image, pixelsAt, slitOffsets=offsets, settings=SETTINGS, axisA="x", axisB="y")


def _curve_estimate(table: pd.DataFrame, pair: str, pixelsAt: object) -> float:
    _, anchor, edgeNm, _ = PAIRS[pair]
    return bright_edge_anchor(table.loc[table["order"] == anchor], edgeNm, pixelsAt, axisA="x")[0]


@pytest.mark.parametrize("pair", list(PAIRS))
def test_a_faint_anchor_is_measured_from_the_edge_profile_when_its_curve_is_biased(pair: str, pixelsAt: object) -> None:
    # ARRANGE: THE CURVE THROUGH THE ANCHOR'S NOISY SAMPLES LANDS ABOUT AN ARCSECOND FROM THE SOURCE AT THE EDGE
    faint, anchor, edgeNm, _ = PAIRS[pair]
    table, image = _scene(pair)
    truth = float(vis_target_slit(edgeNm))
    assert abs(_curve_estimate(table, pair, pixelsAt) - truth) > 0.5

    # ACT
    _, reports = _transfer(table, image, pixelsAt)

    # ASSERT: THE COLLAPSED PROFILE FINDS THE SOURCE
    assert [(report.order, report.anchorOrder) for report in reports] == [(faint, anchor)]
    assert reports[0].anchorSource == "profile"
    assert reports[0].slitArcsec == pytest.approx(truth, abs=0.15)


def test_a_faint_anchor_falls_back_to_the_curve_estimate_when_the_profile_fit_fails(pixelsAt: object) -> None:
    # ARRANGE: THE ANCHOR'S ORDER HOLDS NOTHING BUT NOISE IN THE FRAME, SO NO PROFILE PEAK CAN BE FITTED
    table, image = _scene("u_from_g", anchorPeak=0.0)

    # ACT
    _, reports = _transfer(table, image, pixelsAt)

    # ASSERT: THE CURVE ESTIMATE STANDS, AND THE REPORT SAYS THE PROFILE FAILED
    assert reports[0].anchorSource == "curve after profile failed"
    assert reports[0].slitArcsec == pytest.approx(_curve_estimate(table, "u_from_g", pixelsAt))


def test_an_anchor_that_is_not_faint_keeps_the_curve_estimate(pixelsAt: object) -> None:
    # ARRANGE: G IS CLEAN, SO IT IS NOT FAINT. THE FRAME PUTS ITS SOURCE ARCSECONDS AWAY FROM THE TABLE'S TRACE, WHICH
    # A PROFILE FIT WOULD FOLLOW AND THE CURVE ESTIMATE DOES NOT
    table = vis_transfer_table((2.0,))
    image = vis_transfer_image(
        {1.0: BRIGHT, 2.0: FAINT_PEAK, 3.0: BRIGHT, 4.0: BRIGHT}, dict.fromkeys(range(1, 5), 0.0)
    )
    table = table.assign(cont_x=table["cont_x"].where(table["order"] != 3.0, table["cont_x"] + 2.0))

    # ACT
    _, reports = _transfer(table, image, pixelsAt)

    # ASSERT
    assert reports[0].anchorSource == "curve"
    assert reports[0].slitArcsec == pytest.approx(_curve_estimate(table, "u_from_g", pixelsAt))


def test_the_slit_offsets_move_a_profile_anchor_by_their_difference(pixelsAt: object) -> None:
    # ARRANGE
    table, image = _scene("u_from_g")
    _, plain = _transfer(table, image, pixelsAt)

    # ACT
    _, shifted = _transfer(table, image, pixelsAt, {2.0: 0.5, 3.0: 0.0})

    # ASSERT: ANCHOR = S_BRIGHT - DELTA_BRIGHT + DELTA_FAINT, WITH S_BRIGHT FROM THE PROFILE
    assert shifted[0].anchorSource == "profile"
    assert shifted[0].slitArcsec - plain[0].slitArcsec == pytest.approx(0.5, abs=1e-6)


@pytest.mark.parametrize(
    ("source", "suffix"),
    [
        ("curve", ""),
        ("profile", ", anchor from the edge profile"),
        ("curve after profile failed", ", edge profile fit failed so the anchor is the curve estimate"),
    ],
)
def test_the_report_names_where_the_anchor_came_from_only_when_it_is_not_the_curve(source: str, suffix: str) -> None:
    # ARRANGE
    transferred = transfer_report(2.0, 3.0, 9, 15, 1, 0.3, source)
    skipped = transfer_report(2.0, 3.0, 0, 15, 0, 0.3, source)

    # ACT / ASSERT
    assert (
        transferred.describe("2 (u)", "3 (g)") == f"Order 2 (u) trace transferred from 3 (g): 9/15 bins, deg 1{suffix}"
    )
    assert skipped.describe("2 (u)", "3 (g)") == (
        f"Order 2 (u) from anchor 3 (g): no bin kept, trace left as measured{suffix}"
    )


def _edge_geometry(order: float, edgeNm: float, curveSlit: float) -> bin_geometry:
    return bin_geometry(order, curveSlit, 0.4, 1.5, edgeNm, 3.0)


@pytest.mark.parametrize(
    ("order", "edgeNm", "start"), [(3.0, 600.0, 600.0), (1.0, 1000.0, 800.0)], ids=["blue_edge", "red_edge"]
)
def test_profile_edge_slit_reads_the_stretch_next_to_the_shared_edge(
    order: float, edgeNm: float, start: float, pixelsAt: object
) -> None:
    # ARRANGE: THE TARGET DRIFTS 0.5 ARCSEC ACROSS THE SCENE. THE CURVE ESTIMATE IS 0.4 ARCSEC OFF
    image = vis_transfer_image(dict.fromkeys((1.0, 2.0, 3.0, 4.0), 20.0), dict.fromkeys((1.0, 2.0, 3.0, 4.0), 0.0))
    wavelength = np.linspace(start, start + 200.0, 150)
    truth = float(vis_target_slit(edgeNm))

    # ACT
    slit = profile_edge_slit(
        image,
        pixelsAt,
        _edge_geometry(order, edgeNm, truth + 0.4),
        anchorWavelengths=wavelength,
        widthNm=10.0,
        axisA="x",
    )

    # ASSERT: THE 10 NM STRETCH NEXT TO THE EDGE HOLDS THE SOURCE WITHIN 0.05 ARCSEC OF ITS EDGE POSITION
    assert slit == pytest.approx(truth, abs=0.05)


def test_profile_edge_slit_stops_at_the_end_of_the_anchor_order_when_the_stretch_is_wider_than_it(
    pixelsAt: object,
) -> None:
    # ARRANGE
    image = vis_transfer_image(dict.fromkeys((1.0, 2.0, 3.0, 4.0), 20.0), dict.fromkeys((1.0, 2.0, 3.0, 4.0), 0.0))
    wavelength = np.linspace(600.0, 800.0, 150)

    # ACT
    slit = profile_edge_slit(
        image,
        pixelsAt,
        _edge_geometry(3.0, 600.0, float(vis_target_slit(600.0))),
        anchorWavelengths=wavelength,
        widthNm=5000.0,
        axisA="x",
    )

    # ASSERT: THE WHOLE ORDER IS COLLAPSED, SO THE TRACE'S 0.5 ARCSEC DRIFT BLURS IT BUT THE SOURCE IS STILL FOUND
    assert slit == pytest.approx(float(vis_target_slit(700.0)), abs=0.2)


def test_profile_edge_slit_returns_none_for_a_blank_stretch(pixelsAt: object) -> None:
    # ARRANGE
    image = vis_transfer_image(dict.fromkeys((1.0, 2.0, 3.0, 4.0), 0.0), dict.fromkeys((1.0, 2.0, 3.0, 4.0), 0.0))

    # ACT
    slit = profile_edge_slit(
        image,
        pixelsAt,
        _edge_geometry(3.0, 600.0, 0.5),
        anchorWavelengths=np.linspace(600.0, 800.0, 150),
        widthNm=20.0,
        axisA="x",
    )

    # ASSERT
    assert slit is None
