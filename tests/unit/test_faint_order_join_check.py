"""The same-source check: adjacent VIS orders must put the trace at the same slit position at their join (DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.faint_order_map import make_pixels_at
from soxspipe.commonutils.faint_order_trace import join_slit_mismatches
from tests.unit._continuum_vis_scene import VIS_START_NM, vis_dispersion_map, vis_slit_scale

pytestmark = pytest.mark.unit

LIMIT_ARCSEC = 2.0


class _SilentLog:
    def debug(self, *args: object, **kwargs: object) -> None:
        return None


@pytest.fixture(scope="module", name="pixelsAt")
def vis_pixels_at(tmp_path_factory: pytest.TempPathFactory) -> object:
    mapPath = vis_dispersion_map(tmp_path_factory.mktemp("vis_join_map"))
    return make_pixels_at(_SilentLog(), mapPath, {"x": 1.0, "y": 1.0}, "x", "y")


def _table() -> pd.DataFrame:
    frames = [
        pd.DataFrame({"order": order, "wavelength": np.linspace(start, start + 200.0, 30)})
        for order, start in sorted(VIS_START_NM.items())
    ]
    return pd.concat(frames, ignore_index=True)


def _trace_at_slit(slitByOrder: dict[float, float]):
    """Return a fitted-trace function whose trace sits at `slitByOrder[order]` arcsec in each order."""

    def trace(order: float, axisB: np.ndarray) -> np.ndarray:
        return np.full(np.shape(axisB), -30.0 + 90.0 * order + vis_slit_scale(order) * slitByOrder[order])

    return trace


def _same_source(**overrides: float) -> dict[float, float]:
    return {1.0: 0.4, 2.0: 0.4, 3.0: 0.4, 4.0: 0.4} | {float(key[1:]): value for key, value in overrides.items()}


def _mismatches(pixelsAt: object, slitByOrder: dict, offsets: dict | None = None, limit: float = LIMIT_ARCSEC) -> list:
    return join_slit_mismatches(_table(), _trace_at_slit(slitByOrder), pixelsAt, offsets, limit)


def test_join_slit_mismatches_is_silent_when_every_order_traces_the_same_source(pixelsAt: object) -> None:
    assert _mismatches(pixelsAt, _same_source()) == []


def test_join_slit_mismatches_reports_a_trace_that_jumps_by_more_than_the_limit(pixelsAt: object) -> None:
    # ARRANGE: THE U TRACE SITS 2.5 ARCSEC ABOVE THE G TRACE AT THEIR JOIN
    slits = _same_source(o2=0.4 + 2.5)

    # ACT
    found = _mismatches(pixelsAt, slits)

    # ASSERT: THE DIFFERENCE IS RED MINUS BLUE, SO G IS 2.5 ARCSEC BELOW U
    assert [(m.blueOrder, m.redOrder) for m in found] == [(2.0, 3.0)]
    assert found[0].differenceArcsec == pytest.approx(-2.5, abs=0.01)


def test_join_slit_mismatches_keeps_the_sign_of_the_jump(pixelsAt: object) -> None:
    # ARRANGE
    slits = _same_source(o2=0.4 - 2.5)

    # ACT
    found = _mismatches(pixelsAt, slits)

    # ASSERT
    assert found[0].differenceArcsec == pytest.approx(2.5, abs=0.01)


def test_join_slit_mismatches_is_silent_at_or_below_the_limit(pixelsAt: object) -> None:
    assert _mismatches(pixelsAt, _same_source(o2=0.4 + 1.9)) == []
    assert _mismatches(pixelsAt, _same_source(o2=0.4 + 2.0)) == []


def test_join_slit_mismatches_checks_the_join_across_the_band_groups(pixelsAt: object) -> None:
    # ARRANGE: G AND R ARE THE ORDERS EITHER SIDE OF THE U/G TO R/I BOUNDARY
    slits = _same_source(o1=0.4 + 3.0)

    # ACT
    found = _mismatches(pixelsAt, slits)

    # ASSERT: R IS 3 ARCSEC ABOVE G, AND I, WHICH MATCHES G, IS 3 ARCSEC BELOW R
    assert [(m.blueOrder, m.redOrder) for m in found] == [(3.0, 1.0), (1.0, 4.0)]
    assert found[0].differenceArcsec == pytest.approx(3.0, abs=0.01)
    assert found[1].differenceArcsec == pytest.approx(-3.0, abs=0.01)


def test_join_slit_mismatches_removes_the_slit_offsets_both_orders_report(pixelsAt: object) -> None:
    # ARRANGE: THE U TRACE LOOKS 2.5 ARCSEC OFF, BUT THE SKY SAYS THE WHOLE U IMAGE IS 2.5 ARCSEC OFF
    slits = _same_source(o2=0.4 + 2.5)

    # ACT
    found = _mismatches(pixelsAt, slits, offsets={2: 2.5, 3: 0.0})

    # ASSERT
    assert found == []


def test_join_slit_mismatches_ignores_the_offset_of_a_single_order(pixelsAt: object) -> None:
    # ARRANGE
    slits = _same_source(o2=0.4 + 2.5)

    # ACT
    found = _mismatches(pixelsAt, slits, offsets={2: 2.5})

    # ASSERT: WITH NO OFFSET FOR G THE RAW DIFFERENCE STANDS
    assert len(found) == 1


def test_join_slit_mismatches_skips_a_pair_with_an_order_missing_from_the_table(pixelsAt: object) -> None:
    # ARRANGE: NO G ROWS, SO THE U/G AND G/R JOINS CANNOT BE CHECKED
    table = _table()
    table = table.loc[table["order"] != 3.0]
    slits = _same_source(o2=0.4 + 2.5)

    # ACT
    found = join_slit_mismatches(table, _trace_at_slit(slits), pixelsAt, None, LIMIT_ARCSEC)

    # ASSERT
    assert found == []
