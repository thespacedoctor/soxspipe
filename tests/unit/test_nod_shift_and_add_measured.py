"""A nod pair whose positions are measured by the slit alignment lands on the A position when combined (DY-1674)."""

from __future__ import annotations

import numpy as np
import pytest

from soxspipe.commonutils.nod_shift_and_add import build_shift_and_add_plan, combine_nod_differences
from soxspipe.commonutils.slit_alignment import align_frame_groups_along_slit, frame_group
from tests.factories import pipeline_settings
from tests.factories.slit_alignment import synthetic_alignment_map, synthetic_nod_frame
from tests.unit.test_slit_alignment_frames import DETECT_CONTINUUM_SETTINGS, _use_synthetic_trace_samples

pytestmark = pytest.mark.unit

MIDDLE_ROWS = slice(200, 400)


@pytest.mark.parametrize(("positionA", "positionB"), [(2.0, -2.0), (-1.2, 2.4)])
def test_measured_positions_put_the_positive_trace_of_the_combination_at_the_a_position(
    log, tmp_path, monkeypatch, positionA: float, positionB: float
) -> None:
    # ARRANGE: ONE A FRAME AND ONE B FRAME, WHICH THE ALIGNMENT MEASURES BUT DOES NOT RESAMPLE
    scene = synthetic_alignment_map(tmp_path)
    _use_synthetic_trace_samples(monkeypatch, scene)
    frameA = synthetic_nod_frame(scene, objectSlit=positionA, mjd=60000.0, hasSky=False)
    frameB = synthetic_nod_frame(scene, objectSlit=positionB, mjd=60000.005, offset=-3.0, hasSky=False)
    resultA, resultB = align_frame_groups_along_slit(
        log=log,
        groups=[frame_group([frameA], [frameB], ["A0.fits"]), frame_group([frameB], [frameA], ["B0.fits"])],
        twoDMapPath=scene.path,
        dispersionMap="unused.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"detect-continuum": DETECT_CONTINUUM_SETTINGS},
        recipeName="soxs-nod",
        perOrder=True,
        startNightDate="2024-01-02",
        measureSingleFrames=True,
    )
    abDifference = frameA.subtract(frameB)
    abDifference.header = frameA.header
    baDifference = frameB.subtract(frameA)
    baDifference.header = frameB.header

    # ACT
    plan = build_shift_and_add_plan(log, resultA, resultB)
    combined = combine_nod_differences(abDifference, baDifference, plan.orderShifts, plan.mapImages)

    # ASSERT: THE POSITIVE PEAK OF EVERY ORDER IS AT THE A POSITION, WITH THE HEIGHT OF ONE DIFFERENCE
    assert plan.orderShifts == pytest.approx(dict.fromkeys(scene.orders, positionB - positionA), abs=0.02)
    for order in scene.orders:
        isOrder = (scene.orderMap == order)[MIDDLE_ROWS]
        slit = np.where(isOrder, scene.slit[MIDDLE_ROWS], np.nan)
        counts = np.where(isOrder, combined.data[MIDDLE_ROWS], np.nan)
        columns = np.flatnonzero(isOrder.all(axis=0))
        peakColumn = columns[np.argmax(np.median(counts[:, columns], axis=0))]
        assert np.median(slit[:, peakColumn]) == pytest.approx(positionA, abs=0.15)
        oneDifference = np.nanmedian(np.where(isOrder, abDifference.data[MIDDLE_ROWS], np.nan)[:, columns].max(axis=1))
        assert np.nanmedian(counts[:, peakColumn]) == pytest.approx(oneDifference, rel=0.1)
