"""`soxs-offset` classifies its frames per template run, for both offset templates, before reducing (DY-2354).

These tests run `produce_product` on frames written to disk, with the shared collaborators of the offset
characterization tests replaced at their boundary. They check which frames reach the ON and OFF stacks and extractions.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest

from soxspipe.commonutils.offset_frames import OffsetClassificationError
from soxspipe.recipes import soxs_offset
from tests.factories.offset_frames import (
    FIXED_SKY_TECHNIQUE,
    GENERIC_TECHNIQUE,
    fixed_sky_run,
    generic_offset_run,
)
from tests.integration.test_observing_recipe_orchestration import (
    RouteCollection,
    _configure_nodding_recipe,
    _route,
    _write_prepared_frame,
)
from tests.integration.test_soxs_offset_reduction_characterization import (
    _calibration_routes,
    _extract_each_cycle,
    _patch,
)

pytestmark = pytest.mark.integration


def _recipe(tmpPath: Path, log: Any, frames: list[tuple[str, str, dict[str, object]]]) -> tuple[soxs_offset, list[str]]:
    """Return an unconstructed offset recipe whose SOF holds ``frames``, each a (DPR TYPE, DPR TECH, header) triple.

    **Return:**

    - the recipe, and the file names of the frames in the order given
    """
    routes: dict[frozenset[tuple[str, object]], list[str]] = {}
    names = []
    for index, (dprType, technique, header) in enumerate(frames):
        path = _write_prepared_frame(
            tmpPath / f"frame-{index}_pre.fits",
            seed=700 + index,
            headerOverrides={**header, "DPR_TYPE": dprType, "DPR_TECH": technique, "EXPTIME": 300.0},
        )
        routes.setdefault(_route(DPR_TYPE=dprType, DPR_TECH=technique), []).append(str(path))
        names.append(f"frame-{index}.fits")
    recipe = soxs_offset.__new__(soxs_offset)
    _configure_nodding_recipe(recipe, log=log, tmp_path=tmpPath, objectPaths=[], technique=FIXED_SKY_TECHNIQUE)
    recipe.inputFrames = RouteCollection({**routes, **_calibration_routes(tmpPath)})
    return recipe, names


def _exposures(frames: list[Any]) -> list[int]:
    """Return the template exposure number of each frame."""
    return [frame.header["ESO TPL EXPNO"] for frame in frames]


def test_a_fixed_sky_run_with_sky_frames_labelled_sky_stacks_the_object_and_the_sky_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    headers = fixed_sky_run([(0.0, -0.4), (-2.0, 20.4)])
    recipe, _ = _recipe(
        tmp_path,
        log,
        [("OBJECT", FIXED_SKY_TECHNIQUE, headers[0]), ("SKY", FIXED_SKY_TECHNIQUE, headers[1])],
    )
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    onStack, offStack = captured["stack"]
    assert _exposures(onStack["frames"]) == [1]
    assert _exposures(offStack["frames"]) == [2]


def test_a_generic_offset_sof_stacks_its_object_frames_as_on_and_its_sky_frames_as_off(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    headers = generic_offset_run(
        [(5.4, -7.9), (3.4, 12.1), (5.5, -7.9), (3.4, 12.0)], ["OBJECT", "SKY", "OBJECT", "SKY"]
    )
    recipe, _ = _recipe(tmp_path, log, [(header["ESO DPR TYPE"], GENERIC_TECHNIQUE, header) for header in headers])
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    onStack, offStack = captured["stack"]
    assert _exposures(onStack["frames"]) == [1, 3]
    assert _exposures(offStack["frames"]) == [2, 4]


def test_an_unclassifiable_frame_fails_the_sof_before_any_stacking_or_extraction(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    headers = generic_offset_run(
        [(2.8, -16.3), (4.8, -14.3), (8.8, -10.3), (14.8, -4.2)], ["OBJECT", "SKY", "OBJECT", "SKY"]
    )
    recipe, names = _recipe(tmp_path, log, [(header["ESO DPR TYPE"], GENERIC_TECHNIQUE, header) for header in headers])
    calls, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT / ASSERT
    with pytest.raises(OffsetClassificationError, match=names[2]):
        recipe.produce_product()

    assert calls == []
    assert captured["stack"] == []
    assert captured["extract_cycle"] == []


def test_the_cycles_of_a_sof_with_two_runs_pair_an_on_and_an_off_frame_of_the_same_run(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE SECOND RUN'S FRAMES ARE TAKEN BETWEEN THE FIRST RUN'S, SO A SORT BY TIME ALONE WOULD CROSS THE RUNS
    runA = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T04:00:00")
    runB = fixed_sky_run([(0.0, 0.0), (20.0, 20.0)], start="2025-05-14T05:00:00", acquisition=(50.0, 40.0))
    mjds = [60000.0, 60000.4, 60000.1, 60000.2]
    headers = [{**header, "MJD-OBS": mjd} for header, mjd in zip(runA + runB, mjds, strict=True)]
    recipe, _ = _recipe(tmp_path, log, [("OBJECT", FIXED_SKY_TECHNIQUE, header) for header in headers])
    _extract_each_cycle(recipe)
    _, captured = _patch(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    pairedRuns = [
        (cycle["aFrame"].header["ESO TPL START"], cycle["bFrame"].header["ESO TPL START"])
        for cycle in captured["extract_cycle"]
    ]
    assert pairedRuns == [("2025-05-14T04:00:00",) * 2, ("2025-05-14T05:00:00",) * 2]
