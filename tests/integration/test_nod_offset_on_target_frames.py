"""The nod and offset recipes count the frames that reach their merged spectrum as on target (DY-2356).

A nod A or B frame and an offset ON frame are on target. An offset OFF frame is not. A nod position, or an offset ON
frame, that finds no object trace is not in the spectrum, so its frame leaves the list with a warning that names it. The
other position of the cycle stays (Dave, 2026-11-24).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from tests.factories.offset_frames import GENERIC_TECHNIQUE, generic_offset_run
from tests.integration.test_observing_recipe_orchestration import _route
from tests.integration.test_soxs_nod_reduction_characterization import (
    CALIBRATION_NAMES,
    SINGLE_PAIR,
    TWO_LOCATIONS,
    _extract_each_cycle,
    _nod_recipe,
)
from tests.integration.test_soxs_nod_reduction_characterization import _patch as _patch_nod
from tests.integration.test_soxs_offset_reduction_characterization import (
    OFFSET_TECHNIQUE,
    TWO_CYCLES,
    _offset_headers,
    _offset_recipe,
)
from tests.integration.test_soxs_offset_reduction_characterization import _extract_each_cycle as _each_offset_cycle
from tests.integration.test_soxs_offset_reduction_characterization import _patch as _patch_offset

pytestmark = pytest.mark.integration

SPECTRUM = pd.DataFrame({"WAVE": [500.0], "SNR": [20.0]})


def _with_arcfiles(headers: list[dict[str, object]]) -> list[dict[str, object]]:
    """Give each header an ``ARCFILE`` that differs from its file name."""
    return [{**header, "ARCFILE": f"SOXS.{index}.fits"} for index, header in enumerate(headers)]


def _names(entries: Any) -> list[str]:
    return [entry.name for entry in entries]


def _list_at_each_stack(recipe: Any, monkeypatch: pytest.MonkeyPatch) -> list[Any]:
    """Record the on-target list at each call of `stack_extractions`, which stamps the merged spectrum from it."""
    seen: list[Any] = []
    recorded = recipe.stack_extractions

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(recipe.onTargetFrames)
        return recorded(*args, **kwargs)

    monkeypatch.setattr(recipe, "stack_extractions", spy)
    return seen


def _without_trace(recipe: Any, monkeypatch: pytest.MonkeyPatch, missing: set[tuple[Any, str]]) -> None:
    """Make the extraction of a position find no object trace: ``missing`` holds ``(cycle, position)`` pairs."""

    def extraction(**kwargs: Any) -> tuple[Any, Any, dict[str, int]]:
        cycle = kwargs["locationSetIndex"]
        spectrumA = None if (cycle, "A") in missing else SPECTRUM.copy()
        if "nod" not in recipe.recipeName:
            return spectrumA, False, {"10-11": 11}
        return spectrumA, None if (cycle, "B") in missing else SPECTRUM.copy(), {"10-11": 11}

    monkeypatch.setattr(recipe, "process_single_ab_nodding_cycle", extraction)


def _left_out_warnings(log: Any) -> list[str]:
    return [m for level, m in log.messages if level == "warning" and "left out of the on-target list" in m]


# NOD


def test_every_a_and_b_frame_of_a_stacked_nod_is_on_target(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(SINGLE_PAIR))
    _patch_nod(recipe, monkeypatch, tmp_path)
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _names(seen[0]) == ["SOXS.0.fits", "SOXS.1.fits"]
    assert [entry.exposureSeconds for entry in seen[0]] == [300.0, 300.0]


def test_every_frame_of_cycles_that_all_found_a_trace_is_on_target(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(TWO_LOCATIONS))
    _extract_each_cycle(recipe)
    _patch_nod(recipe, monkeypatch, tmp_path)
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: THE LIST IS IN TIME ORDER
    assert _names(seen[0]) == ["SOXS.0.fits", "SOXS.2.fits", "SOXS.1.fits", "SOXS.3.fits"]
    assert not _left_out_warnings(log)


def test_a_nod_position_with_no_trace_in_one_cycle_leaves_the_list_and_the_other_position_stays(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE B FRAME OF CYCLE 2 IS nod-3, WHICH IS SOXS.3
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(TWO_LOCATIONS))
    _extract_each_cycle(recipe)
    _patch_nod(recipe, monkeypatch, tmp_path)
    _without_trace(recipe, monkeypatch, {(2, "B")})
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _names(seen[0]) == ["SOXS.0.fits", "SOXS.2.fits", "SOXS.1.fits"]
    warnings = _left_out_warnings(log)
    assert len(warnings) == 1 and "SOXS.3.fits" in warnings[0] and "B" in warnings[0]


def test_a_nod_cycle_with_no_trace_at_either_position_leaves_the_list_with_a_warning_for_each(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(TWO_LOCATIONS))
    _extract_each_cycle(recipe)
    _patch_nod(recipe, monkeypatch, tmp_path)
    _without_trace(recipe, monkeypatch, {(1, "A"), (1, "B")})
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: CYCLE 1 IS SOXS.0 AND SOXS.2
    assert _names(seen[0]) == ["SOXS.1.fits", "SOXS.3.fits"]
    assert len(_left_out_warnings(log)) == 2


def test_a_stacked_nod_position_with_no_trace_leaves_the_list_and_the_other_position_stays(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(SINGLE_PAIR))
    _patch_nod(recipe, monkeypatch, tmp_path)
    _without_trace(recipe, monkeypatch, {(1, "A")})
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: THE A FRAME IS SOXS.0
    assert _names(seen[0]) == ["SOXS.1.fits"]
    warnings = _left_out_warnings(log)
    assert len(warnings) == 1 and "SOXS.0.fits" in warnings[0]


def test_the_cycle_inspection_pass_counts_its_own_cycles_and_the_recipe_gets_its_list_back(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE STACKED PASS LEFT SOXS.3 OUT OF THE LIST
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(TWO_LOCATIONS))
    _patch_nod(recipe, monkeypatch, tmp_path)
    recipe.produce_product()
    stackedList = recipe.onTargetFrames[:3]
    recipe.onTargetFrames = stackedList
    seen: list[Any] = []
    monkeypatch.setattr(recipe, "_inspect_each_cycle", lambda *_: seen.append(recipe.onTargetFrames) or pd.DataFrame())
    cycleFrames = ([], [], ["nod-0.fits", "nod-1.fits"], ["nod-2.fits", "nod-3.fits"])

    # ACT
    recipe._extract_each_cycle_as_well(cycleFrames, pd.DataFrame(), "order.fits", False)

    # ASSERT
    assert _names(seen[0]) == ["SOXS.0.fits", "SOXS.2.fits", "SOXS.1.fits", "SOXS.3.fits"]
    assert recipe.onTargetFrames == stackedList


def test_a_nir_nod_frame_spends_its_dit_times_ndit_on_target(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    headers = [{**header, "HIERARCH ESO DET SEQ1 DIT": 900.0, "HIERARCH ESO DET NDIT": 2} for header in SINGLE_PAIR]
    nirRoutes = {
        _route(PRO_CATG=f"{name.removesuffix('_VIS')}_NIR"): [str(tmp_path / f"{name}_NIR.fits")]
        for name in CALIBRATION_NAMES
    }
    recipe, _ = _nod_recipe(log, tmp_path, _with_arcfiles(headers), extraRoutes=nirRoutes)
    recipe.arm = "NIR"
    _patch_nod(recipe, monkeypatch, tmp_path)
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert [entry.exposureSeconds for entry in seen[0]] == [1800.0, 1800.0]


# OFFSET


def test_an_offset_off_frame_is_not_on_target(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: THE SOF HOLDS ONE ON FRAME AND ONE OFF FRAME
    recipe, _ = _offset_recipe(log, tmp_path, _with_arcfiles(_offset_headers([(0.0, 0.0), (20.0, 20.0)], 60000.0)))
    _patch_offset(recipe, monkeypatch, tmp_path)
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: THE ON FRAME IS THE FIRST, SOXS.0. THE OFF FRAME IS STILL AN INPUT OF THE PRODUCT
    assert _names(seen[0]) == ["SOXS.0.fits"]


def test_an_offset_on_frame_of_a_cycle_with_no_trace_leaves_the_list_with_a_warning(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE ON FRAMES ARE SOXS.0 AND SOXS.1, THE OFF FRAMES SOXS.2 AND SOXS.3
    recipe, _ = _offset_recipe(log, tmp_path, _with_arcfiles(TWO_CYCLES))
    _each_offset_cycle(recipe)
    _patch_offset(recipe, monkeypatch, tmp_path)
    _without_trace(recipe, monkeypatch, {(2, "A")})
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _names(seen[0]) == ["SOXS.0.fits"]
    warnings = _left_out_warnings(log)
    assert len(warnings) == 1 and "SOXS.1.fits" in warnings[0]


def test_a_generic_offset_run_counts_its_object_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: GENERICOFFSET TAKES ON AND OFF FROM THE DPR TYPE, AND ITS SKY FRAMES FOLLOW THE OBJECT FRAMES
    headers = _with_arcfiles(
        [
            {**header, "EXPTIME": 300.0}
            for header in generic_offset_run(
                [(0.0, 0.0), (0.0, 0.2), (30.0, 30.0), (30.0, 30.2)], ["OBJECT", "OBJECT", "SKY", "SKY"]
            )
        ]
    )
    recipe, _ = _offset_recipe(
        log,
        tmp_path,
        headers,
        extraRoutes={
            _route(DPR_TYPE="OBJECT", DPR_TECH=OFFSET_TECHNIQUE): [],
            _route(DPR_TYPE="OBJECT", DPR_TECH=GENERIC_TECHNIQUE): [str(path) for path in paths_of(tmp_path, 0, 2)],
            _route(DPR_TYPE="SKY", DPR_TECH=GENERIC_TECHNIQUE): [str(path) for path in paths_of(tmp_path, 2, 4)],
        },
    )
    _patch_offset(recipe, monkeypatch, tmp_path)
    seen = _list_at_each_stack(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert _names(seen[0]) == ["SOXS.0.fits", "SOXS.1.fits"]


def paths_of(tmpPath: Path, first: int, last: int) -> list[Path]:
    """Return the paths `_offset_recipe` writes for frames ``first`` to ``last - 1``."""
    return [tmpPath / f"offset-{index}_pre.fits" for index in range(first, last)]
