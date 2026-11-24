"""A nodding SOF reduced in stare mode counts a nod position that found an object trace as on target (DY-2356).

The A frames and the B frames are extracted apart and the two extractions are combined. A position with no object
trace leaves the combined spectrum, so its frames leave the on-target list, with a warning that names them. The other
position stays. The extraction of a position writes no frame accounting of its own.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from astropy.io import fits

from soxspipe.commonutils.on_target_frames import build_on_target_frames
from soxspipe.recipes import soxs_stare
from soxspipe.recipes.nod_stare_stacking import position_stack
from tests.factories import synthetic_ccd
from tests.factories.on_target_frames import raw_science_header
from tests.integration.test_soxs_stare_nod_mode import MARKER, SOF_NAME, _tail_inputs
from tests.integration.test_soxs_stare_nod_positions import _recipe, _spy_finish, _stacks

pytestmark = pytest.mark.integration

EXPOSURE_SECONDS = 600.0
# THE A FRAMES ARE READ FIRST AND THIRD, THE B FRAMES SECOND AND FOURTH
NAMES_OF = {"A": ("frame-0.fits", "frame-2.fits"), "B": ("frame-1.fits", "frame-3.fits")}


def _listed_recipe(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, fallbacks: tuple[str, ...] = ()
) -> tuple[soxs_stare, list[position_stack]]:
    """Return a recipe that has four on-target frames and two position stacks that know their frame names."""
    recipe, _ = _recipe(log, tmp_path, monkeypatch, fallbacks=fallbacks)
    _spy_finish(recipe, monkeypatch)
    headers = [
        raw_science_header("VIS", mjd=61000.0 + index / 10, arcfile=f"SOXS.{index}.fits", exptime=EXPOSURE_SECONDS)
        for index in range(4)
    ]
    recipe.onTargetFrames = build_on_target_frames(
        headers, [f"frame-{index}.fits" for index in range(4)], recipe.kw, "VIS", log
    )
    stacks = [stack._replace(frameNames=NAMES_OF[stack.position]) for stack in _stacks()]
    return recipe, stacks


def _written_header(tmp_path: Path) -> fits.Header:
    return fits.getheader(tmp_path / f"{SOF_NAME}{MARKER}_SPECTRUM.fits")


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_both_positions_that_found_a_trace_keep_all_their_frames_in_the_combined_spectrum(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, stacks = _listed_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT
    header = _written_header(tmp_path)
    assert header["NCOMBINE"] == 4
    assert header["EXPTIME"] == header["TEXPTIME"] == 4 * EXPOSURE_SECONDS


def test_a_position_with_no_object_trace_leaves_the_list_and_a_warning_names_its_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, stacks = _listed_recipe(log, tmp_path, monkeypatch, fallbacks=("B",))

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT: ONLY THE B FRAMES LEAVE. THE A FRAMES ARE THE SPECTRUM
    header = _written_header(tmp_path)
    assert header["NCOMBINE"] == 2
    assert (header["PROV1"], header["PROV2"]) == ("SOXS.0", "SOXS.2")
    assert header["EXPTIME"] == 2 * EXPOSURE_SECONDS
    leftOut = [message for message in _warnings(log) if "left out of the on-target list" in message]
    assert len(leftOut) == 1
    assert "SOXS.1.fits" in leftOut[0] and "SOXS.3.fits" in leftOut[0] and "SOXS.0.fits" not in leftOut[0]


def test_the_extraction_of_a_position_gets_no_list_and_the_recipe_gets_its_list_back(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe, stacks = _listed_recipe(log, tmp_path, monkeypatch)
    fullList = recipe.onTargetFrames
    seen: list[Any] = []
    realExtract = recipe._extract_stack

    def spy(*args: Any, **kwargs: Any) -> Any:
        seen.append(recipe.onTargetFrames)
        return realExtract(*args, **kwargs)

    monkeypatch.setattr(recipe, "_extract_stack", spy)

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT
    assert seen == [None, None]
    assert recipe.onTargetFrames == fullList


def test_a_stack_of_positions_without_frame_names_leaves_the_list_alone(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: A STACK BUILT WITH NO FRAME NAMES CANNOT SAY WHICH FRAMES TO LEAVE OUT
    recipe, stacks = _listed_recipe(log, tmp_path, monkeypatch, fallbacks=("B",))
    stacks = [stack._replace(frameNames=()) for stack in stacks]

    # ACT
    recipe._reduce_nod_positions(stacks, [], False, _tail_inputs())

    # ASSERT
    assert _written_header(tmp_path)["NCOMBINE"] == 4


def test_a_position_stack_made_without_frame_names_is_still_valid() -> None:
    # ACT
    stack = position_stack("A", [], synthetic_ccd(prepared=True), synthetic_ccd(prepared=True), None)

    # ASSERT
    assert stack.frameNames == ()
