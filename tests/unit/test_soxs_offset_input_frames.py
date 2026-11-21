"""The input frames `soxs-offset` accepts: both offset templates, and the sky frames of either (DY-2354)."""

from __future__ import annotations

from typing import Any

import pytest

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import soxs_nod, soxs_offset
from soxspipe.recipes.base_recipe import base_recipe

pytestmark = pytest.mark.unit

CALIBRATIONS = ["DISP_TAB_VIS", "ORDER_TAB_VIS", "DISP_IMAGE_VIS"]


class InputFramesWithSummary:
    """The input-frame collection, reduced to the summary the error path prints."""

    summary = "SUMMARY TABLE"


def _verified(
    recipeClass: type[Any],
    recipeName: str,
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    types: list[str],
    techs: list[str],
) -> None:
    """Run the recipe's own frame verification on frames of the given types and techniques."""
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = recipeName
    recipe.kw = keyword_lookup(log=log, settings={"instrument": "soxs"}).get
    recipe.inputFrames = InputFramesWithSummary()
    monkeypatch.setattr(
        base_recipe, "_verify_input_frames_basics", lambda self: (list(types), list(techs), list(CALIBRATIONS))
    )
    recipe.verify_input_frames()


@pytest.mark.parametrize("technique", ["ECHELLE,SLIT,OFFSET", "ECHELLE,SLIT,GENOFFSET"])
def test_an_offset_sof_accepts_object_and_sky_frames_of_either_offset_template(
    log: Any, monkeypatch: pytest.MonkeyPatch, technique: str
) -> None:
    # ACT / ASSERT: NO ERROR IS RAISED
    _verified(soxs_offset, "soxs-offset", log, monkeypatch, ["OBJECT", "SKY"], [technique])


def test_a_nod_sof_still_refuses_sky_frames(log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    # ACT / ASSERT
    with pytest.raises(TypeError, match="Found a SKY file"):
        _verified(soxs_nod, "soxs-nod", log, monkeypatch, ["OBJECT", "SKY"], ["ECHELLE,SLIT,NODDING"])
