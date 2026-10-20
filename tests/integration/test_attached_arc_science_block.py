"""The observing recipes give the Horne extraction every raw science frame of the block (DY-1666).

The attached-arc test measures the arc from the first frame's start and the last frame's end. These tests read
real prepared frames through each recipe's own reader, capture the headers the recipe hands to
``horne_extraction``, and run the real attached-arc gate on them.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.recipes import soxs_nod, soxs_offset
from tests.factories import synthetic_ccd
from tests.integration.test_observing_recipe_orchestration import (
    RouteCollection,
    _configure_nodding_recipe,
    _patch_nodding_collaborators,
    _write_prepared_frame,
)
from tests.integration.test_soxs_stare_reduction_characterization import (
    STARE_TECHNIQUE,
    _patch_stare,
    _stare_recipe,
)
from tests.unit._horne_vis_scene import KEYWORDS, build_extractor

pytestmark = pytest.mark.integration

MINUTES_PER_DAY = 1440.0
FIRST_MJD = 61217.19708617
FRAME_EXPTIME = 180.0
CADENCE_MINUTES = 4.0
FRAME_COUNT = 4
FIRST_ABSROT = -0.57699
DEGREES_PER_FRAME = 2.0
# THE ARC STARTS 2 MIN AFTER THE LAST FRAME ENDS: 17 MIN AFTER THE FIRST STARTS, FAR OUTSIDE THE 7 MIN LIMIT
ARC_START_MINUTES = (FRAME_COUNT - 1) * CADENCE_MINUTES + FRAME_EXPTIME / 60.0 + 2.0


def _science_overrides(index: int, technique: str) -> dict[str, object]:
    return {
        "DPR_TYPE": "OBJECT",
        "DPR_TECH": technique,
        KEYWORDS["MJDOBS"]: FIRST_MJD + index * CADENCE_MINUTES / MINUTES_PER_DAY,
        KEYWORDS["EXPTIME"]: FRAME_EXPTIME,
        KEYWORDS["ABSROT"]: FIRST_ABSROT + index * DEGREES_PER_FRAME,
        KEYWORDS["PAR_ANG_START"]: -85.867,
        "ESO SEQ CUMOFF Y": 3.0 if index % 2 == 0 else -3.0,
        "HIERARCH ESO SEQ FIXOFF RA": -2.0 if index % 2 == 0 else 0.0,
        "HIERARCH ESO SEQ FIXOFF DEC": 0.0,
    }


def _write_block(tmp_path: Path, technique: str) -> list[Path]:
    return [
        _write_prepared_frame(
            tmp_path / f"science-{index}.fits", seed=600 + index, headerOverrides=_science_overrides(index, technique)
        )
        for index in range(FRAME_COUNT)
    ]


def _arc_header() -> fits.Header:
    header = fits.Header()
    header[KEYWORDS["MJDOBS"]] = FIRST_MJD + ARC_START_MINUTES / MINUTES_PER_DAY
    header[KEYWORDS["EXPTIME"]] = 30.0
    header[KEYWORDS["ABSROT"]] = FIRST_ABSROT + (FRAME_COUNT - 1) * DEGREES_PER_FRAME + 1.0
    header[KEYWORDS["PAR_ANG_START"]] = -85.867
    return header


def _gate_tunes(log: Any, tmp_path: Path, scienceHeaders: list[fits.Header]) -> bool:
    """Run the real attached-arc gate on the block; True when it tunes to the sky lines, False when it skips."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    extractor = build_extractor(log, tmp_path, dict(scienceHeaders[0]))
    extractor.scienceHeaders = scienceHeaders
    with fits.open(extractor.dispersionMap, mode="update") as hdus:
        hdus[0].header.update(_arc_header())
    calls = []
    extractor.tune_wavelength_calibration_to_skylines = lambda df, arm, byOrder=True: calls.append(arm) or df
    extractor._tune_or_skip(pd.DataFrame({"order": [1], "wavelengthMean": [500.0]}), "VIS")
    return bool(calls)


def _patch_extractor(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, object]]:
    captured: list[dict[str, object]] = []

    class RecordingExtractor:
        noddingSequence = ""
        extractionNotes: list = []

        def __init__(self, **kwargs: object) -> None:
            captured.append(kwargs)

        def extract(self) -> tuple[Any, Any, pd.DataFrame, dict[str, int], str]:
            return False, False, pd.DataFrame({"WAVE": [500.0]}), {"1-2": 1}, "extraction.fits"

    monkeypatch.setattr(import_module("soxspipe.commonutils"), "horne_extraction", RecordingExtractor)
    return captured


def _assert_block(scienceHeaders: list[fits.Header]) -> None:
    assert [header["MJD-OBS"] for header in scienceHeaders] == pytest.approx(
        [FIRST_MJD + index * CADENCE_MINUTES / MINUTES_PER_DAY for index in range(FRAME_COUNT)]
    )


@pytest.mark.parametrize(
    ("recipeClass", "technique"),
    [(soxs_nod, "ECHELLE,SLIT,NODDING"), (soxs_offset, "ECHELLE,SLIT,OFFSET")],
)
def test_nod_and_offset_extractions_get_every_frame_of_the_block(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    recipeClass: type,
    technique: str,
) -> None:
    # ARRANGE
    recipe = recipeClass.__new__(recipeClass)
    _configure_nodding_recipe(
        recipe, log=log, tmp_path=tmp_path, objectPaths=_write_block(tmp_path, technique), technique=technique
    )
    recipe.turnOffMP = True
    recipe.twoDMap = str(tmp_path / "DISP_IMAGE_VIS.fits")
    _patch_nodding_collaborators(
        recipe, monkeypatch=monkeypatch, plotPath=tmp_path / "merged.pdf", extractionPath=tmp_path / "extracted.fits"
    )
    recipe.produce_product()
    firstFrame = synthetic_ccd(prepared=True, headerOverrides=_science_overrides(0, technique))
    captured = _patch_extractor(monkeypatch)

    # ACT
    recipe._extract_cycle_spectra(
        firstFrame, firstFrame, firstFrame, firstFrame, firstFrame, firstFrame, False, 1, False
    )

    # ASSERT
    assert captured
    for index, kwargs in enumerate(captured):
        _assert_block(kwargs["scienceHeaders"])
        assert _gate_tunes(log, tmp_path / f"gate-{index}", kwargs["scienceHeaders"]) is False
    # THE FIRST FRAME'S HEADER ALONE WOULD STILL TUNE: THE ARC STARTS 14 MIN AFTER IT ENDS
    assert _gate_tunes(log, tmp_path / "single", [firstFrame.header]) is True


def test_the_science_headers_are_copies_of_the_raw_frames(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    technique = "ECHELLE,SLIT,NODDING"
    recipe = soxs_nod.__new__(soxs_nod)
    _configure_nodding_recipe(
        recipe, log=log, tmp_path=tmp_path, objectPaths=_write_block(tmp_path, technique), technique=technique
    )
    _, captured, _ = _patch_nodding_collaborators(
        recipe, monkeypatch=monkeypatch, plotPath=tmp_path / "merged.pdf", extractionPath=tmp_path / "extracted.fits"
    )

    # ACT
    recipe.produce_product()
    captured["stack"][0]["frames"][0].header[KEYWORDS["MJDOBS"]] = 0.0

    # ASSERT
    _assert_block(recipe.scienceHeaders)


@pytest.mark.parametrize("stackBeforeExtraction", [True, False])
def test_stare_extractions_get_every_frame_of_the_block(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, stackBeforeExtraction: bool
) -> None:
    # ARRANGE
    paths = [str(path) for path in _write_block(tmp_path, STARE_TECHNIQUE)]
    routes = {frozenset({"DPR_TYPE": "OBJECT", "DPR_TECH": STARE_TECHNIQUE}.items()): paths}
    recipe = _stare_recipe(
        log,
        tmp_path,
        extraRoutes=routes,
        recipeSettings={
            "use_flat": False,
            "sky-subtraction": {"subtract_sky": False},
            "stack_before_extraction": stackBeforeExtraction,
        },
    )
    assert isinstance(recipe.inputFrames, RouteCollection)
    captured = _patch_stare(recipe, monkeypatch, tmp_path)

    # ACT
    recipe.produce_product()

    # ASSERT
    extractions = captured["extractor"]
    assert len(extractions) == (1 if stackBeforeExtraction else FRAME_COUNT)
    for index, kwargs in enumerate(extractions):
        _assert_block(kwargs["scienceHeaders"])
        assert _gate_tunes(log, tmp_path / f"gate-{index}", kwargs["scienceHeaders"]) is False
