"""Directory listings that decide a reduction input must not depend on filesystem order (DY-48)."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from soxspipe.commonutils.set_of_files import set_of_files
from tests.factories import pipeline_settings, raw_fits

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    "listingOrder", [sorted, lambda names: sorted(names, reverse=True)], ids=["forward", "reverse"]
)
def test_directory_sof_picks_the_same_disp_map_whatever_the_listing_order(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, listingOrder
) -> None:
    """When two dispersion maps exist for one arm, the one that sorts last wins on every filesystem."""
    # ARRANGE
    framesDir = tmp_path / "frames"
    framesDir.mkdir()
    raw_fits(framesDir / "a.fits")
    (framesDir / "VIS_DISP_MAP_A.csv").write_text("x\n1\n", encoding="utf-8")
    (framesDir / "VIS_DISP_MAP_B.csv").write_text("x\n1\n", encoding="utf-8")
    settings = pipeline_settings(
        tmp_path, overrides={"summary-keys": {"default": ["DPR_TYPE"], "verbose": [], "nodding_extras": []}}
    )
    realListdir = os.listdir
    monkeypatch.setattr(os, "listdir", lambda path=".": listingOrder(realListdir(path)))

    # ACT
    _collection, supplementary = set_of_files(
        log=log, settings=settings, inputFrames=str(framesDir), verbose=False
    ).get()

    # ASSERT
    assert supplementary["VIS"]["DISP_MAP"] == str(framesDir / "VIS_DISP_MAP_B.csv")
