"""Integration contracts for keeping prepared (`SXSPRE`) frames out of the raw frame index (DY-1703)."""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
import pytest
from astropy.io import fits

from soxspipe.commonutils.data_organiser import _harvest_fits_headers
from soxspipe.commonutils.toolkit import PREPARED_FRAME_KEYWORD
from tests.factories import harvestable_raw_fits
from tests.integration.test_data_organiser_workspace import _organiser_header_keywords

pytestmark = pytest.mark.integration

FILTER_KEYS = ["mjd-obs", "eso dpr catg", "eso dpr tech", "eso dpr type", "exptime"]


def _prepared_copy(rawPath: Path, destination: Path) -> Path:
    """Copy a raw frame and stamp the marker `_write_prepared_frame` sets."""
    with fits.open(rawPath) as hdus:
        hdus[0].header[PREPARED_FRAME_KEYWORD] = ("2024-01-02T03:04:05.678901", "UTC timestamp")
        hdus.writeto(destination)
    return destination


def _harvest(tmp_path: Path, log: logging.Logger, names: list[str]) -> pd.DataFrame:
    kw, keywords = _organiser_header_keywords(log)
    return _harvest_fits_headers(
        batch=names,
        log=log,
        pathToDirectory=str(tmp_path),
        keywords=keywords,
        filterKeys=FILTER_KEYS,
        instrument="SOXS",
        kw=kw,
    )


def test_harvest_keeps_a_raw_frame_that_has_no_sxspre_keyword(tmp_path, log) -> None:
    # ARRANGE
    harvestable_raw_fits(tmp_path / "raw.fits")

    # ACT
    harvested = _harvest(tmp_path, log, ["raw.fits"])

    # ASSERT
    assert [Path(name).name for name in harvested["file"]] == ["raw.fits"]
    assert "sxspre" not in harvested.columns


def test_harvest_leaves_out_a_frame_carrying_sxspre_and_keeps_its_raw_sibling(tmp_path, log) -> None:
    # ARRANGE
    rawPath = harvestable_raw_fits(tmp_path / "raw.fits")
    _prepared_copy(rawPath, tmp_path / "raw_pre.fits")

    # ACT
    harvested = _harvest(tmp_path, log, ["raw.fits", "raw_pre.fits"])

    # ASSERT
    assert [Path(name).name for name in harvested["file"]] == ["raw.fits"]
    assert "sxspre" not in harvested.columns


def test_harvest_returns_no_rows_when_every_frame_in_the_batch_is_prepared(tmp_path, log) -> None:
    # ARRANGE
    rawPath = harvestable_raw_fits(tmp_path / "source.fits")
    _prepared_copy(rawPath, tmp_path / "only_pre.fits")

    # ACT
    harvested = _harvest(tmp_path, log, ["only_pre.fits"])

    # ASSERT
    assert len(harvested.index) == 0
