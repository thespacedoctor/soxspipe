"""The shipped ``soxs-mflat`` ``ordexp90`` QC ranges are the ones Dave ruled on DY-2318."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent


def _mflat_ranges(arm: str) -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings["soxs-mflat"][arm]["qc-acceptable-ranges"]


def test_nir_master_flat_accepts_ordexp90_up_to_70000_electrons():
    # ACT
    ranges = _mflat_ranges("nir")

    # ASSERT
    assert ranges["ordexp90"] == [30000, 70000]


def test_vis_master_flat_ordexp90_range_is_unchanged():
    # ACT
    ranges = _mflat_ranges("vis")

    # ASSERT
    assert ranges["ordexp90"] == [10000, 150000]
