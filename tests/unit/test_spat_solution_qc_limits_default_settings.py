"""The SOXS NIR ``soxs-spat-solution`` block sets the residual QC limits Dave ruled on DY-1609."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent


def _spat_solution_ranges(arm: str) -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings["soxs-spat-solution"][arm]["qc-acceptable-ranges"]


def test_nir_spatial_solution_accepts_residual_scatter_up_to_the_ruled_limits():
    # ACT
    ranges = _spat_solution_ranges("nir")

    # ASSERT
    assert ranges["xy-res-sd"] == [0, 0.5]
    assert ranges["y-res-sd"] == [0, 0.2]


def test_vis_spatial_solution_residual_limits_are_unchanged():
    # ACT
    ranges = _spat_solution_ranges("vis")

    # ASSERT
    assert ranges["xy-res-sd"] == [0, 0.5]
    assert ranges["y-res-sd"] == [0, 0.48]
