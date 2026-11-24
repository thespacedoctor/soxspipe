"""The shipped ``soxs-disp-solution`` ``detlines-frac`` QC ranges are the ones Dave ruled on DY-2207."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent


def _disp_solution_ranges(arm: str) -> dict[str, list[float]]:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return settings["soxs-disp-solution"][arm]["qc-acceptable-ranges"]


def test_nir_dispersion_solution_accepts_detlines_frac_down_to_0_85() -> None:
    # ACT
    ranges = _disp_solution_ranges("nir")

    # ASSERT
    assert ranges["detlines-frac"] == [0.85, 1.01]


def test_vis_dispersion_solution_detlines_frac_range_is_unchanged() -> None:
    # ACT
    ranges = _disp_solution_ranges("vis")

    # ASSERT
    assert ranges["detlines-frac"] == [0.9, 1.01]
