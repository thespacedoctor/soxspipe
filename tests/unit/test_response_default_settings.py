"""The shipped ``soxs-response`` settings configure the spline fit for every arm, and old workspaces still run."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.response_fit import resolve_response_settings

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
REQUIRED_KEYS = {
    "bin_width_nm",
    "knot_spacing_nm",
    "clip_sigma_low",
    "clip_sigma_high",
    "max_iteration",
    "min_bin_samples",
    "telluric_regions",
    "stellar_lines",
    "line_half_width_nm",
    "legacy_poly_order",
    "selfcal_bands",
    "clip_frac_max",
    "max_gap_nm_max",
    "selfcal_res_max",
}
CASES = [
    ("soxs_default_settings.yaml", "vis"),
    ("soxs_default_settings.yaml", "nir"),
    ("xsh_default_settings.yaml", "uvb"),
    ("xsh_default_settings.yaml", "vis"),
    ("xsh_default_settings.yaml", "nir"),
    ("test_settings_xsh.yaml", "uvb"),
    ("test_settings_xsh.yaml", "vis"),
    ("test_settings_xsh.yaml", "nir"),
]


@pytest.mark.parametrize(("fileName", "arm"), CASES)
def test_each_arm_carries_every_response_fit_key(fileName: str, arm: str) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ACT
    armSettings = settings["soxs-response"][arm]

    # ASSERT
    assert set(armSettings) == REQUIRED_KEYS


def test_an_old_workspace_block_is_filled_from_the_packaged_defaults() -> None:
    # ARRANGE: A WORKSPACE YAML WRITTEN BEFORE THE SPLINE FIT HAD ONLY THE POLYNOMIAL KEYS
    oldBlock = {"vis": {"max_iteration": 5, "poly_order": 4}}

    # ACT
    resolved = resolve_response_settings(oldBlock, "VIS", "soxs")

    # ASSERT
    assert set(resolved) == REQUIRED_KEYS
    assert resolved["knot_spacing_nm"] == 40
    assert resolved["max_iteration"] == 5
    assert resolved["legacy_poly_order"] == 4


def test_workspace_values_override_the_packaged_defaults() -> None:
    # ARRANGE
    block = {"nir": {"knot_spacing_nm": 90, "legacy_poly_order": 3}}

    # ACT
    resolved = resolve_response_settings(block, "NIR", "soxs")

    # ASSERT
    assert resolved["knot_spacing_nm"] == 90
    assert resolved["legacy_poly_order"] == 3
    assert resolved["bin_width_nm"] == 5


def test_a_missing_arm_block_resolves_to_the_packaged_defaults() -> None:
    # ACT
    resolved = resolve_response_settings({}, "UVB", "xsh")

    # ASSERT
    assert set(resolved) == REQUIRED_KEYS
