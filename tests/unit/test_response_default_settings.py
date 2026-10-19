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
    "smoothing_penalty",
    "clip_sigma_low",
    "clip_sigma_high",
    "max_iteration",
    "min_bin_samples",
    "telluric_regions",
    "telluric_wing_max_nm",
    "telluric_wing_ref_nm",
    "telluric_wing_threshold",
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


@pytest.mark.parametrize("fileName", ["soxs_default_settings.yaml", "xsh_default_settings.yaml"])
def test_the_smoothing_penalty_is_enabled_for_the_nir_arm_only(fileName: str) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())["soxs-response"]

    # ACT
    penalty = {arm: block["smoothing_penalty"] for arm, block in settings.items()}

    # ASSERT
    assert penalty["nir"] > 0
    assert all(value == 0 for arm, value in penalty.items() if arm != "nir")


@pytest.mark.parametrize("fileName", ["soxs_default_settings.yaml", "xsh_default_settings.yaml"])
def test_wing_growth_is_enabled_for_the_nir_arm_only(fileName: str) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())["soxs-response"]

    # ACT
    wingReach = {arm: block["telluric_wing_max_nm"] for arm, block in settings.items()}

    # ASSERT
    assert wingReach["nir"] > 0
    assert all(reach == 0 for arm, reach in wingReach.items() if arm != "nir")


@pytest.mark.parametrize("fileName", ["soxs_default_settings.yaml", "xsh_default_settings.yaml"])
def test_the_nir_telluric_regions_cover_the_925_to_980_nm_water_band(fileName: str) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())["soxs-response"]

    # ACT
    regions = settings["nir"]["telluric_regions"]

    # ASSERT
    assert any(lower <= 925 and upper >= 980 for lower, upper in regions)
