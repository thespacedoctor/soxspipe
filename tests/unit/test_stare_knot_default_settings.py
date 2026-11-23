"""The shipped sky-fit knot spacing of the stare recipes: finer for the SOXS VIS arm only (DY-1977)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
STARTING_KEY = "starting_points_per_knot"
MINIMUM_KEY = "min_points_per_knot"
# (FILE, RECIPE, ARM, STARTING POINTS PER KNOT, MINIMUM POINTS PER KNOT)
CASES = [
    ("soxs_default_settings.yaml", "soxs-stare", "vis", 500, 20),
    ("soxs_default_settings.yaml", "soxs-stare-std", "vis", 500, 77),
    ("soxs_default_settings.yaml", "soxs-stare", "nir", 1000, 15),
    ("soxs_default_settings.yaml", "soxs-stare-std", "nir", 1000, 15),
    ("xsh_default_settings.yaml", "soxs-stare", "uvb", 300, 11),
    ("xsh_default_settings.yaml", "soxs-stare", "vis", 100, 11),
    ("xsh_default_settings.yaml", "soxs-stare", "nir", 25, 9),
]


@pytest.mark.parametrize(("fileName", "recipe", "arm", "starting", "minimum"), CASES)
def test_each_stare_arm_ships_its_knot_spacing(
    fileName: str, recipe: str, arm: str, starting: int, minimum: int
) -> None:
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    skySettings = settings[recipe][arm]["sky-subtraction"]

    assert (skySettings[STARTING_KEY], skySettings[MINIMUM_KEY]) == (starting, minimum)
