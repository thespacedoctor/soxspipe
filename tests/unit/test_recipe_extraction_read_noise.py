"""The read noise a recipe hands to the Horne extraction's sky-line S/N gate (DY-2323)."""

from __future__ import annotations

import pytest
from astropy import units as u

from soxspipe.recipes.base_recipe import base_recipe

pytestmark = pytest.mark.unit


def _recipe(ron: object) -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.detectorParams = {"ron": ron}
    return recipe


def test_the_extraction_read_noise_is_the_resolved_read_noise_in_electrons() -> None:
    assert _recipe(5.5 * u.electron)._extraction_read_noise() == 5.5


def test_the_extraction_read_noise_of_a_plain_detector_default_is_a_float() -> None:
    assert _recipe(3.2)._extraction_read_noise() == 3.2
