"""The read noise a real recipe hands to the Horne extraction's sky-line S/N gate (DY-2323).

The recipe's own constructor runs for real, with the shipped SOXS settings, so its detector parameters are the ones the
read-noise lookup works on.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pytest
from astropy import units as u

from soxspipe.commonutils import detector_lookup
from tests.integration.test_soxs_stare_flat_slit_gradient_seam import _real_nir_stare_recipe

pytestmark = pytest.mark.integration


def _recipe_with_detector_parameters(log: Any, monkeypatch: pytest.MonkeyPatch, tmpPath: Path) -> Any:
    """The real NIR stare recipe, with the detector parameters that its frame preparation looks up."""
    recipe = _real_nir_stare_recipe(log, monkeypatch, tmpPath)
    recipe.detectorParams = detector_lookup(log=log, settings=recipe.settings).get("NIR")
    return recipe


def test_the_extraction_read_noise_is_the_resolved_read_noise_in_electrons(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE: THE LOOKUP SETS THE READ NOISE AS A QUANTITY
    recipe = _recipe_with_detector_parameters(log, monkeypatch, tmp_path)
    recipe.detectorParams["ron"] = 5.5 * u.electron

    # ACT
    readNoise = recipe._extraction_read_noise()

    # ASSERT
    assert readNoise == 5.5
    assert isinstance(readNoise, float)


def test_the_extraction_read_noise_of_the_detector_default_is_the_default_as_a_float(
    log: Any, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE
    recipe = _recipe_with_detector_parameters(log, monkeypatch, tmp_path)
    default = recipe.detectorParams["ron"]

    # ACT
    readNoise = recipe._extraction_read_noise()

    # ASSERT
    assert readNoise == float(getattr(default, "value", default))
    assert readNoise > 0
