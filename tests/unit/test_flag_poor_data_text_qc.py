"""`flag_poor_data` must range-check numeric QCs when the table also holds a text QC (DY-1608 `UG SHIFT SOURCE`)."""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import base_recipe
from tests.factories import StubInputFrames, noise_frame_pair

pytestmark = pytest.mark.unit


def test_flag_poor_data_range_checks_numeric_qcs_beside_a_text_qc(log: Any) -> None:
    # ARRANGE
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.recipeName = "soxs-stare"
    recipe.dateObs = "2026-07-05T05:26:29"
    recipe.inputFrames = StubInputFrames(noise_frame_pair())
    recipe.qc = pd.DataFrame(
        [
            {"soxspipe_recipe": "soxs-stare", "qc_name": "UG SHIFT SOURCE", "qc_value": "model", "qc_order": None},
        ]
    )
    recipe.inst = "SOXS"
    recipe.detectorTemp = 82.5
    recipe.cptemp = 11.2
    recipe.recipeSettings = {"qc-acceptable-ranges": {"DETECTOR-TEMP": [70.0, 80.0]}}

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    flags = recipe.qc.set_index("qc_name")["qc_flag"]
    assert flags["DETECTOR TEMP"] == "fail"
    assert flags["UG SHIFT SOURCE"] == "pass"
