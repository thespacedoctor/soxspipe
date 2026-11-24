"""The `extraction_stacking` mixin combines per-position extractions in the nod and stare recipes (DY-2290).

The stare recipe reduces the A and B frames of a nodding SOF in stare mode as two stacks, and combines the two
extractions the way the nod recipe combines its cycles. These tests run `stack_extractions` on a stare-marked recipe.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy import units as u
from astropy.table import Table

from soxspipe.recipes import soxs_nod, soxs_stare
from soxspipe.recipes.extraction_stacking import extraction_stacking
from tests.factories import pipeline_settings, qc_table, synthetic_ccd

pytestmark = pytest.mark.unit

SOF_NAME = "20251212T045200_VIS_1X1_1_NOD_OBJ_SOXS"
MARKER = "_STAREMODE"


def _recipe(recipeClass: type, log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, marker: str) -> Any:
    """Return an unconstructed recipe that can run `stack_extractions`, with the SNR QCs passed through."""
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"),
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = "soxs-stare" if recipeClass is soxs_stare else "soxs-nod"
    recipe.sofName = SOF_NAME
    recipe.productNameMarker = marker
    recipe.settings = pipeline_settings(tmpPath)
    recipe.productDir = str(tmpPath)
    _, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.masterHeaderFrame = synthetic_ccd(shape=(3, 3), prepared=True)
    recipe.update_fits_keywords = lambda **_: None
    return recipe


@pytest.fixture
def stare(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> soxs_stare:
    return _recipe(soxs_stare, log, tmp_path, monkeypatch, MARKER)


@pytest.fixture
def nod(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> soxs_nod:
    return _recipe(soxs_nod, log, tmp_path, monkeypatch, "")


def _spectrum(flux: float, sky: float = 1.0, variance: float = 4.0) -> pd.DataFrame:
    """Return an order-merged spectrum of one wavelength, shaped like the extraction's output."""
    return pd.DataFrame(
        {
            "WAVE": [500.0 * u.nm],
            "FLUX_COUNTS": [flux],
            "VARIANCE": [variance],
            "SKY_COUNTS": [sky],
            "SNR": [0.0],
            "FLUX_DENSITY_COUNTS": [flux],
        }
    )


def test_the_stare_and_nod_recipes_share_the_stacking_mixin() -> None:
    assert issubclass(soxs_stare, extraction_stacking)
    assert issubclass(soxs_nod, extraction_stacking)


def test_a_stare_recipe_writes_the_combined_spectrum_under_its_marked_name(stare: soxs_stare) -> None:
    # ACT
    _, filePath = stare.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    assert Path(filePath).name == f"{SOF_NAME}{MARKER}_SPECTRUM.fits"
    assert Path(filePath).exists()
    assert stare.filenameTemplate == f"{SOF_NAME}{MARKER}.fits"
    assert f"{SOF_NAME}{MARKER}_SPECTRUM.fits" in set(stare.products["file_name"])


def test_a_nod_recipe_writes_the_combined_spectrum_under_the_sof_name(nod: soxs_nod) -> None:
    # ACT
    _, filePath = nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    assert Path(filePath).name == f"{SOF_NAME}_SPECTRUM.fits"


def test_two_extractions_combine_to_their_mean(stare: soxs_stare) -> None:
    # ACT
    stacked, filePath = stare.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([15.0])
    assert list(stacked["VARIANCE"]) == pytest.approx([2.0])
    assert list(Table.read(filePath)["FLUX_COUNTS"]) == pytest.approx([15.0])


def test_the_sky_counts_column_is_carried_into_the_combined_spectrum(stare: soxs_stare) -> None:
    # ACT
    stacked, _ = stare.stack_extractions([_spectrum(10.0, sky=2.0), _spectrum(20.0, sky=6.0)])

    # ASSERT
    assert list(stacked["SKY_COUNTS"]) == pytest.approx([4.0])


def test_a_position_with_no_spectrum_is_left_out_of_the_combination(stare: soxs_stare) -> None:
    # ACT
    stacked, _ = stare.stack_extractions([None, _spectrum(20.0)])

    # ASSERT
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([20.0])
    assert list(stacked["VARIANCE"]) == pytest.approx([4.0])


def test_no_spectrum_at_either_position_is_an_error(stare: soxs_stare, log: Any) -> None:
    # ACT / ASSERT
    with pytest.raises(Exception, match="No object trace was found"):
        stare.stack_extractions([None, None])
    assert any(level == "error" for level, _ in log.messages)
