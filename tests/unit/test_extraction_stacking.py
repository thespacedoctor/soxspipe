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
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.on_target_frames import build_on_target_frames
from soxspipe.recipes import soxs_nod, soxs_stare
from soxspipe.recipes.extraction_stacking import extraction_stacking
from tests.factories import pipeline_settings, qc_table, synthetic_ccd
from tests.factories.on_target_frames import raw_science_header

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


def _on_target_frames(recipe: Any, count: int = 2) -> tuple[Any, ...]:
    """Return the on-target list of ``count`` VIS frames of 600 s, named ``SOXS.<n>.fits``."""
    headers = [
        raw_science_header("VIS", mjd=61000.0 + index / 10, arcfile=f"SOXS.{index}.fits", exptime=600.0)
        for index in range(count)
    ]
    names = [f"raw-{index}.fits" for index in range(count)]
    return build_on_target_frames(headers, names, recipe.kw, "VIS", recipe.log)


@pytest.fixture
def listed_nod(nod: soxs_nod) -> soxs_nod:
    """A nod recipe whose master header is a raw frame of 600 s and that has an on-target list of two frames."""
    from soxspipe.commonutils.keyword_lookup import keyword_lookup

    nod.kw = keyword_lookup(log=nod.log, settings=nod.settings).get
    nod.masterHeaderFrame.header["EXPTIME"] = 600.0
    nod.onTargetFrames = _on_target_frames(nod)
    return nod


def test_the_combined_spectrum_carries_the_frame_accounting_cards_of_the_on_target_list(
    listed_nod: soxs_nod,
) -> None:
    # ACT
    _, filePath = listed_nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    header = fits.getheader(filePath)
    assert header["EXPTIME"] == header["TEXPTIME"] == 1200.0
    assert header["NCOMBINE"] == 2
    assert (header["PROV1"], header["PROV2"]) == ("SOXS.0", "SOXS.1")
    assert header["OBID1"] == 42
    assert header["PROG_ID"] == "115.ABCD.001"
    assert header["MJD-OBS"] == 61000.0
    assert header["MJD-END"] == pytest.approx(61000.1 + 600.0 / 86400.0, abs=1e-9)


def test_stamping_the_combined_spectrum_leaves_the_exposure_time_of_the_master_header_alone(
    listed_nod: soxs_nod,
) -> None:
    # ACT
    listed_nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT: FLUX CALIBRATION DIVIDES BY THIS EXPOSURE TIME, WHICH IS THAT OF ONE FRAME
    assert listed_nod.masterHeaderFrame.header["EXPTIME"] == 600.0
    assert "NCOMBINE" not in listed_nod.masterHeaderFrame.header


def test_the_combined_spectrum_has_the_target_position_and_the_uncalibrated_flux_cards(listed_nod: soxs_nod) -> None:
    # ARRANGE: THE MASTER HEADER HAS THE POINTING, AND THE EARLIEST FRAME HAS THE TARGET
    target = {"HIERARCH ESO TEL TARG ALPHA": 3019.8, "HIERARCH ESO TEL TARG DELTA": -331501.08}
    headers = [
        raw_science_header("VIS", mjd=61000.0, arcfile="SOXS.0.fits", exptime=600.0, **target),
        raw_science_header("VIS", mjd=61000.1, arcfile="SOXS.1.fits", exptime=600.0),
    ]
    listed_nod.onTargetFrames = build_on_target_frames(
        headers, ["raw-0.fits", "raw-1.fits"], listed_nod.kw, "VIS", listed_nod.log
    )
    listed_nod.masterHeaderFrame.header["RA"] = 7.9
    listed_nod.masterHeaderFrame.header["DEC"] = -33.0

    # ACT
    _, filePath = listed_nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    header = fits.getheader(filePath)
    assert header["RA"] == pytest.approx(7.5825, abs=1e-9)
    assert header["DEC"] == pytest.approx(-33.2503, abs=1e-9)
    assert (header["RADESYS"], header["EQUINOX"]) == ("FK5", 2000.0)
    assert (header["FLUXCAL"], header["TOT_FLUX"], header["FLUXERR"]) == ("UNCALIBRATED", False, -1.0)


def test_the_not_flattened_stack_is_uncalibrated_too(listed_nod: soxs_nod) -> None:
    # ACT
    _, filePath = listed_nod.stack_extractions([_spectrum(10.0)], notFlattened=True)

    # ASSERT
    header = fits.getheader(filePath)
    assert (header["FLUXCAL"], header["TOT_FLUX"], header["FLUXERR"]) == ("UNCALIBRATED", False, -1.0)


def test_a_recipe_with_no_on_target_list_writes_no_frame_accounting_cards(nod: soxs_nod) -> None:
    # ACT
    _, filePath = nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ASSERT
    header = fits.getheader(filePath)
    assert not {"NCOMBINE", "TEXPTIME", "MJD-END", "PROV1", "OBID1"} & set(header)


def test_the_not_flattened_stack_carries_no_frame_accounting_cards(listed_nod: soxs_nod) -> None:
    # ACT
    _, filePath = listed_nod.stack_extractions([_spectrum(10.0)], notFlattened=True)

    # ASSERT
    assert "NCOMBINE" not in fits.getheader(filePath)


def test_the_merged_spectrum_header_is_the_stamped_header_the_flux_calibration_rewrites(listed_nod: soxs_nod) -> None:
    # ARRANGE
    listed_nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ACT
    header = listed_nod._merged_spectrum_header()

    # ASSERT
    assert header["EXPTIME"] == 1200.0
    assert header is not listed_nod.masterHeaderFrame.header


ACCOUNTING_KEYWORDS = ("NCOMBINE", "PROV", "OBID", "MJD-END", "MJD-OBS", "EXPTIME", "TEXPTIME", "PROG_ID", "PROGID")


def test_the_combined_spectrum_has_no_frame_accounting_finding_from_the_phase_3_checker(listed_nod: soxs_nod) -> None:
    from soxspipe.commonutils.phase3_check import check_phase3

    # ARRANGE
    _, filePath = listed_nod.stack_extractions([_spectrum(10.0), _spectrum(20.0)])

    # ACT
    findings = check_phase3(filePath)

    # ASSERT: THE OTHER KEYWORDS OF THE STANDARD ARE THE BUSINESS OF OTHER ISSUES
    accounting = [
        finding.message for finding in findings if any(keyword in finding.message for keyword in ACCOUNTING_KEYWORDS)
    ]
    assert accounting == []
