"""The header of the flux-calibrated nod and offset spectrum is the header of the merged spectrum (DY-2356).

The flux calibration rewrites ``<SOF stem>_SPECTRUM.fits``. It was given the header of the first frame, which the recipe
had stamped only with the raw frames of its own cycle, or not at all after the frames were aligned and stacked. The
calibrated file then listed no raw frame, or the two frames of the first cycle. It now gets the header the merged
spectrum was written with, which lists every raw input.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import soxs_nod, soxs_offset
from tests.factories import pipeline_settings, qc_table, response_product_fits, synthetic_ccd
from tests.factories.on_target_frames import raw_frame_summary, raw_science_header

pytestmark = pytest.mark.integration

NOD_NAMES = ["nod-A1.fits", "nod-B1.fits", "nod-A2.fits", "nod-B2.fits"]
OFFSET_NAMES = ["off-ON1.fits", "off-OFF1.fits", "off-ON2.fits", "off-OFF2.fits"]


def _frame(index: int) -> Any:
    """Return a prepared frame whose ``ARCFILE`` differs from its file name, with an exposure time."""
    header = raw_science_header("VIS", mjd=61000.0 + index / 10, arcfile=f"SOXS.2026-01-02T03:04:0{index}.fits")
    return synthetic_ccd(prepared=True, headerOverrides={key: header[key] for key in header})


def _recipe(recipeClass: type, tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, names: list[str]) -> Any:
    """Return a recipe that has stacked its extractions with the real writers, ready to flux calibrate."""
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"), "add_snr_efficiency_qcs", lambda **kwargs: kwargs["qcTable"]
    )
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.inst = "SOXS"
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.recipeName = "soxs-offset" if recipeClass is soxs_offset else "soxs-nod"
    recipe.recipeSettings = {}
    recipe.sofName = "OBJECT_VIS"
    recipe.startNightDate = "2024-01-02"
    recipe.productDir = str(tmp_path)
    recipe.imageType = "OBJECT"
    recipe.detectorParams = {"extinction": "extinction.dat"}
    recipe.rawFrames = raw_frame_summary(names)
    recipe.debug = False
    recipe.qc = qc_table().iloc[0:0].copy()
    _, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.masterHeaderFrame = _frame(0)
    return recipe


def _spectrum() -> pd.DataFrame:
    return pd.DataFrame({"WAVE": [500.0 * u.nm], "FLUX_COUNTS": [10.0], "VARIANCE": [4.0], "SNR": [0.0]})


def _capture_flux_calibration(monkeypatch: pytest.MonkeyPatch, recipeModule: str) -> list[dict[str, Any]]:
    captured: list[dict[str, Any]] = []

    class FakeFluxCalibrator:
        def __init__(self, **kwargs: Any) -> None:
            captured.append(kwargs)

        def calibrate(self) -> tuple[str, pd.DataFrame]:
            return "calibrated.fits", pd.DataFrame()

    monkeypatch.setattr(import_module("soxspipe.commonutils.flux_calibration"), "flux_calibration", FakeFluxCalibrator)
    monkeypatch.setattr(import_module(recipeModule), "get_calibrations_path", lambda **_: "/cal")
    return captured


def _raw_names(header: Any) -> list[str]:
    return [
        header[f"ESO PRO REC1 RAW{index} NAME"] for index in range(1, 9) if f"ESO PRO REC1 RAW{index} NAME" in header
    ]


def test_the_calibrated_nod_spectrum_lists_every_raw_frame_after_the_frames_were_stacked(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: ALIGNED FRAMES THAT NO KEYWORD UPDATE EVER TOUCHED, AS IN THE STACKED PATH
    recipe = _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES)
    captured = _capture_flux_calibration(monkeypatch, "soxspipe.recipes.soxs_nod")
    stacked, _ = recipe.stack_extractions([_spectrum(), _spectrum()])
    framesA, framesB = [_frame(1), _frame(3)], [_frame(2), _frame(4)]
    responsePath = str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))

    # ACT
    recipe._flux_calibrate_nod_stack(responsePath, stacked, framesA, framesB)

    # ASSERT
    assert len(captured) == 1
    assert _raw_names(captured[0]["header"]) == NOD_NAMES


def test_the_calibrated_nod_spectrum_lists_every_raw_frame_after_the_cycles_were_extracted_one_by_one(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: THE FIRST A FRAME WAS STAMPED WITH THE TWO FRAMES OF ITS OWN CYCLE ONLY
    recipe = _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES)
    captured = _capture_flux_calibration(monkeypatch, "soxspipe.recipes.soxs_nod")
    framesA, framesB = [_frame(1), _frame(3)], [_frame(2), _frame(4)]
    recipe.update_fits_keywords(frame=framesA[0], rawFrames=["nod-A1.fits", "nod-B1.fits"])
    stacked, _ = recipe.stack_extractions([_spectrum(), _spectrum()])
    responsePath = str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))

    # ACT
    recipe._flux_calibrate_nod_stack(responsePath, stacked, framesA, framesB)

    # ASSERT
    assert _raw_names(captured[0]["header"]) == NOD_NAMES


def test_the_calibrated_offset_spectrum_lists_every_raw_frame_the_off_frames_included(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _recipe(soxs_offset, tmp_path, log, monkeypatch, OFFSET_NAMES)
    captured = _capture_flux_calibration(monkeypatch, "soxspipe.recipes.soxs_offset")
    stacked, _ = recipe.stack_extractions([_spectrum()])
    framesON = [_frame(1), _frame(3)]
    responsePath = str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))

    # ACT
    recipe._flux_calibrate_stack(responsePath, stacked, framesON)

    # ASSERT
    assert _raw_names(captured[0]["header"]) == OFFSET_NAMES


def test_the_merged_nod_spectrum_lists_every_raw_frame(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES)

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum(), _spectrum()])

    # ASSERT
    assert _raw_names(fits.getheader(filePath)) == NOD_NAMES


def test_the_offset_spectrum_has_no_frame_accounting_finding_from_the_phase_3_checker(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    from soxspipe.commonutils.on_target_frames import build_on_target_frames
    from soxspipe.commonutils.phase3_check import check_phase3

    # ARRANGE: TWO ON FRAMES, EACH WITH A DIFFERENT FILE NAME AND ARCFILE
    productDir = tmp_path / "offset"
    productDir.mkdir()
    recipe = _recipe(soxs_offset, productDir, log, monkeypatch, OFFSET_NAMES)
    frames = [_frame(1), _frame(3)]
    recipe.onTargetFrames = build_on_target_frames(
        [frame.header for frame in frames], ["off-ON1.fits", "off-ON2.fits"], recipe.kw, "VIS", log
    )

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])
    findings = check_phase3(filePath)

    # ASSERT: THE OTHER KEYWORDS OF THE STANDARD ARE THE BUSINESS OF OTHER ISSUES
    header = fits.getheader(filePath)
    assert (header["NCOMBINE"], header["PROV1"], header["ESO SEQ ARM"]) == (2, "SOXS.2026-01-02T03:04:01", "VIS")
    keywords = ("NCOMBINE", "PROV", "OBID", "MJD-END", "MJD-OBS", "EXPTIME", "PROG_ID", "PROGID")
    assert [f.message for f in findings if any(keyword in f.message for keyword in keywords)] == []


TARGET_CARDS = {"HIERARCH ESO TEL TARG ALPHA": 3019.8, "HIERARCH ESO TEL TARG DELTA": -331501.08}
TARGET_KEYWORDS = ("FLUXCAL", "FLUXERR", "TOT_FLUX", "RADESYS", "EQUINOX", "RA ", "DEC ", "RA is", "DEC is", "OBJECT")


def _target_frame(index: int) -> Any:
    """Return a prepared frame with the pointing, the target cards and an object name."""
    header = raw_science_header(
        "VIS",
        mjd=61000.0 + index / 10,
        arcfile=f"SOXS.2026-01-02T03:04:0{index}.fits",
        RA=7.9,
        DEC=-33.0,
        OBJECT="SN2026ob",
        **(TARGET_CARDS if index == 1 else {}),
    )
    return synthetic_ccd(prepared=True, headerOverrides={key: header[key] for key in header})


def _with_target_list(recipe: Any, log: Any, names: list[str]) -> Any:
    """Give the recipe the on-target list of two frames, the first of which has the target cards."""
    from soxspipe.commonutils.on_target_frames import build_on_target_frames

    frames = [_target_frame(1), _target_frame(3)]
    recipe.masterHeaderFrame = _target_frame(0)
    recipe.onTargetFrames = build_on_target_frames([f.header for f in frames], names, recipe.kw, "VIS", log)
    return recipe


@pytest.mark.parametrize(
    ("recipeClass", "names", "onNames"),
    [
        (soxs_nod, NOD_NAMES, ["nod-A1.fits", "nod-A2.fits"]),
        (soxs_offset, OFFSET_NAMES, ["off-ON1.fits", "off-ON2.fits"]),
    ],
)
def test_the_merged_nod_and_offset_spectrum_has_the_target_position_and_the_uncalibrated_flux_cards(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch, recipeClass: type, names: list[str], onNames: list[str]
) -> None:
    from soxspipe.commonutils.phase3_check import check_phase3

    # ARRANGE
    recipe = _with_target_list(_recipe(recipeClass, tmp_path, log, monkeypatch, names), log, onNames)

    # ACT
    _, filePath = recipe.stack_extractions([_spectrum()])
    findings = check_phase3(filePath)

    # ASSERT: THE POINTING IS GONE, AND THE OBJECT NAME IS THE RAW ONE
    header = fits.getheader(filePath)
    assert header["RA"] == pytest.approx(7.5825, abs=1e-9)
    assert header["DEC"] == pytest.approx(-33.2503, abs=1e-9)
    assert (header["RADESYS"], header["EQUINOX"], header["OBJECT"]) == ("FK5", 2000.0, "SN2026ob")
    assert (header["FLUXCAL"], header["TOT_FLUX"], header["FLUXERR"]) == ("UNCALIBRATED", False, -1.0)
    failures = [f.message for f in findings if f.severity == "failure"]
    # THE EXTENSION HEADER CARDS ARE THE BUSINESS OF ANOTHER ISSUE
    assert [m for m in failures if "extension" not in m and any(keyword in m for keyword in TARGET_KEYWORDS)] == []


def test_the_header_given_to_the_nod_flux_calibration_has_the_target_position(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    recipe = _with_target_list(
        _recipe(soxs_nod, tmp_path, log, monkeypatch, NOD_NAMES), log, ["nod-A1.fits", "nod-A2.fits"]
    )
    captured = _capture_flux_calibration(monkeypatch, "soxspipe.recipes.soxs_nod")
    stacked, _ = recipe.stack_extractions([_spectrum(), _spectrum()])
    responsePath = str(response_product_fits(tmp_path / "RESP_TAB_VIS.fits"))

    # ACT
    recipe._flux_calibrate_nod_stack(responsePath, stacked, [_frame(1), _frame(3)], [_frame(2), _frame(4)])

    # ASSERT
    header = captured[0]["header"]
    assert header["RA"] == pytest.approx(7.5825, abs=1e-9)
    assert header["OBJECT"] == "SN2026ob"
