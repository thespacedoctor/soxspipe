"""The response and the flux calibration correct for the detector binning, so binned science calibrates correctly.

DY-1622: a VIS science frame binned twice along the dispersion axis collects twice the flat-fielded counts per pixel,
and it calibrated to twice the true flux against a response built from an unbinned standard.
"""

from __future__ import annotations

import importlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.table import Table

from soxspipe.commonutils.flux_calibration import flux_calibration
from soxspipe.commonutils.response_basis import counts_per_unbinned_pixel, dispersion_binning
from soxspipe.commonutils.response_function import response_function

pytestmark = pytest.mark.integration
responseModule = importlib.import_module("soxspipe.commonutils.response_function")
toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")

# THE STAR: CATALOGUE FLUX (ERG/CM2/S/ANGSTROM) AND THE FLAT-FIELDED COUNTS IT GIVES IN ONE UNBINNED PIXEL
CATALOGUE_FLUX = 2.0e-16
COUNTS_PER_UNBINNED_PIXEL = 200.0
EXPOSURE_TIME = 10.0
WAVELENGTHS = np.linspace(400.0, 800.0, 4000)


def _keyword(key: str) -> str:
    """Return a FITS keyword name for a keyword alias, standing in for the keyword lookup."""
    return f"ESO {key.replace('_', ' ')}"


def _spectrum(dispersionBinning: int) -> dict[str, np.ndarray]:
    """Return the WAVE and FLUX_COUNTS columns the star gives with ``dispersionBinning`` pixels binned together."""
    return {
        "WAVE": WAVELENGTHS,
        "FLUX_COUNTS": np.full(WAVELENGTHS.size, COUNTS_PER_UNBINNED_PIXEL * dispersionBinning),
    }


def _header(dispersionBinning: int) -> fits.Header:
    """Return a VIS standard-star header with ``dispersionBinning`` as its ESO DET BINY card."""
    return fits.Header(
        {
            _keyword("SEQ_ARM"): "VIS",
            _keyword("DATE_OBS"): "2024-01-02T03:04:05",
            "DATE-OBS": "2024-01-02T03:04:05",
            _keyword("EXPTIME"): EXPOSURE_TIME,
            _keyword("OBJECT"): "LTT7987",
            _keyword("WIN_BINY"): dispersionBinning,
            "HIERARCH ESO TEL AIRM START": 1.0,
            "HIERARCH ESO TEL AIRM END": 1.0,
        }
    )


def _build_response_product(tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Run the real VIS response_function on an unbinned synthetic standard and return the RESP product path."""
    import soxspipe.commonutils as commonutils

    extractedPath = tmp_path / "standard.fits"
    unflattenedPath = tmp_path / "standard-unflattened.fits"
    header = _header(dispersionBinning=1)
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(Table(_spectrum(1)))]).writeto(extractedPath)
    unflattened = Table({"WAVE": WAVELENGTHS, "FLUX_DENSITY_COUNTS": np.full(WAVELENGTHS.size, 5000.0)})
    fits.HDUList([fits.PrimaryHDU(header=header), fits.BinTableHDU(unflattened)]).writeto(unflattenedPath)
    Table({"WAVE": WAVELENGTHS, "LTT7987": np.full(WAVELENGTHS.size, CATALOGUE_FLUX)}).write(
        tmp_path / "standards.fits"
    )
    (tmp_path / "qc").mkdir()
    (tmp_path / "products").mkdir()

    monkeypatch.setattr(commonutils, "keyword_lookup", lambda **_: SimpleNamespace(get=_keyword))
    monkeypatch.setattr(
        commonutils,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"flux-standards": "standards.fits", "extinction": "ext.fits"}),
    )
    monkeypatch.setattr(toolkitModule, "get_calibrations_path", lambda **_: str(tmp_path))
    monkeypatch.setattr(toolkitModule, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))

    # NO ATMOSPHERE: THE EXTINCTION CORRECTION IS 1 AT EVERY WAVELENGTH
    def no_extinction(wavelengths: object, path: str, airmass: float) -> np.ndarray:
        return np.ones(len(wavelengths))

    monkeypatch.setattr(responseModule, "extinction_correction_factor", no_extinction)
    monkeypatch.setattr(toolkitModule, "extinction_correction_factor", no_extinction)
    response = response_function(
        log=log,
        stdExtractionPath=str(extractedPath),
        recipeName="soxs-response",
        settings={"instrument": "soxs", "soxs-response": {"vis": {}}},
        sofName="SYNTHETIC",
        startNightDate="2024-01-02",
        stdNotFlatExtractionPath=str(unflattenedPath),
        orderJoins={"1-3": 600.0},
    )
    response.get()
    return tmp_path / "products" / "SYNTHETIC_RESP.fits"


def _calibrate(tmp_path: Path, log: object, responsePath: Path, dispersionBinning: int) -> np.ndarray:
    """Flux calibrate the star observed with ``dispersionBinning`` and return the calibrated flux."""
    calibration = object.__new__(flux_calibration)
    calibration.log = log
    calibration.settings = {}
    calibration.responseFunction = str(responsePath)
    calibration.extractedSpectrum = pd.DataFrame(_spectrum(dispersionBinning))
    calibration.exptime = EXPOSURE_TIME
    calibration.airmass = 1.0
    calibration.extinctionPath = "ext.fits"
    calibration.arm = "VIS"
    calibration.debug = False
    calibration.header = _header(dispersionBinning)
    calibration.recipeName = "soxs-stare"
    calibration.sofName = f"synthetic_bin{dispersionBinning}"
    calibration.productDir = str(tmp_path)
    calibration.products = pd.DataFrame()
    calibration.kw = _keyword
    outputPath, _ = calibration.calibrate()
    return Table.read(outputPath)["FLUX_CALIBRATED"].data


def test_science_binned_along_the_dispersion_axis_calibrates_like_unbinned_science(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A 2x2-binned and a 1x1 science spectrum of the standard both calibrate to its catalogue flux (DY-1622)."""
    # ARRANGE
    responsePath = _build_response_product(tmp_path, log, monkeypatch)

    # ACT
    unbinnedFlux = _calibrate(tmp_path, log, responsePath, dispersionBinning=1)
    binnedFlux = _calibrate(tmp_path, log, responsePath, dispersionBinning=2)

    # ASSERT
    np.testing.assert_allclose(unbinnedFlux, CATALOGUE_FLUX, rtol=1e-6)
    np.testing.assert_allclose(binnedFlux, CATALOGUE_FLUX, rtol=1e-6)


def test_response_product_states_its_per_unbinned_pixel_basis(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The RESP product header says its curve scales flat-fielded counts per unbinned pixel (DY-1622)."""
    # ARRANGE AND ACT
    responsePath = _build_response_product(tmp_path, log, monkeypatch)

    # ASSERT
    with fits.open(responsePath) as hdul:
        assert hdul[0].header["ESO PRO RESP BASIS"] == "PER_UNBINNED_PIXEL"


@pytest.mark.parametrize("hasCurve", [True, False], ids=["curve-without-basis", "coefficients-only"])
def test_flux_calibration_refuses_a_response_product_without_the_basis_card(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
    hasCurve: bool,
) -> None:
    """A RESP product from before DY-1622 ignores the binning, so flux calibration refuses it and asks for a rebuild."""
    # ARRANGE
    responsePath = tmp_path / "old_RESP.fits"
    hdus = [fits.PrimaryHDU(), fits.table_to_hdu(Table({"polyOrder": [0], "c0": [2.0]}))]
    if hasCurve:
        curveHDU = fits.table_to_hdu(Table({"WAVE": [400.0, 800.0], "RESPONSE": [2.0, 2.0]}))
        curveHDU.name = "RESPONSE"
        hdus.append(curveHDU)
    fits.HDUList(hdus).writeto(responsePath)
    monkeypatch.setattr(toolkitModule, "extinction_correction_factor", lambda w, p, a: np.ones(len(w)))

    # ACT / ASSERT
    with pytest.raises(ValueError, match="rebuild the response"):
        _calibrate(tmp_path, log, responsePath, dispersionBinning=1)


@pytest.mark.parametrize(
    ("arm", "cards", "expected"),
    [
        ("NIR", {}, 1),
        ("VIS", {"ESO WIN BINY": 2}, 2),
        ("VIS", {"ESO WIN BINY": "2"}, 2),
        ("VIS", {"ESO WIN BINY": 2.0}, 2),
    ],
    ids=["nir", "vis-int", "vis-string", "vis-integral-float"],
)
def test_dispersion_binning_reads_biny_on_vis_and_is_one_on_nir(
    arm: str, cards: dict[str, object], expected: int
) -> None:
    """VIS reads the binning along the dispersion axis from ESO DET BINY; NIR is never binned (DY-1622)."""
    # ARRANGE
    header = fits.Header(cards)

    # ACT
    binning = dispersion_binning(header, arm, _keyword)

    # ASSERT
    assert binning == expected


@pytest.mark.parametrize(
    "cards",
    [{}, {"ESO WIN BINY": 0}, {"ESO WIN BINY": 1.5}, {"ESO WIN BINY": "two"}],
    ids=["missing", "zero", "fractional", "malformed"],
)
def test_dispersion_binning_refuses_a_vis_header_without_a_usable_biny(cards: dict[str, object]) -> None:
    """A VIS header without a whole-number binning of 1 or more cannot be flux calibrated (DY-1622)."""
    # ARRANGE
    header = fits.Header(cards)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="ESO WIN BINY"):
        dispersion_binning(header, "VIS", _keyword)


def test_dispersion_binning_refuses_a_vis_spectrum_with_no_header() -> None:
    """flux_calibration defaults its header to None, which is a missing card, not a TypeError (DY-1622)."""
    # ACT / ASSERT
    with pytest.raises(ValueError, match="ESO WIN BINY"):
        dispersion_binning(None, "VIS", _keyword)


def test_counts_per_unbinned_pixel_divides_the_flat_fielded_counts_by_the_binning() -> None:
    """Counts in a pixel binned twice along the dispersion axis are halved (DY-1622)."""
    # ARRANGE
    spectrum = pd.DataFrame({"FLUX_COUNTS": [400.0, 800.0]})

    # ACT
    counts = counts_per_unbinned_pixel(spectrum, 2)

    # ASSERT
    np.testing.assert_array_equal(counts, [200.0, 400.0])
