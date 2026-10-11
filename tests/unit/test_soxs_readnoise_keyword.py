"""The SOXS read-noise keyword, checked against real SOXS header cards (DY-1457).

SOXS VIS raw frames carry their read noise in `ESO DET OUT1 RON`, and its value
depends on the readout mode. SOXS NIR raw frames carry `ESO DET CHIP RON = 0.0`
and no `ESO DET OUT1 RON`. The `RON` alias used to point at `ESO DET CHIP RON`,
so every VIS readout mode fell back to the 3.8 e- detector-parameter default.

The header cards below are copied from real raw frames:

- VIS mode 1: `SOXS.2026-06-26T04:43:48.245.fits` (STD,FLUX)
- VIS mode 3: `SOXS.2026-07-04T10:00:27.406.fits` (BIAS)
- VIS mode 4: `SOXS.2026-07-03T08:52:04.661.fits` (BIAS)
- VIS mode 1, an older night: `SOXS.2025-05-14T04:02:38.213.fits` (OBJECT)
- NIR: `SOXS.2025-05-14T04:03:11.831.fits` (OBJECT)

The frames are read through the pipeline's own `ImageFileCollection` and
keyword lookup, so a test fails if the alias points at a card the frames lack.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.set_of_files import ImageFileCollection
from soxspipe.commonutils.subtract_sky import subtract_sky
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

LOG = logging.getLogger(__name__)
SETTINGS = {"instrument": "soxs"}
DETECTOR_DEFAULT_RON = 7.3
# THE VIS DETECTOR-PARAMETER DEFAULT EVERY VIS READOUT MODE FELL BACK TO BEFORE DY-1457
VIS_DETECTOR_DEFAULT_RON = 3.8

VIS_MODE_1 = {
    "ESO DET READ CURID": 1,
    "ESO DET READ CURNAME": "Slow High Gain",
    "ESO DET OUT1 RON": 3.3,
}
VIS_MODE_3 = {
    "ESO DET READ CURID": 3,
    "ESO DET READ CURNAME": "Slow Low Gain",
    "ESO DET OUT1 RON": 3.5,
}
VIS_MODE_4 = {
    "ESO DET READ CURID": 4,
    "ESO DET READ CURNAME": "Fast High Gain",
    "ESO DET OUT1 RON": 5.4,
}
VIS_MODE_1_OLDER_NIGHT = {**VIS_MODE_1, "ESO DET OUT1 RON": 3.8}
NIR_UP_THE_RAMP = {
    "ESO DET READ CURID": 3,
    "ESO DET READ CURNAME": "RR_Up-The-Ramp",
    "ESO DET CHIP RON": 0.0,
}


def _kw(name: str) -> str:
    return keyword_lookup(log=LOG, settings=SETTINGS).get(name)


def _write_frame(directory: Path, name: str, arm: str, cards: dict[str, object]) -> None:
    header = fits.Header()
    header["INSTRUME"] = "SOXS"
    header["HIERARCH ESO SEQ ARM"] = arm
    for key, value in cards.items():
        header[f"HIERARCH {key}"] = value
    fits.PrimaryHDU(data=np.zeros((2, 2), dtype=np.float32), header=header).writeto(directory / name)


def _recipe_reading(directory: Path, detectorRon: float = DETECTOR_DEFAULT_RON) -> base_recipe:
    """Return a recipe whose input frames are the FITS files in `directory`.

    The collection loads the same keywords the recipe's verbose frame summary
    loads, so a keyword absent from every frame is a masked column, as in a run.
    """
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = LOG
    recipe.kw = _kw
    recipe.settings = SETTINGS
    recipe.inputFrames = ImageFileCollection(
        location=str(directory),
        keywords=[_kw("SEQ_ARM"), _kw("PRO_CATG"), _kw("RON"), _kw("CHIP_RON"), "file"],
    )
    recipe.detectorParams = {"ron": detectorRon}
    recipe._report_verification_error = lambda **_: None
    return recipe


@pytest.mark.parametrize(
    ("cards", "expectedRon"),
    [(VIS_MODE_1, 3.3), (VIS_MODE_3, 3.5), (VIS_MODE_4, 5.4)],
    ids=["mode-1", "mode-3", "mode-4"],
)
def test_a_vis_recipe_takes_its_read_noise_from_the_readout_mode_header_value(
    tmp_path: Path, cards: dict[str, object], expectedRon: float
) -> None:
    # ARRANGE
    _write_frame(tmp_path, "raw_1.fits", "VIS", cards)
    _write_frame(tmp_path, "raw_2.fits", "VIS", cards)
    recipe = _recipe_reading(tmp_path, detectorRon=VIS_DETECTOR_DEFAULT_RON)

    # ACT
    recipe._verify_single_readnoise()

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(expectedRon)


def test_a_calibration_product_from_another_readout_mode_does_not_trip_the_mix_check(tmp_path: Path) -> None:
    """A product's header read noise describes the frames it was built from, not the frames being reduced.

    Real case: `20260705T023329_VIS_1X1_1_SSOL_MULTPIN_30_0S_SOXS.sof` pairs mode-1 arcs (3.3) with a
    DISP_TAB_VIS built from mode-3 arcs (3.5).
    """
    # ARRANGE
    _write_frame(tmp_path, "raw_arc.fits", "VIS", VIS_MODE_1)
    _write_frame(tmp_path, "disp_tab.fits", "VIS", {**VIS_MODE_3, "ESO PRO CATG": "DISP_TAB_VIS"})
    recipe = _recipe_reading(tmp_path, detectorRon=VIS_DETECTOR_DEFAULT_RON)

    # ACT
    recipe._verify_single_readnoise()

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(3.3)


def test_raw_frames_with_different_read_noise_are_still_rejected(tmp_path: Path) -> None:
    # ARRANGE
    _write_frame(tmp_path, "raw_1.fits", "VIS", VIS_MODE_1)
    _write_frame(tmp_path, "raw_2.fits", "VIS", VIS_MODE_1_OLDER_NIGHT)
    recipe = _recipe_reading(tmp_path, detectorRon=VIS_DETECTOR_DEFAULT_RON)

    # ACT / ASSERT
    with pytest.raises(TypeError, match="mix of readnoise"):
        recipe._verify_single_readnoise()


def test_a_nir_recipe_keeps_the_detector_parameter_read_noise(tmp_path: Path) -> None:
    """NIR frames carry `ESO DET CHIP RON = 0.0` and no `ESO DET OUT1 RON`, so the 7.3 e- default is used."""
    # ARRANGE
    _write_frame(tmp_path, "raw_1.fits", "NIR", NIR_UP_THE_RAMP)
    _write_frame(tmp_path, "raw_2.fits", "NIR", NIR_UP_THE_RAMP)
    recipe = _recipe_reading(tmp_path)

    # ACT
    recipe._verify_single_readnoise()

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(DETECTOR_DEFAULT_RON)


@pytest.mark.parametrize(
    ("arm", "cards", "expectedRon"),
    [("VIS", VIS_MODE_4, 5.4), ("NIR", NIR_UP_THE_RAMP, DETECTOR_DEFAULT_RON)],
    ids=["vis-mode-4", "nir"],
)
def test_the_sky_fit_noise_model_reads_the_frame_header_read_noise(
    tmp_path: Path, arm: str, cards: dict[str, object], expectedRon: float
) -> None:
    # ARRANGE
    _write_frame(tmp_path, "frame.fits", arm, cards)
    header = fits.getheader(tmp_path / "frame.fits")

    # ACT
    ron = subtract_sky._read_noise(header.get(_kw("RON")), DETECTOR_DEFAULT_RON, LOG)

    # ASSERT
    assert ron == pytest.approx(expectedRon)
