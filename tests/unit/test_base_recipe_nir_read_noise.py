"""A NIR recipe takes its read noise from the master dark its SOF names (DY-1581).

NIR raw frames carry no usable read-noise card (`ESO DET CHIP RON` is 0.0 and there is no `ESO DET OUT1 RON`), so the
nominal value is the detector-parameter default. The master dark records `ESO QC RAW RON` and `ESO PRO DATANCOM`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.set_of_files import ImageFileCollection
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

SETTINGS = {"instrument": "soxs"}
YAML_RON = 7.3
MEASURED_RON = 7.6
N_COMBINED = 4
MEASURED_CARDS = {"ESO QC RAW RON": MEASURED_RON, "ESO PRO DATANCOM": N_COMBINED}
VIS_BIAS_CARDS = {"ESO QC RAW RON": 4.0, "ESO PRO DATANCOM": 4}


def _write_frame(directory: Path, name: str, arm: str, cards: dict[str, object]) -> None:
    header = fits.Header()
    header["INSTRUME"] = "SOXS"
    header["HIERARCH ESO SEQ ARM"] = arm
    for key, value in cards.items():
        header[f"HIERARCH {key}"] = value
    fits.PrimaryHDU(data=np.zeros((2, 2), dtype=np.float32), header=header).writeto(directory / name)


def _write_scene(directory: Path, arm: str, rawCards: dict[str, object], masterCards: dict[str, object] | None) -> None:
    """Two raw frames, and a master dark unless `masterCards` is None."""
    for name in ("raw_1.fits", "raw_2.fits"):
        _write_frame(directory, name, arm, rawCards)
    if masterCards is not None:
        _write_frame(directory, "master_dark.fits", arm, {"ESO PRO CATG": f"MASTER_DARK_{arm}", **masterCards})


def _recipe(directory: Path, log: Any, arm: str) -> base_recipe:
    kw = keyword_lookup(log=log, settings=SETTINGS).get
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = kw
    recipe.arm = arm
    recipe.settings = SETTINGS
    recipe.inputFrames = ImageFileCollection(
        location=str(directory),
        keywords=[kw("SEQ_ARM"), kw("PRO_CATG"), kw("RON"), kw("CHIP_RON"), "file"],
    )
    recipe.detectorParams = {"ron": YAML_RON}
    recipe._report_verification_error = lambda **_: None
    return recipe


def _verify_read_noise(recipe: base_recipe) -> None:
    recipe._verify_single_readnoise()
    recipe._apply_measured_read_noise()


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _infos(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "info"]


NIR_RAW = {"ESO DET CHIP RON": 0.0}


def test_a_nir_recipe_takes_the_measured_read_noise_from_its_master_dark(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "NIR", NIR_RAW, MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT: THE PREP ERRS GETS NO 1/N TERM; THE SKY FIT GETS 7.6 * SQRT(1 + 1/4) = 8.4971
    assert recipe.detectorParams["ron"].value == pytest.approx(MEASURED_RON)
    assert recipe.detectorParams["ron"].unit == "electron"
    assert recipe.readNoise.nCombined == N_COMBINED
    assert recipe.readNoise.skyFit == pytest.approx(8.497058, rel=1e-6)
    assert recipe.readNoise.source == "master dark"


def test_a_nir_recipe_without_a_master_dark_keeps_the_yaml_read_noise(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "NIR", NIR_RAW, None)
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(YAML_RON)
    assert recipe.readNoise is None


@pytest.mark.parametrize(
    "masterCards",
    [
        {"ESO QC RAW RON": MEASURED_RON},
        {"ESO PRO DATANCOM": 1},
        {},
    ],
    ids=["old product without DATANCOM", "single raw frame without a RAW RON", "no cards at all"],
)
def test_an_unusable_master_dark_falls_back_to_the_yaml_value_with_a_warning(
    tmp_path: Path, log: Any, masterCards: dict[str, object]
) -> None:
    # ARRANGE
    _write_scene(tmp_path, "NIR", NIR_RAW, masterCards)
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(YAML_RON)
    assert recipe.readNoise.source == "detector default"
    assert any("master dark" in message for message in _warnings(log))


def test_a_single_raw_master_dark_gives_the_sky_fit_its_full_master_noise_term(tmp_path: Path, log: Any) -> None:
    """A master dark from one raw frame carries the full read noise, so the sky fit doubles the variance."""
    # ARRANGE
    _write_scene(tmp_path, "NIR", NIR_RAW, {"ESO PRO DATANCOM": 1})
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT: 7.3 * SQRT(2) = 10.3238
    assert recipe.readNoise.nCombined == 1
    assert recipe.readNoise.skyFit == pytest.approx(10.323759, rel=1e-6)


def test_an_unreadable_master_dark_header_falls_back_to_the_yaml_value_with_a_warning(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _write_scene(tmp_path, "NIR", NIR_RAW, MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "NIR")

    def unreadable(*_: object, **__: object) -> None:
        raise OSError("truncated file")

    monkeypatch.setattr(fits, "getheader", unreadable)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(YAML_RON)
    assert any("master dark" in message and "could not be read" in message for message in _warnings(log))


def test_a_measured_nir_value_over_fifteen_percent_from_the_yaml_value_is_logged_as_info_and_used(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE: 15.0 IS 105% ABOVE 7.3
    _write_scene(tmp_path, "NIR", NIR_RAW, {"ESO QC RAW RON": 15.0, "ESO PRO DATANCOM": N_COMBINED})
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(15.0)
    assert any("nominal 7.30" in message for message in _infos(log))


def test_a_zero_read_noise_card_on_the_nir_raw_frames_is_not_taken_as_the_read_noise(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "NIR", {"ESO DET OUT1 RON": 0.0}, MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    recipe._verify_single_readnoise()

    # ASSERT
    assert recipe.rawHeaderRon is None
    assert recipe.detectorParams["ron"].value == pytest.approx(YAML_RON)


def test_a_vis_recipe_with_a_master_bias_and_a_master_dark_uses_the_master_bias(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE DARK HOLDS A DIFFERENT VALUE, SO A WRONG CHOICE SHOWS
    _write_scene(tmp_path, "VIS", {"ESO DET OUT1 RON": 3.3}, {"ESO QC RAW RON": 9.0, "ESO PRO DATANCOM": 4})
    _write_frame(tmp_path, "master_bias.fits", "VIS", {"ESO PRO CATG": "MASTER_BIAS_VIS", **VIS_BIAS_CARDS})
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(4.0)
    assert recipe.readNoise.source == "master bias"


def test_a_vis_recipe_with_only_a_master_dark_keeps_the_nominal_read_noise(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "VIS", {"ESO DET OUT1 RON": 3.3}, {"ESO QC RAW RON": 9.0, "ESO PRO DATANCOM": 4})
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(3.3)
    assert recipe.readNoise is None


def test_several_nir_master_darks_use_the_last_by_name_and_warn(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE LATER NAME IS WRITTEN FIRST, SO DIRECTORY ORDER CANNOT DECIDE
    _write_scene(tmp_path, "NIR", NIR_RAW, None)
    _write_frame(
        tmp_path,
        "master_dark_b.fits",
        "NIR",
        {"ESO PRO CATG": "MASTER_DARK_NIR", "ESO QC RAW RON": 7.0, "ESO PRO DATANCOM": N_COMBINED},
    )
    _write_frame(tmp_path, "master_dark_a.fits", "NIR", {"ESO PRO CATG": "MASTER_DARK_NIR", **MEASURED_CARDS})
    recipe = _recipe(tmp_path, log, "NIR")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(7.0)
    assert any("2 master darks" in message and "master_dark_b.fits" in message for message in _warnings(log))
