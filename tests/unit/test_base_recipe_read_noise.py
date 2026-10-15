"""A UVB/VIS recipe takes its read noise from the master bias its SOF names (DY-1458)."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.commonutils.set_of_files import ImageFileCollection
from soxspipe.recipes import base_recipe

pytestmark = pytest.mark.unit

baseRecipeModule = importlib.import_module("soxspipe.recipes.base_recipe")

SETTINGS = {"instrument": "soxs"}
NOMINAL_RON = 3.3
DETECTOR_DEFAULT_RON = 3.8
MEASURED_RON = 4.0
N_COMBINED = 4


def _write_frame(directory: Path, name: str, arm: str, cards: dict[str, object]) -> None:
    header = fits.Header()
    header["INSTRUME"] = "SOXS"
    header["HIERARCH ESO SEQ ARM"] = arm
    for key, value in cards.items():
        header[f"HIERARCH {key}"] = value
    fits.PrimaryHDU(data=np.zeros((2, 2), dtype=np.float32), header=header).writeto(directory / name)


def _write_scene(directory: Path, arm: str, masterBiasCards: dict[str, object] | None) -> None:
    """Two raw frames at the nominal read noise, and a master bias unless `masterBiasCards` is None."""
    for name in ("raw_1.fits", "raw_2.fits"):
        _write_frame(directory, name, arm, {"ESO DET OUT1 RON": NOMINAL_RON})
    if masterBiasCards is not None:
        _write_frame(directory, "master_bias.fits", arm, {"ESO PRO CATG": f"MASTER_BIAS_{arm}", **masterBiasCards})


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
    recipe.detectorParams = {"ron": DETECTOR_DEFAULT_RON}
    recipe._report_verification_error = lambda **_: None
    return recipe


def _verify_read_noise(recipe: base_recipe) -> None:
    recipe._verify_single_readnoise()
    recipe._apply_measured_read_noise()


MEASURED_CARDS = {"ESO QC RAW RON": MEASURED_RON, "ESO PRO DATANCOM": N_COMBINED}


def test_a_vis_recipe_takes_the_measured_read_noise_from_its_master_bias(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "VIS", MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(MEASURED_RON)
    assert recipe.detectorParams["ron"].unit == "electron"
    assert recipe.readNoise.nCombined == N_COMBINED
    assert recipe.readNoise.skyFit == pytest.approx(4.472136, rel=1e-6)


def test_a_vis_recipe_without_a_master_bias_keeps_the_nominal_read_noise(tmp_path: Path, log: Any) -> None:
    # ARRANGE
    _write_scene(tmp_path, "VIS", None)
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(NOMINAL_RON)
    assert recipe.readNoise is None


def test_a_nir_recipe_never_looks_for_a_master_bias(tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: A STRAY MASTER BIAS WITH A MEASURED VALUE MUST BE IGNORED
    _write_scene(tmp_path, "NIR", MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "NIR")

    def refuse(*_: object, **__: object) -> None:
        raise AssertionError("a NIR recipe must not read a master bias header")

    monkeypatch.setattr(fits, "getheader", refuse)

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(NOMINAL_RON)
    assert recipe.readNoise is None


def test_an_old_master_bias_falls_back_to_the_nominal_read_noise_with_a_warning(tmp_path: Path, log: Any) -> None:
    # ARRANGE: NO DATANCOM
    _write_scene(tmp_path, "VIS", {"ESO QC RAW RON": MEASURED_RON})
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(NOMINAL_RON)
    assert [level for level, _ in log.messages].count("warning") >= 1


def test_several_master_biases_use_the_last_by_name_and_warn(tmp_path: Path, log: Any) -> None:
    # ARRANGE: THE LATER NAME IS WRITTEN FIRST, SO DIRECTORY ORDER CANNOT DECIDE
    _write_scene(tmp_path, "VIS", None)
    _write_frame(
        tmp_path,
        "master_bias_b.fits",
        "VIS",
        {"ESO PRO CATG": "MASTER_BIAS_VIS", "ESO QC RAW RON": 3.4, "ESO PRO DATANCOM": N_COMBINED},
    )
    _write_frame(tmp_path, "master_bias_a.fits", "VIS", {"ESO PRO CATG": "MASTER_BIAS_VIS", **MEASURED_CARDS})
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(3.4)
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("2 master biases" in message and "master_bias_b.fits" in message for message in warnings)


def test_the_detector_default_is_the_last_resort_when_no_raw_frame_carries_a_read_noise(
    tmp_path: Path, log: Any
) -> None:
    # ARRANGE
    for name in ("raw_1.fits", "raw_2.fits"):
        _write_frame(tmp_path, name, "VIS", {})
    _write_frame(tmp_path, "master_bias.fits", "VIS", {"ESO PRO CATG": "MASTER_BIAS_VIS"})
    recipe = _recipe(tmp_path, log, "VIS")

    # ACT
    _verify_read_noise(recipe)

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(DETECTOR_DEFAULT_RON)
    assert recipe.readNoise.source == "detector default"


def test_the_basic_verification_applies_the_measured_read_noise(
    tmp_path: Path, log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    """`_verify_input_frames_basics` sets the measured noise, so it is in place before frames are prepared."""
    # ARRANGE
    _write_scene(tmp_path, "VIS", MEASURED_CARDS)
    recipe = _recipe(tmp_path, log, "VIS")
    recipe.get_recipe_settings = lambda: {}
    recipe._verify_single_instrument_and_arm = lambda _: None
    recipe._verify_single_binning = lambda _: None
    recipe._verify_single_read_speed = lambda _: None
    recipe._verify_single_gain = lambda: None
    recipe._verify_single_slit_width = lambda: None
    recipe._collect_frame_classifications = lambda: ((), (), ())
    recipe.detectorParams = None
    monkeypatch.setattr(
        baseRecipeModule,
        "detector_lookup",
        lambda **_: type(
            "Lookup", (), {"get": staticmethod(lambda arm: {"ron": DETECTOR_DEFAULT_RON, "dispersion-axis": "x"})}
        )(),
    )

    # ACT
    recipe._verify_input_frames_basics()

    # ASSERT
    assert recipe.detectorParams["ron"].value == pytest.approx(MEASURED_RON)
