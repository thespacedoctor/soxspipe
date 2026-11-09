"""Contracts for the ``tune-wavelength-to-skylines`` recipe setting of the Horne extraction (DY-1154)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
import yaml

import soxspipe
from tests.unit._horne_vis_scene import build_extractor, night_header, qc_by_name

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
SETTING = "tune-wavelength-to-skylines"
TUNING_QC = "SKY TUNING"
SWITCHED_OFF_QC = "skipped: switched off"
DEFAULT_SETTINGS_FILES = ["soxs_default_settings.yaml", "xsh_default_settings.yaml"]
# EVERY RECIPE THAT EXTRACTS WITH horne_extraction IS A STARE, NOD OR OFFSET RECIPE (OR ITS -std VARIANT)
EXTRACTION_RECIPE_PREFIXES = ("soxs-stare", "soxs-nod", "soxs-offset")


def _extractor(log: object, tmp_path: Path, recipeSettings: dict) -> object:
    extractor = build_extractor(log, tmp_path, night_header())
    extractor.recipeSettings = recipeSettings
    extractor.calls = []

    def spy(extractedOrdersDF: pd.DataFrame, arm: str, byOrder: bool = True) -> pd.DataFrame:
        extractor.calls.append((arm, byOrder))
        return extractedOrdersDF.assign(tuned=True)

    extractor.tune_wavelength_calibration_to_skylines = spy
    return extractor


def _spectrum() -> pd.DataFrame:
    return pd.DataFrame({"order": [1, 1], "wavelengthMean": [500.0, 500.05]})


@pytest.mark.parametrize("arm", ["VIS", "NIR"])
def test_tuning_switched_off_passes_the_spectrum_through_untuned(log: object, tmp_path: Path, arm: str) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, {SETTING: False})
    spectrum = _spectrum()

    # ACT
    result = extractor._tune_or_skip(spectrum, arm)

    # ASSERT
    assert extractor.calls == []
    assert result is spectrum
    assert result["wavelengthMean"].tolist() == [500.0, 500.05]
    assert extractor.skyShiftsNm is None
    assert qc_by_name(extractor, TUNING_QC)["qc_value"] == SWITCHED_OFF_QC
    printed = [message for level, message in log.messages if level == "print"]
    assert any("Skyline wavelength tuning is switched off in the settings" in m for m in printed)


def test_tuning_switched_on_tunes_once_order_by_order(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, {SETTING: True})

    # ACT
    result = extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == [("VIS", True)]
    assert result["tuned"].all()
    assert qc_by_name(extractor, TUNING_QC)["qc_value"] == "skylines"


def test_a_settings_file_without_the_key_still_tunes(log: object, tmp_path: Path) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, {"use_flat": True})

    # ACT
    result = extractor._tune_or_skip(_spectrum(), "VIS")

    # ASSERT
    assert extractor.calls == [("VIS", True)]
    assert result["tuned"].all()


@pytest.mark.parametrize("value", ["False", 0, None, "no"])
def test_a_value_that_is_not_true_or_false_is_rejected(log: object, tmp_path: Path, value: object) -> None:
    # ARRANGE
    extractor = _extractor(log, tmp_path, {SETTING: value})

    # ACT / ASSERT
    with pytest.raises(ValueError, match=SETTING):
        extractor._tune_or_skip(_spectrum(), "VIS")
    assert extractor.calls == []


def _extraction_blocks() -> list[tuple[str, str]]:
    blocks = []
    for fileName in DEFAULT_SETTINGS_FILES:
        settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())
        blocks.extend((fileName, name) for name in settings if name.startswith(EXTRACTION_RECIPE_PREFIXES))
    return blocks


def test_the_default_settings_files_have_every_extraction_recipe_block() -> None:
    # ACT
    blocks = set(_extraction_blocks())

    # ASSERT
    assert blocks == {
        ("soxs_default_settings.yaml", "soxs-stare"),
        ("soxs_default_settings.yaml", "soxs-stare-std"),
        ("soxs_default_settings.yaml", "soxs-nod"),
        ("soxs_default_settings.yaml", "soxs-nod-std"),
        ("soxs_default_settings.yaml", "soxs-offset"),
        ("soxs_default_settings.yaml", "soxs-offset-std"),
        ("xsh_default_settings.yaml", "soxs-stare"),
        ("xsh_default_settings.yaml", "soxs-nod"),
    }


@pytest.mark.parametrize(("fileName", "recipe"), _extraction_blocks())
def test_each_extraction_recipe_block_tunes_to_skylines_by_default(fileName: str, recipe: str) -> None:
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / fileName).read_text())

    # ASSERT
    assert settings[recipe][SETTING] is True
