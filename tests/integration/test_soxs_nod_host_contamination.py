"""The nod recipes record host emission at the opposite nod position as an information-only QC (DY-2193).

The detector itself is pinned in ``tests/unit/test_nod_host_contamination.py``. These tests pin the recipe side: the
check reads the flattened A-B difference of each pass, the QC and the warning are recorded before the extractions are
stacked, the warning reaches the merged-spectrum QC plot, and nothing about the check can fail the recipe.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
import yaml

from soxspipe.commonutils.nod_host_contamination import host_feature
from tests.integration.test_soxs_nod_reduction_characterization import (
    SINGLE_PAIR,
    TWO_LOCATIONS,
    _extract_each_cycle,
    _nod_recipe,
    _patch,
)
from tests.integration.test_soxs_nod_shift_and_add import _install_alignment, _real_cycle, _set_aperture, _stacked_run

pytestmark = pytest.mark.integration

QC_NAME = "NOD HOST CONTAM"
SKYLINES_NM = (557.73, 630.03)
HOST_FEATURE = host_feature(order=12, wavelengthMinNm=667.10, wavelengthMaxNm=667.30, significance=11.0)
PACKAGE_PATH = Path(__file__).resolve().parents[2] / "soxspipe"
# THE X-SHOOTER DEFAULTS HAVE NO soxs-nod-std SECTION
NOD_SETTINGS_SECTIONS = [
    ("soxs_default_settings.yaml", "soxs-nod"),
    ("soxs_default_settings.yaml", "soxs-nod-std"),
    ("xsh_default_settings.yaml", "soxs-nod"),
]


def _install_detector(monkeypatch: pytest.MonkeyPatch, features: tuple[host_feature, ...]) -> list[dict[str, Any]]:
    """Replace the detector, the map reader and the skyline catalogue, and return the detector calls."""
    calls: list[dict[str, Any]] = []

    def detect(difference: Any, mapImages: Any, aperturePixels: float, **kwargs: Any) -> tuple[host_feature, ...]:
        calls.append({"difference": difference, "mapImages": mapImages, "aperturePixels": aperturePixels, **kwargs})
        return features

    monkeypatch.setattr(import_module("soxspipe.commonutils.nod_host_contamination"), "find_host_contamination", detect)
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.slit_alignment"),
        "read_two_d_map_images",
        lambda twoDMapPath, binx, biny, kw: ("map", twoDMapPath, binx, biny),
    )
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"),
        "get_skylines_dataframe",
        lambda log, settings, arm: pd.DataFrame({"WAVELENGTH": SKYLINES_NM}),
    )
    return calls


def _qc_rows(recipe: Any) -> pd.DataFrame:
    return recipe.qc.loc[recipe.qc["qc_name"] == QC_NAME]


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def test_a_science_nod_with_a_host_feature_records_the_count_as_a_header_qc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    calls = _install_detector(monkeypatch, (HOST_FEATURE,))

    # ACT
    recipe, _, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT: ONE CHECK OF THE ONE FLATTENED PASS, AND ONE QC ROW WRITTEN TO THE HEADER
    assert len(calls) == 1
    rows = _qc_rows(recipe)
    assert rows["qc_value"].tolist() == [1]
    assert rows["to_header"].tolist() == [True]
    # THE QC IS IN THE TABLE BEFORE THE EXTRACTIONS ARE STACKED, WHICH IS WHEN THE PRODUCT HEADER IS WRITTEN
    assert len(captured["stack_extractions"]) == 1


def test_the_check_reads_the_map_at_the_frame_binning_the_aperture_the_sigma_setting_and_the_skylines(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    calls = _install_detector(monkeypatch, ())
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    _set_aperture(recipe, slitLength=10)
    recipe.recipeSettings = {**recipe.recipeSettings, "host-contamination-sigma": 7}
    _patch(recipe, monkeypatch, tmp_path)
    _install_alignment(monkeypatch, 1.5, -1.5)
    _real_cycle(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT: THE SYNTHETIC FRAMES CARRY NO BINNING, SO THE MAP IS READ 1 BY 1
    call = calls[0]
    assert call["mapImages"] == ("map", recipe.twoDMap, 1, 1)
    assert call["aperturePixels"] == 10
    assert call["sigmaLimit"] == 7
    assert tuple(call["skylineWavelengthsNm"]) == SKYLINES_NM


def test_the_warning_names_the_wavelengths_in_the_log_and_under_the_merged_spectrum_plot(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _install_detector(monkeypatch, (HOST_FEATURE,))

    # ACT
    _, _, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT
    assert any("667.10-667.30 nm" in message for message in _warnings(log))
    notes = captured["plot"][0]["extractionNotes"]
    assert [note for note in notes if "667.10-667.30 nm" in note][0].startswith("WARNING: ")


def test_without_a_host_feature_the_qc_is_zero_and_nothing_is_warned(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _install_detector(monkeypatch, ())

    # ACT
    recipe, _, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT
    assert _qc_rows(recipe)["qc_value"].tolist() == [0]
    assert not any("opposite nod position" in message for message in _warnings(log))
    assert not any("opposite nod position" in note for note in captured["plot"][0]["extractionNotes"])


def test_a_standard_is_checked_on_its_flattened_pass_only_and_records_one_qc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    calls = _install_detector(monkeypatch, (HOST_FEATURE,))

    # ACT
    recipe, constructions, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=True)

    # ASSERT: TWO EXTRACTIONS (FLATTENED AND NOT), ONE CHECK, ONE QC ROW
    assert len(constructions) == 2
    assert len(calls) == 1
    assert _qc_rows(recipe)["qc_value"].tolist() == [1]


def test_the_per_cycle_path_merges_the_features_of_every_cycle_into_one_qc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: EACH CYCLE FINDS THE SAME LINE, SO THE TWO CYCLES COUNT IT ONCE
    calls = _install_detector(monkeypatch, (HOST_FEATURE,))
    recipe, _ = _nod_recipe(log, tmp_path, TWO_LOCATIONS)
    _extract_each_cycle(recipe)
    _set_aperture(recipe)
    _patch(recipe, monkeypatch, tmp_path)
    _real_cycle(recipe, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert len(calls) == 2
    assert _qc_rows(recipe)["qc_value"].tolist() == [1]


def test_a_check_that_cannot_run_is_logged_and_the_reduction_goes_on_without_the_qc(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    _install_detector(monkeypatch, ())

    def broken(*args: Any, **kwargs: Any) -> None:
        raise OSError("no such map")

    monkeypatch.setattr(import_module("soxspipe.commonutils.slit_alignment"), "read_two_d_map_images", broken)

    # ACT
    recipe, _, captured = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)

    # ASSERT
    assert _qc_rows(recipe).empty
    assert any("could not run: no such map" in message for message in _warnings(log))
    assert len(captured["stack_extractions"]) == 1


def test_a_high_count_never_fails_the_recipe(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE: FORTY DISTINCT FEATURES, AND THE RANGES OF THE DEFAULT SETTINGS FILE
    features = tuple(
        host_feature(order=12, wavelengthMinNm=600.0 + step, wavelengthMaxNm=600.1 + step, significance=9.0)
        for step in range(40)
    )
    _install_detector(monkeypatch, features)
    recipe, _, _ = _stacked_run(log, tmp_path, monkeypatch, isStandard=False)
    defaults = yaml.safe_load((PACKAGE_PATH / "soxs_default_settings.yaml").read_text())["soxs-nod"]
    # THE TEMPERATURE QCS OF A SOXS RUN READ HEADER VALUES THE SYNTHETIC RECIPE DOES NOT CARRY
    recipe.detectorTemp = recipe.cptemp = 0.0
    recipe.recipeSettings = {**recipe.recipeSettings, "qc-acceptable-ranges": defaults.get("qc-acceptable-ranges", {})}

    # ACT
    recipe.flag_poor_data()

    # ASSERT
    rows = _qc_rows(recipe)
    assert rows["qc_value"].tolist() == [40]
    assert rows["qc_flag"].tolist() == ["pass"]


@pytest.mark.parametrize(("settingsName", "section"), NOD_SETTINGS_SECTIONS)
def test_the_default_settings_name_the_threshold_and_give_the_qc_no_acceptable_range(
    settingsName: str, section: str
) -> None:
    # ARRANGE
    settings = yaml.safe_load((PACKAGE_PATH / settingsName).read_text())

    # ACT
    recipeSettings = settings[section]

    # ASSERT
    assert recipeSettings["host-contamination-sigma"] == 5
    ranges = {key.lower().replace("-", " ").replace("_", " ") for key in recipeSettings.get("qc-acceptable-ranges", {})}
    assert QC_NAME.lower() not in ranges


def test_an_offset_run_has_no_opposite_nod_position_and_is_not_checked(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    calls = _install_detector(monkeypatch, (HOST_FEATURE,))
    recipe, _ = _nod_recipe(log, tmp_path, SINGLE_PAIR)
    recipe.recipeName = "soxs-offset"

    # ACT
    recipe._check_host_contamination(object())

    # ASSERT
    assert calls == []
    assert recipe.hostFeatures is None
