"""Contracts for the telluric-line shift QC of the Horne extraction (DY-1969)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.telluric_shift import load_telluric_lines
from tests.unit._horne_nir_scene import build_extractor

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
SETTING = "measure-telluric-shift"
EXTRACTION_RECIPES = ["soxs-stare", "soxs-stare-std", "soxs-nod", "soxs-nod-std", "soxs-offset", "soxs-offset-std"]
SHIFT_NM = -0.05
BAND = "H2O_1130"


def _telluric_orders(log: object, shiftNm: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    """One order that covers the whole H2O_1130 band, with every shipped line of the band shifted by ``shiftNm``."""
    lines = load_telluric_lines(log, {"instrument": "soxs"}, "NIR")
    bandLines = lines.loc[lines["BAND"] == BAND]
    wavelength = np.arange(1105.0, 1165.0, 0.02)
    flux = np.ones_like(wavelength)
    for centre in bandLines["WAVELENGTH"]:
        flux -= 0.3 * np.exp(-0.5 * ((wavelength - centre - shiftNm) / 0.1) ** 2)
    orders = pd.DataFrame({"order": 12, "wavelengthMean": wavelength, "extractedFluxOptimal": flux * 1000.0})
    return orders, bandLines


def _extractor(log: object, tmp_path: Path, recipeSettings: dict | None = None, products: object = False) -> object:
    extractor = build_extractor(log, tmp_path)
    extractor.settings = {"instrument": "soxs"}
    extractor.recipeSettings = {} if recipeSettings is None else recipeSettings
    extractor.notFlattened = ""
    extractor.products = products
    extractor.filenameTemplate = "frame.fits"
    extractor.noddingSequence = ""
    extractor.qcDir = str(tmp_path)
    return extractor


def _qc_names(extractor: object) -> list[str]:
    return list(extractor.qc["qc_name"]) if len(extractor.qc.index) else []


def test_the_extraction_reports_the_telluric_shift_of_the_star(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    extractor = _extractor(log, tmp_path)

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    qc = extractor.qc.set_index("qc_name")
    assert qc.loc["TELL SHIFT", "qc_value"] == pytest.approx(SHIFT_NM, abs=0.01)
    assert qc.loc["TELL NLINES", "qc_value"] >= 5
    assert qc.loc[f"TELL SHIFT {BAND}", "qc_value"] == pytest.approx(SHIFT_NM, abs=0.01)
    assert qc.loc["TELL SHIFT H2O_1450", "qc_value"] == "not measured"


def test_the_measurement_leaves_the_wavelengths_of_the_extraction_unchanged(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    before = orders.copy(deep=True)
    extractor = _extractor(log, tmp_path)

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    pd.testing.assert_frame_equal(orders, before)


def test_the_setting_switches_the_measurement_off(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    extractor = _extractor(log, tmp_path, {SETTING: False})

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    assert _qc_names(extractor) == []


def test_a_setting_that_is_not_true_or_false_stops_the_recipe(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    extractor = _extractor(log, tmp_path, {SETTING: "yes"})

    # ACT AND ASSERT
    with pytest.raises(ValueError, match=SETTING):
        extractor._measure_telluric_shift(orders, "NIR")


def test_the_not_flattened_extraction_measures_nothing(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    extractor = _extractor(log, tmp_path)
    extractor.notFlattened = "_NOTFLAT"

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    assert _qc_names(extractor) == []


def test_a_failed_measurement_is_logged_and_never_stops_the_recipe(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    extractor = _extractor(log, tmp_path)

    def broken(*args: object, **kwargs: object) -> None:
        raise RuntimeError("fit exploded")

    monkeypatch.setattr("soxspipe.commonutils.telluric_shift.measure_telluric_shift", broken)

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    assert _qc_names(extractor) == []
    warnings = [message for level, message in log.messages if level == "warning"]
    assert any("fit exploded" in message for message in warnings)


def test_a_band_with_too_little_signal_is_reported_not_measured_and_the_recipe_goes_on(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    flat = orders.assign(extractedFluxOptimal=1000.0)
    extractor = _extractor(log, tmp_path)

    # ACT
    extractor._measure_telluric_shift(flat, "NIR")

    # ASSERT
    qc = extractor.qc.set_index("qc_name")
    assert qc.loc["TELL NLINES", "qc_value"] == 0
    assert qc.loc["TELL SHIFT", "qc_value"] == "not measured"
    assert not qc.loc["TELL SHIFT", "to_header"]
    assert qc.loc[f"TELL SHIFT {BAND}", "qc_value"] == "not measured"


def test_the_qc_plot_is_registered_as_a_qc_product(log: object, tmp_path: Path) -> None:
    # ARRANGE
    orders, _ = _telluric_orders(log, SHIFT_NM)
    products = pd.DataFrame()
    extractor = _extractor(log, tmp_path, products=products)

    # ACT
    extractor._measure_telluric_shift(orders, "NIR")

    # ASSERT
    row = extractor.products.iloc[0]
    assert row["product_label"] == "TELLURIC_SHIFT_QC_PLOT"
    assert row["file_name"] == "frame_TELLURIC_SHIFT_QC_PLOT.pdf"
    assert row["label"] == "QC"
    assert Path(row["file_path"]).read_bytes().startswith(b"%PDF")


@pytest.mark.parametrize("recipe", EXTRACTION_RECIPES)
def test_each_soxs_extraction_recipe_block_measures_the_telluric_shift_by_default(recipe: str) -> None:
    # ACT
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())

    # ASSERT
    assert settings[recipe][SETTING] is True
