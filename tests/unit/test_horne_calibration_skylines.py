"""Contracts for the catalogue lines the Horne-extraction skyline tuning matches against (DY-1645)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

import soxspipe.commonutils as commonutils
import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.horne_extraction import horne_extraction

pytestmark = pytest.mark.unit

DOUBLET = ((785.32587, 8.526), (785.35240, 10.799))
DOUBLET_CENTROID = (DOUBLET[0][0] * DOUBLET[0][1] + DOUBLET[1][0] * DOUBLET[1][1]) / (DOUBLET[0][1] + DOUBLET[1][1])
# THE 0.8 FLUX COMPONENT IS BELOW THE 1.0 CUT OF THE PLOTTING TABLE BUT STILL WEIGHTS THE MERGED CENTROID
FAINT_BLEND_CENTROID = (795.0 * 30.0 + 795.05 * 0.8) / 30.8


def _full_catalogue() -> pd.DataFrame:
    # A DOUBLET WHOSE COMPONENTS CARRY THE ISOLATED FLAG, A CLEAN BRIGHT LINE, AND A FAINT COMPONENT BELOW THE 1.0 CUT
    return pd.DataFrame(
        {
            "WAVELENGTH": [785.32587, 785.35240, 790.0, 795.0, 795.05],
            "FLUX": [8.526, 10.799, 40.0, 30.0, 0.8],
            "ISOLATED": [True, True, True, True, False],
        }
    )


def _tuned_extractor(log: object, arm: str, inst: str = "SOXS") -> horne_extraction:
    extractor = object.__new__(horne_extraction)
    extractor.log = log
    extractor.arm = arm
    extractor.inst = inst
    extractor.settings = {}
    extractor.dispersionMap = "synthetic-dispersion-map.fits"
    catalogue = _full_catalogue()
    extractor.skylinesDF = catalogue.loc[catalogue["FLUX"] > 1.0]
    return extractor


@pytest.fixture
def project_skylines(monkeypatch: pytest.MonkeyPatch) -> None:
    def project(**kwargs: object) -> pd.DataFrame:
        projected = kwargs["orderPixelTable"].copy()
        projected["fit_x"] = projected["wavelength"] - 400.0
        return projected

    monkeypatch.setattr(commonutils, "dispersion_map_to_pixel_arrays", project)


def test_vis_calibration_lines_are_the_merged_isolated_centroids_from_the_unfiltered_catalogue(
    log: object, monkeypatch: pytest.MonkeyPatch, project_skylines: None
) -> None:
    # ARRANGE
    loads: list[dict] = []

    def load_catalogue(log: object, settings: object, arm: str, **kwargs: float) -> pd.DataFrame:
        loads.append({"arm": arm, **kwargs})
        return _full_catalogue()

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", load_catalogue)
    extractor = _tuned_extractor(log, "VIS")

    # ACT
    localSkylines, calibrationSkylines = extractor._get_local_skylines_for_order(785.0, 796.0, 4, "fit_x")

    # ASSERT
    assert loads == [{"arm": "VIS", "minBrightnessVIS": 0.5}]
    assert sorted(calibrationSkylines.tolist()) == pytest.approx(
        [DOUBLET_CENTROID - 400.0, 790.0 - 400.0, FAINT_BLEND_CENTROID - 400.0]
    )
    # THE PLOTTING LINES STILL COME FROM THE CUT TABLE
    assert sorted(localSkylines.tolist()) == pytest.approx([385.32587, 385.3524, 390.0, 395.0])


def test_the_merged_vis_catalogue_is_loaded_once(
    log: object, monkeypatch: pytest.MonkeyPatch, project_skylines: None
) -> None:
    # ARRANGE
    loads: list[int] = []
    monkeypatch.setattr(toolkit, "get_skylines_dataframe", lambda *a, **k: loads.append(1) or _full_catalogue())
    extractor = _tuned_extractor(log, "VIS")

    # ACT
    extractor._get_local_skylines_for_order(785.0, 796.0, 4, "fit_x")
    extractor._get_local_skylines_for_order(785.0, 796.0, 4, "fit_x")

    # ASSERT
    assert loads == [1]


def test_nir_calibration_lines_stay_the_isolated_flagged_lines(
    log: object, monkeypatch: pytest.MonkeyPatch, project_skylines: None
) -> None:
    # ARRANGE
    def refuse_load(*args: object, **kwargs: object) -> pd.DataFrame:
        raise AssertionError("the NIR path must not load the unfiltered catalogue")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", refuse_load)
    extractor = _tuned_extractor(log, "NIR")

    # ACT
    localSkylines, calibrationSkylines = extractor._get_local_skylines_for_order(785.0, 796.0, 4, "fit_x")

    # ASSERT
    assert isinstance(calibrationSkylines, np.ndarray)
    assert sorted(calibrationSkylines.tolist()) == pytest.approx([385.32587, 385.3524, 390.0, 395.0])


def test_xshooter_vis_calibration_lines_stay_the_isolated_flagged_lines(
    log: object, monkeypatch: pytest.MonkeyPatch, project_skylines: None
) -> None:
    # ARRANGE
    def refuse_load(*args: object, **kwargs: object) -> pd.DataFrame:
        raise AssertionError("the X-Shooter VIS path must not load the unfiltered catalogue")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", refuse_load)
    extractor = _tuned_extractor(log, "VIS", inst="XSHOOTER")

    # ACT
    _, calibrationSkylines = extractor._get_local_skylines_for_order(785.0, 796.0, 4, "fit_x")

    # ASSERT
    assert sorted(calibrationSkylines.tolist()) == pytest.approx([385.32587, 385.3524, 390.0, 395.0])


def test_global_tuning_keeps_the_isolated_flagged_lines_in_soxs_vis(
    log: object, monkeypatch: pytest.MonkeyPatch, project_skylines: None
) -> None:
    # ARRANGE
    def refuse_load(*args: object, **kwargs: object) -> pd.DataFrame:
        raise AssertionError("the all-orders lookup must not load the unfiltered catalogue")

    monkeypatch.setattr(toolkit, "get_skylines_dataframe", refuse_load)
    extractor = _tuned_extractor(log, "VIS")

    # ACT
    _, calibrationSkylines = extractor._get_local_skylines_for_order(785.0, 796.0, "all", "wavelength")

    # ASSERT
    assert sorted(calibrationSkylines.tolist()) == pytest.approx([785.32587, 785.3524, 790.0, 795.0])
