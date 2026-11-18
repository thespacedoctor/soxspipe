"""Binned order-centre frames and tables, and flat/table binning ratios below 1, are rejected (DY-1729)."""

from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils.detect_order_edges import detect_order_edges
from soxspipe.recipes.soxs_mflat import soxs_mflat
from soxspipe.recipes.soxs_order_centres import soxs_order_centres
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.unit


def _order_centres_recipe(log: object, *, binx: int, biny: int) -> soxs_order_centres:
    recipe = soxs_order_centres.__new__(soxs_order_centres)
    recipe.log = log
    recipe.kw = lambda keyword: keyword
    recipe.sofName = "single_pinhole_sof"
    recipe.orderFramePath = "/data/night/PINHOLE_VIS.fits"
    recipe.orderFrame = SimpleNamespace(header={"WIN_BINX": binx, "WIN_BINY": biny})
    return recipe


def _order_table(tmp_path: Path, *, binx: int, biny: int) -> Path:
    destination = tmp_path / "orders.fits"
    hdu = fits.PrimaryHDU()
    hdu.header["WIN_BINX"] = binx
    hdu.header["WIN_BINY"] = biny
    hdu.writeto(destination, overwrite=True)
    return destination


def _mflat_recipe(log: object, *, binx: int, biny: int) -> soxs_mflat:
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.kw = lambda keyword: keyword
    flat = CCDData(np.ones((4, 4)), unit=u.electron)
    flat.header["WIN_BINX"] = binx
    flat.header["WIN_BINY"] = biny
    recipe.flats = [flat]
    return recipe


def _edge_detector(log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, **binning: float):
    monkeypatch.setattr("soxspipe.commonutils.toolkit.utility_setup", lambda **_: (str(tmp_path), str(tmp_path)))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    frame = CCDData(np.ones((32, 32)), unit=u.electron, meta=instrument_header(arm="VIS"))
    return detect_order_edges(
        log=log,
        flatFrame=frame,
        orderCentreTable="centres.fits",
        settings=pipeline_settings(tmp_path),
        recipeSettings={"disp-axis-deg": 2, "order-deg": 1},
        startNightDate="2024-01-02",
        **binning,
    )


def test_a_binned_order_centre_frame_is_rejected_naming_the_frame_and_its_binning(log: object) -> None:
    # ARRANGE
    recipe = _order_centres_recipe(log, binx=2, biny=1)

    # ACT
    with pytest.raises(ValueError) as error:
        recipe._read_trace_binning(recipe.kw, "VIS")

    # ASSERT
    message = str(error.value)
    assert "PINHOLE_VIS.fits" in message
    assert "/data/night" not in message
    assert "single_pinhole_sof" in message
    assert "2x1" in message
    assert "must be unbinned" in message


def test_an_unbinned_order_centre_frame_is_read_as_one_by_one(log: object) -> None:
    # ARRANGE
    recipe = _order_centres_recipe(log, binx=1, biny=1)

    # ACT
    binning = recipe._read_trace_binning(recipe.kw, "VIS")

    # ASSERT
    assert binning == (1, 1)


def test_a_binned_order_centre_table_is_rejected_when_the_mflat_reads_it(log: object, tmp_path: Path) -> None:
    # ARRANGE: A 2 X 2 BINNED ORDER-CENTRE TABLE ON AN UNBINNED FLAT, SO THE FLAT/TABLE RATIO IS 0.5
    recipe = _mflat_recipe(log, binx=1, biny=1)
    table = _order_table(tmp_path, binx=2, biny=2)

    # ACT
    with pytest.raises(ValueError) as error:
        recipe._read_binning_ratios(recipe.flats, str(table), recipe.kw)

    # ASSERT
    message = str(error.value)
    assert str(table) in message
    assert "2x2" in message
    assert "must be unbinned" in message


def test_an_unbinned_order_centre_table_gives_the_binned_flat_its_binning_as_the_ratio(
    log: object, tmp_path: Path
) -> None:
    # ARRANGE
    recipe = _mflat_recipe(log, binx=2, biny=2)
    table = _order_table(tmp_path, binx=1, biny=1)

    # ACT
    recipe._read_binning_ratios(recipe.flats, str(table), recipe.kw)

    # ASSERT
    assert (recipe.binRatioX, recipe.binRatioY) == (2, 2)


@pytest.mark.parametrize("binning", [{"binx": 0.5, "biny": 0.5}, {"binx": 1, "biny": 0.5}, {"binx": 0.5, "biny": 1}])
def test_edge_detection_rejects_a_binning_ratio_below_one(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, binning: dict[str, float]
) -> None:
    # ARRANGE, ACT AND ASSERT: THE CONSTRUCTOR REFUSES THE RATIO, AND SHARES THE ORDER-CENTRE RULE WITH THE OTHER CHECKS
    with pytest.raises(ValueError, match="below 1") as error:
        _edge_detector(log, tmp_path, monkeypatch, **binning)

    assert "must be unbinned" in str(error.value)


def test_edge_detection_accepts_a_binned_flat_on_an_unbinned_order_table(
    log: object, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE AND ACT
    detector = _edge_detector(log, tmp_path, monkeypatch, binx=2, biny=2)

    # ASSERT
    assert (detector.axisAbin, detector.axisBbin) == (2, 2)
