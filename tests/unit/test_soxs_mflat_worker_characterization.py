"""Characterization of `soxs_mflat.normalise_flats` and `soxs_mflat.mask_low_sens_pixels`.

These tests pin what the recipe does today, before commit 2 replaces its
inline `datetime.utcnow().strftime(...)` timestamps and inline
`pd.concat([..., pd.DataFrame([{row}])])` QC-row construction with
`toolkit.utcnow_string` and `toolkit.append_qc`.

`tests/unit/test_mflat_helpers.py` already covers `normalise_flats` and
`mask_low_sens_pixels` at a shallow level, on constant data that makes
normalisation trivially `1.0` and never exercises the second pass, the
`axisA == "y"` branch, or the binning fallbacks. These tests use seeded
non-constant data instead, and start `recipe.qc` from the real declared
table (`base_recipe._empty_qc_and_product_tables`), not the bare
`pd.DataFrame()` `test_mflat_helpers.py` uses -- the declared table's columns
are all empty-list float64 columns, and appending the *first* row into it
coerces that row's `to_header=True` into `1.0`, not the Python `True` every
later row keeps. That is pinned below as found.
"""

from __future__ import annotations

import re
from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.recipes.base_recipe import base_recipe
from soxspipe.recipes.soxs_mflat import soxs_mflat

pytestmark = pytest.mark.unit

# IMPORTED VIA `import_module`, NOT A DOTTED `import`: THE PACKAGE `__init__`
# REBINDS THE ATTRIBUTE `soxspipe.recipes.soxs_mflat` TO THE CLASS ITSELF.
mflatModule = import_module("soxspipe.recipes.soxs_mflat")

# THE FORMAT EVERY `datetime.utcnow().strftime(...)` CALL RENDERS.
TIMESTAMP_PATTERN = re.compile(r"\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}")


def _frame(data: np.ndarray, *, binx: int = 1, biny: int = 1, withHeaderBinning: bool = True) -> CCDData:
    """Return a mutable CCD frame with deterministic mask and uncertainty."""
    frame = CCDData(data, unit="electron")
    frame.mask = np.zeros(data.shape, dtype=bool)
    frame.uncertainty = StdDevUncertainty(np.ones(data.shape) * 0.5)
    if withHeaderBinning:
        frame.header["WIN_BINX"] = binx
        frame.header["WIN_BINY"] = biny
    return frame


def _recipe(log: object, *, axisA: str = "x", axisB: str = "y", arm: str = "VIS") -> soxs_mflat:
    """Return the small recipe state required by the normalisation and masking seams.

    `qc` and `products` start from the real declared tables, not a bare
    `pd.DataFrame()`, so the column order and the `to_header` dtype quirk
    pinned here match what a real recipe carries.
    """
    recipe = soxs_mflat.__new__(soxs_mflat)
    recipe.log = log
    recipe.kw = lambda key: key
    recipe.arm = arm
    recipe.axisA = axisA
    recipe.axisB = axisB
    recipe.binx = 1
    recipe.biny = 1
    recipe.recipeName = "soxs-mflat"
    recipe.recipeSettings = {
        "centre-order-window": 2,
        "low-sensitivity-clipping-sigma": 3,
        "scale-d2-to-qth": True,
    }
    recipe.qc, recipe.products = base_recipe._empty_qc_and_product_tables(None, pd)
    recipe.dateObs = "2024-01-01T00:00:00"
    recipe.debug = False
    return recipe


def _order_table(
    tmp_path: Path,
    *,
    binx: int | None = None,
    biny: int | None = None,
    name: str = "orders.fits",
) -> Path:
    """Write a minimal order-table FITS file, optionally carrying its own binning."""
    destination = tmp_path / name
    hdu = fits.PrimaryHDU()
    if binx is not None:
        hdu.header["WIN_BINX"] = binx
    if biny is not None:
        hdu.header["WIN_BINY"] = biny
    hdu.writeto(destination, overwrite=True)
    return destination


def _stub_unpack_order_table(monkeypatch: pytest.MonkeyPatch, pixels: pd.DataFrame) -> None:
    monkeypatch.setattr(mflatModule, "unpack_order_table", lambda **kwargs: (None, pixels, None))


def _centre_pixels(nrows: int, *, xCentre: int = 6) -> pd.DataFrame:
    """Return a vertical order-centre trace covering every row of a frame."""
    return pd.DataFrame({"xcoord_centre": [xCentre] * nrows, "ycoord": list(range(nrows))})


# ---------------------------------------------------------------------------
# normalise_flats
# ---------------------------------------------------------------------------


def test_first_pass_normalises_each_frame_to_its_sigma_clipped_mean(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """On seeded non-constant data, each frame is divided by its own clipped mean flux."""
    # ARRANGE
    recipe = _recipe(log)
    orderTable = _order_table(tmp_path)
    _stub_unpack_order_table(monkeypatch, _centre_pixels(12))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)

    frame1 = _frame(np.random.default_rng(101).normal(loc=1000.0, scale=25.0, size=(12, 12)))
    frame2 = _frame(np.random.default_rng(102).normal(loc=2000.0, scale=40.0, size=(12, 12)))

    # ACT
    normalised = recipe.normalise_flats([frame1, frame2], str(orderTable))

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [0.9754093904226082, 0.9444510799461148, 1.010073854688564],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].uncertainty.array[0, :3],
        [0.0004975328661615877] * 3,
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[1].data[0, :3],
        [1.014338558911565, 1.0451659114864156, 1.0209466647338135],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[1].uncertainty.array[0, :3],
        [0.00025045035670796795] * 3,
        rtol=1e-12,
        atol=0,
    )


def test_first_pass_appends_ordexp_rows_onto_the_declared_qc_table(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """ORDEXP10, ORDEXP50 and ORDEXP90 are appended in that order, with distinct value formats.

    The first row ever appended into the declared-empty QC table lands its
    `to_header=True` as the float `1.0`, not the Python `True` the next two
    rows keep -- the declared table's `to_header` column starts as an
    all-empty `float64` column, and `pd.concat` coerces the first bool it
    receives into that dtype before the column widens to `object`. Pinned as
    found, not as intended.
    """
    # ARRANGE
    recipe = _recipe(log)
    orderTable = _order_table(tmp_path)
    _stub_unpack_order_table(monkeypatch, _centre_pixels(12))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(201).normal(loc=1000.0, scale=25.0, size=(12, 12)))

    # ACT
    recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    assert list(recipe.qc.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_unit",
        "qc_order",
        "qc_comment",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]
    assert list(recipe.qc["qc_name"]) == ["ORDEXP10", "ORDEXP50", "ORDEXP90"]
    assert recipe.qc["qc_value"].tolist() == [
        "958.67519208595945201523136",
        "997.655",
        "1027.085",
    ]
    assert list(recipe.qc["qc_unit"]) == ["electrons"] * 3
    assert list(recipe.qc["qc_comment"]) == [
        "[e-] 10th percentile inter-order flux",
        "[e-] 50th percentile inter-order flux",
        "[e-] 90th percentile inter-order flux",
    ]
    assert list(recipe.qc["soxspipe_recipe"]) == ["soxs-mflat"] * 3
    assert list(recipe.qc["obs_date_utc"]) == [recipe.dateObs] * 3
    for value in recipe.qc["reduction_date_utc"]:
        assert TIMESTAMP_PATTERN.fullmatch(value)
    # THE THREE ORDEXP ROWS ARE STAMPED WITH ONE REDUCTION TIME
    assert len(set(recipe.qc["reduction_date_utc"])) == 1
    toHeaderValues = recipe.qc["to_header"].tolist()
    assert toHeaderValues[0] == 1.0
    assert isinstance(toHeaderValues[0], float)
    assert toHeaderValues[1] is True
    assert toHeaderValues[2] is True


def test_second_pass_divides_by_the_first_pass_master_flat_and_writes_no_qc(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A `firstPassMasterFlat` argument selects the second-pass branch, which writes no QC rows."""
    # ARRANGE
    recipe = _recipe(log)
    orderTable = _order_table(tmp_path)
    _stub_unpack_order_table(monkeypatch, _centre_pixels(12))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)

    frame = _frame(np.random.default_rng(201).normal(loc=1000.0, scale=25.0, size=(12, 12)))
    masterFlat = _frame(1.0 + np.random.default_rng(999).normal(loc=0.0, scale=0.01, size=(12, 12)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable), firstPassMasterFlat=masterFlat)

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [1.0556058037020695, 0.9917432187486499, 0.9548850506963418],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].uncertainty.array[0, :3],
        [0.0005035211340627604] * 3,
        rtol=1e-12,
        atol=0,
    )
    assert len(recipe.qc) == 0


def test_the_axis_a_y_branch_masks_a_horizontal_order_centre_band(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """When `axisA == "y"`, the order-centre mask is applied along the first axis."""
    # ARRANGE
    recipe = _recipe(log, axisA="y", axisB="x")
    orderTable = _order_table(tmp_path)
    pixels = pd.DataFrame({"ycoord_centre": [6] * 12, "xcoord": list(range(12))})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(301).normal(loc=500.0, scale=10.0, size=(12, 12)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [1.0112137819740772, 0.9707859730744216, 0.9663220043672861],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].uncertainty.array[0, :3],
        [0.0010031773543652768] * 3,
        rtol=1e-12,
        atol=0,
    )


def test_binning_ratios_are_read_from_the_frame_and_order_table_headers(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`binRatioX`/`binRatioY` compare the frame's own binning to the order table's."""
    # ARRANGE
    recipe = _recipe(log)
    orderTable = _order_table(tmp_path, binx=2, biny=4)
    _stub_unpack_order_table(monkeypatch, _centre_pixels(12))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(
        np.random.default_rng(101).normal(loc=1000.0, scale=25.0, size=(12, 12)),
        binx=8,
        biny=8,
    )

    # ACT
    recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    assert recipe.binx == 8
    assert recipe.biny == 8
    assert recipe.binRatioX == pytest.approx(4.0)
    assert recipe.binRatioY == pytest.approx(2.0)


def test_a_nir_frame_without_binning_keywords_falls_back_to_unity(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """For the NIR arm, a missing `WIN_BINX`/`WIN_BINY` header falls back to `binx = biny = 1`."""
    # ARRANGE
    recipe = _recipe(log, arm="NIR")
    orderTable = _order_table(tmp_path)
    _stub_unpack_order_table(monkeypatch, _centre_pixels(12))
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(
        np.random.default_rng(101).normal(loc=1000.0, scale=25.0, size=(12, 12)),
        withHeaderBinning=False,
    )

    # ACT
    recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    assert recipe.binx == 1
    assert recipe.biny == 1
    assert recipe.binRatioX == pytest.approx(1.0)
    assert recipe.binRatioY == pytest.approx(1.0)


def test_first_pass_subsamples_chunks_above_ten_thousand_valid_pixels(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A chunk with more than 10000 valid pixels is subsampled with the seeded generator."""
    # ARRANGE
    recipe = _recipe(log)
    recipe.recipeSettings["centre-order-window"] = 1000
    orderTable = _order_table(tmp_path)
    shape = (300, 50)
    # A WIDE ENOUGH `centre-order-window` UNMASKS THE WHOLE ROW, SO THE FIRST
    # 256-ROW CHUNK CARRIES 256 * 50 = 12800 VALID PIXELS, ABOVE THE 10000
    # SUBSAMPLING THRESHOLD; THE SECOND (44-ROW) CHUNK STAYS UNDER IT.
    pixels = pd.DataFrame({"xcoord_centre": [25] * shape[0], "ycoord": list(range(shape[0]))})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(401).normal(loc=1000.0, scale=15.0, size=shape))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [0.9871754659806515, 1.001393804087737, 1.0134972034423777],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].data[150, :3],
        [1.023656861715834, 1.0026766556869509, 1.0051865766192047],
        rtol=1e-12,
        atol=0,
    )
    assert recipe.qc["qc_value"].tolist() == [
        "980.80125993991350696887821",
        "1000.156",
        "1019.592",
    ]


# ---------------------------------------------------------------------------
# mask_low_sens_pixels
# ---------------------------------------------------------------------------


def _order_edges(*, nrows: int = 5, ncols: int = 10, xLow: int = 2, xUp: int = 8, order: int = 10) -> pd.DataFrame:
    """Return a minimal order-edge table spanning every row of a small frame."""
    return pd.DataFrame(
        {
            "order": [order] * nrows,
            "xcoord_edgeup": [xUp] * nrows,
            "xcoord_edgelow": [xLow] * nrows,
            "ycoord": list(range(nrows)),
        }
    )


def test_low_sensitivity_pixels_are_masked_and_recorded_with_write_qc(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A clear low-sensitivity outlier is masked, counted, and logged when `writeQC=True`."""
    # ARRANGE
    recipe = _recipe(log)
    _stub_unpack_order_table(monkeypatch, _order_edges())
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    data = np.random.default_rng(501).normal(loc=100.0, scale=2.0, size=(5, 10))
    data[2, 4] = 0.0
    frame = _frame(data)

    # ACT
    maskedFrame, medianFluxDF = recipe.mask_low_sens_pixels(
        frame, "orders.fits", returnMedianOrderFlux=True, writeQC=True
    )

    # ASSERT
    expectedMask = np.zeros((5, 10), dtype=bool)
    expectedMask[2, 4] = True
    np.testing.assert_array_equal(maskedFrame.mask, expectedMask)

    assert list(recipe.qc.columns) == [
        "soxspipe_recipe",
        "qc_name",
        "qc_value",
        "qc_unit",
        "qc_order",
        "qc_comment",
        "obs_date_utc",
        "reduction_date_utc",
        "to_header",
    ]
    row = recipe.qc.iloc[0]
    assert row["qc_name"] == "N LOW SENS"
    assert row["qc_value"] == 1.0
    assert isinstance(row["qc_value"], float)
    assert row["qc_unit"] == "pixels"
    assert row["qc_comment"] == "Number of low-sensitivity pixels found in master flat"
    assert row["obs_date_utc"] == recipe.dateObs
    assert TIMESTAMP_PATTERN.fullmatch(row["reduction_date_utc"])
    assert pd.isna(row["to_header"])

    assert medianFluxDF.to_dict("records") == pytest.approx([{"order": 10, "medianFlux": 99.36284306903022}])

    printed = [message for level, message in log.messages if level == "print"]
    assert printed == [
        "\n# CLIPPING LOW-SENSITIVITY PIXELS AND SETTING INTER-ORDER AREA TO UNITY",
        "        1 low-sensitivity pixels added to bad-pixel mask",
    ]


def test_write_qc_false_adds_no_row_and_no_log_line(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """With `writeQC=False`, no QC row is appended and the pixel-count line is not printed."""
    # ARRANGE
    recipe = _recipe(log)
    _stub_unpack_order_table(monkeypatch, _order_edges())
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    data = np.random.default_rng(501).normal(loc=100.0, scale=2.0, size=(5, 10))
    data[2, 4] = 0.0
    frame = _frame(data)

    # ACT
    recipe.mask_low_sens_pixels(frame, "orders.fits", writeQC=False)

    # ASSERT
    assert len(recipe.qc) == 0
    printed = [message for level, message in log.messages if level == "print"]
    assert not any("low-sensitivity pixels added" in message for message in printed)


def test_negative_edge_coordinates_are_clamped_to_zero_without_raising(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Negative upper/lower edge coordinates are clamped to zero rather than wrapping the slice."""
    # ARRANGE
    recipe = _recipe(log)
    pixels = pd.DataFrame(
        {
            "order": [10, 10],
            "xcoord_edgeup": [-1, 3],
            "xcoord_edgelow": [-5, -1],
            "ycoord": [0, 1],
        }
    )
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.full((5, 10), 50.0))

    # ACT
    maskedFrame = recipe.mask_low_sens_pixels(frame, "orders.fits", writeQC=False)

    # ASSERT
    np.testing.assert_array_equal(maskedFrame.mask, np.zeros((5, 10), dtype=bool))


def test_the_axis_a_y_branch_masks_a_vertical_inter_order_band(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """When `axisA == "y"`, the inter-order mask and the median-flux sample use the transposed axes.

    The median-flux sample itself is read with `frame.data[b, l:u]` in both
    the `axisA == "x"` and the `axisA == "y"` branches of the source -- for
    the `y` branch this indexes the frame with the coordinates swapped
    relative to `interOrderMask[l:u, b] = 0` a few lines above it, so the
    sampled flux is not the same region that gets unmasked. This looks like
    a defect; it is pinned here, along with the resulting mask, as found.
    """
    # ARRANGE
    recipe = _recipe(log, axisA="y", axisB="x")
    pixels = pd.DataFrame(
        {
            "order": [10] * 5,
            "ycoord_edgeup": [8] * 5,
            "ycoord_edgelow": [2] * 5,
            "xcoord": [0, 1, 2, 3, 4],
        }
    )
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    data = np.random.default_rng(601).normal(loc=100.0, scale=2.0, size=(10, 5))
    data[4, 2] = 0.0
    frame = _frame(data)

    # ACT
    maskedFrame, medianFluxDF = recipe.mask_low_sens_pixels(
        frame, "orders.fits", returnMedianOrderFlux=True, writeQC=False
    )

    # ASSERT
    expectedMask = np.zeros((10, 5), dtype=bool)
    expectedMask[4, 2] = True
    np.testing.assert_array_equal(maskedFrame.mask, expectedMask)
    assert medianFluxDF.to_dict("records") == pytest.approx([{"order": 10, "medianFlux": 99.90835692418564}])


def test_return_median_order_flux_false_returns_the_frame_alone(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`returnMedianOrderFlux=False` returns the masked frame, not a `(frame, table)` pair."""
    # ARRANGE
    recipe = _recipe(log)
    _stub_unpack_order_table(monkeypatch, _order_edges())
    monkeypatch.setattr(mflatModule, "quicklook_image", lambda **kwargs: None)
    data = np.random.default_rng(501).normal(loc=100.0, scale=2.0, size=(5, 10))
    data[2, 4] = 0.0
    frame = _frame(data)

    # ACT
    result = recipe.mask_low_sens_pixels(frame, "orders.fits", returnMedianOrderFlux=False, writeQC=False)

    # ASSERT
    assert isinstance(result, CCDData)
    assert not isinstance(result, tuple)
