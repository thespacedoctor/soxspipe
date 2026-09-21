"""Characterization of `soxs_mflat.normalise_flats` and `find_uvb_overlap_order_and_scale`.

These tests pin what `normalise_flats` does today at its centre-trace mask
guards, its debug-mode plotting, its two "every pixel excluded" edge cases,
and its handling of frames with no uncertainty array, before a later commit
splits `soxs_mflat.py`'s functions into smaller methods. They also pin that
`find_uvb_overlap_order_and_scale` cannot run against the real
`normalise_flats` it calls -- evidence that the method is dead code, not
merely untested. Every defect below is pinned with a docstring noting it is
pinned as found, not as intended.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy.nddata import CCDData

from tests.unit.test_soxs_mflat_characterization import soxs_mflat_module
from tests.unit.test_soxs_mflat_worker_characterization import (
    _frame,
    _order_table,
    _stub_unpack_order_table,
)
from tests.unit.test_soxs_mflat_worker_characterization import (
    _recipe as _worker_recipe,
)

pytestmark = pytest.mark.unit


# ---------------------------------------------------------------------------
# normalise_flats -- mask guards, debug plotting, all-masked edge cases,
# and frames with no uncertainty array
# ---------------------------------------------------------------------------


def test_axis_a_x_mask_guard_skips_out_of_bounds_centre_trace_rows(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A centre-trace row outside the frame is skipped by the `0 <= y < mask.shape[0]` guard.

    Without the guard, `y = -1` would silently wrap around to the frame's
    *last* row under NumPy's negative-index rules, and `y = 10` would raise
    an `IndexError` against the six-row frame -- the guard converts both into
    a silent no-op instead. Only the in-bounds `y = 2` row ends up unmasked.
    """
    # ARRANGE
    recipe = _worker_recipe(log)
    orderTable = _order_table(tmp_path, name="orders_x_guard.fits")
    pixels = pd.DataFrame({"xcoord_centre": [3, 3, 3], "ycoord": [-1, 2, 10]})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(701).normal(loc=1000.0, scale=25.0, size=(6, 8)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[2, :4],
        [0.9835910548425146, 0.9978603888685283, 0.9985020543034777, 1.0014979456965223],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].data[0, :4],
        [0.9963679919091896, 0.9639170142065926, 0.9983043501282696, 0.9432489313744399],
        rtol=1e-12,
        atol=0,
    )


def test_axis_a_y_mask_guard_skips_out_of_bounds_centre_trace_columns(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """When `axisA == "y"`, an out-of-frame column is skipped by the `0 <= x < mask.shape[1]` guard."""
    # ARRANGE
    recipe = _worker_recipe(log, axisA="y", axisB="x")
    orderTable = _order_table(tmp_path, name="orders_y_guard.fits")
    pixels = pd.DataFrame({"ycoord_centre": [3, 3, 3], "xcoord": [-1, 2, 10]})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(702).normal(loc=1000.0, scale=25.0, size=(8, 6)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    np.testing.assert_allclose(
        normalised[0].data[:4, 2],
        [1.0490991135341885, 1.0218096499192484, 1.0044998974449544, 0.9955001025550454],
        rtol=1e-12,
        atol=0,
    )
    np.testing.assert_allclose(
        normalised[0].data[:4, 0],
        [0.9947733823369387, 1.0520341908481134, 1.025253858547396, 1.0083702824919138],
        rtol=1e-12,
        atol=0,
    )


def test_debug_mode_draws_an_extra_quicklook_of_the_masked_first_frame(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`recipe.debug = True` draws one extra quicklook of the first frame with the centre mask applied."""
    # ARRANGE
    recipe = _worker_recipe(log)
    recipe.debug = True
    orderTable = _order_table(tmp_path)
    pixels = pd.DataFrame({"xcoord_centre": [6] * 12, "ycoord": list(range(12))})
    _stub_unpack_order_table(monkeypatch, pixels)
    quicklookCalls: list[dict[str, Any]] = []
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: quicklookCalls.append(kwargs))
    frame = _frame(np.random.default_rng(801).normal(loc=1000.0, scale=25.0, size=(12, 12)))

    # ACT
    recipe.normalise_flats([frame], str(orderTable), lamp="_DLAMP")

    # ASSERT
    assert len(quicklookCalls) == 2
    debugCall = quicklookCalls[0]
    assert debugCall["title"] == "Example input flat frame with order centre mask applied _DLAMP"
    assert debugCall["stdWindow"] == 6
    assert debugCall["show"] is False
    assert debugCall["ext"] is None
    assert debugCall["surfacePlot"] is True

    expectedMask = np.ones((12, 12), dtype=bool)
    expectedMask[:, 5:7] = False
    np.testing.assert_array_equal(debugCall["CCDObject"].mask, expectedMask)


def _empty_pixels() -> pd.DataFrame:
    """Return an order-table pixel table with no centre-trace rows at all."""
    return pd.DataFrame({"xcoord_centre": pd.array([], dtype="int64"), "ycoord": pd.array([], dtype="int64")})


def test_first_pass_with_no_centre_trace_pixels_raises_value_error_on_concatenate(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """With no centre-trace pixels, the mask stays fully `True` and every pixel is excluded.

    `sample_chunks` never receives an array, so `np.concatenate(sample_chunks)`
    is called on an empty list. Pinned as found, not as intended: a
    genuinely empty order table crashes with a NumPy-internal message rather
    than a domain-meaningful error.
    """
    # ARRANGE
    recipe = _worker_recipe(log)
    orderTable = _order_table(tmp_path, name="orders_empty_first.fits")
    _stub_unpack_order_table(monkeypatch, _empty_pixels())
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(801).normal(loc=1000.0, scale=25.0, size=(12, 12)))

    # ACT / ASSERT
    with pytest.raises(ValueError, match="need at least one array to concatenate"):
        recipe.normalise_flats([frame], str(orderTable))


def test_second_pass_with_no_centre_trace_pixels_raises_unbound_local_error_on_norm_level(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """With no centre-trace pixels, the second pass never assigns `norm_level` before reading it.

    `chunk_vals` stays empty, so the `if chunk_vals:` block that assigns
    `norm_level` never runs; the very next line divides by it anyway. Pinned
    as found, not as intended.
    """
    # ARRANGE
    recipe = _worker_recipe(log)
    orderTable = _order_table(tmp_path, name="orders_empty_second.fits")
    _stub_unpack_order_table(monkeypatch, _empty_pixels())
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame(np.random.default_rng(801).normal(loc=1000.0, scale=25.0, size=(12, 12)))
    masterFlat = _frame(1.0 + np.random.default_rng(999).normal(loc=0.0, scale=0.01, size=(12, 12)))

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError, match="norm_level"):
        recipe.normalise_flats([frame], str(orderTable), firstPassMasterFlat=masterFlat)


def _frame_no_uncertainty(data: np.ndarray, *, binx: int = 1, biny: int = 1) -> CCDData:
    """Return a mutable CCD frame carrying no uncertainty array."""
    frame = CCDData(data, unit="electron")
    frame.mask = np.zeros(data.shape, dtype=bool)
    frame.uncertainty = None
    frame.header["WIN_BINX"] = binx
    frame.header["WIN_BINY"] = biny
    return frame


def test_first_pass_normalises_data_and_leaves_a_missing_uncertainty_as_none(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A frame with no uncertainty array normalises its data and keeps `uncertainty` `None`."""
    # ARRANGE
    recipe = _worker_recipe(log)
    orderTable = _order_table(tmp_path, name="orders_no_uncertainty_first.fits")
    pixels = pd.DataFrame({"xcoord_centre": [6] * 12, "ycoord": list(range(12))})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame_no_uncertainty(np.random.default_rng(801).normal(loc=1000.0, scale=25.0, size=(12, 12)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable))

    # ASSERT
    assert normalised[0].uncertainty is None
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [0.9775311708978731, 0.9840835489760086, 1.0434260284288288],
        rtol=1e-12,
        atol=0,
    )


def test_second_pass_normalises_data_and_leaves_a_missing_uncertainty_as_none(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """A frame with no uncertainty array normalises its data in the second pass too, keeping it `None`."""
    # ARRANGE
    recipe = _worker_recipe(log)
    orderTable = _order_table(tmp_path, name="orders_no_uncertainty_second.fits")
    pixels = pd.DataFrame({"xcoord_centre": [6] * 12, "ycoord": list(range(12))})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    frame = _frame_no_uncertainty(np.random.default_rng(801).normal(loc=1000.0, scale=25.0, size=(12, 12)))
    masterFlat = _frame_no_uncertainty(1.0 + np.random.default_rng(999).normal(loc=0.0, scale=0.01, size=(12, 12)))

    # ACT
    normalised = recipe.normalise_flats([frame], str(orderTable), firstPassMasterFlat=masterFlat)

    # ASSERT
    assert normalised[0].uncertainty is None
    np.testing.assert_allclose(
        normalised[0].data[0, :3],
        [0.9766215022207355, 0.9831677828021609, 1.0424550394688552],
        rtol=1e-12,
        atol=0,
    )


# ---------------------------------------------------------------------------
# find_uvb_overlap_order_and_scale
# ---------------------------------------------------------------------------


class _SingleOrderTableCollection:
    """Order-table lookup double returning the same single path for every filter."""

    def __init__(self, path: str) -> None:
        self._path = path

    def filter(self, **_: object) -> _SingleOrderTableCollection:
        return self

    def files_filtered(self, *, include_path: bool) -> list[str]:
        assert include_path is True
        return [self._path]


def test_find_uvb_overlap_order_and_scale_cannot_run_against_the_real_normalise_flats(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    """`find_uvb_overlap_order_and_scale` calls `normalise_flats` expecting a `(frames, table)` pair.

    The real `normalise_flats` returns a single list of normalised frames,
    not a `(frames, table)` pair, so unpacking its result with one input flat
    per lamp raises `ValueError` before any of the per-order flux comparison
    logic runs. Nothing in the package calls this method; this test is the
    evidence that it is dead code, not merely untested. Pinned as found, not
    as intended.
    """
    # ARRANGE
    recipe = _worker_recipe(log, arm="UVB")
    orderTable = _order_table(tmp_path, name="orders_overlap.fits")
    recipe.inputFrames = _SingleOrderTableCollection(str(orderTable))
    pixels = pd.DataFrame({"xcoord_centre": [6] * 12, "ycoord": list(range(12))})
    _stub_unpack_order_table(monkeypatch, pixels)
    monkeypatch.setattr(soxs_mflat_module, "quicklook_image", lambda **kwargs: None)
    dflat = _frame(np.random.default_rng(1).normal(loc=1000.0, scale=25.0, size=(12, 12)))
    qflat = _frame(np.random.default_rng(2).normal(loc=1000.0, scale=25.0, size=(12, 12)))

    # ACT / ASSERT
    with pytest.raises(ValueError, match="not enough values to unpack"):
        recipe.find_uvb_overlap_order_and_scale([dflat], [qflat])
