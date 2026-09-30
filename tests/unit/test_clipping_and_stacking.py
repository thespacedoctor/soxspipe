"""Controlled clipping and stacking contracts for the base recipe."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories import instrument_header, pipeline_settings

pytestmark = pytest.mark.unit


def _recipe(tmp_path: Path, log: object) -> base_recipe:
    recipe = object.__new__(base_recipe)
    recipe.log = log
    recipe.settings = pipeline_settings(tmp_path)
    recipe.kw = keyword_lookup(log=log, settings=recipe.settings).get
    recipe.arm = "VIS"
    recipe.detectorParams = {}
    recipe.imageType = "BIAS"
    recipe.recipeSettings = {
        "stacked-clipping-sigma": 3.0,
        "stacked-clipping-iterations": 2,
        "frame-clipping-sigma": 3.0,
        "frame-clipping-iterations": 2,
    }
    recipe.verbose = False
    recipe.debug = False
    return recipe


def _frame(data: np.ndarray, mask: np.ndarray | None = None) -> CCDData:
    shape = data.shape
    uncertainty = StdDevUncertainty(np.full(shape, 3.0), unit=u.electron)
    return CCDData(
        data.astype(np.float32),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros(shape, dtype=bool) if mask is None else mask.copy(),
        uncertainty=uncertainty,
    )


def test_clip_and_stack_rejects_empty_input(tmp_path: Path, log: object) -> None:
    recipe = _recipe(tmp_path, log)

    with pytest.raises(ValueError, match="No frames were sent"):
        recipe.clip_and_stack([], recipe="soxs_mbias")


def test_clip_and_stack_returns_single_frame_unchanged(
    tmp_path: Path,
    log: object,
) -> None:
    recipe = _recipe(tmp_path, log)
    frame = _frame(np.full((8, 8), 10.0))

    combined = recipe.clip_and_stack([frame], recipe="soxs_mbias")

    assert combined is frame


def test_clip_and_stack_rejects_outlier_and_propagates_mask_and_uncertainty(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None
    )
    commonMask = np.zeros((8, 8), dtype=bool)
    commonMask[0, 0] = True
    firstOnlyMask = commonMask.copy()
    firstOnlyMask[0, 1] = True
    arrays = [np.full((8, 8), 10.0) for _ in range(3)]
    arrays[2][2, 3] = 1000.0
    frames = [
        _frame(arrays[0], firstOnlyMask),
        _frame(arrays[1], commonMask),
        _frame(arrays[2], commonMask),
    ]

    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias")

    assert combined.data.dtype == np.float32
    assert combined.data[2, 3] == pytest.approx(10.0)
    assert combined.mask[0, 0]
    assert not combined.mask[0, 1]
    assert combined.uncertainty.array[1, 1] == pytest.approx(3.0 / np.sqrt(3))
    assert combined.uncertainty.array[2, 3] == pytest.approx(3.0 / np.sqrt(2))


def test_clip_and_stack_can_disable_post_stack_clipping(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    recipe = _recipe(tmp_path, log)
    recipe.recipeSettings = {
        "stacked-clipping-sigma": 3.0,
        "stacked-clipping-iterations": 2,
    }
    monkeypatch.setattr(
        "soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None
    )
    data = np.zeros((8, 8), dtype=np.float32)
    data[4, 4] = 100.0
    frames = [_frame(data) for _ in range(3)]

    combined = recipe.clip_and_stack(
        frames,
        recipe="soxs_mbias",
        post_stack_clipping=False,
    )

    assert combined.data[4, 4] == pytest.approx(100.0)
    assert not combined.mask[4, 4]
