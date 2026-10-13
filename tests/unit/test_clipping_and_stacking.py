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
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
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
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
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


def _frames_masked_in_the_first_only(pixel: tuple[int, int]) -> list[CCDData]:
    """Five frames of value 10 with the values 12, 9, 10, 11, 10.5 at ``pixel``; the first frame flags it as bad."""
    frames = []
    for frameIndex, value in enumerate([12.0, 9.0, 10.0, 11.0, 10.5]):
        data = np.full((8, 8), 10.0)
        data[pixel] = value
        mask = np.zeros((8, 8), dtype=bool)
        mask[pixel] = frameIndex == 0
        frames.append(_frame(data, mask))
    return frames


def test_clip_and_stack_uncertainty_uses_the_frames_that_contributed_at_each_pixel(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: THE PIXEL IS BAD IN ONE OF FIVE FRAMES
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = _frames_masked_in_the_first_only((3, 3))

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", post_stack_clipping=False)

    # ASSERT: FOUR FRAMES CONTRIBUTED AT THE PIXEL, FIVE ELSEWHERE
    assert combined.data[3, 3] == pytest.approx(10.125)
    assert combined.uncertainty.array[3, 3] == pytest.approx(3.0 / np.sqrt(4))
    assert combined.uncertainty.array[0, 0] == pytest.approx(3.0 / np.sqrt(5))
    assert not combined.mask[3, 3]


def test_clip_and_stack_ignoring_input_masks_uses_every_frame_and_leaves_caller_masks_alone(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = _frames_masked_in_the_first_only((3, 3))

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", ignore_input_masks=True, post_stack_clipping=False)

    # ASSERT: ALL FIVE FRAMES CONTRIBUTE, AND THE CALLER'S FRAME KEEPS ITS MASK
    assert combined.data[3, 3] == pytest.approx(10.5)
    assert combined.uncertainty.array[3, 3] == pytest.approx(3.0 / np.sqrt(5))
    assert frames[0].mask[3, 3]


def test_clip_and_stack_keeps_the_value_of_a_pixel_flagged_in_every_frame(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: A STATIC BAD PIXEL, FLAGGED IN ALL FIVE FRAMES, WITH A DIFFERENT VALUE IN EACH
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = []
    for value in [12.0, 9.0, 10.0, 11.0, 10.5]:
        data = np.full((8, 8), 10.0)
        data[3, 3] = value
        mask = np.zeros((8, 8), dtype=bool)
        mask[3, 3] = True
        frames.append(_frame(data, mask))

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", post_stack_clipping=False)

    # ASSERT: THE PIXEL STAYS FLAGGED BUT HOLDS THE MEAN OF ALL FIVE VALUES AND A FINITE UNCERTAINTY
    assert combined.mask[3, 3]
    assert combined.data[3, 3] == pytest.approx(10.5)
    assert combined.uncertainty.array[3, 3] == pytest.approx(3.0 / np.sqrt(5))


def test_clip_and_stack_ignoring_input_masks_still_flags_a_pixel_flagged_in_every_frame(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: PIXEL (3, 3) IS FLAGGED IN EVERY FRAME, PIXEL (4, 4) IN THE FIRST ONLY
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = []
    for frameIndex in range(3):
        mask = np.zeros((8, 8), dtype=bool)
        mask[3, 3] = True
        mask[4, 4] = frameIndex == 0
        frames.append(_frame(np.full((8, 8), 10.0), mask))

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", ignore_input_masks=True, post_stack_clipping=False)

    # ASSERT
    assert combined.mask[3, 3]
    assert not combined.mask[4, 4]


def test_clip_and_stack_leaves_caller_masks_alone_when_honouring_input_masks(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = _frames_masked_in_the_first_only((3, 3))
    masksBefore = [frame.mask.copy() for frame in frames]

    # ACT
    recipe.clip_and_stack(frames, recipe="soxs_mbias", ignore_input_masks=False, post_stack_clipping=False)

    # ASSERT
    for frame, maskBefore in zip(frames, masksBefore, strict=True):
        assert np.array_equal(frame.mask, maskBefore)


def test_clip_and_stack_stacks_a_frame_with_no_mask_alongside_masked_frames(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: THE FIRST FRAME HAS NO MASK; THE SECOND FLAGS THE PIXEL CARRYING 12.0
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = _frames_masked_in_the_first_only((3, 3))
    frames[0], frames[1] = frames[1], frames[0]
    frames[0].mask = None

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", post_stack_clipping=False)

    # ASSERT: THE FLAGGED 12.0 IS LEFT OUT, AND THE RESULT IS A FULL FRAME
    assert combined.shape == (8, 8)
    assert combined.data[3, 3] == pytest.approx(10.125)
    assert combined.uncertainty.array[3, 3] == pytest.approx(3.0 / np.sqrt(4))


def test_clip_and_stack_leaves_masked_values_out_of_the_clipping_statistics(
    tmp_path: Path,
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # ARRANGE: THREE FRAMES FLAG A WILD VALUE (30) AT THE PIXEL; OF THE FOUR UNFLAGGED, ONE IS A TRUE OUTLIER
    recipe = _recipe(tmp_path, log)
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **kwargs: None)
    frames = []
    for value, isFlagged in [(30.0, True)] * 3 + [(9.8, False), (10.0, False), (10.2, False), (30.0, False)]:
        data = np.full((8, 8), 10.0)
        data[3, 3] = value
        mask = np.zeros((8, 8), dtype=bool)
        mask[3, 3] = isFlagged
        frames.append(_frame(data, mask))

    # ACT
    combined = recipe.clip_and_stack(frames, recipe="soxs_mbias", post_stack_clipping=False)

    # ASSERT: THE FLAGGED 30s DO NOT SHIELD THE UNFLAGGED 30, AND DO NOT CLIP THE GOOD VALUES
    assert combined.data[3, 3] == pytest.approx(10.0)
    assert combined.uncertainty.array[3, 3] == pytest.approx(3.0 / np.sqrt(3))
