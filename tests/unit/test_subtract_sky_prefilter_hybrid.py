"""The `hybrid` pre-fit clip and first-fit-only masking of the offset-bin pre-filter flags (DY-2055 tuning)."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import pytest

from tests.unit.test_subtract_sky_bspline_characterization import _subtractor
from tests.unit.test_subtract_sky_prefit_clip_method import (
    COSMIC_ROW,
    METHOD_KEY,
    OBJECT_CENTRE,
    OBJECT_WIDTH,
    _clip,
    _order,
)

pytestmark = pytest.mark.unit

OFFSET_BIN_FLAG = "flagged_offset_bin_clipped"


def test_the_hybrid_method_keeps_the_rolling_clip_flags_for_the_object_search(log: Any) -> None:
    # ARRANGE
    order = _order()
    rolling = _clip(log, order)

    # ACT
    hybrid = _clip(log, order, **{METHOD_KEY: "hybrid"})

    # ASSERT: THE OBJECT SEARCH READS flagged_object_clipped, WHICH STAYS THE ROLLING CLIP'S
    assert np.array_equal(
        hybrid["flagged_object_clipped"].to_numpy(bool), rolling["flagged_object_clipped"].to_numpy(bool)
    )
    assert np.array_equal(hybrid["flagged_all_clipped"].to_numpy(bool), rolling["flagged_all_clipped"].to_numpy(bool))


def test_the_hybrid_method_records_the_offset_bin_flags_beside_the_rolling_flags(log: Any) -> None:
    # ARRANGE
    order = _order()
    offsetBin = _clip(log, order, **{METHOD_KEY: "offset_bin"})

    # ACT
    hybrid = _clip(log, order, **{METHOD_KEY: "hybrid"})

    # ASSERT: THE SAME PIXELS AS THE OFFSET-BIN METHOD FLAGS, MATCHED BY ROW INDEX
    expected = offsetBin["flagged_object_clipped"].reindex(hybrid.index).to_numpy(bool)
    assert hybrid[OFFSET_BIN_FLAG].to_numpy(bool).sum() > 0
    assert np.array_equal(hybrid[OFFSET_BIN_FLAG].to_numpy(bool), expected)


def _clip_and_mask(log: Any, order: pd.DataFrame, **skySettings: Any) -> pd.DataFrame:
    """The pre-fit clip of one order, then the object masking of `subtract`."""
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 3,
            "percentile_clipping_iterations": 15,
            "percentile_rolling_window_size": 35,
            "noise_rolling_window_size": 250,
        }
    )
    subtractor.recipeSettings["sky-subtraction"].update(skySettings)
    clipped = subtractor.get_over_sampled_sky_from_order(order.copy(), clipBPs=True, clipSlitEdge=False)
    return subtractor._mask_object_pixels([clipped], aggressive=True)[0]


def _far_from_object(frame: pd.DataFrame) -> np.ndarray:
    return ((frame["slit_position"] - OBJECT_CENTRE).abs() > 3 * OBJECT_WIDTH).to_numpy()


def test_the_default_object_masking_is_the_develop_object_masking(log: Any) -> None:
    # ARRANGE
    order = _order()
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(
        {
            "percentile_clipping_sigma": 3,
            "percentile_clipping_iterations": 15,
            "percentile_rolling_window_size": 35,
            "noise_rolling_window_size": 250,
        }
    )
    clipped = subtractor.get_over_sampled_sky_from_order(order.copy(), clipBPs=True, clipSlitEdge=False)
    expected = subtractor.clip_object_slit_positions([clipped.copy()], aggressive=True)[0]

    # ACT
    masked = _clip_and_mask(log, order)

    # ASSERT
    pd.testing.assert_frame_equal(masked, expected)


def test_the_hybrid_masking_still_masks_the_object_stripe(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    masked = _clip_and_mask(log, order, **{METHOD_KEY: "hybrid"})

    # ASSERT
    atStripe = (masked["slit_position"] - OBJECT_CENTRE).abs() < 0.3 * OBJECT_WIDTH
    assert masked.loc[atStripe, "flagged_all_clipped"].all()


def test_the_hybrid_masking_takes_the_offset_bin_flags_away_from_the_object(log: Any) -> None:
    # ARRANGE
    order = _order()
    rolling = _clip_and_mask(log, order)

    # ACT
    hybrid = _clip_and_mask(log, order, **{METHOD_KEY: "hybrid"})

    # ASSERT: AWAY FROM THE OBJECT THE FLAGS ARE THE OFFSET-BIN FLAGS, WHICH DIFFER FROM THE ROLLING FLAGS THERE
    far = _far_from_object(hybrid)
    hybridFlags = hybrid["flagged_object_clipped"].to_numpy(bool)[far]
    assert np.array_equal(hybridFlags, hybrid[OFFSET_BIN_FLAG].to_numpy(bool)[far])
    rollingFlags = rolling["flagged_object_clipped"].reindex(hybrid.index).to_numpy(bool)[far]
    assert not np.array_equal(hybridFlags, rollingFlags)
    assert hybrid.loc[COSMIC_ROW, "flagged_all_clipped"]


def test_the_hybrid_masking_keeps_the_edge_and_bad_pixel_clips(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    hybrid = _clip_and_mask(log, order, **{METHOD_KEY: "hybrid"})

    # ASSERT
    isEdgeOrBad = (hybrid["flagged_edge_clipped"] | hybrid["flagged_bad_pixel_clipped"]).to_numpy(bool)
    assert isEdgeOrBad.any()
    assert hybrid["flagged_all_clipped"].to_numpy(bool)[isEdgeOrBad].all()


FIRST_FIT_ONLY = {METHOD_KEY: "offset_bin", "prefilter_mask_first_fit_only": True}
FIRST_FIT_ONLY_FLAG = "flagged_prefilter_first_fit_only"


def _clip_mask_and_fit(log: Any, order: pd.DataFrame, **skySettings: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    """The masked order before the sky fit (a copy) and the fitted order."""
    masked = _clip_and_mask(log, order, **skySettings)
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(skySettings)
    fitted = subtractor.fit_bspline_curve_to_sky(masked.copy())[0]
    return masked, fitted.reindex(masked.index)


def test_first_fit_only_marks_the_offset_bin_flags_away_from_the_object(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    masked = _clip_and_mask(log, order, **FIRST_FIT_ONLY)

    # ASSERT: MARKED PIXELS ARE OFFSET-BIN FLAGS, AND NONE IS ON THE OBJECT STRIPE
    marked = masked[FIRST_FIT_ONLY_FLAG].to_numpy(bool)
    assert marked.sum() > 0
    assert masked[OFFSET_BIN_FLAG].to_numpy(bool)[marked].all()
    atStripe = ((masked["slit_position"] - OBJECT_CENTRE).abs() < 0.3 * OBJECT_WIDTH).to_numpy()
    assert not (marked & atStripe).any()


def test_without_the_setting_no_pixel_is_marked_first_fit_only(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    masked = _clip_and_mask(log, order, **{METHOD_KEY: "offset_bin"})

    # ASSERT
    assert FIRST_FIT_ONLY_FLAG not in masked


def test_first_fit_only_returns_the_marked_pixels_that_fit_the_sky_to_the_fit(log: Any) -> None:
    # ARRANGE
    order = _order()
    _, everyFit = _clip_mask_and_fit(log, order, **{METHOD_KEY: "offset_bin"})

    # ACT
    masked, firstOnly = _clip_mask_and_fit(log, order, **FIRST_FIT_ONLY)

    # ASSERT: THE MARKED PIXELS STAY OUT OF EVERY FIT BY DEFAULT; WITH THE SETTING MOST RETURN, THE COSMIC RAY DOES NOT
    marked = masked[FIRST_FIT_ONLY_FLAG].to_numpy(bool)
    marked[masked.index.get_loc(COSMIC_ROW)] = False
    assert everyFit["flagged_all_clipped"].to_numpy(bool)[marked].all()
    assert (~firstOnly["flagged_all_clipped"].to_numpy(bool)[marked]).mean() > 0.5
    assert firstOnly.loc[COSMIC_ROW, "flagged_all_clipped"]


def test_first_fit_only_keeps_the_object_stripe_out_of_the_fit(log: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT
    _, fitted = _clip_mask_and_fit(log, order, **FIRST_FIT_ONLY)

    # ASSERT
    atStripe = (fitted["slit_position"] - OBJECT_CENTRE).abs() < 0.3 * OBJECT_WIDTH
    assert fitted.loc[atStripe, "flagged_all_clipped"].all()


@pytest.mark.parametrize("value", ["yes", 1, None])
def test_a_bad_first_fit_only_value_raises_a_value_error_naming_the_setting(log: Any, value: Any) -> None:
    # ARRANGE
    order = _order()

    # ACT / ASSERT
    with pytest.raises(ValueError, match="prefilter_mask_first_fit_only"):
        _clip_and_mask(log, order, **{METHOD_KEY: "offset_bin", "prefilter_mask_first_fit_only": value})


def test_first_fit_only_leaves_the_marked_pixels_out_of_the_first_fit_and_in_a_later_one(log: Any) -> None:
    # ARRANGE
    order = _order()
    masked = _clip_and_mask(log, order, **FIRST_FIT_ONLY)
    marked = masked.index[masked[FIRST_FIT_ONLY_FLAG].to_numpy(bool) & (masked.index != COSMIC_ROW)]
    subtractor = _subtractor(log, arm="NIR", iterationLimit=25, pointsPerKnot=1000, noiseSigma=5)
    subtractor.recipeSettings["sky-subtraction"].update(FIRST_FIT_ONLY)
    fittedSamples = []
    original = subtractor._spline_fit_samples

    def recording(imageMapOrder, mask_all_clipped, starterKnots):
        fittedSamples.append(set(imageMapOrder.index[~mask_all_clipped.to_numpy(bool)]))
        return original(imageMapOrder, mask_all_clipped, starterKnots)

    subtractor._spline_fit_samples = recording

    # ACT
    subtractor.fit_bspline_curve_to_sky(masked.copy())

    # ASSERT
    assert len(fittedSamples) > 1
    assert not fittedSamples[0] & set(marked)
    assert len(fittedSamples[1] & set(marked)) > 0.5 * len(marked)


def test_the_hybrid_offset_bin_flags_follow_their_rows_when_the_input_is_shuffled(log: Any) -> None:
    # ARRANGE: SHUFFLED ROWS WITH A NON-MONOTONIC INDEX, SO BOTH WAVELENGTH SORTS MOVE THEM
    order = _order()
    shuffled = order.sample(frac=1.0, random_state=7)
    shuffled.index = shuffled.index[::-1] * 3 + 11
    offsetBin = _clip(log, shuffled, **{METHOD_KEY: "offset_bin"})

    # ACT
    hybrid = _clip(log, shuffled, **{METHOD_KEY: "hybrid"})

    # ASSERT
    expected = offsetBin["flagged_object_clipped"].reindex(hybrid.index).to_numpy(bool)
    assert expected.sum() > 0
    assert np.array_equal(hybrid[OFFSET_BIN_FLAG].to_numpy(bool), expected)
