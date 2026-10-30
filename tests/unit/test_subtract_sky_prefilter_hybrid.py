"""The `hybrid` pre-fit clip of the offset-bin pre-filter prototype (DY-2055 tuning)."""

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
