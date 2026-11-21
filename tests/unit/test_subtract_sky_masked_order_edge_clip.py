"""The slit-edge clip is skipped on orders the slit-illumination exclude mask covers (DY-1280)."""

from __future__ import annotations

import importlib
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

import soxspipe.commonutils.toolkit as toolkit
from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.factories import instrument_header

pytestmark = pytest.mark.unit

subtractSkyModule = importlib.import_module("soxspipe.commonutils.subtract_sky")

# TWO ORDERS OF ELEVEN SLIT POSITIONS (-5 TO 5), ONE DETECTOR ROW PER ORDER
SLIT_POSITIONS = np.arange(-5.0, 6.0)
FRAME_SHAPE = (2, len(SLIT_POSITIONS))
MASKED_ORDER = 13
UNMASKED_ORDER = 14
# A FRACTION OF 0.1 OF THE 10 ARCSEC RANGE CLIPS THE ROWS BEYOND +-4, SO ONLY THE END ROWS GO
EDGE_FRACTION = 0.1
EDGE_CLIPPED_SLIT_POSITIONS = {-5.0, 5.0}
# THE EXCLUDE MASK REMOVES THE SLIT ENDS (+-4 AND BEYOND) OF THE MASKED ORDER, WHICH ARE THEN NOT ITS EDGE ROWS
EXCLUDED_SLIT_POSITIONS = {-5.0, -4.0, 4.0, 5.0}


def _map_table() -> pd.DataFrame:
    frames = []
    for row, order in enumerate([MASKED_ORDER, UNMASKED_ORDER]):
        frames.append(
            pd.DataFrame(
                {
                    "order": order,
                    "x": np.arange(len(SLIT_POSITIONS)),
                    "y": row,
                    "wavelength": 1500.0 + 100.0 * row + np.arange(len(SLIT_POSITIONS)),
                    "slit_position": SLIT_POSITIONS,
                    "mask": False,
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


@pytest.fixture(autouse=True)
def collaborators(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Replace the file-reading and plotting collaborators; the keyword and detector lookups stay real."""
    monkeypatch.setattr(
        subtractSkyModule,
        "twoD_disp_map_image_to_dataframe",
        lambda **_: (_map_table(), np.zeros(FRAME_SHAPE, dtype=bool)),
    )
    monkeypatch.setattr(subtractSkyModule, "quicklook_image", lambda **_: None)
    monkeypatch.setattr(subtractSkyModule, "filenamer", lambda **_: "SYNTHETIC.fits")
    monkeypatch.setattr(toolkit, "utility_setup", lambda **_: (str(tmp_path / "qc"), str(tmp_path / "products")))


def _frame() -> CCDData:
    header = instrument_header(arm="VIS", overrides={"ESO INS VISE NAME": "SLIT1.0"})
    return CCDData(
        np.full(FRAME_SHAPE, 10.0),
        unit=u.electron,
        mask=np.zeros(FRAME_SHAPE, dtype=bool),
        uncertainty=StdDevUncertainty(np.full(FRAME_SHAPE, 2.0)),
        meta=header,
    )


def _exclude_mask() -> np.ndarray:
    """Return a mask over the slit ends of the masked order only."""
    excludeMask = np.zeros(FRAME_SHAPE, dtype=bool)
    isEnd = np.isin(SLIT_POSITIONS, list(EXCLUDED_SLIT_POSITIONS))
    excludeMask[0, isEnd] = True
    return excludeMask


def _fitted_orders(log: Any, **kwargs: Any) -> dict[int, pd.DataFrame]:
    """Run ``subtract`` with the clipping and fitting stages stubbed and return each order as the fit received it."""
    subtractor = subtract_sky(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={
            "sky-subtraction": {
                "bspline_order": 3,
                "clip-slit-edge-fraction": EDGE_FRACTION,
                "aggressive_object_masking": False,
                # THE FIT STAGES ARE STUBBED, SO THE OBJECT-WING REFIT (DY-1279) HAS NO CLIPPING COLUMNS TO READ
                "object_wing_correction": False,
                "sky_model_qc_plot": False,
                "percentile_clipping_sigma": 5,
                "percentile_clipping_iterations": 1,
                "percentile_rolling_window_size": 3,
            }
        },
        objectFrame=_frame(),
        twoDMap="two-d-map.fits",
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        **kwargs,
    )
    fittedOrders: dict[int, pd.DataFrame] = {}

    def record_fit(frame: pd.DataFrame) -> tuple[Any, ...]:
        fittedOrders[int(frame["order"].iloc[0])] = frame
        return (
            frame.assign(sky_model=4.0, sky_subtracted_flux=6.0, error=2.0),
            (None, None, 3),
            np.array([1500.0]),
            np.array([1.0, 2.0]),
            0.5,
        )

    # THE ROLLING-WINDOW CLIP IS NOT UNDER TEST, SO IT LEAVES THE FLAGS AS THE EDGE CLIP SET THEM
    subtractor.rolling_window_clipping = lambda imageMapOrderDF, **_: imageMapOrderDF
    subtractor.clip_object_slit_positions = lambda orders, **_: orders
    subtractor.fit_bspline_curve_to_sky = record_fit
    subtractor.subtract()
    return fittedOrders


def _edge_clipped_slit_positions(order: pd.DataFrame) -> set[float]:
    return set(order.loc[order["flagged_edge_clipped"], "slit_position"].tolist())


def test_an_order_the_exclude_mask_covers_is_not_edge_clipped(log: Any) -> None:
    fittedOrders = _fitted_orders(log, excludeMask=_exclude_mask())

    assert set(fittedOrders[MASKED_ORDER]["slit_position"]) == set(SLIT_POSITIONS) - EXCLUDED_SLIT_POSITIONS
    assert not fittedOrders[MASKED_ORDER]["flagged_edge_clipped"].any()
    assert not fittedOrders[MASKED_ORDER]["flagged_all_clipped"].any()


def test_an_order_the_exclude_mask_does_not_cover_is_still_edge_clipped(log: Any) -> None:
    fittedOrders = _fitted_orders(log, excludeMask=_exclude_mask())

    assert _edge_clipped_slit_positions(fittedOrders[UNMASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS


def test_without_an_exclude_mask_every_order_is_edge_clipped(log: Any) -> None:
    fittedOrders = _fitted_orders(log)

    assert _edge_clipped_slit_positions(fittedOrders[MASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS
    assert _edge_clipped_slit_positions(fittedOrders[UNMASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS


def test_an_all_false_exclude_mask_leaves_every_order_edge_clipped(log: Any) -> None:
    fittedOrders = _fitted_orders(log, excludeMask=np.zeros(FRAME_SHAPE, dtype=bool))

    assert _edge_clipped_slit_positions(fittedOrders[MASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS
    assert _edge_clipped_slit_positions(fittedOrders[UNMASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS


# THE COVERAGE MASK TRIMS THE HIGH END ONLY (DY-2290): POSITIONS 4 AND 5. THE REST OF THE ORDER IS -5 TO 3, A RANGE OF
# 8, SO A FRACTION OF 0.1 CLIPS THE ROWS BEYOND -4.2 AND 2.2
UNCOVERED_SLIT_POSITIONS = {4.0, 5.0}
UNCOVERED_ORDER_EDGE_CLIPPED = {-5.0, 3.0}


def _uncovered_mask(row: int) -> np.ndarray:
    """Return a mask over the high end of the slit of one order."""
    uncoveredMask = np.zeros(FRAME_SHAPE, dtype=bool)
    uncoveredMask[row, np.isin(SLIT_POSITIONS, list(UNCOVERED_SLIT_POSITIONS))] = True
    return uncoveredMask


def test_an_order_the_uncovered_mask_trims_keeps_its_edge_clip(log: Any) -> None:
    # ARRANGE / ACT: THE COVERAGE MASK TRIMS ONE END ONLY, SO THE OTHER END STILL NEEDS ITS EDGE CLIP
    fittedOrders = _fitted_orders(log, uncoveredMask=_uncovered_mask(0))

    # ASSERT
    assert set(fittedOrders[MASKED_ORDER]["slit_position"]) == set(SLIT_POSITIONS) - UNCOVERED_SLIT_POSITIONS
    assert _edge_clipped_slit_positions(fittedOrders[MASKED_ORDER]) == UNCOVERED_ORDER_EDGE_CLIPPED
    assert _edge_clipped_slit_positions(fittedOrders[UNMASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS


def test_the_exclude_mask_still_switches_off_the_edge_clip_of_an_order_the_uncovered_mask_also_trims(log: Any) -> None:
    # ARRANGE / ACT: THE EXCLUDE MASK HAS TRIMMED BOTH ENDS OF THE ORDER, AS IT DID BEFORE THE COVERAGE MASK
    fittedOrders = _fitted_orders(log, excludeMask=_exclude_mask(), uncoveredMask=_uncovered_mask(0))

    # ASSERT
    assert set(fittedOrders[MASKED_ORDER]["slit_position"]) == (
        set(SLIT_POSITIONS) - EXCLUDED_SLIT_POSITIONS - UNCOVERED_SLIT_POSITIONS
    )
    assert not fittedOrders[MASKED_ORDER]["flagged_edge_clipped"].any()
    assert _edge_clipped_slit_positions(fittedOrders[UNMASKED_ORDER]) == EDGE_CLIPPED_SLIT_POSITIONS


def test_an_uncovered_mask_of_the_wrong_shape_raises(log: Any) -> None:
    with pytest.raises(ValueError, match="shape"):
        _fitted_orders(log, uncoveredMask=np.zeros((2, 3), dtype=bool))
