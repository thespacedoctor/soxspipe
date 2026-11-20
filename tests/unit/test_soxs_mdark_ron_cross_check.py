"""The SOXS NIR master dark's `QC RAW RON` is checked against the sample-count table, and sets nothing (DY-2341).

`soxs_mdark` logs a warning when the measured read noise is more than 20% from the table value at the dark's own
`ESO DET NDSAMPLES`.
"""

from __future__ import annotations

from typing import Any

import numpy as np
import pytest
from astropy import units as u
from astropy.io.fits import Header
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.detector_lookup import detector_lookup
from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.recipes import soxs_mdark
from tests.factories import StubInputFrames, instrument_header, qc_table

pytestmark = pytest.mark.unit

SOXS_SETTINGS = {"instrument": "soxs"}
# THE TABLE VALUE AT 10 SAMPLES
TABLE_RON_10 = 13.14
SHAPE = (64, 64)
# A 64 x 64 PAIR GIVES A MAD SIGMA PRECISE TO ABOUT 2%, SO THESE ARE SAFELY INSIDE AND OUTSIDE THE 20% LIMIT
NEAR_TABLE_FACTOR = 1.10
FAR_FROM_TABLE_FACTOR = 1.50


def _recipe(log: Any, arm: str = "NIR", settings: dict[str, str] = SOXS_SETTINGS) -> soxs_mdark:
    recipe = soxs_mdark.__new__(soxs_mdark)
    recipe.log = log
    recipe.kw = keyword_lookup(log=log, settings=settings).get
    recipe.settings = settings
    recipe.recipeName = "soxs-mdark"
    recipe.dateObs = "2024-01-02T03:04:05"
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.arm = arm
    recipe.detectorParams = detector_lookup(log=log, settings=settings).get(arm)
    return recipe


def _header(ndsamples: object = 10) -> Header:
    header = instrument_header(arm="NIR")
    if ndsamples is not None:
        header["HIERARCH ESO DET NDSAMPLES"] = ndsamples
    return header


def _warnings(log: Any) -> list[str]:
    return [message for level, message in log.messages if level == "warning"]


def _ron_warnings(log: Any) -> list[str]:
    return [message for message in _warnings(log) if "QC RAW RON" in message]


@pytest.mark.parametrize(
    ("rawRon", "isWarned"),
    [
        (TABLE_RON_10 * 1.19, False),
        (TABLE_RON_10 * 0.81, False),
        (TABLE_RON_10, False),
        (TABLE_RON_10 * 1.21, True),
        (TABLE_RON_10 * 0.79, True),
    ],
    ids=["19% above", "19% below", "on the table", "21% above", "21% below"],
)
def test_the_nir_dark_ron_warns_only_when_over_twenty_percent_from_the_table_value(
    log: Any, rawRon: float, isWarned: bool
) -> None:
    # ARRANGE
    recipe = _recipe(log)

    # ACT
    recipe._cross_check_nir_dark_ron(rawRon, _header(10))

    # ASSERT
    warnings = _ron_warnings(log)
    assert len(warnings) == int(isWarned)
    if isWarned:
        assert f"{rawRon:0.2f}" in warnings[0]
        assert "13.14" in warnings[0]
        assert "10 samples" in warnings[0]


def test_the_nir_dark_ron_is_checked_at_the_darks_own_sample_count(log: Any) -> None:
    # ARRANGE: 13.14 e- IS RIGHT AT 10 SAMPLES BUT 74% ABOVE THE 128-SAMPLE VALUE
    recipe = _recipe(log)

    # ACT
    recipe._cross_check_nir_dark_ron(TABLE_RON_10, _header(128))

    # ASSERT
    assert len(_ron_warnings(log)) == 1


@pytest.mark.parametrize("ndsamples", [None, 1, 129], ids=["missing card", "below the range", "above the range"])
def test_a_dark_whose_sample_count_the_table_cannot_serve_is_not_cross_checked(log: Any, ndsamples: object) -> None:
    # ARRANGE
    recipe = _recipe(log)

    # ACT
    recipe._cross_check_nir_dark_ron(99.0, _header(ndsamples))

    # ASSERT
    assert _warnings(log) == []


def test_an_xsh_nir_dark_is_not_cross_checked(log: Any) -> None:
    # ARRANGE: THE X-SHOOTER DETECTOR PARAMETERS HOLD NO TABLE
    recipe = _recipe(log, settings={"instrument": "xsh"})

    # ACT
    recipe._cross_check_nir_dark_ron(99.0, _header(10))

    # ASSERT
    assert _warnings(log) == []


def test_a_uvb_or_vis_dark_is_not_cross_checked(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log, arm="VIS")

    # ACT
    recipe._cross_check_nir_dark_ron(99.0, _header(10))

    # ASSERT
    assert _warnings(log) == []


def _record_checks(log: Any, sigma: float, ndsamples: int) -> soxs_mdark:
    """Run the master dark's QC step on a synthetic pair of raw darks with a read noise of `sigma`."""
    rng = np.random.default_rng(5)
    rawFrames = [
        CCDData(
            (100.0 + rng.normal(0.0, sigma, SHAPE)).astype(np.float32),
            unit=u.electron,
            meta=_header(ndsamples),
            mask=np.zeros(SHAPE, dtype=bool),
            uncertainty=StdDevUncertainty(np.ones(SHAPE), unit=u.electron),
        )
        for _ in range(2)
    ]
    master = rawFrames[0].copy()
    recipe = _recipe(log)
    recipe.inputFrames = StubInputFrames(rawFrames)
    recipe._record_dark_quality_checks(
        frame=master, masterMedianFluxLevel=100.0, meanFrameRon=sigma, masterRon=sigma / np.sqrt(2)
    )
    return recipe


def test_the_master_dark_qc_step_warns_when_the_measured_read_noise_is_far_from_the_table(log: Any) -> None:
    # ARRANGE / ACT
    recipe = _record_checks(log, TABLE_RON_10 * FAR_FROM_TABLE_FACTOR, ndsamples=10)

    # ASSERT: THE MEASURED VALUE IS STILL RECORDED
    assert len(_ron_warnings(log)) == 1
    rawRon = recipe.qc.set_index("qc_name").loc["RAW RON", "qc_value"]
    assert rawRon == pytest.approx(TABLE_RON_10 * FAR_FROM_TABLE_FACTOR, rel=0.05)


def test_the_master_dark_qc_step_is_quiet_when_the_measured_read_noise_is_near_the_table(log: Any) -> None:
    # ARRANGE / ACT
    _record_checks(log, TABLE_RON_10 * NEAR_TABLE_FACTOR, ndsamples=10)

    # ASSERT
    assert _ron_warnings(log) == []
