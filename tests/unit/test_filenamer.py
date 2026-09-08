"""Characterization tests for SOXSpipe filename rules."""

from __future__ import annotations

import pytest

from soxspipe.commonutils import filenamer
from tests.factories import synthetic_ccd

pytestmark = pytest.mark.unit


def test_filenamer_names_soxs_bias_from_resolved_headers(log: object) -> None:
    frame = synthetic_ccd(instrument="soxs")

    filename = filenamer(
        log=log,
        frame=frame,
        settings={"instrument": "soxs"},
    )

    assert filename == "2024.01.02T03.04.05.678_VIS_RO2_BIAS.fits"


def test_filenamer_names_prepared_bias_as_master_bias(log: object) -> None:
    frame = synthetic_ccd(instrument="soxs")
    frame.header["SXSPRE"] = True

    filename = filenamer(
        log=log,
        frame=frame,
        settings={"instrument": "soxs"},
    )

    assert filename == "2024.01.02T03.04.05.678_VIS_RO2_MBIAS.fits"


def test_filenamer_returns_none_for_processed_product(log: object) -> None:
    frame = synthetic_ccd(instrument="soxs")
    del frame.header["ESO DPR TYPE"]
    frame.header["ESO PRO TYPE"] = "REDUCED"

    assert filenamer(log=log, frame=frame, settings={"instrument": "soxs"}) is None


def test_filenamer_rejects_unknown_frame_type(log: object) -> None:
    frame = synthetic_ccd(
        instrument="soxs",
        headerOverrides={"DPR_TYPE": "UNKNOWN"},
    )

    with pytest.raises(TypeError, match="Frame type can't be determined"):
        filenamer(log=log, frame=frame, settings={"instrument": "soxs"})
