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


@pytest.mark.parametrize(
    ("instrument", "headerOverrides", "prepared", "expected"),
    (
        ("xsh", {"DET_READ_SPEED": 1}, False, "2024.01.02T03.04.05.678_VIS_ROSPEED1_BIAS.fits"),
        ("xsh", {"DET_READ_SPEED": "100kHz"}, False, "2024.01.02T03.04.05.678_VIS_SLOW_BIAS.fits"),
        ("xsh", {"DET_READ_SPEED": "400kHz"}, False, "2024.01.02T03.04.05.678_VIS_FAST_BIAS.fits"),
        (
            "soxs",
            {"DPR_TYPE": "DARK", "SXSPRE": True},
            True,
            "2024.01.02T03.04.05.678_VIS_RO2_MDARK.fits",
        ),
        (
            "soxs",
            {"DPR_TYPE": "LAMP,FLAT,Q", "DPR_TECH": "ECHELLE,SLIT"},
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_QLAMP_FLAT_SLIT.fits",
        ),
        (
            "soxs",
            {"DPR_TYPE": "LAMP,FLAT,D", "DPR_TECH": "ECHELLE,SLIT"},
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_DLAMP_FLAT_SLIT.fits",
        ),
        (
            "soxs",
            {"DPR_TYPE": "LAMP,FMTCHK", "DPR_TECH": "ECHELLE,PINHOLE"},
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_ARC_ONEPIN.fits",
        ),
        (
            "soxs",
            {"DPR_TYPE": "LAMP,ORDERDEF", "DPR_TECH": "ECHELLE,MULTI-PINHOLE"},
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_FLAT_MULTIPIN.fits",
        ),
        (
            "soxs",
            {
                "DPR_TYPE": "OBJECT",
                "DPR_TECH": "ECHELLE,SLIT,STARE",
                "OBJECT": "HD-123  45",
            },
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_OBJECT_STARE_HD_123_45_SLIT.fits",
        ),
        (
            "soxs",
            {
                "DPR_TYPE": "STD,FLUX",
                "DPR_TECH": "ECHELLE,SLIT,NODDING",
                "OBJECT": "EG 274",
            },
            False,
            "2024.01.02T03.04.05.678_VIS_RO2_STD_FLUX_STARE_EG_274_SLIT.fits",
        ),
    ),
)
def test_filenamer_covers_supported_observation_and_calibration_names(
    log: object,
    instrument: str,
    headerOverrides: dict[str, object],
    prepared: bool,
    expected: str,
) -> None:
    """Resolved headers produce the public filename for each supported class."""
    frame = synthetic_ccd(
        instrument=instrument,
        prepared=prepared,
        headerOverrides=headerOverrides,
    )

    filename = filenamer(log=log, frame=frame, settings={"instrument": instrument})

    assert filename == expected


def test_filenamer_rejects_unknown_xshooter_readout_mode(log: object) -> None:
    """X-shooter filenames reject readout values outside the known conventions."""
    frame = synthetic_ccd(instrument="xsh", headerOverrides={"DET_READ_SPEED": "turbo"})

    with pytest.raises(LookupError, match="Cound not parse readout mode"):
        filenamer(log=log, frame=frame, settings={"instrument": "xsh"})
