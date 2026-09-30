"""Characterization tests for SOXSpipe filename rules."""

from __future__ import annotations

import pytest
from astropy.wcs import WCS

from soxspipe.commonutils import detector_lookup, filenamer, keyword_lookup
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


STAMP = "2024.01.02T03.04.05.678"
SOXS_SETTINGS = {"instrument": "soxs"}


def _name(log: object, frame: object, **kwargs: object) -> str | None:
    return filenamer(log=log, frame=frame, **kwargs)


def _frame_with_cdelt(cdelt: tuple[float, float]) -> object:
    frame = synthetic_ccd()
    wcs = WCS(naxis=2)
    wcs.wcs.cdelt = list(cdelt)
    frame.wcs = wcs
    return frame


@pytest.mark.parametrize(
    ("cdelt", "binning"),
    (
        ((2.0, 2.0), "2X2"),
        ((2.7, 2.9), "2X2"),
        ((-1.9, 3.0), "-1X3"),
        ((0.5, 0.5), "0X0"),
        ((1, 4), "1X4"),
    ),
)
def test_filenamer_adds_binning_from_the_wcs_truncating_towards_zero(
    log: object, cdelt: tuple[float, float], binning: str
) -> None:
    frame = _frame_with_cdelt(cdelt)

    filename = _name(log, frame, settings=SOXS_SETTINGS)

    assert filename == f"{STAMP}_VIS_{binning}_RO2_BIAS.fits"


def test_filenamer_omits_the_readout_fragment_without_a_readout_keyword(
    log: object,
) -> None:
    frame = synthetic_ccd()
    del frame.header["ESO DET READ CURID"]

    assert _name(log, frame, settings=SOXS_SETTINGS) == f"{STAMP}_VIS_BIAS.fits"


def test_filenamer_uses_a_string_readout_value_verbatim_for_soxs(log: object) -> None:
    frame = synthetic_ccd(headerOverrides={"DET_READ_SPEED": "fast"})

    assert _name(log, frame, settings=SOXS_SETTINGS) == f"{STAMP}_VIS_ROFAST_BIAS.fits"


def test_filenamer_reads_a_float_one_as_the_xshooter_speed_one_mode(
    log: object,
) -> None:
    frame = synthetic_ccd(instrument="xsh", headerOverrides={"DET_READ_SPEED": 1.0})

    filename = _name(log, frame, settings={"instrument": "xsh"})

    assert filename == f"{STAMP}_VIS_ROSPEED1_BIAS.fits"


def test_filenamer_lower_cases_the_xshooter_readout_value_before_matching(
    log: object,
) -> None:
    frame = synthetic_ccd(
        instrument="xsh", headerOverrides={"DET_READ_SPEED": "Slow 100K"}
    )

    filename = _name(log, frame, settings={"instrument": "xsh"})

    assert filename == f"{STAMP}_VIS_SLOW_BIAS.fits"


def test_filenamer_prints_the_readout_value_before_rejecting_it(log: object) -> None:
    frame = synthetic_ccd(instrument="xsh", headerOverrides={"DET_READ_SPEED": "turbo"})

    with pytest.raises(LookupError):
        _name(log, frame, settings={"instrument": "xsh"})

    assert ("print", "turbo") in log.messages


@pytest.mark.parametrize(
    ("headerOverrides", "prepared", "expected"),
    (
        ({"DPR_TYPE": "DARK"}, False, f"{STAMP}_VIS_RO2_DARK.fits"),
        ({"DPR_TYPE": "bias"}, False, f"{STAMP}_VIS_RO2_BIAS.fits"),
        (
            {"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,SLIT", "SXSPRE": True},
            True,
            f"{STAMP}_VIS_RO2_MFLAT_SLIT.fits",
        ),
        (
            {"DPR_TYPE": "LAMP,FLAT", "DPR_TECH": "ECHELLE,PINHOLE"},
            False,
            f"{STAMP}_VIS_RO2_FLAT_ONEPIN.fits",
        ),
        (
            {"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,PINHOLE"},
            False,
            f"{STAMP}_VIS_RO2_ARC_ONEPIN.fits",
        ),
        (
            {"DPR_TYPE": "LAMP,FLAT,D,Q", "DPR_TECH": "ECHELLE,SLIT"},
            False,
            f"{STAMP}_VIS_RO2_QLAMP_FLAT_SLIT.fits",
        ),
        (
            {
                "DPR_TYPE": "OBJECT",
                "DPR_TECH": "ECHELLE,SLIT,NODDING",
                "OBJECT": "a-b",
            },
            False,
            f"{STAMP}_VIS_RO2_OBJECT_STARE_A_B_SLIT.fits",
        ),
    ),
)
def test_filenamer_pins_further_type_and_mask_combinations(
    log: object,
    headerOverrides: dict[str, object],
    prepared: bool,
    expected: str,
) -> None:
    frame = synthetic_ccd(prepared=prepared, headerOverrides=headerOverrides)

    assert _name(log, frame, settings=SOXS_SETTINGS) == expected


def test_filenamer_names_the_nir_arm(log: object) -> None:
    frame = synthetic_ccd(arm="NIR")

    assert _name(log, frame, settings=SOXS_SETTINGS) == f"{STAMP}_NIR_RO2_BIAS.fits"


@pytest.mark.parametrize(
    ("headerOverrides", "message"),
    (
        ({"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,SLIT"}, "mask/slit"),
        ({"DPR_TECH": "ECHELLE,SLIT"}, "mask/slit"),
        ({"DPR_TYPE": "OBJECT", "DPR_TECH": "ECHELLE,SLIT"}, "Frame type"),
        ({"DPR_TYPE": "OBJECT", "DPR_TECH": "IMAGE"}, "Frame type"),
        ({"DPR_TYPE": "BIAS,Q"}, "Frame type"),
    ),
)
def test_filenamer_rejects_header_combinations_it_cannot_name(
    log: object, headerOverrides: dict[str, object], message: str
) -> None:
    frame = synthetic_ccd(headerOverrides=headerOverrides)

    with pytest.raises(TypeError, match=message):
        _name(log, frame, settings=SOXS_SETTINGS)


def test_filenamer_logs_the_header_values_before_rejecting_an_unnameable_frame(
    log: object,
    capsys: pytest.CaptureFixture[str],
) -> None:
    frame = synthetic_ccd(
        headerOverrides={"DPR_TYPE": "LAMP,WAVE", "DPR_TECH": "ECHELLE,SLIT"}
    )

    with pytest.raises(TypeError):
        _name(log, frame, settings=SOXS_SETTINGS)

    printed = capsys.readouterr().out
    assert printed.endswith("\nlamp,wave\nechelle,slit\ncalib\n")
    assert "HIERARCH ESO DPR TYPE = 'LAMP,WAVE'" in printed
    assert (
        "error",
        "Frame mask/slit can't be determined - exiting",
    ) in log.messages


@pytest.mark.parametrize(
    ("deleted", "missing"),
    (
        ("ESO OBS ID", "ESO OBS ID"),
        ("ESO SEQ ARM", "ESO SEQ ARM"),
        ("ESO DPR TYPE", "ESO DPR TYPE"),
    ),
)
def test_filenamer_requires_the_obs_id_arm_and_type_keywords(
    log: object, deleted: str, missing: str
) -> None:
    """The observation ID is not in the name but must still be present."""
    frame = synthetic_ccd()
    del frame.header[deleted]

    with pytest.raises(KeyError, match=missing):
        _name(log, frame, settings=SOXS_SETTINGS)


def test_filenamer_looks_up_detector_parameters_only_when_none_are_supplied(
    log: object,
) -> None:
    frame = synthetic_ccd(headerOverrides={"SEQ_ARM": "ZZZ"})

    with pytest.raises(LookupError, match="detector 'ZZZ' cannot be found"):
        _name(log, frame, settings=SOXS_SETTINGS)

    keywordLookup = keyword_lookup(log=log, settings=SOXS_SETTINGS).get
    detectorLookup = detector_lookup(log=log, settings=SOXS_SETTINGS).get("VIS")

    filename = _name(
        log,
        frame,
        keywordLookup=keywordLookup,
        detectorLookup=detectorLookup,
        settings=False,
    )

    assert filename == f"{STAMP}_ZZZ_RO2_BIAS.fits"


@pytest.mark.parametrize("supplied", ("keywords", "detector"))
def test_filenamer_accepts_either_lookup_on_its_own(
    log: object, supplied: str
) -> None:
    frame = synthetic_ccd()
    lookups = {
        "keywords": {
            "keywordLookup": keyword_lookup(log=log, settings=SOXS_SETTINGS).get
        },
        "detector": {
            "detectorLookup": detector_lookup(log=log, settings=SOXS_SETTINGS).get(
                "VIS"
            )
        },
    }

    filename = _name(log, frame, settings=SOXS_SETTINGS, **lookups[supplied])

    assert filename == f"{STAMP}_VIS_RO2_BIAS.fits"


def test_filenamer_emits_only_its_own_start_and_end_debug_events(log: object) -> None:
    frame = synthetic_ccd()

    _name(log, frame, settings=SOXS_SETTINGS)

    filenamerEvents = [m for m in log.messages if "filenamer" in m[1]]
    assert filenamerEvents == [
        ("debug", "starting the ``filenamer`` function"),
        ("debug", "completed the ``filenamer`` function"),
    ]


@pytest.mark.parametrize(
    ("supplyDetector", "firstMissing"),
    ((False, "ESO SEQ ARM"), (True, "ESO OBS ID")),
)
def test_filenamer_reads_the_arm_before_the_obs_id_only_when_it_looks_up_the_detector(
    log: object, supplyDetector: bool, firstMissing: str
) -> None:
    frame = synthetic_ccd()
    del frame.header["ESO OBS ID"]
    del frame.header["ESO SEQ ARM"]
    lookups = {}
    if supplyDetector:
        lookups["detectorLookup"] = detector_lookup(
            log=log, settings=SOXS_SETTINGS
        ).get("VIS")

    with pytest.raises(KeyError, match=firstMissing):
        _name(log, frame, settings=SOXS_SETTINGS, **lookups)
