"""Contracts for the packaged detector-parameters yaml files."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml
from astropy.io import fits

import soxspipe

pytestmark = pytest.mark.unit

RESOURCES_DIR = Path(soxspipe.__file__).parent / "resources"
INSTRUMENTS = ["soxs", "xsh"]

# THE KEYS WHOSE VALUES ARE JOINED ONTO THE STATIC-CALIBRATIONS FOLDER BY THE CODE THAT READS THEM:
# base_recipe ("bad-pixel map"), create_dispersion_map ("predicted pinhole lines", "line-atlas"),
# toolkit ("spectral format table", "skylines"), response_function and the soxs_stare/nod/offset
# recipes ("extinction", "flux-standards")
STATIC_CALIBRATION_KEYS = [
    "bad-pixel map",
    "predicted pinhole lines",
    "spectral format table",
    "line-atlas",
    "skylines",
    "extinction",
    "flux-standards",
]

# NAMED GAPS THAT NEED A DECISION FROM DAVE (NEITHER THE FILE NOR THE KEY CAN BE INVENTED HERE):
# THE SOXS VIS 1x2 BAD-PIXEL MAP IS NAMED BUT NOT PACKAGED
KNOWN_MISSING = {("soxs", "VIS", "bad-pixel map", "BP_MAP_RP_VIS_1x2.fits")}

# THE COLUMNS create_new_static_line_list READS FROM THE LINE ATLAS
LINE_ATLAS_COLUMNS = {"wave", "ion", "source", "amplitude"}


class _UniqueKeyLoader(yaml.SafeLoader):
    """A safe yaml loader that raises when a mapping repeats a key."""


def _construct_unique_mapping(loader: _UniqueKeyLoader, node: yaml.MappingNode, deep: bool = False) -> dict:
    seenKeys = set()
    for keyNode, _ in node.value:
        key = loader.construct_object(keyNode, deep=deep)
        if key in seenKeys:
            raise yaml.constructor.ConstructorError(None, None, f"duplicate key {key!r}", keyNode.start_mark)
        seenKeys.add(key)
    return yaml.SafeLoader.construct_mapping(loader, node, deep=deep)


_UniqueKeyLoader.add_constructor(yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, _construct_unique_mapping)


def _load_rejecting_duplicates(text: str):
    loader = _UniqueKeyLoader(text)
    try:
        return loader.get_single_data()
    finally:
        loader.dispose()


def _column_names(calibrationsDir: Path, armParameters: dict) -> list[str]:
    with fits.open(calibrationsDir / armParameters["line-atlas"]) as hdus:
        return list(hdus[1].columns.names)


def _yaml_path(instrument: str) -> Path:
    return RESOURCES_DIR / f"{instrument}_detector_parameters.yaml"


def _file_names(value) -> list[str]:
    """Collect every file name held by a key's value, however deeply it is nested."""
    if isinstance(value, dict):
        return [name for child in value.values() for name in _file_names(child)]
    return [value]


@pytest.mark.parametrize("instrument", INSTRUMENTS)
def test_detector_parameters_yaml_has_no_duplicate_keys(instrument: str) -> None:
    # ARRANGE
    text = _yaml_path(instrument).read_text()

    # ACT
    parsed = _load_rejecting_duplicates(text)

    # ASSERT
    assert isinstance(parsed, dict)


def test_unique_key_loader_rejects_a_repeated_key() -> None:
    # ARRANGE
    text = "VIS:\n  line-atlas: a.fits\n  line-atlas: b.fits\n"

    # ACT / ASSERT
    with pytest.raises(yaml.constructor.ConstructorError, match="duplicate key 'line-atlas'"):
        _load_rejecting_duplicates(text)


def _missing_static_calibrations(instrument: str) -> set[tuple[str, str, str, str]]:
    parameters = yaml.safe_load(_yaml_path(instrument).read_text())
    calibrationsDir = RESOURCES_DIR / "static_calibrations" / instrument
    named = {
        (instrument, arm, key, name)
        for arm, armParameters in parameters.items()
        for key in STATIC_CALIBRATION_KEYS
        if key in armParameters
        for name in _file_names(armParameters[key])
    }
    assert named, "no static calibration names were found, so the key list is out of date"
    return {item for item in named if not (calibrationsDir / item[3]).is_file()}


@pytest.mark.parametrize("instrument", INSTRUMENTS)
def test_every_static_calibration_named_in_detector_parameters_exists(instrument: str) -> None:
    # ARRANGE / ACT
    missing = _missing_static_calibrations(instrument) - KNOWN_MISSING

    # ASSERT
    assert sorted(missing) == []


@pytest.mark.xfail(
    strict=True,
    reason="soxs_detector_parameters.yaml names BP_MAP_RP_VIS_1x2.fits, which is not packaged; needs a decision",
)
def test_soxs_vis_1x2_bad_pixel_map_is_packaged() -> None:
    # ARRANGE / ACT
    missing = _missing_static_calibrations("soxs")

    # ASSERT
    assert not (missing & KNOWN_MISSING)


@pytest.mark.parametrize(
    ("arm", "atlas"),
    [("VIS", "thar_murphy.fits"), ("NIR", "ArHgNeXe.fits")],
)
def test_soxs_line_atlas_names_the_agreed_file(arm: str, atlas: str) -> None:
    # ARRANGE
    parameters = yaml.safe_load(_yaml_path("soxs").read_text())

    # ACT / ASSERT
    assert parameters[arm]["line-atlas"] == atlas


@pytest.mark.parametrize(
    "instrument",
    [
        "soxs",
        pytest.param(
            "xsh",
            marks=pytest.mark.xfail(
                strict=True,
                reason="xsh/ThAr.fits names its column 'Source', not 'source'; X-Shooter values need a decision",
            ),
        ),
    ],
)
def test_every_line_atlas_has_the_columns_the_static_line_list_reads(instrument: str) -> None:
    # ARRANGE
    parameters = yaml.safe_load(_yaml_path(instrument).read_text())
    calibrationsDir = RESOURCES_DIR / "static_calibrations" / instrument

    # ACT
    lacking = {
        (arm, armParameters["line-atlas"]): sorted(
            LINE_ATLAS_COLUMNS - set(_column_names(calibrationsDir, armParameters))
        )
        for arm, armParameters in parameters.items()
        if "line-atlas" in armParameters
    }
    lacking = {atlas: columns for atlas, columns in lacking.items() if columns}

    # ASSERT
    assert lacking == {}
