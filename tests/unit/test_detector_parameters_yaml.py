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
        # A MERGE KEY (<<) HAS NO CONSTRUCTOR AND MAY BE OVERRIDDEN LOCALLY, SO ONLY EXPLICIT KEYS ARE COUNTED
        if keyNode.tag == "tag:yaml.org,2002:merge":
            continue
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


def _column_names(calibrationsDir: Path, atlas: str) -> list[str]:
    with fits.open(calibrationsDir / atlas) as hdus:
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


def test_unique_key_loader_accepts_merge_keys_and_local_overrides() -> None:
    # ARRANGE
    text = "base: &base\n  gain: 1\n  ron: 2\nVIS:\n  <<: *base\n  ron: 3\n"

    # ACT
    parsed = _load_rejecting_duplicates(text)

    # ASSERT
    assert parsed["VIS"] == {"gain": 1, "ron": 3}


def test_unique_key_loader_rejects_a_repeated_key_beside_a_merge_key() -> None:
    # ARRANGE
    text = "base: &base\n  gain: 1\nVIS:\n  <<: *base\n  ron: 2\n  ron: 3\n"

    # ACT / ASSERT
    with pytest.raises(yaml.constructor.ConstructorError, match="duplicate key 'ron'"):
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
    reason="DY-1604: soxs_detector_parameters.yaml names BP_MAP_RP_VIS_1x2.fits, which is not packaged",
)
def test_soxs_vis_1x2_bad_pixel_map_is_packaged() -> None:
    # ARRANGE / ACT
    missing = _missing_static_calibrations("soxs")

    # ASSERT
    assert not (missing & KNOWN_MISSING)


_XSH_SOURCE_COLUMN = pytest.mark.xfail(
    strict=True,
    reason="DY-1605: xsh/ThAr.fits names its column 'Source', not 'source'",
)


@pytest.mark.parametrize(
    ("instrument", "arm", "atlas"),
    [
        ("soxs", "VIS", "thar_murphy.fits"),
        ("soxs", "NIR", "ArHgNeXe.fits"),
        pytest.param("xsh", "UVB", "ThAr.fits", marks=_XSH_SOURCE_COLUMN),
        pytest.param("xsh", "VIS", "ThAr.fits", marks=_XSH_SOURCE_COLUMN),
        ("xsh", "NIR", "ArHgNeXe.fits"),
    ],
)
def test_line_atlas_of_each_arm_has_the_columns_the_static_line_list_reads(
    instrument: str, arm: str, atlas: str
) -> None:
    # ARRANGE
    parameters = yaml.safe_load(_yaml_path(instrument).read_text())
    calibrationsDir = RESOURCES_DIR / "static_calibrations" / instrument

    # ACT (A MISSING ARM OR A MISSING line-atlas KEY RAISES KeyError AND FAILS THE TEST)
    namedAtlas = parameters[arm]["line-atlas"]
    columns = set(_column_names(calibrationsDir, namedAtlas))

    # ASSERT
    assert namedAtlas == atlas
    assert sorted(LINE_ATLAS_COLUMNS - columns) == []
