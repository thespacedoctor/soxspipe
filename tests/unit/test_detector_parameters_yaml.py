"""Contracts for the packaged detector-parameters yaml files."""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.io import fits

import soxspipe
from soxspipe.commonutils.toolkit import read_spectral_format

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

# THE KEYS WHOSE FILES ARE LINE LISTS: THE ATLAS, THE ARC LINE TABLES AND THE SKY-LINE LIST
LINE_LIST_KEYS = ["line-atlas", "predicted pinhole lines", "skylines"]

# THE LINE LISTS THAT HOLD LAMP LINES, SO EVERY ORDER OF THEIR ARM HAS AT LEAST ONE LINE
LAMP_LINE_LIST_KEYS = ["line-atlas", "predicted pinhole lines"]

# THE NAMES THE LINE LISTS GIVE THEIR WAVELENGTH COLUMN (COMPARED IN LOWER CASE)
WAVELENGTH_COLUMN_NAMES = {"wave", "wavelength"}

# THE IVOA UCD OF AN AIR WAVELENGTH (THE ESO SDP CONVENTION, DY-2255)
AIR_WAVELENGTH_UCD = "em.wl;obs.atmos"

# STRONG Ar I LINES AS (AIR, VACUUM) PAIRS IN nm (DY-2255 AGENT NOTES)
AR_I_AIR_VACUUM_NM = [
    (696.5431, 696.7352),
    (706.7218, 706.9167),
    (750.3869, 750.5935),
    (763.5106, 763.7208),
    (772.3761, 772.5887),
    (794.8176, 795.0362),
    (801.4786, 801.6990),
    (811.5311, 811.7542),
    (826.4522, 826.6794),
    (842.4648, 842.6963),
    (912.2967, 912.5471),
]
LINE_MATCH_TOLERANCE_NM = 0.005

# THE WAVELENGTH SPAN (nm) OF THE X-SHOOTER UVB AND VIS ThAr ATLAS IN AIR (DY-2255 ACCEPTANCE CRITERIA)
XSH_THAR_ATLAS_MIN_NM = 290.0
XSH_THAR_ATLAS_MAX_NM = 1100.0

# THE ion LABELS OF NEUTRAL ARGON ROWS. A LIST WITH AN ion COLUMN IS COMPARED ON THESE ROWS ONLY, BECAUSE A DENSE
# ThAr LIST HOLDS OTHER LINES NEAR A VACUUM VALUE BY CHANCE (IN xsh/ThAr.fits A ThI LINE SITS AT 811.7557 nm, AIR)
AR_I_ION_LABELS = {"ArI", "Ar"}


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


@pytest.mark.parametrize(
    ("instrument", "arm", "atlas"),
    [
        ("soxs", "VIS", "thar_murphy.fits"),
        ("soxs", "NIR", "ArHgNeXe.fits"),
        ("xsh", "UVB", "ThAr.fits"),
        ("xsh", "VIS", "ThAr.fits"),
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


def _configured_line_lists(keys: list[str]) -> list[tuple[str, str, str, str]]:
    """Return (instrument, arm, key, file name) for every line list the yaml files name under ``keys``."""
    lineLists = []
    for instrument in INSTRUMENTS:
        parameters = yaml.safe_load(_yaml_path(instrument).read_text())
        for arm, armParameters in parameters.items():
            if not isinstance(armParameters, dict):
                continue
            for key in keys:
                if key in armParameters:
                    lineLists.extend(
                        (instrument, arm, key, name) for name in sorted(set(_file_names(armParameters[key])))
                    )
    return lineLists


def _line_list_id(lineList: tuple[str, str, str, str]) -> str:
    return "-".join(lineList)


def _line_list_path(instrument: str, name: str) -> Path:
    return RESOURCES_DIR / "static_calibrations" / instrument / name


def _wavelength_column_of(table: fits.BinTableHDU, name: str) -> tuple[int, np.ndarray]:
    columns = table.columns
    matches = [i for i, columnName in enumerate(columns.names) if columnName.lower() in WAVELENGTH_COLUMN_NAMES]
    assert len(matches) == 1, f"{name} has no single wavelength column: {columns.names}"
    index = matches[0]
    return index + 1, np.asarray(table.data[columns.names[index]], dtype=float)


def _wavelength_column(instrument: str, name: str) -> tuple[int, np.ndarray]:
    """Return the 1-based column number and the values of a line list's wavelength column."""
    with fits.open(_line_list_path(instrument, name)) as hdus:
        return _wavelength_column_of(hdus[1], name)


def _ar_i_wavelengths(instrument: str, name: str) -> np.ndarray:
    """Return the wavelengths of the Ar I rows of a line list, or of every row when the list has no ion column."""
    with fits.open(_line_list_path(instrument, name)) as hdus:
        _, waves = _wavelength_column_of(hdus[1], name)
        if "ion" not in hdus[1].columns.names:
            return waves
        ions = np.char.strip(np.asarray(hdus[1].data["ion"], dtype=str))
    return waves[np.isin(ions, list(AR_I_ION_LABELS))]


def _has_line_near(waves: np.ndarray, wavelength: float) -> bool:
    return bool(np.any(np.abs(waves - wavelength) < LINE_MATCH_TOLERANCE_NM))


@pytest.mark.parametrize("lineList", _configured_line_lists(LINE_LIST_KEYS), ids=_line_list_id)
def test_every_configured_line_list_marks_its_wavelength_column_as_air(lineList: tuple[str, str, str, str]) -> None:
    # ARRANGE
    instrument, _, _, name = lineList
    columnNumber, _ = _wavelength_column(instrument, name)

    # ACT
    header = fits.getheader(_line_list_path(instrument, name), 1)

    # ASSERT
    assert header.get(f"TUCD{columnNumber}") == AIR_WAVELENGTH_UCD


@pytest.mark.parametrize("lineList", _configured_line_lists(LINE_LIST_KEYS), ids=_line_list_id)
def test_no_configured_line_list_holds_vacuum_ar_i_wavelengths(lineList: tuple[str, str, str, str]) -> None:
    # ARRANGE
    instrument, _, _, name = lineList
    _, waves = _wavelength_column(instrument, name)
    argonWaves = _ar_i_wavelengths(instrument, name)
    inRange = [vacuum for _, vacuum in AR_I_AIR_VACUUM_NM if waves.min() <= vacuum <= waves.max()]

    # ACT
    vacuumMatches = [vacuum for vacuum in inRange if _has_line_near(argonWaves, vacuum)]

    # ASSERT
    assert vacuumMatches == []


@pytest.mark.parametrize("lineList", _configured_line_lists(LAMP_LINE_LIST_KEYS), ids=_line_list_id)
def test_every_order_of_the_arm_holds_a_line_of_each_lamp_line_list(lineList: tuple[str, str, str, str]) -> None:
    """A list in the wrong unit (for example Angstrom) puts no line in any order, which are in nm."""
    # ARRANGE
    instrument, arm, _, name = lineList
    _, waves = _wavelength_column(instrument, name)
    orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
        log=logging.getLogger(__name__), settings={"instrument": instrument}, arm=arm
    )

    # ACT
    emptyOrders = [
        int(order)
        for order, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax, strict=True)
        if not np.any((waves >= wmin) & (waves <= wmax))
    ]

    # ASSERT
    assert emptyOrders == []


@pytest.mark.parametrize("arm", ["UVB", "VIS"])
def test_xsh_uvb_and_vis_line_atlas_holds_air_wavelengths_in_nm(arm: str) -> None:
    # ARRANGE
    parameters = yaml.safe_load(_yaml_path("xsh").read_text())
    atlas = parameters[arm]["line-atlas"]

    # ACT
    _, waves = _wavelength_column("xsh", atlas)
    argonWaves = _ar_i_wavelengths("xsh", atlas)

    # ASSERT
    assert waves.min() >= XSH_THAR_ATLAS_MIN_NM
    assert waves.max() <= XSH_THAR_ATLAS_MAX_NM
    assert [air for air, _ in AR_I_AIR_VACUUM_NM if not _has_line_near(argonWaves, air)] == []
    assert [vacuum for _, vacuum in AR_I_AIR_VACUUM_NM if _has_line_near(argonWaves, vacuum)] == []
