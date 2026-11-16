"""The SOXS NIR nod SOFs list the master dark as an optional calibration, and nothing else changes (DY-2194)."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

import soxspipe
from soxspipe.commonutils.calibration_matching import CALIBRATION_RULES

pytestmark = pytest.mark.unit

RESOURCES = Path(soxspipe.__file__).parent / "resources"
NOD_ENTRIES = ("science_nodding", "flux_standard_nodding", "telluric_standard_nodding")


def _sof_map(name: str) -> dict:
    return yaml.safe_load((RESOURCES / name).read_text())


@pytest.mark.parametrize("entry", NOD_ENTRIES)
def test_a_nir_nod_sof_lists_the_master_dark_as_optional(entry: str) -> None:
    # ACT
    filters = _sof_map("soxs_sof_map.yaml")[entry]

    # ASSERT
    assert filters["optional calibrations"] == {"NIR": ["mdark"]}


@pytest.mark.parametrize("entry", NOD_ENTRIES)
def test_the_master_dark_is_never_a_required_calibration_of_a_nod_sof(entry: str) -> None:
    # ACT
    required = _sof_map("soxs_sof_map.yaml")[entry]["calibrations"]

    # ASSERT
    assert all("mdark" not in types for types in required.values())
    assert set(required) == {"VIS", "NIR"}


def test_only_the_nir_nod_sofs_have_optional_calibrations() -> None:
    # ACT
    withOptional = sorted(
        name for name, filters in _sof_map("soxs_sof_map.yaml").items() if filters.get("optional calibrations")
    )

    # ASSERT
    assert withOptional == sorted(NOD_ENTRIES)


def test_the_xshooter_sof_map_has_no_optional_calibrations() -> None:
    # ACT
    entries = [name for name, filters in _sof_map("xsh_sof_map.yaml").items() if "optional calibrations" in filters]

    # ASSERT
    assert entries == []


def test_every_optional_calibration_type_has_a_matching_rule() -> None:
    # ACT
    optionalTypes = {
        calibrationType
        for filters in _sof_map("soxs_sof_map.yaml").values()
        for types in (filters.get("optional calibrations") or {}).values()
        for calibrationType in types
    }

    # ASSERT
    assert optionalTypes == {"mdark"}
    assert optionalTypes <= set(CALIBRATION_RULES)
