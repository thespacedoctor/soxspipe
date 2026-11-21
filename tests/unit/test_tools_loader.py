"""*tests for the loader that imports the scripts of tools/ by path*"""

from __future__ import annotations

import sys

import pytest

from tests.unit._tools_loader import PREFIX, load_tool

pytestmark = pytest.mark.unit


def test_load_tool_registers_the_module_under_a_prefixed_name_and_reuses_it() -> None:
    # ARRANGE / ACT
    first = load_tool("measurement_common")
    second = load_tool("measurement_common")

    # ASSERT
    assert first is second
    assert sys.modules[PREFIX + "measurement_common"] is first


def test_load_tool_leaves_no_bare_helper_name_behind() -> None:
    # ARRANGE
    before = {name: sys.modules.get(name) for name in ("measurement_common", "line_fitting", "spec_res_frames")}

    # ACT
    load_tool("measure_spec_res")

    # ASSERT
    assert {name: sys.modules.get(name) for name in before} == before


def test_a_script_shares_the_helper_module_object_the_tests_hold() -> None:
    # ARRANGE / ACT
    script = load_tool("measure_log_wavelength_step")

    # ASSERT
    assert script.common is load_tool("measurement_common")
