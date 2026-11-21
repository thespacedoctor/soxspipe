"""*tests for the loader that imports the scripts of tools/ by path*"""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from tests.unit import _tools_loader
from tests.unit._tools_loader import PREFIX, load_tool

pytestmark = pytest.mark.unit


def test_load_tool_registers_the_module_under_a_prefixed_name_and_reuses_it() -> None:
    # ARRANGE / ACT
    first = load_tool("measurement_common")
    second = load_tool("measurement_common")

    # ASSERT
    assert first is second
    assert sys.modules[PREFIX + "measurement_common"] is first


def test_load_tool_leaves_no_bare_helper_name_behind(monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    for key in [key for key in sys.modules if key.startswith(PREFIX)]:
        monkeypatch.delitem(sys.modules, key)
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


def test_a_load_that_raises_leaves_no_module_registered(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    (tmp_path / "broken_tool.py").write_text("raise RuntimeError('boom')\n")
    monkeypatch.setattr(_tools_loader, "TOOLS_DIR", tmp_path)

    # ACT
    with pytest.raises(RuntimeError, match="boom"):
        load_tool("broken_tool")

    # ASSERT
    assert PREFIX + "broken_tool" not in sys.modules
