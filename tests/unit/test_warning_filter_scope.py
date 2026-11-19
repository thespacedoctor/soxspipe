"""No soxspipe code path may change the process-wide warning filters (DY-1480).

A recipe, the data organiser and the dispersion-map builder each used to install a
filter that stayed in force for the rest of the process, so every later warning
(numpy `RuntimeWarning`, pandas, astropy, `DeprecationWarning`) was dropped before
`logging.captureWarnings` could send it to the recipe log. The contract is that
`warnings.filters` is the same after each of these entry points returns, and that a
warning raised inside a recipe step reaches the log.

pytest wraps every test in `warnings.catch_warnings`, so a filter that leaks inside a
test cannot reach the next one; comparing the filter list before and after the call
is what exposes it.
"""

from __future__ import annotations

import importlib
import logging
import warnings
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData
from photutils.utils import NoDetectionsWarning

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map, measure_line_position
from soxspipe.commonutils.data_organiser import data_organiser
from soxspipe.recipes.base_recipe import base_recipe
from tests.unit.test_base_recipe_prepare_frame_characterization import _bad_pixel_map, _raw_frame, _recipe
from tests.unit.test_soxs_mflat_verification_characterization import (
    _bare_recipe as _bare_mflat_recipe,
)
from tests.unit.test_soxs_mflat_verification_characterization import (
    _stub_basics as _stub_mflat_basics,
)
from tests.unit.test_soxs_mflat_verification_characterization import (
    _VerificationInventory,
)

pytestmark = pytest.mark.unit


@pytest.fixture(autouse=True)
def no_warning_capture_left_over() -> Iterator[None]:
    """Switch `logging.captureWarnings` off before and after each test.

    `logging.captureWarnings(True)` does nothing when an earlier test left it switched on, and a test that
    reaches `_prepare_single_frame` switches it on, so every test in this file starts and ends with it off.
    """
    logging.captureWarnings(False)
    yield
    logging.captureWarnings(False)


@pytest.fixture
def warnings_logged(caplog: pytest.LogCaptureFixture) -> pytest.LogCaptureFixture:
    """Capture what `logging.captureWarnings` sends to the log."""
    caplog.set_level(logging.WARNING)
    return caplog


def _filters_snapshot() -> list[Any]:
    """Return a copy of the process-wide warning filters."""
    return list(warnings.filters)


def test_preparing_a_frame_leaves_the_warning_filters_unchanged(log: Any, tmp_path: Path) -> None:
    """`_prepare_single_frame` installs no process-wide filter."""
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)
    filtersBefore = _filters_snapshot()

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    assert warnings.filters == filtersBefore


def test_a_runtime_warning_during_frame_preparation_reaches_the_log(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    warnings_logged: pytest.LogCaptureFixture,
) -> None:
    """A numpy `RuntimeWarning` raised inside a recipe step is logged, not dropped.

    `logging.captureWarnings` sends it to the `py.warnings` logger, which propagates to the root logger. That is where
    the console and `soxspipe.log` handlers of the `logging settings` block sit, and where `caplog` listens. The
    recipe's own `.log` file does not take it: that handler excludes WARNING and above by design.
    """
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _recipe(log, tmp_path)
    rawPath = _raw_frame(tmp_path)
    badPixelMask = base_recipe._bad_pixel_mask

    def noisy_bad_pixel_mask(self: base_recipe, frame: Any) -> Any:
        np.sqrt(np.array([-1.0]))
        return badPixelMask(self, frame)

    monkeypatch.setattr(base_recipe, "_bad_pixel_mask", noisy_bad_pixel_mask)

    # ACT
    recipe._prepare_single_frame(rawPath)

    # ASSERT
    logged = [record.getMessage() for record in warnings_logged.records if record.name == "py.warnings"]
    assert any("RuntimeWarning" in message and "invalid value encountered in sqrt" in message for message in logged)


def test_verifying_master_flat_inputs_leaves_the_warning_filters_unchanged(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`soxs_mflat.verify_input_frames` installs no process-wide filter."""
    # ARRANGE
    _stub_mflat_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = _bare_mflat_recipe(log, _VerificationInventory())
    filtersBefore = _filters_snapshot()

    # ACT
    recipe.verify_input_frames()

    # ASSERT
    assert warnings.filters == filtersBefore


def test_creating_the_data_organiser_leaves_the_warning_filters_unchanged(log: Any, tmp_path: Path) -> None:
    """`data_organiser.__init__` installs no process-wide filter."""
    # ARRANGE
    filtersBefore = _filters_snapshot()

    # ACT
    data_organiser(log=log, rootDir=str(tmp_path), dbConnect=False)

    # ASSERT
    assert warnings.filters == filtersBefore


def test_creating_the_dispersion_map_builder_leaves_the_warning_filters_unchanged(
    log: Any,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`create_dispersion_map.__init__` installs no process-wide filter."""
    # ARRANGE
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    frame = SimpleNamespace(
        header={"SEQ_ARM": "VIS", "DATE_OBS": "2025-01-02T03:04:05", "INSTRUME": "SOXS", "EXPTIME": 30.0}
    )
    monkeypatch.setattr(module, "keyword_lookup", lambda **_: SimpleNamespace(get=lambda name: name))
    monkeypatch.setattr(toolkitModule, "get_calibration_lamp", lambda **_: "Ar")
    monkeypatch.setattr(
        commonutilsModule, "detector_lookup", lambda **_: SimpleNamespace(get=lambda arm: {"dispersion-axis": "x"})
    )
    monkeypatch.setattr(toolkitModule, "utility_setup", lambda **_: ("qc-dir", "product-dir"))
    filtersBefore = _filters_snapshot()

    # ACT
    create_dispersion_map(
        log=log,
        settings={"instrument": "soxs"},
        recipeSettings={"order-deg": 2},
        pinholeFrame=frame,
        qcTable=pd.DataFrame(),
        productsTable=pd.DataFrame(),
        sofName="synthetic.sof",
        startNightDate="2025-01-02",
        turnOffMP=True,
    )

    # ASSERT
    assert warnings.filters == filtersBefore


def test_a_stamp_with_no_sources_raises_no_photutils_warning_and_leaves_the_filters_unchanged(log: Any) -> None:
    """An empty stamp is routine, so photutils' `NoDetectionsWarning` is suppressed for that call only."""
    # ARRANGE
    noise = np.random.default_rng(1).normal(100.0, 1.0, (11, 11))
    stamp = CCDData(noise, unit=u.adu, mask=np.zeros_like(noise, dtype=bool))
    filtersBefore = _filters_snapshot()

    # ACT
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        filtersInsideProbe = _filters_snapshot()
        detected = measure_line_position(
            (stamp, 10, 21, 20, 31),
            log=log,
            windowHalf=5,
            iraf=False,
            sigmaLimit=3,
            iteration=0,
            returnAll=True,
        )
        filtersAfterCall = _filters_snapshot()

    # ASSERT
    assert detected == []
    assert [warning for warning in raised if issubclass(warning.category, NoDetectionsWarning)] == []
    assert filtersAfterCall == filtersInsideProbe
    assert warnings.filters == filtersBefore
