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

import ast
import importlib
import logging
import warnings
from collections import OrderedDict
from collections.abc import Iterator
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning
from astropy.nddata import CCDData
from astropy.utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning
from astropy.wcs import FITSFixedWarning
from photutils.utils import NoDetectionsWarning

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.create_dispersion_map import create_dispersion_map, measure_line_position
from soxspipe.commonutils.data_organiser import data_organiser
from soxspipe.commonutils.fits_headers import STANDARD_KEYWORD, hierarch, read_ccd, set_header_card
from soxspipe.commonutils.phase3 import add_qc_to_header, sort_keywords
from soxspipe.commonutils.set_of_files import ImageFileCollection
from soxspipe.commonutils.specutils_calls import fit_generic_continuum_quietly, spectrum1d_quietly
from soxspipe.commonutils.toolkit import save_qc_plot
from soxspipe.recipes.base_recipe import base_recipe
from tests.factories import instrument_header, pipeline_settings, raw_fits
from tests.factories.keyword_recipes import keyword_recipe, raw_frame_table_with_calibration
from tests.factories.mflat_verification import VerificationInventory, bare_mflat_recipe, stub_mflat_basics
from tests.factories.prepared_frames import bad_pixel_map_file, prepare_frame_recipe, raw_frame_file

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
    bad_pixel_map_file(tmp_path / "calibrations")
    recipe = prepare_frame_recipe(log, tmp_path)
    rawPath = raw_frame_file(tmp_path)
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
    bad_pixel_map_file(tmp_path / "calibrations")
    recipe = prepare_frame_recipe(log, tmp_path)
    rawPath = raw_frame_file(tmp_path)
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
    stub_mflat_basics(
        monkeypatch,
        imageTypes=["LAMP,FLAT"],
        imageTech=["ECHELLE,SLIT"],
        imageCategories=["MASTER_BIAS_VIS", "ORDER_TAB_VIS"],
        arm="VIS",
    )
    recipe = bare_mflat_recipe(log, VerificationInventory())
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


def _raised_while(call: Any) -> list[warnings.WarningMessage]:
    """Run ``call`` with every warning switched on and return the warnings it raised."""
    with warnings.catch_warnings(record=True) as raised:
        warnings.simplefilter("always")
        call()
    return raised


def _of_category(raised: list[warnings.WarningMessage], category: type[Warning]) -> list[warnings.WarningMessage]:
    return [warning for warning in raised if issubclass(warning.category, category)]


def test_sorting_a_header_of_eso_keywords_raises_no_hierarch_notice(log: Any) -> None:
    """Re-adding long ESO keywords to a header names them `HIERARCH`, so astropy has nothing to warn about."""
    # ARRANGE
    header = fits.Header()
    header["SIMPLE"] = True
    header["HIERARCH ESO DET BINX"] = (2, "binning in x")
    header["HIERARCH ESO DPR TYPE"] = ("BIAS", "frame type")

    # ACT
    raised = _raised_while(lambda: sort_keywords(log=log, header=header))

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert list(header.keys()) == ["ESO DET BINX", "ESO DPR TYPE", "SIMPLE"]
    assert header["ESO DET BINX"] == 2
    assert header.comments["ESO DPR TYPE"] == "frame type"


def test_writing_qc_cards_to_a_header_raises_no_hierarch_notice() -> None:
    """A QC value reaches the header as an `ESO QC <NAME>` card without a `VerifyWarning`."""
    # ARRANGE
    header = fits.Header()
    qc = pd.DataFrame(
        {"qc_name": ["RON"], "qc_value": [3.5], "qc_comment": ["read noise"], "to_header": [True]},
    )

    # ACT
    raised = _raised_while(lambda: add_qc_to_header(header, qc))

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert header["ESO QC RON"] == 3.5


def test_stamping_a_product_header_raises_no_hierarch_notice(log: Any, tmp_path: Path) -> None:
    """`update_fits_keywords` writes its ESO PRO records without a `VerifyWarning`."""
    # ARRANGE
    calibrationPath = tmp_path / "master_bias.fits"
    calibrationPath.write_bytes(b"synthetic calibration payload")
    recipe = keyword_recipe(tmp_path, log, raw_frame_table_with_calibration(tmp_path, calibrationPath))
    recipe.recipeSettings = {"stacked-clipping-sigma": 3.0, "clipping": {"iterations": 2}}
    frame = CCDData(np.zeros((2, 2)), unit=u.electron, meta=instrument_header())

    # ACT
    raised = _raised_while(lambda: recipe.update_fits_keywords(frame))

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert frame.header["ESO PRO REC1 RAW1 NAME"] == "raw_one.fits"
    assert frame.header["ESO PRO REC1 CAL1 CATG"] == "MASTER_BIAS_VIS"
    assert frame.header["ESO PRO REC1 PARAM1 VALUE"] == 3.0


def test_the_dispersion_map_header_raises_no_hierarch_notice(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """`write_map_to_file` builds its ESO PRO header without a `VerifyWarning`."""
    # ARRANGE
    settings = pipeline_settings(tmp_path)
    kw = keyword_lookup(log=log, settings=settings).get
    mapper = object.__new__(create_dispersion_map)
    mapper.log = log
    mapper.arm = "VIS"
    mapper.kw = kw
    mapper.sofName = "SYNTHETIC"
    mapper.pinholeFrame = SimpleNamespace(
        header=fits.Header({kw("DPR_TECH"): "ARC", kw("DPR_CATG"): "CALIB", kw("DPR_TYPE"): "LAMP"})
    )
    mapper.productDir = str(tmp_path)
    mapper.settings = settings
    mapper.qc = pd.DataFrame(columns=["qc_name", "qc_value", "qc_unit", "qc_comment"])
    written: dict[str, Any] = {}
    phase3Module = importlib.import_module("soxspipe.commonutils.phase3")
    monkeypatch.setattr(phase3Module, "write_fits_table_to_disk", lambda **kwargs: written.update(kwargs))

    # ACT
    raised = _raised_while(
        lambda: mapper.write_map_to_file(xcoeff=[1.0, 2.0], ycoeff=[3.0, 4.0], orderDeg=1, wavelengthDeg=0, slitDeg=0)
    )

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert written["header"][kw("PRO_CATG")] == "DISP_TAB_VIS"
    assert written["header"][kw("PRO_TECH")] == "ECHELLE,PINHOLE"


def test_a_filter_that_matches_no_frame_returns_an_empty_collection_without_a_warning(tmp_path: Path) -> None:
    """An empty result is how the recipes probe for an optional calibration, so ccdproc's notice is dropped."""
    # ARRANGE
    framePath = raw_fits(tmp_path / "frame.fits", headerOverrides={"DPR_TYPE": "BIAS"})
    collection = ImageFileCollection(filenames=[str(framePath)], keywords=["*"])
    filtersBefore = _filters_snapshot()

    # ACT
    raised = _raised_while(lambda: collection.filter(**{"ESO DPR TYPE": "NO SUCH TYPE"}))

    # ASSERT
    assert _of_category(raised, AstropyUserWarning) == []
    assert collection.filter(**{"ESO DPR TYPE": "NO SUCH TYPE"}).files == []
    assert warnings.filters == filtersBefore


def test_reading_a_frame_with_a_deprecated_wcs_keyword_raises_no_wcs_fix_notice(tmp_path: Path) -> None:
    """astropy's note about a raw header's old `RADECSYS` card is not for the user, who cannot change the raw frame."""
    # ARRANGE
    framePath = raw_fits(
        tmp_path / "frame.fits",
        headerOverrides={
            "RADECSYS": "FK5",
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRVAL1": 10.0,
            "CRVAL2": 20.0,
            "MJD-OBS": 60311.1278436,
        },
    )

    # ACT
    raised = _raised_while(lambda: read_ccd(str(framePath), hdu=0, unit=u.adu))

    # ASSERT
    assert _of_category(raised, FITSFixedWarning) == []


@pytest.mark.parametrize(
    ("keyword", "expected"),
    [
        ("EXPTIME", "EXPTIME"),
        ("exptime", "exptime"),
        ("DATE-OBS", "DATE-OBS"),
        ("ESO DET BINX", "HIERARCH ESO DET BINX"),
        ("EXPOSURE1", "HIERARCH EXPOSURE1"),
        ("HIERARCH ESO DET BINX", "HIERARCH ESO DET BINX"),
        ("DP1.AXIS.1", "DP1.AXIS.1"),
    ],
)
def test_a_keyword_is_prefixed_with_hierarch_only_when_fits_needs_it(keyword: str, expected: str) -> None:
    """A standard keyword stays as it is; a long or spaced one gets the prefix once."""
    # ACT
    result = hierarch(keyword)

    # ASSERT
    assert result == expected


def test_building_a_spectrum_raises_no_gwcs_deprecation_notice() -> None:
    """gwcs' use of a deprecated astropy function is ignored for the `Spectrum1D` call only."""
    # ARRANGE
    wavelength = np.linspace(500.0, 510.0, 20) * u.nm
    filtersBefore = _filters_snapshot()

    # ACT
    raised = _raised_while(lambda: spectrum1d_quietly(flux=np.ones(20) * u.electron, spectral_axis=wavelength))

    # ASSERT
    assert _of_category(raised, AstropyDeprecationWarning) == []
    assert warnings.filters == filtersBefore


def test_a_continuum_fit_raises_no_linear_model_notice_and_returns_the_continuum() -> None:
    """The fitter note is ignored for the continuum call only, and the fit still comes back."""
    # ARRANGE
    wavelength = np.linspace(500.0, 510.0, 200)
    spectrum = spectrum1d_quietly(flux=(2.0 + 0.01 * wavelength) * u.Unit(""), spectral_axis=wavelength * u.nm)
    fitted: list[Any] = []

    # ACT
    raised = _raised_while(lambda: fitted.append(fit_generic_continuum_quietly(spectrum)))

    # ASSERT
    assert _of_category(raised, AstropyUserWarning) == []
    continuum = fitted[0](wavelength * u.nm)
    assert continuum.value == pytest.approx(2.0 + 0.01 * wavelength, rel=1e-3)


def test_a_header_card_set_through_the_helper_is_the_card_astropy_makes_for_hierarch() -> None:
    """`set_header_card` writes the same bytes as astropy's own `HIERARCH` card, with no warning."""
    # ARRANGE
    header = fits.Header()
    expected = fits.Header()
    expected["HIERARCH ESO DET BINX"] = (2, "binning in x")
    expected["EXPTIME"] = 30.0

    # ACT
    raised = _raised_while(
        lambda: (
            set_header_card(header, "ESO DET BINX", 2, "binning in x"),
            set_header_card(header, "EXPTIME", 30.0),
        )
    )

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert header.tostring() == expected.tostring()


def test_a_record_valued_keyword_is_written_as_astropy_writes_it() -> None:
    """`DP1.AXIS.1` stays a record-valued card, not a `HIERARCH` one."""
    # ARRANGE
    header = fits.Header()
    expected = fits.Header()
    expected["DP1.AXIS.1"] = 1

    # ACT
    set_header_card(header, "DP1.AXIS.1", 1)

    # ASSERT
    assert header.tostring() == expected.tostring()


def test_a_dict_header_gets_the_hierarch_keyword_so_astropy_raises_no_notice_when_it_becomes_a_header() -> None:
    """A CCDData with a plain dict header is turned into a `fits.Header` on write; the long key must not warn then."""
    # ARRANGE
    frame = CCDData(np.zeros((2, 2)), unit=u.adu, meta=OrderedDict())

    # ACT
    set_header_card(frame.header, "ESO PRO CATG", "MASTER_BIAS_VIS")
    raised = _raised_while(lambda: fits.Header(frame.meta))

    # ASSERT
    assert _of_category(raised, VerifyWarning) == []
    assert fits.Header(frame.meta)["ESO PRO CATG"] == "MASTER_BIAS_VIS"
    assert list(frame.meta) == ["HIERARCH ESO PRO CATG"]


def test_reading_a_frame_still_shows_a_wcs_fix_that_is_not_one_of_the_known_notes(tmp_path: Path) -> None:
    """Only the two known WCS notes are ignored; a date fix on the same read is still raised."""
    # ARRANGE
    framePath = raw_fits(
        tmp_path / "frame.fits",
        headerOverrides={
            "CTYPE1": "RA---TAN",
            "CTYPE2": "DEC--TAN",
            "CRVAL1": 10.0,
            "CRVAL2": 20.0,
            "DATE-OBS": "2024-01-02T03:04:05",
            "MJD-OBS": 1.0,
        },
    )

    # ACT
    raised = _raised_while(lambda: read_ccd(str(framePath), hdu=0, unit=u.adu))

    # ASSERT
    messages = [str(warning.message) for warning in _of_category(raised, FITSFixedWarning)]
    assert any("MJD-OBS and DATE-OBS are inconsistent" in message for message in messages)


def test_saving_a_qc_plot_with_collapsed_axes_raises_no_layout_notice(tmp_path: Path) -> None:
    """matplotlib still writes the figure when it cannot apply the constrained layout, so its notice is dropped."""
    # ARRANGE
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(0.3, 0.3), layout="constrained")
    axes = figure.add_subplot(111)
    axes.set_title("a long title")
    axes.set_xlabel("x label")
    axes.set_ylabel("y label")
    plotPath = tmp_path / "plot.pdf"

    # ACT
    try:
        raised = _raised_while(lambda: save_qc_plot(str(plotPath)))
    finally:
        plt.close(figure)

    # ASSERT
    assert [warning for warning in raised if "constrained_layout" in str(warning.message)] == []
    assert plotPath.stat().st_size > 0


PACKAGE_PATH = Path(__file__).resolve().parents[2] / "soxspipe"


def _header_names(node: ast.expr) -> bool:
    """True for ``header`` and for ``<anything>.header``."""
    return (isinstance(node, ast.Name) and node.id == "header") or (
        isinstance(node, ast.Attribute) and node.attr == "header"
    )


def _is_standard_literal(key: ast.expr) -> bool:
    """True for a string literal that FITS stores as a standard card, such as ``"EXTNAME"``."""
    return isinstance(key, ast.Constant) and isinstance(key.value, str) and bool(STANDARD_KEYWORD.fullmatch(key.value))


def test_every_header_card_the_package_writes_goes_through_the_helper() -> None:
    """Only a standard literal keyword is written by subscript, so a new ESO keyword cannot skip `hierarch`."""
    # ARRANGE
    offenders = []

    # ACT
    for sourcePath in sorted(PACKAGE_PATH.rglob("*.py")):
        if "tests" in sourcePath.relative_to(PACKAGE_PATH).parts:
            continue
        tree = ast.parse(sourcePath.read_text())
        helperLines = {
            line
            for function in ast.walk(tree)
            if isinstance(function, ast.FunctionDef) and function.name == "set_header_card"
            for line in range(function.lineno, function.end_lineno + 1)
        }
        for node in ast.walk(tree):
            if getattr(node, "lineno", None) in helperLines:
                continue
            targets = node.targets if isinstance(node, ast.Assign) else [getattr(node, "target", None)]
            for target in targets:
                if (
                    isinstance(target, ast.Subscript)
                    and _header_names(target.value)
                    and not _is_standard_literal(target.slice)
                ):
                    offenders.append(f"{sourcePath.relative_to(PACKAGE_PATH)}:{node.lineno}")

    # ASSERT
    assert offenders == []


@pytest.mark.parametrize("recipeFile", ["soxs_nod.py", "soxs_offset.py"])
def test_the_nodding_recipes_read_their_frames_through_read_ccd(recipeFile: str) -> None:
    """`CCDData.read` called directly would show astropy's WCS notes that `read_ccd` ignores."""
    # ARRANGE
    tree = ast.parse((PACKAGE_PATH / "recipes" / recipeFile).read_text())

    # ACT
    directReads = [
        node.lineno
        for node in ast.walk(tree)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "read"
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "CCDData"
    ]

    # ASSERT
    assert directReads == []
