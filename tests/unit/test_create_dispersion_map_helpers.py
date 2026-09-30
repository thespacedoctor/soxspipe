"""Focused contracts for dispersion-map setup and table helpers."""

from __future__ import annotations

import copy
import importlib
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from soxspipe.commonutils.create_dispersion_map import (
    create_dispersion_map,
    find_largest_cluster_center,
    straighten_mph_sets,
)

pytestmark = pytest.mark.unit


def _mapper(log: object) -> create_dispersion_map:
    """Return a minimal mutable mapper for focused method contracts."""
    mapper = create_dispersion_map.__new__(create_dispersion_map)
    mapper.log = log
    mapper.recipeSettings = {
        "order-deg": 2,
        "wavelength-deg": 3,
        "slit-deg": 1,
        "pixel-window-size": 6,
        "pinhole-detection-thres-sigma": 4,
    }
    mapper.settings = {"bootstrap_dispersion_solution": False}
    mapper.firstGuessMap = False
    mapper.debug = False
    mapper.qc = pd.DataFrame()
    mapper.products = pd.DataFrame()
    mapper.recipeName = "soxs-disp-solution"
    mapper.dateObs = "2024-01-01T00:00:00"
    mapper.arm = "VIS"
    mapper.inst = "SOXS"
    return mapper


def _line_table(*, includeIon: bool = False) -> pd.DataFrame:
    """Return complete fitted-line columns accepted by the table helper."""
    columns = [
        "wavelength",
        "order",
        "slit_index",
        "slit_position",
        "detector_x",
        "detector_y",
        "observed_x",
        "observed_y",
        "x_diff",
        "y_diff",
        "fit_x",
        "fit_y",
        "residuals_x",
        "residuals_y",
        "residuals_xy",
        "sigma_clipped",
        "sharpness",
        "roundness1",
        "roundness2",
        "npix",
        "sky",
        "peak",
        "flux",
        "fwhm_pin_px",
        "R_pin",
        "pixelScaleNm",
        "detector_x_shifted",
        "detector_y_shifted",
        "R_slit",
        "fwhm_slit_px",
    ]
    values: dict[str, list[object]] = {column: [1.0] for column in columns}
    values.update({"wavelength": [500.0], "order": [10], "slit_index": [4]})
    if includeIon:
        values["ion"] = ["Ar"]
    return pd.DataFrame(values)


def test_store_init_params_deep_copies_recipe_settings(log: object) -> None:
    mapper = _mapper(log)
    recipeSettings = {"nested": {"degree": 2}}

    mapper._store_init_params(
        {"workspace-root-dir": "/tmp"},
        recipeSettings,
        "frame",
        "first.fits",
        "orders.fits",
        "qc",
        "products",
        "input.sof",
        True,
        "detections",
        "2024-01-01",
        "arc",
        True,
        copy,
        True,
    )
    recipeSettings["nested"]["degree"] = 9

    assert mapper.recipeSettings == {"nested": {"degree": 2}}
    assert mapper.pinholeFrame == "frame"
    assert mapper.firstGuessMap == "first.fits"
    assert mapper.turnOffMP is True


def test_constructor_preserves_inputs_and_initializes_output_workspace(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Construction resolves frame metadata before creating isolated outputs."""
    module = importlib.import_module("soxspipe.commonutils.create_dispersion_map")
    toolkitModule = importlib.import_module("soxspipe.commonutils.toolkit")
    commonutilsModule = importlib.import_module("soxspipe.commonutils")
    frame = SimpleNamespace(
        header={
            "SEQ_ARM": "VIS",
            "DATE_OBS": "2025-01-02T03:04:05",
            "INSTRUME": "SOXS",
            "EXPTIME": 30.0,
        }
    )
    monkeypatch.setattr(
        module,
        "keyword_lookup",
        lambda **_: SimpleNamespace(get=lambda name: name),
    )
    monkeypatch.setattr(toolkitModule, "get_calibration_lamp", lambda **_: "Ar")
    monkeypatch.setattr(
        commonutilsModule,
        "detector_lookup",
        lambda **_: SimpleNamespace(get=lambda arm: {"dispersion-axis": "x"}),
    )
    monkeypatch.setattr(
        toolkitModule,
        "utility_setup",
        lambda **_: ("/tmp/qc", "/tmp/products"),
    )

    mapper = create_dispersion_map(
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

    assert mapper.recipeName == "soxs-disp-solution"
    assert (mapper.arm, mapper.dateObs, mapper.inst, mapper.exptime) == (
        "VIS",
        "2025-01-02T03:04:05",
        "SOXS",
        30.0,
    )
    assert (mapper.axisA, mapper.axisB) == ("x", "y")
    assert (mapper.qcDir, mapper.productDir) == ("/tmp/qc", "/tmp/products")


@pytest.mark.parametrize(
    ("arm", "lamp", "raises"),
    [("VIS", "Ne", False), ("NIR", "ArHg", False), ("NIR", "Ar", True)],
)
def test_nir_lamp_validation_requires_argon_and_mercury(
    log: object,
    arm: str,
    lamp: str,
    raises: bool,
) -> None:
    mapper = _mapper(log)
    mapper.arm = arm
    mapper.lamp = lamp

    if raises:
        with pytest.raises(Exception, match="requires both Argon"):
            mapper._validate_nir_lamp_requirements()
    else:
        mapper._validate_nir_lamp_requirements()


@pytest.mark.parametrize(
    ("firstGuessMap", "expected"),
    [(False, "soxs-disp-solution"), ("first.fits", "soxs-spat-solution")],
)
def test_recipe_name_tracks_first_guess_mode(
    log: object,
    firstGuessMap: bool | str,
    expected: str,
) -> None:
    mapper = _mapper(log)
    mapper.firstGuessMap = firstGuessMap

    mapper._set_recipe_name()

    assert mapper.recipeName == expected


@pytest.mark.parametrize(
    ("dispersionAxis", "expected"),
    [("x", ("x", "y")), ("y", ("y", "x"))],
)
def test_image_orientation_follows_detector_configuration(
    log: object,
    dispersionAxis: str,
    expected: tuple[str, str],
) -> None:
    mapper = _mapper(log)
    mapper.detectorParams = {"dispersion-axis": dispersionAxis}

    mapper._set_image_orientation()

    assert (mapper.axisA, mapper.axisB) == expected


@pytest.mark.parametrize(
    ("firstGuessMap", "orderDegree", "expectedSlitDegree"),
    [(False, 2, 0), (False, [2, 3], [0, 0]), (True, 2, 1)],
)
def test_recipe_settings_initialize_slit_degree_by_mode(
    log: object,
    firstGuessMap: bool,
    orderDegree: int | list[int],
    expectedSlitDegree: int | list[int],
) -> None:
    mapper = _mapper(log)
    mapper.firstGuessMap = firstGuessMap
    mapper.recipeSettings["order-deg"] = orderDegree

    bootstrap, tightFit, returnedOrder, wavelength, slit = (
        mapper._initialize_recipe_settings()
    )

    assert (bootstrap, tightFit) == (False, False)
    assert returnedOrder == orderDegree
    assert wavelength == 3
    assert slit == expectedSlitDegree
    assert mapper.windowSize == 6


def test_prepare_pinhole_frame_masks_pixels_and_records_statistics(log: object) -> None:
    mapper = _mapper(log)
    mapper.pinholeFrame = SimpleNamespace(
        data=np.arange(9.0).reshape(3, 3),
        mask=np.array(
            [[False, False, False], [False, True, False], [False, False, False]]
        ),
    )

    masked = mapper._prepare_pinhole_frame()

    assert bool(masked.mask[1, 1])
    assert mapper.pinholeFrameMasked is masked
    assert np.isfinite(mapper.meanFrameFlux)
    assert np.isfinite(mapper.stdFrameFlux)


@pytest.mark.parametrize("dispersionAxis", ["x", "y"])
def test_resolution_measurement_fits_a_gaussian_slit_profile(
    log: object,
    dispersionAxis: str,
) -> None:
    """Measured synthetic arc profiles produce finite slit width and resolution."""
    mapper = _mapper(log)
    mapper.arm = "VIS"
    mapper.slitWidth = "1.0"
    mapper.detectorParams = {"dispersion-axis": dispersionAxis}
    coordinates = np.arange(64, dtype=float)
    profile = 100.0 * np.exp(-0.5 * ((coordinates - 32.0) / 2.2) ** 2)
    image = np.zeros((64, 64), dtype=float)
    if dispersionAxis == "x":
        image[:, 32] = profile
    else:
        image[32, :] = profile
    mapper.arcFrame = SimpleNamespace(data=image)
    row = pd.Series(
        {
            "observed_x": 32.0,
            "observed_y": 32.0,
            "pixelScaleNm": 0.1,
            "wavelength": 500.0,
        }
    )

    resolution, fwhm = mapper._calculate_resolution_on_slit(row)

    assert fwhm > 0
    assert resolution == pytest.approx(500.0 / (0.1 * fwhm))


def test_position_differences_add_group_means_once(log: object) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame(
        {
            "wavelength": [500.0, 500.0],
            "order": [10, 10],
            "detector_x": [10.0, 12.0],
            "detector_y": [20.0, 24.0],
            "detector_x_shifted": [11.0, 13.0],
            "detector_y_shifted": [19.0, 23.0],
            "observed_x": [10.0, 10.0],
            "observed_y": [20.0, 20.0],
        }
    )

    result = mapper._calculate_position_differences(table)
    repeated = mapper._calculate_position_differences(result)

    np.testing.assert_allclose(result["xy_diff"], [np.sqrt(2), np.sqrt(18)])
    np.testing.assert_array_equal(result["mph_mean_x"], [11.0, 11.0])
    np.testing.assert_array_equal(repeated["mph_mean_y"], [22.0, 22.0])


@pytest.mark.parametrize(
    ("iteration", "tightFit", "firstGuessMap", "expected"),
    [
        (0, True, False, (3, 4, False, False, True)),
        (0, False, False, (18, 10, True, False, True)),
        (1, False, False, (10, 10, True, False, True)),
        (2, False, False, (3, 4, False, False, True)),
        (0, False, True, (18, 5, False, False, True)),
        (1, False, True, (8, 4, False, False, True)),
    ],
)
def test_detection_parameters_characterize_iteration_modes(
    log: object,
    iteration: int,
    tightFit: bool,
    firstGuessMap: bool,
    expected: tuple[int, int, bool, bool, bool],
) -> None:
    mapper = _mapper(log)
    mapper.windowSize = 6
    mapper.firstGuessMap = firstGuessMap

    assert mapper._get_detection_parameters(iteration, tightFit) == expected


def test_multiple_detection_lists_explode_and_coerce_numeric_values(
    log: object,
) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame(
        {"order": [10], "observed_x": [["1.5", "bad"]], "observed_y": [[2, 3]]}
    )

    result = mapper._explode_multiple_detections(table)

    assert len(result) == 2
    assert result.loc[0, "observed_x"] == pytest.approx(1.5)
    assert np.isnan(result.loc[1, "observed_x"])


@pytest.mark.parametrize("recipeName", ["soxs-disp-solution", "soxs-spat-solution"])
def test_detection_qc_metrics_append_complete_rows(
    log: object,
    recipeName: str,
) -> None:
    mapper = _mapper(log)
    mapper.recipeName = recipeName

    mapper._write_qc_metrics(10, 8, 0.8, 0.75, "2024-01-02T00:00:00")

    assert mapper.qc["qc_name"].tolist() == [
        "DETLINES TOT",
        "DETLINES NUM",
        "DETLINES FRAC",
        "GOODLINES FRAC",
    ]
    assert mapper.qc["qc_value"].tolist() == [10, 8, 0.8, 0.75]
    assert mapper.qc["to_header"].tolist() == [True] * 4


def test_cluster_center_falls_back_after_relaxing_search(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapper = _mapper(log)
    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.create_dispersion_map"
    )
    monkeypatch.setattr(
        dispersionModule,
        "find_largest_cluster_center",
        lambda *args, **kwargs: (None, None),
    )

    assert mapper._find_cluster_center_with_fallback([1], [2], 3.5) == (3.5, None)


def test_cluster_shift_updates_detector_coordinates_when_found(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapper = _mapper(log)
    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.create_dispersion_map"
    )
    monkeypatch.setattr(
        dispersionModule,
        "find_largest_cluster_center",
        lambda *args, **kwargs: (1.0, -2.0),
    )
    table = pd.DataFrame(
        {
            "x_diff": [1.0, 1.1],
            "y_diff": [-2.0, -1.9],
            "detector_x_shifted": [10.0, 11.0],
            "detector_y_shifted": [20.0, 21.0],
        }
    )
    mask = pd.Series([True, True])

    result, centreX, centreY = mapper._find_and_apply_cluster_shift(table, mask)

    assert (centreX, centreY) == (1.0, -2.0)
    np.testing.assert_array_equal(result["detector_x_shifted"], [9.0, 10.0])
    np.testing.assert_array_equal(result["detector_y_shifted"], [22.0, 23.0])


def test_cluster_shift_leaves_empty_selection_unchanged(log: object) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame(
        {
            "x_diff": [1.0],
            "y_diff": [2.0],
            "detector_x_shifted": [3.0],
            "detector_y_shifted": [4.0],
        }
    )

    result, centreX, centreY = mapper._find_and_apply_cluster_shift(
        table,
        pd.Series([False]),
    )

    pd.testing.assert_frame_equal(result, table)
    assert (centreX, centreY) == (None, None)


def test_order_shift_statistics_return_robust_axis_values(log: object) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame(
        {"xy_diff": [1.0, 1.1], "x_diff": [2.0, 2.0], "y_diff": [-1.0, -1.0]}
    )

    medianX, medianY, stdX, stdY, medianCombined, stdCombined = (
        mapper._calculate_order_shift_statistics(table, pd.Series([True, True]))
    )

    assert medianX == pytest.approx(2.0)
    assert medianY == pytest.approx(-1.0)
    assert stdX == pytest.approx(0.0)
    assert stdY == pytest.approx(0.0)
    assert medianCombined == pytest.approx(1.05)
    assert stdCombined == pytest.approx(0.05)


def test_multipinhole_shift_returns_after_late_iteration_and_rejects_short_set(
    log: object,
) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame({"slit_index": [0], "x_diff": [1.0], "y_diff": [1.0]})
    mask = pd.Series([True])

    assert mapper._handle_multipin_hole_big_shift(table, 10, mask, 2) is table
    with pytest.raises(ValueError, match="COULD NOT DETECT ANY PINHOLES"):
        mapper._handle_multipin_hole_big_shift(table, 10, mask, 0)


def test_multipinhole_shift_corrects_the_axis_with_the_larger_end_shift(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A material slit-end displacement is removed from every selected line."""
    mapper = _mapper(log)
    mapper.axisA = "x"
    mapper.axisB = "y"
    mapper.debug = False
    table = pd.DataFrame(
        {
            "slit_index": [0] * 25 + [8] * 25,
            "x_diff": [-0.8] * 25 + [0.1] * 25,
            "y_diff": [0.0] * 50,
            "detector_x_shifted": np.arange(50.0),
        }
    )
    mapper._find_cluster_center_with_fallback = lambda x, y, fallback: (
        float(np.median(x)),
        0.0,
    )

    result = mapper._handle_multipin_hole_big_shift(
        table,
        10,
        pd.Series([True] * len(table)),
        1,
    )

    np.testing.assert_allclose(result["detector_x_shifted"], np.arange(50.0) + 0.8)


def test_straighten_mph_sets_projects_points_without_mutating_input(log: object) -> None:
    """Tilt correction projects observations onto their fitted pinhole line."""
    original = pd.DataFrame(
        {"observed_x": [0.0, 1.0, 2.0], "observed_y": [1.0, 4.0, 5.0]}
    )

    corrected = straighten_mph_sets(original)

    assert "tilt_corrected_x" not in original
    slope, intercept = np.polyfit(original["observed_x"], original["observed_y"], 1)
    np.testing.assert_allclose(
        corrected["tilt_corrected_y"],
        slope * corrected["tilt_corrected_x"] + intercept,
    )


@pytest.mark.parametrize(
    ("poptX", "poptY", "wavelength", "order", "slit", "expected"),
    [
        (None, None, 3, 2, 1, (2, 2, 1)),
        ("xerror", None, [3, 4], [2, 2], [1, 1], ([2, 4], [2, 2], [1, 1])),
        (None, "yerror", [3, 4], [2, 2], [1, 1], ([3, 3], [2, 2], [1, 1])),
    ],
)
def test_polynomial_degree_reduction_targets_failed_axis(
    log: object,
    poptX: str | None,
    poptY: str | None,
    wavelength: int | list[int],
    order: int | list[int],
    slit: int | list[int],
    expected: tuple[object, object, object],
) -> None:
    mapper = _mapper(log)
    mapper.firstGuessMap = True

    result = mapper._reduce_polynomial_degrees(
        poptX,
        poptY,
        wavelength,
        order,
        slit,
    )

    assert result == expected
    assert mapper.recipeSettings["wavelength-deg"] == expected[0]
    assert mapper.recipeSettings["order-deg"] == expected[1]
    assert mapper.recipeSettings["slit-deg"] == expected[2]


@pytest.mark.parametrize(
    ("technique", "expected"),
    [
        ("ECHELLE,PINHOLE", "single"),
        ("ECHELLE,MULTI-PINHOLE", "multi"),
    ],
)
def test_frame_technique_maps_supported_headers(
    log: object,
    technique: str,
    expected: str,
) -> None:
    mapper = _mapper(log)
    mapper.kw = lambda key: key
    mapper.pinholeFrame = SimpleNamespace(header={"DPR_TECH": technique})

    assert mapper._determine_frame_tech() == expected


def test_frame_technique_rejects_unsupported_header(log: object) -> None:
    mapper = _mapper(log)
    mapper.kw = lambda key: key
    mapper.pinholeFrame = SimpleNamespace(header={"DPR_TECH": "IMAGE"})

    with pytest.raises(TypeError, match="single- or multi-pinhole"):
        mapper._determine_frame_tech()


@pytest.mark.parametrize(
    ("arm", "header", "expected"),
    [
        ("VIS", {"WIN_BINX": "2", "WIN_BINY": "4"}, (2, 4)),
        ("VIS", {}, (1, 1)),
        ("NIR", {"WIN_BINX": "2", "WIN_BINY": "4"}, (1, 1)),
    ],
)
def test_binning_parameters_default_for_nir_or_missing_cards(
    log: object,
    arm: str,
    header: dict[str, str],
    expected: tuple[int, int],
) -> None:
    mapper = _mapper(log)
    mapper.arm = arm
    mapper.kw = lambda key: key
    mapper.pinholeFrame = SimpleNamespace(header=header)

    assert mapper._get_binning_params() == expected


def test_line_list_cleaning_standardizes_types_and_removes_flags(log: object) -> None:
    mapper = _mapper(log)
    table = pd.DataFrame(
        {
            "ORDER": ["10", "10", "11"],
            "WAVELENGTH": ["500", "500", "501"],
            "slit_index": ["4", "4", "4"],
            "ion": [1, 1, 2],
            "delete": [0, 0, 1],
        }
    )

    result = mapper._clean_line_list(table)

    assert len(result) == 1
    assert result.iloc[0]["order"] == pytest.approx(10.0)
    assert result.iloc[0]["wavelength"] == pytest.approx(500.0)
    assert result.iloc[0]["slit_index"] == 4
    assert result.iloc[0]["ion"] == "1"


@pytest.mark.parametrize(
    ("instrument", "arm", "expectedX", "expectedY"),
    [("XSHOOTER", "VIS", 9.0, 19.0), ("SOXS", "VIS", 10.0, 20.0)],
)
def test_coordinate_transform_characterizes_instrument_indexing(
    log: object,
    instrument: str,
    arm: str,
    expectedX: float,
    expectedY: float,
) -> None:
    mapper = _mapper(log)
    mapper.inst = instrument
    mapper.arm = arm
    mapper.detectorParams = {
        "science-pixels": {
            "columns": {"start": 0, "end": 100},
            "rows": {"start": 0, "end": 100},
        }
    }
    table = pd.DataFrame({"detector_x": [10.0], "detector_y": [20.0]})

    result = mapper._apply_coordinate_transforms(table)

    assert result.iloc[0]["detector_x"] == expectedX
    assert result.iloc[0]["detector_y"] == expectedY


def test_mid_slit_filter_and_incomplete_set_removal_preserve_complete_groups(
    log: object,
) -> None:
    mapper = _mapper(log)
    mapper.detectorParams = {"mid_slit_index": 4}
    table = pd.DataFrame(
        {
            "wavelength": [500.0, 500.0, 501.0],
            "order": [10, 10, 10],
            "slit_index": [0, 4, 4],
        }
    )

    filtered = mapper._filter_to_mid_slit(table)
    complete = mapper._remove_incomplete_mph_sets(table.copy(deep=True))

    assert filtered["slit_index"].tolist() == [4, 4]
    assert complete["wavelength"].tolist() == [500.0, 500.0]


def test_first_guess_corrections_apply_mid_slit_offset_to_complete_group(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapper = _mapper(log)
    mapper.firstGuessMap = "first-guess.fits"
    mapper.detectorParams = {"mid_slit_index": 4}
    table = pd.DataFrame(
        {
            "wavelength": [500.0, 500.0],
            "order": [10, 10],
            "slit_index": [4, 3],
            "detector_x": [12.0, 15.0],
            "detector_y": [17.0, 20.0],
        }
    )
    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.create_dispersion_map"
    )
    monkeypatch.setattr(
        dispersionModule,
        "dispersion_map_to_pixel_arrays",
        lambda **kwargs: kwargs["orderPixelTable"].assign(
            fit_x=[10.0, 11.0],
            fit_y=[20.0, 22.0],
        ),
    )

    result = mapper._apply_first_guess_corrections(table)

    assert result[["detector_x_shifted", "detector_y_shifted"]].values.tolist() == [
        [10.0, 20.0],
        [13.0, 23.0],
    ]
    assert "fit_x" not in result.columns
    assert "shift_x" not in result.columns


def test_predicted_line_list_runs_mode_specific_collaborators(log: object) -> None:
    mapper = _mapper(log)
    calls: list[str] = []
    initial = pd.DataFrame({"value": [1]})
    mapper._determine_frame_tech = lambda: calls.append("tech") or "single"
    mapper._get_binning_params = lambda: calls.append("binning") or (1, 1)
    mapper._load_predicted_lines = lambda technique, binX, binY: (
        calls.append("load") or initial
    )
    mapper._clean_line_list = lambda table: calls.append("clean") or table
    mapper._apply_coordinate_transforms = lambda table: (
        calls.append("coordinates") or table
    )
    mapper._filter_to_mid_slit = lambda table: calls.append("filter") or table
    mapper._apply_first_guess_corrections = lambda table: (
        calls.append("first_guess") or table
    )

    assert mapper.get_predicted_line_list() is initial
    assert calls == ["tech", "binning", "load", "clean", "coordinates", "filter"]

    calls.clear()
    mapper.firstGuessMap = "first.fits"
    assert mapper.get_predicted_line_list() is initial
    assert calls == [
        "tech",
        "binning",
        "load",
        "clean",
        "coordinates",
        "first_guess",
    ]


def test_output_filenames_and_clipped_line_qc_use_sof_name(log: object) -> None:
    mapper = _mapper(log)
    mapper.sofName = "SCIENCE"

    assert mapper._get_output_filenames() == (
        "SCIENCE_FITTED_LINES.fits",
        "SCIENCE_MISSED_LINES.fits",
    )
    mapper._write_line_list_qc(pd.DataFrame({"line": [1, 2]}), "2024-01-02")

    assert mapper.CLINE == 2
    assert mapper.qc.iloc[-1]["qc_name"] == "DETLINES CLIP NUM"
    assert mapper.qc.iloc[-1]["qc_value"] == 2


def test_output_filenames_derive_line_list_names_from_frame(
    log: object,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    mapper = _mapper(log)
    mapper.sofName = None
    mapper.pinholeFrame = object()
    dispersionModule = importlib.import_module(
        "soxspipe.commonutils.create_dispersion_map"
    )
    monkeypatch.setattr(dispersionModule, "filenamer", lambda **kwargs: "ARC.fits")

    assert mapper._get_output_filenames() == (
        "ARC_FITTED_LINES.fits",
        "ARC_MISSED_LINES.fits",
    )


def test_multipinhole_shift_corrects_larger_first_iteration_offset(
    log: object,
) -> None:
    mapper = _mapper(log)
    mapper.axisA = "x"
    mapper.axisB = "y"
    slitIndexes = [8] * 25 + [0] * 25
    table = pd.DataFrame(
        {
            "slit_index": slitIndexes,
            "x_diff": [1.0] * 25 + [0.2] * 25,
            "y_diff": [0.0] * 50,
            "detector_x_shifted": [20.0] * 50,
        }
    )
    mapper._find_cluster_center_with_fallback = (
        lambda xData, yData, fallbackMedian: (float(xData.iloc[0]), 0.0)
    )

    result = mapper._handle_multipin_hole_big_shift(
        table,
        10,
        pd.Series([True] * len(table)),
        0,
    )

    np.testing.assert_array_equal(result["detector_x_shifted"], [19.0] * 50)


@pytest.mark.parametrize("includeClipped", [False, True])
def test_line_list_columns_sort_and_mark_clipped_rows(
    log: object,
    includeClipped: bool,
) -> None:
    mapper = _mapper(log)
    good = _line_table(includeIon=True)
    clipped = _line_table(includeIon=True) if includeClipped else good.iloc[0:0].copy()

    combined, sortedGood = mapper._prepare_line_list_columns(good, clipped)

    assert combined.iloc[0]["ion"] == "Ar"
    assert sortedGood["order"].tolist() == [10]
    assert len(combined) == (2 if includeClipped else 1)
    if includeClipped:
        assert bool(combined.iloc[0]["sigma_clipped"])


def test_line_list_writers_create_fits_files_and_register_qc_products(
    log: object,
    tmp_path: Path,
) -> None:
    mapper = _mapper(log)
    mapper.qcDir = str(tmp_path)
    mapper.settings = {"tune-pipeline": False}
    fittedLines = _line_table(includeIon=True)
    missingLines = fittedLines[
        [
            "wavelength",
            "order",
            "slit_index",
            "slit_position",
            "detector_x",
            "detector_y",
        ]
    ].copy()

    mapper._write_fitted_lines_file(fittedLines, "fitted.fits", "2024-01-02")
    mapper._write_missing_lines_file(missingLines, "missing.fits", "2024-01-02")

    assert "ion" in Table.read(tmp_path / "fitted.fits").colnames
    assert Table.read(tmp_path / "missing.fits").colnames == list(missingLines.columns)
    assert mapper.products["product_label"].tolist() == [
        "DISP_MAP_LINES",
        "DISP_MAP_LINES_MISSING",
    ]
    assert mapper.products["file_name"].tolist() == ["fitted.fits", "missing.fits"]


def test_calculate_residuals_returns_analytic_pixel_scale_and_resolution(
    log: object,
) -> None:
    """Calculate fitted positions and resolution for a linear two-line solution."""
    mapper = _mapper(log)
    mapper.arcFrame = False
    table = pd.DataFrame(
        {
            "order": [1.0, 1.0],
            "wavelength": [10.0, 20.0],
            "slit_position": [0.0, 0.0],
            "order_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_0": [1.0, 1.0],
            "wavelength_pow_x_1": [10.0, 20.0],
            "slit_position_pow_x_0": [1.0, 1.0],
            "order_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_0": [1.0, 1.0],
            "wavelength_pow_y_1": [10.0, 20.0],
            "slit_position_pow_y_0": [1.0, 1.0],
            "observed_x": [17.0, 40.0],
            "observed_y": [1.0, 5.0],
            "detector_x": [20.0, 40.0],
            "detector_y": [5.0, 5.0],
            "fwhm_pin_px": [2.0, 2.0],
        }
    )

    mean, standardDeviation, median, result = mapper.calculate_residuals(
        table,
        xcoeff=[0.0, 2.0],
        ycoeff=[5.0, 0.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
        pixelRange=True,
        writeQCs=True,
    )

    np.testing.assert_allclose(table["fit_x"], [20.0, 40.0])
    np.testing.assert_allclose(table["fit_y"], [5.0, 5.0])
    assert result is table
    np.testing.assert_allclose(table["residuals_xy"], [5.0, 0.0])
    assert (mean, standardDeviation, median) == pytest.approx((2.5, 2.5, 2.5))
    np.testing.assert_allclose(table["pixelScaleNm"], [0.5, 0.5])
    np.testing.assert_allclose(table["delta_wavelength"], [1.0, 1.0])
    np.testing.assert_allclose(table["R_pin"], [10.0, 20.0])
    assert set(["fit_x_high", "fit_y_high", "fit_x_low", "fit_y_low"]).isdisjoint(
        table.columns
    )
    np.testing.assert_allclose(table["R_slit"], [0.0, 0.0])
    np.testing.assert_allclose(table["fwhm_slit_px"], [0.0, 0.0])
    assert len(mapper.qc.index) == 38
    assert set(mapper.qc["qc_name"]) >= {
        "XY RES MEDIAN",
        "FWHM PIN MEDIAN",
        "R PIN MEDIAN",
    }
    assert set(mapper.qc["qc_order"].dropna()) == {"r"}


def test_largest_cluster_center_ignores_noise_and_returns_none_without_cluster() -> (
    None
):
    xCoordinates = np.array([0.0, 0.1, -0.1, 10.0])
    yCoordinates = np.array([1.0, 1.1, 0.9, 10.0])

    centre = find_largest_cluster_center(
        xCoordinates,
        yCoordinates,
        eps=0.5,
        min_samples=2,
    )
    missing = find_largest_cluster_center(
        xCoordinates,
        yCoordinates,
        eps=0.01,
        min_samples=3,
    )

    assert centre[0] == pytest.approx(0.0)
    assert centre[1] == pytest.approx(1.0)
    assert missing == (None, None)
