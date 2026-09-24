"""Fast behavioral contracts for reusable toolkit functions."""

from __future__ import annotations

import logging
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty
from numpy.testing import assert_allclose, assert_array_equal

from soxspipe.commonutils import toolkit

pytestmark = pytest.mark.unit


def test_cut_image_slice_returns_centered_medians_and_coordinates(log: object) -> None:
    frame = np.arange(49, dtype=float).reshape(7, 7)

    horizontal, xOffset, yCentre = toolkit.cut_image_slice(
        log, frame, width=3, length=4, x=3, y=3, median=True
    )
    vertical, yOffset, xCentre = toolkit.cut_image_slice(
        log, frame, width=3, length=4, x=3, y=3, sliceAxis="y", median=True
    )

    assert_array_equal(horizontal, [22.0, 23.0, 24.0, 25.0])
    assert_array_equal(vertical, [10.0, 17.0, 24.0, 31.0])
    assert (xOffset, yCentre) == (1, 3.5)
    assert (yOffset, xCentre) == (1, 3.5)


def test_cut_image_slice_returns_raw_data_when_median_false(log: object) -> None:
    frame = np.arange(49, dtype=float).reshape(7, 7)

    result, offset, centre = toolkit.cut_image_slice(
        log, frame, width=3, length=4, x=3, y=3, median=False
    )

    assert_array_equal(
        result,
        [
            [15.0, 16.0, 17.0, 18.0],
            [22.0, 23.0, 24.0, 25.0],
            [29.0, 30.0, 31.0, 32.0],
        ],
    )
    assert (offset, centre) == (1, 3.5)


def test_cut_image_slice_rejects_invalid_axis_and_out_of_bounds(log: object) -> None:
    frame = np.zeros((7, 7))

    assert toolkit.cut_image_slice(log, frame, 3, 4, 1, 3) == (None, None, None)
    with pytest.raises(ValueError, match="either 'x' or 'y'"):
        toolkit.cut_image_slice(log, frame, 3, 4, 3, 3, sliceAxis="z")


@pytest.mark.parametrize(
    ("recipeName", "countName", "fractionName", "comment"),
    [
        ("soxs-mdark", "HOTPIX NUM", "HOTPIX FRAC", "hot pixels"),
        ("soxs-mflat", "COLDPIX NUM", "COLDPIX FRAC", "cold pixels"),
        ("soxs-stare", "BADPIX NUM", "BADPIX FRAC", "bad pixels"),
    ],
)
def test_generic_quality_checks_records_mask_count_and_fraction(
    monkeypatch: pytest.MonkeyPatch,
    log: object,
    recipeName: str,
    countName: str,
    fractionName: str,
    comment: str,
) -> None:
    monkeypatch.setattr(
        toolkit, "keyword_lookup", lambda **kwargs: SimpleNamespace(get=lambda key: key)
    )
    frame = SimpleNamespace(
        mask=np.array([[False, True], [True, False]]),
        header={"SEQ_ARM": "VIS", "DATE_OBS": "2024-01-02T03:04:05"},
    )

    result = toolkit.generic_quality_checks(log, frame, {}, recipeName, pd.DataFrame())

    assert result["qc_name"].tolist() == [countName, fractionName]
    assert result["qc_value"].tolist() == [2, 0.5]
    assert result["qc_comment"].tolist() == [
        f"Number of {comment}",
        f"Fraction of {comment}",
    ]
    assert result["to_header"].tolist() == [True, True]


def test_calibration_path_honors_instrument_setting(log: object) -> None:
    defaultPath = toolkit.get_calibrations_path(log, {})
    customPath = toolkit.get_calibrations_path(log, {"instrument": "xshooter"})

    assert defaultPath.endswith("resources/static_calibrations/soxs")
    assert customPath.endswith("resources/static_calibrations/xshooter")


@pytest.mark.parametrize(
    ("arguments", "expectedRecipe"),
    [
        (["soxspipe", "mbias"], "soxs-mbias"),
        (["soxspipe", "-v", "mflat"], "soxs-mflat"),
    ],
)
def test_predict_product_path_uses_observation_night_and_cli_recipe(
    monkeypatch: pytest.MonkeyPatch,
    arguments: list[str],
    expectedRecipe: str,
) -> None:
    class Organiser:
        def __init__(self, **kwargs: object) -> None:
            pass

        def session_list(self, silent: bool) -> tuple[str, list[str]]:
            return "session-01", ["session-01"]

        def close(self) -> None:
            pass

    monkeypatch.setattr("soxspipe.commonutils.data_organiser", Organiser)
    monkeypatch.setattr(sys, "argv", arguments)

    productPath, night = toolkit.predict_product_path("20240102T030405_STARE_STD.sof")

    assert night == "2024-01-01"
    assert productPath == (
        f"./sessions/session-01/reduced/2024-01-01/{expectedRecipe}/"
        "20240102T030405_STARE_STD_RESP.fits"
    )


def test_add_recipe_logger_replaces_handlers_and_separates_messages(
    tmp_path: Path,
) -> None:
    logger = logging.getLogger(f"toolkit-test-{id(tmp_path)}")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    logger.handlers.clear()
    productPath = tmp_path / "nested" / "product.fits"

    toolkit.add_recipe_logger(logger, str(productPath))
    toolkit.add_recipe_logger(logger, str(productPath))
    logger.log(logging.INFO + 1, "recipe progress")
    logger.warning("excluded warning")
    logger.error("recipe failure")
    for handler in logger.handlers:
        handler.flush()

    assert [handler.get_name() for handler in logger.handlers] == [
        "recipeLog",
        "recipeErr",
    ]
    assert (tmp_path / "nested" / "product.log").read_text() == "recipe progress\n"
    errorText = (tmp_path / "nested" / "product_ERROR.log").read_text()
    assert "recipe failure" in errorText
    assert "excluded warning" not in errorText
    for handler in tuple(logger.handlers):
        handler.close()
        logger.removeHandler(handler)


def test_max_filter_accepts_only_records_below_limit() -> None:
    filterObject = toolkit.MaxFilter(logging.WARNING)

    assert (
        filterObject.filter(logging.LogRecord("x", logging.INFO, "", 1, "", (), None))
        is True
    )
    assert (
        filterObject.filter(
            logging.LogRecord("x", logging.WARNING, "", 1, "", (), None)
        )
        is None
    )


def test_calibration_lamp_normalizes_and_concatenates_header_values(
    log: object,
) -> None:
    keywords = {f"LAMP{index}": f"ESO LAMP {index}" for index in range(1, 8)}
    frame = SimpleNamespace(
        header={
            "INSTRUME": "SOXS",
            "ESO LAMP 1": "UVB_High_lamp",
            "ESO LAMP 2": "VIS_Argo_Lamp",
            "ESO LAMP 3": "NIR_Neon_lamp",
        }
    )

    lamp = toolkit.get_calibration_lamp(log, frame, keywords.__getitem__)

    assert lamp == "QTHArNe"
    assert (
        toolkit.get_calibration_lamp(
            log, SimpleNamespace(header={"INSTRUME": "SOXS"}), keywords.__getitem__
        )
        is None
    )


def test_utility_setup_creates_recipe_directories_under_workspace(
    tmp_path: Path, log: object
) -> None:
    workspace = tmp_path / "workspace"

    qcDir, productDir = toolkit.utility_setup(
        log,
        {"workspace-root-dir": str(workspace)},
        "soxs-stare-obj",
        "2024-01-01",
    )

    assert Path(qcDir) == workspace / "qc" / "2024-01-01" / "soxs-stare"
    assert Path(productDir) == workspace / "reduced" / "2024-01-01" / "soxs-stare"
    assert Path(qcDir).is_dir()
    assert Path(productDir).is_dir()


@pytest.mark.parametrize("failingSegment", ["/qc/", "/reduced/"])
def test_utility_setup_propagates_directory_creation_failures(
    tmp_path: Path, log: object, monkeypatch: pytest.MonkeyPatch, failingSegment: str
) -> None:
    """Both directory-creation sites propagate; neither swallows."""
    workspace = tmp_path / "workspace"
    realMakedirs = toolkit.os.makedirs

    def refuse(path, *args, **kwargs):
        # REFUSE ONLY THE DIRECTORY UNDER TEST, SO EACH SITE IS PINNED SEPARATELY
        if failingSegment in str(path):
            raise PermissionError(13, "Permission denied", str(path))
        return realMakedirs(path, *args, **kwargs)

    monkeypatch.setattr(toolkit.os, "makedirs", refuse)

    with pytest.raises(PermissionError):
        toolkit.utility_setup(
            log,
            {"workspace-root-dir": str(workspace)},
            "soxs-stare-obj",
            "2024-01-01",
        )


def test_calculate_rolling_snr_places_finite_values_at_window_centers() -> None:
    source = pd.DataFrame({"flux": [10.0, 11.0, 9.0, 10.0, 10.0, 12.0, 8.0]})

    result = toolkit.calculate_rolling_snr(source.copy(), "flux", 5)
    tooShort = toolkit.calculate_rolling_snr(source.head(4).copy(), "flux", 5)

    firstWindow = 10.0 / (0.6052697 * 2.0)
    remainingWindows = 10.0 / (0.6052697 * 3.0)
    assert_allclose(
        result["SNR"],
        [firstWindow] * 3 + [remainingWindows] * 4,
        rtol=1e-7,
    )
    assert tooShort["SNR"].isna().all()


def test_frame_to_32_converts_data_and_uncertainty_in_place() -> None:
    frame = CCDData(
        np.ones((2, 2), dtype=np.float64),
        unit=u.electron,
        uncertainty=StdDevUncertainty(np.full((2, 2), 2.0, dtype=np.float64)),
    )

    result = toolkit.frame_to_32(frame)

    assert result is frame
    assert result.data.dtype == np.float32
    assert result.uncertainty.array.dtype == np.float32
    sentinel = object()
    assert toolkit.frame_to_32(sentinel) is sentinel


def test_add_snr_efficiency_qcs_records_global_and_per_order_values(
    log: object,
) -> None:
    spectrum = pd.DataFrame(
        {
            "WAVE": [400.0, 450.0, 550.0, 600.0],
            "SNR": [10.0, 20.0, 30.0, 40.0],
            "EFFICIENCY": [0.1, 0.2, 0.3, 0.4],
        }
    )

    result = toolkit.add_snr_efficiency_qcs(
        log,
        spectrum,
        pd.DataFrame(),
        {"1011": 500.0},
        "soxs-stare",
        "2024-01-02T03:04:05",
    )

    assert result["qc_name"].tolist() == [
        "EFF MEDIAN",
        "EFF MEDIAN",
        "EFF MEDIAN",
        "SNR MEDIAN",
        "SNR MEDIAN",
        "SNR MEDIAN",
    ]
    assert result["qc_value"].tolist() == [0.25, 0.35, 0.15, 25.0, 35.0, 15.0]
    assert result["qc_order"].iloc[1:].dropna().tolist() == [10.0, 11.0, 10.0, 11.0]
    assert "ORDER" not in spectrum.columns


def test_qc_settings_plot_tables_limits_qcs_and_hides_axes(log: object) -> None:
    import matplotlib.pyplot as plt

    figure, (qcAxis, settingsAxis) = plt.subplots(1, 2)
    qc = pd.DataFrame(
        {
            "qc_name": [f"QC {index}" for index in range(12)],
            "qc_value": range(12),
            "qc_unit": ["pixel"] * 12,
            "qc_comment": ["Synthetic"] * 12,
            "qc_order": [-1] * 11 + [2],
        }
    )

    result = toolkit.qc_settings_plot_tables(
        log,
        qc,
        qcAxis,
        {"threshold": 3, "vis": {"ignored": True}},
        settingsAxis,
    )

    assert result is None
    assert not qcAxis.axison
    assert not settingsAxis.axison
    assert len(qcAxis.tables) == 1
    assert len(settingsAxis.tables) == 1
    plt.close(figure)
