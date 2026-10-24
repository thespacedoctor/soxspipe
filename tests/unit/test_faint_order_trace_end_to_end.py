"""The faint-order trace through the real sampler and global fit, on the synthetic three-source VIS scene (DY-1941)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import yaml

import soxspipe
from soxspipe.commonutils.detect_continuum import detect_continuum
from soxspipe.commonutils.faint_order_trace import make_pixels_at
from tests.unit._continuum_vis_scene import (
    VIS_ORDER_SAMPLES,
    VIS_START_NM,
    TraceFrame,
    vis_dispersion_map,
    vis_row,
    vis_transfer_image,
    vis_true_centre,
)

pytestmark = pytest.mark.unit

PACKAGE_DIR = Path(soxspipe.__file__).parent
FAINT_QC = "N ORDERS FAINT TRACE"
BRIGHT_AMPLITUDES = {1.0: 100.0, 2.0: 100.0, 3.0: 100.0, 4.0: 100.0}
FLANKS = {1.0: 30.0, 2.0: 0.8, 3.0: 30.0, 4.0: 30.0}
FAINT_U_AMPLITUDES = BRIGHT_AMPLITUDES | {2.0: 0.5}


def _stare_block(isOn: bool) -> dict:
    settings = yaml.safe_load((PACKAGE_DIR / "soxs_default_settings.yaml").read_text())
    return dict(settings["soxs-stare"]["vis"]["detect-continuum"]) | {"faint-order-trace": isOn}


def _sample_grid(pixelsAt: object) -> pd.DataFrame:
    """Return the sample grid of the four orders at slit zero, as `create_pixel_arrays` would."""
    frames = []
    for order in sorted(VIS_START_NM):
        wavelength = VIS_START_NM[order] + np.linspace(0.0, 200.0, VIS_ORDER_SAMPLES, endpoint=False)
        axisA, axisB = pixelsAt(order, wavelength, 0.0)
        frames.append(
            pd.DataFrame(
                {"order": order, "wavelength": wavelength, "slit_position": 0.0, "fit_x": axisA, "fit_y": axisB}
            )
        )
    return pd.concat(frames, ignore_index=True)


def _run(log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path, amplitudes: dict, isOn: bool) -> tuple:
    """Run `get()` on the scene with the real sampler and fit, and return (detector, result)."""
    block = _stare_block(isOn)
    tmp_path.mkdir(parents=True, exist_ok=True)
    detector = object.__new__(detect_continuum)
    detector.log = log
    detector.inst, detector.arm, detector.recipeName = "SOXS", "VIS", "soxs-stare"
    detector.axisA, detector.axisB = "x", "y"
    detector.sliceAxis, detector.sliceAntiAxis = "x", "y"
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.recipeSettings = block
    detector.orderDeg, detector.axisBDeg = block["order-deg"], block["disp-axis-deg"]
    detector.coeff_dict = {"degorder_cent": detector.orderDeg, "degy_cent": detector.axisBDeg}
    detector.settings = {"tune-pipeline": False}
    detector.dateObs = "2024-01-02T00:00:00"
    detector.exptime = 1.0
    detector.qc = pd.DataFrame()
    detector.products = False
    detector.noddingSequence = ""
    detector.notFlattened = ""
    detector.debug = False
    detector.orderPixelTable = False
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.dispersion_map = vis_dispersion_map(tmp_path)
    image = vis_transfer_image(amplitudes, FLANKS)
    detector.traceFrame = np.ma.masked_array(image, mask=np.zeros(image.shape, dtype=bool)).view(TraceFrame)
    pixelsAt = make_pixels_at(log, detector.dispersion_map, {"x": 1.0, "y": 1.0}, "x", "y")
    monkeypatch.setattr(detector, "create_pixel_arrays", lambda: (_sample_grid(pixelsAt), 1, 1))
    monkeypatch.setattr("soxspipe.commonutils.toolkit.quicklook_image", lambda **_: None)
    monkeypatch.setattr(detector, "plot_results", lambda **_: "/qc/continuum.pdf")
    monkeypatch.setattr(detector, "write_order_table_to_file", lambda **_: "/products/orders.fits")
    monkeypatch.setattr(detector, "build_order_metadata", lambda **_: pd.DataFrame({"axis": ["x"]}))
    return detector, detector.get()


def _trace_error(detector: detect_continuum, polyTable: pd.DataFrame, order: float) -> np.ndarray:
    """Return the fitted trace minus the true target position along an order, in pixels."""
    wavelength = VIS_START_NM[order] + np.linspace(0.0, 200.0, VIS_ORDER_SAMPLES, endpoint=False)
    fitted = detector._fitted_trace(polyTable)(order, vis_row(order, wavelength))
    return fitted - vis_true_centre(order, wavelength)


def test_the_fitted_u_trace_lands_within_0_3_px_of_the_target_at_a_snr_of_0_5(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ARRANGE: THE u PEAK IS HALF THE NOISE PER PIXEL, AND TWO SOURCES 3.7 ARCSEC AWAY ARE BRIGHTER THAN IT

    # ACT
    detector, result = _run(log, monkeypatch, tmp_path, FAINT_U_AMPLITUDES, isOn=True)

    # ASSERT
    error = _trace_error(detector, result[3], 2.0)
    assert np.abs(error).max() < 0.3
    assert FAINT_QC in result[1]["qc_name"].tolist()


def test_the_same_scene_without_the_transfer_misplaces_the_u_trace_by_pixels(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ACT
    detector, result = _run(log, monkeypatch, tmp_path, FAINT_U_AMPLITUDES, isOn=False)

    # ASSERT: THE DEFECT THE TRANSFER FIXES
    assert np.abs(_trace_error(detector, result[3], 2.0)).max() > 1.0
    assert FAINT_QC not in result[1]["qc_name"].tolist()


def test_the_bright_orders_barely_move_when_u_is_transferred(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ACT
    onDetector, onResult = _run(log, monkeypatch, tmp_path / "on", FAINT_U_AMPLITUDES, isOn=True)
    offDetector, offResult = _run(log, monkeypatch, tmp_path / "off", FAINT_U_AMPLITUDES, isOn=False)

    # ASSERT: g, r AND i DIFFER BY LESS THAN 0.05 PX BETWEEN THE TWO RUNS
    for order in (1.0, 3.0, 4.0):
        moved = _trace_error(onDetector, onResult[3], order) - _trace_error(offDetector, offResult[3], order)
        assert np.abs(moved).max() < 0.05


def test_a_bright_only_frame_gives_the_same_fit_and_qcs_with_the_trace_on_or_off(
    log: object, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    # ACT
    _, onResult = _run(log, monkeypatch, tmp_path / "on", BRIGHT_AMPLITUDES, isOn=True)
    _, offResult = _run(log, monkeypatch, tmp_path / "off", BRIGHT_AMPLITUDES, isOn=False)

    # ASSERT: BIT-IDENTICAL COEFFICIENTS, QCS (BUT FOR THE WALL-CLOCK STAMP) AND FITTED TABLE, WITH NO NEW COLUMN
    pd.testing.assert_frame_equal(onResult[3], offResult[3], check_exact=True)
    pd.testing.assert_frame_equal(
        onResult[1].drop(columns="reduction_date_utc"),
        offResult[1].drop(columns="reduction_date_utc"),
        check_exact=True,
    )
    pd.testing.assert_frame_equal(onResult[4], offResult[4], check_exact=True)
    assert "faint_transfer" not in onResult[4].columns
    assert FAINT_QC not in onResult[1]["qc_name"].tolist()
