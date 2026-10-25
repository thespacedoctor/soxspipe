"""Shared synthetic VIS scene for the continuum-trace sampling tests (DY-1716, DY-1941)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from soxspipe.commonutils.detect_continuum import detect_continuum
from tests.factories import dispersion_map_fits


class TraceFrame(np.ma.MaskedArray):
    """A masked image that carries the empty `header` the sampler reads, one per instance."""

    header: dict

    def __array_finalize__(self, obj: object) -> None:
        super().__array_finalize__(obj)
        self.header = {}


VIS_ORDER_CENTRES = {1.0: 60.0, 2.0: 150.0, 3.0: 240.0, 4.0: 330.0}
VIS_ORDER_SAMPLE_COUNTS = {1.0: 30, 2.0: 50, 3.0: 100, 4.0: 30}
VIS_TRACE_SIGMA = 1.2
# THE TRANSFER SCENE: FOUR 2000-ROW ORDERS SIDE BY SIDE, U, G, R, I FROM 400 TO 1200 NM IN 200 NM STEPS, WITH THE
# SLIT RUNNING AGAINST THE X AXIS AS IT DOES ON THE FLIPPED VIS DETECTOR. THE DISPERSION MAP IS
# X = -30 + 90 ORDER + (VIS_SLIT_SCALE_BASE + VIS_SLIT_SCALE_PER_ORDER ORDER) SLIT AND
# Y = 10 + VIS_ROWS_PER_NM (WAVELENGTH - VIS_START_NM[ORDER])
VIS_SLIT_SCALE_BASE = -3.6
VIS_SLIT_SCALE_PER_ORDER = -0.1
VIS_START_NM = {1.0: 800.0, 2.0: 400.0, 3.0: 600.0, 4.0: 1000.0}
VIS_ROWS_PER_NM = 10.0
VIS_ORDER_ROWS = 2000
VIS_FRAME_SHAPE = (VIS_ORDER_ROWS + 20, 400)


def vis_slit_scale(order: float) -> float:
    """Return the signed pixels per slit arcsec of an order in the scene's dispersion map."""
    return VIS_SLIT_SCALE_BASE + VIS_SLIT_SCALE_PER_ORDER * order


def vis_row(order: float, wavelength: float | np.ndarray) -> float | np.ndarray:
    """Return the row of a wavelength in an order of the scene's dispersion map."""
    return 10.0 + VIS_ROWS_PER_NM * (wavelength - VIS_START_NM[order])


def vis_wavelength(order: float, row: float | np.ndarray) -> float | np.ndarray:
    """Return the wavelength at a row of an order of the scene's dispersion map."""
    return VIS_START_NM[order] + (row - 10.0) / VIS_ROWS_PER_NM


def vis_dispersion_map(directory: Path) -> str:
    """Write the scene's dispersion map under `directory` and return its path."""
    orders = sorted(VIS_START_NM)
    startCoefficients = np.polyfit(orders, [VIS_START_NM[order] for order in orders], 3)[::-1]
    names = [f"c{i}{j}{k}" for i in range(4) for j in (0, 1) for k in (0, 1)]
    xCoefficients = {"c000": -30.0, "c100": 90.0, "c001": VIS_SLIT_SCALE_BASE, "c101": VIS_SLIT_SCALE_PER_ORDER}
    yCoefficients = {"c010": VIS_ROWS_PER_NM, "c000": 10.0 - VIS_ROWS_PER_NM * startCoefficients[0]} | {
        f"c{power}00": -VIS_ROWS_PER_NM * startCoefficients[power] for power in (1, 2, 3)
    }
    rows = [
        {"axis": axis, "order_deg": 3, "wavelength_deg": 1, "slit_deg": 1}
        | {name: coefficients.get(name, 0.0) for name in names}
        for axis, coefficients in (("x", xCoefficients), ("y", yCoefficients))
    ]
    return str(dispersion_map_fits(directory / "vis_scene_map.fits", coefficients=pd.DataFrame(rows)))


def vis_samples() -> pd.DataFrame:
    """Return the sample grid of four VIS orders, each one a vertical stripe at its own map position."""
    frames = []
    for order, count in VIS_ORDER_SAMPLE_COUNTS.items():
        frames.append(
            pd.DataFrame(
                {
                    "order": order,
                    "wavelength": np.linspace(400.0, 400.0 + count - 1, count),
                    "fit_x": VIS_ORDER_CENTRES[order],
                    "fit_y": np.linspace(10.0, 10.0 + count - 1, count),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)


def vis_trace_frame(offsets: dict[float, float]) -> TraceFrame:
    """Return a noisy frame whose order traces sit `offsets[order]` pixels from the map position."""
    rng = np.random.default_rng(1716)
    columns = np.arange(400, dtype=float)
    image = rng.normal(0.0, 1.0, (120, 400))
    for order, centre in VIS_ORDER_CENTRES.items():
        profile = 100.0 * np.exp(-0.5 * ((columns - centre - offsets[order]) / VIS_TRACE_SIGMA) ** 2)
        image += profile[np.newaxis, :]
    return np.ma.masked_array(image, mask=np.zeros(image.shape, dtype=bool)).view(TraceFrame)


def vis_sampling_detector(detector: detect_continuum, offsets: dict[float, float]) -> detect_continuum:
    """Return `detector` set up as a SOXS VIS detector that samples the synthetic frame with the real slice fitter."""
    detector.inst = "SOXS"
    detector.arm = "VIS"
    detector.recipeName = "soxs-stare"
    detector.traceFrame = vis_trace_frame(offsets)
    detector.kw = lambda key: {"DPR_TYPE": "DPR TYPE", "WIN_BINX": "BINX", "WIN_BINY": "BINY"}[key]
    detector.recipeSettings.update({"slice-length": 70, "peak-sigma-limit": 3.0, "slice-width": 3})
    detector.debug = False
    detector.qc = pd.DataFrame()
    detector.dateObs = "2024-01-02T03:04:05"
    detector.sliceAxis = "x"
    detector.sliceAntiAxis = "y"
    detector.detectorParams = {"dispersion-axis": "x"}
    detector.create_pixel_arrays = lambda: (vis_samples(), 1, 1)
    return detector


# THE TARGET'S SLIT POSITION (ARCSEC) AT 400 NM AND ITS DRIFT TO 1200 NM, TWO FLANKING SOURCES THIS FAR EITHER SIDE
# OF IT, AND THE WIDTH OF EVERY SOURCE IN PIXELS
VIS_TARGET_SLIT_ARCSEC = 0.4
VIS_TARGET_DRIFT_ARCSEC = 0.5
VIS_FLANK_SEPARATION_ARCSEC = 3.7
VIS_SEEING_SIGMA_PX = 1.5
VIS_ORDER_SAMPLES = 150


def vis_target_slit(wavelength: float | np.ndarray) -> float | np.ndarray:
    """Return the slit position (arcsec) of the target at a wavelength."""
    return VIS_TARGET_SLIT_ARCSEC + VIS_TARGET_DRIFT_ARCSEC * (wavelength - 400.0) / 800.0


def vis_true_centre(order: float, wavelength: float | np.ndarray, separationArcsec: float = 0.0) -> np.ndarray:
    """Return the x pixel of a source `separationArcsec` from the target, in an order at a wavelength."""
    return -30.0 + 90.0 * order + vis_slit_scale(order) * (vis_target_slit(wavelength) + separationArcsec)


def vis_transfer_image(
    amplitudes: dict[float, float], flankAmplitudes: dict[float, float], noise: float = 1.0
) -> np.ndarray:
    """Return a frame with the target and two flanking sources in each order, over Gaussian noise of `noise` counts."""
    image = np.random.default_rng(1941).normal(0.0, noise, VIS_FRAME_SHAPE)
    rows = np.arange(10, 10 + VIS_ORDER_ROWS)
    columns = np.arange(VIS_FRAME_SHAPE[1], dtype=float)
    for order, amplitude in amplitudes.items():
        wavelength = vis_wavelength(order, rows)
        for separation, peak in (
            (0.0, amplitude),
            (VIS_FLANK_SEPARATION_ARCSEC, flankAmplitudes[order]),
            (-VIS_FLANK_SEPARATION_ARCSEC, flankAmplitudes[order]),
        ):
            centre = vis_true_centre(order, wavelength, separation)[:, np.newaxis]
            image[rows] += peak * np.exp(-0.5 * ((columns[np.newaxis, :] - centre) / VIS_SEEING_SIGMA_PX) ** 2)
    return image


def vis_transfer_table(faintOrders: tuple[float, ...], fitOffsetArcsec: float = 0.0) -> pd.DataFrame:
    """Return the sampled table of the four orders: clean traces, and noise fits for the orders in `faintOrders`.

    `fitOffsetArcsec` puts the `fit_x` of every order that far from the map's slit zero, as the second sampling pass
    does when it moves `fit_x` onto an order's median trace. `centre_shift` stays `fit_x - cont_x`.
    """
    rng = np.random.default_rng(7)
    frames = []
    for order in sorted(VIS_START_NM):
        wavelength = VIS_START_NM[order] + np.linspace(0.0, 200.0, VIS_ORDER_SAMPLES, endpoint=False)
        fitX = np.full(VIS_ORDER_SAMPLES, -30.0 + 90.0 * order + vis_slit_scale(order) * fitOffsetArcsec)
        isFaint = order in faintOrders
        centre = fitX + rng.normal(0.0, 6.0, VIS_ORDER_SAMPLES) if isFaint else vis_true_centre(order, wavelength)
        centre = centre + (0.0 if isFaint else rng.normal(0.0, 0.1, VIS_ORDER_SAMPLES))
        isDetected = rng.random(VIS_ORDER_SAMPLES) > (0.6 if isFaint else 0.0)
        frames.append(
            pd.DataFrame(
                {
                    "order": order,
                    "wavelength": wavelength,
                    "slit_position": 0.0,
                    "fit_x": fitX,
                    "fit_y": vis_row(order, wavelength),
                    "cont_x": np.where(isDetected, centre, np.nan),
                    "cont_y": np.where(isDetected, vis_row(order, wavelength), np.nan),
                    "gauss_amplitude": np.where(isDetected, 5.0, np.nan),
                    "gauss_stddev": np.where(
                        isDetected, VIS_SEEING_SIGMA_PX + rng.normal(0.0, 0.05, VIS_ORDER_SAMPLES), np.nan
                    ),
                    "gauss_mean": np.where(isDetected, 30.0, np.nan),
                    "pre-clipped": ~isDetected,
                    "centre_shift": np.where(isDetected, fitX - centre, np.nan),
                }
            )
        )
    return pd.concat(frames, ignore_index=True)
