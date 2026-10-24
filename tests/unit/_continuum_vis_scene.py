"""Shared synthetic VIS scene for the continuum-trace sampling tests (DY-1716, DY-1941)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from soxspipe.commonutils.detect_continuum import detect_continuum


class TraceFrame(np.ma.MaskedArray):
    """A masked image that carries the empty `header` the sampler reads, one per instance."""

    header: dict

    def __array_finalize__(self, obj: object) -> None:
        super().__array_finalize__(obj)
        self.header = {}


VIS_ORDER_CENTRES = {1.0: 60.0, 2.0: 150.0, 3.0: 240.0, 4.0: 330.0}
VIS_ORDER_SAMPLE_COUNTS = {1.0: 30, 2.0: 50, 3.0: 100, 4.0: 30}
VIS_TRACE_SIGMA = 1.2


def vis_samples() -> pd.DataFrame:
    """Return the sample grid of four VIS orders, each one a vertical stripe at its own map position."""
    frames = []
    for order, count in VIS_ORDER_SAMPLE_COUNTS.items():
        frames.append(
            pd.DataFrame(
                {
                    "order": order,
                    "wavelength": np.linspace(400.0, 400.0 + count, count),
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
