"""The `WAVE RMS` and `WAVE NLINES` QCs of a dispersion solution (DY-2358).

They are the wavelength accuracy of the arc solution that the Phase 3 keywords `LAMRMS` and `LAMNLIN` come from. The RMS
is over every detection that survived the clipping, along the dispersion axis only, converted to nm with the pixel scale
of each line. The line count is the number of distinct (order, wavelength) pairs: one arc line is detected once per
pinhole, and those detections are not independent measurements (DY-2307).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from tests.unit.test_create_dispersion_map_helpers import _mapper

pytestmark = pytest.mark.unit

# THE LINE TABLE FITS x = 2 * wavelength AND y = 5, SO THE PIXEL SCALE OF EVERY LINE IS 4 nm / 8 px = 0.5 nm/px
PIXEL_SCALE_NM = 0.5


def _line_table() -> pd.DataFrame:
    """Return three detections of two arc lines: the 20 nm line is seen through two pinholes."""
    return pd.DataFrame(
        {
            "order": [1.0, 1.0, 1.0],
            "wavelength": [10.0, 20.0, 20.0],
            "slit_position": [0.0, 0.0, 0.0],
            "order_pow_x_0": 1.0,
            "wavelength_pow_x_0": 1.0,
            "wavelength_pow_x_1": [10.0, 20.0, 20.0],
            "slit_position_pow_x_0": 1.0,
            "order_pow_y_0": 1.0,
            "wavelength_pow_y_0": 1.0,
            "wavelength_pow_y_1": [10.0, 20.0, 20.0],
            "slit_position_pow_y_0": 1.0,
            # x RESIDUALS (FIT - OBSERVED) ARE 3, 2 AND -2 px; y RESIDUALS ARE 4, 0 AND 0 px
            "observed_x": [17.0, 38.0, 42.0],
            "observed_y": [1.0, 5.0, 5.0],
            "detector_x": [20.0, 40.0, 40.0],
            "detector_y": [5.0, 5.0, 5.0],
            "fwhm_pin_px": [2.0, 2.0, 2.0],
        }
    )


def _written_qcs(log: object, dispersionAxis: str) -> pd.DataFrame:
    mapper = _mapper(log)
    mapper.arcFrame = False
    mapper.axisA = dispersionAxis
    mapper.calculate_residuals(
        _line_table(),
        xcoeff=[0.0, 2.0],
        ycoeff=[5.0, 0.0],
        orderDeg=0,
        wavelengthDeg=1,
        slitDeg=0,
        pixelRange=True,
        writeQCs=True,
    )
    return mapper.qc.set_index("qc_name")


def test_wave_rms_is_the_rms_of_the_x_residuals_in_nm_when_the_dispersion_axis_is_x(log: object) -> None:
    # ACT
    qc = _written_qcs(log, "x")

    # ASSERT
    expected = PIXEL_SCALE_NM * np.sqrt((3.0**2 + 2.0**2 + 2.0**2) / 3)
    assert qc.loc["WAVE RMS", "qc_value"] == pytest.approx(expected, abs=1e-5)
    assert qc.loc["WAVE RMS", "qc_unit"] == "nm"
    assert bool(qc.loc["WAVE RMS", "to_header"])


def test_wave_rms_uses_the_y_residuals_when_the_dispersion_axis_is_y(log: object) -> None:
    # ACT
    qc = _written_qcs(log, "y")

    # ASSERT
    expected = PIXEL_SCALE_NM * np.sqrt((4.0**2) / 3)
    assert qc.loc["WAVE RMS", "qc_value"] == pytest.approx(expected, abs=1e-5)


def test_wave_nlines_counts_each_order_and_wavelength_once_however_many_pinholes_saw_it(log: object) -> None:
    # ACT
    qc = _written_qcs(log, "x")

    # ASSERT: THREE DETECTIONS, TWO ARC LINES
    assert qc.loc["WAVE NLINES", "qc_value"] == 2
    assert qc.loc["WAVE NLINES", "qc_unit"] == "lines"
    assert bool(qc.loc["WAVE NLINES", "to_header"])
