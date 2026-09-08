"""Controlled polynomial contracts for dispersion-map conversion."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils import dispersion_map_to_pixel_arrays
from tests.factories import dispersion_map_fits

pytestmark = pytest.mark.unit


def _coordinates() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "order": [1.0, 2.0, -100.0],
            "wavelength": [2.0, 3.0, 1.0],
            "slit_position": [-0.5, 0.0, 0.5],
            "source": ["first", "second", "off detector"],
        }
    )


def test_conversion_evaluates_coefficients_and_preserves_columns(
    tmp_path: Path,
    log: object,
) -> None:
    dispersionPath = dispersion_map_fits(tmp_path / "dispersion.fits")

    converted = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath=str(dispersionPath),
        orderPixelTable=_coordinates(),
        removeOffDetectorLocation=False,
    )

    np.testing.assert_allclose(
        converted["fit_x"],
        [1.212, 1.326, 0.0],
        rtol=1e-6,
        atol=1e-7,
    )
    np.testing.assert_allclose(
        converted["fit_y"],
        [2.424, 2.652, 0.0],
        rtol=1e-6,
        atol=1e-7,
    )
    assert converted["fit_x"].dtype == np.float32
    assert converted["fit_y"].dtype == np.float32
    assert list(converted["source"]) == ["first", "second", "off detector"]


def test_conversion_filters_off_detector_rows_and_trims_columns(
    tmp_path: Path,
    log: object,
) -> None:
    dispersionPath = dispersion_map_fits(tmp_path / "dispersion.fits")

    converted = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath=str(dispersionPath),
        orderPixelTable=_coordinates(),
        removeOffDetectorLocation=True,
        trimColumns=True,
    )

    assert list(converted.columns) == [
        "order",
        "wavelength",
        "slit_position",
        "fit_x",
        "fit_y",
    ]
    assert list(converted["order"]) == [1.0, 2.0]
