"""The arc-line resolution fit is quiet and keeps its limits (DY-2409).

`_calculate_resolution_on_slit` fits a Gaussian across one arc line and returns the resolution and the FWHM, or
``[None, None]`` when the fitted width is outside what the slit allows. It used to print a debug line to stdout for
every line it accepted.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from soxspipe.commonutils.create_dispersion_map import create_dispersion_map
from tests.unit.test_create_dispersion_map_helpers import _mapper

pytestmark = pytest.mark.unit


def _resolution_mapper(log: object, profile: np.ndarray) -> tuple[create_dispersion_map, pd.Series]:
    """Return a mapper and an arc-line row whose 14-sample cut across the line is ``profile``."""
    mapper = _mapper(log)
    mapper.slitWidth = 1.0
    mapper.detectorParams = {"dispersion-axis": "x"}
    image = np.zeros((40, 40))
    # THE CUT RUNS FROM ROW 13 TO ROW 26, SO ROW 20 IS SAMPLE 7
    image[13:27, 20] = profile
    mapper.arcFrame = SimpleNamespace(data=image)
    row = pd.Series({"observed_x": 20.0, "observed_y": 20.0, "wavelength": 500.0, "pixelScaleNm": 0.01})
    return mapper, row


def _gaussian_cut(sigma: float) -> np.ndarray:
    """Return a 14-sample Gaussian profile centred on sample 7."""
    return 100.0 * np.exp(-0.5 * ((np.arange(14) - 7.0) / sigma) ** 2)


def test_the_arc_line_resolution_fit_prints_nothing_to_stdout(log: object, capsys: pytest.CaptureFixture[str]) -> None:
    # ARRANGE
    mapper, row = _resolution_mapper(log, _gaussian_cut(1.5))

    # ACT
    result = mapper._calculate_resolution_on_slit(row)

    # ASSERT
    fwhm = result.iloc[1]
    assert 2.0 < fwhm < 5.0
    assert result.iloc[0] == pytest.approx(500.0 / (0.01 * fwhm))
    assert capsys.readouterr().out == ""


def test_an_arc_line_narrower_than_the_slit_allows_returns_no_resolution(log: object) -> None:
    # ARRANGE
    mapper, row = _resolution_mapper(log, _gaussian_cut(0.3))

    # ACT
    result = mapper._calculate_resolution_on_slit(row)

    # ASSERT
    assert result.tolist() == [None, None]


def test_two_arc_lines_in_one_cut_return_no_resolution_because_the_fit_is_too_wide(log: object) -> None:
    # ARRANGE
    mapper, row = _resolution_mapper(log, np.where(np.isin(np.arange(14), [2, 12]), 100.0, 0.0))

    # ACT
    result = mapper._calculate_resolution_on_slit(row)

    # ASSERT
    assert result.tolist() == [None, None]


def test_a_fit_whose_width_cannot_be_read_returns_no_resolution(log: object, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    from astropy.modeling import fitting

    mapper, row = _resolution_mapper(log, _gaussian_cut(1.5))
    brokenModel = SimpleNamespace(stddev=SimpleNamespace(value="not a number"))
    monkeypatch.setattr(fitting, "LevMarLSQFitter", lambda: lambda *args: brokenModel)

    # ACT
    result = mapper._calculate_resolution_on_slit(row)

    # ASSERT
    assert result.tolist() == [None, None]
