"""Integration tests for real Astropy FITS and CCDData behavior."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData, StdDevUncertainty
from numpy.testing import assert_allclose, assert_array_equal

from tests.factories import prepared_fits, raw_fits, sof_file

pytestmark = pytest.mark.integration


def test_raw_fits_round_trip_uses_primary_hdu(tmp_path: Path) -> None:
    path = raw_fits(tmp_path / "raw.fits", shape=(32, 32), seed=3)

    frame = CCDData.read(path, unit=u.electron)

    assert frame.shape == (32, 32)
    assert frame.header["ESO SEQ ARM"] == "VIS"
    assert frame.uncertainty is None
    assert frame.mask is None


def test_prepared_fits_round_trip_preserves_flux_errors_quality_and_units(
    tmp_path: Path,
) -> None:
    path = prepared_fits(tmp_path / "prepared.fits", shape=(32, 32), seed=3)

    with fits.open(path) as hdus:
        assert [hdu.name for hdu in hdus] == ["FLUX", "QUAL", "ERRS"]
        assert hdus["FLUX"].header["BUNIT"] == "electron"
        assert hdus["ERRS"].header["UTYPE"] == "StdDevUncertainty"

    frame = CCDData.read(
        path,
        hdu=0,
        unit=u.electron,
        hdu_uncertainty="ERRS",
        hdu_mask="QUAL",
        key_uncertainty_type="UTYPE",
    )

    assert frame.shape == (32, 32)
    assert frame.unit == u.electron
    assert frame.header["ESO SEQ ARM"] == "VIS"
    assert isinstance(frame.uncertainty, StdDevUncertainty)
    assert_allclose(frame.uncertainty.array, np.full((32, 32), 1.5))
    expectedMask = np.zeros((32, 32), dtype=bool)
    expectedMask[0, 0] = True
    assert_array_equal(frame.mask, expectedMask)


def test_sof_factory_preserves_member_order(tmp_path: Path) -> None:
    firstPath = tmp_path / "first.fits"
    secondPath = tmp_path / "second.fits"

    path = sof_file(
        tmp_path / "frames.sof",
        [(secondPath, "BIAS_VIS"), (firstPath, "BIAS_VIS")],
    )

    assert path.read_text(encoding="utf-8").splitlines() == [
        f"{secondPath} BIAS_VIS",
        f"{firstPath} BIAS_VIS",
    ]
