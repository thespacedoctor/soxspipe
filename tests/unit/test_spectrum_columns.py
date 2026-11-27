"""The error and flux columns of a merged spectrum, in the SOXS layout and in the X-Shooter layout (DY-2350)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from astropy.table import Table

from soxspipe.commonutils.spectrum_columns import (
    calibrated_flux_column_name,
    error_column_name,
    error_counts,
    is_soxs_instrument,
    require_error_column,
    variance_counts,
)
from tests.factories import merged_spectrum

pytestmark = pytest.mark.unit

WAVES = [500.0, 600.0]
FLUX = [20.0, 40.0]
VARIANCE = [16.0, 64.0]


def test_the_error_column_of_a_soxs_spectrum_is_err_counts() -> None:
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE)

    # ACT / ASSERT
    assert error_column_name(spectrum.columns) == "ERR_COUNTS"


def test_the_error_column_of_an_xshooter_spectrum_is_variance() -> None:
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, hasVarianceColumn=True)

    # ACT / ASSERT
    assert error_column_name(spectrum.columns) == "VARIANCE"


def test_a_spectrum_with_neither_column_has_no_error_column_and_no_error_counts() -> None:
    # ARRANGE
    spectrum = pd.DataFrame({"WAVE": WAVES, "FLUX_COUNTS": FLUX})

    # ACT / ASSERT
    assert error_column_name(spectrum.columns) is None
    with pytest.raises(ValueError, match="ERR_COUNTS"):
        error_counts(spectrum)


@pytest.mark.parametrize("hasVarianceColumn", [False, True], ids=["soxs", "xshooter"])
@pytest.mark.parametrize("asTable", [False, True], ids=["pandas-dataframe", "astropy-table"])
def test_error_counts_and_variance_counts_agree_for_both_layouts(hasVarianceColumn: bool, asTable: bool) -> None:
    # ARRANGE
    spectrum = merged_spectrum(WAVES, FLUX, VARIANCE, asTable=asTable, hasVarianceColumn=hasVarianceColumn)

    # ACT
    error = error_counts(spectrum)
    variance = variance_counts(spectrum)

    # ASSERT: PLAIN FLOAT64 ARRAYS, THE ERROR IS THE ROOT OF THE VARIANCE
    assert error.dtype == np.float64
    assert variance.dtype == np.float64
    np.testing.assert_allclose(error, np.sqrt(VARIANCE), rtol=1e-6)
    np.testing.assert_allclose(variance, VARIANCE, rtol=1e-6)


def test_the_calibrated_flux_column_is_flux_for_soxs_and_flux_calibrated_for_xshooter() -> None:
    # ACT / ASSERT
    assert calibrated_flux_column_name(["WAVE", "FLUX", "ERR", "FLUX_COUNTS", "ERR_COUNTS", "SNR"]) == "FLUX"
    assert (
        calibrated_flux_column_name(Table(names=["WAVE", "FLUX_CALIBRATED", "VARIANCE"]).colnames) == "FLUX_CALIBRATED"
    )


@pytest.mark.parametrize(
    ("instrument", "isSoxs"),
    [("SOXS", True), ("soxs", True), ("Soxs", True), ("XSHOOTER", False), ("xsh", False), ("", False), (None, False)],
)
def test_the_instrument_name_selects_the_soxs_layout_without_regard_to_case(instrument: object, isSoxs: bool) -> None:
    # ACT / ASSERT
    assert is_soxs_instrument(instrument) is isSoxs


@pytest.mark.parametrize(
    ("isSoxs", "columns", "expected"),
    [(True, ["WAVE", "ERR_COUNTS"], "ERR_COUNTS"), (False, ["WAVE", "VARIANCE"], "VARIANCE")],
    ids=["soxs", "xshooter"],
)
def test_the_error_column_the_instrument_needs_is_returned_when_the_spectrum_has_it(
    isSoxs: bool, columns: list[str], expected: str
) -> None:
    # ACT / ASSERT
    assert require_error_column(columns, isSoxs) == expected


@pytest.mark.parametrize(
    ("isSoxs", "columns", "message"),
    [
        (True, ["WAVE", "VARIANCE"], "holds VARIANCE, but a SOXS spectrum holds ERR_COUNTS"),
        (False, ["WAVE", "ERR_COUNTS"], "holds ERR_COUNTS, but a X-Shooter spectrum holds VARIANCE"),
        (True, ["WAVE"], "missing the column ERR_COUNTS"),
        (False, ["WAVE"], "missing the column VARIANCE"),
    ],
    ids=["soxs-stale", "xshooter-soxs-columns", "soxs-none", "xshooter-none"],
)
def test_a_spectrum_of_the_other_layout_or_with_no_error_column_is_refused(
    isSoxs: bool, columns: list[str], message: str
) -> None:
    # ACT / ASSERT
    with pytest.raises(ValueError, match=message):
        require_error_column(columns, isSoxs)
