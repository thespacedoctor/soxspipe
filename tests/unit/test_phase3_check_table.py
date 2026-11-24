"""Unit tests for the table checks of the Phase 3 checker (DY-2345)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

from soxspipe.commonutils.phase3_check import FAILURE, check_phase3
from tests.factories.phase3_products import (
    PRODUCT_NAME,
    SPECTRUM_LENGTH,
    failure_messages,
    finding_messages,
    phase3_spectrum_file,
    rewrite_product,
)

pytestmark = pytest.mark.unit


def test_the_first_three_columns_must_be_wave_flux_and_err_in_order(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, dropColumns=("ERR",))

    # ASSERT
    assert any("WAVE, FLUX, ERR" in message for message in failures)


def test_a_column_of_another_length_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columns={"FLUX_COUNTS": np.ones(SPECTRUM_LENGTH - 1)})

    # ASSERT
    assert any("FLUX_COUNTS" in message and "elements" in message for message in failures)


def test_nelem_that_differs_from_the_array_length_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, extension={"NELEM": SPECTRUM_LENGTH + 1})

    # ASSERT
    assert any("NELEM" in message for message in failures)


@pytest.mark.parametrize("wave", [np.array([600.0, 600.5, 600.5, 601.0, 601.5, 602.0]), np.arange(6.0)[::-1] + 600])
def test_wave_that_does_not_strictly_increase_is_a_failure(tmp_path: Path, log: object, wave: np.ndarray) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columns={"WAVE": wave})

    # ASSERT
    assert any("WAVE" in message and "strictly increas" in message for message in failures)


@pytest.mark.parametrize("value", [np.nan, np.inf])
def test_a_non_finite_value_in_a_column_is_a_failure(tmp_path: Path, log: object, value: float) -> None:
    # ARRANGE
    flux = np.full(SPECTRUM_LENGTH, 1.0e-16)
    flux[2] = value

    # ACT
    failures = failure_messages(tmp_path, log, columns={"FLUX": flux})

    # ASSERT
    assert any("FLUX" in message and "NaN" in message for message in failures)


def test_a_float_qual_column_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columns={"QUAL": np.zeros(SPECTRUM_LENGTH)})

    # ASSERT
    assert any("QUAL" in message and "integer" in message for message in failures)


def test_a_repeated_column_name_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE
    product = phase3_spectrum_file(tmp_path / PRODUCT_NAME, log)

    def repeat_name(hdus: fits.HDUList) -> None:
        hdus[1].header["TTYPE6"] = "FLUX_COUNTS"

    rewrite_product(product, repeat_name)

    # ACT
    failures = finding_messages(check_phase3(product), FAILURE)

    # ASSERT
    assert any("FLUX_COUNTS" in message and "more than once" in message for message in failures)


def test_a_wave_ucd_that_is_neither_air_nor_vacuum_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(
        tmp_path, log, columnCards={"WAVE": ("nm", "em.wl;obs.vacuum", "spec:Data.SpectralAxis.Value")}
    )

    # ASSERT
    assert any("TUCD1" in message and "em.wl" in message for message in failures)


def test_a_vacuum_wave_ucd_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columnCards={"WAVE": ("nm", "em.wl", "spec:Data.SpectralAxis.Value")})

    # ASSERT
    assert failures == []


@pytest.mark.parametrize(
    ("column", "cards", "expected"),
    [
        ("FLUX", ("erg / (s cm2 Angstrom)", "phot.flux.density;em.wl", "spec:Data.FluxAxis.Value"), "meta.main"),
        (
            "ERR",
            ("erg / (s cm2 Angstrom)", "stat.error;phot.flux.density;em.ql;meta.main", "spec:Data.FluxAxis.Value"),
            "TUTYP3",
        ),
        ("WAVE", ("nm", "em.wl;obs.atmos", "Spectrum.Data.SpectralAxis.Value"), "TUTYP1"),
    ],
)
def test_the_main_columns_need_the_v2_utypes_and_the_main_flux_ucds(
    tmp_path: Path, log: object, column: str, cards: tuple[str, str, str], expected: str
) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columnCards={column: cards})

    # ASSERT
    assert any(expected in message for message in failures)


@pytest.mark.parametrize(
    ("cards", "expected"),
    [
        ((None, "phot.flux.density;em.wl", "eso:Data.FluxAxis.Value"), "TUNIT5"),
        (("adu", None, "eso:Data.FluxAxis.Value"), "TUCD5"),
        (("adu", "phot.flux.density;em.wl", None), "TUTYP5"),
    ],
)
def test_a_column_without_unit_ucd_or_utype_is_a_failure(
    tmp_path: Path, log: object, cards: tuple[str | None, str | None, str | None], expected: str
) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, columnCards={"FLUX_COUNTS": cards})

    # ASSERT
    assert any(expected in message and "missing" in message for message in failures)


def test_err_in_another_unit_than_flux_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(
        tmp_path,
        log,
        columnCards={
            "ERR": ("adu", "stat.error;phot.flux.density;em.ql;meta.main", "spec:Data.FluxAxis.Accuracy.StatError")
        },
    )

    # ASSERT
    assert any("ERR" in message and "unit of FLUX" in message for message in failures)


def test_a_bgflux_column_must_match_the_flux_unit_and_ucd(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(
        tmp_path,
        log,
        columns={"BGFLUX": np.ones(SPECTRUM_LENGTH)},
        columnCards={"BGFLUX": ("adu", "phot.flux.density;em.wl", "eso:Data.FluxAxis.Value")},
    )

    # ASSERT
    assert any("BGFLUX" in message and "FLUX" in message for message in failures)


def test_a_missing_tdmin_or_tdmax_of_wave_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, extension={"TDMIN1": None})

    # ASSERT
    assert any("TDMIN1" in message for message in failures)


def test_wave_outside_tdmin_and_tdmax_is_a_failure(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(tmp_path, log, extension={"TDMAX1": 601.0})

    # ASSERT
    assert any("TDMAX1" in message and "WAVE" in message for message in failures)


def test_an_exposure_column_needs_its_ucd_and_an_empty_utype(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(
        tmp_path,
        log,
        columns={"EXPOSURE": np.full(SPECTRUM_LENGTH, 600.0)},
        columnCards={"EXPOSURE": ("s", "time.duration", "eso:Data.Exposure")},
    )

    # ASSERT
    assert any("EXPOSURE" in message and "time.duration;obs.exposure" in message for message in failures)
    assert any("EXPOSURE" in message and "empty" in message for message in failures)


def test_a_good_exposure_column_is_accepted(tmp_path: Path, log: object) -> None:
    # ARRANGE, ACT
    failures = failure_messages(
        tmp_path,
        log,
        columns={"EXPOSURE": np.full(SPECTRUM_LENGTH, 600.0)},
        columnCards={"EXPOSURE": ("s", "time.duration;obs.exposure", "")},
    )

    # ASSERT
    assert failures == []
