"""The variance `soxs_nod.stack_extractions` reports for the per-wavelength median of N spectra (DY-1675).

`stack_extractions` combines N order-merged spectra (the nod per-cycle path, and every `soxs_offset` path) by
taking the median flux at each `WAVE`. The variance it reports is the variance of that median, not the median of
the input variances:

- N <= 2: the median is the mean, so the variance is sum(var) / N**2.
- N > 2: the large-N variance of the median of Gaussian samples, (pi / 2) * mean(var) / N.

N counts only the spectra that have a value at that `WAVE`. The expected values below are worked by hand.
"""

from __future__ import annotations

import math
from importlib import import_module
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from astropy import units as u

from soxspipe.recipes import soxs_nod
from tests.factories import pipeline_settings, qc_table, synthetic_ccd

pytestmark = pytest.mark.unit


@pytest.fixture
def recipe(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> soxs_nod:
    """Return an unconstructed recipe that can run `stack_extractions`, with the SNR QCs passed through."""
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"),
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    nod = soxs_nod.__new__(soxs_nod)
    nod.log = log
    nod.arm = "VIS"
    nod.recipeName = "soxs-nod"
    nod.sofName = "SYNTHETIC"
    nod.settings = pipeline_settings(tmp_path)
    nod.productDir = str(tmp_path)
    _, nod.products = nod._empty_qc_and_product_tables(pd)
    nod.qc = qc_table().iloc[0:0].copy()
    nod.masterHeaderFrame = synthetic_ccd(shape=(3, 3), prepared=True)
    nod.update_fits_keywords = lambda **_: None
    return nod


def _spectrum(waves: list[float], fluxes: list[float], variances: list[float]) -> pd.DataFrame:
    """Return an order-merged spectrum shaped like the extraction's output."""
    return pd.DataFrame(
        {
            "WAVE": [wave * u.nm for wave in waves],
            "FLUX_COUNTS": fluxes,
            "VARIANCE": variances,
            "SKY_COUNTS": [0.0] * len(waves),
            "SNR": [0.0] * len(waves),
            "FLUX_DENSITY_COUNTS": fluxes,
        }
    )


def _row(stacked: Any, wave: float) -> dict[str, float]:
    """Return the stacked row at `wave` as plain floats."""
    index = list(stacked["WAVE"]).index(wave)
    return {name: float(stacked[name][index]) for name in ("FLUX_COUNTS", "VARIANCE", "SNR")}


def test_a_single_spectrum_keeps_its_own_variance(recipe: soxs_nod) -> None:
    """With N = 1 there is nothing to combine, so the variance is not divided by N a second time."""
    stacked, _ = recipe.stack_extractions([_spectrum([500.0], [10.0], [4.0])])

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(10.0)
    assert row["VARIANCE"] == pytest.approx(4.0)
    assert row["SNR"] == pytest.approx(5.0)


def test_two_spectra_report_the_variance_of_their_mean(recipe: soxs_nod) -> None:
    """The median of two values is their mean: flux 15, variance (9 + 16) / 4 = 6.25, SNR 15 / 2.5 = 6."""
    first = _spectrum([500.0], [10.0], [9.0])
    second = _spectrum([500.0], [20.0], [16.0])

    stacked, _ = recipe.stack_extractions([first, second])

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(15.0)
    assert row["VARIANCE"] == pytest.approx(6.25)
    assert row["SNR"] == pytest.approx(6.0)


def test_three_spectra_report_the_variance_of_their_median(recipe: soxs_nod) -> None:
    """N = 3: flux median 2, variance (pi / 2) * mean(3, 4, 5) / 3 = 2 pi / 3, SNR 2 / sqrt(2 pi / 3)."""
    spectra = [
        _spectrum([500.0], [1.0], [3.0]),
        _spectrum([500.0], [2.0], [4.0]),
        _spectrum([500.0], [3.0], [5.0]),
    ]

    stacked, _ = recipe.stack_extractions(spectra)

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(2.0)
    assert row["VARIANCE"] == pytest.approx(2 * math.pi / 3)
    assert row["SNR"] == pytest.approx(round(2.0 / math.sqrt(2 * math.pi / 3), 2))


def test_four_cycles_concatenated_into_one_spectrum_count_as_four(recipe: soxs_nod) -> None:
    """The per-cycle paths concatenate every cycle into one dataframe; each row at a `WAVE` is one input.

    Variances 1, 2, 3, 6 have mean 3, so the variance is (pi / 2) * 3 / 4 = 3 pi / 8.
    """
    cycles = pd.concat([_spectrum([500.0], [flux], [variance]) for flux, variance in [(1, 1), (2, 2), (4, 3), (8, 6)]])

    stacked, _ = recipe.stack_extractions([cycles])

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(3.0)
    assert row["VARIANCE"] == pytest.approx(3 * math.pi / 8)


def test_a_wavelength_only_some_spectra_cover_counts_only_those(recipe: soxs_nod) -> None:
    """At 500 nm three spectra contribute (N = 3); at 510 nm only two do (N = 2); at 520 nm only one (N = 1)."""
    first = _spectrum([500.0, 510.0, 520.0], [1.0, 10.0, 7.0], [3.0, 9.0, 2.0])
    second = _spectrum([500.0, 510.0], [2.0, 20.0], [4.0, 16.0])
    third = _spectrum([500.0], [3.0], [5.0])

    stacked, _ = recipe.stack_extractions([first, second, third])

    assert _row(stacked, 500.0)["VARIANCE"] == pytest.approx(2 * math.pi / 3)
    assert _row(stacked, 510.0)["VARIANCE"] == pytest.approx(6.25)
    assert _row(stacked, 520.0)["VARIANCE"] == pytest.approx(2.0)


def test_a_spectrum_with_no_value_at_a_wavelength_does_not_count_there(recipe: soxs_nod) -> None:
    """A NaN flux is skipped by the median, so its row must not count towards N: N = 2 here, not 3."""
    first = _spectrum([500.0], [10.0], [9.0])
    second = _spectrum([500.0], [20.0], [16.0])
    gap = _spectrum([500.0], [np.nan], [np.nan])

    stacked, _ = recipe.stack_extractions([first, second, gap])

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(15.0)
    assert row["VARIANCE"] == pytest.approx(6.25)


def test_a_flux_without_a_variance_leaves_the_stacked_variance_unknown(recipe: soxs_nod) -> None:
    """The median uses every flux, so a contributing flux whose variance is NaN makes the variance NaN, not N - 1."""
    first = _spectrum([500.0], [10.0], [9.0])
    second = _spectrum([500.0], [20.0], [np.nan])

    stacked, _ = recipe.stack_extractions([first, second])

    row = _row(stacked, 500.0)
    assert row["FLUX_COUNTS"] == pytest.approx(15.0)
    assert math.isnan(row["VARIANCE"])


def test_a_wavelength_no_spectrum_covers_has_no_variance(recipe: soxs_nod) -> None:
    """With N = 0 the variance and SNR are NaN, never infinite."""
    first = _spectrum([500.0], [np.nan], [np.nan])
    second = _spectrum([500.0], [np.nan], [np.nan])

    stacked, _ = recipe.stack_extractions([first, second])

    row = _row(stacked, 500.0)
    assert math.isnan(row["VARIANCE"])
    assert math.isnan(row["SNR"])


def test_the_variance_column_keeps_the_extraction_float32_type(recipe: soxs_nod) -> None:
    """The extraction writes float32 variances; the stacked column keeps that type, as the median of it did."""
    first = _spectrum([500.0], [10.0], [9.0]).astype({"VARIANCE": np.float32})
    second = _spectrum([500.0], [20.0], [16.0]).astype({"VARIANCE": np.float32})

    stacked, _ = recipe.stack_extractions([first, second])

    assert stacked["VARIANCE"].dtype == np.float32
    assert _row(stacked, 500.0)["VARIANCE"] == pytest.approx(6.25)


def _extraction_spectrum(waves: list[float], fluxes: list[float], variances: list[float]) -> pd.DataFrame:
    """Return a spectrum built the way `horne_extraction` builds its merged orders: float32 astropy quantities."""
    spectrum = pd.DataFrame()
    spectrum["WAVE"] = np.array(waves) * u.nm
    spectrum["FLUX_COUNTS"] = np.array(fluxes, dtype=np.float32) * u.electron
    spectrum["VARIANCE"] = np.array(variances, dtype=np.float32) * u.electron
    spectrum["SKY_COUNTS"] = np.zeros(len(waves), dtype=np.float32) * u.electron
    spectrum["SNR"] = spectrum["FLUX_COUNTS"].values / np.sqrt(spectrum["VARIANCE"].values)
    spectrum["FLUX_DENSITY_COUNTS"] = np.array(fluxes, dtype=np.float32) * u.electron / u.nm
    return spectrum


def test_spectra_built_like_the_extraction_output_stack(recipe: soxs_nod) -> None:
    """Real extractions hold float32 quantities in mixed units; stacking two of them gives sum(var) / 4."""
    first = _extraction_spectrum([500.0, 510.0], [10.0, 30.0], [9.0, 4.0])
    second = _extraction_spectrum([500.0, 510.0], [20.0, 50.0], [16.0, 12.0])

    stacked, _ = recipe.stack_extractions([first, second])

    assert _row(stacked, 500.0)["VARIANCE"] == pytest.approx(6.25)
    assert _row(stacked, 510.0)["VARIANCE"] == pytest.approx(4.0)
    assert _row(stacked, 510.0)["FLUX_COUNTS"] == pytest.approx(40.0)
