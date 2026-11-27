"""Stacks group the bins of SOXS spectra by their integer index on the log-wavelength grid (DY-2349).

`stack_extractions` combines the order-merged spectra of nod cycles, offset cycles and the nod positions of a
stare-mode nodding SOF (the stare `--no-stack` grouping). Every SOXS frame sits on the one grid of its arm, so the
bin index is the exact key: it neither splits one bin in two nor merges two bins, and the stacked `WAVE` is the grid
wavelength of the index. X-Shooter keeps its linear grids and its 4-decimal grouping.
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
from astropy.table import Table

from soxspipe.commonutils.log_wavelength_grid import grid_index, grid_step_kms, grid_wavelength
from soxspipe.recipes import soxs_nod, soxs_stare
from soxspipe.recipes.extraction_stacking import median_stack_variance
from tests.factories import pipeline_settings, qc_table, synthetic_ccd

pytestmark = pytest.mark.unit

STEP_PER_PIXEL_KMS = 17.83
BINNING = 2
STEP_KMS = grid_step_kms({"log-wavelength-step-kms": STEP_PER_PIXEL_KMS}, BINNING)
FIRST_BIN = 612


def _recipe(recipeClass: type, log: Any, tmpPath: Path, monkeypatch: pytest.MonkeyPatch, instrument: str) -> Any:
    """Return an unconstructed recipe that can run `stack_extractions`, with the SNR QCs passed through."""
    monkeypatch.setattr(
        import_module("soxspipe.commonutils.toolkit"),
        "add_snr_efficiency_qcs",
        lambda **kwargs: kwargs["qcTable"],
    )
    recipe = recipeClass.__new__(recipeClass)
    recipe.log = log
    recipe.arm = "VIS"
    recipe.recipeName = "soxs-stare" if recipeClass is soxs_stare else "soxs-nod"
    recipe.sofName = "SYNTHETIC"
    recipe.settings = pipeline_settings(tmpPath, instrument=instrument)
    recipe.productDir = str(tmpPath)
    _, recipe.products = recipe._empty_qc_and_product_tables(pd)
    recipe.qc = qc_table().iloc[0:0].copy()
    recipe.masterHeaderFrame = synthetic_ccd(shape=(3, 3), prepared=True)
    recipe.masterHeaderFrame.header["ESO DET BINY"] = BINNING
    recipe.detectorParams = {"log-wavelength-step-kms": STEP_PER_PIXEL_KMS}
    recipe.update_fits_keywords = lambda **_: None
    return recipe


@pytest.fixture
def nod(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> soxs_nod:
    return _recipe(soxs_nod, log, tmp_path, monkeypatch, "soxs")


@pytest.fixture
def stare(log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> soxs_stare:
    return _recipe(soxs_stare, log, tmp_path, monkeypatch, "soxs")


def _spectrum(waves: Any, fluxes: Any, variances: Any = None) -> pd.DataFrame:
    """Return an order-merged spectrum shaped like the extraction's output."""
    waves = np.asarray(waves, dtype=float)
    fluxes = np.asarray(fluxes, dtype=float)
    variances = np.full(len(waves), 4.0) if variances is None else np.asarray(variances, dtype=float)
    return pd.DataFrame(
        {
            "WAVE": [wave * u.nm for wave in waves],
            "FLUX_COUNTS": fluxes,
            "VARIANCE": variances,
            "SKY_COUNTS": np.zeros(len(waves)),
            "SNR": np.zeros(len(waves)),
            "FLUX_DENSITY_COUNTS": fluxes,
        }
    )


def _bins(first: int, count: int) -> np.ndarray:
    return np.arange(first, first + count)


def _grid(first: int, count: int) -> np.ndarray:
    return grid_wavelength(_bins(first, count), STEP_KMS)


def _four_decimal_boundary_around(binIndex: int) -> tuple[float, float]:
    """Return two wavelengths of bin `binIndex` that the old 4-decimal rounding put in different groups."""
    centre = float(grid_wavelength(binIndex, STEP_KMS))
    boundary = (math.floor(centre * 1e4) + 0.5) / 1e4
    below, above = boundary - 1e-6, boundary + 1e-6
    assert round(below, 4) != round(above, 4)
    assert grid_index(np.array([below, above]), STEP_KMS).tolist() == [binIndex, binIndex]
    return below, above


def test_two_spectra_of_one_bin_stack_to_one_row_even_when_four_decimal_rounding_would_split_them(
    nod: soxs_nod,
) -> None:
    # ARRANGE
    below, above = _four_decimal_boundary_around(FIRST_BIN)

    # ACT
    stacked, _ = nod.stack_extractions([_spectrum([below], [10.0]), _spectrum([above], [20.0])])

    # ASSERT
    assert len(stacked) == 1
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([15.0])


def test_the_stacked_wave_is_the_grid_wavelength_of_the_bin_index_exactly(nod: soxs_nod) -> None:
    # ARRANGE: EACH FRAME HOLDS THE WAVELENGTHS A LITTLE OFF THE GRID, AS THE 4-DECIMAL ROUNDING LEFT THEM
    onGrid = _grid(FIRST_BIN, 40)
    first = _spectrum(onGrid + 4e-5, np.arange(40.0))
    second = _spectrum(onGrid - 3e-5, np.arange(40.0))

    # ACT
    stacked, _ = nod.stack_extractions([first, second])

    # ASSERT
    assert np.array_equal(np.asarray(stacked["WAVE"]), onGrid)


def test_the_written_file_holds_the_exact_grid_wavelengths(nod: soxs_nod) -> None:
    # ARRANGE
    onGrid = _grid(FIRST_BIN, 12)

    # ACT
    _, filePath = nod.stack_extractions([_spectrum(onGrid + 2e-5, np.ones(12)), _spectrum(onGrid, np.ones(12))])

    # ASSERT
    assert np.array_equal(np.asarray(Table.read(filePath)["WAVE"]), onGrid)


def test_no_bin_is_lost_or_doubled_when_the_frames_cover_different_ranges(nod: soxs_nod) -> None:
    # ARRANGE: THREE FRAMES THAT START AND END AT DIFFERENT BINS, WITH A JITTER THAT CROSSES A ROUNDING BOUNDARY
    _, above = _four_decimal_boundary_around(FIRST_BIN + 4)
    thirdWaves = np.concatenate([_grid(FIRST_BIN + 2, 2), [above], _grid(FIRST_BIN + 5, 11)])
    first = _spectrum(_grid(FIRST_BIN, 10), np.full(10, 10.0))
    second = _spectrum(_grid(FIRST_BIN + 4, 10) + 4e-5, np.full(10, 20.0))
    third = _spectrum(thirdWaves, np.full(14, 30.0))

    # ACT
    stacked, _ = nod.stack_extractions([first, second, third])

    # ASSERT: EVERY BIN OF THE UNION IS THERE ONCE
    expected = _bins(FIRST_BIN, 16)
    assert grid_index(np.asarray(stacked["WAVE"]), STEP_KMS).tolist() == expected.tolist()


def test_each_bin_is_stacked_from_exactly_the_frames_that_cover_it(nod: soxs_nod) -> None:
    # ARRANGE: FRAME A COVERS BINS 0-5 (FLUX 10), FRAME B COVERS BINS 3-8 (FLUX 20)
    first = _spectrum(_grid(FIRST_BIN, 6), np.full(6, 10.0), np.full(6, 9.0))
    second = _spectrum(_grid(FIRST_BIN + 3, 6) + 4e-5, np.full(6, 20.0), np.full(6, 16.0))

    # ACT
    stacked, _ = nod.stack_extractions([first, second])

    # ASSERT: ONE FRAME GIVES ITS OWN VALUE, TWO FRAMES GIVE THE MEAN WITH THE VARIANCE OF THE MEAN
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([10.0] * 3 + [15.0] * 3 + [20.0] * 3)
    assert list(stacked["VARIANCE"]) == pytest.approx([9.0] * 3 + [6.25] * 3 + [16.0] * 3)


def test_three_frames_that_disagree_by_a_rounding_step_still_stack_as_three(nod: soxs_nod) -> None:
    # ARRANGE
    below, above = _four_decimal_boundary_around(FIRST_BIN)
    centre = float(grid_wavelength(FIRST_BIN, STEP_KMS))
    waves, fluxes, variances = [below, centre, above], [1.0, 2.0, 3.0], [3.0, 4.0, 5.0]
    spectra = [
        _spectrum([wave], [flux], [variance]) for wave, flux, variance in zip(waves, fluxes, variances, strict=True)
    ]

    # ACT
    stacked, _ = nod.stack_extractions(spectra)

    # ASSERT: N = 3, SO THE MEDIAN WITH THE VARIANCE (PI / 2) * MEAN(VAR) / N
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([2.0])
    assert list(stacked["VARIANCE"]) == pytest.approx([2 * math.pi / 3])


def test_a_stare_recipe_stacks_its_nod_positions_by_bin_index(stare: soxs_stare) -> None:
    # ARRANGE
    below, above = _four_decimal_boundary_around(FIRST_BIN)

    # ACT
    stacked, _ = stare.stack_extractions([_spectrum([below], [10.0]), _spectrum([above], [20.0])])

    # ASSERT
    assert len(stacked) == 1
    assert float(stacked["WAVE"][0]) == float(grid_wavelength(FIRST_BIN, STEP_KMS))


@pytest.mark.parametrize("instrument", ["soxs", "SOXS", "Soxs"])
def test_the_instrument_name_is_matched_without_regard_to_case(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, instrument: str
) -> None:
    # ARRANGE: THE BIN MEDIAN IS CALLED DIRECTLY, BECAUSE THE STACK HEADER OPENS <instrument>_detector_parameters.yaml,
    # WHICH ONLY THE LOWER-CASE NAME OPENS ON A CASE-SENSITIVE FILE SYSTEM
    recipe = _recipe(soxs_nod, log, tmp_path, monkeypatch, instrument)
    below, above = _four_decimal_boundary_around(FIRST_BIN)
    frame = pd.DataFrame({"WAVE": [below, above], "FLUX_COUNTS": [10.0, 20.0], "VARIANCE": [4.0, 4.0]})

    # ACT
    medians = recipe._median_by_bin(frame)

    # ASSERT: ONE BIN ON THE GRID, AS FOR "soxs"
    assert list(medians["WAVE"]) == [float(grid_wavelength(FIRST_BIN, STEP_KMS))]


@pytest.mark.parametrize("instrument", ["xsh", "XSH", "Xsh"])
def test_x_shooter_keeps_the_four_decimal_wavelength_grouping(
    log: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, instrument: str
) -> None:
    # ARRANGE: X-SHOOTER HAS NO LOG GRID STEP AND NO BINNING CARD, SO NOTHING OF THE GRID MAY BE READ
    recipe = _recipe(soxs_nod, log, tmp_path, monkeypatch, instrument)
    recipe.detectorParams = {}
    del recipe.masterHeaderFrame.header["ESO DET BINY"]
    first = _spectrum([500.00001, 500.2], [10.0, 1.0])
    second = _spectrum([500.00002, 500.2001], [20.0, 2.0])

    # ACT
    stacked, _ = recipe.stack_extractions([first, second])

    # ASSERT: 500.00001 AND 500.00002 ROUND TO 500.0, BUT 500.2 AND 500.2001 STAY APART (THE FILE'S 2-DECIMAL WAVE IS
    # THE X-SHOOTER OUTPUT ROUNDING, NOT PART OF THE GROUPING)
    assert len(stacked) == 3
    assert list(stacked["FLUX_COUNTS"]) == pytest.approx([15.0, 1.0, 2.0])


def test_median_stack_variance_counts_the_rows_of_each_group_key() -> None:
    # ARRANGE: THREE ROWS IN BIN 5 (N = 3), TWO IN BIN 6 (N = 2), ONE IN BIN 7 (N = 1)
    spectra = _spectrum([500.0] * 6, [1.0] * 6, [3.0, 4.0, 5.0, 9.0, 16.0, 2.0])
    groupKeys = np.array([5, 5, 5, 6, 6, 7])

    # ACT
    variance = median_stack_variance(spectra, groupKeys)

    # ASSERT
    assert list(variance.index) == [5, 6, 7]
    assert list(variance) == pytest.approx([2 * math.pi / 3, 6.25, 2.0])


def test_median_stack_variance_leaves_a_row_with_no_flux_out_of_its_group() -> None:
    # ARRANGE
    spectra = _spectrum([500.0] * 3, [1.0, np.nan, 3.0], [3.0, 100.0, 5.0])

    # ACT
    variance = median_stack_variance(spectra, np.array([5, 5, 5]))

    # ASSERT: N = 2, SO THE VARIANCE OF THE MEAN IS (3 + 5) / 4
    assert list(variance) == pytest.approx([2.0])


def test_a_wavelength_that_is_not_on_the_grid_of_the_stack_is_refused(nod: soxs_nod) -> None:
    # ARRANGE: THE SPECTRA WERE MERGED ON THE GRID OF 1X BINNING, BUT THE MASTER HEADER GIVES 2X, A GRID TWICE AS COARSE
    unbinned = grid_step_kms({"log-wavelength-step-kms": STEP_PER_PIXEL_KMS}, 1)
    waves = grid_wavelength(_bins(FIRST_BIN * BINNING + 1, 6), unbinned)

    # ACT / ASSERT
    with pytest.raises(ValueError, match="not on the 35.660 km/s log-wavelength grid"):
        nod.stack_extractions([_spectrum(waves, np.ones(6))])


def test_a_wavelength_half_a_bin_off_the_grid_is_refused_but_a_rounding_step_is_not(nod: soxs_nod) -> None:
    # ARRANGE
    centre = float(grid_wavelength(FIRST_BIN, STEP_KMS))
    halfBin = centre * (math.exp(0.5 * STEP_KMS / 299792.458) - 1.0)

    # ACT
    stacked, _ = nod.stack_extractions([_spectrum([centre + 5e-5], [1.0])])

    # ASSERT
    assert len(stacked) == 1
    with pytest.raises(ValueError, match="log-wavelength grid"):
        nod.stack_extractions([_spectrum([centre + halfBin], [1.0])])


def test_the_bin_median_leaves_the_frame_it_was_given_unchanged(nod: soxs_nod) -> None:
    # ARRANGE: THE CONCATENATED FRAME, WITH WAVE AS PLAIN FLOATS JUST OFF ITS GRID POINTS, AS THE STACK HANDS IT OVER
    waves = np.asarray(_grid(FIRST_BIN, 3)) * (1.0 + 1.0e-7)
    frame = pd.DataFrame({"WAVE": waves, "FLUX_COUNTS": [1.0, 2.0, 3.0], "VARIANCE": [4.0, 4.0, 4.0]})
    columnsBefore = list(frame.columns)
    wavesBefore = frame["WAVE"].to_numpy().copy()

    # ACT
    medians = nod._median_by_bin(frame)

    # ASSERT: THE RESULT HOLDS THE GRID WAVELENGTHS, BUT THE INPUT FRAME KEEPS ITS OWN
    assert list(medians["WAVE"]) == list(grid_wavelength(np.arange(FIRST_BIN, FIRST_BIN + 3), STEP_KMS))
    assert list(frame.columns) == columnsBefore
    assert np.array_equal(frame["WAVE"].to_numpy(), wavesBefore)


def test_a_wavelength_just_beyond_the_grid_tolerance_is_refused(nod: soxs_nod) -> None:
    # ARRANGE
    centre = float(grid_wavelength(FIRST_BIN, STEP_KMS))

    # ACT
    stacked, _ = nod.stack_extractions([_spectrum([centre * (1.0 + 0.5e-6)], [1.0])])

    # ASSERT: HALF THE TOLERANCE STACKS, TWICE THE TOLERANCE IS REFUSED
    assert len(stacked) == 1
    with pytest.raises(ValueError, match="log-wavelength grid"):
        nod.stack_extractions([_spectrum([centre * (1.0 + 2.0e-6)], [1.0])])
