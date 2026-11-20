"""The per-frame uncertainty map accounts for NDIT (DY-2317).

A raw NIR file is the mean of its NDIT integrations, so its variance is that of one integration over NDIT:
sigma^2 = (max(L, 0) + RON^2) / NDIT, with L in e- and RON the per-integration read noise. Frames with no NDIT card
and NDIT 1 frames keep the error that develop gives them, including every X-Shooter UVB and VIS frame.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import ccdproc
import numpy as np
import pytest
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData

from soxspipe.commonutils.keyword_lookup import keyword_lookup
from soxspipe.recipes import base_recipe
from tests.factories.prepared_frames import TEST_GAIN, TEST_RON
from tests.factories.prepared_frames import bad_pixel_map_file as _bad_pixel_map
from tests.factories.prepared_frames import prepare_frame_recipe as _prepare_recipe
from tests.factories.prepared_frames import raw_frame_file as _raw_frame

pytestmark = pytest.mark.unit

SHAPE = (300, 300)
READ_NOISE = 7.3
SKY_LEVEL = 30.0
Z_SCATTER_TOLERANCE = 0.02


def _recipe(log: Any, instrument: str = "soxs") -> base_recipe:
    recipe = base_recipe.__new__(base_recipe)
    recipe.log = log
    recipe.kw = keyword_lookup(log=log, settings={"instrument": instrument}).get
    recipe.detectorParams = {"ron": READ_NOISE * u.electron}
    return recipe


def _frame(
    recipe: base_recipe,
    data: np.ndarray,
    *,
    ndit: int | None = None,
    dprType: str = "OBJECT",
) -> CCDData:
    """A gain-corrected frame in e-, with the NDIT card written under the recipe's own keyword."""
    header = fits.Header()
    header[recipe.kw("DPR_TYPE")] = dprType
    header[recipe.kw("DPR_TECH")] = "ECHELLE,SLIT"
    header[recipe.kw("DPR_CATG")] = "SCIENCE"
    if ndit is not None:
        header[recipe.kw("DET_NDIT")] = ndit
    return CCDData(data.astype(np.float32), unit=u.electron, meta=header)


def _uniform(level: float) -> np.ndarray:
    return np.full(SHAPE, level, dtype=np.float32)


def _develop_error(recipe: base_recipe, data: np.ndarray) -> np.ndarray:
    """The error develop gives a frame: ccdproc's deviation with the detector read noise."""
    reference = CCDData(data.astype(np.float32), unit=u.electron)
    reference = ccdproc.create_deviation(reference, readnoise=recipe.detectorParams["ron"], disregard_nan=True)
    return reference.uncertainty.array.astype(np.float32)


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
def test_an_ndit_1_frame_and_a_frame_without_the_card_keep_the_develop_error(log: Any, instrument: str) -> None:
    # ARRANGE
    recipe = _recipe(log, instrument)
    data = np.random.default_rng(1).poisson(SKY_LEVEL, SHAPE).astype(np.float32)
    expected = _develop_error(recipe, data)

    # ACT
    withCard = recipe._add_uncertainty_map(_frame(recipe, data, ndit=1))
    withoutCard = recipe._add_uncertainty_map(_frame(recipe, data))

    # ASSERT
    np.testing.assert_array_equal(withCard.uncertainty.array, expected)
    np.testing.assert_array_equal(withoutCard.uncertainty.array, expected)
    assert withCard.uncertainty.array.dtype == np.float32


@pytest.mark.parametrize("instrument", ["soxs", "xsh"])
@pytest.mark.parametrize("ndit", [2, 4])
def test_the_error_of_an_ndit_frame_is_the_develop_error_over_the_square_root_of_ndit(
    log: Any, instrument: str, ndit: int
) -> None:
    # ARRANGE
    recipe = _recipe(log, instrument)
    data = np.random.default_rng(2).poisson(SKY_LEVEL, SHAPE).astype(np.float32)
    expected = _develop_error(recipe, data) / np.sqrt(ndit)

    # ACT
    prepared = recipe._add_uncertainty_map(_frame(recipe, data, ndit=ndit))

    # ASSERT
    np.testing.assert_allclose(prepared.uncertainty.array, expected, rtol=1e-6)
    assert prepared.uncertainty.array.dtype == np.float32


def test_the_error_of_an_ndit_4_frame_is_exactly_half_that_of_an_ndit_1_frame(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log)
    data = np.random.default_rng(3).poisson(SKY_LEVEL, SHAPE).astype(np.float32)

    # ACT
    single = recipe._add_uncertainty_map(_frame(recipe, data, ndit=1))
    quad = recipe._add_uncertainty_map(_frame(recipe, data, ndit=4))

    # ASSERT
    np.testing.assert_array_equal(quad.uncertainty.array, single.uncertainty.array / 2)


@pytest.mark.parametrize("ndit", [1, 2, 4])
def test_the_error_of_a_mean_of_ndit_integrations_has_unit_z_scatter(log: Any, ndit: int) -> None:
    # ARRANGE
    rng = np.random.default_rng(10 + ndit)
    recipe = _recipe(log)
    integrations = [rng.poisson(SKY_LEVEL, SHAPE) + rng.normal(0.0, READ_NOISE, SHAPE) for _ in range(ndit)]
    data = np.mean(integrations, axis=0)

    # ACT
    prepared = recipe._add_uncertainty_map(_frame(recipe, data, ndit=ndit))

    # ASSERT
    z = (prepared.data - SKY_LEVEL) / prepared.uncertainty.array
    assert float(np.std(z)) == pytest.approx(1.0, abs=Z_SCATTER_TOLERANCE)


def test_a_negative_level_gives_the_read_noise_alone_over_the_square_root_of_ndit(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log)
    data = np.full(SHAPE, -50.0, dtype=np.float32)

    # ACT
    prepared = recipe._add_uncertainty_map(_frame(recipe, data, ndit=2))

    # ASSERT
    expected = READ_NOISE / np.sqrt(2)
    np.testing.assert_allclose(prepared.uncertainty.array, expected, rtol=1e-6)


@pytest.mark.parametrize("ndit", [None, 1, 2, 4])
def test_the_bias_error_is_the_read_noise_over_the_square_root_of_ndit(log: Any, ndit: int | None) -> None:
    # ARRANGE
    recipe = _recipe(log)
    data = np.zeros(SHAPE, dtype=np.float32)

    # ACT
    prepared = recipe._add_uncertainty_map(_frame(recipe, data, ndit=ndit, dprType="BIAS"))

    # ASSERT
    expected = np.float32(READ_NOISE) / np.sqrt(np.float32(ndit or 1))
    np.testing.assert_allclose(prepared.uncertainty.array, expected, rtol=1e-6)
    assert prepared.uncertainty.array.dtype == np.float32


def test_an_invalid_ndit_card_leaves_the_develop_error_and_warns_once(log: Any) -> None:
    # ARRANGE
    recipe = _recipe(log)
    data = np.random.default_rng(4).poisson(SKY_LEVEL, SHAPE).astype(np.float32)
    frame = _frame(recipe, data)
    frame.header[recipe.kw("DET_NDIT")] = 0

    # ACT
    prepared = recipe._add_uncertainty_map(frame)

    # ASSERT
    np.testing.assert_array_equal(prepared.uncertainty.array, _develop_error(recipe, data))
    assert len([level for level, _ in log.messages if level == "warning"]) == 1


@pytest.mark.parametrize("value", [1e20, 10**30])
def test_an_out_of_range_ndit_card_leaves_the_develop_error_and_warns_once(log: Any, value: float) -> None:
    # ARRANGE
    recipe = _recipe(log)
    data = np.random.default_rng(5).poisson(SKY_LEVEL, SHAPE).astype(np.float32)
    frame = _frame(recipe, data)
    frame.header[recipe.kw("DET_NDIT")] = value

    # ACT
    prepared = recipe._add_uncertainty_map(frame)

    # ASSERT
    np.testing.assert_array_equal(prepared.uncertainty.array, _develop_error(recipe, data))
    assert len([level for level, _ in log.messages if level == "warning"]) == 1


def test_stacking_ndit_frames_divides_the_per_frame_error_by_root_n_without_scaling_it_again(
    log: Any, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    import soxspipe.commonutils.toolkit as toolkit

    recipe = _recipe(log)
    recipe.imageType = "OBJECT"
    recipe.verbose = False
    recipe.debug = False
    recipe.arm = "NIR"
    recipe.inst = "SOXS"
    recipe.recipeSettings = {"stacked-clipping-sigma": 50.0, "stacked-clipping-iterations": 1}
    monkeypatch.setattr(toolkit, "quicklook_image", lambda **_: None)
    data = _uniform(SKY_LEVEL)
    frames = []
    for _ in range(2):
        frame = recipe._add_uncertainty_map(_frame(recipe, data, ndit=2))
        frame.mask = np.zeros(SHAPE, dtype=bool)
        frame.header["EXPTIME"] = 300.0
        frames.append(frame)
    perFrame = frames[0].uncertainty.array.copy()

    # ACT
    stacked = recipe.clip_and_stack(frames, "soxs_stare", ignore_input_masks=False, post_stack_clipping=False)

    # ASSERT
    np.testing.assert_allclose(stacked.uncertainty.array, perFrame / np.sqrt(2), rtol=1e-5)


def test_a_prepared_nir_frame_with_ndit_2_is_written_with_the_scaled_error(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _prepare_recipe(log, tmp_path, arm="NIR")
    (tmp_path / "ndit2").mkdir()
    singleFrame = _raw_frame(tmp_path, dprType="LAMP,FLAT")
    doubleFrame = _raw_frame(tmp_path / "ndit2", dprType="LAMP,FLAT", DET_NDIT=2)

    # ACT
    single = CCDData.read(recipe._prepare_single_frame(singleFrame), hdu=0, unit=u.electron, hdu_uncertainty="ERRS")
    double = CCDData.read(recipe._prepare_single_frame(doubleFrame), hdu=0, unit=u.electron, hdu_uncertainty="ERRS")

    # ASSERT
    expected = np.sqrt(9.0 * TEST_GAIN + TEST_RON**2)
    assert single.uncertainty.array[1, 1] == pytest.approx(expected, rel=1e-6)
    assert double.uncertainty.array[1, 1] == pytest.approx(expected / np.sqrt(2), rel=1e-6)
    assert double.data[1, 1] == single.data[1, 1]


def test_preparing_an_already_prepared_ndit_2_frame_does_not_scale_the_error_again(log: Any, tmp_path: Path) -> None:
    # ARRANGE
    _bad_pixel_map(tmp_path / "calibrations")
    recipe = _prepare_recipe(log, tmp_path, arm="NIR")
    rawPath = _raw_frame(tmp_path, dprType="LAMP,FLAT", DET_NDIT=2)
    preparedPath = recipe._prepare_single_frame(rawPath)
    first = CCDData.read(preparedPath, hdu=0, unit=u.electron, hdu_uncertainty="ERRS")

    # ACT
    secondPath = recipe._prepare_single_frame(preparedPath)

    # ASSERT
    second = CCDData.read(secondPath, hdu=0, unit=u.electron, hdu_uncertainty="ERRS")
    assert secondPath == preparedPath
    np.testing.assert_array_equal(second.uncertainty.array, first.uncertainty.array)
    assert second.uncertainty.array[1, 1] == pytest.approx(
        np.sqrt(9.0 * TEST_GAIN + TEST_RON**2) / np.sqrt(2), rel=1e-6
    )
