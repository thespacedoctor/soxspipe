"""*tests for the developer script that measures the log-wavelength step of each arm*

The measurement is pure: the tests feed it synthetic wavelength columns on a grid with a known velocity step, so
every expected value is exact to floating point.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
import yaml
from astropy.io import fits

TOOLS_DIR = Path(__file__).resolve().parents[2] / "tools"
sys.path.insert(0, str(TOOLS_DIR))

import measure_log_wavelength_step as step_tool  # noqa: E402

pytestmark = [
    pytest.mark.unit,
    # THE HEADER CARDS OF THE SYNTHETIC FRAMES ARE ESO HIERARCH CARDS, WHICH ASTROPY WARNS ABOUT WHEN IT CREATES THEM
    pytest.mark.filterwarnings("ignore:Keyword name:astropy.io.fits.verify.VerifyWarning"),
]

LIGHT_SPEED_KMS = 299792.458


def _log_grid(startNm: float, stepKms: float, size: int) -> np.ndarray:
    """*wavelengths with a constant velocity step, lambda_k = lambda_0 exp(k dv / c)*"""
    return startNm * np.exp(np.arange(size) * stepKms / LIGHT_SPEED_KMS)


def test_native_step_is_the_velocity_step_of_a_log_grid():
    wavelengths = _log_grid(550.0, 18.0, 200)

    assert step_tool.native_step_kms(wavelengths) == pytest.approx(18.0, rel=1e-9)


def test_native_step_is_the_median_so_one_jump_does_not_move_it():
    wavelengths = _log_grid(550.0, 18.0, 200)
    wavelengths[100:] *= 1.001

    assert step_tool.native_step_kms(wavelengths) == pytest.approx(18.0, rel=1e-6)


def test_native_step_ignores_row_order_and_repeated_wavelengths():
    wavelengths = _log_grid(550.0, 18.0, 50)
    shuffled = np.concatenate([wavelengths, wavelengths])
    np.random.default_rng(1).shuffle(shuffled)

    assert step_tool.native_step_kms(shuffled) == pytest.approx(18.0, rel=1e-9)


def test_native_step_ignores_non_finite_and_non_positive_wavelengths():
    wavelengths = np.concatenate([[np.nan, 0.0, -5.0], _log_grid(550.0, 18.0, 50)])

    assert step_tool.native_step_kms(wavelengths) == pytest.approx(18.0, rel=1e-9)


def test_native_step_raises_with_fewer_than_two_distinct_wavelengths():
    with pytest.raises(step_tool.ToolError, match="two distinct"):
        step_tool.native_step_kms([550.0, 550.0])


def test_per_order_steps_are_the_median_step_of_each_order():
    orders = np.repeat([1, 2], 100)
    wavelengths = np.concatenate([_log_grid(400.0, 20.0, 100), _log_grid(600.0, 18.0, 100)])

    steps = step_tool.per_order_steps_kms(orders, wavelengths)

    assert steps == {1: pytest.approx(20.0, rel=1e-9), 2: pytest.approx(18.0, rel=1e-9)}


def test_per_order_steps_skip_an_order_with_too_few_rows():
    orders = np.array([1] * 100 + [2])
    wavelengths = np.concatenate([_log_grid(400.0, 20.0, 100), [600.0]])

    assert list(step_tool.per_order_steps_kms(orders, wavelengths)) == [1]


def test_unbinned_step_is_the_smallest_order_step_per_unbinned_pixel():
    steps = {1: 20.0, 2: 18.0, 3: 19.5}

    assert step_tool.unbinned_step_kms(steps, binning=1) == 18.0


def test_unbinned_step_divides_the_binned_step_by_the_binning_factor():
    steps = {1: 40.0, 2: 36.0}

    assert step_tool.unbinned_step_kms(steps, binning=2) == 18.0


def test_unbinned_step_raises_when_no_order_was_measured():
    with pytest.raises(step_tool.ToolError, match="no order"):
        step_tool.unbinned_step_kms({}, binning=1)


def _write_extracted_orders(path: Path, arm: str, binning: int, mjd: float, orders, wavelengths) -> Path:
    """*write a table like a `_EXTRACTED_ORDERS` product: cards in the primary header, the rows in extension 1*"""
    primary = fits.PrimaryHDU()
    primary.header["ESO SEQ ARM"] = arm
    primary.header["ESO DET BINX"] = binning if arm == "VIS" else 1
    primary.header["ESO DET BINY"] = binning if arm == "NIR" else 1
    primary.header["MJD-OBS"] = mjd
    table = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="order", format="K", array=np.asarray(orders)),
            fits.Column(name="wavelengthMean", format="D", array=np.asarray(wavelengths)),
        ]
    )
    fits.HDUList([primary, table]).writeto(path)
    return path


def _two_order_table(path: Path, arm: str = "VIS", binning: int = 1, mjd: float = 60000.5, factor: float = 1.0) -> Path:
    orders = np.repeat([1, 2], 100)
    wavelengths = np.concatenate([_log_grid(400.0, 20.0 * factor, 100), _log_grid(600.0, 18.0 * factor, 100)])
    return _write_extracted_orders(path, arm, binning, mjd, orders, wavelengths)


def test_read_product_returns_the_arm_binning_mjd_and_columns(tmp_path):
    path = _two_order_table(tmp_path / "a_EXTRACTED_ORDERS.fits", arm="NIR", binning=1, mjd=60001.25)

    product = step_tool.read_product(path)

    assert (product.arm, product.binning, product.mjd) == ("NIR", 1, 60001.25)
    assert product.orders.size == product.wavelengths.size == 200


def test_measure_arms_pools_the_rows_of_every_file_of_an_arm(tmp_path):
    first = _two_order_table(tmp_path / "a.fits", mjd=60000.5)
    second = _two_order_table(tmp_path / "b.fits", mjd=60002.5)

    results = step_tool.measure_arms([step_tool.read_product(first), step_tool.read_product(second)])

    assert set(results) == {"VIS"}
    assert results["VIS"].stepKms == pytest.approx(18.0, rel=1e-9)
    assert results["VIS"].mjds == [60000.5, 60002.5]
    assert sorted(results["VIS"].perOrderKms) == [1, 2]


def test_measure_arms_does_not_mix_the_bins_of_different_frames(tmp_path):
    # THE BINS OF TWO FRAMES ARE HALF A STEP APART, SO POOLING THEIR WAVELENGTHS WOULD HALVE THE STEP
    orders = np.repeat([1], 100)
    first = _write_extracted_orders(tmp_path / "a.fits", "VIS", 1, 60000.5, orders, _log_grid(400.0, 18.0, 100))
    shifted = _log_grid(400.0, 18.0, 100) * np.exp(9.0 / LIGHT_SPEED_KMS)
    second = _write_extracted_orders(tmp_path / "b.fits", "VIS", 1, 60001.5, orders, shifted)

    results = step_tool.measure_arms([step_tool.read_product(first), step_tool.read_product(second)])

    assert results["VIS"].stepKms == pytest.approx(18.0, rel=1e-9)


def test_measure_arms_keeps_the_arms_apart(tmp_path):
    vis = _two_order_table(tmp_path / "vis.fits", arm="VIS")
    nir = _two_order_table(tmp_path / "nir.fits", arm="NIR", factor=0.6)

    results = step_tool.measure_arms([step_tool.read_product(vis), step_tool.read_product(nir)])

    assert results["VIS"].stepKms == pytest.approx(18.0, rel=1e-9)
    assert results["NIR"].stepKms == pytest.approx(10.8, rel=1e-9)


def test_measure_arms_divides_binned_frames_down_to_the_unbinned_step(tmp_path):
    binned = _two_order_table(tmp_path / "binned.fits", binning=2, factor=2.0)

    results = step_tool.measure_arms([step_tool.read_product(binned)])

    assert results["VIS"].stepKms == pytest.approx(18.0, rel=1e-9)


def test_measure_arms_refuses_to_mix_binnings_within_an_arm(tmp_path):
    unbinned = _two_order_table(tmp_path / "a.fits", binning=1)
    binned = _two_order_table(tmp_path / "b.fits", binning=2, factor=2.0)

    with pytest.raises(step_tool.ToolError, match="binning"):
        step_tool.measure_arms([step_tool.read_product(unbinned), step_tool.read_product(binned)])


def test_read_product_raises_when_the_table_lacks_the_columns(tmp_path):
    path = tmp_path / "bad.fits"
    primary = fits.PrimaryHDU()
    primary.header["ESO SEQ ARM"] = "VIS"
    table = fits.BinTableHDU.from_columns([fits.Column(name="other", format="D", array=np.zeros(3))])
    fits.HDUList([primary, table]).writeto(path)

    with pytest.raises(step_tool.ToolError, match="wavelengthMean"):
        step_tool.read_product(path)


def test_format_yaml_prints_the_provenance_and_the_step_per_arm(tmp_path):
    first = _two_order_table(tmp_path / "a.fits", mjd=60000.5)
    second = _two_order_table(tmp_path / "b.fits", mjd=60002.5)
    results = step_tool.measure_arms([step_tool.read_product(first), step_tool.read_product(second)])

    text = step_tool.format_yaml(results, version="0.19.1", today="2026-11-22")

    assert text == (
        "VIS:\n"
        "  # MEASURED 2026-11-22 WITH soxspipe 0.19.1 FROM 2 FRAMES ON 2 NIGHTS (MJD 60000 TO 60002)\n"
        "  # SMALLEST PER-ORDER MEDIAN STEP PER UNBINNED PIXEL; ORDER 1 20.00, ORDER 2 18.00 km/s\n"
        "  log-wavelength-step-kms: 18.0\n"
    )


def test_format_yaml_block_loads_back_as_the_measured_value(tmp_path):
    path = _two_order_table(tmp_path / "a.fits")
    results = step_tool.measure_arms([step_tool.read_product(path)])

    block = yaml.safe_load(step_tool.format_yaml(results, version="0.19.1", today="2026-11-22"))

    assert block == {"VIS": {"log-wavelength-step-kms": 18.0}}


def test_main_prints_the_yaml_block_for_the_files_it_is_given(tmp_path, capsys, monkeypatch):
    _two_order_table(tmp_path / "a_EXTRACTED_ORDERS.fits")
    monkeypatch.setattr(step_tool, "today_iso", lambda: "2026-11-22")

    status = step_tool.main([str(tmp_path / "*_EXTRACTED_ORDERS.fits")])

    output = capsys.readouterr().out
    assert status == 0
    assert yaml.safe_load(output) == {"VIS": {"log-wavelength-step-kms": 18.0}}


def test_main_raises_for_a_pattern_with_no_match(tmp_path):
    with pytest.raises(step_tool.ToolError, match="no file matches"):
        step_tool.main([str(tmp_path / "nothing*.fits")])
