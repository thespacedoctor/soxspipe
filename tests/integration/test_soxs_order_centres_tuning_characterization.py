"""Characterization of `soxs_order_centres.produce_product`: the tuning branch and the failure points.

Split from `test_soxs_order_centres_reduction_characterization.py` to keep
both files under 800 lines; the builders live there. The `tune-pipeline`
branch was the largest uncovered block at the branch point (lines 331-384).
The failure-point tests pin where a missing dispersion table and a failed
continuum detection raise, because a split could move either earlier.
"""

from __future__ import annotations

from importlib import import_module
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from soxspipe.recipes import soxs_order_centres
from tests.integration.test_mapping_recipe_orchestration import _route
from tests.integration.test_soxs_order_centres_reduction_characterization import (
    DETECTOR_KEYWORDS,
    ORDER_MODULE,
    _fitted_qc,
    _patch_reduction,
    _vis_recipe,
)

pytestmark = pytest.mark.integration

# THE FIVE DIGITS THE TUNING GRID PERMUTES, TWO AT A TIME.
TUNING_DIGITS = [2, 3, 4, 5, 6]
TUNING_PERMUTATIONS = len(TUNING_DIGITS) ** 2

# ---------------------------------------------------------------------------
# THE PARAMETER TUNING BRANCH
# ---------------------------------------------------------------------------


def _tuning_recipe(
    log: Any,
    tmpPath: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> tuple[soxs_order_centres, list[str], dict[str, Any], dict[str, Any]]:
    """Return a recipe configured for the `tune-pipeline` branch, with tuning stubbed."""
    recipe = _vis_recipe(log, tmpPath, settingOverrides={"tune-pipeline": True})
    calls, captured = _patch_reduction(
        recipe,
        monkeypatch,
        productPath=tmpPath / "ORDER_TAB_VIS.fits",
        fittedQc=_fitted_qc(),
        header={"WIN_BINX": 2, "WIN_BINY": 2},
    )

    tuningArguments: dict[str, Any] = {}

    def fake_fmultiprocess(**kwargs: object) -> list[object]:
        tuningArguments.update(kwargs)
        calls.append("fmultiprocess")
        return []

    fundamentals = import_module("fundamentals")
    monkeypatch.setattr(fundamentals, "fmultiprocess", fake_fmultiprocess)
    return recipe, calls, captured, tuningArguments


def test_tuning_samples_the_trace_once_then_sweeps_the_polynomial_grid(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The tuning branch samples the trace, then tunes over 25 permutations."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-run"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    (tuningDirectory / "residuals.txt").write_text("stale residuals\n")
    recipe, calls, captured, tuningArguments = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    returned = recipe.produce_product()

    # ASSERT
    assert returned is None
    assert calls == ["detrend", "keywords", "detect_continuum", "sample_trace", "fmultiprocess"]
    assert not (tuningDirectory / "residuals.txt").exists()

    detector = captured["detectors"][0]
    assert set(detector) == DETECTOR_KEYWORDS
    assert detector["traceFrame"] is recipe.orderFrame
    assert detector["dispersion_map"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert detector["recipeName"] == "soxs-order-centres"
    assert (detector["binx"], detector["biny"]) == (2, 2)

    assert tuningArguments["function"] is ORDER_MODULE.parameterTuning
    assert len(tuningArguments["inputArray"]) == TUNING_PERMUTATIONS
    assert tuningArguments["inputArray"][0] == (2, 2)
    assert tuningArguments["inputArray"][-1] == (6, 6)
    assert tuningArguments["poolSize"] is False
    assert tuningArguments["timeout"] == 360000
    assert tuningArguments["turnOffMP"] is recipe.debug
    assert tuningArguments["mute"] is True
    assert tuningArguments["progressBar"] is True
    assert tuningArguments["log"] is recipe.log
    assert tuningArguments["recipeSettings"] is recipe.recipeSettings
    assert tuningArguments["settings"] is recipe.settings
    assert tuningArguments["orderFrame"] is recipe.orderFrame
    assert tuningArguments["disp_map_table"] == str(tmp_path / "DISP_TAB_VIS.fits")
    assert tuningArguments["orderPixelTable"] == "ORDER-PIXEL-TABLE"
    assert tuningArguments["qc"] is recipe.qc
    assert tuningArguments["products"] is recipe.products
    assert tuningArguments["sofName"] == "synthetic-order-centres"
    assert (tuningArguments["binx"], tuningArguments["biny"]) == (2, 2)


def test_tuning_ignores_poly_orders_and_leaves_the_tables_untouched(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """No product row, no report, no clean-up, and the degrees are not overridden."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-untouched"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, calls, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)
    recipe.polyOrders = 35
    originalProducts = recipe.products.copy(deep=True)
    originalQc = recipe.qc.copy(deep=True)

    # ACT
    recipe.produce_product()

    # ASSERT
    pd.testing.assert_frame_equal(recipe.products, originalProducts)
    pd.testing.assert_frame_equal(recipe.qc, originalQc)
    assert "report" not in calls
    assert "clean_up" not in calls
    assert recipe.polyOrders == 35
    assert recipe.recipeSettings == {"detect-continuum": {}}


def test_tuning_without_a_stale_residuals_file_records_the_failed_removal(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A missing `residuals.txt` is swallowed and logged rather than raised."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-clean"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    debugged = [message for level, message in log.messages if level == "debug"]
    assert debugged[0] == "starting the ``produce_product`` method"
    assert debugged[1].startswith("produce_product: `os.remove('residuals.txt')` failed, continuing: ")
    assert len(debugged) == 2


def test_tuning_announces_itself_on_standard_output(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """The tuning banner is a bare `print`, not a logger call."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-banner"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    recipe, _, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)

    # ACT
    recipe.produce_product()

    # ASSERT
    assert "\n\nTUNING SOXSPIPE\n\n" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# FAILURE POINTS A SPLIT COULD MOVE
# ---------------------------------------------------------------------------


def test_a_missing_dispersion_table_fails_at_the_detector_after_the_degree_override(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The table name is bound only inside its loop, so it fails unbound, late."""
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path, polyOrders=35)
    del recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")]
    calls, _ = _patch_reduction(recipe, monkeypatch, productPath=tmp_path / "ORDER_TAB_VIS.fits", fittedQc=_fitted_qc())

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError, match="disp_map_table"):
        recipe.produce_product()

    assert calls == ["detrend", "keywords"]
    assert recipe.recipeSettings == {"detect-continuum": {"order-deg": 3, "disp-axis-deg": 5}}


def test_a_missing_dispersion_table_fails_after_tuning_clears_the_residuals(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """On the tuning branch the unbound table is read only after `residuals.txt` is removed."""
    # ARRANGE
    tuningDirectory = tmp_path / "tuning-no-table"
    tuningDirectory.mkdir()
    monkeypatch.chdir(tuningDirectory)
    (tuningDirectory / "residuals.txt").write_text("stale residuals\n")
    recipe, calls, _, _ = _tuning_recipe(log, tmp_path, monkeypatch)
    del recipe.inputFrames.filePathsByFilters[_route(PRO_CATG="DISP_TAB_VIS")]

    # ACT / ASSERT
    with pytest.raises(UnboundLocalError, match="disp_map_table"):
        recipe.produce_product()

    assert calls == ["detrend", "keywords"]
    assert not (tuningDirectory / "residuals.txt").exists()


def test_a_failed_continuum_detection_raises_a_named_recipe_failure(
    log: Any,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """A detector that returns no product path fails the recipe with its own message.

    The detector's quality-control rows are merged before the failure, so the
    evidence of the failed attempt survives in `recipe.qc`.
    """
    # ARRANGE
    recipe = _vis_recipe(log, tmp_path)
    originalQc = recipe.qc.copy(deep=True)
    calls, _ = _patch_reduction(recipe, monkeypatch, productPath=None, fittedQc=_fitted_qc())

    # ACT / ASSERT
    with pytest.raises(ArithmeticError) as raised:
        recipe.produce_product()

    assert str(raised.value) == (
        "Could not converge on a good fit to the VIS order-centre continuum. "
        "Please check the quality of your data or adjust your fitting parameters. "
        "No order table was produced."
    )
    assert calls == ["detrend", "keywords", "detect_continuum", "detector_get"]
    pd.testing.assert_frame_equal(recipe.qc, pd.concat([originalQc, _fitted_qc()]))
    assert "ORDER_CENTRES" not in list(recipe.products["product_label"])
