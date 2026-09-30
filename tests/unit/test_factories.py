"""Contract tests for fresh, deterministic synthetic factories."""

from __future__ import annotations

from pathlib import Path

import pytest
from numpy.testing import assert_array_equal

from tests.factories import (
    dispersion_table,
    order_table,
    pipeline_settings,
    product_table,
    qc_table,
    synthetic_ccd,
    synthetic_signal,
)

pytestmark = pytest.mark.unit


def test_settings_factory_returns_fresh_nested_state(tmp_path: Path) -> None:
    first = pipeline_settings(tmp_path)
    second = pipeline_settings(tmp_path)

    updated = {**first, "instrument": "xsh"}

    assert first is not second
    assert updated["instrument"] == "xsh"
    assert first["instrument"] == "soxs"
    assert second["instrument"] == "soxs"
    assert Path(second["workspace-root-dir"]) == tmp_path


def test_dataframe_factories_return_fresh_tables() -> None:
    first = order_table()
    second = order_table()
    updated = first.assign(order=[99.0, 11.0, 12.0])

    assert first is not second
    assert updated.loc[0, "order"] == 99
    assert first.loc[0, "order"] == 10
    assert second.loc[0, "order"] == 10
    assert list(dispersion_table()["axis"]) == ["x", "y"]
    assert list(qc_table()["qc_name"]) == ["RON"]
    assert list(product_table()["file_name"]) == ["MASTER_BIAS_VIS.fits"]


def test_seeded_frame_and_signal_factories_are_deterministic() -> None:
    firstFrame = synthetic_ccd(seed=11)
    secondFrame = synthetic_ccd(seed=11)
    firstSignal = synthetic_signal(seed=11)
    secondSignal = synthetic_signal(seed=11)

    assert firstFrame is not secondFrame
    assert_array_equal(firstFrame.data, secondFrame.data)
    assert_array_equal(firstSignal, secondSignal)
    assert not firstSignal.flags.writeable
