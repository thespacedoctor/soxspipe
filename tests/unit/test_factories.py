"""Contract tests for fresh, deterministic synthetic factories."""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pytest
import uncompresspy
from numpy.testing import assert_array_equal

from tests.factories import (
    dispersion_table,
    lzw_compress,
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


def test_lzw_compress_matches_the_unix_compress_command() -> None:
    # REFERENCE BYTES FROM `printf 'TOBEORNOTTOBEORTOBEORNOT' | compress -c`
    expected = bytes.fromhex("1f9d90549e0829f2448a932754020e2ca890a04184")

    assert lzw_compress(b"TOBEORNOTTOBEORTOBEORNOT") == expected


def test_lzw_compress_matches_the_unix_compress_command_across_code_widths() -> None:
    # 6000 BYTES GROW THE CODES FROM 9 TO 11 BITS BUT STAY UNDER compress(1)'S
    # 10000-BYTE RATIO CHECKPOINT, SO THE REFERENCE HOLDS NO CLEAR CODE.
    # HASH OF `compress -c` ON THE SAME BYTES (macOS, 2091 BYTES OF OUTPUT)
    payload = bytes((i * i * 7 + i * 3) % 251 for i in range(6000))

    archive = lzw_compress(payload)

    assert len(archive) == 2091
    assert hashlib.sha256(archive).hexdigest() == (
        "0eaa1d7e42f4a60e7f291991ea9c32b55dbc6bb2f6377bd374df96efb04ddc2e"
    )


@pytest.mark.parametrize("size", [0, 1, 300, 5000, 400_000])
def test_lzw_compress_round_trips_across_code_width_changes(
    tmp_path: Path, size: int
) -> None:
    # RANDOM BYTES FILL THE DICTIONARY FAST, SO THE LARGER SIZES CROSS EVERY
    # CODE WIDTH FROM 9 TO 16 BITS AND RUN ON PAST A FULL DICTIONARY
    payload = np.random.default_rng(11).integers(0, 256, size, dtype=np.uint8).tobytes()
    archivePath = tmp_path / "payload.Z"
    archivePath.write_bytes(lzw_compress(payload))

    with uncompresspy.open(archivePath) as archive:
        assert archive.read() == payload
