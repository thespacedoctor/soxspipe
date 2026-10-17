"""Out-of-range resampling sources must never make scipy read outside the image (DY-1615)."""

from __future__ import annotations

import ctypes
import dataclasses
import os
import subprocess
import sys
import textwrap
from pathlib import Path

import numpy as np
import pytest
from scipy.ndimage import map_coordinates

from soxspipe.commonutils.slit_alignment import _interpolate, build_two_d_map_images, shift_frame_along_slit
from tests.factories.slit_alignment import synthetic_alignment_map, synthetic_nod_frame

pytestmark = pytest.mark.unit

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PROBE_TIMEOUT_SECONDS = 120
HEIGHT, WIDTH = 8, 16
# A FINITE ROW COORDINATE LARGER THAN THE INT64 RANGE (ABOUT 9.2E18) OVERFLOWS LIKE INFINITY
INT64_OVERFLOW_ROW = 1e19

# THE IMAGE SITS IN AN MMAP WHOSE PAGE BEFORE IT IS UNREADABLE, SO ANY READ JUST BEFORE THE IMAGE IS FATAL
GUARD_PAGE_SCRIPT = textwrap.dedent(
    """
    import ctypes, mmap, sys
    import numpy as np
    from soxspipe.commonutils.slit_alignment import _interpolate

    PAGE = mmap.PAGESIZE
    height, width = 8, 16
    buf = mmap.mmap(-1, PAGE + height * width * 4 + PAGE)
    address = ctypes.addressof(ctypes.c_char.from_buffer(buf))
    libc = ctypes.CDLL(None)
    libc.mprotect.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int]
    assert libc.mprotect(address, PAGE, 0) == 0
    image = np.frombuffer(buf, dtype=np.float32, count=height * width, offset=PAGE).reshape(height, width)
    image[:] = 1.0
    sources = [np.array([float(sys.argv[1])]), np.array([5.0])]
    _interpolate(image, sources, mode=sys.argv[2])
    """
)


def _can_protect_pages() -> bool:
    """Return whether this platform lets the probe make a page unreadable."""
    return os.name == "posix" and hasattr(ctypes.CDLL(None), "mprotect")


@pytest.mark.skipif(not _can_protect_pages(), reason="the guard-page probe needs POSIX mprotect")
@pytest.mark.parametrize("mode", ["nearest", "constant"])
@pytest.mark.parametrize("row", ["inf", str(INT64_OVERFLOW_ROW)])
def test_interpolating_an_overflowing_row_does_not_read_before_the_image(row: str, mode: str) -> None:
    # ARRANGE
    environment = {**os.environ, "PYTHONPATH": str(REPOSITORY_ROOT)}

    # ACT
    completed = subprocess.run(  # noqa: S603
        [sys.executable, "-c", GUARD_PAGE_SCRIPT, row, mode],
        capture_output=True,
        text=True,
        timeout=PROBE_TIMEOUT_SECONDS,
        check=False,
        cwd=REPOSITORY_ROOT,
        env=environment,
    )

    # ASSERT
    assert completed.returncode == 0, f"the probe died with code {completed.returncode}: {completed.stderr}"


@pytest.mark.parametrize("badRow", [np.inf, -np.inf, np.nan])
def test_a_non_finite_source_gives_nan_for_nearest_and_zero_for_constant(badRow: float) -> None:
    # ARRANGE
    image = np.ones((HEIGHT, WIDTH), dtype=np.float32)
    sources = [np.array([badRow, 3.0]), np.array([5.0, 5.0])]

    # ACT
    nearest = _interpolate(image, sources)
    constant = _interpolate(image, sources, mode="constant")

    # ASSERT
    assert np.isnan(nearest[0])
    assert nearest[1] == 1.0
    assert constant[0] == 0.0
    assert constant[1] == 1.0


def test_a_huge_finite_source_gives_the_edge_value_for_nearest_and_zero_for_constant() -> None:
    # ARRANGE: THE SENTINEL ROW SITS JUST BEFORE THE IMAGE, WHERE AN OVERFLOWED READ LANDS
    padded = np.ones((HEIGHT + 1, WIDTH), dtype=np.float32)
    padded[0, :] = 999.0
    image = padded[1:]
    sources = [np.array([INT64_OVERFLOW_ROW, -INT64_OVERFLOW_ROW]), np.array([5.0, 5.0])]

    # ACT
    nearest = _interpolate(image, sources)
    constant = _interpolate(image, sources, mode="constant")

    # ASSERT
    np.testing.assert_array_equal(nearest, [1.0, 1.0])
    np.testing.assert_array_equal(constant, [0.0, 0.0])


@pytest.mark.parametrize("mode", ["nearest", "constant"])
def test_sources_in_range_and_just_outside_it_match_plain_map_coordinates(mode: str) -> None:
    # ARRANGE
    image = np.arange(HEIGHT * WIDTH, dtype=np.float32).reshape(HEIGHT, WIDTH)
    rows = np.array([0.0, 2.5, HEIGHT - 1.0, -0.4, HEIGHT - 0.6, -3.0, HEIGHT + 4.0])
    columns = np.array([1.0, 7.25, WIDTH - 1.0, WIDTH - 0.6, -0.4, 2.0, -9.0])
    sources = [rows, columns]
    expected = map_coordinates(image, [rows.copy(), columns.copy()], order=1, mode=mode, cval=0.0)

    # ACT
    actual = _interpolate(image, sources, mode=mode)

    # ASSERT
    np.testing.assert_array_equal(actual, expected)


def test_interpolating_leaves_the_source_arrays_unchanged() -> None:
    # ARRANGE
    image = np.ones((HEIGHT, WIDTH), dtype=np.float32)
    rows = np.array([np.inf, 1e19, -2.0, 3.0])
    columns = np.array([1.0, np.nan, 99.0, 2.0])
    rowsBefore, columnsBefore = rows.copy(), columns.copy()

    # ACT
    _interpolate(image, [rows, columns])

    # ASSERT
    np.testing.assert_array_equal(rows, rowsBefore)
    np.testing.assert_array_equal(columns, columnsBefore)


def test_a_zero_determinant_pixel_is_flagged_bad_and_no_other_pixel_changes(tmp_path: Path) -> None:
    # ARRANGE: ZERO BOTH SLIT GRADIENTS AT ONE PIXEL SO ITS DETERMINANT IS ZERO AND ITS SOURCE IS INFINITE
    scene = synthetic_alignment_map(tmp_path)
    cleanMap = build_two_d_map_images(scene.wavelength, scene.slit, scene.orderMap)
    badRow, badColumn = np.argwhere(scene.orderMap == 1)[len(np.argwhere(scene.orderMap == 1)) // 2]
    slitRow, slitCol = cleanMap.slitRow.copy(), cleanMap.slitCol.copy()
    slitRow[badRow, badColumn] = slitCol[badRow, badColumn] = 0.0
    degenerateMap = dataclasses.replace(cleanMap, slitRow=slitRow, slitCol=slitCol)
    frame = synthetic_nod_frame(scene, objectSlit=1.0, mjd=60000.0, hasSky=True)
    shifts = {1: 0.6, 2: 0.0, 3: 0.0, 4: 0.0}

    # ACT
    reference = shift_frame_along_slit(frame, shifts, cleanMap)
    shifted = shift_frame_along_slit(frame, shifts, degenerateMap)

    # ASSERT
    assert shifted.mask[badRow, badColumn]
    isOther = np.ones(shifted.data.shape, dtype=bool)
    isOther[badRow, badColumn] = False
    np.testing.assert_array_equal(shifted.data[isOther], reference.data[isOther])
    np.testing.assert_array_equal(shifted.mask[isOther], reference.mask[isOther])
