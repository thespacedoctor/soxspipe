"""Unit tests for the shared FITS frame-name rules."""

from __future__ import annotations

import pytest

from soxspipe.commonutils.fits_frame_names import (
    fits_frame_stem,
    is_fits_frame,
    superseded_frame_names,
)

pytestmark = pytest.mark.unit


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("SOXS.2026-01-01T00:00:00.000.fits", True),
        ("SOXS.2026-01-01T00:00:00.000.fits.Z", True),
        ("./raw/2026-01-01/frame.fits.Z", True),
        ("FRAME.FITS", True),
        ("frame.Fits.Z", True),
        ("frame.fits.z", False),
        ("frame.fits.gz", False),
        ("frame.fz", False),
        ("frame.Z", False),
        ("frame.fits.Z.bak", False),
        ("frame.fitsZ", False),
        ("notes.txt", False),
    ],
)
def test_is_fits_frame_accepts_plain_and_lzw_compressed_fits_only(name: str, expected: bool) -> None:
    assert is_fits_frame(name) is expected


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("SOXS.2026-01-01T00:00:00.000.fits", "SOXS.2026-01-01T00:00:00.000"),
        ("SOXS.2026-01-01T00:00:00.000.fits.Z", "SOXS.2026-01-01T00:00:00.000"),
        ("/data/raw/2026-01-01/frame.fits.Z", "frame"),
        ("FRAME.FITS", "FRAME"),
    ],
)
def test_fits_frame_stem_drops_the_fits_and_lzw_suffixes(name: str, expected: str) -> None:
    assert fits_frame_stem(name) == expected


@pytest.mark.parametrize("name", ["frame.fits.gz", "frame.txt", "frame.Z"])
def test_fits_frame_stem_rejects_a_name_that_is_not_a_fits_frame(name: str) -> None:
    with pytest.raises(ValueError, match=name):
        fits_frame_stem(name)


def test_an_uncompressed_frame_is_superseded_by_its_compressed_twin() -> None:
    names = ["b.fits.Z", "a.fits", "b.fits", "c.fits.Z"]

    assert superseded_frame_names(names) == {"b.fits"}


def test_a_twin_known_only_from_elsewhere_still_supersedes() -> None:
    assert superseded_frame_names(["b.fits", "a.fits"], present=["b.fits.Z"]) == {"b.fits"}


def test_lone_frames_are_never_superseded() -> None:
    assert superseded_frame_names(["a.fits", "b.fits.Z", "notes.txt"]) == set()
