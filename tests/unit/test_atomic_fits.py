"""Tests for the atomic FITS writer that stops an interrupted run leaving a partial product."""

import os
import warnings
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.io.fits.verify import VerifyWarning

from soxspipe.commonutils.atomic_fits import is_complete_fits, write_hdulist_atomically

pytestmark = pytest.mark.unit


@pytest.fixture
def empty_dir(tmp_path: Path) -> Path:
    """An empty directory, since `tmp_path` also holds the suite's home, cache and config folders."""
    directory = tmp_path / "out"
    directory.mkdir()
    return directory


def _names(directory: Path) -> list[str]:
    return sorted(path.name for path in directory.iterdir())


def _hdulist() -> fits.HDUList:
    return fits.HDUList([fits.PrimaryHDU(np.ones((64, 64), dtype="float32"))])


class _FailingHDUList(fits.HDUList):
    """An HDUList whose write stops half way, as a killed or full-disk write does."""

    def writeto(self, name, *args, **kwargs):
        Path(name).write_bytes(b"SIMPLE  =")
        raise OSError("disk full")


def test_writes_the_hdulist_to_the_final_path(tmp_path: Path) -> None:
    # ARRANGE
    finalPath = tmp_path / "product.fits"

    # ACT
    write_hdulist_atomically(_hdulist(), str(finalPath), checksum=True)

    # ASSERT
    with fits.open(finalPath, checksum=True) as written:
        assert written[0].data.shape == (64, 64)


def test_leaves_no_temporary_file_after_a_successful_write(empty_dir: Path) -> None:
    # ACT
    write_hdulist_atomically(_hdulist(), str(empty_dir / "product.fits"))

    # ASSERT
    assert _names(empty_dir) == ["product.fits"]


def test_replaces_an_existing_product_when_overwrite_is_set(tmp_path: Path) -> None:
    # ARRANGE
    finalPath = tmp_path / "product.fits"
    finalPath.write_bytes(b"old")

    # ACT
    write_hdulist_atomically(_hdulist(), str(finalPath), overwrite=True)

    # ASSERT
    assert is_complete_fits(str(finalPath))


def test_refuses_to_replace_an_existing_product_when_overwrite_is_not_set(empty_dir: Path) -> None:
    # ARRANGE
    finalPath = empty_dir / "product.fits"
    finalPath.write_bytes(b"old")

    # ACT / ASSERT
    with pytest.raises(OSError):
        write_hdulist_atomically(_hdulist(), str(finalPath), overwrite=False)
    assert finalPath.read_bytes() == b"old"
    assert _names(empty_dir) == ["product.fits"]


def test_removes_the_temporary_file_when_the_write_raises(empty_dir: Path) -> None:
    # ACT
    with pytest.raises(OSError, match="disk full"):
        write_hdulist_atomically(_FailingHDUList([fits.PrimaryHDU()]), str(empty_dir / "product.fits"))

    # ASSERT
    assert _names(empty_dir) == []


def test_leaves_an_existing_product_untouched_when_the_write_raises(empty_dir: Path) -> None:
    # ARRANGE
    finalPath = empty_dir / "product.fits"
    finalPath.write_bytes(b"previous product")

    # ACT
    with pytest.raises(OSError):
        write_hdulist_atomically(_FailingHDUList([fits.PrimaryHDU()]), str(finalPath), overwrite=True)

    # ASSERT
    assert finalPath.read_bytes() == b"previous product"
    assert _names(empty_dir) == ["product.fits"]


def test_the_written_file_is_readable_by_other_users_like_an_ordinary_write(tmp_path: Path) -> None:
    # ARRANGE
    plainPath = tmp_path / "plain.fits"
    _hdulist().writeto(plainPath)

    # ACT
    write_hdulist_atomically(_hdulist(), str(tmp_path / "atomic.fits"))

    # ASSERT
    assert (tmp_path / "atomic.fits").stat().st_mode == plainPath.stat().st_mode


def test_is_complete_fits_accepts_a_whole_file(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "whole.fits"
    _hdulist().writeto(path)

    # ACT / ASSERT
    assert is_complete_fits(str(path))


@pytest.mark.parametrize("keptFraction", [0.5, 0.99, 0.0])
def test_is_complete_fits_rejects_a_truncated_file(tmp_path: Path, keptFraction: float) -> None:
    # ARRANGE
    path = tmp_path / "cut.fits"
    _hdulist().writeto(path)
    path.write_bytes(path.read_bytes()[: int(path.stat().st_size * keptFraction)])

    # ACT / ASSERT
    assert not is_complete_fits(str(path))


def test_is_complete_fits_rejects_a_file_that_is_not_fits(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "text.fits"
    path.write_text("not a fits file at all " * 200)

    # ACT / ASSERT
    assert not is_complete_fits(str(path))


def test_is_complete_fits_accepts_a_whole_file_whose_header_cards_astropy_warns_about(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE: AN OVER-LONG COMMENT MAKES ASTROPY WARN "COMMENT WILL BE TRUNCATED" ON A PERFECTLY WHOLE FILE
    primary = fits.PrimaryHDU(np.ones((16, 16), dtype="float32"))
    primary.header["HIERARCH ESO QC LONG KEYWORD NAME"] = (1.0, "c" * 80)
    primary.header["LONGCARD"] = ("value", "c" * 90)
    path = tmp_path / "long_cards.fits"
    with pytest.warns(UserWarning):
        fits.HDUList([primary]).writeto(path, checksum=True)
    realOpen = fits.open

    def open_and_warn(*args: object, **kwargs: object) -> fits.HDUList:
        warnings.warn("Card is too long, comment will be truncated.", VerifyWarning, stacklevel=2)
        return realOpen(*args, **kwargs)

    monkeypatch.setattr(fits, "open", open_and_warn)

    # ACT / ASSERT
    assert is_complete_fits(str(path))


def test_is_complete_fits_rejects_a_file_cut_inside_a_later_extension(tmp_path: Path) -> None:
    # ARRANGE
    path = tmp_path / "table.fits"
    table = fits.BinTableHDU.from_columns([fits.Column(name="a", format="D", array=np.zeros(1000))])
    fits.HDUList([fits.PrimaryHDU(np.ones((16, 16), dtype="float32")), table]).writeto(path)
    path.write_bytes(path.read_bytes()[:-1000])

    # ACT / ASSERT
    assert not is_complete_fits(str(path))


def test_is_complete_fits_returns_false_for_any_error_a_corrupt_file_raises(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # ARRANGE
    path = tmp_path / "corrupt.fits"
    _hdulist().writeto(path)

    def fail_to_decode(*args: object, **kwargs: object) -> None:
        raise EOFError("compressed data ended early")

    monkeypatch.setattr(fits, "open", fail_to_decode)

    # ACT / ASSERT
    assert not is_complete_fits(str(path))


def test_writing_through_a_symlink_keeps_the_link_and_replaces_its_target(empty_dir: Path) -> None:
    # ARRANGE
    target = empty_dir / "target.fits"
    target.write_bytes(b"old")
    link = empty_dir / "link.fits"
    link.symlink_to(target)

    # ACT
    write_hdulist_atomically(_hdulist(), str(link), overwrite=True)

    # ASSERT
    assert link.is_symlink()
    assert is_complete_fits(str(target))
    assert _names(empty_dir) == ["link.fits", "target.fits"]


def test_removes_the_temporary_file_when_the_rename_fails(empty_dir: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # ARRANGE
    def fail_to_rename(*args: object, **kwargs: object) -> None:
        raise PermissionError("not allowed")

    monkeypatch.setattr(os, "replace", fail_to_rename)

    # ACT
    with pytest.raises(PermissionError, match="not allowed"):
        write_hdulist_atomically(_hdulist(), str(empty_dir / "product.fits"))

    # ASSERT
    assert _names(empty_dir) == []
