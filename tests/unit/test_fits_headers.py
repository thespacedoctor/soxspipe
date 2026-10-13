"""Unit tests for the fast header reader and the header-reading image collection."""

from __future__ import annotations

import io
import logging
import os
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
import uncompresspy
from astropy import units as u
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection

from soxspipe.commonutils import fits_headers
from soxspipe.commonutils.fits_headers import END_CARD, FITS_CARD_SIZE, HeaderImageFileCollection, read_ccd, read_header
from tests.factories import lzw_compress, lzw_compressed_fits, raw_fits

pytestmark = pytest.mark.unit

FITS_BLOCK_SIZE = 2880
CARDS_PER_BLOCK = 36
LARGE_SHAPE = (256, 256)
BLOCK_FILLING_EXTRA_CARDS = 30
HEADER_CUTTING_BYTES = 100
SYSTEM_COMPRESS = "/usr/bin/compress"
PIXEL_READ_KWARGS = {
    "hdu": 0,
    "unit": u.adu,
    "hdu_uncertainty": "ERRS",
    "hdu_mask": "QUAL",
    "hdu_flags": "FLAGS",
    "key_uncertainty_type": "UTYPE",
}
requiresGzip = pytest.mark.skipif(shutil.which("gzip") is None, reason="the gzip binary is not installed")


def _lzw_frame(tmp_path: Path, name: str = "frame.fits", **rawKwargs: object) -> Path:
    """Write a plain frame and its LZW-compressed twin; return the compressed path."""
    plainPath = raw_fits(tmp_path / name, **rawKwargs)
    return lzw_compressed_fits(plainPath, tmp_path / f"{name}.Z")


def _long_header_frame(tmp_path: Path, cardCount: int) -> Path:
    """Write an LZW frame whose primary header carries ``cardCount`` extra cards."""
    header = fits.Header()
    for index in range(cardCount):
        header[f"KEY{index:04d}"] = index
    plainPath = tmp_path / "long.fits"
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.int16), header=header).writeto(plainPath)
    return lzw_compressed_fits(plainPath, tmp_path / "long.fits.Z")


def _count_blocks_read(monkeypatch: pytest.MonkeyPatch) -> list[int]:
    """Record the size of every chunk read from an LZW stream; return the live list."""
    chunkSizes: list[int] = []
    originalRead = uncompresspy.LZWFile.read

    def _counting_read(self: uncompresspy.LZWFile, size: int = -1) -> bytes:
        chunk = originalRead(self, size)
        chunkSizes.append(len(chunk))
        return chunk

    monkeypatch.setattr(uncompresspy.LZWFile, "read", _counting_read)
    return chunkSizes


def _spy_on_getheader(monkeypatch: pytest.MonkeyPatch) -> list[object]:
    """Record the path of every ``fits.getheader`` call, still running the real function."""
    calledWith: list[object] = []
    originalGetheader = fits.getheader

    def _recording_getheader(path: object, *args: object, **kwargs: object) -> fits.Header:
        calledWith.append(path)
        return originalGetheader(path, *args, **kwargs)

    monkeypatch.setattr(fits, "getheader", _recording_getheader)
    return calledWith


def _forbid_lzw_seek(monkeypatch: pytest.MonkeyPatch) -> None:
    """Make any seek through an LZW stream fail, as a seek decodes the whole file."""

    def _forbidden_seek(*_args: object, **_kwargs: object) -> int:
        raise AssertionError("the fast path must not seek through the data unit")

    monkeypatch.setattr(uncompresspy.LZWFile, "seek", _forbidden_seek)


def test_read_header_matches_getheader_for_an_lzw_frame(tmp_path: Path) -> None:
    framePath = _lzw_frame(tmp_path)

    header = read_header(framePath)

    assert header.tostring() == fits.getheader(framePath, 0).tostring()


def test_read_header_returns_a_header_spanning_several_blocks(tmp_path: Path) -> None:
    framePath = _long_header_frame(tmp_path, cardCount=CARDS_PER_BLOCK * 2)

    header = read_header(framePath)

    assert len(header) > CARDS_PER_BLOCK * 2
    assert header.tostring() == fits.getheader(framePath, 0).tostring()


@pytest.mark.parametrize(
    ("cardCount", "expectedBlocks"),
    [(BLOCK_FILLING_EXTRA_CARDS, 1), (BLOCK_FILLING_EXTRA_CARDS + 1, 2)],
)
def test_read_header_stops_at_the_block_holding_the_end_card(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, cardCount: int, expectedBlocks: int
) -> None:
    # 5 MANDATORY CARDS PLUS 30 EXTRA PUTS END AS THE LAST CARD OF BLOCK 1; ONE MORE PUTS IT FIRST IN BLOCK 2.
    # THE CARDS ALONE CANNOT TELL THE TWO READS APART, SO COUNT THE BLOCKS READ
    framePath = _long_header_frame(tmp_path, cardCount=cardCount)
    chunkSizes = _count_blocks_read(monkeypatch)
    getheaderCalls = _spy_on_getheader(monkeypatch)

    header = read_header(framePath)

    assert chunkSizes == [FITS_BLOCK_SIZE] * expectedBlocks
    assert getheaderCalls == []
    assert header.tostring() == fits.getheader(framePath, 0).tostring()
    assert f"KEY{cardCount - 1:04d}" in header


def test_read_header_is_not_stopped_by_end_text_in_a_value_or_comment(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    header = fits.Header()
    header["ENDKEY"] = "END"
    header.add_comment("END")
    header.add_history("END")
    for index in range(CARDS_PER_BLOCK):
        header[f"KEY{index:04d}"] = index
    plainPath = tmp_path / "endtext.fits"
    fits.PrimaryHDU(data=np.zeros((4, 4), dtype=np.int16), header=header).writeto(plainPath)
    framePath = lzw_compressed_fits(plainPath, tmp_path / "endtext.fits.Z")
    chunkSizes = _count_blocks_read(monkeypatch)

    result = read_header(framePath)

    assert chunkSizes == [FITS_BLOCK_SIZE] * 2
    assert result["ENDKEY"] == "END"
    assert f"KEY{CARDS_PER_BLOCK - 1:04d}" in result
    assert result.tostring() == fits.getheader(framePath, 0).tostring()


def test_read_header_reads_an_lzw_frame_given_as_a_bytes_path(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    framePath = _lzw_frame(tmp_path, shape=LARGE_SHAPE)
    expected = fits.getheader(framePath, 0).tostring()
    _forbid_lzw_seek(monkeypatch)
    getheaderCalls = _spy_on_getheader(monkeypatch)

    header = read_header(os.fsencode(framePath))

    assert getheaderCalls == []
    assert header.tostring() == expected


def test_read_header_falls_back_when_the_end_card_is_malformed(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # THE MATCH IS STRICT: "END" FOLLOWED BY ANYTHING BUT SPACES IS NOT AN END CARD, SO THE FAST PATH GIVES UP
    plainPath = raw_fits(tmp_path / "malformed.fits")
    rawBytes = plainPath.read_bytes()
    endAt = rawBytes.index(END_CARD)
    malformedCard = (b"END" + b" " * 7 + b"JUNK").ljust(FITS_CARD_SIZE)
    malformedPath = tmp_path / "malformed.fits.Z"
    malformedPath.write_bytes(lzw_compress(rawBytes[:endAt] + malformedCard + rawBytes[endAt + FITS_CARD_SIZE :]))
    getheaderCalls = _spy_on_getheader(monkeypatch)

    header = read_header(malformedPath)

    assert getheaderCalls == [malformedPath]
    assert header.tostring() == fits.getheader(malformedPath, 0).tostring()


def test_read_header_delegates_a_plain_frame_to_getheader(tmp_path: Path) -> None:
    framePath = raw_fits(tmp_path / "plain.fits")

    header = read_header(framePath)

    assert header.tostring() == fits.getheader(framePath, 0).tostring()


def test_read_header_returns_the_requested_extension_of_an_lzw_frame(tmp_path: Path) -> None:
    extensionHeader = fits.Header({"EXTNAME": "SECOND", "MARKER": 42})
    hdus = fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data=np.zeros((4, 4)), header=extensionHeader)])
    plainPath = tmp_path / "multi.fits"
    hdus.writeto(plainPath)
    framePath = lzw_compressed_fits(plainPath, tmp_path / "multi.fits.Z")

    header = read_header(framePath, 1)

    assert header["MARKER"] == 42
    assert header.tostring() == fits.getheader(framePath, 1).tostring()


def test_read_header_does_not_decode_the_data_unit(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    framePath = _lzw_frame(tmp_path, shape=LARGE_SHAPE)
    decodedSize = fits.getdata(framePath).nbytes
    bytesRead = []
    originalRead = uncompresspy.LZWFile.read

    def _counting_read(self: uncompresspy.LZWFile, size: int = -1) -> bytes:
        chunk = originalRead(self, size)
        bytesRead.append(len(chunk))
        return chunk

    monkeypatch.setattr(uncompresspy.LZWFile, "read", _counting_read)
    _forbid_lzw_seek(monkeypatch)

    header = read_header(framePath)

    assert header["NAXIS1"] == LARGE_SHAPE[1]
    assert sum(bytesRead) < decodedSize / 4


def test_read_header_raises_the_getheader_error_for_a_corrupt_lzw_frame(tmp_path: Path) -> None:
    framePath = _lzw_frame(tmp_path)
    corruptPath = tmp_path / "corrupt.fits.Z"
    corruptPath.write_bytes(framePath.read_bytes()[:3])
    with pytest.raises(Exception) as expected:  # noqa: PT011
        fits.getheader(corruptPath, 0)

    with pytest.raises(expected.type):
        read_header(corruptPath)


def test_read_header_falls_back_when_the_stream_ends_before_the_end_card(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # THE CUT FALLS INSIDE THE HEADER BLOCK, SO THE STREAM ENDS BEFORE ANY END CARD
    framePath = _lzw_frame(tmp_path)
    truncatedPath = tmp_path / "truncated.fits.Z"
    truncatedPath.write_bytes(framePath.read_bytes()[:HEADER_CUTTING_BYTES])
    with pytest.raises(OSError, match="Empty or corrupt FITS file|Header size is not multiple") as expected:
        fits.getheader(truncatedPath, 0)
    getheaderCalls = _spy_on_getheader(monkeypatch)

    with pytest.raises(type(expected.value)):
        read_header(truncatedPath)

    # ONLY THE FALLBACK CAN HAVE RAISED THIS, SO THE FAST PATH FAILED FIRST
    assert getheaderCalls == [truncatedPath]


def test_read_header_falls_back_and_logs_when_the_fast_path_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture
) -> None:
    framePath = _lzw_frame(tmp_path)
    expected = fits.getheader(framePath, 0).tostring()

    def _broken_fast_read(*_args: object, **_kwargs: object) -> fits.Header:
        raise RuntimeError("decoder failure")

    monkeypatch.setattr(fits_headers, "_read_lzw_primary_header", _broken_fast_read)

    with caplog.at_level(logging.DEBUG, logger="soxspipe.commonutils.fits_headers"):
        header = read_header(framePath)

    assert header.tostring() == expected
    assert "decoder failure" in caplog.text


def test_read_header_does_not_swallow_the_fallback_error(tmp_path: Path) -> None:
    missingPath = tmp_path / "missing.fits.Z"

    with pytest.raises(FileNotFoundError):
        read_header(missingPath)


def test_collection_summary_matches_ccdproc_for_plain_and_lzw_frames(tmp_path: Path) -> None:
    plainPath = raw_fits(tmp_path / "a_plain.fits", seed=1, headerOverrides={"EXTRA": 5})
    lzwPath = _lzw_frame(tmp_path, "b_lzw.fits", seed=2)
    secondLzwPath = _long_header_frame(tmp_path, cardCount=3)
    fileNames = [str(plainPath), str(lzwPath), str(secondLzwPath)]

    ours = HeaderImageFileCollection(location=str(tmp_path), filenames=[Path(f).name for f in fileNames])
    upstream = ImageFileCollection(location=str(tmp_path), filenames=[Path(f).name for f in fileNames])

    assert ours.summary.colnames == upstream.summary.colnames
    assert ours.summary.pformat(max_lines=-1, max_width=-1) == upstream.summary.pformat(max_lines=-1, max_width=-1)


def test_collection_lowercases_keywords_and_joins_comment_cards(tmp_path: Path) -> None:
    header = fits.Header({"DPR_TYPE": "BIAS"})
    header.add_comment("first")
    header.add_comment("second")
    plainPath = tmp_path / "one.fits"
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(plainPath)
    framePath = lzw_compressed_fits(plainPath, tmp_path / "one.fits.Z")

    collection = HeaderImageFileCollection(location=str(tmp_path), filenames=[framePath.name], keywords=["*"])

    assert list(collection.summary["dpr_type"]) == ["BIAS"]
    assert list(collection.summary["comment"]) == ["first,second"]


def test_collection_keeps_the_first_value_and_warns_on_a_duplicate_keyword(tmp_path: Path) -> None:
    header = fits.Header({"DPR_TYPE": "BIAS"})
    header.add_blank("free text")
    header.append(fits.Card("DUPKEY", "first"), useblanks=False, bottom=True)
    header.append(fits.Card("DUPKEY", "second"), useblanks=False, bottom=True)
    plainPath = tmp_path / "one.fits"
    fits.PrimaryHDU(data=np.zeros((2, 2)), header=header).writeto(plainPath)
    framePath = lzw_compressed_fits(plainPath, tmp_path / "one.fits.Z")

    with pytest.warns(UserWarning, match="multiple entries for .dupkey."):
        collection = HeaderImageFileCollection(location=str(tmp_path), filenames=[framePath.name], keywords=["*"])

    assert list(collection.summary["dupkey"]) == ["first"]


def test_read_header_delegates_an_open_file_object_to_getheader(tmp_path: Path) -> None:
    framePath = raw_fits(tmp_path / "plain.fits")

    with framePath.open("rb") as handle:
        header = read_header(handle)

    assert header.tostring() == fits.getheader(framePath, 0).tostring()


def _spy_on_decompress(monkeypatch: pytest.MonkeyPatch) -> list[object]:
    """Record the path of every ``_decompress_lzw`` call, still running the real function."""
    calledWith: list[object] = []
    realDecompress = fits_headers._decompress_lzw

    def _recording(filePath: object) -> bytes:
        calledWith.append(filePath)
        return realDecompress(filePath)

    monkeypatch.setattr(fits_headers, "_decompress_lzw", _recording)
    return calledWith


def _spy_on_ccddata_read(monkeypatch: pytest.MonkeyPatch, result: object = "sentinel") -> list[tuple[object, dict]]:
    """Replace ``CCDData.read`` with a recorder returning ``result``; return the live call list."""
    calls: list[tuple[object, dict]] = []

    def _recording(source: object, **kwargs: object) -> object:
        calls.append((source, kwargs))
        return result

    monkeypatch.setattr(fits_headers.CCDData, "read", _recording)
    return calls


@requiresGzip
def test_read_ccd_returns_the_same_frame_as_ccddata_read_for_an_lzw_frame(tmp_path: Path) -> None:
    framePath = _lzw_frame(tmp_path)

    frame = read_ccd(framePath, **PIXEL_READ_KWARGS)

    expected = CCDData.read(framePath, **PIXEL_READ_KWARGS)
    assert np.array_equal(frame.data, expected.data)
    assert frame.data.dtype == expected.data.dtype
    assert frame.header.tostring() == expected.header.tostring()
    assert frame.unit == expected.unit


@requiresGzip
def test_read_ccd_decompresses_an_lzw_frame_with_the_gzip_subprocess(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    framePath = _lzw_frame(tmp_path)
    decompressed = _spy_on_decompress(monkeypatch)
    pixelReads = _spy_on_ccddata_read(monkeypatch)

    read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert decompressed == [os.fspath(framePath)]
    source, kwargs = pixelReads[0]
    assert isinstance(source, io.BytesIO)
    assert kwargs == {"format": "fits", **PIXEL_READ_KWARGS}


@requiresGzip
@pytest.mark.skipif(not os.access(SYSTEM_COMPRESS, os.X_OK), reason="the system compress binary is not installed")
def test_read_ccd_reads_a_frame_compressed_by_the_system_compress_binary(tmp_path: Path) -> None:
    plainPath = raw_fits(tmp_path / "frame.fits")
    expectedData = fits.getdata(plainPath)
    subprocess.run([SYSTEM_COMPRESS, "-f", str(plainPath)], check=True)  # noqa: S603
    framePath = tmp_path / "frame.fits.Z"

    frame = read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert np.array_equal(frame.data, expectedData)


def test_read_ccd_falls_back_to_ccddata_read_when_gzip_is_not_installed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    framePath = _lzw_frame(tmp_path)
    monkeypatch.setattr(fits_headers.shutil, "which", lambda name: None)
    pixelReads = _spy_on_ccddata_read(monkeypatch)

    result = read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert result == "sentinel"
    assert pixelReads == [(framePath, PIXEL_READ_KWARGS)]


@pytest.mark.parametrize(
    "error",
    [
        subprocess.CalledProcessError(1, ["gzip"]),
        FileNotFoundError("gzip"),
        PermissionError("gzip"),
    ],
    ids=["nonzero-exit", "binary-vanished", "os-error"],
)
def test_read_ccd_falls_back_to_ccddata_read_when_the_subprocess_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, caplog: pytest.LogCaptureFixture, error: Exception
) -> None:
    framePath = _lzw_frame(tmp_path)

    def _failing(filePath: object) -> bytes:
        raise error

    monkeypatch.setattr(fits_headers, "_decompress_lzw", _failing)
    pixelReads = _spy_on_ccddata_read(monkeypatch)

    with caplog.at_level(logging.DEBUG, logger=fits_headers.log.name):
        read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert pixelReads == [(framePath, PIXEL_READ_KWARGS)]
    assert any("falling back" in record.getMessage() for record in caplog.records)


@requiresGzip
def test_read_ccd_falls_back_when_gzip_rejects_the_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    framePath = tmp_path / "garbage.fits.Z"
    framePath.write_bytes(b"this is not compressed data")
    pixelReads = _spy_on_ccddata_read(monkeypatch)

    read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert pixelReads == [(framePath, PIXEL_READ_KWARGS)]


def test_read_ccd_starts_no_subprocess_for_a_plain_frame(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    framePath = raw_fits(tmp_path / "plain.fits")

    def _forbidden(*args: object, **kwargs: object) -> None:
        raise AssertionError("a subprocess was started")

    monkeypatch.setattr(fits_headers.subprocess, "run", _forbidden)

    frame = read_ccd(framePath, **PIXEL_READ_KWARGS)

    assert np.array_equal(frame.data, fits.getdata(framePath))


@requiresGzip
def test_read_ccd_propagates_a_parse_error_from_the_decompressed_bytes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    framePath = _lzw_frame(tmp_path)

    def _corrupt_read(source: object, **kwargs: object) -> None:
        raise TypeError("buffer is too small for requested array")

    monkeypatch.setattr(fits_headers.CCDData, "read", _corrupt_read)

    with pytest.raises(TypeError, match="buffer is too small"):
        read_ccd(framePath, **PIXEL_READ_KWARGS)
