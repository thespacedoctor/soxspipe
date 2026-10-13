#!/usr/bin/env python
"""
*Fast FITS header and pixel reading for plain and LZW-compressed frames*

``fits.getheader`` on a ``.fits.Z`` file parses the primary header and then seeks
past the data unit. ``uncompresspy`` has no random access, so that seek decodes
the whole file. For the primary header, ``read_header`` instead streams only the
blocks up to the ``END`` card, which is about 100 times faster. Any other request,
and any failure of the fast path, is handed to ``fits.getheader`` so results and
errors stay the same as before.

``CCDData.read`` on a ``.fits.Z`` file decodes the pixels with the same pure-Python
decoder, which takes about 1.2 s per real NIR frame. ``read_ccd`` instead decompresses
the file with the system ``gzip`` (about 0.06 s), which reads LZW natively, and gives
the bytes to astropy. If ``gzip`` is missing or fails, it falls back to ``CCDData.read``
on the path.

This module imports nothing from soxspipe except ``fits_frame_names``, which
itself imports nothing from soxspipe, so it cannot create an import cycle.

Author
: David Young

Date Created
: October 12, 2026
"""

import io
import logging
import os
import shutil
import subprocess
import warnings
from collections import OrderedDict

import uncompresspy
from astropy.io import fits
from astropy.nddata import CCDData
from ccdproc import ImageFileCollection

from soxspipe.commonutils.fits_frame_names import LZW_SUFFIX, is_fits_frame

FITS_BLOCK_SIZE = 2880
GZIP_TIMEOUT_SECONDS = 120
GZIP_STDERR_LOG_LIMIT = 500
FITS_CARD_SIZE = 80
END_CARD = b"END" + b" " * (FITS_CARD_SIZE - 3)

log = logging.getLogger(__name__)

# SET ONCE THE MISSING-GZIP WARNING HAS BEEN LOGGED, SO A WORKSPACE OF FRAMES WARNS ONLY ONCE
_hasWarnedGzipMissing = False


def _block_holds_end_card(block):
    """*return True when an 80-byte-aligned card in the block is the END card*

    **Key Arguments:**

    - ``block`` -- one 2880-byte FITS header block

    **Return:**

    - ``hasEnd`` -- True when the block contains the END card
    """
    return any(block[start : start + FITS_CARD_SIZE] == END_CARD for start in range(0, len(block), FITS_CARD_SIZE))


def _is_lzw_frame_path(path):
    """*return True when ``path`` is a file-system path naming a ``.fits.Z`` frame*

    **Key Arguments:**

    - ``path`` -- the argument given to ``read_header``; a ``str``, ``bytes`` or path-like value is a path,
      and a file object is not

    **Return:**

    - ``isLzw`` -- True for a ``str``, ``bytes`` or path-like value ending ``.fits.Z``
    """
    if not isinstance(path, (str, bytes, os.PathLike)):
        return False
    pathText = os.fsdecode(path)
    return pathText.endswith(LZW_SUFFIX) and is_fits_frame(pathText)


def _read_lzw_primary_header(filePath):
    """*stream an LZW-compressed frame up to its END card and parse the header*

    **Key Arguments:**

    - ``filePath`` -- the path of a ``.fits.Z`` file

    **Return:**

    - ``header`` -- the primary header, as an ``astropy.io.fits.Header``

    **Raises:**

    - ``ValueError`` -- when the stream ends before an END card

    The END card match is strict: ``END`` followed by 77 spaces at an 80-byte-aligned
    offset. A malformed END card (for example, ``END`` followed by other characters)
    is not matched, so the stream ends first and the caller falls back to ``fits.getheader``.
    """
    headerBytes = bytearray()
    with uncompresspy.LZWFile(filePath, "rb") as stream:
        while True:
            block = stream.read(FITS_BLOCK_SIZE)
            if len(block) != FITS_BLOCK_SIZE:
                raise ValueError("stream ended before the END card")
            headerBytes += block
            if _block_holds_end_card(block):
                break
    # HEADER.FROMSTRING REJECTS A BYTEARRAY (TESTED ON ASTROPY 7.2.2), SO THE COPY TO bytes STAYS
    return fits.Header.fromstring(bytes(headerBytes), sep="")


def read_header(path, ext=0):
    """*return a FITS header, streaming only the header blocks of an LZW frame*

    The fast path applies to the primary header (``ext=0``) of a ``.fits.Z``
    file. Everything else, and any failure of the fast path, goes to
    ``fits.getheader``, whose result or exception is returned unchanged.

    **Key Arguments:**

    - ``path`` -- the path of the FITS file
    - ``ext`` -- the extension to read. Default *0*

    **Return:**

    - ``header`` -- the header, as an ``astropy.io.fits.Header``
    """
    if ext == 0 and _is_lzw_frame_path(path):
        try:
            return _read_lzw_primary_header(os.fspath(path))
        except Exception as e:
            # UNCOMPRESSPY DECODER ERRORS ARE NOT A STABLE SET, SO CATCH BROADLY.
            # THE FALLBACK BELOW RE-RAISES THE AUTHORITATIVE ERROR
            log.debug(f"read_header: fast LZW read of `{path}` failed, falling back to fits.getheader: {e}")
    return fits.getheader(path, ext)


def _decompress_lzw(filePath):
    """*decompress a ``.fits.Z`` file with the system ``gzip``, which reads LZW natively*

    **Key Arguments:**

    - ``filePath`` -- the path of a ``.fits.Z`` file

    **Return:**

    - ``data`` -- the decompressed bytes of the FITS file

    **Raises:**

    - ``FileNotFoundError`` -- when no ``gzip`` binary is on the ``PATH``
    - ``subprocess.CalledProcessError`` -- when ``gzip`` exits with an error
    - ``subprocess.TimeoutExpired`` -- when ``gzip`` runs longer than ``GZIP_TIMEOUT_SECONDS``
    - ``OSError`` -- when ``gzip`` cannot be started
    """
    gzipPath = shutil.which("gzip")
    if gzipPath is None:
        raise FileNotFoundError("gzip is not on the PATH")
    # LIST ARGUMENTS AND NO SHELL: THE PATH IS NEVER PARSED BY A SHELL. `--` ENDS OPTION PARSING,
    # SO A FILE NAME STARTING WITH A DASH IS NOT TAKEN FOR A GZIP OPTION
    result = subprocess.run(  # noqa: S603
        [gzipPath, "-dc", "--", filePath],
        capture_output=True,
        check=True,
        timeout=GZIP_TIMEOUT_SECONDS,
    )
    return result.stdout


def _describe_failure(error):
    """*return a log-ready description of a failed gzip run, with its bounded error output*

    **Key Arguments:**

    - ``error`` -- the exception raised by ``_decompress_lzw``

    **Return:**

    - ``description`` -- the exception text, followed by the first ``GZIP_STDERR_LOG_LIMIT``
      characters of ``gzip``'s error output when there is any
    """
    stderr = getattr(error, "stderr", None)
    if not stderr:
        return str(error)
    if isinstance(stderr, bytes):
        stderr = stderr.decode(errors="replace")
    return f"{error}: {stderr.strip()[:GZIP_STDERR_LOG_LIMIT]}"


def _warn_once_gzip_missing():
    """*log a warning, once per process, that gzip is missing and LZW reads will be slow*"""
    global _hasWarnedGzipMissing
    if _hasWarnedGzipMissing:
        return
    _hasWarnedGzipMissing = True
    log.warning(
        "gzip was not found, so .fits.Z frames are read with the much slower uncompresspy; "
        "install gzip to speed up reading"
    )


def read_ccd(path, **kwargs):
    """*return a CCDData frame, decompressing an LZW frame with the system* ``gzip``

    The pure-Python ``uncompresspy`` decoder takes about 1.2 s per real NIR frame. ``gzip`` decodes the
    same file in about 0.06 s. Only the decompression step falls back to ``CCDData.read`` on the path,
    when ``gzip`` is missing, fails or times out. A missing ``gzip`` is logged once as a warning. An error
    raised by astropy while it parses the decompressed bytes is not caught, so the caller sees the same
    exception as before. The decompressed bytes are read from memory, so ``memmap`` does not apply on the
    fast path.

    **Key Arguments:**

    - ``path`` -- the path of the FITS file
    - ``kwargs`` -- the keyword arguments given to ``CCDData.read``. A ``format`` given here is
      replaced by ``"fits"`` on the fast path

    **Return:**

    - ``frame`` -- the frame, as an ``astropy.nddata.CCDData``
    """
    if not _is_lzw_frame_path(path):
        return CCDData.read(path, **kwargs)
    try:
        data = _decompress_lzw(os.fsdecode(path))
    except FileNotFoundError as e:
        _warn_once_gzip_missing()
        log.debug(f"read_ccd: gzip decompression of `{path}` failed, falling back to CCDData.read: {e}")
        return CCDData.read(path, **kwargs)
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        log.debug(
            f"read_ccd: gzip decompression of `{path}` failed, falling back to CCDData.read: {_describe_failure(e)}"
        )
        return CCDData.read(path, **kwargs)
    return CCDData.read(io.BytesIO(data), **{**kwargs, "format": "fits"})


class HeaderImageFileCollection(ImageFileCollection):
    """*a ccdproc ImageFileCollection that reads headers with* ``read_header``

    ``_dict_from_fits_header`` mirrors ccdproc 2.4.2 (pinned in ``pyproject.toml``),
    with ``fits.getheader`` replaced by ``read_header``. Re-check it against
    upstream whenever the ccdproc pin changes.

    Unlike the subclass in ``set_of_files``, it deliberately has no fallback to
    extension 0 when the requested extension cannot be read. It mirrors plain
    ccdproc, which is what ``data_organiser`` used before.
    """

    def _dict_from_fits_header(self, file_name, input_summary=None, missing_marker=None):
        """*summarise one file's header, lower-casing keywords (as ccdproc 2.4.2 does)*

        **Key Arguments:**

        - ``file_name`` -- the path of the FITS file
        - ``input_summary`` -- the summary built so far. Default *None*
        - ``missing_marker`` -- the value recorded for a keyword the file lacks. Default *None*

        **Return:**

        - ``summary`` -- the summary, with this file's row appended
        """

        def _add_val_to_dict(key, value, tbl_dict, n_previous, missing_marker):
            try:
                tbl_dict[key].append(value)
            except KeyError:
                tbl_dict[key] = [missing_marker] * n_previous
                tbl_dict[key].append(value)

        if input_summary is None:
            summary = OrderedDict()
            n_previous = 0
        else:
            summary = input_summary
            n_previous = len(summary["file"])

        h = read_header(file_name, self.ext)

        # KEEP THE BARE ASSERT TO MIRROR UPSTREAM: AN EXPLICIT RAISE WOULD ALSO FIRE UNDER `python -O`
        assert "file" not in h  # noqa: S101

        # WITH A LOCATION WE CAN RECONSTRUCT THE PATH USING IT; WITHOUT ONE, USE WHATEVER PATH THE USER PASSED IN
        name_for_file_column = os.path.basename(file_name) if self.location else file_name

        # OPEN THE HEADER BEFORE THIS SO THAT THE FILE NAME IS ONLY ADDED IF THE FILE IS VALID FITS
        try:
            summary["file"].append(name_for_file_column)
        except KeyError:
            summary["file"] = [name_for_file_column]

        missing_in_this_file = [k for k in summary if (k not in h and k != "file")]

        multi_entry_keys = {"comment": [], "history": []}

        alreadyencountered = set()
        for k, v in h.items():
            if k == "":
                continue

            k = k.lower()

            if k in ["comment", "history"]:
                # ACCUMULATE THESE SEPARATELY UNTIL THE END TO AVOID ADDING MULTIPLE ENTRIES TO THE SUMMARY
                multi_entry_keys[k].append(str(v))
                continue
            if k in alreadyencountered:
                # THE NORMAL MULTI-ENTRY KEYS ARE ALREADY HANDLED, SO ANY OTHER DUPLICATE IS PROBABLY A MISTAKE
                warnings.warn(
                    f'Header from file "{file_name}" contains multiple entries for '
                    f'"{k}", the pair "{k}={v}" will be ignored.',
                    UserWarning,
                    stacklevel=1,
                )
                continue
            alreadyencountered.add(k)

            _add_val_to_dict(k, v, summary, n_previous, missing_marker)

        for k, v in multi_entry_keys.items():
            if v:
                joined = ",".join(v)
                _add_val_to_dict(k, joined, summary, n_previous, missing_marker)

        for missing in missing_in_this_file:
            summary[missing].append(missing_marker)

        return summary
