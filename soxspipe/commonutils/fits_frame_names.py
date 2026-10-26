#!/usr/bin/env python
"""
*Shared rules for recognising FITS frames and naming files derived from them*

A raw frame arrives either as a plain ``.fits`` file or as an LZW-compressed
``.fits.Z`` file from the ESO archive. astropy reads both directly, so the
pipeline treats them as the same kind of frame and never uncompresses on disk.

This module imports nothing from soxspipe, so any module can import it without
creating an import cycle.

Author
: David Young

Date Created
: October 7, 2026
"""

import os

FITS_SUFFIX = ".fits"
LZW_SUFFIX = ".Z"


def is_fits_frame(path):
    """*return True when the file name is a plain or LZW-compressed FITS frame*

    The ``.fits`` suffix matches in any letter case. The ``.Z`` suffix must be
    an upper-case ``Z``, as written by the Unix ``compress`` command.

    **Key Arguments:**

    - ``path`` -- a file name or path

    **Return:**

    - ``isFrame`` -- True for names ending ``.fits`` or ``.fits.Z``
    """
    name = os.path.basename(path)
    if name.endswith(LZW_SUFFIX):
        name = name[: -len(LZW_SUFFIX)]
    return name.lower().endswith(FITS_SUFFIX)


def fits_frame_stem(path):
    """*return the frame's file name without its ``.fits`` or ``.fits.Z`` suffix*

    Use this wherever a new file name is built from a raw frame name, so no
    output name carries the ``.Z`` of a compressed input.

    **Key Arguments:**

    - ``path`` -- a file name or path that passes ``is_fits_frame``

    **Return:**

    - ``stem`` -- the base name with the FITS suffixes removed

    **Raises:**

    - ``ValueError`` -- when ``path`` is not a FITS frame name
    """
    if not is_fits_frame(path):
        raise ValueError(f"`{path}` is not a .fits or .fits.Z frame")
    name = os.path.basename(path)
    if name.endswith(LZW_SUFFIX):
        name = name[: -len(LZW_SUFFIX)]
    return name[: -len(FITS_SUFFIX)]


def superseded_frame_names(names, present=()):
    """*return the uncompressed frame names that have a compressed twin*

    When ``X.fits`` and ``X.fits.Z`` both exist, the compressed file is the
    one the pipeline keeps. The result is a set, so it does not depend on the
    order in which a directory was listed.

    **Key Arguments:**

    - ``names`` -- the candidate file names or paths
    - ``present`` -- further names known to exist (for example, frames already in the database). Default *()*

    **Return:**

    - ``superseded`` -- the members of ``names`` whose ``.Z`` twin is in ``names`` or ``present``
    """
    available = set(names) | set(present)
    return {
        name
        for name in names
        if not name.endswith(LZW_SUFFIX) and is_fits_frame(name) and name + LZW_SUFFIX in available
    }
