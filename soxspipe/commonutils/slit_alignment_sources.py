"""Where a slit-shifted pixel takes its value from: Newton steps through the 2D map, and safe bilinear reads.

The slit alignment resamples a frame so that a pixel of the target takes its value from a source with the same
wavelength and a slit coordinate shifted by a given amount. This module finds those sources from the map images and
reads an image at them.
"""

import numpy as np
from scipy.ndimage import map_coordinates

# HOW MANY TIMES A SLIT-SHIFT SOURCE IS RE-STEPPED FROM ITSELF, FOR A MAP WHOSE PLATE SCALE CHANGES ACROSS THE SHIFT
SOURCE_REFINEMENTS = 2
# THE WAVELENGTH RESIDUAL IS TRUSTED IN A SOURCE STEP ONLY WHERE THE SINE OF THE ANGLE BETWEEN THE WAVELENGTH AND SLIT
# GRADIENTS IS AT LEAST THIS. THE REAL MAPS HAVE A SINE NEAR 1
MIN_JACOBIAN_SINE = 1e-3


def _gradient_determinant(wavelengthRow, wavelengthCol, slitRow, slitCol):
    """*the determinant of the wavelength and slit gradients: the signed (wavelength, slit) area of a pixel*"""
    return wavelengthRow * slitCol - wavelengthCol * slitRow


def _interpolate(image, sources, mode="nearest"):
    """*the bilinear interpolation of a float32 image at (row, column) sources*

    Beyond the image edge, ``mode`` ``"nearest"`` takes the edge value and ``"constant"`` takes zero.
    """
    return _map_coordinates(image, sources, mode=mode)


def _map_coordinates(image, sources, mode):
    """*bilinear ``map_coordinates`` of an image, with every source made safe to read first*"""
    return map_coordinates(image, _safe_sources(sources, image.shape), order=1, mode=mode, cval=0.0)


def _safe_sources(sources, shape):
    """*the (row, column) sources, made safe for ``map_coordinates``, as new arrays*

    scipy's ``map_coordinates`` converts each coordinate to a 64-bit integer index without a range check in
    ``mode="nearest"``. A coordinate of ``+inf``, or a finite one above about ``9.2e18``, overflows that index, and
    the read lands just before the image: a crash when the page there is unmapped, and a garbage value otherwise.

    A non-finite coordinate becomes NaN, which gives NaN for ``"nearest"`` and ``cval`` for ``"constant"``. A finite
    coordinate is clipped to one pixel beyond the image, which keeps the edge value for ``"nearest"`` and ``cval``
    for ``"constant"``, since bilinear weight at -1 or ``n`` lies wholly outside the image.

    **Key Arguments:**

    - ``sources`` -- the (rows, columns) coordinate arrays
    - ``shape`` -- the (height, width) of the image they index
    """
    safe = []
    for coordinates, size in zip(sources, shape, strict=True):
        coordinates = np.asarray(coordinates, dtype=float)
        safe.append(np.where(np.isfinite(coordinates), np.clip(coordinates, -1, size), np.nan))
    return safe


def _order_sources(mapImages, order, shift):
    """*where each pixel of an order takes its value from, for a slit shift in arcsec*

    The source is the position that has the target's wavelength and a slit coordinate ``shift`` arcsec higher. It is
    found with Newton steps on the two residuals, wavelength and slit, using the 2x2 inverse of the local gradients
    (see ``_newton_step``). The first step starts from the target, where the wavelength residual is zero. The plate
    scale can change across a nod throw, by about 8% at the red end of the SOXS NIR order 10, and the map is tilted
    and curved, so one step overshoots in the slit coordinate and drifts in wavelength, by 0.013 nm (about a tenth of
    a pixel) at the red end of the SOXS NIR order 10. The step is then repeated from the source, along the gradients
    there, for the wavelength and slit residuals still left, ``SOURCE_REFINEMENTS`` times.

    **Return:**

    - ``targets`` -- the (rows, columns) index arrays of the order's pixels
    - ``sources`` -- the (rows, columns) source coordinates of those pixels
    """
    targets = np.nonzero(mapImages.order == order)
    targetWavelength = mapImages.wavelength[targets].astype(np.float64)
    targetSlit = mapImages.slit[targets].astype(np.float64) + shift
    gradients = [image[targets] for image in _gradient_images(mapImages)]
    rows, columns = _newton_step(targets[0], targets[1], gradients, 0.0, shift)
    for _ in range(SOURCE_REFINEMENTS):
        sources = [rows, columns]
        wavelengthMissing = targetWavelength - _interpolate(mapImages.wavelength, sources)
        slitMissing = targetSlit - _interpolate(mapImages.slit, sources)
        gradients = [_interpolate(image, sources) for image in _gradient_images(mapImages)]
        newRows, newColumns = _newton_step(rows, columns, gradients, wavelengthMissing, slitMissing)
        # A STEP THAT IS NOT FINITE, AS THE SOURCE STENCIL TOUCHES AN OFF-ORDER PIXEL (NaN WAVELENGTH), KEEPS THE
        # PREVIOUS SOURCE. A PIXEL WITH A ZERO DETERMINANT AT THE TARGET IS NON-FINITE FROM THE FIRST STEP AND STAYS
        # NON-FINITE. THE CALLER FLAGS BOTH BAD
        isRefined = np.isfinite(newRows) & np.isfinite(newColumns)
        rows, columns = np.where(isRefined, newRows, rows), np.where(isRefined, newColumns, columns)
    return targets, [rows, columns]


def _gradient_images(mapImages):
    """*the wavelength row, wavelength column, slit row and slit column gradient images, in that order*"""
    return (mapImages.wavelengthRow, mapImages.wavelengthCol, mapImages.slitRow, mapImages.slitCol)


def _newton_step(rows, columns, gradients, wavelengthResidual, slitResidual):
    """*the pixel coordinates a Newton step away that closes a wavelength and a slit residual*

    The step solves ``J @ (dRow, dColumn) = (wavelengthResidual, slitResidual)``, where ``J`` is the 2x2 matrix of the
    wavelength and slit gradients, by its inverse: ``dRow = (slitColumn * wavelengthResidual - wavelengthColumn *
    slitResidual) / det`` and ``dColumn = (wavelengthRow * slitResidual - slitRow * wavelengthResidual) / det``.

    The wavelength residual is trusted only where ``J`` is well conditioned, where the sine of the angle between
    the wavelength and slit gradients is at least ``MIN_JACOBIAN_SINE``. Elsewhere the step is the slit-only step,
    the wavelength residual taken as zero, so a noisy wavelength residual is not amplified by a near-zero determinant. A
    zero determinant, or non-finite gradients, still give a non-finite coordinate, never a finite wrong one. The
    caller flags that pixel bad, and ``_safe_sources`` keeps the non-finite value from being read as a memory offset.

    **Key Arguments:**

    - ``rows``, ``columns`` -- the coordinates to step from
    - ``gradients`` -- the wavelength row, wavelength column, slit row and slit column gradients at those coordinates
    - ``wavelengthResidual`` -- the wavelength still missing in nm, a scalar or one value per coordinate
    - ``slitResidual`` -- the slit coordinate still missing in arcsec, a scalar or one value per coordinate

    **Return:**

    - ``rows``, ``columns`` -- the new coordinates
    """
    wavelengthRow, wavelengthCol, slitRow, slitCol = gradients
    with np.errstate(divide="ignore", invalid="ignore"):
        determinant = _gradient_determinant(*gradients)
        sine = np.abs(determinant) / (np.hypot(wavelengthRow, wavelengthCol) * np.hypot(slitRow, slitCol))
        isWellConditioned = sine >= MIN_JACOBIAN_SINE
        wavelengthResidual = np.where(isWellConditioned, wavelengthResidual, 0.0)
        newRows = rows + (slitCol * wavelengthResidual - wavelengthCol * slitResidual) / determinant
        newColumns = columns + (wavelengthRow * slitResidual - slitRow * wavelengthResidual) / determinant
        # AN INFINITE GRADIENT WOULD GIVE AN INFINITE DETERMINANT AND A ZERO STEP IN ONE AXIS, A FINITE WRONG SOURCE
        isFinite = np.isfinite(determinant)
        newRows, newColumns = np.where(isFinite, newRows, np.nan), np.where(isFinite, newColumns, np.nan)
    return newRows, newColumns
