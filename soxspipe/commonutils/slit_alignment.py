"""Slit alignment: shift nodding and offset frames along the slit so that their object traces stack on each other.

Telescope jitter and drift move the object along the slit from frame to frame. This module measures that movement
and removes it, so that a group of frames can be stacked before extraction.

The work has three steps:

1. The object trace of each frame is found on that frame's cycle difference, and its position along the slit is
   read from the 2D dispersion map.
2. Each frame's slit shift is its position minus the position in the first frame in time.
3. Each frame is resampled through the 2D map to apply the shift in slit coordinates.
"""

import copy
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.ndimage import distance_transform_edt, map_coordinates

# A TRACE NEEDS AT LEAST THIS MANY SAMPLES IN AN ORDER TO GIVE THE ORDER A SLIT POSITION
MIN_TRACE_SAMPLES = 10
# THE LARGEST SCATTER (ROBUST STANDARD DEVIATION, IN ARCSEC) OF THE SLIT COORDINATES OF AN ORDER'S TRACE SAMPLES
MAX_TRACE_SCATTER_ARCSEC = 0.5
# THE SCALE FROM A MEDIAN ABSOLUTE DEVIATION TO A STANDARD DEVIATION
MAD_TO_SIGMA = 1.4826
# A TRACE FOUND IN FEWER THAN THIS PERCENTAGE OF THE SAMPLES IS NOT TRUSTED, AS IN detect_continuum.get
MIN_DETECTION_PERCENT = 10
# A MASK OR COVERAGE VALUE THIS CLOSE TO ITS CLEAN VALUE COUNTS AS CLEAN
INTERPOLATION_TOLERANCE = 1e-6
# HOW MANY TIMES A SLIT-SHIFT SOURCE IS RE-STEPPED FROM ITSELF, FOR A MAP WHOSE PLATE SCALE CHANGES ACROSS THE SHIFT
SOURCE_REFINEMENTS = 2
# THE ORDER VALUE OF A SHIFT THAT BELONGS TO A WHOLE FRAME
WHOLE_FRAME = -1

POSITION_COLUMNS = ["order", "position_arcsec", "n_samples", "scatter_arcsec"]
SHIFT_COLUMNS = ["frame", "mjd", "order", "position_arcsec", "shift_arcsec", "is_fallback", "is_reference"]
TRACE_COLUMNS = ["order", "cont_x", "cont_y", "pre-clipped"]


@dataclass(frozen=True)
class two_d_map_images:
    """The SLIT and ORDER images of a 2D map, with the local gradients of the WAVELENGTH and SLIT images.

    The gradient images (``wavelengthRow``, ``wavelengthCol``, ``slitRow``, ``slitCol``) hold the change of the
    wavelength and of the slit coordinate per pixel step along the rows and the columns. They are measured inside
    each order only, and are valid in the pixels of an order. ``orders`` lists the order numbers the map holds.
    """

    slit: np.ndarray
    order: np.ndarray
    wavelengthRow: np.ndarray
    wavelengthCol: np.ndarray
    slitRow: np.ndarray
    slitCol: np.ndarray
    orders: tuple


def build_two_d_map_images(wavelength, slit, orderMap):
    """*the map images with their per-order gradients*

    **Key Arguments:**

    - ``wavelength`` -- the WAVELENGTH image
    - ``slit`` -- the SLIT image, in arcsec
    - ``orderMap`` -- the ORDER image, zero or NaN outside the orders

    **Return:**

    - ``mapImages`` -- a frozen ``two_d_map_images``, holding float32 copies of the slit and order images and the
      gradients of the wavelength and slit images
    """
    images = [np.asarray(image, dtype=np.float32) for image in (wavelength, slit, orderMap)]
    terms = [np.zeros(images[0].shape, dtype=np.float32) for _ in range(4)]
    wavelengthImage, slitImage, orderImage = images
    isOrder = np.isfinite(orderImage) & (orderImage > 0) & (orderImage == np.round(orderImage))
    orders = np.unique(orderImage[isOrder])
    for order in orders:
        onOrder = orderImage == order
        box = _bounding_box(onOrder)
        local = onOrder[box]
        gradients = [
            *_filled_gradients(np.where(local, wavelengthImage[box], np.nan)),
            *_filled_gradients(np.where(local, slitImage[box], np.nan)),
        ]
        for term, gradient in zip(terms, gradients, strict=True):
            term[box][local] = gradient[local]
    return two_d_map_images(slitImage, orderImage, *terms, tuple(int(order) for order in orders))


def _bounding_box(mask):
    """*the pair of slices that bound the True pixels of a mask*"""
    rows = np.flatnonzero(mask.any(axis=1))
    columns = np.flatnonzero(mask.any(axis=0))
    return slice(rows[0], rows[-1] + 1), slice(columns[0], columns[-1] + 1)


def _filled_gradients(image):
    """*the row and column gradients of an image with NaN outside an order, the NaN rim filled from inside*"""
    gradients = []
    for gradient in np.gradient(image):
        isMissing = ~np.isfinite(gradient)
        if isMissing.any() and not isMissing.all():
            nearest = distance_transform_edt(isMissing, return_distances=False, return_indices=True)
            gradient = gradient[tuple(nearest)]
        gradients.append(gradient)
    return gradients


def order_slit_positions(orderPixelTable, mapImages):
    """*the slit position of the object trace in each order of the 2D map*

    A trace sample counts when it was not clipped and all four pixels around it belong to its order. An order has a
    position when it has at least ``MIN_TRACE_SAMPLES`` samples whose slit coordinates scatter by no more than
    ``MAX_TRACE_SCATTER_ARCSEC``. The position is then their median.

    **Key Arguments:**

    - ``orderPixelTable`` -- the trace samples ``sample_trace`` returned, with the columns ``order``, ``cont_x``,
      ``cont_y`` and ``pre-clipped``
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning

    **Return:**

    - ``positions`` -- one row per order of the map, with the columns ``order``, ``position_arcsec`` (NaN when the
      order has no reliable trace), ``n_samples`` and ``scatter_arcsec``
    """
    isUsable = (
        ~orderPixelTable["pre-clipped"].astype(bool)
        & np.isfinite(orderPixelTable["cont_x"])
        & np.isfinite(orderPixelTable["cont_y"])
    )
    samples = orderPixelTable.loc[isUsable]
    records = []
    for order in mapImages.orders:
        ofOrder = samples.loc[samples["order"] == order]
        slitCoordinates = _slit_coordinates(
            mapImages, order, ofOrder["cont_y"].to_numpy(), ofOrder["cont_x"].to_numpy()
        )
        slitCoordinates = slitCoordinates[np.isfinite(slitCoordinates)]
        position, scatter = np.nan, np.nan
        if len(slitCoordinates):
            position = float(np.median(slitCoordinates))
            scatter = float(MAD_TO_SIGMA * np.median(np.abs(slitCoordinates - position)))
        isReliable = len(slitCoordinates) >= MIN_TRACE_SAMPLES and scatter <= MAX_TRACE_SCATTER_ARCSEC
        records.append(
            {
                "order": order,
                "position_arcsec": position if isReliable else np.nan,
                "n_samples": len(slitCoordinates),
                "scatter_arcsec": scatter,
            }
        )
    return pd.DataFrame(records, columns=POSITION_COLUMNS)


def _slit_coordinates(mapImages, order, rows, columns):
    """*the slit coordinate at each (row, column) inside an order, and NaN where a neighbouring pixel is off-order*"""
    height, width = mapImages.order.shape
    rowBelow, columnBelow = np.floor(rows).astype(int), np.floor(columns).astype(int)
    isInside = (rowBelow >= 0) & (columnBelow >= 0) & (rowBelow + 1 < height) & (columnBelow + 1 < width)
    rowBelow, columnBelow = np.clip(rowBelow, 0, height - 2), np.clip(columnBelow, 0, width - 2)
    for rowStep in (0, 1):
        for columnStep in (0, 1):
            isInside &= mapImages.order[rowBelow + rowStep, columnBelow + columnStep] == order
    values = _map_coordinates(mapImages.slit, [rows, columns], mode="nearest")
    return np.where(isInside, values, np.nan)


def slit_shifts(log, positions, perOrder):
    """*the slit shift of every frame and order, relative to the first frame in time*

    **Key Arguments:**

    - ``log`` -- logger
    - ``positions`` -- one row per frame and order, with the columns ``frame``, ``mjd``, ``order`` and
      ``position_arcsec`` (NaN where the order has no reliable trace)
    - ``perOrder`` -- if True, every order of every frame has its own shift. If False, every frame has one shift

    **Return:**

    - ``shifts`` -- a dataframe with the columns ``frame``, ``mjd``, ``order`` (``-1`` for a whole-frame shift),
      ``position_arcsec``, ``shift_arcsec``, ``is_fallback`` and ``is_reference``
    - ``fallbackCount`` -- the number of shifts that were not measured, and so were taken from other orders or frames
    """
    frames = positions[["frame", "mjd"]].drop_duplicates().sort_values("mjd", kind="stable")
    referenceFrame = frames["frame"].iloc[0]
    measure = _per_order_shifts if perOrder else _per_frame_shifts
    shifts = measure(positions, referenceFrame)
    shifts["is_reference"] = shifts["frame"] == referenceFrame
    shifts.loc[shifts["is_reference"], "shift_arcsec"] = 0.0
    shifts["is_fallback"] = shifts["is_reference"] & shifts["position_arcsec"].isna()
    _warn_unmeasured_reference(log, shifts)

    shifts = _fill_fallback_shifts(log, shifts)
    return shifts[SHIFT_COLUMNS], int(shifts["is_fallback"].sum())


def _per_order_shifts(positions, referenceFrame):
    """*the position of every order of every frame minus the position of the same order in the reference frame*"""
    referencePositions = positions.loc[positions["frame"] == referenceFrame].set_index("order")["position_arcsec"]
    shifts = positions.copy()
    shifts["shift_arcsec"] = shifts["position_arcsec"] - shifts["order"].map(referencePositions)
    return shifts


def _per_frame_shifts(positions, referenceFrame):
    """*one row per frame: the median over orders of the position differences from the reference frame*"""
    referencePositions = positions.loc[positions["frame"] == referenceFrame].set_index("order")["position_arcsec"]
    differences = positions.assign(
        shift_arcsec=positions["position_arcsec"] - positions["order"].map(referencePositions)
    )
    grouped = differences.groupby("frame", sort=False)
    shifts = grouped.agg(
        mjd=("mjd", "first"), position_arcsec=("position_arcsec", "median"), shift_arcsec=("shift_arcsec", "median")
    )
    shifts["order"] = WHOLE_FRAME
    return shifts.reset_index()


def _warn_unmeasured_reference(log, shifts):
    """*warn about each order (or the whole frame) of the reference frame that has no reliable trace*

    The other frames are measured against the reference position, so they have no measured shift there either.
    """
    for _, row in shifts.loc[shifts["is_fallback"]].iterrows():
        where = f"order {int(row['order'])} of " if row["order"] != WHOLE_FRAME else ""
        log.info(
            f"No reliable object trace in {where}the reference frame {row['frame']}, so no other frame can be "
            "measured against it there."
        )


def _fill_fallback_shifts(log, shifts):
    """*give every unmeasured shift the median of the measured shifts of the closest group, and warn about it*

    The closest group is the other orders of the same frame, then the same order in the other frames, then every
    measured shift. A group with no measured shift at all is left unshifted.
    """
    measured = shifts.loc[shifts["shift_arcsec"].notna() & ~shifts["is_fallback"]]
    filled = shifts.copy()
    for index, row in shifts.loc[shifts["shift_arcsec"].isna()].iterrows():
        groups = (
            ("its other orders", measured.loc[measured["frame"] == row["frame"], "shift_arcsec"]),
            ("the same order in the other frames", measured.loc[measured["order"] == row["order"], "shift_arcsec"]),
            ("every measured shift", measured["shift_arcsec"]),
        )
        source, values = next(((name, group) for name, group in groups if len(group)), ("nothing", pd.Series([0.0])))
        filled.loc[index, "shift_arcsec"] = float(values.median())
        filled.loc[index, "is_fallback"] = True
        where = f"order {int(row['order'])} of {row['frame']}" if row["order"] != WHOLE_FRAME else row["frame"]
        log.warning(
            f"No reliable object trace in {where}, so its slit shift is the median of {source} "
            f"({filled.loc[index, 'shift_arcsec']:+.3f} arcsec)."
        )
    return filled


def shift_frame_along_slit(frame, orderShifts, mapImages):
    """*resample a frame through the 2D map, so its sources move along the slit at fixed wavelength*

    A pixel of order ``o`` takes its value from where the slit coordinate is ``orderShifts[o]`` arcsec higher at the
    same wavelength, so an object that sits ``orderShifts[o]`` above the reference position moves onto it. The data
    are interpolated bilinearly, with the edge value beyond the detector edge. A pixel is flagged bad when any pixel
    that contributes to it is bad, when its interpolated value or variance is not finite, when its source lies
    outside the order, or when its source lies off the detector. Orders with no shift, and pixels outside every
    order, are copied unchanged. A shifted value is multiplied by the ratio of the target pixel's area in the
    (wavelength, slit) plane to the source's, and its variance by the square of that ratio, so the object keeps its
    flux where the plate scale changes along the slit.

    The variance is interpolated bilinearly too, so a resampled pixel carries the weighted mean of the contributing
    variances, the sum of ``w * sigma**2``. Independent noise would give the smaller ``sum(w**2 * sigma**2)``. The
    resampled pixels are correlated with their neighbours, so the larger value is the conservative one, and it is
    the convention elsewhere in the pipeline.

    **Key Arguments:**

    - ``frame`` -- the CCDData frame to shift. It needs a standard-deviation uncertainty, or none
    - ``orderShifts`` -- the slit shift in arcsec of each order, as a dictionary from order number
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning

    **Return:**

    - ``shifted`` -- a new float32 CCDData frame. The input frame is not changed. Flags are not carried over, as
      nothing downstream of the alignment reads them
    """
    from astropy.nddata import CCDData, StdDevUncertainty

    if frame.data.shape != mapImages.order.shape:
        raise ValueError(f"the frame shape {frame.data.shape} differs from the 2D map shape {mapImages.order.shape}")
    hasUncertainty = frame.uncertainty is not None
    if hasUncertainty and frame.uncertainty.uncertainty_type != "std":
        raise ValueError(f"the frame uncertainty is {frame.uncertainty.uncertainty_type}, not std")

    # EVERY SOURCE IMAGE IS FLOAT32 AND MADE ONCE. THE NEW IMAGES ARE COPIES, SO THE SOURCES STAY UNCHANGED
    data = np.asarray(frame.data, dtype=np.float32)
    newData = data.copy()
    newVariance = variance = None
    if hasUncertainty:
        variance = np.asarray(frame.uncertainty.array, dtype=np.float32) ** 2
        newVariance = variance.copy()
    newIsBad = np.zeros(data.shape, dtype=bool) if frame.mask is None else np.asarray(frame.mask) != 0
    badness = newIsBad.astype(np.float32)

    height, width = data.shape
    pixelAreas = None
    for order, shift in orderShifts.items():
        if shift == 0:
            continue
        if pixelAreas is None:
            pixelAreas = _pixel_areas(mapImages)
        targets, (rows, columns) = _order_sources(mapImages, order, shift)
        sources = [rows, columns]
        isInBounds = (rows >= 0) & (rows <= height - 1) & (columns >= 0) & (columns <= width - 1)
        coverage = _interpolate((mapImages.order == order).astype(np.float32), sources, mode="constant")
        # A PIXEL HOLDS THE FLUX THAT FALLS ON ITS AREA OF THE (WAVELENGTH, SLIT) PLANE. THAT AREA CHANGES ALONG THE
        # SLIT, BY ABOUT 8% OVER A NOD THROW AT THE RED END OF THE SOXS NIR ORDER 10, SO THE VALUE IS RESCALED TO THE
        # TARGET'S AREA TO KEEP THE OBJECT FLUX (DY-1679)
        # A SOURCE OFF THE ORDER HAS NO AREA. ITS PIXEL IS FLAGGED BAD AND KEEPS THE UNSCALED VALUE, NOT A NaN.
        # NEXT TO THE ORDER EDGE THE INTERPOLATED SOURCE AREA MIXES IN ZERO AREA AND THE RATIO IS TOO LARGE; THOSE
        # PIXELS STAY SAFE ONLY BECAUSE THE COVERAGE TEST BELOW FLAGS THEM BAD
        with np.errstate(divide="ignore", invalid="ignore"):
            areaRatio = pixelAreas[targets] / _interpolate(pixelAreas, sources)
        hasArea = np.isfinite(areaRatio) & (areaRatio > 0)
        areaRatio = np.where(hasArea, areaRatio, np.float32(1.0))
        newData[targets] = _interpolate(data, sources) * areaRatio
        isBad = (_interpolate(badness, sources) > INTERPOLATION_TOLERANCE) | (coverage < 1.0 - INTERPOLATION_TOLERANCE)
        isBad |= ~isInBounds | ~np.isfinite(newData[targets]) | ~hasArea
        if variance is not None:
            newVariance[targets] = _interpolate(variance, sources) * areaRatio**2
            isBad |= ~np.isfinite(newVariance[targets])
        newIsBad[targets] = isBad

    uncertainty = None
    if hasUncertainty:
        uncertainty = StdDevUncertainty(np.sqrt(newVariance), unit=frame.uncertainty.unit)
    return CCDData(newData, unit=frame.unit, meta=copy.deepcopy(frame.meta), mask=newIsBad, uncertainty=uncertainty)


def _pixel_areas(mapImages):
    """*the area of each pixel in the (wavelength, slit) plane, in nm arcsec, from the map gradients*

    The area is the absolute determinant of the gradients. It is zero outside the orders, where the gradients are.
    """
    return np.abs(_gradient_determinant(*_gradient_images(mapImages))).astype(np.float32)


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

    The source is a step of ``shift`` in slit coordinate at fixed wavelength, taken along the local gradients of the
    wavelength and slit images. The plate scale can change across a nod throw, by about 8% at the red end of the SOXS
    NIR order 10, so one step from the target overshoots. The step is then repeated from the source, along the
    gradients there, for the slit coordinate still missing, ``SOURCE_REFINEMENTS`` times. Only the slit coordinate is
    refined: on a tilted, curved map the source can drift in wavelength by a second-order amount, 0.013 nm (about a
    tenth of a pixel) at the red end of the SOXS NIR order 10.

    **Return:**

    - ``targets`` -- the (rows, columns) index arrays of the order's pixels
    - ``sources`` -- the (rows, columns) source coordinates of those pixels
    """
    targets = np.nonzero(mapImages.order == order)
    gradients = [image[targets] for image in _gradient_images(mapImages)]
    rows, columns = _slit_step(targets[0], targets[1], gradients, shift)
    targetSlit = mapImages.slit[targets] + shift
    for _ in range(SOURCE_REFINEMENTS):
        sources = [rows, columns]
        missing = targetSlit - _interpolate(mapImages.slit, sources)
        gradients = [_interpolate(image, sources) for image in _gradient_images(mapImages)]
        newRows, newColumns = _slit_step(rows, columns, gradients, missing)
        # A SOURCE THAT CANNOT BE REFINED, AS IT LIES OFF THE ORDER, KEEPS ITS FIRST STEP. THE CALLER FLAGS IT BAD
        isRefined = np.isfinite(newRows) & np.isfinite(newColumns)
        rows, columns = np.where(isRefined, newRows, rows), np.where(isRefined, newColumns, columns)
    return targets, [rows, columns]


def _gradient_images(mapImages):
    """*the wavelength row, wavelength column, slit row and slit column gradient images, in that order*"""
    return (mapImages.wavelengthRow, mapImages.wavelengthCol, mapImages.slitRow, mapImages.slitCol)


def _slit_step(rows, columns, gradients, slitStep):
    """*the pixel coordinates a step of ``slitStep`` arcsec along the slit away, at fixed wavelength*

    **Key Arguments:**

    - ``rows``, ``columns`` -- the coordinates to step from
    - ``gradients`` -- the wavelength row, wavelength column, slit row and slit column gradients at those coordinates
    - ``slitStep`` -- the slit step in arcsec, a scalar or one value per coordinate

    **Return:**

    - ``rows``, ``columns`` -- the new coordinates
    """
    wavelengthRow, wavelengthCol, _, _ = gradients
    determinant = _gradient_determinant(*gradients)
    # A ZERO DETERMINANT GIVES AN INFINITE OR NaN SOURCE. THE CALLER FLAGS THAT PIXEL BAD, AND
    # _safe_sources KEEPS THE NON-FINITE VALUE FROM BEING READ AS A MEMORY OFFSET
    with np.errstate(divide="ignore", invalid="ignore"):
        return rows - wavelengthCol * slitStep / determinant, columns + wavelengthRow * slitStep / determinant


def read_two_d_map_images(twoDMapPath, binx, biny, kw):
    """*read a 2D dispersion map at a frame's binning and build its slit-alignment images*

    **Key Arguments:**

    - ``twoDMapPath`` -- the path to the 2D dispersion map image
    - ``binx`` -- the frame's binning in the x-axis
    - ``biny`` -- the frame's binning in the y-axis
    - ``kw`` -- the fits keyword lookup, used to read the map's own binning

    **Return:**

    - ``mapImages`` -- the ``two_d_map_images`` of the map
    """
    from soxspipe.commonutils.toolkit import _open_disp_map_hdus

    hdul, _, _ = _open_disp_map_hdus(twoDMapPath, binx, biny, kw=kw)
    try:
        return build_two_d_map_images(*(hdul[name].data for name in ("WAVELENGTH", "SLIT", "ORDER")))
    finally:
        hdul.close()


def measure_trace_slit_positions(
    log, differenceFrame, mapImages, dispersionMap, settings, recipeSettings, recipeName, startNightDate=""
):
    """*find the object trace on a cycle difference and read its slit position in each order*

    The trace is found with ``detect_continuum.sample_trace`` only. No order table is fitted, and no product or QC is
    written. A trace that cannot be found at all leaves every order without a position.

    **Key Arguments:**

    - ``log`` -- logger
    - ``differenceFrame`` -- the cycle difference (this frame minus its partner), with the frame's header
    - ``mapImages`` -- the ``two_d_map_images`` of the frame's binning
    - ``dispersionMap`` -- the path to the dispersion map table
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe settings, with the ``detect-continuum`` block
    - ``recipeName`` -- the recipe name
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night, which names the QC and product folders

    **Return:**

    - ``positions`` -- the output of ``order_slit_positions``
    """
    from soxspipe.commonutils.detect_continuum import detect_continuum

    detector = detect_continuum(
        log=log,
        traceFrame=differenceFrame,
        dispersion_map=dispersionMap,
        settings=settings,
        recipeSettings=recipeSettings,
        recipeName=recipeName,
        qcTable=pd.DataFrame(),
        productsTable=False,
        startNightDate=startNightDate,
    )
    try:
        orderPixelTable, detectionPercent = detector.sample_trace()
    except ValueError as error:
        log.warning(f"The object trace could not be sampled ({error}).")
        orderPixelTable, detectionPercent = None, 0.0
    if orderPixelTable is None or detectionPercent < MIN_DETECTION_PERCENT:
        return order_slit_positions(_no_trace_samples(), mapImages)
    return order_slit_positions(orderPixelTable, mapImages)


def _no_trace_samples():
    """*an empty trace-sample table with numeric columns, which ``order_slit_positions`` can test for finite values*"""
    return pd.DataFrame({column: pd.Series(dtype=float) for column in TRACE_COLUMNS})


@dataclass(frozen=True)
class slit_alignment_result:
    """The frames of a group after slit alignment, and the shifts that were applied to them."""

    frames: tuple
    shifts: pd.DataFrame
    fallbackCount: int
    mapImages: two_d_map_images | None = None


@dataclass(frozen=True)
class frame_group:
    """A group of frames to align to each other, with the partner frames that are subtracted to trace them.

    ``frames`` is a list, and ``align_frame_groups_along_slit`` replaces its items with the aligned frames.
    """

    frames: list
    partnerFrames: list
    frameNames: list | None = None


def align_frames_along_slit(
    log,
    frames,
    partnerFrames,
    twoDMapPath,
    dispersionMap,
    settings,
    recipeSettings,
    recipeName,
    perOrder,
    frameNames=None,
    startNightDate="",
):
    """*shift each frame of one group along the slit so that its object trace lies where the first frame's does*

    This is ``align_frame_groups_along_slit`` for a single group. It leaves the ``frames`` list as it is and returns
    the aligned frames, so the caller holds both the original and the aligned frames.

    **Key Arguments:**

    - ``log`` -- logger
    - ``frames`` -- the CCDData frames to align, with an observation-time keyword
    - ``partnerFrames`` -- the frames to subtract from them, as many as ``frames``
    - ``twoDMapPath`` -- the path to the 2D dispersion map image
    - ``dispersionMap`` -- the path to the dispersion map table
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe settings, with the ``detect-continuum`` block
    - ``recipeName`` -- the recipe name
    - ``perOrder`` -- if True, every order of a frame has its own shift. If False, a frame has one shift
    - ``frameNames`` -- the name of each frame, for the log and the shift table. Default *None*, for the ``ARCFILE``
      or ``ORIGFILE`` keyword, or a number
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""

    **Return:**

    - ``result`` -- a ``slit_alignment_result``: the aligned frames in the order of ``frames``, the shift table of
      ``slit_shifts`` and the number of shifts that were not measured
    """
    group = frame_group(list(frames), partnerFrames, frameNames)
    return align_frame_groups_along_slit(
        log, [group], twoDMapPath, dispersionMap, settings, recipeSettings, recipeName, perOrder, startNightDate
    )[0]


def align_frame_groups_along_slit(
    log,
    groups,
    twoDMapPath,
    dispersionMap,
    settings,
    recipeSettings,
    recipeName,
    perOrder,
    startNightDate="",
    measureSingleFrames=False,
):
    """*align each group of frames along the slit to the first frame of the group in time*

    Each frame's trace is measured on its cycle difference, the frame minus its partner in the other nod position.
    Frames and partners are paired in time order. One difference is held in memory at a time. All the traces are
    measured before any frame is resampled, because a partner frame must be unshifted when it is subtracted. The
    frames are then resampled one at a time, and each is replaced in its group's ``frames`` list by its aligned
    frame, so an original is released as soon as its replacement exists. A group of one frame is left as it is,
    with a zero shift.

    **Key Arguments:**

    - ``log`` -- logger
    - ``groups`` -- the ``frame_group`` objects to align. Their ``frames`` lists are changed in place
    - ``twoDMapPath`` -- the path to the 2D dispersion map image
    - ``dispersionMap`` -- the path to the dispersion map table
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe settings, with the ``detect-continuum`` block
    - ``recipeName`` -- the recipe name
    - ``perOrder`` -- if True, every order of a frame has its own shift. If False, a frame has one shift
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
    - ``measureSingleFrames`` -- if True, a group of one frame has its trace position measured too, so its shift
      table holds the position of its one frame. The frame is not resampled and no shift counts as a fallback.
      Default *False*

    **Return:**

    - ``results`` -- a ``slit_alignment_result`` for each group: its aligned frames, the shift table of
      ``slit_shifts``, the number of shifts that were not measured, and the 2D map images when the group was
      measured (else *None*)
    """
    from soxspipe.commonutils.keyword_lookup import keyword_lookup

    kw = keyword_lookup(log=log, settings=settings).get
    mjdKey = kw("MJDOBS")
    for group in groups:
        if len(group.frames) != len(group.partnerFrames):
            raise ValueError(f"{len(group.frames)} frames cannot be paired with {len(group.partnerFrames)} partners")
    names = [_frame_names(group.frames, group.frameNames) for group in groups]
    isAligned = [len(group.frames) > 1 for group in groups]
    isMeasured = [True] * len(groups) if measureSingleFrames else isAligned

    mapImages, measured = _measure_groups(
        log, groups, names, isMeasured, (twoDMapPath, dispersionMap, settings, recipeSettings, recipeName),
        perOrder, startNightDate, kw,
    )  # fmt: skip
    return _resample_groups(groups, names, isAligned, measured, mapImages, mjdKey)


def _measure_groups(log, groups, names, isMeasured, inputs, perOrder, startNightDate, kw):
    """*measure the trace position and the slit shifts of every group that is to be measured*

    **Key Arguments:**

    - ``groups`` -- the ``frame_group`` objects
    - ``names`` -- the frame names of each group
    - ``isMeasured`` -- True for each group to measure
    - ``inputs`` -- the 2D map path, the dispersion map path, the settings, the recipe settings and the recipe name
    - ``perOrder`` -- if True, every order of a frame has its own shift
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night
    - ``kw`` -- the fits keyword lookup

    **Return:**

    - ``mapImages`` -- the ``two_d_map_images`` read for the measurement, or *None* when no group is measured
    - ``measured`` -- for each group, ``(shifts, fallbackCount)``, or *None* for a group that was not measured
    """
    twoDMapPath, dispersionMap, settings, recipeSettings, recipeName = inputs
    mapImages, measured = None, [None] * len(groups)
    if not any(isMeasured):
        return mapImages, measured
    firstGroup = groups[isMeasured.index(True)]
    binx, biny = _frame_binning(log, firstGroup.frames[0], kw)
    mapImages = read_two_d_map_images(twoDMapPath, binx, biny, kw)
    for index, group in enumerate(groups):
        if not isMeasured[index]:
            continue
        positions = _measure_group(
            log, group, names[index], kw("MJDOBS"), mapImages, dispersionMap, settings, recipeSettings, recipeName,
            startNightDate,
        )  # fmt: skip
        shifts, fallbackCount = slit_shifts(log, positions, perOrder)
        _log_shifts(log, shifts)
        measured[index] = (shifts, fallbackCount)
    return mapImages, measured


def _resample_groups(groups, names, isAligned, measured, mapImages, mjdKey):
    """*resample each group of more than one frame by its shifts, and build the result of every group*

    Each aligned frame replaces its original in the group's ``frames`` list.

    **Return:**

    - ``results`` -- a ``slit_alignment_result`` for each group
    """
    results = []
    for index, group in enumerate(groups):
        if not isAligned[index]:
            results.append(_unaligned_result(group, names[index], mjdKey, measured[index], mapImages))
            continue
        shifts, fallbackCount = measured[index]
        for frameIndex in range(len(group.frames)):
            orderShifts = _order_shifts(shifts, names[index][frameIndex], mapImages.orders)
            group.frames[frameIndex] = shift_frame_along_slit(group.frames[frameIndex], orderShifts, mapImages)
        results.append(slit_alignment_result(tuple(group.frames), shifts, fallbackCount, mapImages))
    return results


def _unaligned_result(group, names, mjdKey, measurement, mapImages):
    """*the result for a group of one frame: the frame as it is, with its measured positions when it was measured*

    A group of one frame has nothing to align to, so none of its shifts count as a fallback.
    """
    if measurement is None:
        return slit_alignment_result(tuple(group.frames), _unshifted(group.frames, names, mjdKey), 0)
    return slit_alignment_result(tuple(group.frames), measurement[0], 0, mapImages)


def _measure_group(
    log, group, names, mjdKey, mapImages, dispersionMap, settings, recipeSettings, recipeName, startNightDate
):
    """*the slit position of the trace in each order of each frame of a group, from the frame minus its partner*"""
    partnerOf = _partner_indices(group.frames, group.partnerFrames, mjdKey)
    measured = []
    for index, frame in enumerate(group.frames):
        difference = frame.subtract(group.partnerFrames[partnerOf[index]])
        difference.header = frame.header
        positions = measure_trace_slit_positions(
            log, difference, mapImages, dispersionMap, settings, recipeSettings, recipeName, startNightDate
        )
        del difference
        measured.append(positions.assign(frame=names[index], mjd=float(frame.header[mjdKey])))
    return pd.concat(measured, ignore_index=True)


def _frame_names(frames, frameNames):
    """*one distinct name per frame: the given name, else the ARCFILE or ORIGFILE keyword, else a number*"""
    names = []
    for index, frame in enumerate(frames):
        name = frameNames[index] if frameNames is not None else None
        name = name or frame.header.get("ARCFILE") or frame.header.get("ORIGFILE") or f"frame {index + 1}"
        names.append(name if names.count(name) == 0 else f"{name} ({index + 1})")
    return names


def _unshifted(frames, names, mjdKey):
    """*the shift table of a group that needs no alignment: one zero whole-frame shift per frame*"""
    return pd.DataFrame(
        {
            "frame": names,
            "mjd": [float(frame.header[mjdKey]) for frame in frames],
            "order": WHOLE_FRAME,
            "position_arcsec": np.nan,
            "shift_arcsec": 0.0,
            "is_fallback": False,
            "is_reference": True,
        }
    )[SHIFT_COLUMNS]


def _frame_binning(log, frame, kw):
    """*the frame's binning from its header: 1 by 1 on the NIR arm, or with a warning if either key is unusable*"""
    if frame.header[kw("SEQ_ARM")] == "NIR":
        return 1, 1
    try:
        return int(frame.header[kw("WIN_BINX")]), int(frame.header[kw("WIN_BINY")])
    except (KeyError, TypeError, ValueError) as error:
        log.warning(
            f"The frame header has no usable integer binning keywords ({type(error).__name__}: {error}), "
            "so a binning of 1 by 1 is used to read the 2D map."
        )
        return 1, 1


def _partner_indices(frames, partnerFrames, mjdKey):
    """*for each frame, the index of the partner frame that is in the same place in time order*"""
    frameOrder = sorted(range(len(frames)), key=lambda index: frames[index].header[mjdKey])
    partnerOrder = sorted(range(len(partnerFrames)), key=lambda index: partnerFrames[index].header[mjdKey])
    return dict(zip(frameOrder, partnerOrder, strict=True))


def _order_shifts(shifts, name, orders):
    """*the shift of each order of a frame, from its rows of the shift table*"""
    rows = shifts.loc[shifts["frame"] == name]
    if (rows["order"] == WHOLE_FRAME).all():
        return {order: float(rows["shift_arcsec"].iloc[0]) for order in orders}
    return dict(zip(rows["order"].astype(int), rows["shift_arcsec"], strict=True))


def _log_shifts(log, shifts):
    """*print the slit shift of every frame, in time order*"""
    for name, rows in shifts.sort_values("mjd", kind="stable").groupby("frame", sort=False):
        parts = [
            ("" if order == WHOLE_FRAME else f"order {order}: ") + f"{shift:+.3f}"
            for order, shift in zip(rows["order"], rows["shift_arcsec"], strict=True)
        ]
        log.print(f"Slit shift of {name} (arcsec): " + ", ".join(parts))
