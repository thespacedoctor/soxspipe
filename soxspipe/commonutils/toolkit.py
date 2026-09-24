#!/usr/bin/env python
"""
*small reusable functions used throughout soxspipe*

Author
: David Young

Date Created
: September 18, 2020
"""

import logging
import os
import sys
from datetime import UTC, datetime
from os.path import expanduser

from soxspipe.commonutils import detector_lookup, keyword_lookup
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import (
    dispersion_map_to_pixel_arrays,
)
from soxspipe.commonutils.polynomials import (
    chebyshev_order_xy_polynomials,
)

os.environ["TERM"] = "vt100"


def cut_image_slice(log, frame, width, length, x, y, sliceAxis="x", median=False, debug=False):
    """*cut and return an N-pixel wide and M-pixels long slice, centred on a given coordinate from an image frame*

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the data array to cut the slice from (masked array)
    - ``width`` -- width of the slice (odd number)
    - ``length`` -- length of the slice
    - ``x`` -- x-coordinate
    - ``y`` -- y-coordinate
    - ``sliceAxis`` -- the axis along which slice is to be taken. Default *x*
    - ``median`` -- collapse the slice to a median value across its width
    - ``debug`` -- generate a plot of slice. Useful for debugging.

    **Return:**

    - ``slice`` -- the median-collapsed slice when ``median`` is True
    - ``slice_length_offset`` -- the pixel offset of the slice start along its length
    - ``slice_width_centre`` -- the pixel coordinate of the slice centre across its width

    All three are *None* when the slice would fall outside the frame.

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import cut_image_slice
    slice, slice_length_offset, slice_width_centre = cut_image_slice(
        log=self.log, frame=self.pinholeFlat.data, width=1, length=sliceLength, x=x_fit, y=y_fit, median=True
    )
    if slice is None:
        return None
    ```
    """
    log.debug("starting the ``cut_image_slice`` function")

    import random

    import numpy as np
    import numpy.ma as ma

    halfSlice = length / 2
    # NEED AN EVEN PIXEL SIZE
    halfwidth = (width - 1) / 2 if width % 2 != 0 else width / 2

    if sliceAxis == "x":
        axisA = x
        axisB = y
        axisALen = frame.shape[1]
        axisBLen = frame.shape[0]
    elif sliceAxis == "y":
        axisB = x
        axisA = y
        axisALen = frame.shape[0]
        axisBLen = frame.shape[1]
    else:
        raise ValueError("sliceAxis needs to be either 'x' or 'y'")

    # CHECK WE ARE NOT GOING BEYOND BOUNDS OF FRAME
    if (axisA > axisALen - halfSlice) or (axisB > axisBLen - halfwidth) or (axisA < halfSlice) or (axisB < halfwidth):
        return None, None, None

    slice_length_offset = int(axisA - halfSlice)
    if sliceAxis == "x":
        sliceFull = frame[
            int(axisB - halfwidth) : int(axisB + halfwidth + 1),
            slice_length_offset : int(axisA + halfSlice),
        ]
    else:
        sliceFull = frame[
            slice_length_offset : int(axisA + halfSlice),
            int(axisB - halfwidth) : int(axisB + halfwidth + 1),
        ]
    slice_width_centre = (int(axisB + halfwidth + 1) + int(axisB - halfwidth)) / 2

    # # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
    # maskedDataArray = np.ma.array(sliceFull.data, mask=sliceFull.mask)
    # try:
    #     sliceFull.data=sliceFull.data - np.percentile(maskedDataArray, 30)
    # except:
    #     pass

    if median:
        slice = ma.median(sliceFull, axis=1) if sliceAxis == "y" else ma.median(sliceFull, axis=0)

    # DELIBERATELY DISABLED DEBUG PLOT: THE LEADING `False` SHORT-CIRCUITS, SO `random` NEVER RUNS
    if False and debug and random.randint(1, 101) < 5:  # noqa: SIM223, S311
        import matplotlib.pyplot as plt

        # CHECK THE SLICE POINTS IF NEEDED
        sliceImg = np.rot90(sliceFull, 1) if sliceAxis == "y" else sliceFull
        plt.imshow(sliceImg)
        plt.show()
        xx = np.arange(0, len(slice))
        plt.figure(figsize=(8, 5))
        if sliceAxis == "y":
            plt.plot(xx, slice, "ko", label=f"x={axisB}, y={axisA}, sliceAxis={sliceAxis}")
        if sliceAxis == "x":
            plt.plot(xx, slice, "ko", label=f"x={axisA}, y={axisB}, sliceAxis={sliceAxis}")
        plt.xlabel("Position")
        plt.ylabel("Flux")
        plt.legend()
        plt.show()

    log.debug("completed the ``cut_image_slice`` function")
    return slice, slice_length_offset, slice_width_centre


def quicklook_image(
    log,
    CCDObject,
    show=True,
    ext="data",
    stdWindow=3,
    title=False,
    surfacePlot=False,
    dispMap=False,
    dispMapImage=False,
    inst=False,
    settings=False,
    skylines=False,
    saveToPath=False,
):
    """*generate a quicklook image of a CCDObject - useful for development/debugging*

    **Key Arguments:**

    - ``log`` -- logger
    - ``CCDObject`` -- the CCDObject to plot
    - ``show`` -- show the image. Set to False to skip
    - ``ext`` -- the name of the the extension to show. Can be "data", "mask" or "err". Default "data".
    - ``stdWindow`` -- the width of the colour scale in standard deviations, centred on the median. Default *3*
    - ``title`` -- give a title for the plot
    - ``surfacePlot`` -- plot as a 3D surface plot
    - ``dispMap`` -- path to dispersion map. Default *False*
    - ``dispMapImage`` -- the 2D dispersion map image
    - ``inst`` -- provide instrument name if no header exists
    - ``settings`` -- the soxspipe settings dictionary, used to look up the arm and skylines. Default *False*
    - ``skylines`` -- mark skylines on image
    - ``saveToPath`` -- path to save the plot to. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import quicklook_image
    quicklook_image(
        log=self.log, CCDObject=myframe, show=True)
    ```
    """
    log.debug("starting the ``quicklook_image`` function")

    if not show and not saveToPath:
        return

    from copy import copy

    import matplotlib as mpl

    originalRC = dict(mpl.rcParams)
    import matplotlib.pyplot as plt

    if settings:
        kw, arm = _quicklook_header_lookups(log, CCDObject, settings)

    frame = _quicklook_frame_array(CCDObject, ext)

    if inst is False:
        inst = _quicklook_instrument(log, CCDObject)

    skylinesDF = get_skylines_dataframe(log, settings, arm) if skylines else False

    # COMBINE MASK WITH THE BAD PIXEL MASK
    if not isinstance(dispMapImage, bool):
        gridLinePixelTable = _apply_inter_order_mask(
            log=log,
            frame=frame,
            CCDObject=CCDObject,
            dispMap=dispMap,
            dispMapImage=dispMapImage,
            kw=kw,
            skylinesDF=skylinesDF,
        )

    rotatedImg = _rotate_for_display(frame, inst)

    from astropy.stats import sigma_clipped_stats

    mean, median, std = sigma_clipped_stats(frame, sigma=50.0, stdfunc="mad_std", cenfunc="median", maxiters=3)

    # std = np.nanstd(frame)
    # mean = np.nanmean(frame)
    # median = np.nanmean(frame)
    palette = copy(plt.cm.viridis)
    palette.set_bad("#dc322f", 1.0)
    vmax = median + stdWindow * 0.5 * std
    vmin = median - stdWindow * 0.5 * std

    if surfacePlot:
        fig, ax2 = _draw_surface_plot(rotatedImg, frame, inst, vmin, vmax)
    else:
        if rotatedImg.shape[0] - rotatedImg.shape[1] > 1000:
            fig = plt.figure(figsize=(5, 12))
        else:
            fig = plt.figure(figsize=(12, 5))
        # palette.set_over('r', 1.0)
        # palette.set_under('g', 1.0)
        ax2 = fig.add_subplot(111)

    if not isinstance(dispMapImage, bool):
        _draw_dispersion_grid_lines(ax2, gridLinePixelTable, inst)

    _draw_detector_image(
        fig=fig,
        ax2=ax2,
        rotatedImg=rotatedImg,
        vmin=vmin,
        vmax=vmax,
        palette=palette,
        mean=mean,
        surfacePlot=surfacePlot,
        title=title,
        inst=inst,
    )

    if show:
        # plt.pause(0.1)
        plt.show()

    if saveToPath:
        save_qc_plot(saveToPath)
        plt.clf()  # CLEAR FIGURE
    mpl.rcParams.update(originalRC)
    plt.close("all")

    log.debug("completed the ``quicklook_image`` function")
    return


def _quicklook_header_lookups(log, CCDObject, settings):
    """*build the keyword lookup for a quicklook frame and read its arm, validating the lookups on the way*

    **Key Arguments:**

    - ``log`` -- logger
    - ``CCDObject`` -- the CCDObject to plot
    - ``settings`` -- the soxspipe settings dictionary

    **Return:**

    - ``kw`` -- the fits keyword lookup function
    - ``arm`` -- the spectrograph arm the frame was taken with
    """
    # RESOLVED PER CALL, AS `quicklook_image` DID: THE MODULE-LEVEL NAMES ARE BOUND AT IMPORT AND WOULD NOT
    # SEE A LOOKUP PATCHED ON THE PACKAGE
    from soxspipe.commonutils import detector_lookup, keyword_lookup

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(log=log, settings=settings).get
    arm = CCDObject.header[kw("SEQ_ARM")]
    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING KEYWORD; DELETING IT REMOVES THAT FAILURE
    dateObs = CCDObject.header[kw("DATE_OBS")]  # noqa: F841

    # DETECTOR PARAMETERS LOOKUP OBJECT
    detectorParams = detector_lookup(log=log, settings=settings).get(arm)

    # USE THIS ELSEWHERE IN THE OBJECT METHODS
    dp = detectorParams
    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING PARAMETER; DELETING IT REMOVES THAT FAILURE
    science_pixels = dp["science-pixels"]  # noqa: F841

    return kw, arm


def _apply_inter_order_mask(log, frame, CCDObject, dispMap, dispMapImage, kw, skylinesDF):
    """*mask a quicklook frame's inter-order pixels and return the dispersion grid lines to draw*

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the array being plotted, masked in place where it supports a mask
    - ``CCDObject`` -- the CCDObject the frame came from
    - ``dispMap`` -- path to dispersion map
    - ``dispMapImage`` -- the 2D dispersion map image
    - ``kw`` -- the fits keyword lookup function
    - ``skylinesDF`` -- the skylines dataframe, or *False* for none

    **Return:**

    - ``gridLinePixelTable`` -- the pixel coordinates of the dispersion solution grid lines
    """
    gridLinePixelTable, interOrderMask = create_dispersion_solution_grid_lines_for_plot(
        log=log,
        dispMap=dispMap,
        dispMapImage=dispMapImage,
        associatedFrame=CCDObject,
        kw=kw,
        skylines=skylinesDF,
    )

    try:
        mask = (frame.mask == 1) | (interOrderMask == 1)
    except (AttributeError, ValueError) as e:
        log.debug(f"quicklook_image: `mask = (frame.mask == 1) | (interOrderMask == 1)` failed, continuing: {e}")
        mask = interOrderMask == 1
    try:
        frame.mask = mask
    except AttributeError as e:
        log.debug(f"quicklook_image: `frame.mask = mask` failed, continuing: {e}")

    return gridLinePixelTable


def _draw_detector_image(fig, ax2, rotatedImg, vmin, vmax, palette, mean, surfacePlot, title, inst):
    """*draw the detector image, its colour bar and its title onto the quicklook figure*

    **Key Arguments:**

    - ``fig`` -- the quicklook figure
    - ``ax2`` -- the axes to draw the image on
    - ``rotatedImg`` -- the image array, rotated for display
    - ``vmin`` -- the lower limit of the colour scale
    - ``vmax`` -- the upper limit of the colour scale
    - ``palette`` -- the colour map
    - ``mean`` -- the frame's sigma-clipped mean, which sets the colour-bar number format
    - ``surfacePlot`` -- is the figure a 3D surface plot?
    - ``title`` -- the figure title, or *False* for none
    - ``inst`` -- the instrument name
    """
    import matplotlib.pyplot as plt

    if rotatedImg.shape[0] - rotatedImg.shape[1] > 1000:
        ax2.set_box_aspect(2.0)
    else:
        ax2.set_box_aspect(0.5)
    detectorPlot = plt.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap=palette, alpha=1, aspect="auto")

    shrink = 0.5 if surfacePlot else 1.0

    if mean > 10:
        fmt = "%1.0f"
        fig.colorbar(detectorPlot, shrink=shrink, format=fmt)
    else:
        fig.colorbar(detectorPlot, shrink=shrink)
        # plt.colorbar()
    if title:
        fig.suptitle(title, fontsize=20)
    if inst == "XSHOOTER":
        ax2.invert_yaxis()
    # cbar.ticklabel_format(useOffset=False)
    _label_detector_axes(inst)


def _quicklook_frame_array(CCDObject, ext):
    """*return the array ``quicklook_image`` plots for the requested extension*

    **Key Arguments:**

    - ``CCDObject`` -- the CCDObject (or plain array) to plot
    - ``ext`` -- the extension name: "data", "mask" or "uncertainty". Anything else treats
      ``CCDObject`` as a plain array

    **Return:**

    - ``frame`` -- the array to plot
    """
    if ext == "data":
        frame = CCDObject.data
    elif ext == "mask":
        frame = CCDObject.mask
    elif ext == "uncertainty":
        frame = CCDObject.uncertainty.array
    else:
        # ASSUME ONLY NDARRAY
        frame = CCDObject
    return frame


def _quicklook_instrument(log, CCDObject):
    """*read the instrument name from the frame header, falling back to XSHOOTER*

    **Key Arguments:**

    - ``log`` -- logger
    - ``CCDObject`` -- the CCDObject (or plain array) being plotted

    **Return:**

    - ``inst`` -- the instrument name
    """
    try:
        inst = CCDObject.header["INSTRUME"]
    except (KeyError, AttributeError) as e:
        log.debug(f"quicklook_image: `inst = CCDObject.header['INSTRUME']` failed, continuing: {e}")
        inst = "XSHOOTER"
    return inst


def _rotate_for_display(frame, inst):
    """*orient a detector frame so it displays with the instrument's axis convention*

    **Key Arguments:**

    - ``frame`` -- the array to orient
    - ``inst`` -- the instrument name

    **Return:**

    - ``rotatedImg`` -- the oriented array
    """
    import numpy as np

    if inst == "SOXS":
        rotatedImg = np.flipud(frame)
    elif inst == "XSHOOTER":
        rotatedImg = np.rot90(frame, 1)
    else:
        rotatedImg = frame
    return np.flipud(rotatedImg)


def _label_detector_axes(inst):
    """*label the current axes with the detector axis names for the instrument*

    **Key Arguments:**

    - ``inst`` -- the instrument name
    """
    import matplotlib.pyplot as plt

    if inst == "SOXS":
        plt.xlabel("x-axis", fontsize=16)
        plt.ylabel("y-axis", fontsize=16)
    else:
        plt.xlabel("y-axis", fontsize=16)
        plt.ylabel("x-axis", fontsize=16)


def _draw_surface_plot(rotatedImg, frame, inst, vmin, vmax):
    """*draw the 3D surface panel of a ``quicklook_image`` surface plot*

    **Key Arguments:**

    - ``rotatedImg`` -- the oriented array to draw
    - ``frame`` -- the unrotated array, used for the z-axis limit
    - ``inst`` -- the instrument name
    - ``vmin`` -- the lower colour and z-axis limit
    - ``vmax`` -- the upper colour limit

    **Return:**

    - ``fig`` -- the new figure
    - ``ax2`` -- the empty 2D axes beside the surface, for the detector image
    """
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import rc

    axisColour = "#002b36"
    rc("axes", edgecolor=axisColour, labelcolor=axisColour, linewidth=0.6)
    rc("xtick", color=axisColour)
    rc("ytick", color=axisColour)
    rc("grid", color=axisColour)
    rc("text", color=axisColour)

    fig = plt.figure(figsize=(20, 8))
    ax = fig.add_subplot(121, projection="3d")
    if inst == "XSHOOTER":
        plt.gca().invert_yaxis()
    ax.set_box_aspect(aspect=(2, 1, 1))
    # REMOVE GRAY PANES AND AXIS GRID
    ax.xaxis.pane.fill = False
    ax.zaxis.pane.set_facecolor("#dc322f")
    ax.zaxis.pane.set_alpha(1.0)
    ax.yaxis.pane.fill = False

    ax.grid(False)
    # REMOVE Z-AXIS
    # ax.w_zaxis.line.set_lw(0.)
    # ax.set_zticks([])

    X, Y = np.meshgrid(
        np.linspace(0, rotatedImg.shape[1], rotatedImg.shape[1]),
        np.linspace(0, rotatedImg.shape[0], rotatedImg.shape[0]),
    )
    ax.plot_surface(
        X=X,
        Y=Y,
        Z=rotatedImg,
        cmap="viridis",
        antialiased=True,
        vmin=vmin,
        vmax=vmax,
    )

    if inst == "SOXS":
        ax.azim = 70
    else:
        ax.azim = -120
    ax.elev = 30

    ax.set_xlim(0, rotatedImg.shape[1])
    ax.set_ylim(0, rotatedImg.shape[0])
    ax.set_zlim(vmin, min(np.nanmax(frame), vmax * 1.2))

    if inst == "SOXS":
        ax.invert_yaxis()
    backgroundColour = "white"
    fig.set_facecolor(backgroundColour)
    ax.set_facecolor(backgroundColour)
    ax.xaxis.pane.set_edgecolor(backgroundColour)
    ax.yaxis.pane.set_edgecolor(backgroundColour)
    ax.zaxis.pane.set_edgecolor(backgroundColour)

    _label_detector_axes(inst)

    ax2 = fig.add_subplot(122)
    return fig, ax2


def _draw_dispersion_grid_lines(ax2, gridLinePixelTable, inst):
    """*draw the dispersion-solution grid lines over a ``quicklook_image`` detector image*

    **Key Arguments:**

    - ``ax2`` -- the axes holding the detector image
    - ``gridLinePixelTable`` -- the grid-line pixel table from
      ``create_dispersion_solution_grid_lines_for_plot``
    - ``inst`` -- the instrument name
    """
    for lineIndex in range(int(gridLinePixelTable["line"].max())):
        mask = gridLinePixelTable["line"] == lineIndex
        if inst == "SOXS":
            ax2.plot(
                gridLinePixelTable.loc[mask]["fit_x"],
                gridLinePixelTable.loc[mask]["fit_y"],
                "w-",
                linewidth=0.5,
                alpha=0.8,
                color="black",
            )
        else:
            ax2.plot(
                gridLinePixelTable.loc[mask]["fit_y"],
                gridLinePixelTable.loc[mask]["fit_x"],
                "w-",
                linewidth=0.5,
                alpha=0.8,
                color="black",
            )


def unpack_order_table(
    log,
    orderTablePath,
    extend=0.0,
    pixelDelta=1,
    binx=1,
    biny=1,
    prebinned=False,
    order=False,
    limitToDetectorFormat=False,
):
    """*Unpack an order location table into polynomial, pixel and metadata dataframes.*

    Return an ``orderPolyTable`` dataframe containing the polynomial coefficients for the order centres and edges, an
    ``orderPixelTable`` dataframe containing the pixel-coordinates for each order centre and edges, and finally, an
    ``orderMetaTable`` dataframe giving metadata about the frame binning and format.

    **Key Arguments:**

    - ``log`` -- logger
    - ``orderTablePath`` -- path to the order table
    - ``extend`` -- fractional increase to the order area in the y-axis (needed for masking)
    - ``pixelDelta`` -- space between returned data points. Default *1* (sampled at every pixel)
    - ``binx`` -- binning in the x-axis (from FITS header). Default *1*
    - ``biny`` -- binning in the y-axis (from FITS header). Default *1*
    - ``prebinned`` -- was the order-table measured on a pre-binned frame (typically only for mflats). Default *False*
    - ``order`` -- unpack only a single order
    - ``limitToDetectorFormat`` -- limit the pixels return to those limited by the detector format static
      calibration table

    **Return:**

    - ``orderPolyTable`` -- the polynomial coefficients for the order centres and edges
    - ``orderPixelTable`` -- the pixel-coordinates for each order centre and edges
    - ``orderMetaTable`` -- metadata about the frame binning and format

    **Usage:**

    ```python
    # UNPACK THE ORDER TABLE
    from soxspipe.commonutils.toolkit import unpack_order_table
    orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
        log=self.log, orderTablePath=orderTablePath, extend=0.)
    ```
    """
    log.debug("starting the ``functionName`` function")
    from astropy.table import Table

    # PIXEL DELTA NEEDS TO BE ODD .. ELSE MASKING ON BINNED DATA GETS MESSED UP
    if pixelDelta % 2 == 0:
        pixelDelta += 1  # RETURN THE NEAREST ODD NUMBER ABOVE IF IT'S EVEN

    # MAKE RELATIVE HOME PATH ABSOLUTE

    home = expanduser("~")
    orderTablePath = orderTablePath.replace("~", home)

    dat = Table.read(orderTablePath, format="fits", hdu=1)
    orderPolyTable = dat.to_pandas()

    dat = Table.read(orderTablePath, format="fits", hdu=2)
    orderMetaTable = dat.to_pandas()

    if order:
        mask = orderMetaTable["order"] == order
        orderMetaTable = orderMetaTable.loc[mask]

    axisA, axisB, axisAbin, axisBbin = _order_table_axes(orderPolyTable, binx, biny)

    orderPixelTable = _order_table_pixel_grid(
        orderMetaTable=orderMetaTable,
        axisB=axisB,
        ratio=axisBbin if prebinned else 1,
        extend=extend,
        pixelDelta=pixelDelta,
    )

    orderPixelTable[f"{axisA}coord_centre"] = _evaluate_order_polynomial(
        log=log,
        orderPolyTable=orderPolyTable,
        orderPixelTable=orderPixelTable,
        axisB=axisB,
        degreeSuffix="cent",
        coefficients=_order_polynomial_coefficients(orderPolyTable, "cent_"),
    )

    std_coeff = _order_polynomial_coefficients(orderPolyTable, "std_")
    if len(std_coeff):
        # THE STANDARD-DEVIATION POLYNOMIAL DELIBERATELY REUSES THE CENTRE-TRACE DEGREES
        orderPixelTable["std"] = _evaluate_order_polynomial(
            log=log,
            orderPolyTable=orderPolyTable,
            orderPixelTable=orderPixelTable,
            axisB=axisB,
            degreeSuffix="cent",
            coefficients=std_coeff,
        )

    for edge in ["edgeup", "edgelow"]:
        if f"deg{axisB}_{edge}" in orderPolyTable.columns:
            orderPixelTable[f"{axisA}coord_{edge}"] = _evaluate_order_polynomial(
                log=log,
                orderPolyTable=orderPolyTable,
                orderPixelTable=orderPixelTable,
                axisB=axisB,
                degreeSuffix=edge,
                coefficients=_order_polynomial_coefficients(orderPolyTable, f"{edge}_"),
            )

    orderPixelTable, orderMetaTable = _rescale_order_tables_for_binning(
        orderPixelTable=orderPixelTable,
        orderMetaTable=orderMetaTable,
        axisA=axisA,
        axisB=axisB,
        axisAbin=axisAbin,
        axisBbin=axisBbin,
    )

    log.debug("completed the ``functionName`` function")
    return orderPolyTable, orderPixelTable, orderMetaTable


def _order_table_axes(orderPolyTable, binx, biny):
    """*name the order table's dispersion and cross-dispersion axes, with the binning that applies to each*

    **Key Arguments:**

    - ``orderPolyTable`` -- the order table's polynomial coefficient dataframe
    - ``binx`` -- binning in the x-axis (from FITS header)
    - ``biny`` -- binning in the y-axis (from FITS header)

    **Return:**

    - ``axisA`` -- the axis the order traces are solved for
    - ``axisB`` -- the axis the order traces are sampled along
    - ``axisAbin`` -- the binning of ``axisA``
    - ``axisBbin`` -- the binning of ``axisB``
    """
    if "degy_cent" in orderPolyTable.columns:
        return "x", "y", binx, biny
    return "y", "x", biny, binx


def _order_table_pixel_grid(orderMetaTable, axisB, ratio, extend, pixelDelta):
    """*build the order-pixel dataframe holding one sampled axis-B coordinate per row*

    **Key Arguments:**

    - ``orderMetaTable`` -- the order table's metadata dataframe
    - ``axisB`` -- the axis the order traces are sampled along
    - ``ratio`` -- the factor the metadata limits are scaled by before sampling
    - ``extend`` -- fractional increase to the order area along axis B
    - ``pixelDelta`` -- space between sampled coordinates

    **Return:**

    - ``orderPixelTable`` -- a dataframe of axis-B coordinates and their order numbers
    """
    import math

    import numpy as np
    import pandas as pd

    blower = orderMetaTable[f"{axisB}min"].values * ratio
    bupper = orderMetaTable[f"{axisB}max"].values * ratio
    brange = orderMetaTable[f"{axisB}max"].values * ratio - orderMetaTable[f"{axisB}min"].values * ratio

    axisBcoords = [
        np.arange(
            (0 if (math.floor(lower) - int(r * extend)) < 0 else (math.floor(lower) - int(r * extend))),
            (4200 if (math.ceil(u) + int(r * extend)) > 4200 else (math.ceil(u) + int(r * extend))),
            pixelDelta,
        )
        for lower, u, r in zip(blower, bupper, brange, strict=False)
    ]

    orders = [np.full_like(a, o) for a, o in zip(axisBcoords, orderMetaTable["order"].values, strict=False)]

    # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
    myDict = {
        f"{axisB}coord": np.concatenate(axisBcoords),
        "order": np.concatenate(orders),
    }
    return pd.DataFrame(myDict)


def _order_polynomial_coefficients(orderPolyTable, coefficientPrefix):
    """*collect one order-table polynomial's coefficients, in column order*

    **Key Arguments:**

    - ``orderPolyTable`` -- the order table's polynomial coefficient dataframe
    - ``coefficientPrefix`` -- the column-name prefix of the wanted coefficients, such as ``cent_``

    **Return:**

    - ``coefficients`` -- the matching coefficients as floats, empty when the table carries none
    """
    return [float(v) for k, v in orderPolyTable.iloc[0].items() if coefficientPrefix in k]


def _evaluate_order_polynomial(log, orderPolyTable, orderPixelTable, axisB, degreeSuffix, coefficients):
    """*evaluate one of the order table's Chebyshev polynomials over the order-pixel grid*

    **Key Arguments:**

    - ``log`` -- logger
    - ``orderPolyTable`` -- the order table's polynomial coefficient dataframe
    - ``orderPixelTable`` -- the order-pixel dataframe to evaluate the polynomial over
    - ``axisB`` -- the axis the order traces are sampled along
    - ``degreeSuffix`` -- the column-name suffix giving the polynomial degrees, such as ``cent``
    - ``coefficients`` -- the polynomial coefficients

    **Return:**

    - ``values`` -- the evaluated axis-A coordinate for each row of ``orderPixelTable``
    """
    poly = chebyshev_order_xy_polynomials(
        log=log,
        axisBCol=f"{axisB}coord",
        orderCol="order",
        orderDeg=int(orderPolyTable.iloc[0][f"degorder_{degreeSuffix}"]),
        axisBDeg=int(orderPolyTable.iloc[0][f"deg{axisB}_{degreeSuffix}"]),
    ).poly
    return poly(orderPixelTable, *coefficients)


def _rescale_order_tables_for_binning(orderPixelTable, orderMetaTable, axisA, axisB, axisAbin, axisBbin):
    """*scale the unpacked order tables from unbinned pixels to the frame's binning*

    **Key Arguments:**

    - ``orderPixelTable`` -- the order-pixel dataframe
    - ``orderMetaTable`` -- the order table's metadata dataframe
    - ``axisA`` -- the axis the order traces are solved for
    - ``axisB`` -- the axis the order traces are sampled along
    - ``axisAbin`` -- the binning of ``axisA``
    - ``axisBbin`` -- the binning of ``axisB``

    **Return:**

    - ``orderPixelTable`` -- the order-pixel dataframe in binned pixels
    - ``orderMetaTable`` -- the metadata dataframe in binned pixels
    """
    if axisAbin != 1:
        for c in ["coord_centre", "coord_edgeup", "coord_edgelow"]:
            if f"{axisA}{c}" in orderPixelTable.columns:
                orderPixelTable[f"{axisA}{c}"] /= axisAbin
    if axisBbin != 1:
        orderMetaTable[f"{axisB}min"] /= axisBbin
        orderMetaTable[f"{axisB}max"] /= axisBbin
        orderPixelTable[f"{axisB}coord"] /= axisBbin
        orderPixelTable["std"] /= axisBbin
        mask = orderPixelTable[f"{axisB}coord"].mod(1) > 0
        orderPixelTable = orderPixelTable.loc[~mask]
        orderPixelTable[f"{axisB}coord"] = orderPixelTable[f"{axisB}coord"].round().astype("int")

    return orderPixelTable, orderMetaTable


def generic_quality_checks(log, frame, settings, recipeName, qcTable):
    """*measure very basic quality checks on a frame and return the QC table with results appended*

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- CCDData object
    - ``settings`` -- soxspipe settings
    - ``recipeName`` -- the name of the recipe
    - ``qcTable`` -- the QC pandas data-frame to save the QC measurements

    **Return:**

    - ``qcTable`` -- the QC table with the new measurements appended

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import generic_quality_checks
    qcTable = generic_quality_checks(log=log, frame=myFrame, settings=settings, recipeName="my recipe", qcTable=qcTable)
    ```
    """
    log.debug("starting the ``functionName`` function")

    import numpy as np
    import pandas as pd

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(log=log, settings=settings).get
    kw = kw
    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING KEYWORD; DELETING IT REMOVES THAT FAILURE
    arm = frame.header[kw("SEQ_ARM")]  # noqa: F841
    dateObs = frame.header[kw("DATE_OBS")]

    # nanCount = np.count_nonzero(np.isnan(frame.data))

    utcnow = utcnow_string()

    # COUNT BAD-PIXELS
    badCount = frame.mask.sum()
    totalPixels = np.size(frame.mask)
    percent = float(badCount) / float(totalPixels)
    percent = float(f"{percent:.6f}")

    if "dark" in recipeName.lower():
        qcComment = "Number of hot pixels"
        qcName = "HOTPIX NUM"
    elif "flat" in recipeName.lower():
        qcComment = "Number of cold pixels"
        qcName = "COLDPIX NUM"
    else:
        qcComment = "Number of bad pixels"
        qcName = "BADPIX NUM"

    qcTable = pd.concat(
        [
            qcTable,
            pd.DataFrame([
                {
                    "soxspipe_recipe": recipeName,
                    "qc_name": qcName,
                    "qc_value": int(badCount),
                    "qc_comment": qcComment,
                    "qc_unit": "",
                    "obs_date_utc": dateObs,
                    "reduction_date_utc": utcnow,
                    "to_header": True,
                }
            ]),
        ],
        ignore_index=True,
    )

    if "dark" in recipeName.lower():
        qcComment = "Fraction of hot pixels"
        qcName = "HOTPIX FRAC"
    elif "flat" in recipeName.lower():
        qcComment = "Fraction of cold pixels"
        qcName = "COLDPIX FRAC"
    else:
        qcComment = "Fraction of bad pixels"
        qcName = "BADPIX FRAC"

    qcTable = pd.concat(
        [
            qcTable,
            pd.DataFrame([
                {
                    "soxspipe_recipe": recipeName,
                    "qc_name": qcName,
                    "qc_value": percent,
                    "qc_comment": qcComment,
                    "qc_unit": "",
                    "obs_date_utc": dateObs,
                    "reduction_date_utc": utcnow,
                    "to_header": True,
                }
            ]),
        ],
        ignore_index=True,
    )

    log.debug("completed the ``functionName`` function")
    return qcTable


def spectroscopic_image_quality_checks(log, frame, orderTablePath, settings, recipeName, qcTable):
    """*measure and record spectroscopic image quailty checks*

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- CCDData object
    - ``orderTablePath`` -- path to the order table
    - ``settings`` -- soxspipe settings
    - ``recipeName`` -- the name of the recipe
    - ``qcTable`` -- the QC pandas data-frame to save the QC measurements

    **Return:**

    - ``qcTable`` -- the QC table with the new measurements appended

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import spectroscopic_image_quality_checks
    qcTable = spectroscopic_image_quality_checks(
        log=log,
        frame=myFrame,
        orderTablePath=orderTablePath,
        settings=settings,
        recipeName="this recipe",
        qcTable=qcTable,
    )
    ```
    """
    log.debug("starting the ``functionName`` function")

    import numpy as np
    import numpy.ma as ma
    import pandas as pd

    from soxspipe.commonutils import detector_lookup

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(log=log, settings=settings).get

    arm = frame.header[kw("SEQ_ARM")]
    dateObs = frame.header[kw("DATE_OBS")]

    try:
        binx = frame.header[kw("WIN_BINX")]
        biny = frame.header[kw("WIN_BINY")]
    except KeyError as e:
        log.debug(f"spectroscopic_image_quality_checks: `binx = frame.header[kw('WIN_BINX')]` failed, continuing: {e}")
        if arm.lower() == "nir":
            binx = 1
            biny = 1

    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING KEYWORD; DELETING IT REMOVES THAT FAILURE
    inst = frame.header[kw("INSTRUME")]  # noqa: F841

    # DETECTOR PARAMETERS LOOKUP OBJECT
    detectorParams = detector_lookup(log=log, settings=settings).get(arm)

    if detectorParams["dispersion-axis"] == "x":
        axisA = "x"
        axisB = "y"
    else:
        axisA = "y"
        axisB = "x"

    # UNPACK THE ORDER TABLE
    orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
        log=log, orderTablePath=orderTablePath, binx=binx, biny=biny, prebinned=True
    )

    mask = _inner_order_mask(frame, orderTablePixels, axisA, axisB)

    # COMBINE MASK WITH THE BAD PIXEL MASK
    mask = (mask == 1) | (frame.mask == 1)

    # PLOT ONE OF THE MASKED FRAMES TO CHECK
    maskedFrame = ma.array(frame.data, mask=mask)
    quicklook_image(log=log, CCDObject=maskedFrame, show=False, ext=None)

    mean = np.ma.mean(maskedFrame)
    flux = np.ma.sum(maskedFrame)

    utcnow = utcnow_string()

    # A FULLY MASKED FRAME GIVES np.ma.masked, WHICH `%0.*f` FORMATS AS "nan" BUT `:.3f` FORMATS AS "--"
    mean = "%0.*f" % (3, mean)  # noqa: UP031
    flux = "%0.*f" % (3, flux)  # noqa: UP031

    qcTable = pd.concat(
        [
            qcTable,
            pd.DataFrame([
                {
                    "soxspipe_recipe": recipeName,
                    "qc_name": "INNER ORDER PIX MEAN",
                    "qc_value": mean,
                    "qc_comment": "[e-] Mean inner-order pixel value",
                    "qc_unit": "electrons",
                    "obs_date_utc": dateObs,
                    "reduction_date_utc": utcnow,
                    "to_header": True,
                }
            ]),
        ],
        ignore_index=True,
    )

    qcTable = pd.concat(
        [
            qcTable,
            pd.DataFrame([
                {
                    "soxspipe_recipe": recipeName,
                    "qc_name": "INNER ORDER PIX SUM",
                    "qc_value": flux,
                    "qc_comment": "[e-] Sum of all inner-order pixel values",
                    "qc_unit": "electrons",
                    "obs_date_utc": dateObs,
                    "reduction_date_utc": utcnow,
                    "to_header": True,
                }
            ]),
        ],
        ignore_index=True,
    )

    log.debug("completed the ``functionName`` function")
    return qcTable


def _inner_order_mask(frame, orderTablePixels, axisA, axisB):
    """*build a mask that keeps only the pixels lying between the order edges*

    **Key Arguments:**

    - ``frame`` -- CCDData object
    - ``orderTablePixels`` -- the unpacked order-pixel dataframe carrying the order edge coordinates
    - ``axisA`` -- the axis the order edges are given on
    - ``axisB`` -- the axis the order edges are sampled along

    **Return:**

    - ``mask`` -- an array of ones, zeroed at every inner-order pixel
    """
    import numpy as np

    mask = np.ones_like(frame.data)

    axisACoords_up = orderTablePixels[f"{axisA}coord_edgeup"].values
    axisACoords_low = orderTablePixels[f"{axisA}coord_edgelow"].values
    axisBCoords = orderTablePixels[f"{axisB}coord"].values
    axisACoords_up = axisACoords_up.astype(int)
    axisACoords_low = axisACoords_low.astype(int)

    if axisA == "x":
        for u, lower, y in zip(axisACoords_up, axisACoords_low, axisBCoords, strict=False):
            y = int(y)
            lower = int(max(0, lower))
            u = int(min(mask.shape[1], u))
            if 0 <= y < mask.shape[0] and lower < u:
                mask[y, lower:u] = 0
    else:
        for u, lower, x in zip(axisACoords_up, axisACoords_low, axisBCoords, strict=False):
            x = int(x)
            lower = int(max(0, lower))
            u = int(min(mask.shape[0], u))
            if 0 <= x < mask.shape[1] and lower < u:
                mask[lower:u, x] = 0

    return mask


def read_spectral_format(log, settings, arm, dispersionMap=False, extended=True, binx=1, biny=1):
    """*read the spectral format table to get some key parameters*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- soxspipe settings
    - ``arm`` -- arm to retrieve format for
    - ``dispersionMap`` -- if a dispersion map is given, the minimum and maximum dispersion axis pixel limits are
      computed
    - ``extended`` -- the spectral format table can provide WLMIN/WLMAX (extended=False) or WLMINFUL/WLMAXFUL
      (extended=True)
    - ``binx`` -- binning in the x-axis (from FITS header). Default *1*
    - ``biny`` -- binning in the y-axis (from FITS header). Default *1*

    **Return:**

    - ``orderNums`` -- an array of the order numbers
    - ``waveLengthMin`` -- an array of the minimum wavelengths reached by each order
    - ``waveLengthMax`` -- an array of the maximum wavelengths reached by each order
    - ``amins`` -- the minimum dispersion-axis pixel limit of each order (only when ``dispersionMap`` is given)
    - ``amaxs`` -- the maximum dispersion-axis pixel limit of each order (only when ``dispersionMap`` is given)

    Three values are returned without a ``dispersionMap``, and five with one.

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import read_spectral_format
    # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
    orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=arm)
    ```
    """
    log.debug("starting the ``read_spectral_format`` function")

    # DETECTOR PARAMETERS LOOKUP OBJECT
    dp = detector_lookup(log=log, settings=settings).get(arm)

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    # UNUSED, BUT BUILDING THE LOOKUP READS AND VALIDATES THE KEYWORD MAP; DELETING IT REMOVES THAT FAILURE
    kw = keyword_lookup(log=log, settings=settings).get  # noqa: F841

    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING PARAMETER; DELETING IT REMOVES THAT FAILURE
    science_pixels = dp["science-pixels"]  # noqa: F841

    # READ THE SPECTRAL FORMAT TABLE FILE
    calibrationRootPath = get_calibrations_path(log=log, settings=settings)
    spectralFormatFile = calibrationRootPath + "/" + dp["spectral format table"]

    # SPEC FORMAT TO PANDAS DATAFRAME
    from astropy.table import Table

    dat = Table.read(spectralFormatFile, format="fits")
    specFormatTable = dat.to_pandas()

    # EXTRACT REQUIRED PARAMETERS
    orderNums = specFormatTable["ORDER"].values

    if extended or "WLMIN" not in specFormatTable.columns:
        waveLengthMin = specFormatTable["WLMINFUL"].values
        waveLengthMax = specFormatTable["WLMAXFUL"].values
    else:
        waveLengthMin = specFormatTable["WLMIN"].values
        waveLengthMax = specFormatTable["WLMAX"].values

    # USE DISPERSION MAP TO FIND X-Y LIMITS OF THE SPECTRAL FORMAT FOR EACH ORDER
    # WE WANT TO LIMIT THE EXTRACTION TO THESE REGIONS
    if not isinstance(dispersionMap, bool):
        amins, amaxs = _dispersion_axis_pixel_limits(
            log=log,
            dp=dp,
            dispersionMap=dispersionMap,
            orderNums=orderNums,
            waveLengthMin=waveLengthMin,
            waveLengthMax=waveLengthMax,
            binx=binx,
            biny=biny,
        )

        log.debug("completed the ``read_spectral_format`` function")
        return orderNums, waveLengthMin, waveLengthMax, amins, amaxs

    log.debug("completed the ``read_spectral_format`` function")
    return orderNums, waveLengthMin, waveLengthMax


def _dispersion_axis_pixel_limits(log, dp, dispersionMap, orderNums, waveLengthMin, waveLengthMax, binx, biny):
    """*find each order's dispersion-axis pixel limits, clipped to the detector and scaled to the binning*

    **Key Arguments:**

    - ``log`` -- logger
    - ``dp`` -- the detector parameters dictionary
    - ``dispersionMap`` -- path to the dispersion map
    - ``orderNums`` -- an array of the order numbers
    - ``waveLengthMin`` -- an array of the minimum wavelengths reached by each order
    - ``waveLengthMax`` -- an array of the maximum wavelengths reached by each order
    - ``binx`` -- binning in the x-axis (from FITS header)
    - ``biny`` -- binning in the y-axis (from FITS header)

    **Return:**

    - ``amins`` -- the minimum dispersion-axis pixel limit of each order
    - ``amaxs`` -- the maximum dispersion-axis pixel limit of each order
    """
    import numpy as np
    import pandas as pd

    myDict = {
        "order": np.asarray([]),
        "wavelength": np.asarray([]),
        "slit_position": np.asarray([]),
    }
    for o, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax, strict=False):
        wlArray = np.array([wmin, wmax])
        myDict["wavelength"] = np.append(myDict["wavelength"], wlArray)
        myDict["order"] = np.append(myDict["order"], np.ones(len(wlArray)) * o)
        myDict["slit_position"] = np.append(myDict["slit_position"], np.zeros(len(wlArray)))
    orderPixelTable = pd.DataFrame(myDict)
    orderPixelTable = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath=dispersionMap,
        orderPixelTable=orderPixelTable,
        removeOffDetectorLocation=False,
    )

    if dp["dispersion-axis"] == "x":
        axis = "y"
        rowCol = "rows"
        abinFactor = biny
    else:
        axis = "x"
        rowCol = "columns"
        abinFactor = binx

    amins = []
    amaxs = []
    for o in orderNums:
        amin = orderPixelTable.loc[orderPixelTable["order"] == o, f"fit_{axis}"].min()
        amax = orderPixelTable.loc[orderPixelTable["order"] == o, f"fit_{axis}"].max()
        if amin < 0:
            amin = 0
        if amax > dp["science-pixels"][rowCol]["end"]:
            amax = dp["science-pixels"][rowCol]["end"]
        amins.append(amin / abinFactor)
        amaxs.append(amax / abinFactor)

    return amins, amaxs


def get_calibrations_path(log, settings):
    """*return the root path to the static calibrations*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary

    **Return:**

    - ``calibrationRootPath`` -- the root path to the instrument's static calibrations

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import get_calibrations_path
    calibrationRootPath = get_calibrations_path(log=log, settings=settings)
    ```
    """
    log.debug("starting the ``get_calibrations_path`` function")

    # GENERATE PATH TO STATIC CALIBRATION DATA
    # KEPT AS `in` PLUS INDEXING: FOR A NON-DICT `settings` (E.G. THE `False` DEFAULT ELSEWHERE IN THIS
    # MODULE), `.get` WOULD RAISE AttributeError WHERE THIS RAISES TypeError
    instrument = settings["instrument"] if "instrument" in settings else "soxs"  # noqa: SIM401
    calibrationRootPath = os.path.dirname(os.path.dirname(__file__)) + "/resources/static_calibrations/" + instrument

    log.debug("completed the ``get_calibrations_path`` function")
    return calibrationRootPath


# THE NAME IS PUBLIC AND IMPORTED ACROSS THE PACKAGE, SO RENAMING IT WOULD CHANGE THE PUBLIC SURFACE.
def twoD_disp_map_image_to_dataframe(  # noqa: N802
    log,
    slit_length,
    twoDMapPath,
    kw=False,
    associatedFrame=False,
    removeMaskedPixels=False,
    dispAxis="y",
):
    """*convert the 2D dispersion image map to a pandas dataframe*

    **Key Arguments:**

    - ``log`` -- logger
    - ``slit_length`` -- length of the slit; pixels with a slit position beyond half this length either side of
      the centre are removed
    - ``twoDMapPath`` -- 2D dispersion map image path
    - ``kw`` -- fits keyword lookup dictionary
    - ``associatedFrame`` -- include a flux column in returned dataframe from a frame associated with the
      dispersion map. Default *False*
    - ``removeMaskedPixels`` -- remove the masked pixels from the associated image? Default *False*
    - ``dispAxis`` -- x or y. Needed for pixel scale calculation

    **Return:**

    - ``mapDF`` -- the dispersion map as a dataframe, one row per pixel
    - ``interOrderMask`` -- mask array flagging the inter-order pixels

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
    mapDF, interOrderMask = twoD_disp_map_image_to_dataframe(
        log=log, slit_length=11, twoDMapPath=twoDMap, associatedFrame=objectFrame, kw=kw
    )
    ```
    """
    log.debug("starting the ``twoD_disp_map_image_to_dataframe`` function")

    # MAKE RELATIVE HOME PATH ABSOLUTE
    from os.path import expanduser

    import numpy as np
    import pandas as pd

    home = expanduser("~")
    if twoDMapPath[0] == "~":
        twoDMapPath = twoDMapPath.replace("~", home)

    binx, biny = _associated_frame_binning(associatedFrame, kw)

    hdul, minimumBinnedPixelValue, binned = _open_disp_map_hdus(twoDMapPath, binx, biny)

    mapDF = pd.DataFrame.from_dict(
        _disp_map_pixel_columns(hdul, minimumBinnedPixelValue, binned, associatedFrame)
    )
    if removeMaskedPixels:
        mask = mapDF["mask"].eq(False)
        mapDF = mapDF.loc[mask]

    # REMOVE ZEROS
    mask = mapDF["wavelength"] == 0
    mapDF = mapDF.loc[~mask]

    interOrderMask = hdul["ORDER"].data.copy()
    interOrderMask = np.where(interOrderMask > 0, 0, interOrderMask)
    interOrderMask = np.where(np.isnan(interOrderMask), 1, interOrderMask)

    mapDF.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

    # REMOVE FILTERED ROWS FROM DATA FRAME
    mask = (mapDF["slit_position"] < -slit_length / 2) | (mapDF["slit_position"] > slit_length / 2)
    mapDF = mapDF.loc[~mask]
    mask = mapDF["min"] == 0
    mapDF = mapDF.loc[~mask]

    # SORT BY COLUMN NAME
    if not mapDF.empty:
        _add_pixel_scale_column(mapDF, dispAxis)

    log.debug("completed the ``twoD_disp_map_image_to_dataframe`` function")
    return mapDF, interOrderMask


def _associated_frame_binning(associatedFrame, kw):
    """*read the x and y binning of the frame associated with a 2D dispersion map*

    **Key Arguments:**

    - ``associatedFrame`` -- the frame associated with the dispersion map, or *False* for none
    - ``kw`` -- fits keyword lookup dictionary

    **Return:**

    - ``binx`` -- binning in the x-axis, *1* when the frame does not report one
    - ``biny`` -- binning in the y-axis, *1* when the frame does not report one
    """
    # FIND THE APPROPRIATE PREDICTED LINE-LIST
    if associatedFrame:
        arm = associatedFrame.header[kw("SEQ_ARM")]
        if arm != "NIR" and kw("WIN_BINX") in associatedFrame.header:
            return int(associatedFrame.header[kw("WIN_BINX")]), int(associatedFrame.header[kw("WIN_BINY")])
    return 1, 1


def _open_disp_map_hdus(twoDMapPath, binx, biny):
    """*open a 2D dispersion map and block-reduce its extensions to the frame's binning*

    **Key Arguments:**

    - ``twoDMapPath`` -- 2D dispersion map image path
    - ``binx`` -- binning in the x-axis
    - ``biny`` -- binning in the y-axis

    **Return:**

    - ``hdul`` -- the open dispersion-map HDU list, block-reduced when the frame is binned
    - ``minimumBinnedPixelValue`` -- the minimum unbinned wavelength within each binned pixel, or *None* when
      the frame is unbinned
    - ``binned`` -- was the map block-reduced?
    """
    import numpy as np
    from astropy.io import fits

    hdul = fits.open(twoDMapPath)

    hdul["WAVELENGTH"].data = hdul["WAVELENGTH"].data.astype("float32")
    hdul["SLIT"].data = hdul["SLIT"].data.astype("float32")
    hdul["ORDER"].data = hdul["ORDER"].data.astype("float32")

    binned = False
    minimumBinnedPixelValue = None
    if binx > 1 or biny > 1:
        binned = True
        from astropy.nddata import block_reduce

        minimumBinnedPixelValue = hdul["WAVELENGTH"].data.copy()
        hdul["WAVELENGTH"].data = block_reduce(hdul["WAVELENGTH"].data, (biny, binx), func=np.mean)
        hdul["SLIT"].data = block_reduce(hdul["SLIT"].data, (biny, binx), func=np.mean)
        hdul["ORDER"].data = block_reduce(hdul["ORDER"].data, (biny, binx), func=np.mean)
        minimumBinnedPixelValue = block_reduce(minimumBinnedPixelValue, (biny, binx), func=np.min)
        minimumBinnedPixelValue = minimumBinnedPixelValue.flatten()

    return hdul, minimumBinnedPixelValue, binned


def _disp_map_pixel_columns(hdul, minimumBinnedPixelValue, binned, associatedFrame):
    """*build the one-row-per-pixel columns of the dispersion-map dataframe*

    **Key Arguments:**

    - ``hdul`` -- the open dispersion-map HDU list
    - ``minimumBinnedPixelValue`` -- the minimum unbinned wavelength within each binned pixel
    - ``binned`` -- was the map block-reduced?
    - ``associatedFrame`` -- the frame associated with the dispersion map, or *False* for none

    **Return:**

    - ``thisDict`` -- a dictionary of equal-length pixel columns
    """
    import numpy as np

    # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
    xdim = hdul[0].data.shape[1]
    ydim = hdul[0].data.shape[0]
    xarray = np.tile(np.arange(0, xdim), ydim)
    yarray = np.repeat(np.arange(0, ydim), xdim)

    if not binned:
        minimumBinnedPixelValue = np.ones_like(yarray)

    thisDict = {
        "x": xarray,
        "y": yarray,
        "wavelength": hdul["WAVELENGTH"].data.flatten().astype("float32"),
        "slit_position": hdul["SLIT"].data.flatten().astype("float32"),
        "order": hdul["ORDER"].data.flatten().astype("float32"),
        "min": minimumBinnedPixelValue.astype("float32"),
    }

    if associatedFrame:
        thisDict["flux"] = associatedFrame.data.flatten().astype("float32")
        thisDict["mask"] = associatedFrame.mask.flatten().astype(bool)
        thisDict["error"] = associatedFrame.uncertainty.array.flatten().astype("float32")

    # REMOVE IF ABOVE .astype(float) IS WORKING
    # try:
    #     if associatedFrame:
    #         thisDict["flux"] = associatedFrame.data.flatten()
    #         thisDict["mask"] = associatedFrame.mask.flatten()
    #         thisDict["error"] = associatedFrame.uncertainty.array.flatten()

    # except Exception as e:

    #     if binned:
    #         minimumBinnedPixelValue = minimumBinnedPixelValue.byteswap().newbyteorder()

    #     thisDict = {
    #         "x": xarray,
    #         "y": yarray,
    #         "wavelength": hdul["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
    #         "slit_position": hdul["SLIT"].data.flatten().byteswap().newbyteorder(),
    #         "order": hdul["ORDER"].data.flatten().byteswap().newbyteorder(),
    #         "min": minimumBinnedPixelValue
    #     }
    #     if associatedFrame:
    #         thisDict["flux"] = associatedFrame.data.flatten().byteswap().newbyteorder()
    #         thisDict["mask"] = associatedFrame.mask.flatten().byteswap().newbyteorder()
    #         thisDict["error"] = associatedFrame.uncertainty.array.flatten().byteswap().newbyteorder()

    return thisDict


def _add_pixel_scale_column(mapDF, dispAxis):
    """*add the per-pixel wavelength step to the dispersion-map dataframe, in place*

    The dataframe is left sorted by wavelength.

    **Key Arguments:**

    - ``mapDF`` -- the dispersion-map dataframe, modified in place
    - ``dispAxis`` -- x or y. Sets the order the pixels are stepped through in
    """
    import numpy as np

    mapDF.sort_values(["wavelength"], inplace=True, kind="stable")

    # CALCULATE PIXEL SCALE
    if dispAxis == "y":
        mapDF.sort_values(["x", "y"], inplace=True, kind="stable")
    else:
        mapDF.sort_values(["y", "x"], inplace=True, kind="stable")
    shiftedWlArray = list(mapDF["wavelength"].values)[1:]
    shiftedWlArray.append(np.nan)
    mapDF["pixelScale"] = mapDF["wavelength"] - shiftedWlArray
    mask = (mapDF["pixelScale"] > 2) | (mapDF["pixelScale"] < -2)
    mapDF.loc[mask, "pixelScale"] = 0.0
    mapDF["pixelScale"] = mapDF["pixelScale"].abs()

    # SORT BY COLUMN NAME
    mapDF.sort_values(["wavelength"], inplace=True, kind="stable")


def predict_product_path(sofName, recipeName=False):
    """*predict the path of the recipe product from a given SOF name*

    **Key Arguments:**

    - ``sofName`` -- name or full path to the sof file
    - ``recipeName`` -- name of the recipe being considered. Default *False*.

    **Return:**

    - ``productPath`` -- the predicted path of the recipe product
    - ``startNightDate`` -- the start-of-night date of the observations

    **Usage:**

    ```python
    from soxspipe.commonutils import toolkit
    productPath, startNightDate = toolkit.predict_product_path(sofFilePath)
    ```
    """

    from astropy.time import Time, TimeDelta

    startNightDate = False

    # TRY AND READ startNightDate FROM RAW FRAME DIRECTORY
    # with codecs.open(sofName, encoding="utf-8", mode="r") as readFile:
    #     thisData = readFile.read()
    #     for l in thisData.split("\n"):
    #         if "raw/" in l:
    #             startNightDate = l.split("raw/")[1].split("/")[0]
    #             break

    try:
        sofName = os.path.basename(sofName)
    except (TypeError, AttributeError) as e:
        # NO LOGGER EXISTS YET AT THIS POINT IN THE FUNCTION, SO REPORT THE SAME WAY THE OBSDATE HANDLER BELOW DOES
        print(f"predict_product_path: no sof filename in {sofName!r}, the product path will not resolve: {e}")

    if not recipeName:
        recipeName = sys.argv[1]
        if recipeName[0] == "-":
            recipeName = sys.argv[2]
        recipeName = "soxs-" + recipeName

    sofName = sofName.replace(".sof", "")

    if not startNightDate:

        obsDate = sofName.split("_")[0]

        startNightDate = ""

        try:
            obsDate = Time.strptime(obsDate, "%Y%m%dT%H%M%S", scale="utc")
            night_start_offset = TimeDelta(15.0 * 60 * 60, format="sec")
            startNightDate = obsDate - night_start_offset
            startNightDate = startNightDate.strftime("%Y-%m-%d")
        except (ValueError, TypeError):
            print("Could not determine OBSDATE from sof filename")
            pass

    from fundamentals.logs import emptyLogger

    from soxspipe.commonutils import data_organiser

    log = emptyLogger()
    do = data_organiser(log=log, rootDir=".", dbConnect=False)
    currentSession, allSessions = do.session_list(silent=True)
    do.close()

    if "_STARE_STD" in sofName or "_NOD_STD" in sofName:
        sofName += "_RESP"
    elif "_STARE" in sofName or "_NOD" in sofName:
        sofName += "_EXTRACTED_MERGED"
    productPath = (
        f"./sessions/{currentSession}/reduced/{startNightDate}/"
        + recipeName.replace("_", "-")
        + "/"
        + sofName
        + ".fits"
    )
    productPath = productPath.replace("//", "/")

    return productPath, startNightDate


def add_recipe_logger(log, productPath):
    """*add a recipe-specific handler to the default logger that writes the recipe's logs beside the product*

    **Key Arguments:**

    - ``log`` -- original logger
    - ``productPath`` -- path to the recipe product

    **Return:**

    - ``log`` -- the logger with the recipe-specific handlers attached

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import add_recipe_logger
    log = add_recipe_logger(log, productPath="/path/to/product")
    ```

    """
    import logging
    import os

    i = 0
    while i < 3:
        for handler in log.handlers:
            if handler.get_name() == "recipeLog":
                log.removeHandler(handler)
            if handler.get_name() == "recipeErr":
                log.removeHandler(handler)
        i += 1

    # GET THE EXTENSION (WITH DOT PREFIX)
    loggingPath = os.path.splitext(productPath)[0] + ".log"
    loggingErrorPath = os.path.splitext(productPath)[0] + "_ERROR.log"
    try:
        os.remove(loggingPath)
        os.remove(loggingErrorPath)
    except OSError as e:
        log.debug(f"add_recipe_logger: `os.remove(loggingPath)` failed, continuing: {e}")

    # PARENT DIRECTORY PATH NEEDS TO EXIST FOR LOGGER TO WRITE
    parentDirectory = os.path.dirname(loggingPath)
    if not os.path.exists(parentDirectory):
        try:
            os.makedirs(parentDirectory)
        except OSError as e:
            log.debug(f"add_recipe_logger: `os.makedirs(parentDirectory)` failed, continuing: {e}")

    recipeLog = logging.FileHandler(loggingPath, mode="a", encoding=None, delay=False)
    recipeLogFormatter = logging.Formatter("%(message)s")
    recipeLog.set_name("recipeLog")
    recipeLog.setLevel(logging.INFO + 1)
    recipeLog.setFormatter(recipeLogFormatter)
    recipeLog.addFilter(MaxFilter(logging.WARNING))
    log.addHandler(recipeLog)

    recipeErr = logging.FileHandler(loggingErrorPath, mode="a", encoding=None, delay=True)
    recipeErrFormatter = logging.Formatter(
        '%(asctime)s %(levelname)s: "%(pathname)s", line %(lineno)d, in %(funcName)s > %(message)s',
        "%Y-%m-%d %H:%M:%S",
    )
    recipeErr.set_name("recipeErr")
    recipeErr.setLevel(logging.ERROR)
    recipeErr.setFormatter(recipeErrFormatter)
    log.addHandler(recipeErr)

    return log


class MaxFilter:
    def __init__(self, max_level):
        self.max_level = max_level

    def filter(self, record):
        if record.levelno < self.max_level:
            return True
        return None


def create_dispersion_solution_grid_lines_for_plot(
    log,
    dispMap,
    dispMapImage,
    associatedFrame,
    kw,
    skylines=False,
    slitPositions=False,
    slit_length=11,
):
    """*generate the grid lines to add to QC plots from a dispersion solution and its 2D dispersion map image*

    **Key Arguments:**

    - ``log`` -- logger
    - ``dispMap`` -- path to dispersion map. Default *False*
    - ``dispMapImage`` -- the 2D dispersion map image
    - ``associatedFrame`` -- a frame associated with the reduction (to read arm and binning info).
    - ``kw`` -- fits header kw dictionary
    - ``skylines`` -- a list of skylines to use as the grid. Default *False*
    - ``slitPositions`` -- slit positions to plot (else plot min and max)
    - ``slit_length`` -- length of the slit to use for the dispersion map dataframe (default 11)

    **Return:**

    - ``orderPixelTable`` -- DataFrame containing the pixel coordinates for grid lines to plot.
    - ``interOrderMask`` -- Mask array indicating inter-order regions.

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import create_dispersion_solution_grid_lines_for_plot
    gridLinePixelTable, interOrderMask = create_dispersion_solution_grid_lines_for_plot(
        log=log,
        dispMap=dispMap,
        dispMapImage=dispMapImage,
        associatedFrame=CCDObject,
        kw=kw,
        skylines=skylines
    )

    for l in range(int(gridLinePixelTable['line'].max())):
        mask = (gridLinePixelTable['line'] == l)
        ax.plot(
            gridLinePixelTable.loc[mask]["fit_y"],
            gridLinePixelTable.loc[mask]["fit_x"],
            "w-",
            linewidth=0.5,
            alpha=0.8,
            color="black",
        )
    ```
    """
    log.debug("starting the ``create_dispersion_solution_grid_lines_for_plot`` function")

    import numpy as np
    import pandas as pd

    dispMapDF, interOrderMask = twoD_disp_map_image_to_dataframe(
        log=log,
        slit_length=slit_length,
        twoDMapPath=dispMapImage,
        associatedFrame=associatedFrame,
        kw=kw,
    )

    uniqueOrders = dispMapDF["order"].unique()
    wlLims = []
    sPos = []

    for o in uniqueOrders:
        filDF = dispMapDF.loc[dispMapDF["order"] == o]
        wlRange = filDF["wavelength"].max() - filDF["wavelength"].min()
        wlLims.append((filDF["wavelength"].min(), filDF["wavelength"].max()))
        if isinstance(slitPositions, bool):
            sPos.append((filDF["slit_position"].min(), filDF["slit_position"].max()))
        else:
            sPos.append(slitPositions)

    lineNumber = 0
    orderPixelTable_list = []
    for o, wlLim, spLim in zip(uniqueOrders, wlLims, sPos, strict=False):
        wlRange = np.arange(wlLim[0], wlLim[1], 1)
        wlRange = np.append(wlRange, [wlLim[1]])
        for e in spLim:
            myDict = {
                "line": np.full_like(wlRange, lineNumber),
                "order": np.full_like(wlRange, o),
                "wavelength": wlRange,
                "slit_position": np.full_like(wlRange, e),
            }
            orderPixelTable_list.append(pd.DataFrame(myDict))
            lineNumber += 1

        spRange = np.arange(min(spLim), max(spLim), 1)
        spRange = np.append(spRange, [max(spLim)])
        if not isinstance(skylines, bool):
            mask = skylines["WAVELENGTH"].between(wlLim[0], wlLim[1])
            wlRange = skylines.loc[mask]["WAVELENGTH"].values
        else:
            step = max(1, int((wlLim[1] - wlLim[0]) // 400))
            wlRange = np.arange(wlLim[0], wlLim[1], step)
        wlRange = np.append(wlRange, [wlLim[1]])

        for wl in wlRange:
            myDict = {
                "line": np.full_like(spRange, lineNumber),
                "order": np.full_like(spRange, o),
                "wavelength": np.full_like(spRange, wl),
                "slit_position": spRange,
            }
            orderPixelTable_list.append(pd.DataFrame(myDict))
            lineNumber += 1

    orderPixelTable = pd.concat(orderPixelTable_list, ignore_index=True)

    orderPixelTable = dispersion_map_to_pixel_arrays(
        log=log, dispersionMapPath=dispMap, orderPixelTable=orderPixelTable
    )

    log.debug("completed the ``create_dispersion_solution_grid_lines_for_plot`` function")
    return orderPixelTable, interOrderMask


def get_calibration_lamp(log, frame, kw):
    """*given a frame, determine which calibration lamp is being used*

    **Key Arguments:**

    - ``log`` -- logger
    - ``frame`` -- the frame to determine the calibration lamp for
    - ``kw`` -- the FITS header keyword dictionary

    **Return:**

    - ``lamp`` -- the calibration lamp names found in the header, or *None* if there are none

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import get_calibration_lamp
    lamp = get_calibration_lamp(log=log, frame=frame, kw=kw)
    ```
    """
    log.debug("starting the ``read_calibration_lamp`` function")

    # UNUSED, BUT THE LOOKUP RAISES KeyError FOR A MISSING KEYWORD; DELETING IT REMOVES THAT FAILURE
    inst = frame.header["INSTRUME"]  # noqa: F841
    lamp = None

    for lampKeyword in [
        kw("LAMP1"),
        kw("LAMP2"),
        kw("LAMP3"),
        kw("LAMP4"),
        kw("LAMP5"),
        kw("LAMP6"),
        kw("LAMP7"),
    ]:
        if lampKeyword in frame.header:
            newLamp = frame.header[lampKeyword]
            newLamp = (
                newLamp.replace("UVB_High", "QTH")
                .replace("UVB_Low_", "")
                .replace("NIR_", "")
                .replace("VIS_", "")
                .replace("UVB_", "")
                .replace("_lamp", "")
                .replace("_Lamp", "")
                .replace("Argo", "Ar")
                .replace("Neon", "Ne")
                .replace("Merc", "Hg")
                .replace("Xeno", "Xe")
            )
            if lamp:
                lamp += newLamp
            else:
                lamp = newLamp

    log.debug("completed the ``read_calibration_lamp`` function")
    return lamp


def qc_settings_plot_tables(log, qc, qcAx, settings, settingsAx):
    """*generate QC and settings table to be placed at the bottom of the QC plots*

    **Key Arguments:**

    - ``log`` -- logger
    - ``qc`` -- date frame of collected QCs
    - ``qcAx`` -- the axis to add the QC table to
    - ``settings`` -- settings to report in settings table
    - ``settingsAx`` -- the axis to add the settings table to

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import qc_settings_plot_tables
    qc_settings_plot_tables(log=log,qc=self.qc,qcAx=qcAx, settings=settings,settingsAx=settingsAx)
    ```
    """
    log.debug("starting the ``qc_settings_plot_tables`` function")

    import matplotlib as plt
    import numpy as np
    import pandas as pd

    tables = []
    cols = []

    qcCopy = qc.copy()
    if "qc_order" in qcCopy.columns:
        qcCopy = qcCopy.loc[((qcCopy["qc_order"] == "-1") | (qcCopy["qc_order"].isna()) | (qcCopy["qc_order"] == -1))]
    qcCopy["value"] = qcCopy["qc_value"].astype(str) + " " + qcCopy["qc_unit"]
    qcCopy.loc[qcCopy["value"].isnull(), "value"] = qcCopy.loc[qcCopy["value"].isnull(), "qc_value"]

    if len(qcCopy.index) > 10:
        qcCopy = qcCopy.head(10)

    columns1 = ["value", "qc_comment"]
    colColours = plt.cm.Greys(np.full(len(columns1), 0.1))
    rowColours = plt.cm.Greys(np.full(len(qcCopy.index), 0.1))
    rowLabels = qcCopy["qc_name"].values

    if len(qcCopy[columns1].values):
        qcTable = qcAx.table(
            cellText=qcCopy[columns1].values,
            colLabels=columns1,
            loc="center",
            cellLoc="left",
            rowColours=rowColours,
            colColours=colColours,
            rowLabels=rowLabels,
            rowLoc="right",
            fontsize=14,
        )
        tables.append(qcTable)
        cols.append(columns1)
    # qcAx.set_title(
    #     "QC Table", fontsize=9)

    settingsCopy = {k: v for k, v in settings.items() if k not in ["nir", "vis", "uvb", "qc-acceptable-ranges"]}

    settingsCopy = {"setting": settingsCopy.keys(), "value": settingsCopy.values()}

    settingsDF = pd.DataFrame(settingsCopy)

    columns2 = ["value"]
    colColours = plt.cm.Greys(np.full(len(columns2), 0.1))
    rowColours = plt.cm.Greys(np.full(len(settingsDF.index), 0.1))
    rowLabels = settingsDF["setting"].values
    settingsTable = settingsAx.table(
        cellText=settingsDF[columns2].values,
        colLabels=columns2,
        loc="center",
        cellLoc="left",
        rowColours=rowColours,
        colColours=colColours,
        rowLabels=rowLabels,
        rowLoc="right",
        fontsize=14,
    )
    tables.append(settingsTable)
    cols.append(columns2)
    # settingsAx.set_title(
    #     "Parameters", fontsize=9, loc='left')
    settingsAx.margins(x=0, y=0)

    for t, c in zip(tables, cols, strict=False):
        t.scale(1, 1.5)
        t.auto_set_font_size(False)
        t.set_fontsize(4)
        table_cells = t.properties()["children"]
        for cell in table_cells:
            cell.set_linewidth(0.3)
        t.auto_set_column_width(list(range(len(c))))

    for a in [qcAx, settingsAx]:

        # HIDE AXES
        a.get_xaxis().set_visible(False)
        a.get_yaxis().set_visible(False)
        a.axis("off")

    log.debug("completed the ``qc_settings_plot_tables`` function")
    return


def utility_setup(log, settings, recipeName, startNightDate):
    """*setup some tools needed by most soxspipe utils*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night

    **Return:**

    - ``qcDir`` -- the QC directory (created if missing)
    - ``productDir`` -- the product directory (created if missing)

    **Usage:**

    ```python
    # Example usage with timezone-aware UTC datetime
    from datetime import datetime, UTC
    startNightDate = datetime.now(UTC).strftime("%Y-%m-%d")
    qcDir, productDir = utility_setup(log=log, settings=settings, recipeName="my_recipe", startNightDate=startNightDate)
    ```
    """
    log.debug("starting the ``utility_setup`` function")

    from os.path import expanduser

    home = expanduser("~")

    recipeName = recipeName.replace("-obj", "")
    # QC DIR
    qcDir = settings["workspace-root-dir"].replace("~", home) + f"/qc/{startNightDate}/{recipeName}/"
    qcDir = qcDir.replace("//", "/")
    # RECURSIVELY CREATE MISSING DIRECTORIES. A CONCURRENT RUN MAY WIN THE RACE,
    # SO exist_ok ABSORBS THAT; ANY OTHER FAILURE PROPAGATES.
    os.makedirs(qcDir, exist_ok=True)

    # PRODUCT DIR
    productDir = settings["workspace-root-dir"].replace("~", home) + f"/reduced/{startNightDate}/{recipeName}/"
    productDir = productDir.replace("//", "/")
    # RECURSIVELY CREATE MISSING DIRECTORIES. A CONCURRENT RUN MAY WIN THE RACE,
    # SO exist_ok ABSORBS THAT; ANY OTHER FAILURE PROPAGATES.
    os.makedirs(productDir, exist_ok=True)

    log.debug("completed the ``utility_setup`` function")
    return qcDir, productDir


def plot_merged_spectrum_qc(
    merged_orders,
    products,
    log,
    qcDir,
    filenameTemplate,
    noddingSequence,
    dateObs,
    arm,
    recipeName,
    orderJoins=False,
    debug=False,
    fluxCalibrated=False,
    qcTable=False,
    settings=False,
):
    """*plot the order-merged spectrum QC plot, save it to the QC directory and record it in the products table*

    **Key Arguments:**

    - ``merged_orders`` -- the order-merged spectrum, with ``WAVE`` and ``FLUX_COUNTS`` columns
    - ``products`` -- the products table. Nothing is plotted if this is *False*
    - ``log`` -- logger
    - ``qcDir`` -- the directory to save the QC plot in
    - ``filenameTemplate`` -- the product filename the QC plot filename is built from
    - ``noddingSequence`` -- suffix for the QC plot filename and product label. *False* for none
    - ``dateObs`` -- the observation date recorded in the products table
    - ``arm`` -- the spectrograph arm
    - ``recipeName`` -- the name of the recipe
    - ``orderJoins`` -- a dictionary of order-join wavelengths to mark on the plot. Default *False*
    - ``debug`` -- show the plot before saving it. Default *False*
    - ``fluxCalibrated`` -- is the spectrum flux calibrated? Default *False*
    - ``qcTable`` -- the QC table holding the SNR values to plot. Default *False*
    - ``settings`` -- the soxspipe settings dictionary, used to look up the skylines. Default *False*

    **Return:**

    - ``products`` -- the products table with the QC plot appended
    - ``filePath`` -- the path to the saved QC plot PDF, or *None* if nothing was plotted

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc
    products, filePath = plot_merged_spectrum_qc(
        merged_orders=mergedSpectrum,
        products=products,
        log=log,
        qcDir=qcDir,
        filenameTemplate=filenameTemplate,
        noddingSequence=False,
        dateObs=dateObs,
        arm=arm,
        recipeName=recipeName,
        qcTable=qcTable,
        settings=settings,
    )
    ```
    """
    log.debug("starting the ``plot_merged_spectrum_qc`` function")

    # DETECTOR PARAMETERS LOOKUP OBJECT

    # DO NOT PLOT IF PRODUCT TABLE HAS NOT BEEN PASSED
    if isinstance(products, bool) and not products:
        return products, None

    import matplotlib.pyplot as plt
    import pandas as pd

    if not noddingSequence:
        noddingSequence = ""

    skylinesDF = get_skylines_dataframe(log, settings, arm)

    fig = plt.figure(figsize=(14, 10), constrained_layout=True, dpi=180)
    # ADJUSTED HEIGHT RATIOS
    gs = fig.add_gridspec(5, 1, height_ratios=[3, 1, 1, 1, 0])

    top_panel, middle_panel = _draw_merged_flux_panels(
        fig=fig,
        gs=gs,
        merged_orders=merged_orders,
        arm=arm,
        filenameTemplate=filenameTemplate,
        fluxCalibrated=fluxCalibrated,
    )

    bottom_panel = _draw_merged_snr_panel(
        log=log,
        fig=fig,
        gs=gs,
        merged_orders=merged_orders,
        qcTable=qcTable,
    )

    # SKY PANEL WITH SKY FLUX
    sky_panel = _draw_sky_panel(fig, gs, merged_orders, fluxCalibrated)

    if orderJoins:
        _mark_order_joins(orderJoins, [top_panel, middle_panel, bottom_panel])

    # PLOT SKY LINES AS VERTICAL LINES ON SKY PANEL
    _mark_skylines(skylinesDF, [top_panel, middle_panel, bottom_panel, sky_panel], top_panel)

    if fluxCalibrated:
        filename = filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_FLUXCALIBRATED_QC_PLOT{noddingSequence}.pdf")
    else:
        filename = filenameTemplate.replace(".fits", f"_EXTRACTED_MERGED_QC_PLOT{noddingSequence}.pdf")
    filePath = f"{qcDir}/{filename}"
    if debug:
        plt.show()
    save_qc_plot(filePath)
    plt.close("all")

    utcnow = utcnow_string()
    products = pd.concat(
        [
            products,
            pd.DataFrame([
                {
                    "soxspipe_recipe": recipeName,
                    "product_label": (
                        f"EXTRACTED_MERGED_FLUXCALIBRATED_QC_PLOT{noddingSequence}"
                        if fluxCalibrated
                        else f"EXTRACTED_MERGED_QC_PLOT{noddingSequence}"
                    ),
                    "file_name": filename,
                    "file_type": "PDF",
                    "obs_date_utc": dateObs,
                    "reduction_date_utc": utcnow,
                    "product_desc": "QC plot of extracted order-merged source",
                    "file_path": filePath,
                    "label": "QC",
                }
            ]),
        ],
        ignore_index=True,
    )

    log.debug("completed the ``plot_merged_spectrum_qc`` function")
    return products, filePath


def _draw_merged_flux_panels(fig, gs, merged_orders, arm, filenameTemplate, fluxCalibrated):
    """*draw the linear-scale and log-scale flux panels of the merged-spectrum QC plot*

    **Key Arguments:**

    - ``fig`` -- the QC figure
    - ``gs`` -- the figure's gridspec
    - ``merged_orders`` -- the order-merged spectrum
    - ``arm`` -- the spectrograph arm
    - ``filenameTemplate`` -- the product filename the plot title is built from
    - ``fluxCalibrated`` -- is the spectrum flux calibrated?

    **Return:**

    - ``top_panel`` -- the linear-scale flux panel
    - ``middle_panel`` -- the log-scale flux panel
    """
    from astropy.stats import sigma_clip, sigma_clipped_stats

    # TOP PANEL WITH LINEAR SCALE
    top_panel = fig.add_subplot(gs[0, :])
    if fluxCalibrated:
        top_panel.set_ylabel("flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)", fontsize=10)
    else:
        top_panel.set_ylabel("flux ($e^{-}$)", fontsize=10)

    top_panel.set_title(
        f"Optimally Extracted Order-Merged Object Spectrum ({arm.upper()})\n{filenameTemplate.replace('.fits', '')}",
        fontsize=11,
        linespacing=2.0,
    )

    top_panel.plot(
        merged_orders["WAVE"],
        merged_orders["FLUX_COUNTS"],
        linewidth=0.3,
        color="#dc322f" if not fluxCalibrated else "#2aa198",
        zorder=1,
    )

    _set_wavelength_xlim(top_panel, merged_orders)

    # MIDDLE PANEL WITH LOG SCALE
    middle_panel = fig.add_subplot(gs[1, :])
    if not fluxCalibrated:
        middle_panel.set_ylabel("flux ($e^{-}$)", fontsize=10)
    else:
        middle_panel.set_ylabel("flux (erg s$^{-1}$ cm$^{-2}$ $\\AA^{-1}$)", fontsize=10)
    middle_panel.set_xlabel("wavelength (nm)", fontsize=10)
    middle_panel.set_yscale("log")

    middle_panel.plot(
        merged_orders["WAVE"],
        merged_orders["FLUX_COUNTS"],
        linewidth=0.3,
        color="#dc322f" if not fluxCalibrated else "#2aa198",
        zorder=1,
    )

    # SIGMA-CLIP THE DATA
    arrayMask = sigma_clip(
        merged_orders["FLUX_COUNTS"],
        sigma_lower=3,
        sigma_upper=15.0,
        maxiters=1,
        cenfunc="mean",
        stdfunc="std",
    )
    mean, median, std = sigma_clipped_stats(
        merged_orders["FLUX_COUNTS"],
        sigma=5.0,
        stdfunc="std",
        cenfunc="mean",
        maxiters=3,
    )
    maxFlux = arrayMask.max() + 3 * std
    minFlux = arrayMask.min() - 3 * std

    top_panel.set_ylim(minFlux, maxFlux)

    middle_panel.set_ylim(max(arrayMask.min() * 0.5, 0), arrayMask.max() * 2)
    _set_wavelength_xlim(middle_panel, merged_orders)

    return top_panel, middle_panel


def _draw_merged_snr_panel(log, fig, gs, merged_orders, qcTable):
    """*draw the SNR panel of the merged-spectrum QC plot, annotated with the per-order SNR values*

    **Key Arguments:**

    - ``log`` -- logger
    - ``fig`` -- the QC figure
    - ``gs`` -- the figure's gridspec
    - ``merged_orders`` -- the order-merged spectrum
    - ``qcTable`` -- the QC table holding the SNR values to annotate, or *False* for none

    **Return:**

    - ``bottom_panel`` -- the SNR panel
    """
    from astropy.stats import sigma_clipped_stats

    # BOTTOM PANEL WITH LINEAR SCALE FOR SNR
    bottom_panel = fig.add_subplot(gs[2, :])
    bottom_panel.set_ylabel("SNR", fontsize=10)
    bottom_panel.set_xlabel("wavelength (nm)", fontsize=10)

    if not isinstance(qcTable, bool):
        orderValue, snrValue = _snr_orders_and_values(log, qcTable)

    bottom_panel.plot(
        merged_orders["WAVE"],
        merged_orders["SNR"],
        linewidth=0.4,
        color="black",
        zorder=1,
    )
    _set_wavelength_xlim(bottom_panel, merged_orders)

    mean, median, std = sigma_clipped_stats(merged_orders["SNR"], sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)

    bottom_panel.set_ylim(0, mean + 4 * std)

    # ADD SNR VALUES TO BOTTOM PANEL
    if not isinstance(qcTable, bool) and len(orderValue):
        _annotate_snr_values(bottom_panel, orderValue, snrValue)

    return bottom_panel


def _set_wavelength_xlim(panel, merged_orders):
    """*set a merged-spectrum QC panel's x-axis limits to the spectrum's wavelength range*

    **Key Arguments:**

    - ``panel`` -- the axes to set the limits on
    - ``merged_orders`` -- the merged spectrum, with a ``WAVE`` column
    """
    try:
        panel.set_xlim(merged_orders["WAVE"].min().value, merged_orders["WAVE"].max().value)
    except Exception:
        panel.set_xlim(merged_orders["WAVE"].min(), merged_orders["WAVE"].max())


def _snr_orders_and_values(log, qcTable):
    """*pull the per-order median SNR values out of a QC table for the merged-spectrum QC plot*

    **Key Arguments:**

    - ``log`` -- logger
    - ``qcTable`` -- the QC table holding the ``SNR MEDIAN`` rows

    **Return:**

    - ``orderValue`` -- the order of each SNR value, with ``GLOBAL`` for the whole-spectrum value
    - ``snrValue`` -- the SNR values
    """
    import numpy as np
    import pandas as pd

    qcTable = qcTable.drop_duplicates(subset=["qc_name", "qc_order"], keep="last")
    snrValue = qcTable.loc[qcTable["qc_name"] == "SNR MEDIAN"]
    orderValue = snrValue["qc_order"].values
    for i in range(len(orderValue)):
        try:
            orderValue[i] = int(orderValue[i])
        except (ValueError, TypeError) as e:
            log.debug(f"plot_merged_spectrum_qc: `orderValue[i] = int(orderValue[i])` failed, continuing: {e}")

    orderValue = np.array(["GLOBAL" if pd.isna(v) else v for v in orderValue])
    snrValue = snrValue["qc_value"].values
    return orderValue, snrValue


def _annotate_snr_values(bottom_panel, orderValue, snrValue):
    """*write the per-order median SNR values in the corner of the merged-spectrum SNR panel*

    **Key Arguments:**

    - ``bottom_panel`` -- the SNR panel
    - ``orderValue`` -- the order of each SNR value
    - ``snrValue`` -- the SNR values
    """
    _vis_rank = {"GLOBAL": 0, "u": 1, "g": 2, "r": 3, "i": 4}

    def _order_key(pair):
        o = pair[0]
        if o in _vis_rank:
            return (0, _vis_rank[o])
        try:
            return (1, float(o))
        except (ValueError, TypeError):
            return (2, str(o))

    pairs = sorted(zip(orderValue, snrValue, strict=False), key=_order_key)
    snr_text = "\n".join(f"{o}: {v:.0f}" for o, v in pairs)
    bottom_panel.text(
        0.99,
        0.98,
        snr_text,
        transform=bottom_panel.transAxes,
        ha="right",
        va="top",
        fontsize=5,
        family="monospace",
        color="black",
        zorder=1,
    )


def _draw_sky_panel(fig, gs, merged_orders, fluxCalibrated):
    """*add the sky-flux panel to the merged-spectrum QC plot*

    **Key Arguments:**

    - ``fig`` -- the QC figure
    - ``gs`` -- the figure's grid spec
    - ``merged_orders`` -- the merged spectrum
    - ``fluxCalibrated`` -- whether the spectrum is flux calibrated

    **Return:**

    - ``sky_panel`` -- the new sky panel
    """
    sky_panel = fig.add_subplot(gs[3, :])
    sky_panel.set_xlabel("wavelength (nm)", fontsize=10)

    sky_panel.set_ylabel("sky flux ($e^{-}$)", fontsize=10)

    sky_panel.set_yscale("log")

    if "SKY_COUNTS" in merged_orders.columns and merged_orders["SKY_COUNTS"].max()  > 0:
        sky_panel.plot(
            merged_orders["WAVE"],
            merged_orders["SKY_COUNTS"],
            linewidth=0.3,
            color="#859900" if not fluxCalibrated else "#2aa198",
            zorder=1,
        )

        sky_panel.set_ylim(10, merged_orders["SKY_COUNTS"].max() * 1.1)

    _set_wavelength_xlim(sky_panel, merged_orders)
    return sky_panel


def _mark_order_joins(orderJoins, panels):
    """*mark each order join on the merged-spectrum QC panels*

    **Key Arguments:**

    - ``orderJoins`` -- a dictionary of order-join wavelengths
    - ``panels`` -- the panels to mark
    """
    for v in orderJoins.values():
        for panel in panels:
            panel.axvline(v, color="black", linestyle="--", linewidth=0.5, alpha=0.5)
            panel.text(
                v + 5,
                0.9 * panel.get_ylim()[1],
                "ORDER JOIN",
                rotation=90,
                verticalalignment="top",
                fontsize=6,
                color="black",
                alpha=0.5,
            )


def _mark_skylines(skylinesDF, panels, labelPanel):
    """*mark the sky lines as faint vertical lines on the merged-spectrum QC panels*

    **Key Arguments:**

    - ``skylinesDF`` -- the sky-line table, with ``WAVELENGTH`` and ``ISOLATED`` columns
    - ``panels`` -- the panels to mark
    - ``labelPanel`` -- the one panel whose lines carry a legend label
    """
    import pandas as pd

    # `== True` IS KEPT: A BARE TRUTH MASK WOULD TREAT A NON-BOOLEAN ISOLATED COLUMN DIFFERENTLY
    mask = skylinesDF["ISOLATED"] == True  # noqa: E712
    calibrationSkylines = pd.to_numeric(skylinesDF.loc[mask, "WAVELENGTH"], errors="coerce").dropna().to_numpy()
    otherSkylines = pd.to_numeric(skylinesDF.loc[~mask, "WAVELENGTH"], errors="coerce").dropna().to_numpy()

    for ww in calibrationSkylines:
        for panel in panels:
            panel.axvline(
                ww,
                color="blue",
                linestyle="-",
                linewidth=0.2,
                alpha=0.08,
                zorder=0,
                label="calibration skyline" if panel == labelPanel else "",
            )
    for ww in otherSkylines:
        for panel in panels:
            panel.axvline(
                ww,
                color="grey",
                linestyle="-",
                linewidth=0.2,
                alpha=0.08,
                zorder=0,
                label="skyline" if panel == labelPanel else "",
            )


def calculate_rolling_snr(dataframe, flux_column, window_size):
    """*calculate the rolling signal-to-noise ratio (SNR) for a given column in a pandas dataframe*

    The SNR is calculated as the ratio of the median signal to the noise, where the noise is estimated using a robust
    statistical method.

    **Key Arguments:**

    - ``dataframe`` -- the input pandas dataframe containing the data
    - ``flux_column`` -- the name of the column in the dataframe for which the rolling SNR will be calculated
    - ``window_size`` -- the size of the rolling window to use for the calculation

    **Return:**

    - ``dataframe`` -- the input dataframe with an additional column 'SNR' containing the calculated rolling SNR
      values

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import calculate_rolling_snr
    df_with_snr = calculate_rolling_snr(dataframe=df, flux_column='flux', window_size=5)
    ```
    """
    import numpy as np
    from numpy.lib.stride_tricks import sliding_window_view

    values = dataframe[flux_column].to_numpy(dtype=np.float64, copy=False)
    n = values.size
    snr_full = np.full(n, np.nan, dtype=np.float64)

    # NEED AT LEAST 5 POINTS FOR THE NOISE ESTIMATOR AND ENOUGH DATA FOR ONE WINDOW
    if window_size >= 5 and n >= window_size:
        windows = sliding_window_view(values, window_shape=window_size)

        # SIGNAL: ROLLING MEDIAN
        signal = np.median(windows, axis=1)

        # NOISE: ROBUST ESTIMATOR BASED ON 5-POINT SECOND-DIFFERENCE PATTERN
        diff = np.abs(2.0 * windows[:, 2:-2] - windows[:, :-4] - windows[:, 4:])
        noise = 0.6052697 * np.median(diff, axis=1)

        snr = np.divide(signal, noise, out=np.full_like(signal, np.nan), where=noise > 0)

        # MATCH `center=True` PLACEMENT
        left = (window_size - 1) // 2
        snr_full[left : left + snr.size] = snr

    dataframe["SNR"] = snr_full
    dataframe["SNR"] = dataframe["SNR"].replace([np.inf, -np.inf], np.nan).bfill().ffill()
    return dataframe


def extinction_correction_factor(wave, extinctionTablePath, airmass):
    """*calculate the atmospheric extinction correction factor at each wavelength*

    **Key Arguments:**

    - ``wave`` -- the wavelengths to calculate the factor at, in nm
    - ``extinctionTablePath`` -- path to the observatory extinction curve FITS table (wavelength in Angstrom)
    - ``airmass`` -- the airmass of the observation

    **Return:**

    - ``extCorrectionFactor`` -- the multiplicative correction factor at each wavelength

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import extinction_correction_factor
    extCorrectionFactor = extinction_correction_factor(wave=wave, extinctionTablePath=extinctionTablePath, airmass=1.2)
    ```
    """
    import numpy as np
    from astropy.table import Table
    from scipy.interpolate import interp1d

    # READ THE EXTINCTION CURVE FOR THE OBSERVATORY
    # DATA IS ORGANIZED AS FOLLOWS:
    # FIRST COLUMN, WAVELENGTH (IN ANGSTROM), SECOND COLUMN MAG/AIRMASS
    extinctionData = Table.read(extinctionTablePath, format="fits")
    extinctionData = extinctionData.to_pandas()

    # CONVERT ANG TO NM
    wave_ext = extinctionData["WAVE"] / 10
    # INTERPOLATING ON THE REQUIRED WAVE SCALE
    refitted_ext = interp1d(
        np.array(wave_ext),
        np.array(extinctionData["MAG_AIRMASS"]),
        kind="next",
        fill_value="extrapolate",
    )

    return 10 ** (0.4 * refitted_ext(wave) * airmass)


def frame_to_32(frame):
    """*convert a given frame to 32-bit float format*

    **Key Arguments:**

    - ``frame`` -- the input frame

    **Return:**

    - ``frame`` -- the converted frame

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import frame_to_32
    frame = frame_to_32(frame)
    ```
    """
    import numpy as np

    try:
        if frame.data.dtype != np.float32:
            frame.data = frame.data.astype(np.float32, copy=False)
    except AttributeError as e:
        # NO SOXSPIPE LOGGER IS IN SCOPE IN THIS MODULE-LEVEL HELPER
        logging.getLogger(__name__).debug(f"frame_to_32: could not cast the frame data to float32, continuing: {e}")

    try:
        frame.uncertainty.array = frame.uncertainty.array.astype(np.float32, copy=False)
    except AttributeError as e:
        # NO SOXSPIPE LOGGER IS IN SCOPE IN THIS MODULE-LEVEL HELPER
        logging.getLogger(__name__).debug(
            f"frame_to_32: could not cast the frame uncertainty array to float32, continuing: {e}"
        )

    return frame


def add_snr_efficiency_qcs(log, spectrumDF, qcTable, orderJoins, recipeName, dateObs):
    """*add quality checks to the qc table*

    **Key Arguments:**

    - ``log`` -- logger
    - ``spectrumDF`` -- the dataframe containing the extracted spectrum with a column named 'SNR' or 'EFFICIENCY'
      to calculate the checks from.
    - ``qcTable`` -- the qc table to which the SNR checks will be added
    - ``orderJoins`` -- a dictionary containing the wavelengths of order joins (if any) to be added as QCs.
    - ``recipeName`` -- name of the recipe to add to the QC entries
    - ``dateObs`` -- the observation date to add to the QC entries (in ISO format, e.g., "2024-06-01T00:00:00")

    **Return:**

    - ``qcTable`` -- the updated qc table with SNR checks added

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import add_snr_efficiency_qcs
    qcTable = add_snr_efficiency_qcs(log, spectrumDF, qcTable, orderJoins, recipeName, dateObs)
    ```
    """
    import numpy as np
    import pandas as pd

    spectrumDF = spectrumDF.copy(deep=False)
    spectrumDF["ORDER"] = np.nan
    try:
        spectrumDF["WAVE"] = spectrumDF["WAVE"].values.value
    except AttributeError as e:
        log.debug(f"add_snr_efficiency_qcs: `spectrumDF['WAVE'] = spectrumDF['WAVE'].value...` failed, continuing: {e}")

    ## REVERSE DICTIONARY KEYS SO FIRST KEY IS LAST
    orderJoins = dict(reversed(list(orderJoins.items())))

    for i, (orders, wl) in enumerate(orderJoins.items()):
        if len(orders) == 2:
            visOrders = ["i", "r", "g", "u"]
            lastOrder = visOrders[int(orders[0])]
            order = visOrders[int(orders[0]) - 1]

        elif len(orders) == 3:
            lastOrder = int(orders[:1])
            order = int(orders[1:])
        else:
            order = int(orders[:2])
            lastOrder = int(orders[2::])
        if i == 0:
            spectrumDF.loc[(spectrumDF["WAVE"] <= wl), "ORDER"] = lastOrder
        spectrumDF.loc[
            spectrumDF["WAVE"] >= wl,
            "ORDER",
        ] = order

    utcnow = utcnow_string()

    # CALCULATE THE MEDIAN EFFICIENCY ACROSS ALL ORDERS
    if "EFFICIENCY" in spectrumDF.columns:
        medianEfficiency = spectrumDF["EFFICIENCY"].median()
        qcTable = pd.concat(
            [
                qcTable,
                pd.DataFrame([
                    {
                        "soxspipe_recipe": recipeName,
                        "qc_name": "EFF MEDIAN",
                        "qc_value": float(f"{medianEfficiency:0.4f}"),
                        "qc_comment": "Median efficiency across all orders",
                        "qc_unit": None,
                        "obs_date_utc": dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                ]),
            ],
            ignore_index=True,
        )
        # CALCULATE THE MEDIAN EFFICIENCY IN EACH ORDER
        orderEfficiency = spectrumDF.groupby("ORDER")["EFFICIENCY"].median().reset_index()
        orderEfficiency.columns = ["ORDER", "MEDIAN_EFFICIENCY"]

        for _, row in orderEfficiency.iterrows():
            qcTable = pd.concat(
                [
                    qcTable,
                    pd.DataFrame([
                        {
                            "soxspipe_recipe": recipeName,
                            "qc_name": "EFF MEDIAN",
                            "qc_value": float(f"{row['MEDIAN_EFFICIENCY']:0.4f}"),
                            "qc_comment": f"Median efficiency in order {row['ORDER']}",
                            "qc_order": row["ORDER"],
                            "qc_unit": None,
                            "obs_date_utc": dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]),
                ],
                ignore_index=True,
            )

    if "SNR" in spectrumDF.columns:
        medianSNR = spectrumDF["SNR"].median()
        qcTable = pd.concat(
            [
                qcTable,
                pd.DataFrame([
                    {
                        "soxspipe_recipe": recipeName,
                        "qc_name": "SNR MEDIAN",
                        "qc_value": float(f"{medianSNR:0.3f}"),
                        "qc_comment": "Median SNR across all orders",
                        "qc_unit": None,
                        "obs_date_utc": dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                ]),
            ],
            ignore_index=True,
        )
        orderSNR = spectrumDF.groupby("ORDER")["SNR"].median().reset_index()
        orderSNR.columns = ["ORDER", "MEDIAN_SNR"]

        for _, row in orderSNR.iterrows():
            qcTable = pd.concat(
                [
                    qcTable,
                    pd.DataFrame([
                        {
                            "soxspipe_recipe": recipeName,
                            "qc_name": "SNR MEDIAN",
                            "qc_value": float(f"{row['MEDIAN_SNR']:0.3f}"),
                            "qc_comment": f"Median SNR in order {row['ORDER']}",
                            "qc_order": row["ORDER"],
                            "qc_unit": None,
                            "obs_date_utc": dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]),
                ],
                ignore_index=True,
            )

    return qcTable


################# SHARED QC / PRODUCT / PLOT HELPERS ####################
# THESE HELPERS FACTOR OUT THE ROW-BUILDING BOILERPLATE CURRENTLY DUPLICATED
# AT DOZENS OF CALL SITES ACROSS THE RECIPES. THEY ADD CODE ONLY -- NO
# EXISTING CALL SITE IS CONVERTED TO USE THEM HERE.


def utcnow_string(microseconds=False):
    """*return the current UTC timestamp formatted the way the codebase's
    40 existing ``datetime.utcnow().strftime(...)`` call sites already do*

    The rendered string is byte-identical to those call sites because the
    format string contains no `%z`/`%Z` offset marker -- swapping the naive
    ``datetime.utcnow()`` for the timezone-aware ``datetime.now(UTC)`` (this
    module already imports ``UTC`` at the top) changes nothing about the
    rendered text, only about whether the underlying `datetime` object
    carries timezone info.

    **Key Arguments:**

    - ``microseconds`` -- also render the `%f` fractional-second component. Default *False*

    **Return:**

    - ``timestamp`` -- the formatted UTC timestamp string

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import utcnow_string
    utcnow = utcnow_string()
    ```
    """
    timestampFormat = "%Y-%m-%dT%H:%M:%S"
    if microseconds:
        timestampFormat += ".%f"
    return datetime.now(UTC).strftime(timestampFormat)


class _omitted:
    def __repr__(self):
        return "<omitted>"


# SENTINEL USED TO DISTINGUISH "CALLER DID NOT PASS THIS KEYWORD" FROM
# "CALLER PASSED None" -- None IS A REAL VALUE IN THESE QC/PRODUCT COLUMNS
# (E.G. `soxs_mflat.py:1225` BUILDS A QC ROW WITH NO `to_header` KEY AT ALL;
# A HELPER THAT DEFAULTED `toHeader=True` WOULD FLIP THAT ROW FROM NaN TO
# True, WHICH IS A REAL BEHAVIOUR CHANGE), SO "NOT PASSED" AND "PASSED AS
# None" MUST STAY DISTINGUISHABLE.
#
# THE HELPERS BELOW TEST FOR IT BY IDENTITY (`value is not OMITTED`), NOT BY
# `isinstance(value, _omitted)`: THIS ONE OBJECT IS THE SENTINEL, AND A SECOND
# `_omitted()` ARRIVING FROM ANYWHERE -- A deepcopy, A pickle ROUND TRIP, A
# FUTURE REFACTOR -- MUST NOT READ AS "OMITTED".
OMITTED = _omitted()


def append_qc(
    qcTable,
    recipeName,
    qcName,
    qcValue,
    qcComment,
    obsDateUtc,
    reductionDateUtc,
    qcUnit=OMITTED,
    toHeader=OMITTED,
    qcOrder=OMITTED,
):
    """*append a single QC row to a QC table without mutating the input table*

    ``reductionDateUtc`` is a required argument and this helper never mints
    its own timestamp: call sites deliberately share one timestamp across
    several rows (e.g. `soxs_mbias.py:329` covers two QC rows,
    `soxs_stare.py:510` covers three product rows), and a per-call
    timestamp would split those apart.

    The optional keyword arguments (``qcUnit``, ``toHeader``, ``qcOrder``)
    only appear in the built row when the caller actually passes them --
    passing ``None`` explicitly still adds the column (with a `None`/NaN
    value); not passing it at all omits the column entirely.

    **Key Arguments:**

    - ``qcTable`` -- the QC table to append to (not mutated)
    - ``recipeName`` -- the recipe name to record against the QC row
    - ``qcName`` -- the QC metric name
    - ``qcValue`` -- the QC metric value
    - ``qcComment`` -- the QC metric comment
    - ``obsDateUtc`` -- the observation date (UTC) to record against the QC row
    - ``reductionDateUtc`` -- the reduction date (UTC) to record against the QC row. Shared
      across several rows by the caller, never generated here
    - ``qcUnit`` -- the QC metric unit. Omit to leave the column absent for this row
    - ``toHeader`` -- whether the QC metric should be written to the FITS header. Omit to
      leave the column absent for this row
    - ``qcOrder`` -- the echelle order the QC metric applies to. Omit to leave the column absent for this row

    **Return:**

    - ``qcTable`` -- a new QC table with the row appended

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import append_qc
    self.qc = append_qc(
        self.qc,
        recipeName=self.recipeName,
        qcName="RON",
        qcValue=1.2,
        qcComment="[e-] RON in single BIAS",
        obsDateUtc=self.dateObs,
        reductionDateUtc=utcnow,
        qcUnit="electron",
        toHeader=True,
    )
    ```
    """
    import pandas as pd

    row = {
        "soxspipe_recipe": recipeName,
        "qc_name": qcName,
        "qc_value": qcValue,
    }
    if qcUnit is not OMITTED:
        row["qc_unit"] = qcUnit
    if qcOrder is not OMITTED:
        row["qc_order"] = qcOrder
    row["qc_comment"] = qcComment
    row["obs_date_utc"] = obsDateUtc
    row["reduction_date_utc"] = reductionDateUtc
    if toHeader is not OMITTED:
        row["to_header"] = toHeader

    return pd.concat([qcTable, pd.DataFrame([row])], ignore_index=True)


def append_product(
    productsTable,
    recipeName,
    productLabel,
    fileName,
    filePath,
    productDesc,
    obsDateUtc,
    reductionDateUtc,
    fileType=OMITTED,
    label=OMITTED,
):
    """*append a single product row to a products table without mutating the input table*

    ``recipeName`` is required and deliberately NOT read from `self`: eight
    rows across the codebase hardcode the literal `"soxs-stare"` regardless
    of the running recipe (`soxs_stare.py:516`, `:547`, `:579`,
    `subtract_sky.py:307`, `:339`, `horne_extraction.py:444`, `:532`,
    `:1366`), and `base_recipe.py` rewrites `self.recipeName` to a `-std`
    variant for standard-star input.

    As with `append_qc`, ``reductionDateUtc`` is required and never minted
    here, and the optional keyword arguments (``fileType``, ``label``) only
    appear in the built row when the caller actually passes them.

    **Key Arguments:**

    - ``productsTable`` -- the products table to append to (not mutated)
    - ``recipeName`` -- the recipe name to record against the product row
    - ``productLabel`` -- the product label
    - ``fileName`` -- the product file name
    - ``filePath`` -- the product file path
    - ``productDesc`` -- the product description
    - ``obsDateUtc`` -- the observation date (UTC) to record against the product row
    - ``reductionDateUtc`` -- the reduction date (UTC) to record against the product row. Shared
      across several rows by the caller, never generated here
    - ``fileType`` -- the product file type. Omit to leave the column absent for this row
    - ``label`` -- the product label category (e.g. `PROD`). Omit to leave the column absent for this row

    **Return:**

    - ``productsTable`` -- a new products table with the row appended

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import append_product
    self.products = append_product(
        self.products,
        recipeName=self.recipeName,
        productLabel="MBIAS",
        fileName=filename,
        filePath=productPath,
        productDesc=f"{self.arm} Master bias frame",
        obsDateUtc=self.dateObs,
        reductionDateUtc=utcnow,
        fileType="FITS",
        label="PROD",
    )
    ```
    """
    import pandas as pd

    row = {
        "soxspipe_recipe": recipeName,
        "product_label": productLabel,
        "file_name": fileName,
    }
    if fileType is not OMITTED:
        row["file_type"] = fileType
    row["obs_date_utc"] = obsDateUtc
    row["reduction_date_utc"] = reductionDateUtc
    row["product_desc"] = productDesc
    row["file_path"] = filePath
    if label is not OMITTED:
        row["label"] = label

    return pd.concat([productsTable, pd.DataFrame([row])], ignore_index=True)


def save_qc_plot(filePath, dpi=120, bboxInches="tight"):
    """*save the current matplotlib figure as a QC PDF plot*

    ``bboxInches=None`` means OMIT the `bbox_inches` keyword from the
    `plt.savefig` call entirely, not pass `bbox_inches=None` -- that is what
    reproduces the existing `subtract_sky.py` and `create_dispersion_map.py`
    call sites exactly, where `bbox_inches` is never passed at all.

    **Key Arguments:**

    - ``filePath`` -- the path to save the plot to
    - ``dpi`` -- the plot resolution in dots per inch. Default *120*
    - ``bboxInches`` -- the `bbox_inches` value to forward to `plt.savefig`, or `None` to omit
      the keyword entirely. Default *"tight"*

    **Return:**

    - ``filePath`` -- the path the plot was saved to

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import save_qc_plot
    save_qc_plot(filePath)
    ```
    """
    import matplotlib.pyplot as plt

    kwargs = {"dpi": dpi, "format": "pdf"}
    if bboxInches is not None:
        kwargs["bbox_inches"] = bboxInches
    plt.savefig(filePath, **kwargs)
    return filePath


def get_skylines_dataframe(log, settings, arm, minBrightnessVIS=5, minBrightnessNIR=100):
    """*load the static skyline table for an arm and keep only the strong skylines, for QC plotting*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the soxspipe settings dictionary
    - ``arm`` -- the spectrograph arm
    - ``minBrightnessVIS`` -- in the VIS arm, keep only skylines with flux above this. Default *5*
    - ``minBrightnessNIR`` -- in every other arm, keep only skylines with flux above this. Default *100*

    **Return:**

    - ``skylinesDF`` -- the strong skylines as a dataframe

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import get_skylines_dataframe
    skylinesDF = get_skylines_dataframe(log, settings, arm)
    ```
    """
    from astropy.table import Table

    from soxspipe.commonutils import detector_lookup
    from soxspipe.commonutils.toolkit import get_calibrations_path

    dp = detector_lookup(log=log, settings=settings).get(arm)
    calibrationRootPath = get_calibrations_path(log=log, settings=settings)
    skylines = calibrationRootPath + "/" + dp["skylines"]

    dat = Table.read(skylines, format="fits")
    skylinesDF = dat.to_pandas()
    skylinesDF["WAVELENGTH"] = skylinesDF["WAVELENGTH"].astype(float)
    skylinesDF["FLUX"] = skylinesDF["FLUX"].astype(float)
    if "ISOLATED" in skylinesDF.columns:
        skylinesDF["ISOLATED"] = skylinesDF["ISOLATED"].astype(bool)

    mask = skylinesDF["FLUX"] > minBrightnessVIS if arm == "VIS" else skylinesDF["FLUX"] > minBrightnessNIR

    return skylinesDF.loc[mask]
