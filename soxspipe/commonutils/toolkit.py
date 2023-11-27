#!/usr/bin/env python
# encoding: utf-8
"""
*small reusable functions used throughout soxspipe*

:Author:
    David Young

:Date Created:
    September 18, 2020
"""


from os.path import expanduser
from soxspipe.commonutils import detector_lookup
from datetime import datetime
from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial, chebyshev_order_xy_polynomials
from fundamentals import tools
from builtins import object
import sys
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
import os


os.environ['TERM'] = 'vt100'


def cut_image_slice(
        log,
        frame,
        width,
        length,
        x,
        y,
        sliceAxis="x",
        median=False,
        plot=False):
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
    - ``plot`` -- generate a plot of slice. Useful for debugging.

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import cut_image_slice
    slice = cut_image_slice(log=self.log, frame=self.pinholeFlat.data,
                                    width=1, length=sliceLength, x=x_fit, y=y_fit, plot=False)
    if slice is None:
        return None
    ```
    """
    log.debug('starting the ``cut_image_slice`` function')

    import numpy.ma as ma
    import numpy as np
    import random

    halfSlice = length / 2
    # NEED AN EVEN PIXEL SIZE
    if (width % 2) != 0:
        halfwidth = (width - 1) / 2
    else:
        halfwidth = width / 2

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
        raise ValueError("sliceAxis needs to be eith 'x' or 'y'")

    # CHECK WE ARE NOT GOING BEYOND BOUNDS OF FRAME
    if (axisA > axisALen - halfSlice) or (axisB > axisBLen - halfwidth) or (axisA < halfSlice) or (axisB < halfwidth):
        return None, None, None

    slice_length_offset = int(axisA - halfSlice)
    if sliceAxis == "x":
        sliceFull = frame[int(axisB - halfwidth):int(axisB + halfwidth + 1),
                          slice_length_offset:int(axisA + halfSlice)]
    else:
        sliceFull = frame[slice_length_offset:int(axisA + halfSlice), int(axisB - halfwidth):int(axisB + halfwidth + 1)]
    slice_width_centre = (int(axisB + halfwidth + 1) + int(axisB - halfwidth)) / 2

    if median:
        if sliceAxis == "y":
            slice = ma.median(sliceFull, axis=1)
        else:
            slice = ma.median(sliceFull, axis=0)

    if plot and random.randint(1, 101) < 10000:
        import matplotlib.pyplot as plt
        # CHECK THE SLICE POINTS IF NEEDED
        if sliceAxis == "y":
            sliceImg = np.rot90(sliceFull, 1)
        else:
            sliceImg = sliceFull
        plt.imshow(sliceImg)
        plt.show()
        xx = np.arange(0, len(slice))
        plt.figure(figsize=(8, 5))
        if sliceAxis == "y":
            plt.plot(xx, slice, 'ko', label=f"x={axisB}, y={axisA}, sliceAxis={sliceAxis}")
        if sliceAxis == "x":
            plt.plot(xx, slice, 'ko', label=f"x={axisA}, y={axisB}, sliceAxis={sliceAxis}")
        plt.xlabel('Position')
        plt.ylabel('Flux')
        plt.legend()
        plt.show()

    log.debug('completed the ``cut_image_slice`` function')
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
        skylines=False):
    """*generate a quicklook image of a CCDObject - useful for development/debugging*

    **Key Arguments:**

    - ``log`` -- logger
    - ``CCDObject`` -- the CCDObject to plot
    - ``show`` -- show the image. Set to False to skip
    - ``ext`` -- the name of the the extension to show. Can be "data", "mask" or "err". Default "data".
    - ``title`` -- give a title for the plot
    - ``surfacePlot`` -- plot as a 3D surface plot
    - ``dispMap`` -- path to dispersion map. Default *False*
    - ``dispMapImage`` -- the 2D dispersion map image
    - ``inst`` -- provide instrument name if no header exists
    - ``skylines`` -- mark skylines on image

    ```python
    from soxspipe.commonutils.toolkit import quicklook_image
    quicklook_image(
        log=self.log, CCDObject=myframe, show=True)
    ```
    """
    log.debug('starting the ``quicklook_image`` function')

    import pandas as pd
    import matplotlib as mpl
    import numpy as np
    from copy import copy
    from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
    from soxspipe.commonutils import keyword_lookup
    from soxspipe.commonutils import detector_lookup

    if settings:
        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=log,
            settings=settings
        ).get
        arm = CCDObject.header[kw("SEQ_ARM")]
        dateObs = CCDObject.header[kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(arm)

        # USE THIS ELSEWHERE IN THE OBJECT METHODS
        dp = detectorParams
        science_pixels = dp["science-pixels"]

    originalRC = dict(mpl.rcParams)

    if not show:
        return
    import matplotlib.pyplot as plt

    if ext == "data":
        frame = CCDObject.data
    elif ext == "mask":
        frame = CCDObject.mask
    elif ext == "uncertainty":
        frame = CCDObject.uncertainty.array
    else:
        # ASSUME ONLY NDARRAY
        frame = CCDObject

    if inst is False:
        try:
            inst = CCDObject.header["INSTRUME"]
        except:
            inst = "XSHOOTER"

    if inst == "SOXS":
        rotatedImg = np.flipud(frame)
    elif inst == "XSHOOTER":
        rotatedImg = np.rot90(frame, 1)
    else:
        rotatedImg = frame
    rotatedImg = np.flipud(rotatedImg)

    std = np.nanstd(frame)
    mean = np.nanmean(frame)
    median = np.nanmean(frame)
    palette = copy(plt.cm.viridis)
    palette.set_bad("#dc322f", 1.0)
    vmax = median + stdWindow * 0.5 * std
    vmin = median - stdWindow * 0.5 * std

    if surfacePlot:

        from matplotlib import rc

        axisColour = '#dddddd'
        rc('axes', edgecolor=axisColour, labelcolor=axisColour, linewidth=0.6)
        rc('xtick', color=axisColour)
        rc('ytick', color=axisColour)
        rc('grid', color=axisColour)
        rc('text', color=axisColour)

        fig = plt.figure(figsize=(40, 10))
        ax = fig.add_subplot(121, projection='3d')
        if inst == "XSHOOTER":
            plt.gca().invert_yaxis()
        ax.set_box_aspect(aspect=(2, 1, 1))
        # Remove gray panes and axis grid
        ax.xaxis.pane.fill = False
        ax.zaxis.pane.set_facecolor("#dc322f")
        ax.zaxis.pane.set_alpha(1.)
        ax.yaxis.pane.fill = False

        ax.grid(False)
        # Remove z-axis
        # ax.w_zaxis.line.set_lw(0.)
        # ax.set_zticks([])

        X, Y = np.meshgrid(np.linspace(0, rotatedImg.shape[1], rotatedImg.shape[1]), np.linspace(0, rotatedImg.shape[0], rotatedImg.shape[0]))
        surface = ax.plot_surface(X=X, Y=Y, Z=rotatedImg, cmap='viridis', antialiased=True, vmin=vmin, vmax=vmax)

        if inst == "SOXS":
            ax.azim = 70
        else:
            ax.azim = -120
        ax.elev = 30

        ax.set_xlim(0, rotatedImg.shape[1])
        ax.set_ylim(0, rotatedImg.shape[0])
        # ax.set_zlim(0, min(np.nanmax(frame), mean + stdWindow * 10 * std))

        if inst == "SOXS":
            ax.invert_yaxis()
        backgroundColour = '#404040'
        fig.set_facecolor(backgroundColour)
        ax.set_facecolor(backgroundColour)
        ax.xaxis.pane.set_edgecolor(backgroundColour)
        ax.yaxis.pane.set_edgecolor(backgroundColour)
        ax.zaxis.pane.set_edgecolor(backgroundColour)

        if inst == "SOXS":
            plt.xlabel(
                "x-axis", fontsize=10)
            plt.ylabel(
                "y-axis", fontsize=10)
        else:
            plt.xlabel(
                "y-axis", fontsize=10)
            plt.ylabel(
                "x-axis", fontsize=10)

        ax2 = fig.add_subplot(122)
    else:
        fig = plt.figure(figsize=(12, 5))

        # palette.set_over('r', 1.0)
        # palette.set_under('g', 1.0)
        ax2 = fig.add_subplot(111)

    if skylines:
        calibrationRootPath = get_calibrations_path(log=log, settings=settings)
        skylines = calibrationRootPath + "/" + dp["skylines"]
        # SPEC FORMAT TO PANDAS DATAFRAME
        from astropy.table import Table
        dat = Table.read(skylines, format='fits')
        skylinesDF = dat.to_pandas()

    if not isinstance(dispMapImage, bool):

        gridLinePixelTable = create_dispersion_solution_grid_lines_for_plot(
            log=log,
            dispMap=dispMap,
            dispMapImage=dispMapImage,
            associatedFrame=CCDObject,
            kw=kw,
            skylines=skylines
        )

        for l in range(int(gridLinePixelTable['line'].max())):
            mask = (gridLinePixelTable['line'] == l)
            ax2.plot(gridLinePixelTable.loc[mask]["fit_y"], gridLinePixelTable.loc[mask]["fit_x"], "w-", linewidth=0.5, alpha=0.8, color="black")

    ax2.set_box_aspect(0.5)
    detectorPlot = plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                              cmap=palette, alpha=1, aspect='auto')

    if surfacePlot:
        shrink = 0.5
    else:
        shrink = 1.0

    if mean > 10:
        fmt = '%1.0f'
        fig.colorbar(detectorPlot, shrink=shrink, format=fmt)
    else:
        fig.colorbar(detectorPlot, shrink=shrink)
        # plt.colorbar()
    if title:
        fig.suptitle(title, fontsize=16)
    if inst == "XSHOOTER":
        ax2.invert_yaxis()
    # cbar.ticklabel_format(useOffset=False)
    if inst == "SOXS":
        plt.xlabel(
            "x-axis", fontsize=10)
        plt.ylabel(
            "y-axis", fontsize=10)
    else:
        plt.xlabel(
            "y-axis", fontsize=10)
        plt.ylabel(
            "x-axis", fontsize=10)

    plt.show()
    mpl.rcParams.update(originalRC)

    log.debug('completed the ``quicklook_image`` function')
    return None


def unpack_order_table(
        log,
        orderTablePath,
        extend=0.,
        pixelDelta=1,
        binx=1,
        biny=1,
        prebinned=False,
        order=False):
    """*unpack an order table and return a top-level `orderPolyTable` data-frame and a second `orderPixelTable` data-frame with the central-trace coordinates of each order given

    **Key Arguments:**

    - ``orderTablePath`` -- path to the order table
    - ``extend`` -- fractional increase to the order area in the y-axis (needed for masking)
    - ``pixelDelta`` -- space between returned data points. Default *1* (sampled at every pixel)
    - ``binx`` -- binning in the x-axis (from FITS header). Default *1*
    - ``biny`` -- binning in the y-axis (from FITS header). Default *1*
    - ``prebinned`` -- was the order-table measured on a pre-binned frame (typically only for mflats). Default *False*
    - ``order`` -- unpack only a single order

    **Usage:**

    ```python
    # UNPACK THE ORDER TABLE
    from soxspipe.commonutils.toolkit import unpack_order_table
    orderPolyTable, orderPixelTable = unpack_order_table(
        log=self.log, orderTablePath=orderTablePath, extend=0.)
    ```
    """
    log.debug('starting the ``functionName`` function')
    from astropy.table import Table
    import pandas as pd
    import numpy as np
    import math
    # MAKE RELATIVE HOME PATH ABSOLUTE

    home = expanduser("~")
    orderTablePath = orderTablePath.replace("~", home)

    dat = Table.read(orderTablePath, format='fits', hdu=1)
    orderPolyTable = dat.to_pandas()

    dat = Table.read(orderTablePath, format='fits', hdu=2)
    orderMetaTable = dat.to_pandas()

    if order:
        mask = (orderMetaTable["order"] == order)
        orderMetaTable = orderMetaTable.loc[mask]

    if "degy_cent" in orderPolyTable.columns:
        axisA = "x"
        axisB = "y"
        axisAbin = binx
        axisBbin = biny
    else:
        axisA = "y"
        axisB = "x"
        axisAbin = biny
        axisBbin = binx

    # ADD AXIS B COORD LIST
    if prebinned:
        ratio = axisBbin
    else:
        ratio = 1
    axisBcoords = [np.arange(0 if (math.floor(l) - int(r * extend)) < 0 else (math.floor(l) - int(r * extend)), 4200 if (math.ceil(u) + int(r * extend)) > 4200 else (math.ceil(u) + int(r * extend)), pixelDelta) for l, u, r in zip(
        orderMetaTable[f"{axisB}min"].values * ratio, orderMetaTable[f"{axisB}max"].values * ratio, orderMetaTable[f"{axisB}max"].values * ratio - orderMetaTable[f"{axisB}min"].values * ratio)]
    orders = [np.full_like(a, o) for a, o in zip(
        axisBcoords, orderMetaTable["order"].values)]

    # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
    myDict = {
        f"{axisB}coord": np.concatenate(axisBcoords),
        "order": np.concatenate(orders)
    }
    orderPixelTable = pd.DataFrame(myDict)

    cent_coeff = [float(v) for k, v in orderPolyTable.iloc[0].items() if "cent_" in k]
    poly = chebyshev_order_xy_polynomials(log=log, axisBCol=f"{axisB}coord", orderCol="order", orderDeg=int(orderPolyTable.iloc[0]["degorder_cent"]), axisBDeg=int(orderPolyTable.iloc[0][f"deg{axisB}_cent"])).poly
    orderPixelTable[f"{axisA}coord_centre"] = poly(orderPixelTable, *cent_coeff)

    std_coeff = [float(v) for k, v in orderPolyTable.iloc[0].items() if "std_" in k]
    if len(std_coeff):
        poly = chebyshev_order_xy_polynomials(log=log, axisBDeg=int(orderPolyTable.iloc[0][f"deg{axisB}_cent"]), orderDeg=int(orderPolyTable.iloc[0]["degorder_cent"]), orderCol="order", axisBCol=f"{axisB}coord").poly
        orderPixelTable["std"] = poly(orderPixelTable, *std_coeff)

    if f"deg{axisB}_edgeup" in orderPolyTable.columns:
        upper_coeff = [float(v) for k, v in orderPolyTable.iloc[0].items() if "edgeup_" in k]
        poly = chebyshev_order_xy_polynomials(log=log, axisBDeg=int(orderPolyTable.iloc[0][f"deg{axisB}_edgeup"]), orderDeg=int(orderPolyTable.iloc[0]["degorder_edgeup"]), orderCol="order", axisBCol=f"{axisB}coord").poly
        orderPixelTable[f"{axisA}coord_edgeup"] = poly(orderPixelTable, *upper_coeff)

    if f"deg{axisB}_edgelow" in orderPolyTable.columns:
        upper_coeff = [float(v) for k, v in orderPolyTable.iloc[0].items() if "edgelow_" in k]
        poly = chebyshev_order_xy_polynomials(log=log, axisBDeg=int(orderPolyTable.iloc[0][f"deg{axisB}_edgelow"]), orderDeg=int(orderPolyTable.iloc[0]["degorder_edgelow"]), orderCol="order", axisBCol=f"{axisB}coord").poly
        orderPixelTable[f"{axisA}coord_edgelow"] = poly(orderPixelTable, *upper_coeff)

    if axisAbin != 1:
        orderPixelTable[f"{axisA}coord_centre"] /= axisAbin
        orderPixelTable[f"{axisA}coord_edgeup"] /= axisAbin
        orderPixelTable[f"{axisA}coord_edgelow"] /= axisAbin
    if axisBbin != 1:
        orderMetaTable[f"{axisB}min"] /= axisBbin
        orderMetaTable[f"{axisB}max"] /= axisBbin
        orderPixelTable[f"{axisB}coord"] /= axisBbin
        orderPixelTable["std"] /= axisBbin
        mask = (orderPixelTable[f"{axisB}coord"].mod(1) > 0)
        orderPixelTable = orderPixelTable.loc[~mask]
        orderPixelTable[f"{axisB}coord"] = orderPixelTable[f"{axisB}coord"].round().astype('int')

    log.debug('completed the ``functionName`` function')
    return orderPolyTable, orderPixelTable, orderMetaTable


def generic_quality_checks(
        log,
        frame,
        settings,
        recipeName,
        qcTable):
    """*measure very basic quality checks on a frame and return the QC table with results appended*

    **Key Arguments:**

    - `log` -- logger
    - `frame` -- CCDData object
    - `settings` -- soxspipe settings
    - `recipeName` -- the name of the recipe
    - `qcTable` -- the QC pandas data-frame to save the QC measurements

    **Usage:**

    ```eval_rst
    .. todo::

            add usage info
            create a sublime snippet for usage
    ```

    ```python
    usage code
    ```
    """
    log.debug('starting the ``functionName`` function')

    import numpy as np
    import pandas as pd

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(
        log=log,
        settings=settings
    ).get
    kw = kw
    arm = frame.header[kw("SEQ_ARM")]
    dateObs = frame.header[kw("DATE_OBS")]

    nanCount = np.count_nonzero(np.isnan(frame.data))

    utcnow = datetime.utcnow()
    utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

    qcTable = pd.concat([qcTable, pd.Series({
        "soxspipe_recipe": recipeName,
        "qc_name": "N NAN PIXELS",
        "qc_value": nanCount,
        "qc_comment": "Number of NaN pixels",
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow,
        "to_header": False
    }).to_frame().T], ignore_index=True)

    # COUNT BAD-PIXELS
    badCount = frame.mask.sum()
    totalPixels = np.size(frame.mask)
    percent = (float(badCount) / float(totalPixels))
    percent = float("{:.6f}".format(percent))

    qcTable = pd.concat([qcTable, pd.Series({
        "soxspipe_recipe": recipeName,
        "qc_name": "N BAD PIXELS",
        "qc_value": int(badCount),
        "qc_comment": "Number of bad pixels",
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow,
        "to_header": True
    }).to_frame().T], ignore_index=True)

    qcTable = pd.concat([qcTable, pd.Series({
        "soxspipe_recipe": recipeName,
        "qc_name": "FRAC BAD PIXELS",
        "qc_value": percent,
        "qc_comment": "Fraction of bad pixels",
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow,
        "to_header": True
    }).to_frame().T], ignore_index=True)

    log.debug('completed the ``functionName`` function')
    return qcTable


def spectroscopic_image_quality_checks(
        log,
        frame,
        orderTablePath,
        settings,
        recipeName,
        qcTable):
    """*measure and record spectroscopic image quailty checks*

    **Key Arguments:**

    - `log` -- logger
    - `frame` -- CCDData object
    - ``orderTablePath`` -- path to the order table
    - `settings` -- soxspipe settings
    - `recipeName` -- the name of the recipe
    - `qcTable` -- the QC pandas data-frame to save the QC measurements

    **Usage:**

    ```eval_rst
    .. todo::

            add usage info
            create a sublime snippet for usage
    ```

    ```python
    usage code
    ```
    """
    log.debug('starting the ``functionName`` function')

    import numpy.ma as ma
    import numpy as np
    import pandas as pd

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(
        log=log,
        settings=settings
    ).get
    kw = kw
    arm = frame.header[kw("SEQ_ARM")]
    dateObs = frame.header[kw("DATE_OBS")]

    try:
        binx = frame.header[kw("WIN_BINX")]
        biny = frame.header[kw("WIN_BINY")]
    except:
        if arm.lower() == "nir":
            binx = 1
            biny = 1

    inst = frame.header[kw("INSTRUME")]
    if inst == "SOXS":
        axisA = "y"
        axisB = "x"
    elif inst == "XSHOOTER":
        axisA = "x"
        axisB = "y"

    # UNPACK THE ORDER TABLE
    orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
        log=log, orderTablePath=orderTablePath, binx=binx, biny=biny)

    mask = np.ones_like(frame.data)

    axisACoords_up = orderTablePixels[f"{axisA}coord_edgeup"].values
    axisACoords_low = orderTablePixels[f"{axisA}coord_edgelow"].values
    axisBCoords = orderTablePixels[f"{axisB}coord"].values
    axisACoords_up = axisACoords_up.astype(int)
    axisACoords_low = axisACoords_low.astype(int)

    # UPDATE THE MASK
    for u, l, y in zip(axisACoords_up, axisACoords_low, axisBCoords):
        mask[y][l:u] = 0

    # COMBINE MASK WITH THE BAD PIXEL MASK
    mask = (mask == 1) | (frame.mask == 1)

    # PLOT ONE OF THE MASKED FRAMES TO CHECK
    maskedFrame = ma.array(frame.data, mask=mask)
    quicklook_image(log=log, CCDObject=np.copy(mask),
                    show=False, ext=None)

    mean = np.ma.mean(maskedFrame)
    flux = np.ma.sum(maskedFrame)

    utcnow = datetime.utcnow()
    utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

    qcTable = pd.concat([qcTable, pd.Series({
        "soxspipe_recipe": recipeName,
        "qc_name": "INNER ORDER PIX MEAN",
        "qc_value": mean,
        "qc_comment": "[e-] Mean inner-order pixel value",
        "qc_unit": "electrons",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow,
        "to_header": True
    }).to_frame().T], ignore_index=True)

    qcTable = pd.concat([qcTable, pd.Series({
        "soxspipe_recipe": recipeName,
        "qc_name": "INNER ORDER PIX SUM",
        "qc_value": flux,
        "qc_comment": "[e-] Sum of all inner-order pixel values",
        "qc_unit": "electrons",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow,
        "to_header": True
    }).to_frame().T], ignore_index=True)

    log.debug('completed the ``functionName`` function')
    return qcTable


def read_spectral_format(
        log,
        settings,
        arm):
    """*read the spectral format table to get some key parameters*

    **Key Arguments:**

    - `log` -- logger
    - `settings` -- soxspipe settings
    - `arm` -- arm to retrieve format for

    **Return:**
        - ``orderNums`` -- a list of the order numbers
        - ``waveLengthMin`` -- a list of the maximum wavelengths reached by each order
        - ``waveLengthMax`` -- a list of the minimum wavelengths reached by each order

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import read_spectral_format
    # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
    orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=arm)
    ```
    """
    log.debug('starting the ``read_spectral_format`` function')

    # DETECTOR PARAMETERS LOOKUP OBJECT
    dp = detector_lookup(
        log=log,
        settings=settings
    ).get(arm)

    science_pixels = dp["science-pixels"]

    # READ THE SPECTRAL FORMAT TABLE FILE
    home = expanduser("~")

    calibrationRootPath = get_calibrations_path(log=log, settings=settings)
    spectralFormatFile = calibrationRootPath + \
        "/" + dp["spectral format table"]

    # SPEC FORMAT TO PANDAS DATAFRAME
    from astropy.table import Table
    dat = Table.read(spectralFormatFile, format='fits')
    specFormatTable = dat.to_pandas()

    # EXTRACT REQUIRED PARAMETERS
    orderNums = specFormatTable["ORDER"].values
    waveLengthMin = specFormatTable["WLMINFUL"].values - 10
    waveLengthMax = specFormatTable["WLMAXFUL"].values + 10

    log.debug('completed the ``read_spectral_format`` function')
    return orderNums, waveLengthMin, waveLengthMax


def get_calibrations_path(
        log,
        settings):
    """*return the root path to the static calibrations*

    **Key Arguments:**

    - `log` -- logger
    - ``settings`` -- the settings dictionary

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import get_calibrations_path
    calibrationRootPath = get_calibrations_path(log=log, settings=settings)
    ```
    """
    log.debug('starting the ``get_calibrations_path`` function')

    # GENERATE PATH TO STATIC CALIBRATION DATA
    if "instrument" in settings:
        instrument = settings["instrument"]
    else:
        instrument = "soxs"
    calibrationRootPath = os.path.dirname(os.path.dirname(
        __file__)) + "/resources/static_calibrations/" + instrument

    log.debug('completed the ``get_calibrations_path`` function')
    return calibrationRootPath


def twoD_disp_map_image_to_dataframe(
        log,
        slit_length,
        twoDMapPath,
        kw=False,
        associatedFrame=False,
        removeMaskedPixels=False):
    """*convert the 2D dispersion image map to a pandas dataframe*

    **Key Arguments:**

    - `log` -- logger
    - `twoDMapPath` -- 2D dispersion map image path
    - `kw` -- fits keyword lookup dictionary
    - `associatedFrame` -- include a flux column in returned dataframe from a frame assosiated with the dispersion map. Default *False*
    - `removeMaskedPixels` -- remove the masked pixels from the assosicated image? Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe
    mapDF = twoD_disp_map_image_to_dataframe(log=log, twoDMapPath=twoDMap, associatedFrame=objectFrame, kw=kw)
    ```           
    """
    log.debug('starting the ``twoD_disp_map_image_to_dataframe`` function')

    import pandas as pd
    import numpy as np
    from astropy.io import fits

    # MAKE RELATIVE HOME PATH ABSOLUTE
    from os.path import expanduser
    home = expanduser("~")
    if twoDMapPath[0] == "~":
        twoDMapPath = twoDMapPath.replace("~", home)

    binx = 1
    biny = 1

    # FIND THE APPROPRIATE PREDICTED LINE-LIST
    if associatedFrame:
        arm = associatedFrame.header[kw("SEQ_ARM")]
        if arm != "NIR" and kw('WIN_BINX') in associatedFrame.header:
            binx = int(associatedFrame.header[kw('WIN_BINX')])
            biny = int(associatedFrame.header[kw('WIN_BINY')])

    hdul = fits.open(twoDMapPath)

    if binx > 1 or biny > 1:
        from astropy.nddata import block_reduce
        hdul["WAVELENGTH"].data = block_reduce(hdul["WAVELENGTH"].data, (biny, binx), func=np.mean)
        hdul["SLIT"].data = block_reduce(hdul["SLIT"].data, (biny, binx), func=np.mean)
        hdul["ORDER"].data = block_reduce(hdul["ORDER"].data, (biny, binx), func=np.mean)

    # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
    xdim = hdul[0].data.shape[1]
    ydim = hdul[0].data.shape[0]
    xarray = np.tile(np.arange(0, xdim), ydim)
    yarray = np.repeat(np.arange(0, ydim), xdim)

    thisDict = {
        "x": xarray,
        "y": yarray,
        "wavelength": hdul["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
        "slit_position": hdul["SLIT"].data.flatten().byteswap().newbyteorder(),
        "order": hdul["ORDER"].data.flatten().byteswap().newbyteorder(),
    }

    try:
        if associatedFrame:
            thisDict["flux"] = associatedFrame.data.flatten()
            thisDict["mask"] = associatedFrame.mask.flatten()
            thisDict["error"] = associatedFrame.uncertainty.flatten()

    except:
        if associatedFrame:
            thisDict["flux"] = associatedFrame.data.flatten().byteswap().newbyteorder()
            thisDict["mask"] = associatedFrame.mask.flatten().byteswap().newbyteorder()
            thisDict["error"] = associatedFrame.uncertainty.array.flatten().byteswap().newbyteorder()

    mapDF = pd.DataFrame.from_dict(thisDict)
    if removeMaskedPixels:
        mask = (mapDF["mask"] == False)
        mapDF = mapDF.loc[mask]

    # REMOVE ZEROS
    mask = (mapDF['wavelength'] == 0)
    mapDF = mapDF.loc[~mask]

    mapDF.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

    # REMOVE FILTERED ROWS FROM DATA FRAME
    mask = ((mapDF['slit_position'] < -slit_length / 2) | (mapDF['slit_position'] > slit_length / 2))
    mapDF = mapDF.loc[~mask]

    # SORT BY COLUMN NAME
    mapDF.sort_values(['wavelength'], inplace=True)

    log.debug('completed the ``twoD_disp_map_image_to_dataframe`` function')
    return mapDF


def predict_product_path(
        sofName):
    """*predict the path of the recipe product from a given SOF name*

    **Key Arguments:**

    - `log` -- logger,
    - `sofName` -- name or full path to the sof file

    **Usage:**

    ```python
    from soxspipe.commonutils import toolkit
    productPath = toolkit.predict_product_path(sofFilePath)
    ```           
    """
    try:
        sofName = os.path.basename(sofName)
    except:
        pass

    from soxspipe.commonutils import data_organiser
    from fundamentals.logs import emptyLogger
    log = emptyLogger()
    do = data_organiser(
        log=log,
        rootDir="."
    )
    currentSession, allSessions = do.session_list(silent=True)

    recipeName = sys.argv[1]
    if recipeName[0] == "-":
        recipeName = sys.argv[2]

    sofName = sofName.replace(".sof", "")
    if "_STARE_" in sofName:
        sofName += "_SKYSUB"
    productPath = f"./sessions/{currentSession}/product/soxs-" + recipeName.replace("_", "-").replace("sol", "solution").replace("centres", "centre").replace("spat", "spatial") + "/" + sofName + ".fits"
    productPath = productPath.replace("//", "/")

    return productPath


def add_recipe_logger(
        log,
        productPath):
    """*add a recipe-specific handler to the default logger that writes the recipe's logs adjacent to the recipe project*
    """
    import logging
    import os

    for handler in log.handlers:
        if handler.get_name() == "recipeLog":
            log.removeHandler(handler)
        if handler.get_name() == "recipeErr":
            log.removeHandler(handler)

    # GET THE EXTENSION (WITH DOT PREFIX)
    loggingPath = os.path.splitext(productPath)[0] + ".log"
    loggingErrorPath = os.path.splitext(productPath)[0] + "_ERROR.log"
    try:
        os.remove(loggingPath)
        os.remove(loggingErrorPath)
    except:
        pass

    # PARENT DIRECTORY PATH NEEDS TO EXIST FOR LOGGER TO WRITE
    parentDirectory = os.path.dirname(loggingPath)
    if not os.path.exists(parentDirectory):
        os.makedirs(parentDirectory)

    recipeLog = logging.FileHandler(loggingPath, mode='a', encoding=None, delay=False)
    recipeLogFormatter = logging.Formatter("%(message)s")
    recipeLog.set_name("recipeLog")
    recipeLog.setLevel(logging.INFO + 1)
    recipeLog.setFormatter(recipeLogFormatter)
    recipeLog.addFilter(MaxFilter(logging.WARNING))
    log.addHandler(recipeLog)

    recipeErr = logging.FileHandler(loggingErrorPath, mode='a', encoding=None, delay=True)
    recipeErrFormatter = logging.Formatter('%(asctime)s %(levelname)s: "%(pathname)s", line %(lineno)d, in %(funcName)s > %(message)s', '%Y-%m-%d %H:%M:%S')
    recipeErr.set_name("recipeErr")
    recipeErr.setLevel(logging.WARNING)
    recipeErr.setFormatter(recipeErrFormatter)
    log.addHandler(recipeErr)

    return log


class MaxFilter:
    def __init__(self, max_level):
        self.max_level = max_level

    def filter(self, record):
        if record.levelno < self.max_level:
            return True


def create_dispersion_solution_grid_lines_for_plot(
        log,
        dispMap,
        dispMapImage,
        associatedFrame,
        kw,
        skylines=False,
        slitPositions=False):
    """*give a dispersion solution and accompanying 2D dispersion map image, generate the grid lines to add to QC plots*

    **Key Arguments:**

    - `log` -- logger
    - ``dispMap`` -- path to dispersion map. Default *False*
    - ``dispMapImage`` -- the 2D dispersion map image
    - `associatedFrame` -- a frame associated with the reduction (to read arm and binning info).
    - `kw` -- fits header kw dictionary
    - `skylines` -- a list of skylines to use as the grid. Default *False*
    - `slitPositions` -- slit positions to plot (else plot min and max)

    **Usage:**

    ```python
    from soxspipe.commonutils.toolkit import create_dispersion_solution_grid_lines_for_plot
    gridLinePixelTable = create_dispersion_solution_grid_lines_for_plot(
        log=log,
        dispMap=dispMap,
        dispMapImage=dispMapImage,
        associatedFrame=CCDObject,
        kw=kw,
        skylines=skylines
    )

    for l in range(int(gridLinePixelTable['line'].max())):
        mask = (gridLinePixelTable['line'] == l)
        ax.plot(gridLinePixelTable.loc[mask]["fit_y"], gridLinePixelTable.loc[mask]["fit_x"], "w-", linewidth=0.5, alpha=0.8, color="black")
    ```           
    """
    log.debug('starting the ``create_dispersion_solution_grid_lines_for_plot`` function')

    import numpy as np
    import pandas as pd

    dispMapDF = twoD_disp_map_image_to_dataframe(log=log, slit_length=11, twoDMapPath=dispMapImage, associatedFrame=associatedFrame, kw=kw)
    uniqueOrders = dispMapDF['order'].unique()
    wlLims = []
    sPos = []

    for o in uniqueOrders:
        filDF = dispMapDF.loc[dispMapDF["order"] == o]
        wlLims.append((filDF['wavelength'].min() - 200, filDF['wavelength'].max() + 200))
        if isinstance(slitPositions, bool):
            sPos.append((filDF['slit_position'].min(), filDF['slit_position'].max()))
        else:
            sPos.append(slitPositions)

    lineNumber = 0
    for o, wlLim, spLim in zip(uniqueOrders, wlLims, sPos):
        wlRange = np.arange(wlLim[0], wlLim[1], 1)
        wlRange = np.append(wlRange, [wlLim[1]])
        for e in spLim:
            myDict = {
                "line": np.full_like(wlRange, lineNumber),
                "order": np.full_like(wlRange, o),
                "wavelength": wlRange,
                "slit_position": np.full_like(wlRange, e)
            }
            if lineNumber == 0:
                orderPixelTable = pd.DataFrame(myDict)
            else:
                orderPixelTableNew = pd.DataFrame(myDict)
                orderPixelTable = pd.concat([orderPixelTable, orderPixelTableNew], ignore_index=True)
            lineNumber += 1

        spRange = np.arange(spLim[0], spLim[-1], 1)
        spRange = np.append(spRange, [spLim[-1]])
        if skylines:
            mask = skylinesDF['WAVELENGTH'].between(wlLim[0], wlLim[1])
            wlRange = skylinesDF.loc[mask]['WAVELENGTH'].values
        else:
            step = int(wlLim[1] - wlLim[0]) / 400
            wlRange = np.arange(wlLim[0], wlLim[1], step)
        wlRange = np.append(wlRange, [wlLim[1]])
        for l in wlRange:
            myDict = {
                "line": np.full_like(spRange, lineNumber),
                "order": np.full_like(spRange, o),
                "wavelength": np.full_like(spRange, l),
                "slit_position": spRange
            }
            orderPixelTableNew = pd.DataFrame(myDict)
            orderPixelTable = pd.concat([orderPixelTable, orderPixelTableNew], ignore_index=True)
            lineNumber += 1

    orderPixelTable = dispersion_map_to_pixel_arrays(
        log=log,
        dispersionMapPath=dispMap,
        orderPixelTable=orderPixelTable
    )

    log.debug('completed the ``create_dispersion_solution_grid_lines_for_plot`` function')
    return orderPixelTable

# use the tab-trigger below for new function
# xt-def-function
