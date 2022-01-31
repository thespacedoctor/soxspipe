#!/usr/bin/env python
# encoding: utf-8
"""
*small reusable functions used throughout soxspipe*

:Author:
    David Young

:Date Created:
    September 18, 2020
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
import unicodecsv as csv
import numpy as np
import matplotlib.pyplot as plt
import random
import numpy.ma as ma
import pandas as pd
from soxspipe.commonutils import detector_lookup
import math
from soxspipe.commonutils import keyword_lookup
from datetime import datetime
from copy import copy
from os.path import expanduser


def cut_image_slice(
        log,
        frame,
        width,
        length,
        x,
        y,
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

    halfSlice = length / 2
    # NEED AN EVEN PIXEL SIZE
    if (width % 2) != 0:
        halfwidth = (width - 1) / 2
    else:
        halfwidth = width / 2

    # CHECK WE ARE NOT GOING BEYOND BOUNDS OF FRAME
    if (x > frame.shape[1] - halfSlice) or (y > frame.shape[0] - halfwidth) or (x < halfSlice) or (y < halfwidth):
        return None

    slice = frame[int(y - halfwidth):int(y + halfwidth + 1),
                  int(x - halfSlice):int(x + halfSlice)]

    if median:
        slice = ma.median(slice, axis=0)

    if plot and random.randint(1, 101) < 10000:
        # CHECK THE SLICE POINTS IF NEEDED
        x = np.arange(0, len(slice))
        plt.figure(figsize=(8, 5))
        plt.plot(x, slice, 'ko')
        plt.xlabel('Position')
        plt.ylabel('Flux')
        plt.show()

    log.debug('completed the ``cut_image_slice`` function')
    return slice


def quicklook_image(
        log,
        CCDObject,
        show=True,
        ext="data",
        stdWindow=3):
    """*generate a quicklook image of a CCDObject - useful for development/debugging*

    **Key Arguments:**

    - ``log`` -- logger
    - ``CCDObject`` -- the CCDObject to plot
    - ``show`` -- show the image. Set to False to skip
    - ``ext`` -- the name of the the extension to show. Can be "data", "mask" or "err". Default "data". 

    ```python
    from soxspipe.commonutils.toolkit import quicklook_image
    quicklook_image(
        log=self.log, CCDObject=myframe, show=True)
    ```           
    """
    log.debug('starting the ``quicklook_image`` function')

    if not show:
        return

    if ext == "data":
        frame = CCDObject.data
    elif ext == "mask":
        frame = CCDObject.mask
    elif ext == "err":
        frame = CCDObject.err
    else:
        # ASSUME ONLY NDARRAY
        frame = CCDObject

    rotatedImg = np.rot90(frame, 1)
    rotatedImg = np.flipud(np.rot90(frame, 1))

    std = np.nanstd(frame)
    mean = np.nanmean(frame)
    vmax = mean + stdWindow * std
    vmin = mean - stdWindow * std
    palette = copy(plt.cm.gray)
    # palette.set_over('r', 1.0)
    # palette.set_under('g', 1.0)
    palette.set_bad("#dc322f", 1.0)
    plt.figure(figsize=(12, 5))
    plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
               cmap=palette, alpha=1, aspect='auto')
    if mean > 10:
        fmt = '%1.0f'
        cbar = plt.colorbar(format=fmt)
    else:
        plt.colorbar()
    # cbar.ticklabel_format(useOffset=False)
    plt.xlabel(
        "y-axis", fontsize=10)
    plt.ylabel(
        "x-axis", fontsize=10)
    plt.gca().invert_yaxis()
    plt.show()

    log.debug('completed the ``quicklook_image`` function')
    return None


def unpack_order_table(
        log,
        orderTablePath,
        extend=0.):
    """*unpack an order table and return a top-level `orderPolyTable` data-frame and a second `orderPixelTable` data-frame with the central-trace coordinates of each order given

    **Key Arguments:**

    - ``orderTablePath`` -- path to the order table
    - ``extend`` -- fractional increase to the order area in the y-axis (needed for masking)

    **Usage:**

    ```python
    # UNPACK THE ORDER TABLE
    from soxspipe.commonutils.toolkit import unpack_order_table
    orderPolyTable, orderPixelTable = unpack_order_table(
        log=self.log, orderTablePath=orderTablePath, extend=0.)
    ```           
    """
    log.debug('starting the ``functionName`` function')

    # READ CSV FILE TO PANDAS DATAFRAME
    orderPolyTable = pd.read_csv(orderTablePath, index_col=False,
                                 na_values=['NA', 'MISSING'])

    # ADD Y-COORD LIST
    ycoords = [np.arange(math.floor(l) - int(r * extend), math.ceil(u) + int(r * extend), 1) for l, u, r in zip(
        orderPolyTable["ymin"].values, orderPolyTable["ymax"].values, orderPolyTable["ymax"].values - orderPolyTable["ymin"].values)]
    orders = [np.full_like(a, o) for a, o in zip(
        ycoords, orderPolyTable["order"].values)]

    # RUN COORDINATES THROUGH POLYNOMIALS TO GET X-COORDS
    xcoords_centre = []
    xcoords_edgeup = []
    xcoords_edgelow = []
    for index, row in orderPolyTable.iterrows():
        cent_coeff = [float(v) for k, v in row.items() if "cent_" in k]
        poly = chebyshev_xy_polynomial(log=log, deg=int(row["degy_cent"])).poly
        xcoords_centre.append(np.array(poly(ycoords[index], *cent_coeff)))
        if "degy_edgeup" in row:
            degy_edge = int(row["degy_edgeup"])
            edgeup_coeff = [float(v)
                            for k, v in row.items() if "edgeup_" in k]
            poly = chebyshev_xy_polynomial(
                log=log, deg=degy_edge).poly
            xcoords_edgeup.append(
                np.array(poly(ycoords[index], *edgeup_coeff)))

            edgelow_coeff = [float(v)
                             for k, v in row.items() if "edgelow_" in k]
            poly = chebyshev_xy_polynomial(
                log=log, deg=degy_edge).poly
            xcoords_edgelow.append(
                np.array(poly(ycoords[index], *edgelow_coeff)))

    # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
    myDict = {
        "order": np.concatenate(orders),
        "ycoord": np.concatenate(ycoords),
        'xcoord_centre': np.concatenate(xcoords_centre)
    }
    # ADD ODER EDGES IF NEEDED
    if len(xcoords_edgeup):
        myDict['xcoord_edgeup'] = np.concatenate(xcoords_edgeup)
        myDict['xcoord_edgelow'] = np.concatenate(xcoords_edgelow)

    # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
    orderPixelTable = pd.DataFrame(myDict)

    log.debug('completed the ``functionName`` function')
    return orderPolyTable, orderPixelTable


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

    qcTable = qcTable.append({
        "soxspipe_recipe": recipeName,
        "qc_name": "NaN count",
        "qc_value": nanCount,
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow
    }, ignore_index=True)

    # COUNT BAD-PIXELS
    badCount = frame.mask.sum()
    totalPixels = np.size(frame.mask)
    percent = (float(badCount) / float(totalPixels)) * 100.

    qcTable = qcTable.append({
        "soxspipe_recipe": recipeName,
        "qc_name": "bad-pixel count",
        "qc_value": badCount,
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow
    }, ignore_index=True)

    qcTable = qcTable.append({
        "soxspipe_recipe": recipeName,
        "qc_name": "bad-pixel percentage",
        "qc_value": percent,
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow
    }, ignore_index=True)

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

    # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
    # FOLDER
    kw = keyword_lookup(
        log=log,
        settings=settings
    ).get
    kw = kw
    arm = frame.header[kw("SEQ_ARM")]
    dateObs = frame.header[kw("DATE_OBS")]

    # UNPACK THE ORDER TABLE
    orderTableMeta, orderTablePixels = unpack_order_table(
        log=log, orderTablePath=orderTablePath)

    mask = np.ones_like(frame.data)

    xcoords_up = orderTablePixels["xcoord_edgeup"].values
    xcoords_low = orderTablePixels["xcoord_edgelow"].values
    ycoords = orderTablePixels["ycoord"].values
    xcoords_up = xcoords_up.astype(int)
    xcoords_low = xcoords_low.astype(int)

    # UPDATE THE MASK
    for u, l, y in zip(xcoords_up, xcoords_low, ycoords):
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

    qcTable = qcTable.append({
        "soxspipe_recipe": recipeName,
        "qc_name": "Mean inner-order pixel value",
        "qc_value": mean,
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow
    }, ignore_index=True)

    qcTable = qcTable.append({
        "soxspipe_recipe": recipeName,
        "qc_name": "Sum of all inner-order pixel values",
        "qc_value": flux,
        "qc_unit": "",
        "obs_date_utc": dateObs,
        "reduction_date_utc": utcnow
    }, ignore_index=True)

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
    calibrationRootPath = settings[
        "calibration-data-root"].replace("~", home)
    spectralFormatFile = calibrationRootPath + \
        "/" + dp["spectral format table"]

    # READ CSV FILE TO PANDAS DATAFRAME
    specFormatTable = pd.read_csv(
        spectralFormatFile, index_col=False, na_values=['NA', 'MISSING'])
    # # DATAFRAME INFO

    # EXTRACT REQUIRED PARAMETERS
    orderNums = specFormatTable["ORDER"].values
    waveLengthMin = specFormatTable["WLMINFUL"].values
    waveLengthMax = specFormatTable["WLMAXFUL"].values

    log.debug('completed the ``read_spectral_format`` function')
    return orderNums, waveLengthMin, waveLengthMax
