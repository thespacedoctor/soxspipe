#!/usr/bin/env python
# encoding: utf-8
"""
*fit and subtract background flux from scattered light from frame*

:Author:
    David Young

:Date Created:
    June  3, 2021
"""
from soxspipe.commonutils import keyword_lookup
from os.path import expanduser
import random
from scipy.ndimage.filters import median_filter
from scipy.signal import medfilt2d
from soxspipe.commonutils.toolkit import quicklook_image
from scipy.interpolate import BSpline, splrep
from astropy import units as u
from astropy.nddata import CCDData
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
from soxspipe.commonutils.toolkit import unpack_order_table
from fundamentals import tools
from builtins import object
import sys
import math
import pandas as pd
import os
os.environ['TERM'] = 'vt100'


class subtract_background(object):
    """
    *fit and subtract background flux from scattered light from frame*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``frame`` -- the frame to subtract background light from
        - ``orderTable`` -- the order geometry table
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``surfaceFit`` -- fit a bspline surface to the background pixel data rather than series of linear detector-column bspline fits

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To fit and subtract the background from an image:

    ```python
    from soxspipe.commonutils import subtract_background
    background = subtract_background(
        log=log,
        frame=myCCDDataObject,
        orderTable="/path/to/orderTable",
        settings=settings,
        surfaceFit=True
    )
    backgroundFrame, backgroundSubtractedFrame = background.subtract()
    ```

    """
    # Initialisation

    def __init__(
            self,
            log,
            frame,
            orderTable,
            settings=False,
            qcTable=False,
            productsTable=False,
            surfaceFit=False
    ):
        self.log = log
        log.debug("instantiating a new 'subtract_background' object")
        self.settings = settings
        self.frame = frame
        self.orderTable = orderTable
        self.qc = qcTable
        self.products = productsTable
        self.surfaceFit = surfaceFit

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = frame.header[kw("SEQ_ARM")]
        self.dateObs = frame.header[kw("DATE_OBS")]

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext='data', stdWindow=3, surfacePlot=True, title="Initial input frame needing scattered light subtraction")

        return None

    def subtract(self):
        """
        *fit and subtract background light from frame*

        **Return:**
            - ``backgroundSubtractedFrame`` -- a CCDData object of the original input frame with fitted background light subtracted
        """
        self.log.debug('starting the ``subtract`` method')

        kw = self.kw
        imageType = self.frame.header[kw("DPR_TYPE").lower()].replace(",", "-")
        imageTech = self.frame.header[kw("DPR_TECH").lower()].replace(",", "-")
        imageCat = self.frame.header[kw("DPR_CATG").lower()].replace(",", "-")

        print(f"\n# FITTING AND SUBTRACTING SCATTERED LIGHT BACKGROUND FROM {self.arm} {imageCat} {imageTech} {imageType} FRAME")

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=1.5)

        originalMask = np.copy(self.frame.mask)

        # MASK THE INNER ORDER AREA (AND BAD PIXELS)
        self.mask_order_locations(orderPixelTable)

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext=None, surfacePlot=True, title="Initial input frame with order pixels masked")

        backgroundFrame = self.create_background_image(
            rowFitOrder=self.settings['background-subtraction']['bsline-deg'], medianFilterSize=self.settings['background-subtraction']['median-filter-pixels'])

        # REPLACE MASK
        self.frame.mask = originalMask

        backgroundSubtractedFrame = self.frame.subtract(backgroundFrame)
        backgroundSubtractedFrame.header = self.frame.header

        quicklook_image(
            log=self.log, CCDObject=backgroundSubtractedFrame, show=False, ext='data')

        self.log.debug('completed the ``subtract`` method')
        return backgroundFrame, backgroundSubtractedFrame

    def mask_order_locations(
            self,
            orderPixelTable):
        """*mask the order locations and return the masked frame*

        **Key Arguments:**
            - ``orderPixelTable`` -- the order location in a pandas datafrmae.
        """
        self.log.debug('starting the ``mask_order_locations`` method')

        # MASK DATA INSIDE OF ORDERS (EXPAND THE INNER-ORDER AREA IF NEEDED)
        uniqueOrders = orderPixelTable['order'].unique()
        expandEdges = 2
        for o in uniqueOrders:
            ycoord = orderPixelTable.loc[
                (orderPixelTable["order"] == o)]["ycoord"]
            xcoord_edgeup = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                "xcoord_edgeup"] + expandEdges
            xcoord_edgelow = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                "xcoord_edgelow"] - expandEdges
            xcoord_edgelow, xcoord_edgeup, ycoord = zip(*[(x1, x2, y) for x1, x2, y in zip(xcoord_edgelow, xcoord_edgeup, ycoord) if x1 > 0 and x1 < self.frame.data.shape[
                                                        1] and x2 > 0 and x2 < self.frame.data.shape[1] and y > 0 and y < self.frame.data.shape[0]])
            for y, u, l in zip(ycoord, np.ceil(xcoord_edgeup).astype(int), np.floor(xcoord_edgelow).astype(int)):
                self.frame.mask[y, l:u] = 1

        self.log.debug('completed the ``mask_order_locations`` method')
        return None

    def create_background_image(
            self,
            rowFitOrder,
            medianFilterSize):
        """*model the background image from intra-order flux detected*

        **Key Arguments:**
            - ``rowFitOrder`` -- order of the polynomial fit to flux in each row
            - ``medianFilterSize`` -- the length of the line to median filter in y-direction (after bspline fitting)
        """
        self.log.debug('starting the ``create_background_image`` method')

        from astropy.stats import sigma_clip, mad_std
        maskedImage = np.ma.array(self.frame.data, mask=self.frame.mask)
        # SIGMA-CLIP THE DATA
        clippedDataMask = sigma_clip(
            maskedImage, sigma_lower=2, sigma_upper=5, maxiters=3, cenfunc='median', stdfunc=mad_std)
        # COMBINE MASK WITH THE BAD PIXEL MASK
        mask = (clippedDataMask.mask == 1) | (self.frame.mask == 1)
        maskedImage = np.ma.array(self.frame.data, mask=mask)

        quicklook_image(
            log=self.log, CCDObject=maskedImage, show=False, ext=True, stdWindow=3, surfacePlot=True, title="Sigma clipped masked image")

        # PLACEHOLDER ARRAY FOR BACKGROUND IMAGE
        backgroundMap = np.zeros_like(self.frame)

        if self.surfaceFit:

            fig = plt.figure(figsize=(40, 10))
            from scipy.interpolate import griddata
            flatMask = mask.flatten()
            X, Y = np.meshgrid(np.linspace(0, self.frame.data.shape[1], self.frame.data.shape[1]), np.linspace(0, self.frame.data.shape[0], self.frame.data.shape[0]))
            Z = self.frame.data.flatten()
            Xflat = X.flatten()
            Yflat = Y.flatten()
            Xdata = Xflat[~flatMask]
            Ydata = Yflat[~flatMask]
            Zdata = Z[~flatMask]

            backgroundMap = griddata((Xdata, Ydata), Zdata, (X, Y), method='cubic')

            quicklook_image(
                log=self.log, CCDObject=backgroundMap, show=False, ext=None, surfacePlot=True, title="Scattered light background image - fitted with Cubic")

        else:
            for idx, row in enumerate(maskedImage):

                # SET X TO A MASKED RANGE
                xunmasked = ma.masked_array(np.linspace(
                    0, len(row), len(row), dtype=int), mask=row.mask)

                # fit = ma.polyfit(xunmasked, row, deg=rowFitOrder)
                xfit = np.linspace(0, len(row), len(row), dtype=int)

                # yfit = np.polyval(fit, xfit)

                xmasked = xunmasked[~xunmasked.mask]
                xmin = xmasked.min()
                xmax = xmasked.max()
                rowmasked = row[~row.mask].byteswap().newbyteorder()

                window = 9
                hw = math.floor(window / 2)
                # rowmaskedSmoothed = pd.Series(rowmasked).rolling(window=window, center=True).quantile(.1)
                rowmaskedSmoothed = pd.Series(rowmasked).rolling(window=window, center=True).median()
                rowmaskedSmoothed[:hw] = rowmaskedSmoothed.iloc[hw + 1]
                rowmaskedSmoothed[-hw:] = rowmaskedSmoothed.iloc[-hw - 1]
                rowmasked[:hw] = rowmasked[hw + 1]
                rowmasked[-hw:] = rowmasked[-hw - 1]

                seedKnots = xmasked[1:-1:25]
                t, c, k = splrep(xmasked, rowmaskedSmoothed, t=seedKnots, s=10.0, k=rowFitOrder)

                spline = BSpline(t, c, k, extrapolate=True)
                yfit = spline(xfit)

                # ADD FITTED ROW TO BACKGROUND IMAGE
                backgroundMap[idx, :] = yfit

                if random.randint(1, 501) == 42 and 1 == 0:
                    fig, (ax1) = plt.subplots(1, 1, figsize=(30, 15))
                    plt.scatter(xmasked, rowmasked)
                    plt.scatter(xmasked, rowmaskedSmoothed)
                    plt.plot(xfit, yfit, 'b')
                    plt.ylabel("flux")
                    plt.xlabel("pixels along row")
                    plt.show()

            quicklook_image(
                log=self.log, CCDObject=backgroundMap, show=False, ext=None, surfacePlot=True, title="Scattered light background image")

        backgroundMap = medfilt2d(
            backgroundMap, medianFilterSize)

        quicklook_image(
            log=self.log, CCDObject=backgroundMap, show=False, ext=None, surfacePlot=True, title="Scattered light background image with median filtering")

        backgroundMap = CCDData(backgroundMap, unit=self.frame.unit)

        self.log.debug('completed the ``create_background_image`` method')
        return backgroundMap
