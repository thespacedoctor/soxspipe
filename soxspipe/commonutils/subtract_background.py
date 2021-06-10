#!/usr/bin/env python
# encoding: utf-8
"""
*fit and subtract background flux from scattered light from frame*

:Author:
    David Young

:Date Created:
    June  3, 2021
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
from soxspipe.commonutils.toolkit import unpack_order_table
import numpy.ma as ma
import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
from astropy import units as u
from scipy.interpolate import BSpline, splrep
from soxspipe.commonutils.toolkit import quicklook_image
from scipy.ndimage.filters import median_filter
import random
from os.path import expanduser


class subtract_background(object):
    """
    *fit and subtract background flux from scattered light from frame*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``frame`` -- the frame to subtract background light from
        - ``orderTable`` -- the order geometry table

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To fit and subtract the background from an image:

    ```python
    from soxspipe.commonutils import subtract_background
    background = subtract_background(
        log=log,
        frame=myCCDDataObject,
        orderTable="/path/to/orderTable",
        settings=settings
    )
    backgroundSubtractedFrame = background.subtract()
    ```

    """
    # Initialisation

    def __init__(
            self,
            log,
            frame,
            orderTable,
            settings=False
    ):
        self.log = log
        log.debug("instantiating a new 'subtract_background' object")
        self.settings = settings
        self.frame = frame
        self.orderTable = orderTable

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext='data')

        return None

    def subtract(self):
        """
        *fit and subtract background light from frame*

        **Return:**
            - ``backgroundSubtractedFrame`` -- a CCDData object of the original input frame with fitted background light subtracted
        """
        self.log.debug('starting the ``subtract`` method')

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=self.settings['background-subtraction']['order-extension-fraction-for-background-subtraction'])

        # MASK THE INNER ORDER AREA (AND BAD PIXELS)
        self.mask_order_locations(orderPixelTable)

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext=None)

        backgroundFrame = self.create_background_image(
            rowFitOrder=self.settings['background-subtraction']['bsline-deg'], medianFilterSize=self.settings['background-subtraction']['median-filter-pixels'])

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
        expandEdges = 3
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

        maskedImage = np.ma.array(self.frame.data, mask=self.frame.mask)

        # PLACEHOLDER ARRAY FOR BACKGROUND IMAGE
        backgroundMap = np.zeros_like(self.frame)
        for idx, row in enumerate(maskedImage):

            # SET X TO A MASKED RANGE
            xunmasked = ma.masked_array(np.linspace(
                0, len(row), len(row), dtype=int), mask=row.mask)
            # fit = ma.polyfit(xunmasked, row, deg=rowFitOrder)
            xfit = np.linspace(0, len(row), len(row), dtype=int)
            # yfit = np.polyval(fit, xfit)

            t, c, k = splrep(xunmasked[~xunmasked.mask], row[
                             ~row.mask], s=0.0006, k=rowFitOrder)
            spline = BSpline(t, c, k, extrapolate=False)
            yfit = spline(xfit)

            # ADD FITTED ROW TO BACKGROUND IMAGE
            backgroundMap[idx, :] = yfit

            if random.randint(1, 101) == 42 and 1 == 0:
                print(idx)
                fig, (ax1) = plt.subplots(1, 1, figsize=(30, 15))
                plt.scatter(xunmasked, row)
                plt.plot(xfit, yfit, 'b')
                plt.ylabel("flux")
                plt.xlabel("pixels along row")
                plt.show()

        quicklook_image(
            log=self.log, CCDObject=backgroundMap, show=False, ext=None)

        backgroundMap = median_filter(
            backgroundMap, size=(medianFilterSize, medianFilterSize))

        quicklook_image(
            log=self.log, CCDObject=backgroundMap, show=False, ext=None)

        backgroundMap = CCDData(backgroundMap, unit=u.electron)

        self.log.debug('completed the ``create_background_image`` method')
        return backgroundMap
