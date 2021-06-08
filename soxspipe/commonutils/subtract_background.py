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
from astropy.convolution import Gaussian2DKernel, Box2DKernel
from astropy.convolution import convolve
from scipy.ndimage.filters import median_filter
from astropy.nddata import CCDData
from astropy import units as u


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
    backgroundSubtractedFrame = background.subtracted()
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

        backgroundFrame = self.create_background_image(
            rowFitOrder=self.settings['background-subtraction']['poly-deg'], convolutionRadius=self.settings['background-subtraction']['smoothing-radius-pixels'])

        backgroundSubtractedFrame = self.frame.subtract(backgroundFrame)

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=backgroundSubtractedFrame, show=False)

        self.log.debug('completed the ``subtract`` method')
        return backgroundSubtractedFrame

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
            convolutionRadius):
        """*model the background image from intra-order flux detected*

        **Key Arguments:**
            - ``rowFitOrder`` -- order of the polynomial fit to flux in each row
            - ``convolutionRadius`` -- radius of the kernel used to convolve the fitted background
        """
        self.log.debug('starting the ``create_background_image`` method')

        maskedImage = np.ma.array(self.frame.data, mask=self.frame.mask)

        # PLACEHOLDER ARRAY FOR BACKGROUND IMAGE
        backgroundMap = np.zeros_like(self.frame)
        for idx, row in enumerate(maskedImage):

            # SET X TO A MASKED RANGE
            xunmasked = ma.masked_array(np.linspace(
                0, len(row), len(row), dtype=int), mask=row.mask)
            fit = ma.polyfit(xunmasked, row, deg=rowFitOrder)
            xfit = np.linspace(0, len(row), len(row), dtype=int)
            yfit = np.polyval(fit, xfit)

            # ADD FITTED ROW TO BACKGROUND IMAGE
            backgroundMap[idx, :] = yfit

            if idx == 1000 and 1 == 0:
                fig, (ax1) = plt.subplots(1, 1, figsize=(30, 15))
                plt.scatter(xunmasked, row)
                plt.plot(xfit, yfit, 'b')
                plt.ylabel("flux")
                plt.xlabel("pixels along row")
                plt.show()

        kernel = Gaussian2DKernel(convolutionRadius)
        backgroundMap = convolve(backgroundMap, kernel)

        backgroundFrame = CCDData(backgroundFrame, unit=u.electron)

        self.log.debug('completed the ``create_background_image`` method')
        return backgroundMap
