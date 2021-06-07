#!/usr/bin/env python
# encoding: utf-8
"""
*Fit and subtract background flux from scattered light from frame*

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


class subtract_background(object):
    """
    *The worker class for the subtract_background module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``frame`` -- the frame to subtract background light from
        - ``orderTable`` -- the order geometry table
        - ``orderExt`` -- the order geometry table of inner-order masked along the y-axis

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a subtract_background object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``subtract_background`` to documentation
        - create a blog post about what ``subtract_background`` does
    ```

    ```python
    usage code 
    ```

    """
    # Initialisation

    def __init__(
            self,
            log,
            frame,
            orderTable,
            orderExt=0.0,
            settings=False
    ):
        self.log = log
        log.debug("instansiating a new 'subtract_background' object")
        self.settings = settings
        self.frame = frame
        self.orderTable = orderTable
        self.orderExt = orderExt

        return None

    def get(self):
        """
        *get the subtract_background object*

        **Return:**
            - ``subtract_background``

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code 
        ```
        """
        self.log.debug('starting the ``get`` method')

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=self.orderExt)

        self.mask_order_locations(orderPixelTable)

        # SUPPRESS MATPLOTLIB WARNINGS
        # import warnings
        # warnings.filterwarnings("ignore")
        # import matplotlib.pyplot as plt
        # fig, (ax1) = plt.subplots(1, 1, figsize=(30, 15))
        # plt.imshow(self.frame.mask, 'gray', interpolation='none')
        # # plt.imshow(maskImage, 'Reds_r', interpolation='none', alpha=0.3)
        # plt.show()
        # plt.imshow(self.frame, 'gray', interpolation='none')
        # plt.show()

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=self.frame, show=True, ext=None)

        subtract_background = None

        self.log.debug('completed the ``get`` method')
        return subtract_background

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

    # use the tab-trigger below for new method
    # xt-class-method
