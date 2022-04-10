#!/usr/bin/env python
# encoding: utf-8
"""
*definition of polynomial functions needed throughout code*

:Author:
    David Young

:Date Created:
    September 10, 2020
"""
################# GLOBAL IMPORTS ####################
import numpy as np
import math
from fundamentals import tools
from builtins import object
import sys
import os
import pandas as pd
os.environ['TERM'] = 'vt100'


class chebyshev_order_wavelength_polynomials():
    """*the chebyshev polynomial fits for the single pinhole frames; to be iteratively fitted to minimise errors*

    **Key Arguments:**
        - ``log`` -- logger
        - ``order_deg`` -- degree of the order polynomial components
        - ``wavelength_deg`` -- degree of wavelength polynomial components
        - ``slit_deg`` -- degree of the slit polynomial components
        - ``exponents_included`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
    poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg).poly
    ```
    """

    def __init__(
            self,
            log,
            order_deg,
            wavelength_deg,
            slit_deg,
            exponents_included=False
    ):
        self.log = log
        self.order_deg = order_deg
        self.wavelength_deg = wavelength_deg
        self.slit_deg = slit_deg
        self.exponents_included = exponents_included

        return None

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**
        - ``orderPixelTable`` -- a pandas dataframe containing wavelengths, orders and slit positions
        - ``*coeff`` -- a list of the initial coefficients

        **Return:**
        - ``lhsVals`` -- the left-hand-side vals of the fitted polynomials
        """
        self.log.debug('starting the ``poly`` method')

        # UNPACK TUPLE INPUT
        order_deg = self.order_deg
        wavelength_deg = self.wavelength_deg
        slit_deg = self.slit_deg

        n_coeff = 0
        lhsVals = np.zeros(len(orderPixelTable.index))

        # FOR LOOPS ARE THE RIGHT TOOL TO PERFORM COMPUTATIONS OR RUN FUNCTIONS. LIST COMPREHENSION IS SLOW IN THESE CASES

        if self.exponents_included == False:
            orderVals = orderPixelTable["order"].values
            wlVals = orderPixelTable["wavelength"].values
            spVals = orderPixelTable["slit_position"].values

            for i in range(0, order_deg + 1):
                for j in range(0, wavelength_deg + 1):
                    for k in range(0, slit_deg + 1):
                        lhsVals += coeff[n_coeff] * orderVals**i * \
                            wlVals**j * \
                            spVals**k
                        n_coeff += 1
        else:
            for i in range(0, order_deg + 1):
                for j in range(0, wavelength_deg + 1):
                    for k in range(0, slit_deg + 1):
                        lhsVals += coeff[n_coeff] * orderPixelTable[f"order_pow_{i}"].values * orderPixelTable[f"wavelength_pow_{j}"].values * orderPixelTable[f"slit_position_pow_{k}"].values
                        n_coeff += 1

        self.log.debug('completed the ``poly`` method')

        return lhsVals


class chebyshev_xy_polynomial():
    """*the chebyshev polynomial fits for the pinhole flat frame order tracing; to be iteratively fitted to minimise errors*

    **Key Arguments:**
        - ``log`` -- logger
        - ``yCol`` -- name of the yCol
        - ``y_deg`` -- y degree of the polynomial components
        - ``exponents_included`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
    poly = chebyshev_xy_polynomial(
            log=self.log, deg=deg).poly
    ```
    """

    def __init__(
            self,
            log,
            y_deg,
            yCol=False,
            exponents_included=False
    ):
        self.log = log
        self.y_deg = y_deg
        self.yCol = yCol
        self.exponents_included = exponents_included

        return None

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**
        - ``orderPixelTable`` -- data frame with all pixel data arrays
        - ``*coeff`` -- a list of the initial coefficients

        **Return:**
        - ``xvals`` -- the x values of the fitted polynomial
        """
        self.log.info('starting the ``poly`` method')

        n_coeff = 0
        if not isinstance(orderPixelTable, pd.core.frame.DataFrame):
            yarray = orderPixelTable
            lhsVals = np.zeros(len(orderPixelTable))
        else:
            yarray = orderPixelTable[self.yCol].values
            lhsVals = np.zeros(len(orderPixelTable.index))

        if not self.exponents_included:
            # POLYNOMIALS SUMS
            for i in range(0, self.y_deg + 1):
                lhsVals += coeff[n_coeff] * yarray**i
                n_coeff += 1
        else:
            for i in range(0, self.y_deg + 1):
                lhsVals += coeff[n_coeff] * orderPixelTable[f"y_pow_{i}"].values
                n_coeff += 1

        self.log.info('completed the ``poly`` method')

        return lhsVals


class chebyshev_order_xy_polynomials():
    """*the chebyshev polynomial fits FIX ME*

    **Key Arguments:**
        - ``log`` -- logger
        - ``order_deg`` -- degree of the order polynomial components
        - ``y_deg`` -- degree of y polynomial components
        - ``yCol`` -- name of the y column (if needed). Default *False*
        - ``orderCol`` -- name of the order column (if needed). Default *False*
        - ``exponents_included`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
    poly = chebyshev_order_wavelength_polynomials(
            log=self.log, order_deg=order_deg, wavelength_deg=wavelength_deg, slit_deg=slit_deg).poly
    ```
    """

    def __init__(
            self,
            log,
            order_deg,
            y_deg,
            yCol=False,
            orderCol=False,
            exponents_included=False
    ):
        self.log = log
        self.order_deg = order_deg
        self.y_deg = y_deg
        self.yCol = yCol
        self.orderCol = orderCol
        self.exponents_included = exponents_included

        return None

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**
        - ``orderPixelTable`` -- a pandas dataframe containing x, y, order
        - ``*coeff`` -- a list of the initial coefficients

        **Return:**
        - ``lhsVals`` -- the left-hand-side vals of the fitted polynomials
        """
        self.log.debug('starting the ``poly`` method')

        # UNPACK TUPLE INPUT
        order_deg = self.order_deg
        y_deg = self.y_deg

        n_coeff = 0
        lhsVals = np.zeros(len(orderPixelTable.index))

        # FOR LOOPS ARE THE RIGHT TOOL TO PERFORM COMPUTATIONS OR RUN FUNCTIONS. LIST COMPREHENSION IS SLOW IN THESE CASES

        if self.exponents_included == False:
            orderVals = orderPixelTable[self.orderCol].values
            yVals = orderPixelTable[self.yCol].values

            for i in range(0, order_deg + 1):
                for j in range(0, y_deg + 1):
                    lhsVals += coeff[n_coeff] * orderVals**i * \
                        yVals**j
                    n_coeff += 1
        else:
            for i in range(0, order_deg + 1):
                for j in range(0, y_deg + 1):
                    lhsVals += coeff[n_coeff] * orderPixelTable[f"order_pow_{i}"].values * orderPixelTable[f"y_pow_{j}"].values
                    n_coeff += 1

        self.log.debug('completed the ``poly`` method')

        return lhsVals
