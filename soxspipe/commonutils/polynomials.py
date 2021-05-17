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
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools
import math
import numpy as np


class chebyshev_order_wavelength_polynomials():
    """*the chebyshev polynomial fits for the single pinhole frames; to be iteratively fitted to minimise errors*

    **Key Arguments:**
        - ``log`` -- logger
        - ``order_deg`` -- degree of the order polynomial components
        - ``wavelength_deg`` -- degree of wavelength polynomial components
        - ``slit_deg`` -- degree of the slit polynomial components

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
            slit_deg
    ):
        self.log = log
        self.order_deg = order_deg
        self.wavelength_deg = wavelength_deg
        self.slit_deg = slit_deg

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

        # lhsVals = np.sum([v * c for v, c in zip([orderPixelTable["order"].values**i * orderPixelTable["wavelength"].values**j * orderPixelTable["slit_position"].values**k for i in range(0, order_deg + 1)
        # for j in range(0, wavelength_deg + 1) for k in range(0, slit_deg +
        # 1)], coeff)], axis=0)

        # THE LIST COMPREHENSION ABOVE DID NOT SPEED UP THE NEST LOOPS BELOW -
        # KEEPING LOOPS!
        n_coeff = 0
        lhsVals = np.zeros(len(orderPixelTable.index))
        for i in range(0, order_deg + 1):
            for j in range(0, wavelength_deg + 1):
                for k in range(0, slit_deg + 1):
                    lhsVals += coeff[n_coeff] * orderPixelTable["order"].values**i * \
                        orderPixelTable["wavelength"].values**j * \
                        orderPixelTable["slit_position"].values**k
                    n_coeff += 1

        self.log.debug('completed the ``poly`` method')

        return lhsVals


class chebyshev_xy_polynomial():
    """*the chebyshev polynomial fits for the pinhole flat frame order tracing; to be iteratively fitted to minimise errors*

    **Key Arguments:**
        - ``log`` -- logger
        - ``deg`` -- degree of the polynomial components

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
            deg
    ):
        self.log = log
        self.deg = deg

        return None

    def poly(self, yarray, *coeff):
        """the polynomial definition

        **Key Arguments:**
        - ``yarray`` -- the y coordinates
        - ``*coeff`` -- a list of the initial coefficients

        **Return:**
        - ``xvals`` -- the x values of the fitted polynomial
        """
        self.log.info('starting the ``poly`` method')

        xvals = []

        # POLYNOMIALS SUMS
        for y in yarray:
            n_coeff = 0
            val = 0
            for i in range(0, self.deg + 1):
                val += coeff[n_coeff] * \
                    math.pow(y, i)
                n_coeff += 1
            xvals.append(val)

        self.log.info('completed the ``poly`` method')

        return xvals
