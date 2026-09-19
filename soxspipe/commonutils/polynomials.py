#!/usr/bin/env python
"""
*definition of polynomial functions needed throughout code*

Author
: David Young

Date Created
: September 10, 2020
"""

################# GLOBAL IMPORTS ####################

import os

os.environ["TERM"] = "vt100"


class chebyshev_order_wavelength_polynomials:
    """*the chebyshev polynomial fits for the single frames; to be iteratively fitted to minimise errors*

    **Key Arguments:**

    - ``log`` -- logger
    - ``orderDeg`` -- degree of the order polynomial components
    - ``wavelengthDeg`` -- degree of wavelength polynomial components
    - ``slitDeg`` -- degree of the slit polynomial components
    - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to
      regenerate. Default *False*
    - ``axis`` -- x, y or False. Default *False*.

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
    poly = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDeg, wavelengthDeg=wavelengthDeg, slitDeg=slitDeg).poly
    ```
    """

    def __init__(
        self, log, orderDeg, wavelengthDeg, slitDeg, exponentsIncluded=False, axis=False
    ):
        self.log = log
        self.orderDeg = orderDeg
        self.wavelengthDeg = wavelengthDeg
        self.slitDeg = slitDeg
        self.exponentsIncluded = exponentsIncluded

        if axis:
            self.axis = f"{axis.lower()}_"
        else:
            self.axis = ""

        return

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**

        - ``orderPixelTable`` -- a pandas dataframe containing wavelengths, orders and slit positions
        - ``coeff`` -- a list of the initial coefficients

        **Return:**

        - ``lhsVals`` -- the left-hand-side vals of the fitted polynomials
        """
        self.log.debug("starting the ``poly`` method")

        import numpy as np

        # UNPACK TUPLE INPUT
        orderDeg = self.orderDeg
        wavelengthDeg = self.wavelengthDeg
        slitDeg = self.slitDeg

        n = len(orderPixelTable.index)
        if n == 0:
            self.log.debug("completed the ``poly`` method")
            return np.zeros(0)

        # BUILD EACH VARIABLE'S POWER MATRIX ONCE (SHAPE (N, DEG+1)) INSTEAD
        # OF RECOMPUTING/RE-FETCHING POWERS PER POLYNOMIAL TERM
        if self.exponentsIncluded == False:
            orderVals = orderPixelTable["order"].to_numpy(dtype=float)
            wlVals = orderPixelTable["wavelength"].to_numpy(dtype=float)
            spVals = orderPixelTable["slit_position"].to_numpy(dtype=float)

            orderPow = np.power.outer(orderVals, np.arange(orderDeg + 1))
            wlPow = np.power.outer(wlVals, np.arange(wavelengthDeg + 1))
            spPow = np.power.outer(spVals, np.arange(slitDeg + 1))
        else:
            axis = self.axis
            orderPow = orderPixelTable[
                [f"order_pow_{axis}{i}" for i in range(orderDeg + 1)]
            ].to_numpy(dtype=float)
            wlPow = orderPixelTable[
                [f"wavelength_pow_{axis}{j}" for j in range(wavelengthDeg + 1)]
            ].to_numpy(dtype=float)
            spPow = orderPixelTable[
                [f"slit_position_pow_{axis}{k}" for k in range(slitDeg + 1)]
            ].to_numpy(dtype=float)

        # RESHAPE THE FLAT COEFF TUPLE INTO (I, J, K); ROW-MAJOR MATCHES THE
        # I-OUTER/J-MIDDLE/K-INNER ORDER THE COEFFS WERE ORIGINALLY WRITTEN IN.
        # NOTE: THIS RELIES ON ALL POLY DEGREES BEING SINGLE-DIGIT (TRUE FOR
        # ALL CURRENT SETTINGS FILES), A PRE-EXISTING ASSUMPTION, NOT
        # INTRODUCED HERE.
        coeffArr = np.asarray(coeff, dtype=float).reshape(
            orderDeg + 1, wavelengthDeg + 1, slitDeg + 1
        )

        # CONTRACT ONE AXIS AT A TIME SO WE NEVER MATERIALISE AN
        # (N, I, J, K) INTERMEDIATE ARRAY
        step1 = np.tensordot(orderPow, coeffArr, axes=([1], [0]))
        step2 = np.einsum("njk,nj->nk", step1, wlPow, optimize=True)
        lhsVals = np.einsum("nk,nk->n", step2, spPow, optimize=True)

        self.log.debug("completed the ``poly`` method")

        return lhsVals


class chebyshev_xy_polynomial:
    """*the chebyshev polynomial fits for the pinhole flat frame order tracing; to be iteratively fitted to minimise
    errors*

    **Key Arguments:**

    - ``log`` -- logger
    - ``y_deg`` -- y degree of the polynomial components
    - ``yCol`` -- name of the yCol. Default *False*
    - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to
      regenerate. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial
    poly = chebyshev_xy_polynomial(
            log=self.log, y_deg=y_deg).poly
    ```
    """

    def __init__(self, log, y_deg, yCol=False, exponentsIncluded=False):
        self.log = log
        self.y_deg = y_deg
        self.yCol = yCol
        self.exponentsIncluded = exponentsIncluded

        return

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**

        - ``orderPixelTable`` -- data frame with all pixel data arrays
        - ``coeff`` -- a list of the initial coefficients

        **Return:**

        - ``xvals`` -- the x values of the fitted polynomial
        """
        self.log.info("starting the ``poly`` method")

        import numpy as np
        import pandas as pd

        n_coeff = 0
        if not isinstance(orderPixelTable, pd.core.frame.DataFrame):
            yarray = np.array(orderPixelTable).astype("float")
            lhsVals = np.zeros(len(orderPixelTable))
        else:
            yarray = orderPixelTable[self.yCol].values.astype("float")
            lhsVals = np.zeros(len(orderPixelTable.index))

        if not self.exponentsIncluded:
            # POLYNOMIALS SUMS
            for i in range(0, self.y_deg + 1):
                lhsVals += coeff[n_coeff] * yarray**i
                n_coeff += 1
        else:
            for i in range(0, self.y_deg + 1):
                lhsVals += coeff[n_coeff] * orderPixelTable[f"y_pow_{i}"].values.astype(
                    "float"
                )
                n_coeff += 1

        self.log.info("completed the ``poly`` method")

        return lhsVals


class chebyshev_order_xy_polynomials:
    """*the chebyshev polynomial fits FIX ME*

    **Key Arguments:**

    - ``log`` -- logger
    - ``orderDeg`` -- degree of the order polynomial components
    - ``axisBDeg`` -- degree for polynomial to fit free axis-values
    - ``axisB`` -- the free axis related to `axisBDeg`. Default *'y'*. ['x'|'y']
    - ``axisBCol`` -- name of the free axis column (if needed). Default *False*
    - ``orderCol`` -- name of the order column (if needed). Default *False*
    - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to
      regenerate. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils.polynomials import chebyshev_order_xy_polynomials
    poly = chebyshev_order_xy_polynomials(
            log=self.log, orderDeg=orderDeg, axisBDeg=axisBDeg, axisB="y", axisBCol="y", orderCol="order").poly
    ```
    """

    def __init__(
        self,
        log,
        orderDeg,
        axisBDeg,
        axisB="y",
        axisBCol=False,
        orderCol=False,
        exponentsIncluded=False,
    ):
        self.log = log
        self.orderDeg = orderDeg
        self.axisBDeg = axisBDeg
        self.axisB = axisB
        self.axisBCol = axisBCol
        self.orderCol = orderCol
        self.exponentsIncluded = exponentsIncluded

        return

    def poly(self, orderPixelTable, *coeff):
        """the polynomial definition

        **Key Arguments:**

        - ``orderPixelTable`` -- a pandas dataframe containing x, y, order
        - ``coeff`` -- a list of the initial coefficients

        **Return:**

        - ``lhsVals`` -- the left-hand-side vals of the fitted polynomials
        """
        self.log.debug("starting the ``poly`` method")

        import numpy as np

        # UNPACK TUPLE INPUT
        orderDeg = int(self.orderDeg)
        axisBDeg = int(self.axisBDeg)
        axisB = self.axisB

        n_coeff = 0
        lhsVals = np.zeros(len(orderPixelTable.index))

        # FOR LOOPS ARE THE RIGHT TOOL TO PERFORM COMPUTATIONS OR RUN FUNCTIONS. LIST COMPREHENSION IS SLOW IN THESE
        # CASES

        if self.exponentsIncluded == False:
            orderVals = orderPixelTable[self.orderCol].values.astype("float")
            bVals = orderPixelTable[self.axisBCol].values.astype("float")

            for i in range(0, orderDeg + 1):
                for j in range(0, axisBDeg + 1):
                    lhsVals += coeff[n_coeff] * orderVals**i * bVals**j
                    n_coeff += 1
        else:
            for i in range(0, orderDeg + 1):
                for j in range(0, axisBDeg + 1):
                    lhsVals += (
                        float(coeff[n_coeff])
                        * orderPixelTable[f"order_pow_{i}"].values.astype("float")
                        * orderPixelTable[f"{axisB}_pow_{j}"].values.astype("float")
                    )
                    n_coeff += 1

        self.log.debug("completed the ``poly`` method")

        return lhsVals
