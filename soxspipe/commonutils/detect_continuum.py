#!/usr/bin/env python
# encoding: utf-8
"""
*find and fit the continuum trace across all echelle orders with low-order polynomials.*

Author
: David Young & Marco Landoni

Date Created
: September 10, 2020
"""
################# GLOBAL IMPORTS ####################


from soxspipe.commonutils.toolkit import read_spectral_format
from soxspipe.commonutils.toolkit import cut_image_slice
from soxspipe.commonutils.toolkit import get_calibration_lamp
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import chebyshev_xy_polynomial, chebyshev_order_xy_polynomials
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
import collections
from random import random
from os.path import expanduser
from fundamentals import tools
from builtins import object
import sys
import math
import os
from io import StringIO
import copy
from contextlib import suppress
from datetime import datetime
os.environ['TERM'] = 'vt100'


class _base_detect(object):

    def fit_order_polynomial(
            self,
            pixelList,
            order,
            axisBDeg,
            axisACol,
            axisBCol,
            exponentsIncluded=False):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**

        - ``pixelList`` -- data-frame group containing x,y pixel array
        - ``order`` -- the order to fit
        - ``axisBDeg`` -- degree for polynomial to fit
        - ``axisACol`` -- name of columns containing axis to be fitted
        - ``axisBCol`` -- name of columns containing free axis (values known)
        - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**

        - ``coeffs`` -- the coefficients of the polynomial fit
        - ``pixelList`` -- the pixel list but now with fits and residuals included
        """
        self.log.debug('starting the ``fit_order_polynomial`` method')

        import numpy as np
        from astropy.stats import sigma_clip
        from scipy.optimize import curve_fit

        arm = self.arm
        self.axisBDeg = axisBDeg

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log, axisBCol=axisBCol, axisBDeg=self.axisBDeg, exponentsIncluded=exponentsIncluded).poly

        clippingSigma = self.recipeSettings[
            "poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings[
            "poly-clipping-iteration-limit"]

        iteration = 0
        mask = (pixelList['order'] == order)
        pixelListFiltered = pixelList.loc[mask]

        coeff = np.ones((self.axisBDeg + 1))
        while clippedCount > 0 and iteration < clippingIterationLimit:
            pixelListFiltered = pixelList.loc[mask]

            startCount = len(pixelListFiltered.index)
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY
            # NOTE X AND Y COLUMN ARE CORRECLY IN xdata AND ydata - WANT TO
            # FIND X (UNKNOWN) WRT Y (KNOWNN)
            try:
                coeff, pcov_x = curve_fit(
                    poly, xdata=pixelListFiltered, ydata=pixelListFiltered[axisACol].values, p0=coeff, maxfev=30000)
            except TypeError as e:
                # REMOVE THIS ORDER FROM PIXEL LIST
                pixelList = pixelList.loc[~mask]
                coeff = None
                return coeff, pixelList
            except Exception as e:
                raise e

            res, res_mean, res_std, res_median, xfit, poly = self.calculate_residuals(
                orderPixelTable=pixelListFiltered,
                coeff=coeff,
                axisACol=axisACol,
                axisBCol=axisBCol)

            pixelList.loc[mask, "x_fit_res"] = res
            pixelList.loc[mask, "x_fit"] = xfit

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                res, sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
            pixelList.loc[mask, "mask"] = masked_residuals.mask

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = (pixelList["mask"] == True)
            pixelList = pixelList.loc[~removeMask]
            pixelListFiltered = pixelList.loc[mask]
            clippedCount = startCount - len(pixelListFiltered.index)

            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print(f'\t\tORDER {order:0.0f}: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial')

        self.log.debug('completed the ``fit_order_polynomial`` method')
        return coeff, pixelList

    def fit_global_polynomial(
            self,
            pixelList,
            axisACol="cont_x",
            axisBCol="cont_y",
            orderCol="order",
            exponentsIncluded=False,
            writeQCs=False):
        """*iteratively fit the global polynomial to the data, fitting axisA as a function of axisB, clipping residuals with each iteration*

        **Key Arguments:**

        - ``pixelList`` -- data-frame group containing x,y pixel array
        - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**

        - ``coeffs`` -- the coefficients of the polynomial fit
        - ``pixelList`` -- the pixel list but now with fits and residuals included
        - ``allClipped`` -- data that was sigma-clipped
        """
        self.log.debug('starting the ``fit_global_polynomial`` method')

        import numpy as np
        from astropy.stats import sigma_clip
        from scipy.optimize import curve_fit
        import pandas as pd

        arm = self.arm

        clippedCount = 1

        poly = chebyshev_order_xy_polynomials(log=self.log, axisBCol=axisBCol, orderCol=orderCol, orderDeg=self.orderDeg, axisBDeg=self.axisBDeg, axisB=self.axisB, exponentsIncluded=exponentsIncluded).poly

        clippingSigmaHigh = self.recipeSettings[
            "poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings[
            "poly-clipping-iteration-limit"]

        if axisACol == "stddev":
            clippingSigmaLow = clippingSigmaHigh
            clippingSigmaHigh = clippingSigmaLow * 2
        else:
            clippingSigmaLow = clippingSigmaHigh

        iteration = 0

        allClipped = []
        coeff = np.ones((self.axisBDeg + 1) * (self.orderDeg + 1))
        while clippedCount > 0 and iteration < clippingIterationLimit:
            startCount = len(pixelList.index)
            iteration += 1

            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY

            if len(pixelList.index) == 0:
                # REMOVE THIS ORDER FROM PIXEL LIST
                coeff = None
                return coeff, pixelList, pixelList

            if iteration < 3:
                coeff = np.ones((self.axisBDeg + 1) * (self.orderDeg + 1))

            try:
                coeff, pcov_x = curve_fit(
                    poly, xdata=pixelList, ydata=pixelList[axisACol].values, p0=coeff, maxfev=30000)
            except TypeError as e:
                # REMOVE THIS ORDER FROM PIXEL LIST
                coeff = None
                return coeff, pixelList, pixelList
            except Exception as e:
                raise e

            res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
                orderPixelTable=pixelList,
                coeff=coeff,
                orderCol=orderCol,
                axisACol=axisACol,
                axisBCol=axisBCol)

            pixelList[f"{axisACol}_fit_res"] = res
            pixelList[f"{axisACol}_fit"] = xfit

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                res, sigma_lower=clippingSigmaLow, sigma_upper=clippingSigmaHigh, maxiters=1, cenfunc='mean', stdfunc='std')
            pixelList["mask"] = masked_residuals.mask

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = (pixelList["mask"] == True)
            allClipped.append(pixelList.loc[removeMask])
            pixelList = pixelList.loc[~removeMask]
            clippedCount = startCount - len(pixelList.index)

            if iteration > 1:
                sys.stdout.flush()
                sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print(f'\t\tGLOBAL FIT: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial')

        if len(allClipped):
            allClipped = pd.concat(allClipped, ignore_index=True)

        res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
            orderPixelTable=pixelList,
            coeff=coeff,
            orderCol=orderCol,
            axisACol=axisACol,
            axisBCol=axisBCol,
            writeQCs=writeQCs)

        self.log.debug('completed the ``fit_global_polynomial`` method')
        return coeff, pixelList, allClipped

    def calculate_residuals(
            self,
            orderPixelTable,
            coeff,
            axisACol,
            axisBCol,
            orderCol=False,
            writeQCs=False):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**

        - ``orderPixelTable`` -- data-frame containing pixel list for given order
        - ``coeff`` -- the coefficients of the fitted polynomial
        - ``axisACol`` -- name of x-pixel column
        - ``axisBCol`` -- name of y-pixel column
        - ``orderCol`` -- name of the order column (global fits only)
        - ``writeQCs`` -- write the QCs to dataframe? Default *False*

        **Return:**

        - ``res`` -- x residuals
        - ``mean`` -- the mean of the residuals
        - ``std`` -- the stdev of the residuals
        - ``median`` -- the median of the residuals
        - ``xfit`` -- fitted x values
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        import numpy as np
        import pandas as pd

        arm = self.arm

        poly = chebyshev_order_xy_polynomials(
            log=self.log, axisBCol=axisBCol, orderCol=orderCol, orderDeg=self.orderDeg, axisBDeg=self.axisBDeg).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        thisFit = poly(
            orderPixelTable, *coeff)
        res = thisFit - orderPixelTable[axisACol].values

        # GET UNIQUE VALUES IN COLUMN
        uniqueorders = len(orderPixelTable['order'].unique())

        # CALCULATE COMBINED RESIDUALS AND STATS
        res_mean = np.ma.mean(res)
        res_std = np.ma.std(res)
        res_median = np.ma.median(res)

        if writeQCs:
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

            tag = "continuum"
            if "order-centre" in self.recipeName.lower():
                tag = "order centre"

            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": f"{self.axisA.upper()}RESMIN",
                "qc_value": f"{res.min():0.2f}",
                "qc_comment": f"[px] Minimum residual in {tag} fit along {self.axisA}-axis",
                "qc_unit": "px",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": f"{self.axisA.upper()}RESMAX",
                "qc_value": f"{res.max():0.2f}",
                "qc_comment": f"[px] Maximum residual in {tag} fit along {self.axisA}-axis",
                "qc_unit": "px",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": f"{self.axisA.upper()}RESRMS",
                "qc_value": f"{res_std:0.2f}",
                "qc_comment": f"[px] Std-dev of residual {tag} fit along {self.axisA}-axis",
                "qc_unit": "px",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

            if "order" in self.recipeName.lower():
                c = f"Number of order centre traces found"
            else:
                c = f"Number of orders containing an object trace"

            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "NORDERS",
                "qc_value": uniqueorders,
                "qc_comment": c,
                "qc_unit": None,
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``calculate_residuals`` method')
        return res, res_mean, res_std, res_median, thisFit

    def write_order_table_to_file(
            self,
            frame,
            orderPolyTable,
            orderMetaTable):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**

        - ``frame`` -- the calibration frame used to generate order location data
        - ``orderPolyTable`` -- data-frames containing centre location coefficients (and possibly also order edge coeffs)
        - ``orderMetaTable`` -- extra order meta data to be added in an extra FITS extension

        **Return:**

        - ``order_table_path`` -- path to the order table file
        """
        from astropy.table import Table
        from astropy.io import fits
        self.log.debug('starting the ``write_order_table_to_file`` method')

        arm = self.arm
        kw = self.kw

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        if (False and (self.binx > 1 or self.biny > 1)) or (isinstance(self.products, bool) and self.products == False):
            outDir = self.settings["workspace-root-dir"] + "/tmp"
        else:
            outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
            outDir = outDir.replace("//", "/")
        # Recursively create missing directories
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        if not self.sofName:
            filename = filenamer(
                log=self.log,
                frame=frame,
                settings=self.settings
            )

        else:
            filename = self.sofName + ".fits"
        filename = filename.replace("MFLAT", "FLAT")

        # if self.inst.upper() == "SOXS":
        #     filename = filename.replace("_DLAMP", "")
        #     filename = filename.replace("_QLAMP", "")

        if "mflat" in self.recipeName.lower():
            filename = filename.upper().split("FLAT")[0] + "ORDER_LOCATIONS.fits"
        elif "stare" in self.recipeName.lower():
            filename = filename.upper().split(".FITS")[0] + "_OBJECT_TRACE.fits"
        elif "nod" in self.recipeName.lower():
            # sequence = "A" if int(frame.header['HIERARCH ESO SEQ CUMOFF Y'] > 0) else "B"
            filename = filename.upper().split(".FITS")[0] + "_OBJECT_TRACE" + self.noddingSequence + ".fits"

        if self.lampTag and self.inst.upper() != "SOXS":
            filename = filename.replace(".fits", f"{self.lampTag}.fits")

        order_table_path = f"{outDir}/{filename}"

        header = copy.deepcopy(frame.header)
        with suppress(KeyError):
            header.pop(kw("DPR_TECH"))
        with suppress(KeyError):
            header.pop(kw("DPR_CATG"))
        with suppress(KeyError):
            header.pop(kw("DPR_TYPE"))
        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        header["HIERARCH " + kw("PRO_TECH")] = "ECHELLE,SLIT"

        orderPolyTable = Table.from_pandas(orderPolyTable)
        BinTableHDU = fits.table_to_hdu(orderPolyTable)
        orderMetaTable = Table.from_pandas(orderMetaTable)
        BinTableHDU2 = fits.table_to_hdu(orderMetaTable)

        header[kw("SEQ_ARM")] = arm
        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        if "stare" not in self.recipeName.lower() and "nod" not in self.recipeName.lower() and "offset" not in self.recipeName.lower():
            header["HIERARCH " + kw("PRO_CATG")] = f"ORDER_TAB_{arm}".upper()
        else:
            header["HIERARCH " + kw("PRO_CATG")] = f"OBJECT_TAB_{arm}".upper()
        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU, BinTableHDU2])
        hduList.writeto(order_table_path, checksum=True, overwrite=True)

        self.log.debug('completed the ``write_order_table_to_file`` method')
        return order_table_path


class detect_continuum(_base_detect):
    """
    *find and fit the continuum trace across all echelle orders with low-order polynomials.*

    **Key Arguments:**

    - ``log`` -- logger
    - ``traceFrame`` -- calibrated frame containing a source trace (CCDObject)
    - ``dispersion_map`` -- path to dispersion map file containing polynomial fits of the dispersion solution for the frame
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``recipeName`` -- the recipe name as given in the settings dictionary
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products
    - ``sofName`` ---- name of the originating SOF file
    - ``binx`` -- binning in x-axis
    - ``biny`` -- binning in y-axis
    - ``lampTag`` -- add this tag to the end of the product filename. Default *False*
    - ``locationSetIndex`` -- the index of the AB cycle locations (nodding mode only). Default *False*
    - ``orderPixelTable`` -- this is used for tuning the pipeline.  Default *False*

    **Usage:**

    To use the ``detect_continuum`` object, use the following:

    ```python
    from soxspipe.commonutils import detect_continuum
    detector = detect_continuum(
        log=log,
        traceFrame=traceFrame,
        dispersion_map=dispersion_map,
        settings=settings,
        recipeName="soxs-order-centre"
    )
    order_table_path = detector.get()
    ```
    """

    def __init__(
            self,
            log,
            traceFrame,
            dispersion_map,
            settings=False,
            recipeSettings=False,
            recipeName=False,
            qcTable=False,
            productsTable=False,
            sofName=False,
            binx=1,
            biny=1,
            lampTag=False,
            locationSetIndex=False,
            orderPixelTable=False
    ):
        self.log = log
        log.debug("instantiating a new 'detect_continuum' object")

        import copy

        self.settings = settings
        try:
            self.noddingSequence = "_A" if int(traceFrame.header['HIERARCH ESO SEQ CUMOFF Y'] > 0) else "_B"
            if locationSetIndex:
                self.noddingSequence += str(locationSetIndex)
        except:
            self.noddingSequence = ""

        self.recipeName = recipeName
        self.traceFrame = traceFrame
        self.dispersion_map = dispersion_map
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.binx = binx
        self.biny = biny
        self.recipeSettings = copy.deepcopy(recipeSettings["detect-continuum"])
        self.lampTag = lampTag
        self.orderPixelTable = orderPixelTable

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.arm = traceFrame.header[self.kw("SEQ_ARM")]
        self.dateObs = traceFrame.header[self.kw("DATE_OBS")]
        self.inst = traceFrame.header[self.kw("INSTRUME")]
        self.exptime = traceFrame.header[self.kw("EXPTIME")]

        # if self.exptime < 59 and recipeName != "soxs-stare":
        #     raise Exception("too low")

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.axisBDeg = self.recipeSettings["disp-axis-deg"]
        self.orderDeg = self.recipeSettings["order-deg"]

        self.lamp = get_calibration_lamp(log=log, frame=traceFrame, kw=self.kw)

        home = expanduser("~")
        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.qcDir = self.qcDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.qcDir):
            os.makedirs(self.qcDir)

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
            self.coeff_dict = {"degorder_cent": self.orderDeg,
                               "degy_cent": self.axisBDeg}
        else:
            self.axisA = "y"
            self.axisB = "x"
            self.coeff_dict = {"degorder_cent": self.orderDeg,
                               "degx_cent": self.axisBDeg}

        return None

    def get(self):
        """
        *return the order centre table filepath*

        **Return:**

        - ``order_table_path`` -- file path to the order centre table giving polynomial coeffs to each order fit
        """
        self.log.debug('starting the ``get`` method')

        import numpy as np
        import pandas as pd
        from datetime import datetime

        arm = self.arm
        coeff_dict = self.coeff_dict

        if isinstance(self.orderPixelTable, bool):
            orderPixelTable = self.sample_trace()
        else:
            orderPixelTable = self.orderPixelTable

        self.log.print("\n\t## FINDING GLOBAL POLYNOMIAL SOLUTION FOR CONTINUUM TRACES\n")

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = orderPixelTable['order'].unique()

        orderLocations = {}
        orderPixelTable[f'cont_{self.axisA}_fit'] = np.nan
        orderPixelTable[f'cont_{self.axisA}_fit_res'] = np.nan

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        fitFound = False
        tryCount = 0
        backupOrderPixelTable = orderPixelTable.copy()
        while not fitFound and tryCount < 5:
            # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
            for i in range(0, self.axisBDeg + 1):
                orderPixelTable[f"{self.axisB}_pow_{i}"] = orderPixelTable[f"cont_{self.axisB}"].pow(i)
            for i in range(0, self.orderDeg + 1):
                orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)
            try:

                coeff, orderPixelTable, clippedDataCentre = self.fit_global_polynomial(
                    pixelList=orderPixelTable,
                    axisACol=f"cont_{self.axisA}",
                    axisBCol=f"cont_{self.axisB}",
                    exponentsIncluded=True,
                    writeQCs=True
                )
                mean_res = np.mean(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))

                if "order" in self.recipeName.lower() and mean_res > 1:
                    # BAD FIT ... FORCE A FAIL
                    orderPixelTable = backupOrderPixelTable
                    raise AttributeError("Failed to continuum trace")

                n_coeff = 0
                for i in range(0, self.orderDeg + 1):
                    for j in range(0, self.axisBDeg + 1):
                        coeff_dict[f'cent_{i}{j}'] = coeff[n_coeff]
                        n_coeff += 1

                # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
                coeff, orderPixelTable, clippedData = self.fit_global_polynomial(
                    pixelList=orderPixelTable,
                    axisACol="stddev",
                    axisBCol=f"cont_{self.axisB}",
                    exponentsIncluded=True
                )

                if len(clippedData) and len(clippedData.index):
                    clippedData = pd.concat([clippedDataCentre, clippedData], ignore_index=True)
                else:
                    clippedData = clippedDataCentre
                fitFound = True

            except Exception as e:
                degList = [self.axisBDeg, self.orderDeg]
                degList[degList.index(max(degList))] -= 1
                self.axisBDeg, self.orderDeg = degList
                coeff_dict["degorder_cent"] = self.orderDeg
                coeff_dict[f"deg{self.axisB}_cent"] = self.axisBDeg
                self.log.print(f"{self.axisB} and Order fitting orders reduced to {self.axisBDeg}, {self.orderDeg} to try and successfully fit the continuum.")
                self.recipeSettings["disp-axis-deg"] = self.axisBDeg
                self.recipeSettings["order-deg"] = self.orderDeg
                tryCount += 1
                if tryCount == 5:
                    self.log.print(f"Could not converge on a good fit to the continuum. Please check the quality of your data or adjust your fitting parameters.")
                    raise e
                if self.settings["tune-pipeline"]:
                    raise e

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        try:
            CSAMP = len(clippedData.index)
        except:
            CSAMP = 0
        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "CSAMP",
            "qc_value": CSAMP,
            "qc_comment": "Number of continuum sample clipped during solution fitting",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        # orderLocations[o] = coeff
        coeff_dict["degorder_std"] = self.orderDeg
        coeff_dict[f"deg{self.axisB}_std"] = self.axisBDeg
        n_coeff = 0
        for i in range(0, self.orderDeg + 1):
            for j in range(0, self.axisBDeg + 1):
                coeff_dict[f'std_{i}{j}'] = coeff[n_coeff]
                n_coeff += 1
        coeffColumns = coeff_dict.keys()

        orderPolyTable = pd.DataFrame([coeff_dict])

        # HERE IS THE LINE LIST IF NEEDED FOR QC
        orderPixelTable.drop(columns=['mask'], inplace=True)

        plotPath, orderMetaTable = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable,
            clippedData=clippedData
        )

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        basename = os.path.basename(plotPath)

        if "order" in self.recipeName.lower():
            label = "ORDER_CENTRES_RES"
            product_desc = f"Residuals of the order centre polynomial fit"
        else:
            label = "OBJECT_TRACE_RES"
            product_desc = f"Residuals of the object trace polynomial fit"

        if not isinstance(self.products, bool):
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": label + self.noddingSequence,
                "file_name": basename,
                "file_type": "PDF",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": product_desc,
                "file_path": plotPath,
                "label": "QC"
            }).to_frame().T], ignore_index=True)
            # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
            order_table_path = self.write_order_table_to_file(
                frame=self.traceFrame, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)
        else:
            order_table_path = self.write_order_table_to_file(
                frame=self.traceFrame, orderPolyTable=orderPolyTable, orderMetaTable=orderMetaTable)

        # mean_res = np.mean(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))
        # std_res = np.std(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))

        self.log.debug('completed the ``get`` method')
        return order_table_path, self.qc, self.products, orderPolyTable, orderPixelTable, orderMetaTable

    def create_pixel_arrays(
            self):
        """*create a pixel array for the approximate centre of each order*

        **Return:**

        - ``orderPixelTable`` -- a data-frame containing lines and associated pixel locations
        """
        self.log.debug('starting the ``create_pixel_arrays`` method')

        import numpy as np
        import pandas as pd

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax, amins, amaxs = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm, dispersionMap=self.dispersion_map)

        # READ ORDER SAMPLING RESOLUTION FROM SETTINGS
        sampleCount = self.recipeSettings["order-sample-count"]

        orderPixelRanges = []
        for o, amin, amax in zip(orderNums, amins, amaxs):
            arange = amax - amin
            orderPixelRanges.append(arange)

        smallestRange = min(orderPixelRanges)
        samplePixelSep = int(smallestRange / sampleCount)

        # CREATE THE WAVELENGTH/ORDER ARRAYS TO BE CONVERTED TO PIXELS
        myDict = {
            "order": np.asarray([]),
            "wavelength": np.asarray([]),
            "slit_position": np.asarray([])
        }

        for o, wmin, wmax, pixelRange in zip(orderNums, waveLengthMin, waveLengthMax, orderPixelRanges):

            orderSampleCount = int(pixelRange / samplePixelSep)
            wrange = wmax - wmin

            wlArray = np.arange(
                wmin, wmax, (wmax - wmin) / orderSampleCount)
            myDict["wavelength"] = np.append(myDict["wavelength"], wlArray)
            myDict["order"] = np.append(
                myDict["order"], np.ones(len(wlArray)) * o)
            myDict["slit_position"] = np.append(
                myDict["slit_position"], np.zeros(len(wlArray)))

        orderPixelTable = pd.DataFrame(myDict)
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersion_map,
            orderPixelTable=orderPixelTable
        )

        self.log.debug('completed the ``create_pixel_arrays`` method')
        return orderPixelTable

    def fit_1d_gaussian_to_slice(
            self,
            pixelPostion):
        """*cut a slice from the pinhole flat along the cross-dispersion direction centred on pixel position, fit 1D gaussian and return the peak pixel position*

        **Key Arguments:**

        - ``pixelPostion`` -- the x,y pixel coordinate from orderPixelTable data-frame (series)

        **Return:**

        - ``pixelPostion`` -- now including gaussian fit peak xy position
        """
        self.log.debug('starting the ``fit_1d_gaussian_to_slice`` method')

        import numpy as np
        from astropy.stats import mad_std
        from astropy.modeling import models, fitting
        from scipy.signal import find_peaks

        # CLIP OUT A SLICE TO INSPECT CENTRED AT POSITION
        halfSlice = self.sliceLength / 2

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            sliceAxis = "x"
            sliceAntiAxis = "y"
        else:
            sliceAxis = "y"
            sliceAntiAxis = "x"

        slice, slice_length_offset, slice_width_centre = cut_image_slice(log=self.log, frame=self.traceFrame,
                                                                         width=self.sliceWidth, length=self.sliceLength, x=pixelPostion["fit_x"], y=pixelPostion["fit_y"], sliceAxis=sliceAxis, median=True, plot=False)

        if slice is None:
            pixelPostion[f"cont_{self.axisA}"] = np.nan
            pixelPostion[f"cont_{self.axisB}"] = np.nan
            return pixelPostion

        # CHECK THE SLICE POINTS IF NEEDED
        if 1 == 0:
            import matplotlib.pyplot as plt
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            plt.show()

        origSlice = slice.copy()

        slice[slice < -500] = np.nan

        # EVALUATING THE MEAN AND STD-DEV FOR PEAK FINDING - REMOVES SLICE
        # CONTAINING JUST NOISE
        try:
            median_r = np.nanmedian(slice)
            std_r = mad_std(slice, ignore_nan=True)
        except:
            median_r = None

        # print(median_r)
        # print()

        if not median_r:
            pixelPostion[f"cont_{self.axisA}"] = np.nan
            pixelPostion[f"cont_{self.axisB}"] = np.nan
            return pixelPostion

        peaks, _ = find_peaks(slice, height=median_r +
                              self.peakSigmaLimit * std_r, width=1)

        # CHECK PEAK HAS BEEN FOUND
        if peaks is None or len(peaks) <= 0:
            # CHECK THE SLICE POINTS IF NEEDED
            if False:
                import matplotlib.pyplot as plt
                x = np.arange(0, len(origSlice))
                plt.figure(figsize=(8, 5))
                plt.plot(x, slice, 'ko')
                plt.xlabel('Position')
                plt.ylabel('Flux')
                plt.show()
                print(median_r, std_r)
                import matplotlib.pyplot as plt
                x = np.arange(0, len(slice))
                plt.figure(figsize=(8, 5))
                plt.plot(x, slice, 'ko')
                plt.xlabel('Position')
                plt.ylabel('Flux')
                plt.show()
            pixelPostion[f"cont_{self.axisA}"] = np.nan
            pixelPostion[f"cont_{self.axisB}"] = np.nan
            return pixelPostion

        if len(peaks) > 1:
            closest = peaks[0]
            smallest_diff = abs(self.sliceLength / 2. - closest)  # Initialize the smallest difference
            for num in peaks:
                current_diff = abs(self.sliceLength / 2. - num)  # Calculate the difference
                if current_diff < smallest_diff:  # Check if current is closer
                    smallest_diff = current_diff
                    closest = num  # Update closest integer
            peaks = [closest]

        # FIT THE DATA USING A 1D GAUSSIAN - USING astropy.modeling
        # CENTRE THE GAUSSIAN ON THE PEAK
        g_init = models.Gaussian1D(
            amplitude=1000., mean=peaks[0], stddev=1.)
        # self.log.print(f"g_init: {g_init}")
        fit_g = fitting.LevMarLSQFitter()

        # NOW FIT
        try:
            # MASK OUT NAN VALUES
            mask = np.isfinite(slice)
            g = fit_g(g_init, np.arange(0, len(slice))[mask], slice[mask])
        except:
            pixelPostion[f"cont_{self.axisA}"] = np.nan
            pixelPostion[f"cont_{self.axisB}"] = np.nan
            return pixelPostion
        pixelPostion[f"cont_{sliceAxis}"] = g.mean.value + \
            max(0, slice_length_offset)
        pixelPostion[f"cont_{sliceAntiAxis}"] = slice_width_centre
        pixelPostion["amplitude"] = g.amplitude.value
        pixelPostion["stddev"] = g.stddev.value

        # PRINT A FEW PLOTS IF NEEDED - GAUSSIAN FIT OVERLAID
        if False and random() < 0.02:
            import matplotlib.pyplot as plt
            x = np.arange(0, len(slice))
            plt.figure(figsize=(8, 5))
            plt.plot(x, slice, 'ko')
            plt.xlabel('Position')
            plt.ylabel('Flux')
            gaussx = np.arange(0, max(x), 0.05)
            plt.plot(gaussx, g(gaussx), label=f'Mean = {g.mean.value:0.2f}, std = {g.stddev.value:0.2f}, cont_{self.axisA} = {pixelPostion[f"cont_{self.axisA}"]:0.2f}, fit_{self.axisA} = {pixelPostion[f"fit_{self.axisA}"]:0.2f},cont_{self.axisB} = {pixelPostion[f"cont_{self.axisB}"]:0.2f},fit_y = {pixelPostion[f"fit_{self.axisB}"]:0.2f}')
            plt.legend()
            plt.show()

        self.log.debug('completed the ``fit_1d_gaussian_to_slice`` method')

        return pixelPostion

    def plot_results(
            self,
            orderPixelTable,
            orderPolyTable,
            clippedData):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**

        - ``orderPixelTable`` -- the pixel table with residuals of fits
        - ``orderPolyTable`` -- data-frame of order-location polynomial coeff
        - ``clippedData`` -- the sigma-clipped data

        **Return:**

        - ``filePath`` -- path to the plot pdf
        - ``orderMetaTable`` -- dataframe of useful order fit metadata
        """
        self.log.debug('starting the ``plot_results`` method')

        import numpy as np
        import pandas as pd
        import matplotlib.pyplot as plt
        from soxspipe.commonutils.toolkit import qc_settings_plot_tables

        arm = self.arm

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotateImage = self.detectorParams["rotate-qc-plot"]
        flipImage = self.detectorParams["flip-qc-plot"]

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = self.traceFrame.data
        if rotateImage:
            rotatedImg = np.rot90(rotatedImg, rotateImage / 90)
        if flipImage:
            rotatedImg = np.flipud(rotatedImg)
            if not rotateImage:
                aLen = rotatedImg.shape[0]
                orderPixelTable[f"cont_{self.axisA}"] = aLen - orderPixelTable[f"cont_{self.axisA}"]
                clippedData[f"cont_{self.axisA}"] = aLen - clippedData[f"cont_{self.axisA}"]

        if rotatedImg.shape[0] / rotatedImg.shape[1] > 0.8:
            fig = plt.figure(figsize=(6, 13.5), constrained_layout=True)
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(7, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:5, 0:2])
            bottomright = fig.add_subplot(gs[4:5, 2:])
            fwhmaxis = fig.add_subplot(gs[5:6, :])
            settingsAx = fig.add_subplot(gs[6:, 2:])
            qcAx = fig.add_subplot(gs[6:, 0:2])
        else:
            fig = plt.figure(figsize=(6, 19), constrained_layout=True)
            # CREATE THE GID OF AXES
            gs = fig.add_gridspec(10, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:6, 0:2])
            bottomright = fig.add_subplot(gs[4:6, 2:])
            fwhmaxis = fig.add_subplot(gs[6:7, :])
            settingsAx = fig.add_subplot(gs[7:, 2:])
            qcAx = fig.add_subplot(gs[7:, 0:2])

        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "1D guassian peak positions (post-clipping)", fontsize=10)
        toprow.scatter(orderPixelTable[f"cont_{self.axisB}"], orderPixelTable[f"cont_{self.axisA}"], marker='o', c='green', s=0.3, alpha=0.6, label="cross-dispersion 1D gaussian peak-position")
        toprow.scatter(clippedData[f"cont_{self.axisB}"], clippedData[f"cont_{self.axisA}"], marker='x', c='red', s=5, alpha=0.6, linewidths=0.5, label="peaks clipped during continuum fitting")
        # Put a legend below current axis

        toprow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.1),
                      fontsize=4)

        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])

        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        toprow.tick_params(axis='both', which='major', labelsize=9)
        toprow.set_xlim([0, rotatedImg.shape[1]])

        if self.detectorParams["dispersion-axis"] == "x":
            toprow.invert_yaxis()
            midrow.invert_yaxis()
        toprow.set_ylim([0, rotatedImg.shape[0]])
        midrow.set_ylim([0, rotatedImg.shape[0]])

        if "order" in self.recipeName.lower():
            midrow.set_title(
                "order-location fit solutions", fontsize=10)
        else:
            midrow.set_title(
                "global polynomal fit of order-centres", fontsize=10)
        if self.axisB == "y":
            axisALength = self.traceFrame.data.shape[1]
            axisBLength = self.traceFrame.data.shape[0]
        elif self.axisB == "x":
            axisALength = self.traceFrame.data.shape[0]
            axisBLength = self.traceFrame.data.shape[1]

        axisBlinelist = np.arange(0, axisBLength, 3)

        poly = chebyshev_order_xy_polynomials(
            log=self.log, axisBCol=self.axisB, orderCol="order", orderDeg=self.orderDeg, axisBDeg=self.axisBDeg).poly
        for index, row in orderPolyTable.iterrows():
            cent_coeff = [float(v) for k, v in row.items() if "cent_" in k]
            std_coeff = [float(v) for k, v in row.items() if "std_" in k]

        uniqueOrders = orderPixelTable['order'].unique()
        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {f"{self.axisB}": axisBlinelist}
        df = pd.DataFrame(myDict)
        ymin = []
        ymax = []
        xmin = []
        xmax = []
        foundOrders = []
        colors = []
        labelAdded = None

        for o in uniqueOrders:
            df["order"] = o
            xfit = poly(df, *cent_coeff)
            stdfit = poly(df, *std_coeff)

            try:
                xfit, yfit, stdfit, lower, upper = zip(
                    *[(x, y, std, x - 3 * std, x + 3 * std) for x, y, std in zip(xfit, axisBlinelist, stdfit) if x > 0 and x < (axisALength) - 10])
            except:
                continue
            if flipImage and not rotateImage:
                xfit = aLen - np.array(xfit)
                lower = aLen - np.array(lower)
                upper = aLen - np.array(upper)
            foundOrders.append(o)
            # lower = xfit - 3 * stdfit
            # upper = xfit + 3 * stdfit
            if labelAdded == None:
                label1 = "$3\sigma$ deviation"
                label2 = "polynomial fit"
                labelAdded = True
            else:
                label1 = None
                label2 = None
            c = midrow.plot(yfit, xfit, linewidth=0.7, label=label2)

            midrow.fill_between(yfit, lower, upper, color=c[0].get_color(), alpha=0.3, label=label1)
            colors.append(c[0].get_color())
            ymin.append(min(yfit))
            ymax.append(max(yfit))
            xmin.append(axisALength - max(xfit))
            xmax.append(axisALength - min(xfit))
            try:
                midrow.text(yfit[10], xfit[10] + 10, int(o), fontsize=6, c=c[0].get_color(), verticalalignment='bottom')
            except:
                pass

        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        orderMetaTable = {
            "order": foundOrders,
            f"{self.axisB}min": ymin,
            f"{self.axisB}max": ymax,
            f"{self.axisA}min": xmin,
            f"{self.axisA}max": xmax,
        }
        orderMetaTable = pd.DataFrame(orderMetaTable)

        # xfit = np.ones(len(xfit)) * \
        #     self.pinholeFrame.data.shape[1] - xfit
        # midrow.scatter(yfit, xfit, marker='x', c='blue', s=4)
        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis='both', which='major', labelsize=9)

        midrow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.1),
                      fontsize=4)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        for o, c in zip(uniqueOrders, colors):
            mask = (orderPixelTable['order'] == o)
            bottomleft.scatter(orderPixelTable.loc[mask][f'cont_{self.axisA}'].values, orderPixelTable.loc[mask][
                f'cont_{self.axisA}_fit_res'].values, alpha=0.2, s=1, c=c)
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                bottomleft.text(orderPixelTable.loc[mask][f'cont_{self.axisA}'].values[labelIndex], orderPixelTable.loc[mask][f'cont_{self.axisA}_fit_res'].values[labelIndex], int(o), fontsize=8, c=c, verticalalignment='bottom')
            except:
                pass

        bottomleft.set_xlabel(f'{self.axisA} pixel position', fontsize=10)
        bottomleft.set_ylabel(f'{self.axisA} residual', fontsize=10)
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        for o, c in zip(uniqueOrders, colors):
            mask = (orderPixelTable['order'] == o)
            bottomright.scatter(orderPixelTable.loc[mask][f'cont_{self.axisB}'].values, orderPixelTable.loc[mask][
                f'cont_{self.axisA}_fit_res'].values, alpha=0.2, s=1, c=c)
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                bottomright.text(orderPixelTable.loc[mask][f'cont_{self.axisB}'].values[labelIndex], orderPixelTable.loc[mask][
                    f'cont_{self.axisA}_fit_res'].values[labelIndex], int(o), fontsize=8, c=c, verticalalignment='bottom')
            except:
                pass

        bottomright.set_xlabel(f'{self.axisB} pixel position', fontsize=10)
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        # bottomright.set_ylabel('x residual')
        bottomright.set_yticklabels([])

        stdToFwhm = 2 * (2 * math.log(2))**0.5
        for o, c in zip(uniqueOrders, colors):
            mask = (orderPixelTable['order'] == o)
            fwhmaxis.scatter(orderPixelTable.loc[mask]['wavelength'].values, orderPixelTable.loc[mask]['stddev_fit'].values * stdToFwhm, alpha=0.2, s=1, c=c)
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                fwhmaxis.text(orderPixelTable.loc[mask]['wavelength'].values[labelIndex], orderPixelTable.loc[mask]['stddev_fit'].values[labelIndex] * stdToFwhm, int(o), fontsize=8, c=c, verticalalignment='bottom')
            except:
                pass
        fwhmaxis.set_xlabel('wavelength (nm)', fontsize=10)
        fwhmaxis.set_ylabel('Cross-dispersion\nFWHM (pixels)', fontsize=10)

        fwhmaxis.set_ylim(orderPixelTable['stddev_fit'].min() * stdToFwhm * 0.5, orderPixelTable['stddev_fit'].max() * stdToFwhm * 1.2)

        # REMOVE DUPLICATE ENTRIES IN COLUMN 'qc_name' AND KEEP THE LAST ENTRY
        self.qc = self.qc.drop_duplicates(subset=['qc_name'], keep='last')
        qc_settings_plot_tables(log=self.log, qc=self.qc, qcAx=qcAx, settings={**self.recipeSettings, **{"exptime": self.exptime}}, settingsAx=settingsAx)

        mean_res = np.mean(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))
        std_res = np.std(orderPixelTable[f'cont_{self.axisA}_fit_res'].values)
        res_min = np.min(orderPixelTable[f'cont_{self.axisA}_fit_res'].values)
        res_max = np.max(orderPixelTable[f'cont_{self.axisA}_fit_res'].values)
        res_range = res_max - res_min

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        lamp = ""
        if self.lamp:
            lamp = f" {self.lamp} lamp"
        if "order" in self.recipeName.lower():
            fig.suptitle(f"traces of order-centre locations - {arm}{lamp} pinhole flat-frame\n{subtitle}", fontsize=12, y=0.99)
        else:
            fig.suptitle(f"{arm} object trace locations\n{subtitle}", fontsize=12)

        polyOrders = [self.orderDeg, self.axisBDeg]
        polyOrders[:] = [str(l) for l in polyOrders]
        polyOrders = "".join(polyOrders)

        # plt.show()
        if not self.sofName:
            filename = filenamer(
                log=self.log,
                frame=self.traceFrame,
                settings=self.settings
            )
            filename = filename.split("FLAT")[0] + "ORDER_CENTRES_residuals.pdf"
        elif "order" in self.recipeName.lower():
            filename = self.sofName + f"_residuals_{polyOrders}.pdf"
        elif "nod" in self.recipeName.lower():
            filename = self.sofName + "_OBJECT_TRACE_residuals" + self.noddingSequence + f"_{polyOrders}.pdf"
        else:
            filename = self.sofName + f"_OBJECT_TRACE_residuals_{polyOrders}.pdf"

        filePath = f"{self.qcDir}/{filename}"
        plt.tight_layout()
        if not self.settings["tune-pipeline"]:
            plt.savefig(filePath, dpi=720)
        plt.close()

        if self.settings["tune-pipeline"]:
            import codecs
            filePath = f"residuals.txt"
            exists = os.path.exists(filePath)
            if not exists:
                with codecs.open(filePath, encoding='utf-8', mode='w') as writeFile:
                    writeFile.write(f"polyOrders,mean_res,std_res,res_min,res_max,res_range \n")
            with codecs.open(filePath, encoding='utf-8', mode='a') as writeFile:
                writeFile.write(f"{polyOrders},{mean_res:2.4f},{std_res:2.4f},{res_min:2.4f},{res_max:2.4f},{res_range:2.4f}\n")

        self.log.debug('completed the ``plot_results`` method')
        return filePath, orderMetaTable

    def sample_trace(
            self):
        """*take many cross-dispersion samples across each order to try and find an object trace*

        **Return:**

        - ``orderPixelTable`` -- the detector locations at which a trace was found
        """
        self.log.debug('starting the ``sample_trace`` method')

        import pandas as pd
        import numpy as np

        # CONVERT WAVELENGTH TO PIXEL POSITIONS AND RETURN ARRAY OF POSITIONS TO
        # SAMPLE THE TRACES
        orderPixelTable = self.create_pixel_arrays()

        binx = 1
        biny = 1
        try:
            binx = self.traceFrame.header[self.kw("WIN_BINX")]
            biny = self.traceFrame.header[self.kw("WIN_BINY")]
        except:
            pass

        # FIT_X AND FIT_Y FROM DISP-SOLUTION
        orderPixelTable['fit_x'] = orderPixelTable['fit_x'] / binx
        orderPixelTable['fit_y'] = orderPixelTable['fit_y'] / biny

        # SLICE LENGTH TO SAMPLE TRACES IN THE CROSS-DISPERSION DIRECTION
        self.sliceLength = self.recipeSettings["slice-length"]
        self.peakSigmaLimit = self.recipeSettings["peak-sigma-limit"]
        self.sliceWidth = self.recipeSettings["slice-width"]

        if self.kw("DPR_TYPE").upper() in self.traceFrame.header and "STD" in self.traceFrame.header[self.kw("DPR_TYPE")].upper() and self.sliceWidth > 3:
            self.sliceWidth = 3

        # PREP LISTS WITH NAN VALUE IN CONT_X AND CONT_Y BEFORE FITTING
        orderPixelTable[f'cont_{self.axisA}'] = np.nan
        orderPixelTable[f'cont_{self.axisB}'] = np.nan

        # FOR EACH ORDER, FOR EACH PIXEL POSITION SAMPLE, FIT A 1D GAUSSIAN IN
        # CROSS-DISPERSION DIRECTION. RETURN PEAK POSITIONS
        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=self.traceFrame, show=False, ext='data', stdWindow=3, title=False, surfacePlot=True)

        if "order" in self.recipeName.lower():
            self.log.print("\n# FINDING & FITTING ORDER-CENTRE CONTINUUM TRACES\n")
        else:
            self.log.print("\n# FINDING & FITTING OBJECT CONTINUUM TRACES\n")

        orderPixelTable = orderPixelTable.apply(
            self.fit_1d_gaussian_to_slice, axis=1)
        allLines = len(orderPixelTable.index)
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (orderPixelTable['cont_x'] < 0)
        orderPixelTable.loc[mask, 'cont_x'] = np.nan
        mask = (orderPixelTable['cont_y'] < 0)
        orderPixelTable.loc[mask, 'cont_y'] = np.nan

        # xpd-update-filter-dataframe-column-values

        # DROP ROWS WITH NAN VALUES
        orderPixelTable.dropna(axis='index', how='any',
                               subset=['cont_x'], inplace=True)
        orderPixelTable.dropna(axis='index', how='any',
                               subset=['cont_y'], inplace=True)

        foundLines = len(orderPixelTable.index)
        percent = 100 * foundLines / allLines

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "TSAMP",
            "qc_value": allLines,
            "qc_comment": "Total number of samples along orders",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "NSAMP",
            "qc_value": foundLines,
            "qc_comment": "Number of samples where a continuum is detected",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        self.psamp = percent / 100
        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "PSAMP",
            "qc_value": f"{self.psamp:0.3f}",
            "qc_comment": "Proportion of samples where a continuum is detected",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        self.log.print(f"\tContinuum found in {foundLines} out of {allLines} order slices ({percent:2.0f}%)")

        self.log.debug('completed the ``sample_trace`` method')
        return orderPixelTable

    # use the tab-trigger below for new method
    # xt-class-method
