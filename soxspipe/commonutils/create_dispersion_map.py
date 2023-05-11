#!/usr/bin/env python
# encoding: utf-8
"""
*detect arc-lines on a pinhole frame to generate a dispersion solution*

:Author:
    Marco Landoni & David Young

:Date Created:
    September  1, 2020
"""
################# GLOBAL IMPORTS ####################
from fundamentals import fmultiprocess
from tabulate import tabulate
from soxspipe.commonutils.toolkit import unpack_order_table, read_spectral_format
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials
import warnings
import math
from soxspipe.commonutils.toolkit import get_calibrations_path
from os.path import expanduser
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from fundamentals import tools
from builtins import object
import sys
import os
from io import StringIO
from contextlib import suppress
import copy
from datetime import datetime

os.environ['TERM'] = 'vt100'


class create_dispersion_map(object):
    """
    *detect arc-lines on a pinhole frame to generate a dispersion solution*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``pinholeFrame`` -- the calibrated pinhole frame (single or multi)
        - ``firstGuessMap`` -- the first guess dispersion map from the `soxs_disp_solution` recipe (needed in `soxs_spat_solution` recipe). Default *False*.
        - ``orderTable`` -- the order geometry table
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``sofName`` -- name of the originating SOF file
        - ``create2DMap`` -- create the 2D image map of wavelength, slit-position and order from disp solution.

    **Usage:**

    ```python
    from soxspipe.commonutils import create_dispersion_map
    mapPath, mapImagePath, res_plots, qcTable, productsTable = create_dispersion_map(
        log=log,
        settings=settings,
        pinholeFrame=frame,
        firstGuessMap=False,
        qcTable=self.qc,
        productsTable=self.products
    ).get()
    ```
    """

    def __init__(
            self,
            log,
            settings,
            pinholeFrame,
            firstGuessMap=False,
            orderTable=False,
            qcTable=False,
            productsTable=False,
            sofName=False,
            create2DMap=True
    ):
        self.log = log
        log.debug("instantiating a new 'create_dispersion_map' object")

        from photutils.utils import NoDetectionsWarning

        self.settings = settings
        self.pinholeFrame = pinholeFrame
        self.firstGuessMap = firstGuessMap
        self.orderTable = orderTable
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.create2DMap = create2DMap

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        self.kw = kw
        self.arm = pinholeFrame.header[kw("SEQ_ARM")]
        self.dateObs = pinholeFrame.header[kw("DATE_OBS")]
        self.inst = pinholeFrame.header[kw("INSTRUME")]

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            self.recipeName = "soxs-spatial-solution"
            self.recipeSettings = self.settings["soxs-spatial-solution"]
        else:
            self.recipeName = "soxs-disp-solution"
            self.recipeSettings = self.settings["soxs-disp-solution"]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        warnings.simplefilter('ignore', NoDetectionsWarning)

        # SET IMAGE ORIENTATION
        if self.inst == "SOXS":
            self.axisA = "y"
            self.axisB = "x"
        elif self.inst == "XSHOOTER":
            self.axisA = "x"
            self.axisB = "y"

        return None

    def get(self):
        """
        *generate the dispersion map*

        **Return:**
            - ``mapPath`` -- path to the file containing the coefficients of the x,y polynomials of the global dispersion map fit
        """
        self.log.debug('starting the ``get`` method')

        import pandas as pd

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            slitDeg = self.recipeSettings["slit-deg"]
        else:
            slitDeg = 0

        # READ PREDICTED LINE POSITIONS FROM FILE - RETURNED AS DATAFRAME
        orderPixelTable = self.get_predicted_line_list()
        totalLines = len(orderPixelTable.index)

        # GET THE WINDOW SIZE FOR ATTEMPTING TO DETECT LINES ON FRAME
        windowSize = self.recipeSettings["pixel-window-size"]
        self.windowHalf = int(windowSize / 2)

        # DETECT THE LINES ON THE PINHOLE FRAME AND
        # ADD OBSERVED LINES TO DATAFRAME
        iteration = 0
        print(f"FINDING PINHOLE ARC-LINES ON IMAGE")
        while iteration < 2:
            iteration += 1
            orderPixelTable = orderPixelTable.apply(
                self.detect_pinhole_arc_line, axis=1)

            if 'detector_x_shifted' in orderPixelTable.columns:
                orderPixelTable['x_diff'] = orderPixelTable['detector_x_shifted'] - orderPixelTable['observed_x']
                orderPixelTable['y_diff'] = orderPixelTable['detector_y_shifted'] - orderPixelTable['observed_y']
                orderPixelTable['detector_x_shifted'] = orderPixelTable['detector_x_shifted'] - orderPixelTable['x_diff'].mean()
                orderPixelTable['detector_y_shifted'] = orderPixelTable['detector_y_shifted'] - orderPixelTable['y_diff'].mean()
            else:
                orderPixelTable['x_diff'] = orderPixelTable['detector_x'] - orderPixelTable['observed_x']
                orderPixelTable['y_diff'] = orderPixelTable['detector_y'] - orderPixelTable['observed_y']
                orderPixelTable['detector_x_shifted'] = orderPixelTable['detector_x'] - orderPixelTable['x_diff'].mean()
                orderPixelTable['detector_y_shifted'] = orderPixelTable['detector_y'] - orderPixelTable['y_diff'].mean()

            print(f"\t ITERATION {iteration}: Mean X Y difference between predicted and measured positions: {orderPixelTable['x_diff'].mean():0.3f},{orderPixelTable['y_diff'].mean():0.3f}")

        # COLLECT MISSING LINES
        mask = (orderPixelTable['observed_x'].isnull())
        missingLines = orderPixelTable.loc[mask]

        # GROUP RESULTS BY WAVELENGTH
        lineGroups = missingLines.groupby(['wavelength', 'order'])
        lineGroups = lineGroups.size().to_frame(name='count').reset_index()

        # CREATE THE LIST OF INCOMPLETE MULTIPINHOLE WAVELENGTHS & ORDER SETS TO DROP
        missingLineThreshold = 1
        mask = (lineGroups['count'] > missingLineThreshold)
        orderPixelTable["dropped"] = False
        lineGroups = lineGroups.loc[mask]
        setsToDrop = lineGroups[['wavelength', 'order']]
        s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["dropped"] = s["dropped"].values

        # DROP MISSING VALUES
        orderPixelTable.dropna(axis='index', how='any', subset=[
            'observed_x'], inplace=True)

        detectedLines = len(orderPixelTable.index)
        percentageDetectedLines = (float(detectedLines) / float(totalLines))
        percentageDetectedLines = float("{:.6f}".format(percentageDetectedLines))

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        if "DISP" in self.recipeName.upper():
            tag = "single"
        else:
            tag = "multi"

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "NLINE",
            "qc_value": detectedLines,
            "qc_comment": f"[lines] Number of lines detected in {tag} pinhole frame",
            "qc_unit": "lines",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "PLINE",
            "qc_value": percentageDetectedLines,
            "qc_comment": f"Proportion of lines detected in {tag} pinhole frame",
            "qc_unit": None,
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        orderDeg = self.recipeSettings["order-deg"]
        wavelengthDeg = self.recipeSettings["wavelength-deg"]

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        fitFound = False
        tryCount = 0
        while not fitFound and tryCount < 5:
            try:
                popt_x, popt_y, res_plots = self.fit_polynomials(
                    orderPixelTable=orderPixelTable,
                    wavelengthDeg=wavelengthDeg,
                    orderDeg=orderDeg,
                    slitDeg=slitDeg,
                    missingLines=missingLines
                )
                fitFound = True
            except Exception as e:
                degList = [wavelengthDeg, orderDeg, slitDeg]
                degList[degList.index(max(degList))] -= 1
                wavelengthDeg, orderDeg, slitDeg = degList
                print(f"Wavelength, Order and Slit fitting orders reduced to {wavelengthDeg}, {orderDeg}, {slitDeg} to try and achieve a dispersion solution.")
                tryCount += 1
                if tryCount == 5:
                    print(f"Could not converge on a good fit to the dispersion solution. Please check the quality of your data or adjust your fitting parameters.")
                    raise e

        # WRITE THE MAP TO FILE
        mapPath = self.write_map_to_file(
            popt_x, popt_y, orderDeg, wavelengthDeg, slitDeg)

        if self.firstGuessMap and self.orderTable and self.create2DMap:
            mapImagePath = self.map_to_image(dispersionMapPath=mapPath)
            return mapPath, mapImagePath, res_plots, self.qc, self.products

        self.log.debug('completed the ``get`` method')
        return mapPath, None, res_plots, self.qc, self.products

    def get_predicted_line_list(
            self):
        """*lift the predicted line list from the static calibrations*

        **Return:**
            - ``orderPixelTable`` -- a panda's data-frame containing wavelength,order,slit_index,slit_position,detector_x,detector_y
        """
        self.log.debug('starting the ``get_predicted_line_list`` method')

        from astropy.table import Table

        kw = self.kw
        pinholeFrame = self.pinholeFrame
        dp = self.detectorParams

        # WHICH TYPE OF PINHOLE FRAME DO WE HAVE - SINGLE OR MULTI
        if self.pinholeFrame.header[kw("DPR_TECH")] == "ECHELLE,PINHOLE":
            frameTech = "single"
        elif self.pinholeFrame.header[kw("DPR_TECH")] == "ECHELLE,MULTI-PINHOLE":
            frameTech = "multi"
        else:
            raise TypeError(
                "The input frame needs to be a calibrated single- or multi-pinhole arc lamp frame")

        # FIND THE APPROPRIATE PREDICTED LINE-LIST
        arm = self.arm
        if arm != "NIR" and kw('WIN_BINX') in pinholeFrame.header:
            binx = int(self.pinholeFrame.header[kw('WIN_BINX')])
            biny = int(self.pinholeFrame.header[kw('WIN_BINY')])
        else:
            binx = 1
            biny = 1

        # READ THE FILE
        home = expanduser("~")

        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        predictedLinesFile = calibrationRootPath + "/" + dp["predicted pinhole lines"][frameTech][f"{binx}x{biny}"]

        # LINE LIST TO PANDAS DATAFRAME
        dat = Table.read(predictedLinesFile, format='fits')
        orderPixelTable = dat.to_pandas()

        # FITS TO PYTHON INDEXING
        # PHOTUTILS CENTRE OF BOTTOM LEFT PIXEL IS (0,0) BUT FOR WCS IT IS (1,1)
        # AND FOR PYTHON IT IS ALSO (0,0)
        orderPixelTable["detector_x"] -= 1.0
        orderPixelTable["detector_y"] -= 1.0

        # RENAME ALL COLUMNS FOR CONSISTENCY
        listName = []
        listName[:] = [l if l else l for l in listName]
        orderPixelTable.columns = [d.lower() if d.lower() in [
            "order", "wavelength"] else d for d in orderPixelTable.columns]

        if not self.firstGuessMap:
            slitIndex = int(dp["mid_slit_index"])
            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (orderPixelTable['slit_index'] != slitIndex)
            orderPixelTable.drop(index=orderPixelTable[mask].index, inplace=True)

        # WANT TO DETERMINE SYSTEMATIC SHIFT IF FIRST GUESS SOLUTION PRESENT
        if self.firstGuessMap:
            # ADD SOME EXTRA COLUMNS TO DATAFRAME

            # FILTER THE PREDICTED LINES TO ONLY SLIT POSITION INCLUDED IN
            # SINGLE PINHOLE FRAMES
            slitIndex = int(dp["mid_slit_index"])

            # GET THE OBSERVED PIXELS VALUES
            orderPixelTable = dispersion_map_to_pixel_arrays(
                log=self.log,
                dispersionMapPath=self.firstGuessMap,
                orderPixelTable=orderPixelTable
            )

            # CREATE A COPY OF THE DATA-FRAME TO DETERMINE SHIFTS
            tmpList = orderPixelTable.copy()

            mask = (tmpList['slit_index'] == slitIndex)
            tmpList.loc[mask, 'shift_x'] = tmpList.loc[
                mask, 'detector_x'].values - tmpList.loc[mask, 'fit_x'].values
            tmpList.loc[mask, 'shift_y'] = tmpList.loc[
                mask, 'detector_y'].values - tmpList.loc[mask, 'fit_y'].values

            # MERGING SHIFTS INTO MAIN DATAFRAME
            tmpList = tmpList.loc[tmpList['shift_x'].notnull(
            ), ['wavelength', 'order', 'shift_x', 'shift_y']]
            orderPixelTable = orderPixelTable.merge(tmpList, on=[
                'wavelength', 'order'], how='outer')

            # DROP ROWS WITH MISSING SHIFTS
            orderPixelTable.dropna(axis='index', how='any', subset=[
                'shift_x'], inplace=True)

            # SHIFT DETECTOR LINE PIXEL POSITIONS BY SHIFTS
            # UPDATE FILTERED VALUES
            orderPixelTable.loc[
                :, 'detector_x'] -= orderPixelTable.loc[:, 'shift_x']
            orderPixelTable.loc[
                :, 'detector_y'] -= orderPixelTable.loc[:, 'shift_y']

            # DROP HELPER COLUMNS
            orderPixelTable.drop(columns=['fit_x', 'fit_y',
                                          'shift_x', 'shift_y'], inplace=True)

        self.log.debug('completed the ``get_predicted_line_list`` method')
        return orderPixelTable

    def detect_pinhole_arc_line(
            self,
            predictedLine):
        """*detect the observed position of an arc-line given the predicted pixel positions*

        **Key Arguments:**
            - ``predictedLine`` -- single predicted line coordinates from predicted line-list

        **Return:**
            - ``predictedLine`` -- the line with the observed pixel coordinates appended (if detected, otherwise nan)
        """
        self.log.debug('starting the ``detect_pinhole_arc_line`` method')

        import numpy as np
        from photutils import DAOStarFinder
        from astropy.stats import sigma_clipped_stats

        pinholeFrame = self.pinholeFrame
        windowHalf = self.windowHalf
        if 'detector_x_shifted' in predictedLine:
            x = predictedLine['detector_x_shifted']
            y = predictedLine['detector_y_shifted']
        else:
            x = predictedLine['detector_x']
            y = predictedLine['detector_y']

        # CLIP A STAMP FROM IMAGE AROUNDS PREDICTED POSITION
        xlow = int(np.max([x - windowHalf, 0]))
        xup = int(np.min([x + windowHalf, pinholeFrame.shape[1]]))
        ylow = int(np.max([y - windowHalf, 0]))
        yup = int(np.min([y + windowHalf, pinholeFrame.shape[0]]))
        stamp = pinholeFrame[ylow:yup, xlow:xup]

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        stamp = np.ma.array(stamp.data, mask=stamp.mask)
        # USE DAOStarFinder TO FIND LINES WITH 2D GUASSIAN FITTING
        mean, median, std = sigma_clipped_stats(stamp, sigma=3.0)

        try:
            daofind = DAOStarFinder(
                fwhm=2.5, threshold=1.3 * std, roundlo=-5.0, roundhi=5.0, sharplo=0.0, sharphi=2.0, exclude_border=True)
            # SUBTRACTING MEDIAN MAKES LITTLE TO NO DIFFERENCE
            # sources = daofind(stamp - median)
            sources = daofind(stamp)
        except Exception as e:
            sources = None

        if 1 == 0:
            import random
            ran = random.randint(1, 300)
            # if ran == 200:
            if ran:
                import matplotlib.pyplot as plt
                plt.clf()
                plt.imshow(stamp)

                # plt.show()
        old_resid = windowHalf * 4
        if sources:
            # FIND SOURCE CLOSEST TO CENTRE
            if len(sources) > 1:
                for source in sources:
                    tmp_x = source['xcentroid']
                    tmp_y = source['ycentroid']
                    new_resid = ((windowHalf - tmp_x)**2 +
                                 (windowHalf - tmp_y)**2)**0.5
                    if new_resid < old_resid:
                        observed_x = tmp_x + xlow
                        observed_y = tmp_y + ylow
                        old_resid = new_resid
            else:
                observed_x = sources[0]['xcentroid'] + xlow
                observed_y = sources[0]['ycentroid'] + ylow

            if 1 == 0 and ran:
                plt.scatter(0, 0, marker='x', s=30)
                plt.scatter(observed_x - xlow, observed_y -
                            ylow, s=30)
                plt.text(windowHalf - 2, windowHalf - 2, f"{observed_x-xlow:0.2f},{observed_y -ylow:0.2f}", fontsize=16, c="black", verticalalignment='bottom')
                plt.show()

        else:
            observed_x = np.nan
            observed_y = np.nan
        # plt.show()

        predictedLine['observed_x'] = observed_x
        predictedLine['observed_y'] = observed_y

        self.log.debug('completed the ``detect_pinhole_arc_line`` method')
        return predictedLine

    def write_map_to_file(
            self,
            xcoeff,
            ycoeff,
            orderDeg,
            wavelengthDeg,
            slitDeg):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**
            - ``xcoeff`` -- the x-coefficients
            - ``ycoeff`` -- the y-coefficients
            - ``orderDeg`` -- degree of the order fitting
            - ``wavelengthDeg`` -- degree of wavelength fitting
            - ``slitDeg`` -- degree of the slit fitting (False for single pinhole)

        **Return:**
            - ``disp_map_path`` -- path to the saved file
        """
        self.log.debug('starting the ``write_map_to_file`` method')

        import pandas as pd
        from astropy.table import Table
        from astropy.io import fits

        arm = self.arm
        kw = self.kw

        # SORT X COEFFICIENT OUTPUT TO WRITE TO FILE
        coeff_dict_x = {}
        coeff_dict_x["axis"] = "x"
        coeff_dict_x["order-deg"] = orderDeg
        coeff_dict_x["wavelength-deg"] = wavelengthDeg
        coeff_dict_x["slit-deg"] = slitDeg
        n_coeff = 0
        for i in range(0, orderDeg + 1):
            for j in range(0, wavelengthDeg + 1):
                for k in range(0, slitDeg + 1):
                    coeff_dict_x[f'c{i}{j}{k}'] = xcoeff[n_coeff]
                    n_coeff += 1

        # SORT Y COEFFICIENT OUTPUT TO WRITE TO FILE
        coeff_dict_y = {}
        coeff_dict_y["axis"] = "y"
        coeff_dict_y["order-deg"] = orderDeg
        coeff_dict_y["wavelength-deg"] = wavelengthDeg
        coeff_dict_y["slit-deg"] = slitDeg
        n_coeff = 0
        for i in range(0, orderDeg + 1):
            for j in range(0, wavelengthDeg + 1):
                for k in range(0, slitDeg + 1):
                    coeff_dict_y[f'c{i}{j}{k}'] = ycoeff[n_coeff]
                    n_coeff += 1

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home)
        outDir += f"product/{self.recipeName}/"
        outDir = outDir.replace("//", "/")
        # Recursively create missing directories
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        if not self.sofName:
            filename = filenamer(
                log=self.log,
                frame=self.pinholeFrame,
                settings=self.settings
            )
        else:
            filename = self.sofName + ".fits"

        header = copy.deepcopy(self.pinholeFrame.header)
        header.pop(kw("DPR_TECH"))
        header.pop(kw("DPR_CATG"))
        header.pop(kw("DPR_TYPE"))

        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError, LookupError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        if slitDeg == 0:
            # filename = filename.split("ARC")[0] + "DISP_MAP.fits"
            header[kw("PRO_TECH").upper()] = "ECHELLE,PINHOLE"
        else:
            # filename = filename.split("ARC")[0] + "2D_MAP.fits"
            header[kw("PRO_TECH").upper()] = "ECHELLE,MULTI-PINHOLE"
        filePath = f"{outDir}/{filename}"

        df = pd.DataFrame([coeff_dict_x, coeff_dict_y])
        t = Table.from_pandas(df)
        BinTableHDU = fits.table_to_hdu(t)

        header[kw("SEQ_ARM").upper()] = arm
        header[kw("PRO_TYPE").upper()] = "REDUCED"
        header[kw("PRO_CATG").upper()] = f"DISP_TAB_{arm}".upper()
        self.dispMapHeader = header

        # WRITE QCs TO HEADERS
        for n, v, c, h in zip(self.qc["qc_name"].values, self.qc["qc_value"].values, self.qc["qc_comment"].values, self.qc["to_header"].values):
            if h:
                header[f"ESO QC {n}".upper()] = (v, c)

        priHDU = fits.PrimaryHDU(header=header)

        hduList = fits.HDUList([priHDU, BinTableHDU])
        hduList.writeto(filePath, checksum=True, overwrite=True)

        self.log.debug('completed the ``write_map_to_file`` method')
        return filePath

    def calculate_residuals(
            self,
            orderPixelTable,
            xcoeff,
            ycoeff,
            orderDeg,
            wavelengthDeg,
            slitDeg,
            writeQCs=False):
        """*calculate residuals of the polynomial fits against the observed line positions*

        **Key Arguments:**

            - ``orderPixelTable`` -- the predicted line list as a data frame
            - ``xcoeff`` -- the x-coefficients
            - ``ycoeff`` -- the y-coefficients
            - ``orderDeg`` -- degree of the order fitting
            - ``wavelengthDeg`` -- degree of wavelength fitting
            - ``slitDeg`` -- degree of the slit fitting (False for single pinhole)
            - ``writeQCs`` -- write the QCs to dataframe? Default *False*

        **Return:**
            - ``residuals`` -- combined x-y residuals
            - ``mean`` -- the mean of the combine residuals
            - ``std`` -- the stdev of the combine residuals
            - ``median`` -- the median of the combine residuals
        """
        self.log.debug('starting the ``calculate_residuals`` method')

        import numpy as np
        import pandas as pd

        arm = self.arm

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # ADD EXPONENTS TO ORDERTABLE UP-FRONT
        if isinstance(orderDeg, list):
            orderDegx = orderDeg[0]
            orderDegy = orderDeg[1]
            wavelengthDegx = wavelengthDeg[0]
            wavelengthDegy = wavelengthDeg[1]
            slitDegx = slitDeg[0]
            slitDegy = slitDeg[1]
        else:
            orderDegx = orderDeg
            orderDegy = orderDeg
            wavelengthDegx = wavelengthDeg
            wavelengthDegy = wavelengthDeg
            slitDegx = slitDeg
            slitDegy = slitDeg

        polyx = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDegx, wavelengthDeg=wavelengthDegx, slitDeg=slitDegx, exponentsIncluded=True, axis="x").poly
        polyy = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDegy, wavelengthDeg=wavelengthDegy, slitDeg=slitDegy, exponentsIncluded=True, axis="y").poly

        # CALCULATE X & Y RESIDUALS BETWEEN OBSERVED LINE POSITIONS AND POLY
        # FITTED POSITIONS
        orderPixelTable["fit_x"] = polyx(orderPixelTable, *xcoeff)
        orderPixelTable["fit_y"] = polyy(orderPixelTable, *ycoeff)
        orderPixelTable["residuals_x"] = orderPixelTable[
            "fit_x"] - orderPixelTable["observed_x"]
        orderPixelTable["residuals_y"] = orderPixelTable[
            "fit_y"] - orderPixelTable["observed_y"]

        # CALCULATE COMBINED RESIDUALS AND STATS
        orderPixelTable["residuals_xy"] = np.sqrt(np.square(
            orderPixelTable["residuals_x"]) + np.square(orderPixelTable["residuals_y"]))
        combined_res_mean = np.mean(orderPixelTable["residuals_xy"])
        combined_res_std = np.std(orderPixelTable["residuals_xy"])
        combined_res_median = np.median(orderPixelTable["residuals_xy"])

        if writeQCs:
            absx = abs(orderPixelTable["residuals_x"])
            absy = abs(orderPixelTable["residuals_y"])
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMIN",
                "qc_value": absx.min(),
                "qc_comment": "[px] Minimum residual in dispersion solution fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESMAX",
                "qc_value": absx.max(),
                "qc_comment": "[px] Maximum residual in dispersion solution fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XRESRMS",
                "qc_value": absx.std(),
                "qc_comment": "[px] Std-dev of residual in dispersion solution fit along x-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "YRESMIN",
                "qc_value": absy.min(),
                "qc_comment": "[px] Minimum residual in dispersion solution fit along y-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "YRESMAX",
                "qc_value": absy.max(),
                "qc_comment": "[px] Maximum residual in dispersion solution fit along y-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "YRESRMS",
                "qc_value": absy.std(),
                "qc_comment": "[px] Std-dev of residual in dispersion solution fit along y-axis",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XYRESMIN",
                "qc_value": orderPixelTable["residuals_xy"].min(),
                "qc_comment": "[px] Minimum residual in dispersion solution fit (XY combined)",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XYRESMAX",
                "qc_value": orderPixelTable["residuals_xy"].max(),
                "qc_comment": "[px] Maximum residual in dispersion solution fit (XY combined)",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)
            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "XYRESRMS",
                "qc_value": orderPixelTable["residuals_xy"].std(),
                "qc_comment": "[px] Std-dev of residual in dispersion solution (XY combined)",
                "qc_unit": "pixels",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``calculate_residuals`` method')
        return combined_res_mean, combined_res_std, combined_res_median, orderPixelTable

    def fit_polynomials(
            self,
            orderPixelTable,
            wavelengthDeg,
            orderDeg,
            slitDeg,
            missingLines=False):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**
            - ``orderPixelTable`` -- data frame containing order, wavelengths, slit positions and observed pixel positions
            - ``wavelengthDeg`` -- degree of wavelength fitting
            - ``orderDeg`` -- degree of the order fitting
            - ``slitDeg`` -- degree of the slit fitting (0 for single pinhole)
            - ``missingLines`` -- lines not detected on the image

        **Return:**
            - ``xcoeff`` -- the x-coefficients post clipping
            - ``ycoeff`` -- the y-coefficients post clipping
            - ``res_plots`` -- plot of fit residuals
        """
        self.log.debug('starting the ``fit_polynomials`` method')

        # FIRST REMOVE DROPPED LINES FILTERED ROWS FROM DATA FRAME
        allClippedLines = []
        mask = (orderPixelTable['dropped'] == True)
        allClippedLines.append(orderPixelTable.loc[mask])
        orderPixelTable.drop(index=orderPixelTable[
                             mask].index, inplace=True)

        import numpy as np
        from astropy.stats import sigma_clip
        from scipy.optimize import curve_fit
        import pandas as pd

        arm = self.arm
        dp = self.detectorParams

        # XSH
        if self.inst == "XSHOOTER":
            rotateImage = 90
            flipImage = 1
        else:
            # SOXS
            rotateImage = 0
            flipImage = 0

        if self.firstGuessMap:
            recipe = "soxs-spatial-solution"
        else:
            recipe = "soxs-disp-solution"

        clippedCount = 1

        # ADD EXPONENTS TO ORDERTABLE UP-FRONT
        if isinstance(orderDeg, list):
            orderDegx = orderDeg[0]
            orderDegy = orderDeg[1]
            wavelengthDegx = wavelengthDeg[0]
            wavelengthDegy = wavelengthDeg[1]
            slitDegx = slitDeg[0]
            slitDegy = slitDeg[1]
        else:
            orderDegx = orderDeg
            orderDegy = orderDeg
            wavelengthDegx = wavelengthDeg
            wavelengthDegy = wavelengthDeg
            slitDegx = slitDeg
            slitDegy = slitDeg

        for i in range(0, orderDegx + 1):
            orderPixelTable[f"order_pow_x_{i}"] = orderPixelTable["order"].pow(i)
        for j in range(0, wavelengthDegx + 1):
            orderPixelTable[f"wavelength_pow_x_{j}"] = orderPixelTable["wavelength"].pow(j)
        for k in range(0, slitDegx + 1):
            orderPixelTable[f"slit_position_pow_x_{k}"] = orderPixelTable["slit_position"].pow(k)
        for i in range(0, orderDegy + 1):
            orderPixelTable[f"order_pow_y_{i}"] = orderPixelTable["order"].pow(i)
        for j in range(0, wavelengthDegy + 1):
            orderPixelTable[f"wavelength_pow_y_{j}"] = orderPixelTable["wavelength"].pow(j)
        for k in range(0, slitDegy + 1):
            orderPixelTable[f"slit_position_pow_y_{k}"] = orderPixelTable["slit_position"].pow(k)

        polyx = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDegx, wavelengthDeg=wavelengthDegx, slitDeg=slitDegx, exponentsIncluded=True, axis="x").poly
        polyy = chebyshev_order_wavelength_polynomials(
            log=self.log, orderDeg=orderDegy, wavelengthDeg=wavelengthDegy, slitDeg=slitDegy, exponentsIncluded=True, axis="y").poly

        clippingSigma = self.settings[
            recipe]["poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.settings[
            recipe]["poly-clipping-iteration-limit"]

        print("\n# FINDING DISPERSION SOLUTION\n")

        iteration = 0

        xcoeff = np.ones((orderDegx + 1) *
                         (wavelengthDegx + 1) * (slitDegx + 1))
        ycoeff = np.ones((orderDegy + 1) *
                         (wavelengthDegy + 1) * (slitDegy + 1))
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            observed_x = orderPixelTable["observed_x"].to_numpy()
            observed_y = orderPixelTable["observed_y"].to_numpy()
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            self.log.info("""curvefit x""" % locals())

            xcoeff, pcov_x = curve_fit(
                polyx, xdata=orderPixelTable, ydata=observed_x, p0=xcoeff)

            # NOW Y
            self.log.info("""curvefit y""" % locals())

            ycoeff, pcov_y = curve_fit(
                polyy, xdata=orderPixelTable, ydata=observed_y, p0=ycoeff)

            self.log.info("""calculate_residuals""" % locals())
            mean_res, std_res, median_res, orderPixelTable = self.calculate_residuals(
                orderPixelTable=orderPixelTable,
                xcoeff=xcoeff,
                ycoeff=ycoeff,
                orderDeg=orderDeg,
                wavelengthDeg=wavelengthDeg,
                slitDeg=slitDeg,
                writeQCs=False)

            # SIGMA-CLIP THE DATA
            self.log.info("""sigma_clip""" % locals())

            masked_residuals = sigma_clip(
                orderPixelTable["residuals_xy"], sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
            orderPixelTable["residuals_masked"] = masked_residuals.mask

            # COUNT THE CLIPPED LINES
            mask = (orderPixelTable['residuals_masked'] == True)
            allClippedLines.append(orderPixelTable.loc[mask])
            totalAllClippedLines = pd.concat(allClippedLines, ignore_index=True)

            # GROUP MASKED RESIDUALS BY WAVELENGTH
            lineGroups = totalAllClippedLines.groupby(['wavelength', 'order'])
            lineGroups = lineGroups.size().to_frame(name='count').reset_index()
            # CREATE THE LIST OF INCOMPLETE MULTIPINHOLE WAVELENGTHS & ORDER SETS TO DROP
            missingLineThreshold = 3
            mask = (lineGroups['count'] > missingLineThreshold)
            orderPixelTable["dropped"] = False
            lineGroups = lineGroups.loc[mask]
            setsToDrop = lineGroups[['wavelength', 'order']]
            s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
            s["dropped"] = False
            s.loc[(s["_merge"] == "both"), "dropped"] = True
            orderPixelTable["dropped"] = s["dropped"].values

            # ADD NEWLY CLIPPED LINES FROM SETS THEN FLAG AS MASKED
            mask = ((orderPixelTable['dropped'] == True) & (orderPixelTable['residuals_masked'] == False))
            allClippedLines.append(orderPixelTable.loc[mask])

            orderPixelTable.loc[(orderPixelTable['dropped'] == True), 'residuals_masked'] = True

            # RETURN BREAKDOWN OF COLUMN VALUE COUNT
            valCounts = orderPixelTable[
                'residuals_masked'].value_counts(normalize=False)
            if True in valCounts:
                clippedCount = valCounts[True]
            else:
                clippedCount = 0

            if iteration > 1:
                # Cursor up one line and clear line
                sys.stdout.write("\x1b[1A\x1b[2K")

            print(f'ITERATION {iteration:02d}: {clippedCount} arc lines where clipped in this iteration of fitting a global dispersion map')

            mask = (orderPixelTable['residuals_masked'] == True)
            orderPixelTable.drop(index=orderPixelTable[
                                 mask].index, inplace=True)

        if len(allClippedLines):
            allClippedLines = pd.concat(allClippedLines, ignore_index=True)

        mean_res, std_res, median_res, orderPixelTable = self.calculate_residuals(
            orderPixelTable=orderPixelTable,
            xcoeff=xcoeff,
            ycoeff=ycoeff,
            orderDeg=orderDeg,
            wavelengthDeg=wavelengthDeg,
            slitDeg=slitDeg,
            writeQCs=True)

        mean_x_res = abs(orderPixelTable["residuals_x"]).mean()
        mean_y_res = abs(orderPixelTable["residuals_y"]).mean()

        # a = plt.figure(figsize=(40, 15))
        import matplotlib.pyplot as plt
        if arm == "UVB" or self.settings["instrument"].lower() == "soxs":
            fig = plt.figure(figsize=(6, 13.5), constrained_layout=True)
            # CREATE THE GRID OF AXES
            gs = fig.add_gridspec(5, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:, 0:2])
            bottomright = fig.add_subplot(gs[4:, 2:])
        else:
            fig = plt.figure(figsize=(6, 11), constrained_layout=True)
            # CREATE THE GRID OF AXES
            gs = fig.add_gridspec(6, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:, 0:2])
            bottomright = fig.add_subplot(gs[4:, 2:])

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = self.pinholeFrame.data
        if self.axisA == "x":
            rotatedImg = np.rot90(rotatedImg, rotateImage / 90)
            rotatedImg = np.flipud(rotatedImg)
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        toprow.set_title(
            "observed arc-line positions (post-clipping)", fontsize=10)

        if isinstance(missingLines, pd.core.frame.DataFrame):
            toprow.scatter(missingLines[f"detector_{self.axisB}"], missingLines[f"detector_{self.axisA}"], marker='o', c='red', s=5, alpha=0.2, linewidths=0.5)
        if len(allClippedLines):
            toprow.scatter(allClippedLines[f"observed_{self.axisB}"], allClippedLines[f"observed_{self.axisA}"], marker='x', c='red', s=5, alpha=0.6, linewidths=0.5)
        toprow.scatter(orderPixelTable[f"observed_{self.axisB}"], orderPixelTable[f"observed_{self.axisA}"], marker='o', c='green', s=1, alpha=0.6)
        # SHOW MID-SLIT PINHOLE - CHECK WE ARE LATCHING ONTO THE CORRECT PINHOLE POSITION
        mask = (orderPixelTable['slit_index'] == int(dp["mid_slit_index"]))
        toprow.scatter(orderPixelTable.loc[mask][f"observed_{self.axisB}"], orderPixelTable.loc[mask][f"observed_{self.axisA}"], marker='o', c='green', s=5, alpha=0.6, linewidths=0.5)
        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)

        toprow.tick_params(axis='both', which='major', labelsize=9)

        if self.axisA == "x":
            toprow.invert_yaxis()

        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap='gray', alpha=0.5)
        midrow.set_title(
            "global dispersion solution", fontsize=10)

        xfit = orderPixelTable[f"fit_{self.axisA}"]
        midrow.scatter(orderPixelTable[f"observed_{self.axisB}"], orderPixelTable[f"observed_{self.axisA}"], marker='o', c='red', s=1, alpha=0.9)
        midrow.scatter(orderPixelTable[f"fit_{self.axisB}"],
                       orderPixelTable[f"fit_{self.axisA}"], marker='o', c='blue', s=orderPixelTable[f"residuals_xy"] * 30, alpha=0.6)

        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis='both', which='major', labelsize=9)
        if self.axisA == "x":
            midrow.invert_yaxis()

        # PLOT THE FINAL RESULTS:
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(orderPixelTable[f"residuals_{self.axisA}"], orderPixelTable[
            f"residuals_{self.axisB}"], alpha=0.1)
        bottomleft.set_xlabel(f'{self.axisA} residual (mean: {mean_x_res:2.2f} pix)')
        bottomleft.set_ylabel(f'{self.axisB} residual (mean: {mean_y_res:2.2f} pix)')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)

        from astropy.visualization import hist
        hist(orderPixelTable["residuals_xy"], bins='scott', ax=bottomright, histtype='stepfilled',
             alpha=0.7, density=True)
        bottomright.set_xlabel('xy residual')
        bottomright.tick_params(axis='both', which='major', labelsize=9)
        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        if self.firstGuessMap:
            fig.suptitle(f"residuals of global dispersion solution fitting - multi-pinhole\n{subtitle}", fontsize=12)
        else:
            fig.suptitle(f"residuals of global dispersion solution fitting - single pinhole\n{subtitle}", fontsize=12)

        utcnow = datetime.utcnow()
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # GET FILENAME FOR THE RESIDUAL PLOT
        if not self.sofName:
            res_plots = filenamer(
                log=self.log,
                frame=self.pinholeFrame,
                settings=self.settings
            )
            res_plots = res_plots.replace(".fits", ".pdf")
        else:
            res_plots = self.sofName + "_RESIDUALS.pdf"
        # plt.show()

        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + "/qc/pdf/"
        outDir = outDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIes
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        if self.firstGuessMap:
            filePath = f"{outDir}/{res_plots}"
        else:
            filePath = f"{outDir}/{res_plots}"
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "DISP_MAP_RES",
            "file_name": res_plots,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} dispersion solution QC plots",
            "file_path": filePath
        }).to_frame().T], ignore_index=True)

        plt.tight_layout()
        # plt.show()
        plt.savefig(filePath, dpi=720)

        self.log.debug('completed the ``fit_polynomials`` method')
        return xcoeff, ycoeff, res_plots

    def create_placeholder_images(
            self,
            order=False,
            plot=False,
            reverse=False):
        """*create CCDData objects as placeholders to host the 2D images of the wavelength and spatial solutions from dispersion solution map*

        **Key Arguments:**
            - ``order`` -- specific order to generate the placeholder pixels for. Inner-order pixels set to NaN, else set to 0. Default *False* (generate all orders)
            - ``plot`` -- generate plots of placeholder images (for debugging). Default *False*.
            - ``reverse`` -- Inner-order pixels set to 0, else set to NaN (reverse of default output).

        **Return:**
            - ``slitMap`` -- placeholder image to add pixel slit positions to
            - ``wlMap`` -- placeholder image to add pixel wavelength values to

        **Usage:**

        ```python
        slitMap, wlMap, orderMap = self._create_placeholder_images(order=order)
        ```
        """
        self.log.debug('starting the ``create_placeholder_images`` method')

        import numpy as np
        from astropy.nddata import CCDData

        kw = self.kw
        dp = self.detectorParams

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=0.)

        # CREATE THE IMAGE SAME SIZE AS DETECTOR - NAN INSIDE ORDERS, 0 OUTSIDE
        science_pixels = dp["science-pixels"]
        xlen = science_pixels["columns"]["end"] - \
            science_pixels["columns"]["start"]
        ylen = science_pixels["rows"]["end"] - science_pixels["rows"]["start"]
        xlen, ylen
        if reverse:
            seedArray = np.empty((ylen, xlen))
            seedArray[:] = np.nan
        else:
            seedArray = np.zeros((ylen, xlen))
        wlMap = CCDData(seedArray, unit="adu")
        orderMap = wlMap.copy()
        uniqueOrders = orderPixelTable['order'].unique()
        expandEdges = 0

        if self.axisA == "x":
            axisALen = wlMap.data.shape[1]
            axisBLen = wlMap.data.shape[0]
        else:
            axisALen = wlMap.data.shape[0]
            axisBLen = wlMap.data.shape[1]

        for o in uniqueOrders:
            if order and o != order:
                continue
            axisBcoord = orderPixelTable.loc[
                (orderPixelTable["order"] == o)][f"{self.axisB}coord"]
            axisACoord_edgeup = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                f"{self.axisA}coord_edgeup"] + expandEdges
            axisACoord_edgelow = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                f"{self.axisA}coord_edgelow"] - expandEdges
            axisACoord_edgelow, axisACoord_edgeup, axisBcoord = zip(*[(l, u, b) for l, u, b in zip(axisACoord_edgelow, axisACoord_edgeup, axisBcoord) if l > 0 and l < axisALen and u > 0 and u < axisALen and b > 0 and b < axisBLen])
            if reverse:
                for b, u, l in zip(axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)):
                    if self.axisA == "x":
                        wlMap.data[b, l:u] = 0
                        orderMap.data[b, l:u] = o
                    else:
                        wlMap.data[l:u, b] = 0
                        orderMap.data[l:u, b] = o
            else:
                for b, u, l in zip(axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)):
                    if self.axisA == "x":
                        wlMap.data[b, l:u] = np.NaN
                        orderMap.data[b, l:u] = np.NaN
                    else:
                        wlMap.data[l:u, b] = np.NaN
                        orderMap.data[l:u, b] = np.NaN

        # SLIT MAP PLACEHOLDER SAME AS WAVELENGTH MAP PLACEHOLDER
        slitMap = wlMap.copy()

        # PLOT CCDDATA OBJECT
        if plot:
            import matplotlib.pyplot as plt
            rotatedImg = np.rot90(slitMap.data, 1)
            std = np.nanstd(slitMap.data)
            mean = np.nanmean(slitMap.data)
            vmax = mean + 3 * std
            vmin = mean - 3 * std
            plt.figure(figsize=(12, 5))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                       cmap='gray', alpha=1)
            plt.colorbar()
            plt.xlabel(
                "y-axis", fontsize=10)
            plt.ylabel(
                "x-axis", fontsize=10)
            plt.show()

        self.log.debug('completed the ``create_placeholder_images`` method')
        return slitMap, wlMap, orderMap

    def map_to_image(
            self,
            dispersionMapPath):
        """*convert the dispersion map to images in the detector format showing pixel wavelength values and slit positions*

        **Key Arguments:**
            - ``dispersionMapPath`` -- path to the full dispersion map to convert to images

        **Return:**
            - ``dispersion_image_filePath`` -- path to the FITS image with an extension for wavelength values and another for slit positions

        **Usage:**

        ```python
        mapImagePath = self.map_to_image(dispersionMapPath=mapPath)
        ```
        """
        self.log.debug('starting the ``map_to_image`` method')

        from soxspipe.commonutils.combiner import Combiner
        import numpy as np
        from astropy.io import fits

        print("\n# CREATING 2D IMAGE MAP FROM DISPERSION SOLUTION\n\n")

        self.dispersionMapPath = dispersionMapPath
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        self.map_to_image_displacement_threshold = 0.0001
        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        combinedSlitImage = False
        combinedWlImage = False

        # DEFINE AN INPUT ARRAY
        inputArray = [(order, minWl, maxWl) for order, minWl,
                      maxWl in zip(orderNums, waveLengthMin, waveLengthMax)]

        # NUMPY CAN BE TRICKY WITH MP
        numThreads = '1'
        os.environ['OPENBLAS_NUM_THREADS'] = numThreads
        os.environ['OMP_NUM_THREADS'] = numThreads
        os.environ['BLAS_NUM_THREADS'] = numThreads
        results = fmultiprocess(log=self.log, function=self.order_to_image,
                                inputArray=inputArray, poolSize=6, timeout=3600, turnOffMP=False)
        del os.environ['OPENBLAS_NUM_THREADS']
        del os.environ['OMP_NUM_THREADS']
        del os.environ['BLAS_NUM_THREADS']

        slitImages = [r[0] for r in results]
        wlImages = [r[1] for r in results]

        slitMap, wlMap, orderMap = self.create_placeholder_images(reverse=True)

        combinedSlitImage = Combiner(slitImages)
        combinedSlitImage = combinedSlitImage.sum_combine()
        combinedWlImage = Combiner(wlImages)
        combinedWlImage = combinedWlImage.sum_combine()

        combinedWlImage.data += wlMap.data
        combinedSlitImage.data += wlMap.data

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home) + f"/product/{self.recipeName}"
        outDir = outDir.replace("//", "/")
        # Recursively create missing directories
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        #

        # GET THE EXTENSION (WITH DOT PREFIX)
        extension = os.path.splitext(dispersionMapPath)[1]
        filename = os.path.basename(
            dispersionMapPath).replace(extension, "_IMAGE.fits")

        dispersion_image_filePath = f"{outDir}/{filename}"
        # WRITE CCDDATA OBJECT TO FILE
        # from astropy.io import fits
        header = copy.deepcopy(self.dispMapHeader)
        header[kw("PRO_CATG")] = f"DISP_IMAGE_{arm}".upper()
        primary_hdu = fits.PrimaryHDU(combinedWlImage.data, header=header)
        primary_hdu.header['EXTNAME'] = 'WAVELENGTH'
        image_hdu = fits.ImageHDU(combinedSlitImage.data)
        image_hdu.header['EXTNAME'] = 'SLIT'
        image_hdu2 = fits.ImageHDU(orderMap.data)
        image_hdu2.header['EXTNAME'] = 'ORDER'
        hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2])
        hdul.writeto(dispersion_image_filePath, output_verify='exception',
                     overwrite=True, checksum=False)

        self.log.debug('completed the ``map_to_image`` method')
        return dispersion_image_filePath

    def order_to_image(
            self,
            orderInfo):
        """*convert a single order in the dispersion map to wavelength and slit position images*

        **Key Arguments:**
            - ``orderInfo`` -- tuple containing the order number to generate the images for, the minimum wavelength to consider (from format table) and maximum wavelength to consider (from format table).

        **Return:**
            - ``slitMap`` -- the slit map with order values filled
            - ``wlMap`` -- the wavelengths map with order values filled

        **Usage:**

        ```python
        slitMap, wlMap = self.order_to_image(order=order,minWl=minWl, maxWl=maxWl)
        ```
        """
        self.log.debug('starting the ``order_to_image`` method')

        import numpy as np
        from scipy.interpolate import griddata

        (order, minWl, maxWl) = orderInfo
        slitMap, wlMap, orderMap = self.create_placeholder_images(order=order)

        # FIRST GENERATE A WAVELENGTH SURFACE - FINE WL, CHUNKY SLIT-POSTION
        wlRange = maxWl - minWl
        if self.arm.lower == "nir":
            grid_res_wavelength = wlRange / 3500
        else:
            grid_res_wavelength = wlRange / 6000
        slitLength = self.detectorParams["slit_length"]
        grid_res_slit = slitLength / 100

        halfGrid = (slitLength / 2) * 1.2
        slitArray = np.arange(-halfGrid, halfGrid +
                              grid_res_slit, grid_res_slit)
        wlArray = np.arange(minWl - int(wlRange / 10), maxWl + int(wlRange / 10), grid_res_wavelength)
        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate(
            [np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
        # NOW THE BIG WAVELEGTH ARRAY
        bigWlArray = np.tile(wlArray, np.shape(slitArray)[0])

        iteration = 0
        remainingPixels = 1
        iterationLimit = 20
        remainingCount = 1
        while remainingPixels and remainingCount and iteration < iterationLimit:
            iteration += 1

            # GENERATE THE ORDERPIXEL TABLE FROM WL AND SLIT-POSTION GRID .. IF WITHIN THRESHOLD OF CENTRE OF DETECTOR PIXEL THEN INJECT INTO MAPS
            orderPixelTable, remainingCount = self.convert_and_fit(
                order=order, bigWlArray=bigWlArray, bigSlitArray=bigSlitArray, slitMap=slitMap, wlMap=wlMap, iteration=iteration, plots=False)

            if remainingCount < 3:
                break

            orderPixelTable = orderPixelTable.drop_duplicates(subset=['pixel_x', 'pixel_y', 'order'])
            train_wlx = orderPixelTable["fit_x"].values
            train_wly = orderPixelTable["fit_y"].values
            train_wl = orderPixelTable["wavelength"].values
            train_sp = orderPixelTable["slit_position"].values

            # USE CUBIC SPLINE NEIGHEST NEIGHBOUR TO SEED RESULTS
            bigWlArray = griddata((train_wlx, train_wly), train_wl, (orderPixelTable['pixel_x'].values, orderPixelTable['pixel_y'].values), method="cubic")
            bigSlitArray = griddata((train_wlx, train_wly), train_sp, (orderPixelTable['pixel_x'].values, orderPixelTable['pixel_y'].values), method="cubic")

        self.log.debug('completed the ``order_to_image`` method')
        return slitMap, wlMap

    def convert_and_fit(
            self,
            order,
            bigWlArray,
            bigSlitArray,
            slitMap,
            wlMap,
            iteration,
            plots=False):
        """*convert wavelength and slit position grids to pixels*

        **Key Arguments:**
            - ``order`` -- the order being considered
            - ``bigWlArray`` -- 1D array of all wavelengths to be converted
            - ``bigSlitArray`` -- 1D array of all split-positions to be converted (same length as `bigWlArray`)
            - ``slitMap`` -- place-holder image hosting fitted pixel slit-position values
            - ``wlMap`` -- place-holder image hosting fitted pixel wavelength values
            - ``iteration`` -- the iteration index (used for CL reporting)

        **Return:**
            - ``orderPixelTable`` -- dataframe containing unfitted pixel info
            - ``remainingCount`` -- number of remaining pixels in orderTable

        **Usage:**

        ```python
        orderPixelTable = self.convert_and_fit(
                order=order, bigWlArray=bigWlArray, bigSlitArray=bigSlitArray, slitMap=slitMap, wlMap=wlMap)
        ```
        """
        self.log.debug('starting the ``convert_and_fit`` method')

        import pandas as pd
        import numpy as np

        # CREATE PANDAS DATAFRAME WITH LARGE ARRAYS - ONE ROW PER
        # WAVELENGTH-SLIT GRID CELL
        myDict = {
            "order": np.ones(bigWlArray.shape[0]) * order,
            "wavelength": bigWlArray,
            "slit_position": bigSlitArray
        }
        orderPixelTable = pd.DataFrame(myDict)

        # GET DETECTOR PIXEL POSITIONS FOR ALL WAVELENGTH-SLIT GRID CELLS
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersionMapPath,
            orderPixelTable=orderPixelTable
        )

        # INTEGER PIXEL VALUES & FIT DISPLACEMENTS FROM PIXEL CENTRES
        orderPixelTable["pixel_x"] = np.round(orderPixelTable["fit_x"].values)
        orderPixelTable["pixel_y"] = np.round(orderPixelTable["fit_y"].values)
        orderPixelTable["residual_x"] = orderPixelTable[
            "fit_x"] - orderPixelTable["pixel_x"]
        orderPixelTable["residual_y"] = orderPixelTable[
            "fit_y"] - orderPixelTable["pixel_y"]
        orderPixelTable["residual_xy"] = np.sqrt(np.square(
            orderPixelTable["residual_x"]) + np.square(orderPixelTable["residual_y"]))

        # ADD A COUNT COLUMN FOR THE NUMBER OF SMALL SLIT/WL PIXELS FALLING IN
        # LARGE DETECTOR PIXELS
        count = orderPixelTable.groupby(
            ['pixel_x', 'pixel_y']).size().reset_index(name='count')
        orderPixelTable = pd.merge(orderPixelTable, count, how='left', left_on=[
            'pixel_x', 'pixel_y'], right_on=['pixel_x', 'pixel_y'])
        orderPixelTable = orderPixelTable.sort_values(
            ['order', 'pixel_x', 'pixel_y', 'residual_xy'])

        # FILTER TO WL/SLIT POSITION CLOSE ENOUGH TO CENTRE OF PIXEL
        mask = (orderPixelTable['residual_xy'] <
                self.map_to_image_displacement_threshold)
        # KEEP ONLY VALUES CLOSEST TO CENTRE OF PIXEL
        newPixelValue = orderPixelTable.loc[mask].drop_duplicates(
            subset=['pixel_x', 'pixel_y'], keep="first")
        # REMOVE PIXELS FOUND IN newPixelValue FROM orderPixelTable
        orderPixelTable = newPixelValue[['pixel_x', 'pixel_y']].merge(orderPixelTable, on=[
            'pixel_x', 'pixel_y'], how='right', indicator=True).query('_merge == "right_only"').drop(columns=['_merge'])

        remainingCount = orderPixelTable.drop_duplicates(
            subset=['pixel_x', 'pixel_y'], keep="first")

        # ADD FITTED PIXELS TO PLACE HOLDER IMAGES
        for xx, yy, wavelength, slit_position in zip(newPixelValue["pixel_x"].values.astype(int), newPixelValue["pixel_y"].values.astype(int), newPixelValue["wavelength"].values, newPixelValue["slit_position"].values):
            try:
                wlMap.data[yy, xx] = np.where(
                    np.isnan(wlMap.data[yy, xx]), wavelength, wlMap.data[yy, xx])
                slitMap.data[yy, xx] = np.where(
                    np.isnan(slitMap.data[yy, xx]), slit_position, slitMap.data[yy, xx])
            except (IndexError):
                # PIXELS OUTSIDE OF DETECTOR EDGES - IGNORE
                pass

        sys.stdout.write("\x1b[1A\x1b[2K")
        percentageFound = (1 - (np.count_nonzero(np.isnan(wlMap.data)) / np.count_nonzero(wlMap.data))) * 100
        print(f"ORDER {order:02d}, iteration {iteration:02d}. {percentageFound:0.2f}% order pixels now fitted.")
        # print(f"ORDER {order:02d}, iteration {iteration:02d}. {percentageFound:0.2f}% order pixels now fitted. Fit found for {len(newPixelValue.index)} new pixels, {len(remainingCount.index)} image pixel remain to be constrained ({np.count_nonzero(np.isnan(wlMap.data))} nans in place-holder image)")

        if plots:
            from matplotlib import cm
            import matplotlib.pyplot as plt
            # PLOT CCDDATA OBJECT
            rotatedImg = slitMap.data
            if self.axisA == "x":
                rotatedImg = np.rot90(rotatedImg, 1)
                rotatedImg = np.flipud(rotatedImg)
            std = np.nanstd(rotatedImg)
            mean = np.nanmean(rotatedImg)
            cmap = cm.gray
            cmap.set_bad(color='#ADD8E6')
            vmax = np.nanmax(rotatedImg)
            vmin = np.nanmin(rotatedImg)
            plt.figure(figsize=(24, 10))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                       cmap=cmap, alpha=1)
            if self.axisA == "x":
                plt.gca().invert_yaxis()
            plt.colorbar()
            plt.ylabel(f"{self.axisA}-axis", fontsize=12)
            plt.xlabel(f"{self.axisB}-axis", fontsize=12)
            plt.show()

        remainingCount = len(remainingCount.index)

        self.log.debug('completed the ``convert_and_fit`` method')
        return orderPixelTable, remainingCount
