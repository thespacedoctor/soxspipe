#!/usr/bin/env python
# encoding: utf-8
"""
*detect arc-lines on a pinhole frame to generate a dispersion solution*

Author
: Marco Landoni & David Young

Date Created
: September  1, 2020
"""
################# GLOBAL IMPORTS ####################
from fundamentals import fmultiprocess
from soxspipe.commonutils.toolkit import unpack_order_table, read_spectral_format
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import chebyshev_order_wavelength_polynomials

from soxspipe.commonutils.toolkit import get_calibrations_path
from os.path import expanduser
from soxspipe.commonutils import detector_lookup
from soxspipe.commonutils import keyword_lookup
from fundamentals import tools
from builtins import object
import sys
import os
from datetime import datetime

os.environ['TERM'] = 'vt100'


class create_dispersion_map(object):
    """
    *detect arc-lines on a pinhole frame to generate a dispersion solution*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``pinholeFrame`` -- the calibrated pinhole frame (single or multi)
    - ``firstGuessMap`` -- the first guess dispersion map from the `soxs_disp_solution` recipe (needed in `soxs_spat_solution` recipe). Default *False*.
    - ``orderTable`` -- the order geometry table
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products
    - ``sofName`` -- name of the originating SOF file
    - ``create2DMap`` -- create the 2D image map of wavelength, slit-position and order from disp solution.
        - ``lineDetectionTable`` -- the list of arc-lines detected on the pinhole frame (used only for pipeline tuning)

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
            recipeSettings,
            pinholeFrame,
            firstGuessMap=False,
            orderTable=False,
            qcTable=False,
            productsTable=False,
            sofName=False,
            create2DMap=True,
            lineDetectionTable=False
    ):
        self.log = log
        log.debug("instantiating a new 'create_dispersion_map' object")

        import warnings
        import copy

        self.settings = settings
        self.pinholeFrame = pinholeFrame
        self.firstGuessMap = firstGuessMap
        self.orderTable = orderTable
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.create2DMap = create2DMap
        self.recipeSettings = copy.deepcopy(recipeSettings)
        self.lineDetectionTable = lineDetectionTable

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
        self.exptime = pinholeFrame.header[kw("EXPTIME")]

        from soxspipe.commonutils.toolkit import get_calibration_lamp
        self.lamp = get_calibration_lamp(log=self.log, frame=pinholeFrame, kw=kw)

        if self.arm.upper() == "NIR" and not ("ar" in self.lamp.lower() and "hg" in self.lamp.lower()):
            raise Exception("wrong lamp")

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            self.recipeName = "soxs-spatial-solution"
        else:
            self.recipeName = "soxs-disp-solution"

        # DETECTOR PARAMETERS LOOKUP OBJECToff
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        from photutils.utils import NoDetectionsWarning
        warnings.simplefilter('ignore', NoDetectionsWarning)
        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging
        logging.getLogger().setLevel(logging.INFO + 5)

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        home = expanduser("~")
        self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
        self.qcDir = self.qcDir.replace("//", "/")
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.qcDir):
            os.makedirs(self.qcDir)

        return None

    def get(self):
        """
        *generate the dispersion map*

        **Return:**

        - ``mapPath`` -- path to the file containing the coefficients of the x,y polynomials of the global dispersion map fit
        """
        self.log.debug('starting the ``get`` method')

        import pandas as pd
        from astropy.table import Table
        import numpy as np
        from astropy.stats import sigma_clipped_stats
        from astropy.stats import sigma_clip

        bootstrap_dispersion_solution = self.settings["bootstrap_dispersion_solution"]

        orderDeg = self.recipeSettings["order-deg"]
        wavelengthDeg = self.recipeSettings["wavelength-deg"]

        # WHICH RECIPE ARE WE WORKING WITH?
        if self.firstGuessMap:
            slitDeg = self.recipeSettings["slit-deg"]
        else:
            if isinstance(orderDeg, list):
                slitDeg = [0, 0]
            else:
                slitDeg = 0

        # READ PREDICTED LINE POSITIONS FROM FILE - RETURNED AS DATAFRAME
        orderPixelTable = self.get_predicted_line_list()
        totalLines = len(orderPixelTable.index)
        self.uniqueSlitPos = orderPixelTable['slit_position'].unique()

        # GET THE WINDOW SIZE FOR ATTEMPTING TO DETECT LINES ON FRAME
        windowSize = self.recipeSettings["pixel-window-size"]

        # MASK THE PINHOLE FRAME
        pinholeFrame = self.pinholeFrame
        pinholeFrameMasked = np.ma.array(pinholeFrame.data, mask=pinholeFrame.mask)
        mean, median, std = sigma_clipped_stats(pinholeFrameMasked, sigma=5.0, stdfunc="mad_std", cenfunc="median", maxiters=3)

        # pinholeFrameMasked = pinholeFrameMasked - median
        self.pinholeFrameMasked = pinholeFrameMasked

        self.mean = mean
        self.std = std

        from soxspipe.commonutils.toolkit import quicklook_image
        quicklook_image(
            log=self.log, CCDObject=pinholeFrame, show=False, ext='data', stdWindow=3, title=False, surfacePlot=True)

        if False:
            import matplotlib.pyplot as plt
            from matplotlib.patches import Rectangle
            fig, ax = plt.subplots()
            ax.imshow(pinholeFrameMasked, cmap='gray', origin='lower', vmin=self.mean, vmax=self.mean + 5 * self.std)
            for index, row in orderPixelTable.iterrows():
                x = row['detector_x']
                y = row['detector_y']
                ax.add_patch(Rectangle((x - windowSize / 2, y - windowSize / 2), windowSize, windowSize, fill=None, edgecolor='red'))
            plt.show()

        boost = True
        while boost:

            # SORT BY COLUMN NAME
            orderPixelTable.sort_values(['wavelength'], inplace=True)

            boost = False

            if not isinstance(self.lineDetectionTable, bool):
                orderPixelTable = self.lineDetectionTable.copy()
            else:
                # DETECT THE LINES ON THE PINHOLE FRAME AND
                # ADD OBSERVED LINES TO DATAFRAME
                iteration = 0
                self.log.print(f"\n# FINDING PINHOLE ARC-LINES ON IMAGE")
                iraf = False
                tmpDF = orderPixelTable.copy()
                while iteration < 3:

                    if iteration == 0:
                        sigmaLimit = 20
                        self.windowHalf = round(windowSize * 1.5)
                    elif iteration == 1:
                        sigmaLimit = 10
                        self.windowHalf = windowSize
                    else:
                        self.windowHalf = round(windowSize / 2)
                        sigmaLimit = self.recipeSettings['pinhole-detection-thres-sigma']

                    orderPixelTable = orderPixelTable.apply(self.detect_pinhole_arc_line, axis=1, iraf=iraf, sigmaLimit=sigmaLimit, iteration=iteration)

                    iteration += 1

                    # AFTER FIRST ITERATION WE CAN NARROW DOWN THE STAMP WINDOW TO REDUCE THE RISK OF IDENTIFYING THE WRONG PINHOLE
                    # self.windowHalf = 5
                    iraf = True

                    if 'detector_x_shifted' not in orderPixelTable.columns:
                        orderPixelTable['detector_x_shifted'] = orderPixelTable['detector_x']
                        orderPixelTable['detector_y_shifted'] = orderPixelTable['detector_y']

                    orderPixelTable['x_diff'] = orderPixelTable['detector_x_shifted'] - orderPixelTable['observed_x']
                    orderPixelTable['y_diff'] = orderPixelTable['detector_y_shifted'] - orderPixelTable['observed_y']

                    # FIND THE GLOBAL SHIFT SEEN BETWEEN PREDICTED AND OBSERVED PIXEL POSITIONS IN EACH ORDER ... THE MOVE TO A NEW SET OF PREDICTED LINE POSITIONS USING THAT SHIFT
                    uniqueorders = orderPixelTable['order'].unique()
                    for o in uniqueorders:
                        mask = (orderPixelTable['order'] == o)

                        meanx, medianx, stdx = sigma_clipped_stats(orderPixelTable.loc[mask]['x_diff'], sigma=3.0, stdfunc="mad_std", cenfunc="median")
                        meany, mediany, stdy = sigma_clipped_stats(orderPixelTable.loc[mask]['y_diff'], sigma=3.0, stdfunc="mad_std", cenfunc="median")
                        if np.isnan(medianx) or np.isnan(mediany):
                            self.log.warning(f"Could not find any arc lines in order {o}.")
                            continue

                        orderPixelTable.loc[mask, 'detector_x_shifted'] = orderPixelTable.loc[mask]['detector_x_shifted'] - medianx
                        orderPixelTable.loc[mask, 'detector_y_shifted'] = orderPixelTable.loc[mask]['detector_y_shifted'] - mediany
                        orderPixelTable.loc[mask, 'xy_diff'] = np.sqrt(np.square(orderPixelTable.loc[mask]["x_diff"]) + np.square(orderPixelTable.loc[mask]["y_diff"]))
                        sys.stdout.flush()
                        sys.stdout.write("\x1b[1A\x1b[2K")
                        self.log.print(f"\t ITERATION {iteration}: Median X Y difference between predicted and measured positions: {medianx:0.5f},{mediany:0.5f} (order {o})")

            # MAKE A CLEAN COPY OF THE DETECTION TABLE ... USED FOR PIPELINE TUNING ONLY
            lineDetectionTable = orderPixelTable.copy()

            # COLLECT MISSING LINES
            mask = (orderPixelTable['observed_x'].isnull())
            missingLines = orderPixelTable.loc[mask]
            # GROUP RESULTS BY WAVELENGTH
            lineGroups = missingLines.groupby(['wavelength', 'order'])
            lineGroups = lineGroups.size().to_frame(name='count').reset_index()

            # CREATE THE LIST OF INCOMPLETE MULTIPINHOLE WAVELENGTHS & ORDER SETS TO DROP
            orderPixelTable["dropped"] = False
            if "SPAT" in self.recipeName.upper():
                missingLineThreshold = 9 - self.recipeSettings['mph_line_set_min']
                mask = (lineGroups['count'] > missingLineThreshold)
                lineGroups = lineGroups.loc[mask]
                setsToDrop = lineGroups[['wavelength', 'order']]
                s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
                s["dropped"] = False
                s.loc[(s["_merge"] == "both"), "dropped"] = True
                orderPixelTable["droppedOnMissing"] = s["dropped"].values
                orderPixelTable.loc[(orderPixelTable["droppedOnMissing"] == True), "dropped"] = True

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
                "qc_name": "TLINE",
                "qc_value": totalLines,
                "qc_comment": f"Total number of line in {tag} line-list",
                "qc_unit": "lines",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "NLINE",
                "qc_value": detectedLines,
                "qc_comment": f"Number of lines detected in {tag} pinhole frame",
                "qc_unit": "lines",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

            self.qc = pd.concat([self.qc, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "qc_name": "PLINE",
                "qc_value": percentageDetectedLines,
                "qc_comment": f"Proportion of input line-list lines detected on {tag} pinhole frame",
                "qc_unit": None,
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "to_header": True
            }).to_frame().T], ignore_index=True)

            # GROUP FOUND LINES INTO SETS AND CLIP ON MEAN XY SHIFT RESIDUALS
            if True:
                orderPixelTable = self._clip_on_measured_line_metrics(orderPixelTable)

            # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
            fitFound = False
            tryCount = 0
            while not fitFound and tryCount < 5:

                popt_x, popt_y, goodLinesTable, clippedLinesTable = self.fit_polynomials(
                    orderPixelTable=orderPixelTable,
                    wavelengthDeg=wavelengthDeg,
                    orderDeg=orderDeg,
                    slitDeg=slitDeg,
                    missingLines=missingLines
                )

                if isinstance(popt_x, np.ndarray) and isinstance(popt_y, np.ndarray):
                    fitFound = True
                else:
                    if self.settings["tune-pipeline"]:
                        raise ArithmeticError("Could not converge on a good fit to the dispersion solution. Please check the quality of your data or adjust your fitting parameters.")

                    if not isinstance(wavelengthDeg, list):
                        degList = [wavelengthDeg, orderDeg, slitDeg]
                        degList[degList.index(max(degList))] -= 1
                        wavelengthDeg, orderDeg, slitDeg = degList
                    elif popt_x == "xerror":
                        degList = [wavelengthDeg[0], orderDeg[0], slitDeg[0]]
                        degList[degList.index(max(degList))] -= 1
                        wavelengthDeg[0], orderDeg[0], slitDeg[0] = degList
                    elif popt_y == "yerror":
                        degList = [wavelengthDeg[1], orderDeg[1], slitDeg[1]]
                        degList[degList.index(max(degList))] -= 1
                        wavelengthDeg[1], orderDeg[1], slitDeg[1] = degList

                    self.recipeSettings["order-deg"] = orderDeg
                    self.recipeSettings["wavelength-deg"] = wavelengthDeg

                    # WHICH RECIPE ARE WE WORKING WITH?
                    if self.firstGuessMap:
                        self.recipeSettings["slit-deg"] = slitDeg

                    self.log.print(f"Wavelength, Order and Slit fitting orders reduced to {wavelengthDeg}, {orderDeg}, {slitDeg} to try and achieve a dispersion solution.")
                    tryCount += 1
                    if tryCount == 5:
                        self.log.error(f"Could not converge on a good fit to the dispersion solution. Please check the quality of your data or adjust your fitting parameters.")
                        raise ArithmeticError("Could not converge on a good fit to the dispersion solution. Please check the quality of your data or adjust your fitting parameters.")

            if bootstrap_dispersion_solution:
                # WRITE THE MAP TO FILE
                mapPath = self.write_map_to_file(
                    popt_x, popt_y, orderDeg, wavelengthDeg, slitDeg)
                orderPixelTable = self.create_new_static_line_list(dispersionMapPath=mapPath)

                boost = True
                bootstrap_dispersion_solution = False
                from tabulate import tabulate
                print(tabulate(orderPixelTable, headers='keys', tablefmt='psql'))

        # GET FILENAME FOR THE LINE LISTS
        if not self.sofName:
            filename = filenamer(
                log=self.log,
                frame=self.pinholeFrame,
                settings=self.settings
            )
            goodLinesFN = filename.replace(".fits", "_FITTED_LINES.fits")
            missingLinesFN = filename.replace(".fits", "_MISSED_LINES.fits")
        else:
            goodLinesFN = self.sofName + "_FITTED_LINES.fits"
            missingLinesFN = self.sofName + "_MISSED_LINES.fits"

        # WRITE CLIPPED LINE LIST TO FILE
        keepColumns = ['wavelength', 'order', 'slit_index', 'slit_position', 'detector_x', 'detector_y', 'observed_x', 'observed_y', 'x_diff', 'y_diff', 'fit_x', 'fit_y', 'residuals_x', 'residuals_y', 'residuals_xy', 'sigma_clipped', "sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux", 'fwhm_px', 'R', 'pixelScaleNm', 'detector_x_shifted', 'detector_y_shifted', ]
        if "ion" in goodLinesTable.columns:
            keepColumns.insert(0, 'ion')
        clippedLinesTable['sigma_clipped'] = True
        clippedLinesTable['R'] = np.nan
        clippedLinesTable['pixelScaleNm'] = np.nan

        self.CLINE = len(clippedLinesTable.index)
        self.qc = pd.concat([self.qc, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "qc_name": "CLINE",
            "qc_value": self.CLINE,
            "qc_comment": f"Total number of detected lines clipped during solution fitting",
            "qc_unit": "lines",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "to_header": True
        }).to_frame().T], ignore_index=True)

        if False:
            # SAVE OUT NEW STATIC LINE LIST
            from astropy.table import Table
            freshTable = goodLinesTable[["ion", 'wavelength', 'order', 'slit_index', 'slit_position', 'observed_x', 'observed_y']]
            freshTable.rename(columns={"observed_x": "detector_x", "observed_y": "detector_y"}, inplace=True)
            t = Table.from_pandas(freshTable)
            t.write("Hg_clean.fits", overwrite=True)

        if len(clippedLinesTable.index):
            goodAndClippedLines = pd.concat([clippedLinesTable[keepColumns], goodLinesTable[keepColumns]], ignore_index=True)
        else:
            goodAndClippedLines = goodLinesTable[keepColumns]
        goodLinesTable = goodLinesTable[keepColumns]

        # SORT BY COLUMN NAME
        goodAndClippedLines.sort_values(['order', 'wavelength', 'slit_index'], inplace=True)
        goodLinesTable.sort_values(['order', 'wavelength', 'slit_index'], inplace=True)

        # WRITE GOOD LINE LIST TO FILE
        t = Table.from_pandas(goodAndClippedLines)
        filePath = f"{self.qcDir}/{goodLinesFN}"
        if not self.settings["tune-pipeline"]:
            t.write(filePath, overwrite=True)
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "DISP_MAP_LINES",
            "file_name": goodLinesFN,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} dispersion solution fitted lines",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        # WRITE MISSING LINE LIST TO FILE
        keepColumns = ['wavelength', 'order', 'slit_index', 'slit_position', 'detector_x', 'detector_y']
        # SORT BY COLUMN NAME
        missingLines.sort_values(['order', 'wavelength', 'slit_index'], inplace=True)
        t = Table.from_pandas(missingLines[keepColumns])
        filePath = f"{self.qcDir}/{missingLinesFN}"
        if not self.settings["tune-pipeline"]:
            t.write(filePath, overwrite=True)
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "DISP_MAP_LINES_MISSING",
            "file_name": missingLinesFN,
            "file_type": "FITS",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} undetected arc lines",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        # WRITE THE MAP TO FILE
        if not self.settings["tune-pipeline"]:
            mapPath = self.write_map_to_file(
                popt_x, popt_y, orderDeg, wavelengthDeg, slitDeg)
        else:
            mapPath = None
            self.create2DMap = False

        if self.firstGuessMap and self.orderTable and self.create2DMap:
            mapImagePath = self.map_to_image(dispersionMapPath=mapPath, orders=list(goodLinesTable['order'].unique()))
            res_plots = self._create_dispersion_map_qc_plot(
                xcoeff=popt_x,
                ycoeff=popt_y,
                orderDeg=orderDeg,
                wavelengthDeg=wavelengthDeg,
                slitDeg=slitDeg,
                orderPixelTable=goodLinesTable,
                missingLines=missingLines,
                allClippedLines=clippedLinesTable,
                dispMap=mapPath,
                dispMapImage=mapImagePath
            )
            return mapPath, mapImagePath, res_plots, self.qc, self.products, lineDetectionTable

        res_plots = self._create_dispersion_map_qc_plot(
            xcoeff=popt_x,
            ycoeff=popt_y,
            orderDeg=orderDeg,
            wavelengthDeg=wavelengthDeg,
            slitDeg=slitDeg,
            orderPixelTable=goodLinesTable,
            missingLines=missingLines,
            allClippedLines=clippedLinesTable
        )

        self.log.debug('completed the ``get`` method')
        return mapPath, None, res_plots, self.qc, self.products, lineDetectionTable

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

        # RENAME ALL COLUMNS FOR CONSISTENCY
        listName = []
        listName[:] = [l if l else l for l in listName]
        orderPixelTable.columns = [d.lower() if d.lower() in [
            "order", "wavelength"] else d for d in orderPixelTable.columns]

        floats = ["wavelength", 'slit_position', 'detector_x', 'detector_y']
        ints = ['order', 'slit_index']
        strings = ['ion']
        for f in floats:
            try:
                orderPixelTable[f] = orderPixelTable[f].astype(float)
            except:
                pass
        for i in ints:
            try:
                orderPixelTable[i] = orderPixelTable[i].astype(int)
            except:
                pass
        for s in strings:
            try:
                orderPixelTable[s] = orderPixelTable[s].astype(str)
            except:
                pass

        # THERE ARE DUPLICATES IN XSHOOTER LINE LIST
        orderPixelTable.drop_duplicates(inplace=True)

        # FITS TO PYTHON INDEXING
        # PHOTUTILS CENTRE OF BOTTOM LEFT PIXEL IS (0,0) BUT FOR WCS IT IS (1,1)
        # AND FOR PYTHON IT IS ALSO (0,0)
        if self.inst != "SOXS":
            orderPixelTable["detector_x"] -= 1.0
            orderPixelTable["detector_y"] -= 1.0
        elif True:

            if arm == "NIR" or arm == "VIS":
                dp = self.detectorParams
                science_pixels = dp["science-pixels"]
                xlen = science_pixels["columns"]["end"] - \
                    science_pixels["columns"]["start"]
                ylen = science_pixels["rows"]["end"] - science_pixels["rows"]["start"]
            else:
                orderPixelTable["detector_x"] -= science_pixels["columns"]["start"]
                orderPixelTable["detector_y"] -= science_pixels["rows"]["start"]
            if arm == "NIR":
                pass
                # THIS IS TEMPORARALLY NEEDED TO ADJUST OLD LINE LIST
                # orderPixelTable["detector_x"] += 7
                # orderPixelTable["detector_y"] -= 9
                # mask = (orderPixelTable['order'] == 24)
                # orderPixelTable.loc[mask, "detector_x"] += 1
                # orderPixelTable.loc[mask, "detector_y"] += 1
                # mask = (orderPixelTable['order'] == 23)
                # orderPixelTable.loc[mask, "detector_x"] += 1
                # orderPixelTable.loc[mask, "detector_y"] += 1
                # mask = (orderPixelTable['order'] == 22)
                # orderPixelTable.loc[mask, "detector_x"] += 1
                # orderPixelTable.loc[mask, "detector_y"] += 1
                # mask = (orderPixelTable['order'] == 21)
                # orderPixelTable.loc[mask, "detector_x"] += 1
                # orderPixelTable.loc[mask, "detector_y"] += 2
                # mask = (orderPixelTable['order'] == 20)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 2
                # mask = (orderPixelTable['order'] == 19)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 2
                # mask = (orderPixelTable['order'] == 18)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 2
                # mask = (orderPixelTable['order'] == 17)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 16)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 15)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 14)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 13)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 12)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 11)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 10)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # mask = (orderPixelTable['order'] == 9)
                # orderPixelTable.loc[mask, "detector_x"] += 2
                # orderPixelTable.loc[mask, "detector_y"] += 3
                # # mask = (orderPixelTable['order'] == 2)
                # # orderPixelTable.loc[mask, "detector_x"] -= 23
                # # orderPixelTable.loc[mask, "detector_y"] -= 0
                # # mask = (orderPixelTable['order'] == 3)
                # # orderPixelTable.loc[mask, "detector_x"] -= 25
                # # orderPixelTable.loc[mask, "detector_y"] -= 0
                # # mask = (orderPixelTable['order'] == 4)
                # # orderPixelTable.loc[mask, "detector_x"] -= 10
                # # orderPixelTable.loc[mask, "detector_y"] -= 2
                # from astropy.table import Table
                # t = Table.from_pandas(orderPixelTable)
                # # t.write("Xe.fits", overwrite=True)
                # t.write("ArHgNeXe_clean_within_2.0pixel.fits", overwrite=True)

        if not self.firstGuessMap:
            slitIndex = int(dp["mid_slit_index"])
            # REMOVE FILTERED ROWS FROM DATA FRAME
            mask = (orderPixelTable['slit_index'] == slitIndex)
            orderPixelTable = orderPixelTable.loc[mask]

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

            # DROP MPH SETS THAT DON'T CONTAIN AlL SLIT-POSITIONS
            lineGroups = orderPixelTable.groupby(['wavelength', 'order'])
            lineGroups = lineGroups.size().to_frame(name='count').reset_index()
            fullSet = lineGroups['count'].max()
            mask = (lineGroups['count'] < fullSet)
            lineGroups = lineGroups.loc[mask]
            setsToDrop = lineGroups[['wavelength', 'order']]
            s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
            s["dropped"] = False
            s.loc[(s["_merge"] == "both"), "dropped"] = True
            orderPixelTable["droppedOnMissing"] = s["dropped"].values
            orderPixelTable = orderPixelTable.loc[(orderPixelTable["droppedOnMissing"] == False)]
            orderPixelTable.drop(columns=['droppedOnMissing'], inplace=True)

        self.log.debug('completed the ``get_predicted_line_list`` method')
        return orderPixelTable

    def detect_pinhole_arc_line(
            self,
            predictedLine,
            iraf=True,
            sigmaLimit=3,
            iteration=False):
        """*detect the observed position of an arc-line given the predicted pixel positions*

        **Key Arguments:**

        - ``predictedLine`` -- single predicted line coordinates from predicted line-list
        - ``iraf`` -- use IRAF star finder to generate a FWHM
        - ``sigmaLimit`` -- the lower sigma limit for arc line to be considered detected
        - ``iteration`` -- which detect and shift iteration are we on?

        **Return:**

        - ``predictedLine`` -- the line with the observed pixel coordinates appended (if detected, otherwise nan)
        """
        self.log.debug('starting the ``detect_pinhole_arc_line`` method')

        import numpy as np
        from photutils import DAOStarFinder, IRAFStarFinder
        from astropy.stats import sigma_clipped_stats
        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging
        logging.getLogger().setLevel(logging.INFO + 5)

        if iteration is False:
            iterationText = ""
        else:
            iterationText = f", iter #{iteration}"

        pinholeFrame = self.pinholeFrameMasked
        windowHalf = self.windowHalf
        if 'detector_x_shifted' in predictedLine:
            x = predictedLine['detector_x_shifted']
            y = predictedLine['detector_y_shifted']
        else:
            x = predictedLine['detector_x']
            y = predictedLine['detector_y']

        order = predictedLine['order']
        wl = predictedLine['wavelength']

        # CLIP A STAMP FROM IMAGE AROUNDS PREDICTED POSITION
        xlow = round(np.max([x - windowHalf, 0]))
        xup = round(np.min([x + windowHalf, pinholeFrame.shape[1]]))
        ylow = round(np.max([y - windowHalf, 0]))
        yup = round(np.min([y + windowHalf, pinholeFrame.shape[0]]))
        stamp = pinholeFrame[ylow:yup, xlow:xup]

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        stamp = np.ma.array(stamp.data, mask=stamp.mask)
        # USE DAOStarFinder TO FIND LINES WITH 2D GAUSSIAN FITTING
        mean, median, std = sigma_clipped_stats(stamp, sigma=3.0, stdfunc="mad_std", cenfunc="median")

        try:
            # LET AS MANY LINES BE DETECTED AS POSSIBLE ... WE WILL CLEAN UP LATER
            daofind = DAOStarFinder(
                fwhm=3., threshold=sigmaLimit * std, roundlo=-2.0, roundhi=2.0, sharplo=-1, sharphi=3.0, exclude_border=False)
            # SUBTRACTING MEDIAN MAKES LITTLE TO NO DIFFERENCE .. EXCEPT FOR LOW SIGNAL IMAGES
            # sources = daofind(stamp.data, mask=stamp.mask)
            sources = daofind(stamp.data - median, mask=stamp.mask)
        except Exception as e:
            sources = None

        old_resid = windowHalf * 4
        selectedSource = None
        observed_x = np.nan
        observed_y = np.nan
        fwhm = np.nan

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
                        stamp_x = tmp_x
                        stamp_y = tmp_y
                        old_resid = new_resid
                        selectedSource = source
                        detectionSigma = source['peak'] / std

            else:
                observed_x = sources[0]['xcentroid'] + xlow
                observed_y = sources[0]['ycentroid'] + ylow
                stamp_x = sources[0]['xcentroid']
                stamp_y = sources[0]['ycentroid']
                selectedSource = sources[0]
                detectionSigma = sources[0]['peak'] / std

            if iraf:
                try:
                    # Rerun detection with IRAFStarFinder
                    iraf_find = IRAFStarFinder(
                        fwhm=3., threshold=1 * std, roundlo=-2.0, roundhi=2.0, sharplo=-1, sharphi=3.0, exclude_border=True, xycoords=[(stamp_x, stamp_y)])
                    iraf_sources = iraf_find(stamp)
                    fwhm = iraf_sources['fwhm'][0]
                except Exception as e:
                    fwhm = np.nan
                    # print(np.shape(stamp))
                    # print(stamp_x,stamp_y)
                    # print(e)
                    pass

            if False and iteration == 2:
                import random
                ran = random.randint(1, 3)
                # if int(wl) == 7125:
                # if True or order == 2.:
                if True:
                    print(x, y)
                    # if ran == 1:
                    import matplotlib.pyplot as plt
                    plt.clf()
                    plt.imshow(stamp)
                    # plt.show()
                    plt.scatter(0, 0, marker='x', s=30)
                    plt.scatter(observed_x - xlow, observed_y -
                                ylow, s=30)
                    plt.text(windowHalf - 2, windowHalf - 2, f"{observed_x-xlow:0.2f},{observed_y -ylow:0.2f}", fontsize=16, c="black", verticalalignment='bottom')
                    plt.text(2, 2, f"{sigmaLimit}$\\sigma$, order {order}{iterationText}\n{wl:0.1f},{detectionSigma:0.1f}", fontsize=16, c="white", verticalalignment='bottom')
                    plt.show()

        # plt.show()

        keepValues = ["sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux"]
        if not selectedSource:
            for k in keepValues:
                predictedLine[k] = np.nan
        else:
            for k in keepValues:
                predictedLine[k] = selectedSource[k]

        predictedLine['observed_x'] = observed_x
        predictedLine['observed_y'] = observed_y
        if iraf:
            predictedLine['fwhm_px'] = fwhm

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
        from contextlib import suppress
        import copy
        import math
        arm = self.arm
        kw = self.kw

        # SORT X COEFFICIENT OUTPUT TO WRITE TO FILE
        if isinstance(orderDeg, list):
            coeff_dict_x = {}
            coeff_dict_x["axis"] = "x"
            coeff_dict_x["order_deg"] = orderDeg[0]
            coeff_dict_x["wavelength_deg"] = wavelengthDeg[0]
            coeff_dict_x["slit_deg"] = slitDeg[0]
            n_coeff = 0
            for i in range(0, orderDeg[0] + 1):
                for j in range(0, wavelengthDeg[0] + 1):
                    for k in range(0, slitDeg[0] + 1):
                        coeff_dict_x[f'c{i}{j}{k}'] = xcoeff[n_coeff]
                        n_coeff += 1
        else:
            coeff_dict_x = {}
            coeff_dict_x["axis"] = "x"
            coeff_dict_x["order_deg"] = orderDeg
            coeff_dict_x["wavelength_deg"] = wavelengthDeg
            coeff_dict_x["slit_deg"] = slitDeg
            n_coeff = 0
            for i in range(0, orderDeg + 1):
                for j in range(0, wavelengthDeg + 1):
                    for k in range(0, slitDeg + 1):
                        coeff_dict_x[f'c{i}{j}{k}'] = xcoeff[n_coeff]
                        n_coeff += 1

        # SORT Y COEFFICIENT OUTPUT TO WRITE TO FILE
        if isinstance(orderDeg, list):
            coeff_dict_y = {}
            coeff_dict_y["axis"] = "y"
            coeff_dict_y["order_deg"] = orderDeg[1]
            coeff_dict_y["wavelength_deg"] = wavelengthDeg[1]
            coeff_dict_y["slit_deg"] = slitDeg[1]
            n_coeff = 0
            for i in range(0, orderDeg[1] + 1):
                for j in range(0, wavelengthDeg[1] + 1):
                    for k in range(0, slitDeg[1] + 1):
                        coeff_dict_y[f'c{i}{j}{k}'] = ycoeff[n_coeff]
                        n_coeff += 1
        else:
            coeff_dict_y = {}
            coeff_dict_y["axis"] = "y"
            coeff_dict_y["order_deg"] = orderDeg
            coeff_dict_y["wavelength_deg"] = wavelengthDeg
            coeff_dict_y["slit_deg"] = slitDeg
            n_coeff = 0
            for i in range(0, orderDeg + 1):
                for j in range(0, wavelengthDeg + 1):
                    for k in range(0, slitDeg + 1):
                        coeff_dict_y[f'c{i}{j}{k}'] = ycoeff[n_coeff]
                        n_coeff += 1

        # DETERMINE WHERE TO WRITE THE FILE
        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home)
        outDir += f"/product/{self.recipeName}/"
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

        if slitDeg == 0 or slitDeg == [0, 0]:
            # filename = filename.split("ARC")[0] + "DISP_MAP.fits"
            header[kw("PRO_TECH").upper()] = "ECHELLE,PINHOLE"
        else:
            # filename = filename.split("ARC")[0] + "2D_MAP.fits"
            header[kw("PRO_TECH").upper()] = "ECHELLE,MULTI-PINHOLE"
        filePath = f"{outDir}/{filename}"

        df = pd.DataFrame([coeff_dict_x, coeff_dict_y])
        t = Table.from_pandas(df)

        # SORT COLUMNS
        cols = list(t.columns)
        startList = ['axis', 'order_deg', 'wavelength_deg', 'slit_deg']
        cols = [c for c in cols if c not in startList]
        cols.sort()
        cols = startList + cols
        t = t[cols]

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

        cache = self.settings["workspace-root-dir"].replace("~", home) + "/.cache"
        # Recursively create missing directories
        if not os.path.exists(cache):
            os.makedirs(cache)
        polyOrders = [orderDeg, wavelengthDeg, slitDeg]
        if isinstance(orderDeg, list):
            merged_list = []
            for sublist in polyOrders:
                merged_list.extend(sublist)
            polyOrders = merged_list
        polyOrders[:] = [str(l) for l in polyOrders]
        polyOrders = "".join(polyOrders)
        filename = f"{self.recipeName}_{self.arm}_{polyOrders}.fits"
        cacheFilePath = f"{cache}/{filename}"
        hduList.writeto(cacheFilePath, checksum=True, overwrite=True)

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
            writeQCs=False,
            pixelRange=False):
        """*calculate residuals of the polynomial fits against the observed line positions*

        **Key Arguments:**

        - ``orderPixelTable`` -- the predicted line list as a data frame
        - ``xcoeff`` -- the x-coefficients
        - ``ycoeff`` -- the y-coefficients
        - ``orderDeg`` -- degree of the order fitting
        - ``wavelengthDeg`` -- degree of wavelength fitting
        - ``slitDeg`` -- degree of the slit fitting (False for single pinhole)
        - ``writeQCs`` -- write the QCs to dataframe? Default *False*
        - ``pixelRange`` -- return centre pixel *and* +- 2nm from the centre pixel (to measure the pixel scale)

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

        if pixelRange == True:
            polyx = chebyshev_order_wavelength_polynomials(
                log=self.log, orderDeg=orderDegx, wavelengthDeg=wavelengthDegx, slitDeg=slitDegx, exponentsIncluded=False, axis="x").poly
            polyy = chebyshev_order_wavelength_polynomials(
                log=self.log, orderDeg=orderDegy, wavelengthDeg=wavelengthDegy, slitDeg=slitDegy, exponentsIncluded=False, axis="y").poly
            # GET THE PIXEL SCALE
            orderPixelTableHigh = orderPixelTable.copy()
            nmRange = 4.
            orderPixelTableHigh["wavelength"] = orderPixelTableHigh["wavelength"] + nmRange / 2.
            orderPixelTableHigh["fit_x"] = polyx(orderPixelTableHigh, *xcoeff)
            orderPixelTableHigh["fit_y"] = polyy(orderPixelTableHigh, *ycoeff)

            orderPixelTableLow = orderPixelTable.copy()
            orderPixelTableLow["wavelength"] = orderPixelTableLow["wavelength"] - nmRange / 2.
            orderPixelTableLow["fit_x"] = polyx(orderPixelTableLow, *xcoeff)
            orderPixelTableLow["fit_y"] = polyy(orderPixelTableLow, *ycoeff)

            orderPixelTable["fit_x_high"] = orderPixelTableHigh["fit_x"]
            orderPixelTable["fit_y_high"] = orderPixelTableHigh["fit_y"]
            orderPixelTable["fit_x_low"] = orderPixelTableLow["fit_x"]
            orderPixelTable["fit_y_low"] = orderPixelTableLow["fit_y"]

            orderPixelTable["pixelScaleNm"] = nmRange / np.power(np.power(orderPixelTable["fit_x_high"] - orderPixelTable["fit_x_low"], 2) + np.power(orderPixelTable["fit_y_high"] - orderPixelTable["fit_y_low"], 2), 0.5)
            orderPixelTable["delta_wavelength"] = orderPixelTable["pixelScaleNm"] * orderPixelTable["fwhm_px"]
            orderPixelTable["R"] = orderPixelTable["wavelength"] / orderPixelTable["delta_wavelength"]

            # REMOVE COLUMN FROM DATA FRAME
            orderPixelTable.drop(columns=['fit_x_high', 'fit_y_high', 'fit_x_low', 'fit_y_low'], inplace=True)

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
        - ``goodLinesTable`` -- the fitted line-list with metrics
        - ``clippedLinesTable`` -- the lines that were sigma-clipped during polynomial fitting
        """
        self.log.debug('starting the ``fit_polynomials`` method')

        # FIRST REMOVE DROPPED LINES FILTERED ROWS FROM DATA FRAME
        allClippedLines = []
        mask = (orderPixelTable['dropped'] == True)
        allClippedLines.append(orderPixelTable.loc[mask])
        orderPixelTable = orderPixelTable.loc[~mask]

        import numpy as np
        from astropy.stats import sigma_clip
        from scipy.optimize import curve_fit
        import pandas as pd
        from soxspipe.commonutils import get_cached_coeffs

        arm = self.arm
        dp = self.detectorParams

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

        clippingSigma = self.recipeSettings["poly-fitting-residual-clipping-sigma"]
        clippingSigmaX = clippingSigma
        clippingSigmaY = clippingSigma
        clippingIterationLimit = self.recipeSettings["poly-clipping-iteration-limit"]
        if "poly-clipping-pinhole-sets" in self.recipeSettings and self.recipeSettings["poly-clipping-pinhole-sets"]:
            clipOnMphSets = True
        else:
            clipOnMphSets = False

        self.log.print("\n# FINDING DISPERSION SOLUTION")

        iteration = 0

        # CREATE A QUICK FIRST GUESS AT COEFFS ... SPEEDS UP FIRST ITERATION OF FULL SET
        tmpDF = orderPixelTable.copy()
        mean_res = 100.

        orderPixelTable['sigma_clipped'] = False
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            observed_x = orderPixelTable["observed_x"].to_numpy()
            observed_y = orderPixelTable["observed_y"].to_numpy()

            # IF mean_res < 10 WE WANT TO START FROM SCRATCH AGAIN SO NOT TO INFLUENCE THE FINAL RESULT
            if True and mean_res > 10:
                # FIND CACHED COEFF ELSE RETURN ARRAYS OF 1s
                xcoeff, ycoeff = get_cached_coeffs(
                    log=self.log,
                    arm=arm,
                    settings=self.settings,
                    recipeName=self.recipeName,
                    orderDeg=orderDeg,
                    wavelengthDeg=wavelengthDeg,
                    slitDeg=slitDeg
                )

            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            self.log.info("""curvefit x""" % locals())

            try:
                xcoeff, pcov_x = curve_fit(
                    polyx, xdata=orderPixelTable, ydata=observed_x, p0=xcoeff, maxfev=30000)
            except:
                return "xerror", None, None, None

            # NOW Y
            self.log.info("""curvefit y""" % locals())

            try:
                ycoeff, pcov_y = curve_fit(
                    polyy, xdata=orderPixelTable, ydata=observed_y, p0=ycoeff, maxfev=30000)
            except:
                return None, "yerror", None, None

            self.log.info("""calculate_residuals""" % locals())
            mean_res, std_res, median_res, orderPixelTable = self.calculate_residuals(
                orderPixelTable=orderPixelTable,
                xcoeff=xcoeff,
                ycoeff=ycoeff,
                orderDeg=orderDeg,
                wavelengthDeg=wavelengthDeg,
                slitDeg=slitDeg,
                pixelRange=True,
                writeQCs=False)

            # DO SOME CLIPPING ON THE PROFILES OF THE DETECTED LINES
            if iteration == -1:
                orderPixelTable = self._clip_on_measured_line_metrics(orderPixelTable)

            # COUNT THE CLIPPED LINES
            mask = (orderPixelTable['sigma_clipped'] == True)
            allClippedLines.append(orderPixelTable.loc[mask])
            mask = (orderPixelTable['sigma_clipped'] == True)
            orderPixelTable = orderPixelTable.loc[~mask]

            # SIGMA-CLIP THE DATA
            self.log.info("""sigma_clip""" % locals())

            if clipOnMphSets:

                columnsNoStrings = list(orderPixelTable.columns)
                try:
                    columnsNoStrings.remove("ion")
                except:
                    pass

                # GROUP BY ARC LINES (MPH SETS)
                lineGroups = orderPixelTable[columnsNoStrings].groupby(['wavelength', 'order']).mean()
                lineGroups = lineGroups.reset_index()

                # SIGMA-CLIP THE DATA ON SCATTER
                masked_residuals = sigma_clip(
                    lineGroups["residuals_x"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaX, maxiters=1, cenfunc='median', stdfunc='mad_std')
                lineGroups["sigma_clipped_x"] = masked_residuals.mask
                masked_residuals = sigma_clip(
                    lineGroups["residuals_y"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaY, maxiters=1, cenfunc='median', stdfunc='mad_std')
                lineGroups["sigma_clipped_y"] = masked_residuals.mask
                lineGroups.loc[((lineGroups["sigma_clipped_y"] == True) | (lineGroups["sigma_clipped_x"] == True)), "sigma_clipped"] = True

                if True:
                    # CLIP ALSO ON COMBINED RESIDUALS
                    masked_residuals = sigma_clip(
                        lineGroups["residuals_xy"], sigma_lower=5000, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
                    lineGroups["sigma_clipped_xy"] = masked_residuals.mask
                    lineGroups.loc[((lineGroups["sigma_clipped_y"] == True) | (lineGroups["sigma_clipped_x"] == True) | (lineGroups["sigma_clipped_xy"] == True)), "sigma_clipped"] = True

                # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FLUX
                mask = (lineGroups["sigma_clipped"] == True)
                clippedGroups = lineGroups.loc[mask]
                clippedGroups = clippedGroups[['wavelength', 'order']]
                s = orderPixelTable[['wavelength', 'order']].merge(clippedGroups, indicator=True, how='left')
                s["clipped"] = False
                s.loc[(s["_merge"] == "both"), "clipped"] = True

                orderPixelTable["sigma_clipped"] = s["clipped"].values

                # CLIP THE MOST DEVIATE SINGLE PINHOLES
                if iteration == 1 and True:
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_x"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaX * 3, maxiters=1, cenfunc='median', stdfunc='mad_std')
                    orderPixelTable["sigma_clipped_x"] = masked_residuals.mask
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_y"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaY * 3, maxiters=1, cenfunc='median', stdfunc='mad_std')
                    orderPixelTable["sigma_clipped_y"] = masked_residuals.mask
                    masked_residuals = sigma_clip(
                        orderPixelTable["R"], sigma_lower=3000, sigma_upper=10, maxiters=1, cenfunc='median', stdfunc='mad_std')
                    orderPixelTable["sigma_clipped_R"] = masked_residuals.mask
                    orderPixelTable.loc[((orderPixelTable["sigma_clipped_y"] == True) | (orderPixelTable["sigma_clipped_x"] == True) | (orderPixelTable["sigma_clipped_R"] == True)), "sigma_clipped"] = True

            else:
                masked_residuals = sigma_clip(
                    orderPixelTable["residuals_x"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaX, maxiters=1, cenfunc='median', stdfunc='mad_std')
                orderPixelTable["sigma_clipped_x"] = masked_residuals.mask
                masked_residuals = sigma_clip(
                    orderPixelTable["residuals_y"].abs(), sigma_lower=3000, sigma_upper=clippingSigmaY, maxiters=1, cenfunc='median', stdfunc='mad_std')
                orderPixelTable["sigma_clipped_y"] = masked_residuals.mask
                orderPixelTable.loc[((orderPixelTable["sigma_clipped_y"] == True) | (orderPixelTable["sigma_clipped_x"] == True)), "sigma_clipped"] = True

                if False:
                    # CLIP ALSO ON COMBINED RESIDUALS
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_xy"], sigma_lower=clippingSigma, sigma_upper=clippingSigma, maxiters=1, cenfunc='median', stdfunc='mad_std')
                    orderPixelTable["sigma_clipped_xy"] = masked_residuals.mask
                    orderPixelTable.loc[((orderPixelTable["sigma_clipped_y"] == True) | (orderPixelTable["sigma_clipped_x"] == True) | (orderPixelTable["sigma_clipped_xy"] == True)), "sigma_clipped"] = True

            # COUNT THE CLIPPED LINES
            mask = (orderPixelTable['sigma_clipped'] == True)
            allClippedLines.append(orderPixelTable.loc[mask])
            totalAllClippedLines = pd.concat(allClippedLines, ignore_index=True)

            # RETURN BREAKDOWN OF COLUMN VALUE COUNT
            valCounts = orderPixelTable[
                'sigma_clipped'].value_counts(normalize=False)
            if True in valCounts:
                clippedCount = valCounts[True]
            else:
                clippedCount = 0

            if iteration > 1:
                # Cursor up one line and clear line
                sys.stdout.flush()
                sys.stdout.write("\x1b[1A\x1b[2K")

            try:
                self.log.print(f'\tITERATION {iteration:02d}: {clippedCount} arc lines where clipped in this iteration of fitting a global dispersion map')
            except:
                pass

            mask = (orderPixelTable['sigma_clipped'] == True)
            orderPixelTable = orderPixelTable.loc[~mask]

        if len(allClippedLines):
            allClippedLines = pd.concat(allClippedLines, ignore_index=True)

        mean_res, std_res, median_res, orderPixelTable = self.calculate_residuals(
            orderPixelTable=orderPixelTable,
            xcoeff=xcoeff,
            ycoeff=ycoeff,
            orderDeg=orderDeg,
            wavelengthDeg=wavelengthDeg,
            slitDeg=slitDeg,
            writeQCs=False,
            pixelRange=True)

        self.log.debug('completed the ``fit_polynomials`` method')
        return xcoeff, ycoeff, orderPixelTable, allClippedLines

    def create_placeholder_images(
            self,
            order=False,
            plot=False,
            reverse=False):
        """*create CCDData objects as placeholders to host the 2D images of the wavelength and spatial solutions from dispersion solution map*

        **Key Arguments:**

        - ``order`` -- specific order to generate the placeholder pixels for. Inner-order pixels set to nan, else set to 0. Default *False* (generate all orders)
        - ``plot`` -- generate plots of placeholder images (for debugging). Default *False*.
        - ``reverse`` -- Inner-order pixels set to 0, else set to nan (reverse of default output).

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
            log=self.log, orderTablePath=self.orderTable, extend=0., order=order)

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

        if self.detectorParams["dispersion-axis"] == "x":
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
            axisACoord_edgeup[axisACoord_edgeup > axisALen] = axisALen
            axisACoord_edgeup[axisACoord_edgeup < 0] = 0
            axisACoord_edgelow[axisACoord_edgelow > axisALen] = axisALen
            axisACoord_edgelow[axisACoord_edgelow < 0] = 0
            axisACoord_edgelow, axisACoord_edgeup, axisBcoord = zip(*[(l, u, b) for l, u, b in zip(axisACoord_edgelow, axisACoord_edgeup, axisBcoord) if l >= 0 and l <= axisALen and u >= 0 and u <= axisALen and b >= 0 and b < axisBLen])
            if reverse:
                for b, u, l in zip(axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)):
                    if self.axisA == "x":
                        wlMap.data[b, l: u] = 0
                        orderMap.data[b, l: u] = o
                    else:
                        wlMap.data[l: u, b] = 0
                        orderMap.data[l: u, b] = o
            else:
                for b, u, l in zip(axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)):
                    if self.axisA == "x":
                        wlMap.data[b, l: u] = np.nan
                        orderMap.data[b, l: u] = np.nan
                    else:
                        wlMap.data[l: u, b] = np.nan
                        orderMap.data[l: u, b] = np.nan

        # SLIT MAP PLACEHOLDER SAME AS WAVELENGTH MAP PLACEHOLDER
        slitMap = wlMap.copy()

        # PLOT CCDDATA OBJECT
        if False:
            import matplotlib.pyplot as plt
            rotatedImg = np.rot90(slitMap.data, 1)
            rotatedImg = slitMap.data
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
            plt.close("all")

        self.log.debug('completed the ``create_placeholder_images`` method')
        return slitMap, wlMap, orderMap

    def map_to_image(
            self,
            dispersionMapPath,
            orders=False):
        """*convert the dispersion map to images in the detector format showing pixel wavelength values and slit positions*

        **Key Arguments:**

        - ``dispersionMapPath`` -- path to the full dispersion map to convert to images
        - ``orders`` -- orders to add to the image map. List. Default *False* (add all orders)

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
        import copy

        self.log.print("\n# CREATING 2D IMAGE MAP FROM DISPERSION SOLUTION\n\n")

        self.dispersionMapPath = dispersionMapPath
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        self.map_to_image_displacement_threshold = self.recipeSettings["map_to_image_displacement_threshold"]
        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        combinedSlitImage = False
        combinedWlImage = False

        if orders:
            theseOrders = orders
        else:
            theseOrders = orderNums

        # DEFINE AN INPUT ARRAY
        inputArray = [(order, minWl, maxWl) for order, minWl,
                      maxWl in zip(orderNums, waveLengthMin, waveLengthMax) if order in theseOrders]

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
        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outDir):
            os.makedirs(outDir)

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
                     overwrite=True, checksum=True)

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

        if np.count_nonzero(wlMap.data) == 0:
            return slitMap, wlMap

        # FIRST GENERATE A WAVELENGTH SURFACE - FINE WL, CHUNKY SLIT-POSITION
        wlRange = maxWl - minWl
        if self.arm.lower == "nir":
            grid_res_wavelength = wlRange / 200
        else:
            grid_res_wavelength = wlRange / 200
        slitLength = self.detectorParams["slit_length"]
        grid_res_slit = slitLength / 20

        halfGrid = (slitLength / 2) * 1.2
        slitArray = np.arange(-halfGrid, halfGrid +
                              grid_res_slit, grid_res_slit)
        wlArray = np.arange(minWl, maxWl, grid_res_wavelength)

        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate(
            [np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
        # NOW THE BIG WAVELEGTH ARRAY
        bigWlArray = np.tile(wlArray, np.shape(slitArray)[0])

        iteration = 0
        remainingPixels = 1
        iterationLimit = 20
        remainingCount = 1

        # GET A COMPLETE LIST OF THE PIXEL WE NEED
        nan_indexes = np.argwhere(np.isnan(slitMap.data))
        ally, allx = nan_indexes[:, 0], nan_indexes[:, 1]

        while remainingPixels and remainingCount and iteration < iterationLimit:
            iteration += 1

            # GENERATE THE ORDER PIXEL TABLE FROM WL AND SLIT-POSITION GRID .. IF WITHIN THRESHOLD OF CENTRE OF DETECTOR PIXEL THEN INJECT INTO MAPS
            orderPixelTable, remainingCount = self.convert_and_fit(
                order=order, bigWlArray=bigWlArray, bigSlitArray=bigSlitArray, slitMap=slitMap, wlMap=wlMap, iteration=iteration, plots=False)

            if remainingCount < 3:
                break

            orderPixelTable = orderPixelTable.drop_duplicates(subset=['pixel_x', 'pixel_y', 'order'])
            train_wlx = orderPixelTable["fit_x"].values
            train_wly = orderPixelTable["fit_y"].values
            train_wl = orderPixelTable["wavelength"].values
            train_sp = orderPixelTable["slit_position"].values

            targetX = orderPixelTable['pixel_x'].values
            targetY = orderPixelTable['pixel_y'].values

            if iteration == 1:
                # ADD MISSING PIXELS
                targetX = np.concatenate([targetX, allx])
                targetY = np.concatenate([targetY, ally])

            # USE CUBIC SPLINE NEAREST NEIGHBOUR TO SEED RESULTS
            bigWlArray = griddata((train_wlx, train_wly), train_wl, (targetX, targetY), method="cubic")
            bigSlitArray = griddata((train_wlx, train_wly), train_sp, (targetX, targetY), method="cubic")

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
        - ``plots`` -- show plot of the slit-map

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

        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        percentageFound = (1 - (np.count_nonzero(np.isnan(wlMap.data)) / np.count_nonzero(wlMap.data))) * 100
        try:
            self.log.print(f"ORDER {order:02d}, iteration {iteration:02d}. {percentageFound:0.2f}% order pixels now fitted.")
        except:
            pass
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
            plt.close("all")

        remainingCount = len(remainingCount.index)

        self.log.debug('completed the ``convert_and_fit`` method')
        return orderPixelTable, remainingCount

    def _create_dispersion_map_qc_plot(
            self,
            xcoeff,
            ycoeff,
            orderDeg,
            wavelengthDeg,
            slitDeg,
            orderPixelTable,
            missingLines,
            allClippedLines,
            dispMap=False,
            dispMapImage=False):
        """*create the QC plot for the dispersion map solution*

        **Key Arguments:**

        - ``xcoeff`` -- the x-coefficients
        - ``ycoeff`` -- the y-coefficients
        - ``orderDeg`` -- degree of the order fitting
        - ``wavelengthDeg`` -- degree of wavelength fitting
        - ``slitDeg`` -- degree of the slit fitting (False for single pinhole)
        - ``orderPixelTable`` -- a panda's data-frame containing wavelength,order,slit_index,slit_position,detector_x,detector_y
        - ``missingLines`` -- lines not detected on the image
        - `allClippedLines` -- lines clipped during dispersion solution fitting
        - `dispMap` -- path to dispersion map. Default *False*
        - `dispMapImage` -- the 2D dispersion map image

        **Return:**

        - ``res_plots`` -- path the the output QC plot
        ```
        """
        self.log.debug('starting the ``create_dispersion_map_qc_plot`` method')
        # print('Creating the QC plot:')
        # print(orderPixelTable)

        import numpy as np
        from astropy.visualization import hist
        import matplotlib.pyplot as plt
        import pandas as pd
        from soxspipe.commonutils.toolkit import qc_settings_plot_tables

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        rotateImage = dp["rotate-qc-plot"]
        flipImage = dp["flip-qc-plot"]

        gridLinePixelTable = False

        # CREATE THE GRID LINES FOR THE FULL SPATIAL SOLUTION PLOTS (MPH)
        if not isinstance(dispMapImage, bool):
            from soxspipe.commonutils.toolkit import create_dispersion_solution_grid_lines_for_plot
            gridLinePixelTable, interOrderMask = create_dispersion_solution_grid_lines_for_plot(
                log=self.log,
                dispMap=dispMap,
                dispMapImage=dispMapImage,
                associatedFrame=self.pinholeFrame,
                kw=kw,
                skylines=False,
                slitPositions=self.uniqueSlitPos
            )

        # DROP MISSING VALUES
        orderPixelTable.dropna(axis='index', how='any', subset=[
            'residuals_x'], inplace=True)
        orderPixelTable["residuals_xy"] = np.sqrt(np.square(
            orderPixelTable["residuals_x"]) + np.square(orderPixelTable["residuals_y"]))
        mean_res = np.mean(orderPixelTable["residuals_xy"])
        std_res = np.std(orderPixelTable["residuals_xy"])
        median_res = np.median(orderPixelTable["residuals_xy"])

        # mean_x_res = abs(orderPixelTable["residuals_x"]).mean()
        # mean_y_res = abs(orderPixelTable["residuals_y"]).mean()

        mean_x_res = abs(orderPixelTable["residuals_x"]).mean()
        mean_y_res = abs(orderPixelTable["residuals_y"]).mean()
        median_x_res = np.median(abs(orderPixelTable["residuals_x"]))
        median_y_res = np.median(abs(orderPixelTable["residuals_y"]))

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = np.ma.array(self.pinholeFrameMasked.data, mask=self.pinholeFrameMasked.mask, fill_value=0).filled()
        if rotateImage:
            rotatedImg = np.rot90(rotatedImg, rotateImage / 90)
        if flipImage:
            rotatedImg = np.flipud(rotatedImg)
            if not rotateImage:
                aLen = rotatedImg.shape[0]
                orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]
                missingLines[f"detector_{self.axisA}_shifted"] = aLen - missingLines[f"detector_{self.axisA}_shifted"]
                allClippedLines[f"detector_{self.axisA}_shifted"] = aLen - allClippedLines[f"detector_{self.axisA}_shifted"]
                orderPixelTable[f"detector_{self.axisA}_shifted"] = aLen - orderPixelTable[f"detector_{self.axisA}_shifted"]
                allClippedLines[f"detector_{self.axisA}"] = aLen - allClippedLines[f"detector_{self.axisA}"]
                missingLines[f"detector_{self.axisA}"] = aLen - missingLines[f"detector_{self.axisA}"]
                orderPixelTable[f"fit_{self.axisA}"] = aLen - orderPixelTable[f"fit_{self.axisA}"]
                allClippedLines[f"observed_{self.axisA}"] = aLen - allClippedLines[f"observed_{self.axisA}"]
                if not isinstance(gridLinePixelTable, bool):
                    gridLinePixelTable[f"fit_{self.axisA}"] = aLen - gridLinePixelTable[f"fit_{self.axisA}"]
                # orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]
                # orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]

        # a = plt.figure(figsize=(40, 15))

        if rotatedImg.shape[0] / rotatedImg.shape[1] > 0.8:
            fig = plt.figure(figsize=(6, 20), constrained_layout=True)
            # CREATE THE GRID OF AXES
            gs = fig.add_gridspec(10, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:6, 0:2])
            bottomright = fig.add_subplot(gs[4:6, 2:])
            resAx = fig.add_subplot(gs[6:8, :])
            settingsAx = fig.add_subplot(gs[8:, 2:])
            qcAx = fig.add_subplot(gs[8:, 0:2])
        else:
            fig = plt.figure(figsize=(6, 20), constrained_layout=True)
            # CREATE THE GRID OF AXES
            gs = fig.add_gridspec(10, 4)
            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:6, 0:2])
            bottomright = fig.add_subplot(gs[4:6, 2:])
            resAx = fig.add_subplot(gs[6:8, :])
            settingsAx = fig.add_subplot(gs[8:, 2:])
            qcAx = fig.add_subplot(gs[8:, 0:2])

        std = self.std
        mean = self.mean

        vmax = mean + 25 * std
        vmin = mean
        toprow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=0.5)
        toprow.set_title(
            "observed arc-line positions (post-clipping)", fontsize=10)

        alphaBoost = 1.0
        if not self.firstGuessMap:
            alphaBoost = 1.7

        if False:
            # PLOT WHERE LINES GET SHIFTED TO
            toprow.scatter(missingLines[f"detector_{self.axisB}_shifted"], missingLines[f"detector_{self.axisA}_shifted"], marker='o', c='black', s=20, alpha=0.1 * alphaBoost, linewidths=0.5, label="missed lines")
            toprow.scatter(allClippedLines[f"detector_{self.axisB}_shifted"], allClippedLines[f"detector_{self.axisA}_shifted"], marker='o', c='black', s=20, alpha=0.1 * alphaBoost, linewidths=0.5, label="missed lines")
            toprow.scatter(orderPixelTable[f"detector_{self.axisB}_shifted"], orderPixelTable[f"detector_{self.axisA}_shifted"], marker='o', c='black', s=20, alpha=0.1 * alphaBoost, linewidths=0.5, label="missed lines")

        else:

            if isinstance(missingLines, pd.core.frame.DataFrame):
                toprow.scatter(missingLines[f"detector_{self.axisB}"], missingLines[f"detector_{self.axisA}"], marker='o', c='black', s=20, alpha=0.1 * alphaBoost, linewidths=0.5, label="undetected line location")
            if len(allClippedLines.index):
                mask = (allClippedLines['dropped'] == True)
                if self.firstGuessMap:
                    toprow.scatter(allClippedLines.loc[mask][f"observed_{self.axisB}"], allClippedLines.loc[mask][f"detector_{self.axisA}"], marker='o', c='blue', s=5, alpha=0.3 * alphaBoost, linewidths=0.5, label="dropped multi-pinhole set")
                else:
                    toprow.scatter(allClippedLines.loc[mask][f"observed_{self.axisB}"], allClippedLines.loc[mask][f"detector_{self.axisA}"], marker='o', c='blue', s=5, alpha=0.3 * alphaBoost, linewidths=0.5, label="dropped pinhole")
                toprow.scatter(allClippedLines.loc[~mask][f"observed_{self.axisB}"], allClippedLines.loc[~mask][f"observed_{self.axisA}"], marker='o', c='green', s=5, alpha=0.3 * alphaBoost, linewidths=0.5 * alphaBoost, )
                toprow.scatter(allClippedLines.loc[~mask][f"observed_{self.axisB}"], allClippedLines.loc[~mask][f"observed_{self.axisA}"], marker='x', c='red', s=5, alpha=0.3 * alphaBoost, linewidths=0.5, label="clipped during dispersion solution fitting")
            if len(orderPixelTable.index):
                toprow.scatter(orderPixelTable[f"observed_{self.axisB}"], orderPixelTable[f"observed_{self.axisA}"], marker='o', c='green', s=5, alpha=0.3 * alphaBoost, linewidths=0.5 * alphaBoost, label="detected line location")

        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        toprow.tick_params(axis='both', which='major', labelsize=9)
        toprow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.1), fontsize=4)

        toprow.set_xlim([0, rotatedImg.shape[1]])
        if self.axisA == "x":
            toprow.invert_yaxis()
        toprow.set_ylim([0, rotatedImg.shape[0]])

        midrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap='gray', alpha=0.5)
        midrow.set_title(
            "global dispersion solution", fontsize=10)

        # ADD FULL DISPERSION SOLUTION GRID-LINES TO PLOT
        if not isinstance(gridLinePixelTable, bool):
            for l in range(int(gridLinePixelTable['line'].max())):
                mask = (gridLinePixelTable['line'] == l)
                if l == 1:
                    midrow.plot(gridLinePixelTable.loc[mask][f"fit_{self.axisB}"], gridLinePixelTable.loc[mask][f"fit_{self.axisA}"], "w-", linewidth=0.2, alpha=0.5 * alphaBoost, color="blue", label="dispersion solution")
                else:
                    midrow.plot(gridLinePixelTable.loc[mask][f"fit_{self.axisB}"], gridLinePixelTable.loc[mask][f"fit_{self.axisA}"], "w-", linewidth=0.2, alpha=0.5 * alphaBoost, color="blue")
        else:
            midrow.scatter(orderPixelTable[f"fit_{self.axisB}"],
                           orderPixelTable[f"fit_{self.axisA}"], marker='o', c='blue', s=orderPixelTable[f"residuals_xy"] * 30, alpha=0.1 * alphaBoost, label="fitted line (size proportional to line-fit residual)")

        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis='both', which='major', labelsize=9)
        midrow.set_xlim([0, rotatedImg.shape[1]])
        if self.axisA == "x":
            midrow.invert_yaxis()
        midrow.set_ylim([0, rotatedImg.shape[0]])

        midrow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.1), fontsize=4)

        # PLOT THE RESIDUALS
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(orderPixelTable[f"residuals_{self.axisA}"], orderPixelTable[
            f"residuals_{self.axisB}"], alpha=0.1)
        bottomleft.set_xlabel(f'{self.axisA} residual (mean: {mean_x_res:2.2f} pix)')
        bottomleft.set_ylabel(f'{self.axisB} residual (mean: {mean_y_res:2.2f} pix)')
        bottomleft.tick_params(axis='both', which='major', labelsize=9)
        hist(orderPixelTable["residuals_xy"], bins='scott', ax=bottomright, histtype='stepfilled',
             alpha=0.7, density=True)
        bottomright.set_xlabel('xy residual')
        bottomright.tick_params(axis='both', which='major', labelsize=9)

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        if self.firstGuessMap:
            fig.suptitle(f"residuals of global dispersion solution fitting - {arm} multi-pinhole\n{subtitle}", fontsize=10, y=0.99)
        else:
            fig.suptitle(f"residuals of global dispersion solution fitting - {arm} single pinhole\n{subtitle}", fontsize=10, y=0.99)
        orderPixelTable_groups = orderPixelTable.groupby(['order'])
        for name, group in orderPixelTable_groups:
            resAx.scatter(group["wavelength"], group["R"], alpha=0.1)
            # CALCULATE THE MEAN AND STD DEV OF THE GROUP AND ADD TO THE PLOT
            mean_resol = group["R"].mean()
            std_resol = group["R"].std()

            mean_wavelength = group["wavelength"].mean()
            # ADD THIS POINT TO THE PLOT USING STD_RES AS ERROR BAR
            # make marker big

            # ADD TO THE POINT THE ERROR BAR CONTAINED IN STD_RED
            resAx.errorbar(mean_wavelength, mean_resol, yerr=std_resol, fmt='o', color='black', alpha=1.0)

        # resAx.scatter(orderPixelTable["wavelength"], orderPixelTable["R"], alpha=0.1)
        resAx.set_xlabel("Wavelength (nm)", fontsize=10)
        resAx.set_ylabel("Resolution", fontsize=10)

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
            polyOrders = [orderDeg, wavelengthDeg, slitDeg]
            if isinstance(orderDeg, list):
                merged_list = []
                for sublist in polyOrders:
                    merged_list.extend(sublist)
                polyOrders = merged_list
            polyOrders[:] = [str(l) for l in polyOrders]
            polyOrders = "".join(polyOrders)
            res_plots = self.sofName + f"_RESIDUALS_{polyOrders}.pdf"

        if self.firstGuessMap:
            filePath = f"{self.qcDir}/{res_plots}"
        else:
            filePath = f"{self.qcDir}/{res_plots}"
        self.products = pd.concat([self.products, pd.Series({
            "soxspipe_recipe": self.recipeName,
            "product_label": "DISP_MAP_RES",
            "file_name": res_plots,
            "file_type": "PDF",
            "obs_date_utc": self.dateObs,
            "reduction_date_utc": utcnow,
            "product_desc": f"{self.arm} dispersion solution QC plots",
            "file_path": filePath,
            "label": "QC"
        }).to_frame().T], ignore_index=True)

        qc_settings_plot_tables(log=self.log, qc=self.qc, qcAx=qcAx, settings={**self.recipeSettings, **{"exptime": self.exptime}}, settingsAx=settingsAx)

        plt.tight_layout()
        # plt.show()
        if not self.settings["tune-pipeline"]:
            plt.savefig(filePath, dpi=720)

        plt.clf()
        plt.close(fig)

        if self.settings["tune-pipeline"]:
            import codecs
            filePath = f"residuals.txt"
            exists = os.path.exists(filePath)
            if not exists:
                with codecs.open(filePath, encoding='utf-8', mode='w') as writeFile:
                    writeFile.write(f"polyOrders,mean_x_res,mean_y_res,mean_res,std_res,median_res,median_x_res,median_y_res,CLINE \n")
            with codecs.open(filePath, encoding='utf-8', mode='a') as writeFile:
                writeFile.write(f"{polyOrders},{mean_x_res:0.4f},{mean_y_res:0.4f},{mean_res:2.4f},{std_res:2.4f},{median_res:2.4f},{median_x_res:2.4f},{median_y_res:2.4f},{self.CLINE}\n")

        self.log.debug('completed the ``create_dispersion_map_qc_plot`` method')
        return res_plots

    def _clip_on_measured_line_metrics(
            self,
            orderPixelTable):
        """*clip lines & sets of lines based on measured line metrics (from daostarfinder etc)*

        **Key Arguments:**

        - `orderPixelTable` -- panda's data-frame containing measure line metrics

        **Return:**

        - `orderPixelTable` -- the data-frame with clipped lines indicated in the "sigma_clipped" column
        """
        self.log.debug('starting the ``_clip_on_measured_line_metrics`` method')

        import matplotlib.pyplot as plt
        from astropy.stats import sigma_clip, sigma_clipped_stats
        from astropy.visualization import hist

        # LAYOUT THE FIGURE
        fig = plt.figure(figsize=(6, 12), constrained_layout=True)
        gs = fig.add_gridspec(8, 4)
        toprow = fig.add_subplot(gs[0:2, :])
        midrow = fig.add_subplot(gs[2:4, :])
        midrow2 = fig.add_subplot(gs[4:6, :])
        if self.firstGuessMap:
            bottomleft = fig.add_subplot(gs[6:8, 0:2])
            bottomright = fig.add_subplot(gs[6:8, 2:])

        toprow.set_title(
            "Pre-fitting QC Clipping on Pinhole-Sets", fontsize=10)

        columnsNoStrings = list(orderPixelTable.columns)
        try:
            columnsNoStrings.remove("ion")
        except:
            pass

        if self.firstGuessMap:

            # GROUP BY ARC LINES (MPH SETS)

            lineGroups = orderPixelTable.loc[(orderPixelTable["dropped"] == False)][columnsNoStrings].groupby(['wavelength', 'order']).std()
            lineGroups = lineGroups.reset_index()

            # SIGMA-CLIP THE DATA ON SCATTER
            # masked_residuals = sigma_clip(
            #     lineGroups["xy_diff"], sigma_lower=5000, sigma_upper=7, maxiters=5, cenfunc='median', stdfunc='mad_std')
            # lineGroups["sigma_clipped_scatter"] = masked_residuals.mask

            # SIGMA-CLIP THE DATA ON SCATTER
            lineGroups["sigma_clipped_scatter"] = False
            masked_residuals = sigma_clip(
                lineGroups["x_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc='median', stdfunc='mad_std')
            lineGroups["sigma_clipped_x"] = masked_residuals.mask
            masked_residuals = sigma_clip(
                lineGroups["y_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc='median', stdfunc='mad_std')
            lineGroups["sigma_clipped_y"] = masked_residuals.mask
            masked_residuals = sigma_clip(
                lineGroups["xy_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc='median', stdfunc='mad_std')
            lineGroups["sigma_clipped_xy"] = masked_residuals.mask
            lineGroups.loc[((lineGroups["sigma_clipped_y"] == True) | (lineGroups["sigma_clipped_x"] == True) | (lineGroups["sigma_clipped_xy"] == True)), "sigma_clipped_scatter"] = True

            bottomleft.scatter(
                x=lineGroups["x_diff"],  # numpy array of x-points
                y=lineGroups["y_diff"],  # numpy array of y-points
                s=1,    # 1 number or array of areas for each datapoint (i.e. point size)
                c="black",    # color or sequence of color, optional, default
                marker='x',
                alpha=0.6,
                label="arc line shift from predicted position")

            bottomleft.scatter(
                x=lineGroups.loc[(lineGroups["sigma_clipped_scatter"] == True)]["x_diff"],  # numpy array of x-points
                y=lineGroups.loc[(lineGroups["sigma_clipped_scatter"] == True)]["y_diff"],  # numpy array of y-points
                s=5,    # 1 number or array of areas for each datapoint (i.e. point size)
                c="red",    # color or sequence of color, optional, default
                marker='x',
                alpha=0.6,
                label="clipped arc lines")

            bottomleft.set_ylabel(f"y-shift rms (px)", fontsize=12)
            bottomleft.set_xlabel(f"x-shift rms (px)", fontsize=12)
            bottomleft.tick_params(axis='both', which='major', labelsize=9)
            bottomleft.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

            hist(lineGroups[(lineGroups["sigma_clipped_scatter"] == False)]["xy_diff"], bins='scott', ax=bottomright, histtype='stepfilled',
                 alpha=0.7, density=True)
            bottomright.set_xlabel('xy residual')
            bottomright.tick_params(axis='both', which='major', labelsize=9)

            # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FWHM
            mask = (lineGroups["sigma_clipped_scatter"] == True)
            dropGroups = lineGroups.loc[mask]
            setsToDrop = dropGroups[['wavelength', 'order']]
            s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
            s["dropped"] = False
            s.loc[(s["_merge"] == "both"), "dropped"] = True
            orderPixelTable["droppedOnScatter"] = s["dropped"].values
            orderPixelTable.loc[(orderPixelTable["droppedOnScatter"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FWHM
        lineGroups = orderPixelTable.loc[(orderPixelTable["dropped"] == False)][columnsNoStrings].groupby(['wavelength', 'order']).mean()
        lineGroups = lineGroups.reset_index()
        masked_residuals = sigma_clip(
            lineGroups["fwhm_px"], sigma_lower=2.5, sigma_upper=5, maxiters=3, cenfunc='median', stdfunc='mad_std')
        lineGroups["sigma_clipped_fwhm"] = masked_residuals.mask
        lineGroups["sigma_clipped"] = masked_residuals.mask

        toprow.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["fwhm_px"],  # numpy array of y-points
            s=1,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="black",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="measured pinhole lines")

        toprow.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_fwhm"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_fwhm"] == True)]["fwhm_px"],  # numpy array of y-points
            s=5,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="red",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="clipped pinhole lines")

        toprow.set_ylabel(f"fwhm (px)", fontsize=12)
        toprow.set_xlabel(f"wavelength (nm)", fontsize=12)
        toprow.tick_params(axis='both', which='major', labelsize=9)
        toprow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FLUX
        mask = (lineGroups["sigma_clipped_fwhm"] == True)
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[['wavelength', 'order']]
        s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnFWHM"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnFWHM"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FLUX
        lineGroups = lineGroups.loc[~mask]
        masked_residuals = sigma_clip(
            lineGroups["flux"], sigma_lower=5, sigma_upper=5, maxiters=5, cenfunc='mean', stdfunc='std')
        lineGroups["sigma_clipped_flux"] = masked_residuals.mask

        midrow.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["flux"],  # numpy array of y-points
            s=1,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="black",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="pinhole flux")

        midrow.set_ylabel(f"pinhole flux", fontsize=12)
        midrow.set_xlabel(f"wavelength (nm)", fontsize=12)
        midrow.tick_params(axis='both', which='major', labelsize=9)
        midrow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

        midrow.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_flux"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_flux"] == True)]["flux"],  # numpy array of y-points
            s=5,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="red",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="clipped pinhole lines")

        midrow.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON PEAK
        mask = (lineGroups["sigma_clipped_flux"] == True)
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[['wavelength', 'order']]
        s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnFlux"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnFlux"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FLUX
        lineGroups = lineGroups.loc[~mask]
        lineGroups["peak"] = lineGroups["peak"] / lineGroups["flux"]
        masked_residuals = sigma_clip(
            lineGroups["peak"], sigma_lower=5000, sigma_upper=7, maxiters=5, cenfunc='mean', stdfunc='std')
        lineGroups["sigma_clipped_peak"] = masked_residuals.mask

        midrow2.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["peak"],  # numpy array of y-points
            s=1,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="black",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="pinhole peak flux / mean flux")

        midrow2.set_ylabel(f"pinhole peak flux", fontsize=12)
        midrow2.set_xlabel(f"wavelength (nm)", fontsize=12)
        midrow2.tick_params(axis='both', which='major', labelsize=9)
        midrow2.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

        midrow2.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_peak"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_peak"] == True)]["peak"],  # numpy array of y-points
            s=5,    # 1 number or array of areas for each datapoint (i.e. point size)
            c="red",    # color or sequence of color, optional, default
            marker='x',
            alpha=0.6,
            label="clipped pinhole lines")

        midrow2.legend(loc='upper right', bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA
        mask = (lineGroups["sigma_clipped_peak"] == True)
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[['wavelength', 'order']]
        s = orderPixelTable[['wavelength', 'order']].merge(setsToDrop, indicator=True, how='left')
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnPeak"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnPeak"] == True), "dropped"] = True

        # plt.show()
        plt.close("all")

        self.log.debug('completed the ``_clip_on_measured_line_metrics`` method')
        return orderPixelTable

    def create_new_static_line_list(
            self,
            dispersionMapPath):
        """*using a first pass dispersion solution, use a line atlas to generate a more accurate and more complete static line list*

        **Key Arguments:**

        - `dispersionMapPath` -- path to the first pass dispersion solution

        **Return:**

        - `newPredictedLineList` -- a new predicted line list (to replace the static calibration line-list)
        """
        self.log.debug('starting the ``create_new_static_line_list`` method')

        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        from soxspipe.commonutils.toolkit import read_spectral_format
        import pandas as pd
        from astropy.table import Table

        dp = self.detectorParams

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm)

        # FIND THE LINE ATLAS
        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        lineAtlas = calibrationRootPath + "/" + dp["line-atlas"]
        # LINE LIST TO PANDAS DATAFRAME
        lineAtlas = Table.read(lineAtlas, format='fits')
        lineAtlas = lineAtlas.to_pandas()

        # CLEAN UP LINE LIST DATA
        lineAtlas['ion'] = lineAtlas['ion'].str.decode('ascii')
        lineAtlas['source'] = lineAtlas['source'].str.decode('ascii')

        lineAtlas['wave'] = lineAtlas['wave'].round(5)

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK

        # GET UNIQUE VALUES IN COLUMN
        uniqueions = lineAtlas['ion'].unique()

        # # mask = (lineAtlas['ion'].isin(["XeI", 'UNK']))
        # mask = (lineAtlas['ion'].isin(["HgI", 'UNK']))
        # # mask = (lineAtlas['ion'].isin(["NeI", "NeII", 'UNK']))
        # # mask = (lineAtlas['ion'].isin(["Ar", "ArI", "ArII", "ArIII", 'UNK']))
        # lineAtlas = lineAtlas.loc[mask]

        dfCollection = []
        for o, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax):

            wrange = wmax - wmin
            wmin -= wrange / 5
            wmax += wrange / 5

            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            mask = (lineAtlas['wave'].between(wmin, wmax))
            filteredDf = lineAtlas.loc[mask]
            filteredDf = filteredDf.drop_duplicates(subset=["wave", "ion"])
            wave = filteredDf["wave"]
            ion = filteredDf["ion"]

            if not self.firstGuessMap:
                slit_positions = [0.0]
                slit_indexes = [4]
            else:
                slit_positions = self.uniqueSlitPos
                slit_indexes = list(range(0, len(self.uniqueSlitPos), 1))

            # xpd-update-filter-dataframe-column-values

            for si, sp in zip(slit_indexes, slit_positions):
                myDict = {
                    "order": [o] * len(wave),
                    "wavelength": wave,
                    "slit_position": [sp] * len(wave),
                    "slit_index": [si] * len(wave),
                    "ion": ion}

                orderPixelTable = pd.DataFrame(myDict)
                floats = ["wavelength", 'slit_position']
                ints = ['order', 'slit_index']
                for f in floats:
                    orderPixelTable[f] = orderPixelTable[f].astype(float)
                for i in ints:
                    orderPixelTable[i] = orderPixelTable[i].astype(int)

                orderPixelTable = dispersion_map_to_pixel_arrays(
                    log=self.log,
                    dispersionMapPath=dispersionMapPath,
                    orderPixelTable=orderPixelTable,
                    removeOffDetectorLocation=False
                )

                orderPixelTable.rename(columns={"fit_x": "detector_x", "fit_y": "detector_y"}, inplace=True)
                orderPixelTable = orderPixelTable[['ion', 'wavelength', 'order', 'slit_index', 'slit_position', 'detector_x', 'detector_y']]

                dfCollection.append(orderPixelTable)

        newPredictedLineList = pd.concat(dfCollection, ignore_index=True)

        from astropy.table import Table
        t = Table.from_pandas(newPredictedLineList)
        # t.write("Xe.fits", overwrite=True)
        t.write("Hg.fits", overwrite=True)
        # t.write("Ne.fits", overwrite=True)
        # t.write("Ar.fits", overwrite=True)

        # sys.exit(0)

        self.log.debug('completed the ``create_new_static_line_list`` method')
        return newPredictedLineList

    # use the tab-trigger below for new method
    # xt-class-method
