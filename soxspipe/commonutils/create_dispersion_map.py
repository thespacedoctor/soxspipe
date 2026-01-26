#!/usr/bin/env python
# encoding: utf-8
"""
*detect arc-lines on a pinhole frame to generate a dispersion solution*

Author
: Marco Landoni & David Young

Date Created
: September  1, 2020

Module Structure
: This module contains the main `create_dispersion_map` class and helper functions:

  **Main Class:**
  - `create_dispersion_map` - Generates dispersion solutions from pinhole arc frames

  **Helper Functions:**
  - `measure_line_position()` - Detect and measure arc line positions on image stamps
  - `straighten_mph_sets()` - Straighten multi-pinhole sets by projecting onto fitted line
  - `find_largest_cluster_center()` - Find center of largest cluster using DBSCAN
  - `_plot_slit_index_comparisons()` - Debug visualization of slit position residuals

  **Workflow:**
  1. Load predicted line positions from static calibrations
  2. Detect arc lines on pinhole frame iteratively
  3. Calculate shifts between predicted and observed positions
  4. Fit polynomial dispersion solution (wavelength, order, slit)
  5. Generate QC plots and write output files
"""
################# GLOBAL IMPORTS ####################
from soxspipe.commonutils.toolkit import unpack_order_table, read_spectral_format, twoD_disp_map_image_to_dataframe
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
from datetime import datetime, timezone

os.environ["TERM"] = "vt100"


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
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
    - ``arcFrame`` -- the calibrated arc frame used to determine spectral resolution. Default *None*
    - ``debug`` -- debug mode. Default *False*

    **Usage:**

    ```python
    from soxspipe.commonutils import create_dispersion_map
    mapPath, mapImagePath, res_plots, qcTable, productsTable = create_dispersion_map(
        log=log,
        settings=settings,
        pinholeFrame=frame,
        firstGuessMap=False,
        qcTable=self.qc,
        productsTable=self.products,
        sofName=sofName,
        create2DMap=True,
        recipeSettings=recipeSettings,
        startNightDate=startNightDate,
        arcFrame=arcFrame,
        debug=debug
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
        lineDetectionTable=False,
        startNightDate="",
        arcFrame=None,
        debug=False,
    ):
        self.log = log
        log.debug("instantiating a new 'create_dispersion_map' object")

        import warnings
        import copy

        # STORE INITIALIZATION PARAMETERS
        self._store_init_params(
            settings,
            recipeSettings,
            pinholeFrame,
            firstGuessMap,
            orderTable,
            qcTable,
            productsTable,
            sofName,
            create2DMap,
            lineDetectionTable,
            startNightDate,
            arcFrame,
            debug,
            copy,
        )

        # SETUP KEYWORD LOOKUP AND EXTRACT FRAME METADATA
        self._setup_keywords_and_metadata()

        # VALIDATE NIR ARM LAMP REQUIREMENTS
        self._validate_nir_lamp_requirements()

        # DETERMINE RECIPE NAME BASED ON firstGuessMap FLAG
        self._set_recipe_name()

        # INITIALIZE DETECTOR PARAMETERS
        self._setup_detector_params()

        # CONFIGURE WARNINGS AND LOGGING LEVELS
        self._configure_warnings_and_logging(warnings)

        # SET IMAGE ORIENTATION AXES
        self._set_image_orientation()

        # CREATE OUTPUT DIRECTORIES FOR QC AND PRODUCTS
        self._setup_output_directories()

        return None

    def _store_init_params(
        self,
        settings,
        recipeSettings,
        pinholeFrame,
        firstGuessMap,
        orderTable,
        qcTable,
        productsTable,
        sofName,
        create2DMap,
        lineDetectionTable,
        startNightDate,
        arcFrame,
        debug,
        copy,
    ):
        """*Store all initialization parameters as instance variables*"""
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
        self.startNightDate = startNightDate
        self.arcFrame = arcFrame
        self.debug = debug

    def _setup_keywords_and_metadata(self):
        """*Initialize keyword lookup and extract frame header metadata*"""
        from soxspipe.commonutils.toolkit import get_calibration_lamp

        # SETUP KEYWORD LOOKUP OBJECT
        kw = keyword_lookup(log=self.log, settings=self.settings).get
        self.kw = kw

        # EXTRACT HEADER METADATA
        self.arm = self.pinholeFrame.header[kw("SEQ_ARM")]
        self.dateObs = self.pinholeFrame.header[kw("DATE_OBS")]
        self.inst = self.pinholeFrame.header[kw("INSTRUME")]
        self.exptime = self.pinholeFrame.header[kw("EXPTIME")]

        # DETERMINE CALIBRATION LAMP TYPE
        self.lamp = get_calibration_lamp(log=self.log, frame=self.pinholeFrame, kw=kw)

    def _validate_nir_lamp_requirements(self):
        """*Validate that NIR arm has required Argon and Mercury lamps*"""
        if self.arm.upper() == "NIR":
            has_ar = "ar" in self.lamp.lower()
            has_hg = "hg" in self.lamp.lower()
            if not (has_ar and has_hg):
                raise Exception("NIR arm requires both Argon (Ar) and Mercury (Hg) lamps")

    def _set_recipe_name(self):
        """*Determine recipe name based on whether this is first guess or spatial solution*"""
        if self.firstGuessMap:
            self.recipeName = "soxs-spat-solution"
        else:
            self.recipeName = "soxs-disp-solution"

    def _setup_detector_params(self):
        """*Initialize detector parameters lookup for the current arm*"""
        self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

    def _configure_warnings_and_logging(self, warnings):
        """*Configure warning filters and reset logging levels*"""
        from photutils.utils import NoDetectionsWarning
        import logging

        # SUPPRESS PHOTUTILS NO DETECTIONS WARNINGS
        warnings.simplefilter("ignore", NoDetectionsWarning)

        # FIX ASTROPY LOGGING LEVEL RESET
        logging.getLogger().setLevel(logging.INFO + 5)

    def _set_image_orientation(self):
        """*Set dispersion and spatial axes based on detector configuration*"""
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

    def _setup_output_directories(self):
        """*Create QC and product output directories*"""
        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log, settings=self.settings, recipeName=self.recipeName, startNightDate=self.startNightDate
        )

    def _initialize_recipe_settings(self):
        """*INITIALIZE POLYNOMIAL DEGREES AND DETECTION WINDOW SETTINGS*"""
        # READ BOOTSTRAP AND FITTING PARAMETERS FROM SETTINGS
        bootstrap_dispersion_solution = self.settings["bootstrap_dispersion_solution"]
        tightFit = self.settings["bootstrap_dispersion_solution"]
        orderDeg = self.recipeSettings["order-deg"]
        wavelengthDeg = self.recipeSettings["wavelength-deg"]
        self.windowSize = self.recipeSettings["pixel-window-size"]

        # DETERMINE SLIT POLYNOMIAL DEGREE BASED ON FRAME TYPE
        if self.firstGuessMap:
            slitDeg = self.recipeSettings["slit-deg"]
        else:
            # SINGLE PINHOLE: NO SLIT VARIATION
            slitDeg = [0, 0] if isinstance(orderDeg, list) else 0

        return bootstrap_dispersion_solution, tightFit, orderDeg, wavelengthDeg, slitDeg

    def _prepare_pinhole_frame(self):
        """*MASK FRAME AND CALCULATE STATISTICS FOR LINE DETECTION*"""
        import numpy as np
        from astropy.stats import sigma_clipped_stats

        # CREATE MASKED ARRAY FROM FRAME DATA
        pinholeFrameMasked = np.ma.array(self.pinholeFrame.data, mask=self.pinholeFrame.mask)

        # CALCULATE ROBUST STATISTICS FOR BACKGROUND ESTIMATION
        mean, median, std = sigma_clipped_stats(
            pinholeFrameMasked, sigma=5.0, stdfunc="mad_std", cenfunc="median", maxiters=3
        )

        # STORE FOR LATER USE IN LINE DETECTION
        self.meanFrameFlux = mean
        self.stdFrameFlux = std
        self.pinholeFrameMasked = pinholeFrameMasked

        return pinholeFrameMasked

    def _generate_quicklook_and_debug_plots(self, orderPixelTable, pinholeFrameMasked):
        """*GENERATE QUICKLOOK IMAGE AND OPTIONAL DEBUG PLOTS*"""
        from soxspipe.commonutils.toolkit import quicklook_image

        # GENERATE QUICKLOOK IMAGE WITH SURFACE PLOT
        quicklook_image(
            log=self.log,
            CCDObject=self.pinholeFrame,
            show=self.debug,
            ext=False,
            stdWindow=3,
            surfacePlot=True,
            title="Pinhole Frame",
        )

        # DISPLAY DEBUG PLOT OF PREDICTED LINE WINDOWS IF IN DEBUG MODE
        if self.debug:
            self._plot_predicted_line_windows(orderPixelTable, pinholeFrameMasked)

    def _calculate_position_differences(self, orderPixelTable):
        """*CALCULATE DIFFERENCES BETWEEN PREDICTED AND OBSERVED POSITIONS*"""
        import numpy as np

        # CALCULATE X AND Y DIFFERENCES
        orderPixelTable["x_diff"] = orderPixelTable["detector_x_shifted"] - orderPixelTable["observed_x"]
        orderPixelTable["y_diff"] = orderPixelTable["detector_y_shifted"] - orderPixelTable["observed_y"]

        # CALCULATE EUCLIDEAN DISTANCE
        orderPixelTable["xy_diff"] = np.sqrt(
            np.square(orderPixelTable["x_diff"]) + np.square(orderPixelTable["y_diff"])
        )

        if "mph_mean_x" not in orderPixelTable.columns or "mph_mean_y" not in orderPixelTable.columns:
            # FOR EACH UNIQUE (wavelength, order), ADD mph_mean_x AND mph_mean_y COLUMNS
            # GROUP BY (wavelength, order) AND CALCULATE MEAN detector_x AND detector_y
            means = (
                orderPixelTable.groupby(["wavelength", "order"])[["detector_x", "detector_y"]]
                .mean()
                .rename(columns={"detector_x": "mph_mean_x", "detector_y": "mph_mean_y"})
            )
            # MERGE THE MEANS BACK TO THE ORIGINAL DATAFRAME
            orderPixelTable = orderPixelTable.merge(means, on=["wavelength", "order"], how="left")

        return orderPixelTable

    def _get_detection_parameters(self, iteration, tightFit):
        """*DETERMINE DETECTION PARAMETERS BASED ON ITERATION NUMBER*"""
        # INITIALIZE DEFAULT PARAMETERS
        returnAll = False
        brightest = False
        exclude_border = True

        if tightFit:
            # USE STRICT WINDOW SIZE FOR TIGHT FITTING
            windowHalf = round(self.windowSize / 2)
            sigmaLimit = self.recipeSettings["pinhole-detection-thres-sigma"]
        else:
            # PROGRESSIVELY REFINE SEARCH WINDOW AND DETECTION THRESHOLD
            if iteration == 0:
                # ITERATION 0: WIDE SEARCH TO CATCH ALL POTENTIAL LINES
                windowHalf = min(round(self.windowSize * 3), 25)
                sigmaLimit = (
                    10 if not self.firstGuessMap else max(self.recipeSettings["pinhole-detection-thres-sigma"], 5)
                )
                returnAll = not self.firstGuessMap
            elif iteration == 1:
                # ITERATION 1: INTERMEDIATE REFINEMENT
                windowHalf = min(round(self.windowSize * 2), 10 if not self.firstGuessMap else 8)
                sigmaLimit = (
                    10 if not self.firstGuessMap else max(self.recipeSettings["pinhole-detection-thres-sigma"], 2)
                )
                returnAll = not self.firstGuessMap
            else:
                # ITERATION 2+: FINAL PRECISE SEARCH
                windowHalf = round(self.windowSize / 2)
                sigmaLimit = self.recipeSettings["pinhole-detection-thres-sigma"]

        return windowHalf, sigmaLimit, returnAll, brightest, exclude_border

    def _explode_multiple_detections(self, orderPixelTable):
        """*EXPLODE DATAFRAME WHEN MULTIPLE LINES DETECTED FOR SINGLE PREDICTION*"""
        import pandas as pd

        # FIND COLUMNS CONTAINING LISTS (MULTIPLE DETECTIONS)
        exploded_columns = [col for col in orderPixelTable.columns if isinstance(orderPixelTable[col].iloc[0], list)]

        # EXPLODE LISTS INTO SEPARATE ROWS
        orderPixelTable = orderPixelTable.explode(exploded_columns, ignore_index=True)

        # CONVERT STRING VALUES BACK TO NUMERIC
        for col in exploded_columns:
            orderPixelTable[col] = pd.to_numeric(orderPixelTable[col], errors="coerce")

        return orderPixelTable

    def _write_qc_metrics(self, totalLines, detectedLines, percentageDetectedLines, utcnow):
        """*WRITE LINE DETECTION QC METRICS TO TABLE*"""
        import pandas as pd

        # DETERMINE TAG BASED ON RECIPE TYPE
        tag = "single" if "DISP" in self.recipeName.upper() else "multi"

        # WRITE TOTAL LINES QC
        self.qc = pd.concat(
            [
                self.qc,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "qc_name": "TLINE",
                        "qc_value": totalLines,
                        "qc_comment": f"Total number of line in {tag} line-list",
                        "qc_unit": "lines",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

        # WRITE DETECTED LINES QC
        self.qc = pd.concat(
            [
                self.qc,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "qc_name": "NLINE",
                        "qc_value": detectedLines,
                        "qc_comment": f"Number of lines detected in {tag} pinhole frame",
                        "qc_unit": "lines",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

        # WRITE PERCENTAGE DETECTED QC
        self.qc = pd.concat(
            [
                self.qc,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "qc_name": "PLINE",
                        "qc_value": percentageDetectedLines,
                        "qc_comment": f"Proportion of input line-list lines detected on {tag} pinhole frame",
                        "qc_unit": None,
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

    def _find_cluster_center_with_fallback(self, x_data, y_data, fallback_median):
        """*FIND CLUSTER CENTER USING DBSCAN WITH AUTOMATIC PARAMETER ADJUSTMENT*"""
        # START WITH MODERATE CLUSTERING PARAMETERS
        centre = None
        eps = 0.7

        # PROGRESSIVELY RELAX PARAMETERS IF NO CLUSTER FOUND
        while not centre and eps < 2.0:
            min_samples = 25
            while not centre and min_samples > 5:
                centre_x, centre_y = find_largest_cluster_center(x_data, y_data, eps=eps, min_samples=min_samples)
                if centre_x:
                    centre = (centre_x, centre_y)
                min_samples -= 1
            if not centre:
                eps += 0.1

        # FALLBACK TO MEDIAN IF NO CLUSTER FOUND
        if not centre:
            centre = (fallback_median, None)

        return centre

    def _calculate_order_shift_statistics(self, orderPixelTable, mask):
        """*CALCULATE SHIFT STATISTICS FOR A SINGLE ORDER*"""
        from astropy.stats import sigma_clipped_stats
        import numpy as np

        # CALCULATE XY DISTANCE STATISTICS
        meanxy, medianxy, stdxy = sigma_clipped_stats(
            orderPixelTable.loc[mask]["xy_diff"], sigma=1.0, stdfunc="std", cenfunc="mean", maxiters=20
        )

        # CALCULATE X AND Y SHIFT STATISTICS SEPARATELY
        updatedMask = mask
        meanx, medianx, stdx = sigma_clipped_stats(
            orderPixelTable.loc[updatedMask]["x_diff"], sigma=1.5, stdfunc="std", cenfunc="mean", maxiters=7
        )
        meany, mediany, stdy = sigma_clipped_stats(
            orderPixelTable.loc[updatedMask]["y_diff"], sigma=1.5, stdfunc="std", cenfunc="mean", maxiters=7
        )

        return medianx, mediany, stdx, stdy, medianxy, stdxy

    def _find_and_apply_cluster_shift(self, orderPixelTable, mask):
        """*FIND CLUSTER CENTER IN SHIFT DISTRIBUTION AND APPLY CORRECTION*"""
        # FIND CLUSTER CENTER IN SHIFT SPACE
        centrex = None
        centrey = None
        eps = 0.7
        orderLen = len(orderPixelTable.loc[mask].index)

        # PROGRESSIVELY RELAX CLUSTERING PARAMETERS
        while not centrex and eps < 1.2 and orderLen:
            min_samples = 25
            while not centrex and min_samples > 10:
                centrex, centrey = find_largest_cluster_center(
                    orderPixelTable.loc[mask]["x_diff"],
                    orderPixelTable.loc[mask]["y_diff"],
                    eps=eps,
                    min_samples=min_samples,
                )
                min_samples -= 1
            if not centrex:
                eps += 0.1

        # APPLY SHIFT CORRECTION IF CLUSTER FOUND
        if centrex and centrey:
            orderPixelTable.loc[mask, "detector_x_shifted"] = orderPixelTable.loc[mask]["detector_x_shifted"] - centrex
            orderPixelTable.loc[mask, "detector_y_shifted"] = orderPixelTable.loc[mask]["detector_y_shifted"] - centrey

        return orderPixelTable, centrex, centrey

    def _handle_multipin_hole_big_shift(self, orderPixelTable, order_num, mask, iteration):
        """*DETECT AND CORRECT LARGE SHIFTS BETWEEN SINGLE/MULTI-PINHOLE FRAMES*"""
        from astropy.stats import sigma_clipped_stats

        # ONLY CHECK ON FIRST ITERATION
        if iteration != 0:
            return orderPixelTable

        # CALCULATE MEDIAN SHIFTS FOR TOP AND BOTTOM SLIT POSITIONS
        _, medTop, _ = sigma_clipped_stats(
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 8))][f"{self.axisA}_diff"],
            sigma=1.0,
            stdfunc="mad_std",
            cenfunc="median",
            maxiters=3,
        )
        _, medBottom, _ = sigma_clipped_stats(
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 0))][f"{self.axisA}_diff"],
            sigma=1.0,
            stdfunc="mad_std",
            cenfunc="median",
            maxiters=3,
        )

        # FIND CLUSTER CENTERS FOR TOP AND BOTTOM SLITS
        centreTop = self._find_cluster_center_with_fallback(
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 8))][f"{self.axisA}_diff"],
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 8))][f"{self.axisB}_diff"],
            medTop,
        )
        centreBottom = self._find_cluster_center_with_fallback(
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 0))][f"{self.axisA}_diff"],
            orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 0))][f"{self.axisB}_diff"],
            medBottom,
        )

        centreATop, centreBTop = centreTop
        centreABottom, centreBBottom = centreBottom

        # CHECK IF SHIFT IS SIGNIFICANT (> 0.4 PIXELS)
        if abs(centreATop - centreABottom) > 0.4:
            # USE LARGER SHIFT
            shift = centreATop if abs(centreATop) > abs(centreABottom) else centreABottom
            print(f"{order_num} APPLYING BIG SHIFT {shift}")

            # APPLY CORRECTION
            orderPixelTable.loc[mask, f"detector_{self.axisA}_shifted"] = (
                orderPixelTable.loc[mask][f"detector_{self.axisA}_shifted"] - shift
            )

            # DEBUG PLOTS IF ENABLED
            if self.debug:
                _plot_slit_index_comparisons(orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 8))])
                _plot_slit_index_comparisons(orderPixelTable.loc[(mask & (orderPixelTable["slit_index"] == 0))])
        elif self.debug:
            _plot_slit_index_comparisons(orderPixelTable.loc[mask])

        return orderPixelTable

    def _reduce_polynomial_degrees(self, popt_x, popt_y, wavelengthDeg, orderDeg, slitDeg):
        """*REDUCE POLYNOMIAL DEGREES WHEN FIT FAILS*"""
        # HANDLE DIFFERENT POLYNOMIAL DEGREE FORMATS
        if not isinstance(wavelengthDeg, list):
            # SINGLE DEGREE FOR BOTH AXES
            degList = [wavelengthDeg, orderDeg, slitDeg]
            degList[degList.index(max(degList))] -= 1
            wavelengthDeg, orderDeg, slitDeg = degList
        elif popt_x == "xerror":
            # REDUCE X-AXIS DEGREES
            degList = [wavelengthDeg[0], orderDeg[0], slitDeg[0]]
            degList[degList.index(max(degList))] -= 1
            wavelengthDeg[0], orderDeg[0], slitDeg[0] = degList
        elif popt_y == "yerror":
            # REDUCE Y-AXIS DEGREES
            degList = [wavelengthDeg[1], orderDeg[1], slitDeg[1]]
            degList[degList.index(max(degList))] -= 1
            wavelengthDeg[1], orderDeg[1], slitDeg[1] = degList

        # UPDATE RECIPE SETTINGS
        self.recipeSettings["order-deg"] = orderDeg
        self.recipeSettings["wavelength-deg"] = wavelengthDeg
        if self.firstGuessMap:
            self.recipeSettings["slit-deg"] = slitDeg

        return wavelengthDeg, orderDeg, slitDeg

    def get(self):
        """
        *generate the dispersion map*

        **Return:**

        - ``mapPath`` -- path to the file containing the coefficients of the x,y polynomials of the global dispersion map fit

        **Workflow:**

        1. LOAD PREDICTED LINE POSITIONS FROM CALIBRATION FILES
        2. PREPARE PINHOLE FRAME (MASKING, STATISTICS)
        3. ITERATIVELY DETECT ARC LINES ON FRAME
        4. FIT POLYNOMIAL DISPERSION SOLUTION
        5. WRITE OUTPUTS (LINE LISTS, MAP FILES, QC PLOTS)
        """
        self.log.debug("starting the ``get`` method")

        import pandas as pd
        from astropy.table import Table
        import numpy as np
        from astropy.stats import sigma_clipped_stats
        from astropy.stats import sigma_clip

        # STEP 1: INITIALIZE RECIPE SETTINGS AND POLYNOMIAL DEGREES
        bootstrap_dispersion_solution, tightFit, orderDeg, wavelengthDeg, slitDeg = self._initialize_recipe_settings()

        # STEP 2: LOAD PREDICTED LINE POSITIONS FROM CALIBRATIONS
        orderPixelTable = self.get_predicted_line_list()
        originalOrderPixelTable = orderPixelTable.copy()
        totalLines = len(orderPixelTable.index)
        self.uniqueSlitPos = orderPixelTable["slit_position"].unique()

        # STEP 3: PREPARE PINHOLE FRAME FOR LINE DETECTION
        pinholeFrameMasked = self._prepare_pinhole_frame()

        # STEP 4: GENERATE QUICKLOOK IMAGE AND DEBUG PLOTS
        self._generate_quicklook_and_debug_plots(orderPixelTable, pinholeFrameMasked)

        boost = True
        while boost:
            # SORT BY COLUMN NAME
            orderPixelTable.sort_values(["wavelength"], inplace=True)

            # BOOST WILL BE SET TO TRUE LATER IF FOUND TO BE TRUE IN THE SETTINGS FILE
            boost = False
            bigShift = False

            if not isinstance(self.lineDetectionTable, bool):
                # USE THE PROVIDED LINE DETECTION TABLE (FOR PIPELINE TUNING ONLY)
                orderPixelTable = self.lineDetectionTable.copy()
            else:
                # DETECT THE LINES ON THE PINHOLE FRAME AND
                # ADD OBSERVED LINES TO DATAFRAME
                iteration = 0
                self.log.print(f"\n# FINDING PINHOLE ARC-LINES ON IMAGE\n")
                iraf = False
                while iteration < 3:

                    # GET DETECTION PARAMETERS FOR CURRENT ITERATION
                    self.windowHalf, sigmaLimit, returnAll, brightest, exclude_border = self._get_detection_parameters(
                        iteration, tightFit
                    )

                    # DETECT LINES ON PINHOLE FRAME
                    orderPixelTable = self.detect_pinhole_arc_lines(
                        orderPixelTable=orderPixelTable,
                        iraf=iraf,
                        sigmaLimit=sigmaLimit,
                        iteration=iteration,
                        brightest=brightest,
                        exclude_border=exclude_border,
                        returnAll=returnAll,
                    )

                    iteration += 1

                    # USE IRAF STAR FINDER FOR SUBSEQUENT ITERATIONS
                    iraf = True

                    # INITIALIZE SHIFTED DETECTOR COLUMNS ON FIRST PASS
                    if "detector_x_shifted" not in orderPixelTable.columns:
                        orderPixelTable["detector_x_shifted"] = orderPixelTable["detector_x"]
                        orderPixelTable["detector_y_shifted"] = orderPixelTable["detector_y"]

                    # EXPLODE MULTIPLE DETECTIONS INTO SEPARATE ROWS
                    orderPixelTable = self._explode_multiple_detections(orderPixelTable)

                    # CALCULATE POSITION DIFFERENCES BETWEEN PREDICTED AND OBSERVED
                    orderPixelTable = self._calculate_position_differences(orderPixelTable)

                    # STEP 5: APPLY SHIFT CORRECTIONS ORDER-BY-ORDER
                    if self.arm.upper() == "VIS" and self.inst.upper() == "SOXS":
                        # VIS ARM OF SOXS: GROUP BY SHIFT GROUPS
                        orderPixelTable["shift_group"] = orderPixelTable["order"]
                    else:
                        # NIR ARM OF SOXS OR OTHER INSTRUMENTS
                        # ASSIGN SHIFT GROUP BASED ON NxN GRID OF mph_mean_x AND mph_mean_y
                        # SET GRID SIZE (N=3 MEANS 3x3 GRID)
                        grid_size = 2

                        x = orderPixelTable["mph_mean_x"]
                        y = orderPixelTable["mph_mean_y"]
                        # COMPUTE BIN EDGES
                        x_bins = np.linspace(x.min(), x.max(), grid_size + 1)
                        y_bins = np.linspace(y.min(), y.max(), grid_size + 1)
                        # DIGITIZE TO GET GRID INDICES (1 to grid_size)
                        x_idx = np.digitize(x, x_bins, right=False)
                        y_idx = np.digitize(y, y_bins, right=False)
                        # CLIP TO 1-grid_size (IN CASE OF EDGE CASES)
                        x_idx = np.clip(x_idx, 1, grid_size)
                        y_idx = np.clip(y_idx, 1, grid_size)
                        # GRID CELL NUMBER: (row-major, 1 to grid_size^2)
                        orderPixelTable["shift_group"] = (y_idx - 1) * grid_size + x_idx

                    shiftGroups = orderPixelTable["shift_group"].unique()

                    for sg in shiftGroups:
                        mask = (orderPixelTable["shift_group"] == sg) & (orderPixelTable["observed_x"].notnull())

                        # CHECK FOR LARGE SHIFTS IN MULTI-PINHOLE FRAMES
                        if self.firstGuessMap and bigShift is False:
                            orderPixelTable = self._handle_multipin_hole_big_shift(orderPixelTable, sg, mask, iteration)
                            continue

                        # FIND AND APPLY CLUSTER-BASED SHIFT CORRECTION
                        orderPixelTable, centrex, centrey = self._find_and_apply_cluster_shift(orderPixelTable, mask)

                        # CALCULATE SHIFT STATISTICS FOR QC REPORTING
                        medianx, mediany, stdx, stdy, medianxy, stdxy = self._calculate_order_shift_statistics(
                            orderPixelTable, mask
                        )

                        # PRINT DEBUG INFORMATION IF ENABLED
                        if self.debug:
                            print(f"# SHIFT GROUP {sg}")
                            print(f"StdXY: {stdxy:.3f}, MedianXY: {medianxy:.3f}")
                            print(f"StdX: {stdx:.3f}, StdY: {stdy:.3f}")
                            print(f"MedianX: {medianx:.3f}, MedianY: {mediany:.3f}")
                            if centrex:
                                print(f"centrex: {centrex:.3f}, centrey: {centrey:.3f}")
                            else:
                                print("no centre found")
                            print(f"Length shift group {sg}:", len(orderPixelTable.loc[mask].index))
                            _plot_slit_index_comparisons(orderPixelTable.loc[mask])
                            print()

                        # SKIP ORDER IF NO VALID LINES FOUND
                        if np.isnan(medianx) or np.isnan(mediany):
                            self.log.warning(f"Could not find any arc lines in shift group {sg}.")
                            continue

                        # PRINT PROGRESS UPDATE
                        sys.stdout.flush()
                        sys.stdout.write("\x1b[1A\x1b[2K")
                        self.log.print(
                            f"\t ITERATION {iteration}: Median X Y difference between predicted "
                            f"and measured positions: {medianx:0.5f},{mediany:0.5f} "
                            f"(stdx,stdy = {stdx:0.2f},{stdy:0.2f}) (shift group {sg})"
                        )

                    # RETURN TO DISCRETE PREDICTED LINE POSITIONS BY REMOVING DUPLICATE DETECTOR_X/Y VALUES
                    orderPixelTable.drop_duplicates(
                        subset=["detector_x", "detector_y"], keep="first", inplace=True, ignore_index=True
                    )
                    bigShift = True

            # MAKE A CLEAN COPY OF THE DETECTION TABLE ... USED FOR PIPELINE TUNING ONLY
            lineDetectionTable = orderPixelTable.copy()

            # COLLECT MISSING LINES
            mask = orderPixelTable["observed_x"].isnull()
            missingLines = orderPixelTable.loc[mask]
            # GROUP RESULTS BY WAVELENGTH
            lineGroups = missingLines.groupby(["wavelength", "order"])
            lineGroups = lineGroups.size().to_frame(name="count").reset_index()

            # CREATE THE LIST OF INCOMPLETE MULTIPINHOLE WAVELENGTHS & ORDER SETS TO DROP
            orderPixelTable["dropped"] = False
            if "SPAT" in self.recipeName.upper():
                missingLineThreshold = 9 - self.recipeSettings["mph_line_set_min"]
                mask = lineGroups["count"] > missingLineThreshold
                lineGroups = lineGroups.loc[mask]
                setsToDrop = lineGroups[["wavelength", "order"]]
                s = orderPixelTable[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
                s["dropped"] = False
                s.loc[(s["_merge"] == "both"), "dropped"] = True
                orderPixelTable["droppedOnMissing"] = s["dropped"].values
                orderPixelTable.loc[(orderPixelTable["droppedOnMissing"] == True), "dropped"] = True

            # DROP MISSING VALUES
            orderPixelTable.dropna(axis="index", how="any", subset=["observed_x"], inplace=True)

            detectedLines = len(orderPixelTable.index)

            if self.firstGuessMap and False:
                orderPixelTable = (
                    orderPixelTable.groupby(["wavelength", "order"]).apply(straighten_mph_sets).reset_index(drop=True)
                )

            # CALCULATE LINE DETECTION STATISTICS
            percentageDetectedLines = float("{:.6f}".format(float(detectedLines) / float(totalLines)))

            if percentageDetectedLines < 0.5:
                message = f"Only {percentageDetectedLines*100:.2f}% of the input line-list lines were detected on the pinhole frame. Please check the quality of the raw frames."
                raise ValueError(message)

            # GET CURRENT UTC TIMESTAMP FOR QC RECORDS
            utcnow = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%S")

            # WRITE LINE DETECTION QC METRICS
            self._write_qc_metrics(totalLines, detectedLines, percentageDetectedLines, utcnow)

            # GROUP FOUND LINES INTO SETS AND CLIP ON MEAN XY SHIFT RESIDUALS
            if True:
                orderPixelTable = self._clip_on_measured_line_metrics(orderPixelTable)
            else:
                orderPixelTable["fwhm_pin_px"] = 3.0

            if self.firstGuessMap:
                detectedLineGroups = orderPixelTable.groupby(["wavelength", "order"]).size().reset_index(name="count")
                # FIND THE MODE COUNT FOR EACH ORDER

                modeCounts = (
                    detectedLineGroups.groupby("order")["count"]
                    .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else 0)
                    .reset_index(name="pinhole count")
                )
                if self.inst.upper() == "SOXS":
                    # REMOVE ORDER = 10
                    modeCounts = modeCounts[~modeCounts["order"].isin([10, 11, 12])]
                # IF ANY MODE IS LESS THAN 9, FAIL GRACEFULLY
                if (modeCounts["pinhole count"] < 9).any():
                    from tabulate import tabulate

                    self.log.print("\n")
                    self.log.print(tabulate(modeCounts, headers="keys", tablefmt="psql"))
                    self.log.print("\n")
                    raise AttributeError(
                        f"All 9 pinholes cannot be seen in the data. A dispersion solution cannot be created."
                    )

            # STEP 6: ITERATIVELY FIT POLYNOMIAL SOLUTIONS WITH AUTOMATIC DEGREE REDUCTION
            fitFound = False
            tryCount = 0
            while not fitFound and tryCount < 5:

                # ATTEMPT TO FIT POLYNOMIAL DISPERSION SOLUTION
                popt_x, popt_y, goodLinesTable, clippedLinesTable = self.fit_polynomials(
                    orderPixelTable=orderPixelTable,
                    wavelengthDeg=wavelengthDeg,
                    orderDeg=orderDeg,
                    slitDeg=slitDeg,
                    missingLines=missingLines,
                )

                # CHECK IF FIT CONVERGED SUCCESSFULLY
                if isinstance(popt_x, np.ndarray) and isinstance(popt_y, np.ndarray):
                    fitFound = True
                else:
                    # IN TUNING MODE, FAIL IMMEDIATELY WITHOUT RETRY
                    if self.settings["tune-pipeline"]:
                        raise ArithmeticError(
                            "Could not converge on a good fit to the dispersion solution. "
                            "Please check the quality of your data or adjust your fitting parameters."
                        )

                    # REDUCE POLYNOMIAL DEGREES AND RETRY
                    wavelengthDeg, orderDeg, slitDeg = self._reduce_polynomial_degrees(
                        popt_x, popt_y, wavelengthDeg, orderDeg, slitDeg
                    )

                    self.log.print(
                        f"Wavelength, Order and Slit fitting orders reduced to "
                        f"{wavelengthDeg}, {orderDeg}, {slitDeg} to try and achieve a dispersion solution."
                    )

                    tryCount += 1
                    if tryCount == 5:
                        self.log.error(
                            "Could not converge on a good fit to the dispersion solution after 5 attempts. "
                            "Please check the quality of your data or adjust your fitting parameters."
                        )
                        raise ArithmeticError(
                            "Could not converge on a good fit to the dispersion solution. "
                            "Please check the quality of your data or adjust your fitting parameters."
                        )

            if bootstrap_dispersion_solution:
                # WRITE THE MAP TO FILE
                mapPath = self.write_map_to_file(popt_x, popt_y, orderDeg, wavelengthDeg, slitDeg)
                # orderPixelTable = self.update_static_line_list_detector_positions(originalOrderPixelTable, mapPath)
                orderPixelTable = self.create_new_static_line_list(dispersionMapPath=mapPath)

                boost = True
                bootstrap_dispersion_solution = False

        # STEP 7: GENERATE OUTPUT FILENAMES
        goodLinesFN, missingLinesFN = self._get_output_filenames()

        # STEP 8: PREPARE LINE LISTS WITH PROPER COLUMNS AND QC METRICS
        self._write_line_list_qc(clippedLinesTable, utcnow)
        goodAndClippedLines, goodLinesTable = self._prepare_line_list_columns(goodLinesTable, clippedLinesTable)

        # STEP 9: WRITE FITTED LINES TO FILE
        self._write_fitted_lines_file(goodAndClippedLines, goodLinesFN, utcnow)

        # STEP 10: WRITE MISSING LINES TO FILE
        self._write_missing_lines_file(missingLines, missingLinesFN, utcnow)

        # WRITE THE MAP TO FILE
        if not self.settings["tune-pipeline"]:
            mapPath = self.write_map_to_file(popt_x, popt_y, orderDeg, wavelengthDeg, slitDeg)
        else:
            mapPath = None
            self.create2DMap = False

        if self.firstGuessMap and self.orderTable and self.create2DMap:
            mapImagePath = self.map_to_image(dispersionMapPath=mapPath, orders=list(goodLinesTable["order"].unique()))
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
                dispMapImage=mapImagePath,
            )

            # if self.CLINE > 160 and self.arm == "VIS":
            #     raise ValueError(f"Total number of clipped lines {self.CLINE} is above threshold.")
            # elif self.CLINE > 2000 and self.arm == "NIR":
            #     raise ValueError(f"Total number of clipped lines {self.CLINE} is above threshold.")
            return mapPath, mapImagePath, res_plots, self.qc, self.products, lineDetectionTable

        res_plots = self._create_dispersion_map_qc_plot(
            xcoeff=popt_x,
            ycoeff=popt_y,
            orderDeg=orderDeg,
            wavelengthDeg=wavelengthDeg,
            slitDeg=slitDeg,
            orderPixelTable=goodLinesTable,
            missingLines=missingLines,
            allClippedLines=clippedLinesTable,
        )

        self.log.debug("completed the ``get`` method")
        return mapPath, None, res_plots, self.qc, self.products, lineDetectionTable

    def get_predicted_line_list(self):
        """*lift the predicted line list from the static calibrations*

        **Return:**

        - ``orderPixelTable`` -- a panda's data-frame containing wavelength,order,slit_index,slit_position,detector_x,detector_y
        """
        self.log.debug("starting the ``get_predicted_line_list`` method")

        from astropy.table import Table
        from astropy.stats import sigma_clipped_stats
        import numpy as np

        # DETERMINE FRAME TYPE (SINGLE OR MULTI-PINHOLE)
        frameTech = self._determine_frame_tech()

        # GET BINNING PARAMETERS
        binx, biny = self._get_binning_params()

        # LOAD PREDICTED LINE LIST FROM CALIBRATIONS
        orderPixelTable = self._load_predicted_lines(frameTech, binx, biny)

        # CLEAN AND PREPARE LINE LIST
        orderPixelTable = self._clean_line_list(orderPixelTable)

        # APPLY COORDINATE TRANSFORMATIONS
        orderPixelTable = self._apply_coordinate_transforms(orderPixelTable)

        # FILTER TO SPECIFIC SLIT INDEX IF NOT FIRST GUESS MAP
        if not self.firstGuessMap:
            orderPixelTable = self._filter_to_mid_slit(orderPixelTable)

        # APPLY FIRST GUESS CORRECTIONS IF AVAILABLE
        if self.firstGuessMap:
            orderPixelTable = self._apply_first_guess_corrections(orderPixelTable)

        self.log.debug("completed the ``get_predicted_line_list`` method")
        return orderPixelTable

    def _determine_frame_tech(self):
        """*Determine if frame is single or multi-pinhole*"""
        kw = self.kw
        tech = self.pinholeFrame.header[kw("DPR_TECH")]

        if tech == "ECHELLE,PINHOLE":
            return "single"
        elif tech == "ECHELLE,MULTI-PINHOLE":
            return "multi"
        else:
            raise TypeError("The input frame needs to be a calibrated single- or multi-pinhole arc lamp frame")

    def _get_binning_params(self):
        """*Extract binning parameters from frame header*"""
        kw = self.kw

        if self.arm != "NIR" and kw("WIN_BINX") in self.pinholeFrame.header:
            binx = int(self.pinholeFrame.header[kw("WIN_BINX")])
            biny = int(self.pinholeFrame.header[kw("WIN_BINY")])
        else:
            binx, biny = 1, 1

        return binx, biny

    def _load_predicted_lines(self, frameTech, binx, biny):
        """*Load predicted line list from calibration files*"""
        from astropy.table import Table

        calibrationPath = get_calibrations_path(log=self.log, settings=self.settings)
        predictedFile = self.detectorParams["predicted pinhole lines"][frameTech][f"{binx}x{biny}"]
        fullPath = f"{calibrationPath}/{predictedFile}"

        # LOAD AND CONVERT TO PANDAS
        dat = Table.read(fullPath, format="fits")
        return dat.to_pandas()

    def _clean_line_list(self, df):
        """*Clean and standardize line list dataframe*"""
        # STANDARDIZE COLUMN NAMES
        df.columns = [c.lower() if c.lower() in ["order", "wavelength"] else c for c in df.columns]

        # SET COLUMN DTYPES
        floatCols = ["wavelength", "slit_position", "detector_x", "detector_y", "order"]
        intCols = ["slit_index"]
        stringCols = ["ion"]

        for col in floatCols:
            if col in df.columns:
                try:
                    df[col] = df[col].astype(float)
                except:
                    pass

        for col in intCols:
            if col in df.columns:
                try:
                    df[col] = df[col].astype(int)
                except:
                    pass

        for col in stringCols:
            if col in df.columns:
                try:
                    df[col] = df[col].astype(str)
                except:
                    pass

        # REMOVE FLAGGED LINES
        if "delete" in df.columns:
            mask = df["delete"] == 1
            df.drop(index=df[mask].index, inplace=True)

        # REMOVE DUPLICATES
        df.drop_duplicates(inplace=True)

        return df

    def _apply_coordinate_transforms(self, df):
        """*Apply coordinate system transformations (FITS to Python indexing)*"""
        # PHOTUTILS: BOTTOM LEFT PIXEL CENTER IS (0,0)
        # FITS WCS: BOTTOM LEFT PIXEL CENTER IS (1,1)
        if self.inst != "SOXS":
            df["detector_x"] -= 1.0
            df["detector_y"] -= 1.0
        elif self.arm in ["NIR", "VIS"]:
            # CALCULATE SCIENCE PIXEL DIMENSIONS FOR SOXS
            dp = self.detectorParams
            science_pixels = dp["science-pixels"]
            xlen = science_pixels["columns"]["end"] - science_pixels["columns"]["start"]
            ylen = science_pixels["rows"]["end"] - science_pixels["rows"]["start"]

        return df

    def _filter_to_mid_slit(self, df):
        """*Filter line list to only include mid-slit position*"""
        slitIndex = int(self.detectorParams["mid_slit_index"])
        mask = df["slit_index"] == slitIndex
        return df.loc[mask]

    def _apply_first_guess_corrections(self, df):
        """*Apply systematic shifts from first guess dispersion solution*"""
        from astropy.stats import sigma_clipped_stats
        import numpy as np

        slitIndex = int(self.detectorParams["mid_slit_index"])

        # GET PREDICTED PIXEL VALUES FROM FIRST GUESS MAP
        df = dispersion_map_to_pixel_arrays(log=self.log, dispersionMapPath=self.firstGuessMap, orderPixelTable=df)

        # CALCULATE SHIFTS BETWEEN PREDICTED AND ACTUAL POSITIONS
        tmpList = df.copy()
        mask = tmpList["slit_index"] == slitIndex
        tmpList.loc[mask, "shift_x"] = tmpList.loc[mask, "detector_x"].values - tmpList.loc[mask, "fit_x"].values
        tmpList.loc[mask, "shift_y"] = tmpList.loc[mask, "detector_y"].values - tmpList.loc[mask, "fit_y"].values
        tmpList.loc[mask, "shift_xy"] = np.sqrt(
            tmpList.loc[mask, "shift_x"].values ** 2 + tmpList.loc[mask, "shift_y"].values ** 2
        )

        # EXTRACT SHIFT COLUMNS
        tmpList = tmpList.loc[tmpList["shift_x"].notnull(), ["wavelength", "order", "shift_x", "shift_y", "shift_xy"]]

        # CALCULATE MEDIAN SHIFTS PER ORDER, REPLACING OUTLIERS
        uniqueOrders = df["order"].unique()
        for o in uniqueOrders:
            mask = tmpList["order"] == o
            meanxy, medianxy, stdxy = sigma_clipped_stats(
                tmpList.loc[mask, "shift_xy"], sigma=10.0, stdfunc="mad_std", cenfunc="median", maxiters=3
            )
            meanx, medianx, stdx = sigma_clipped_stats(
                tmpList.loc[mask, "shift_x"], sigma=10.0, stdfunc="mad_std", cenfunc="median", maxiters=1
            )
            meany, mediany, stdy = sigma_clipped_stats(
                tmpList.loc[mask, "shift_y"], sigma=10.0, stdfunc="mad_std", cenfunc="median", maxiters=1
            )

            # USE MEDIAN FOR ORDERS WITH HIGH SCATTER
            if stdxy > 1.5:
                tmpList.loc[mask, "shift_y"] = mediany
                tmpList.loc[mask, "shift_x"] = medianx

        # MERGE SHIFTS BACK INTO MAIN DATAFRAME
        df = df.merge(tmpList, on=["wavelength", "order"], how="outer")

        # DROP ROWS WITHOUT VALID SHIFTS
        df.dropna(axis="index", how="any", subset=["shift_x"], inplace=True)

        # APPLY SHIFTS TO DETECTOR POSITIONS
        df.loc[:, "detector_x"] -= df.loc[:, "shift_x"]
        df.loc[:, "detector_y"] -= df.loc[:, "shift_y"]

        # DROP TEMPORARY SHIFT COLUMNS
        df.drop(columns=["fit_x", "fit_y", "shift_x", "shift_y"], inplace=True)

        # REMOVE INCOMPLETE MULTI-PINHOLE SETS
        df = self._remove_incomplete_mph_sets(df)

        return df

    def _remove_incomplete_mph_sets(self, df):
        """*Remove multi-pinhole line sets that don't contain all slit positions*"""
        # GROUP BY WAVELENGTH AND ORDER
        lineGroups = df.groupby(["wavelength", "order"]).size().to_frame(name="count").reset_index()
        fullSet = lineGroups["count"].max()

        # IDENTIFY INCOMPLETE SETS
        mask = lineGroups["count"] < fullSet
        setsToDrop = lineGroups.loc[mask, ["wavelength", "order"]]

        # MERGE AND FLAG INCOMPLETE SETS
        s = df[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
        s["dropped"] = s["_merge"] == "both"
        df["droppedOnMissing"] = s["dropped"].values

        # FILTER OUT DROPPED SETS
        df = df.loc[df["droppedOnMissing"] == False]
        df.drop(columns=["droppedOnMissing"], inplace=True)

        return df

    def _plot_predicted_line_windows(self, orderPixelTable, pinholeFrameMasked):
        """*Plot predicted line positions with detection windows overlaid*"""
        import matplotlib.pyplot as plt
        from matplotlib.patches import Rectangle

        fig, ax = plt.subplots()
        ax.imshow(
            pinholeFrameMasked,
            cmap="gray",
            origin="lower",
            vmin=self.meanFrameFlux,
            vmax=self.meanFrameFlux + 5 * self.stdFrameFlux,
        )

        # OVERLAY DETECTION WINDOWS
        for index, row in orderPixelTable.iterrows():
            x, y = row["detector_x"], row["detector_y"]
            ax.add_patch(
                Rectangle(
                    (x - self.windowSize / 2, y - self.windowSize / 2),
                    self.windowSize,
                    self.windowSize,
                    fill=None,
                    edgecolor="red",
                )
            )
        plt.show()

    def _get_output_filenames(self):
        """*Generate output filenames for line lists*"""
        if not self.sofName:
            filename = filenamer(log=self.log, frame=self.pinholeFrame, settings=self.settings)
            goodLinesFN = filename.replace(".fits", "_FITTED_LINES.fits")
            missingLinesFN = filename.replace(".fits", "_MISSED_LINES.fits")
        else:
            goodLinesFN = self.sofName + "_FITTED_LINES.fits"
            missingLinesFN = self.sofName + "_MISSED_LINES.fits"

        return goodLinesFN, missingLinesFN

    def _write_line_list_qc(self, clippedLinesTable, utcnow):
        """*Write QC metric for number of clipped lines*"""
        import pandas as pd

        self.CLINE = len(clippedLinesTable.index)

        self.qc = pd.concat(
            [
                self.qc,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "qc_name": "CLINE",
                        "qc_value": self.CLINE,
                        "qc_comment": "Total number of detected lines clipped during solution fitting",
                        "qc_unit": "lines",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "to_header": True,
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

    def _prepare_line_list_columns(self, goodLinesTable, clippedLinesTable):
        """*Prepare and combine good and clipped line lists with proper columns*"""
        import pandas as pd
        import numpy as np

        # DEFINE COLUMNS TO KEEP
        keepColumns = [
            "wavelength",
            "order",
            "slit_index",
            "slit_position",
            "detector_x",
            "detector_y",
            "observed_x",
            "observed_y",
            "x_diff",
            "y_diff",
            "fit_x",
            "fit_y",
            "residuals_x",
            "residuals_y",
            "residuals_xy",
            "sigma_clipped",
            "sharpness",
            "roundness1",
            "roundness2",
            "npix",
            "sky",
            "peak",
            "flux",
            "fwhm_pin_px",
            "R_pin",
            "pixelScaleNm",
            "detector_x_shifted",
            "detector_y_shifted",
            "R_slit",
            "fwhm_slit_px",
        ]

        # ADD ION COLUMN IF PRESENT
        if "ion" in goodLinesTable.columns:
            keepColumns.insert(0, "ion")

        # MARK CLIPPED LINES
        clippedLinesTable["sigma_clipped"] = True
        clippedLinesTable["R"] = np.nan
        clippedLinesTable["pixelScaleNm"] = np.nan

        # COMBINE GOOD AND CLIPPED LINES
        if len(clippedLinesTable.index):
            goodAndClippedLines = pd.concat(
                [clippedLinesTable[keepColumns], goodLinesTable[keepColumns]], ignore_index=True
            )
        else:
            goodAndClippedLines = goodLinesTable[keepColumns]

        # SORT BOTH DATAFRAMES
        goodAndClippedLines.sort_values(["order", "wavelength", "slit_index"], inplace=True)
        goodLinesTable = goodLinesTable[keepColumns]
        goodLinesTable.sort_values(["order", "wavelength", "slit_index"], inplace=True)

        return goodAndClippedLines, goodLinesTable

    def _write_fitted_lines_file(self, goodAndClippedLines, goodLinesFN, utcnow):
        """*Write fitted lines (good + clipped) to FITS file*"""
        from astropy.table import Table
        import pandas as pd

        t = Table.from_pandas(goodAndClippedLines)
        filePath = f"{self.qcDir}/{goodLinesFN}"

        if not self.settings["tune-pipeline"]:
            t.write(filePath, overwrite=True)

        self.products = pd.concat(
            [
                self.products,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "product_label": "DISP_MAP_LINES",
                        "file_name": goodLinesFN,
                        "file_type": "FITS",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "product_desc": f"{self.arm} dispersion solution fitted lines",
                        "file_path": filePath,
                        "label": "QC",
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

    def _write_missing_lines_file(self, missingLines, missingLinesFN, utcnow):
        """*Write missing lines to FITS file*"""
        from astropy.table import Table
        import pandas as pd

        keepColumns = ["wavelength", "order", "slit_index", "slit_position", "detector_x", "detector_y"]

        # SORT AND EXTRACT RELEVANT COLUMNS
        missingLines.sort_values(["order", "wavelength", "slit_index"], inplace=True)
        t = Table.from_pandas(missingLines[keepColumns])
        filePath = f"{self.qcDir}/{missingLinesFN}"

        if not self.settings["tune-pipeline"]:
            t.write(filePath, overwrite=True)

        self.products = pd.concat(
            [
                self.products,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "product_label": "DISP_MAP_LINES_MISSING",
                        "file_name": missingLinesFN,
                        "file_type": "FITS",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "product_desc": f"{self.arm} undetected arc lines",
                        "file_path": filePath,
                        "label": "QC",
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

    def detect_pinhole_arc_lines(
        self,
        orderPixelTable,
        iraf=True,
        sigmaLimit=3,
        iteration=False,
        brightest=False,
        exclude_border=False,
        returnAll=False,
    ):
        """*detect the observed position of an arc-line given the predicted pixel positions*

        **Key Arguments:**

        - ``orderPixelTable`` -- the initial line list
        - ``iraf`` -- use IRAF star finder to generate a FWHM
        - ``sigmaLimit`` -- the lower sigma limit for arc line to be considered detected
        - ``iteration`` -- which detect and shift iteration are we on?
        - ``brightest`` -- find the brightest source?
        - ``exclude_border`` -- exclude border pixels from detection?

        **Return:**

        - ``predictedLine`` -- the line with the observed pixel coordinates appended (if detected, otherwise nan)
        """
        self.log.debug("starting the ``detect_pinhole_arc_lines`` method")

        import numpy as np
        from fundamentals import fmultiprocess
        import logging

        # FIX ASTROPY LOGGING LEVEL RESET ISSUE
        logging.getLogger().setLevel(logging.INFO + 5)

        # GET FRAME AND WINDOW PARAMETERS
        pinholeFrame = self.pinholeFrameMasked
        windowHalf = self.windowHalf

        # USE SHIFTED POSITIONS IF AVAILABLE, OTHERWISE USE ORIGINAL PREDICTIONS
        if "detector_x_shifted" in orderPixelTable.columns:
            xArray = orderPixelTable["detector_x_shifted"]
            yArray = orderPixelTable["detector_y_shifted"]
        else:
            xArray = orderPixelTable["detector_x"]
            yArray = orderPixelTable["detector_y"]

        # CALCULATE STAMP BOUNDARIES AROUND EACH PREDICTED POSITION
        # ENSURE BOUNDARIES DON'T EXCEED FRAME DIMENSIONS
        xlows = np.clip(np.round(xArray - windowHalf).astype(int), 0, None)
        xups = np.clip(np.round(xArray + windowHalf).astype(int), None, pinholeFrame.shape[1])
        ylows = np.clip(np.round(yArray - windowHalf).astype(int), 0, None)
        yups = np.clip(np.round(yArray + windowHalf).astype(int), None, pinholeFrame.shape[0])

        # CREATE IMAGE STAMPS (CUTOUTS) FOR EACH PREDICTED LINE POSITION
        stamps = [
            (pinholeFrame[ylow:yup, xlow:xup], xlow, xup, ylow, yup)
            for xlow, xup, ylow, yup in zip(xlows, xups, ylows, yups)
        ]

        # INITIALIZE OUTPUT DICTIONARY FOR LINE MEASUREMENTS
        predictedLines = {
            "sharpness": [],  # POINT SOURCE SHARPNESS METRIC
            "roundness1": [],  # ROUNDNESS METRIC 1
            "roundness2": [],  # ROUNDNESS METRIC 2
            "npix": [],  # NUMBER OF PIXELS IN DETECTION
            "sky": [],  # LOCAL SKY BACKGROUND LEVEL
            "peak": [],  # PEAK PIXEL VALUE
            "flux": [],  # INTEGRATED FLUX
            "observed_x": [],  # MEASURED X POSITION
            "observed_y": [],  # MEASURED Y POSITION
            "fwhm_pin_px": [],  # MEASURED FWHM IN PIXELS
        }

        # PROCESS ALL STAMPS USING MULTIPROCESSING (OR SERIAL IN DEBUG MODE)
        results = fmultiprocess(
            log=self.log,
            function=measure_line_position,
            inputArray=stamps,
            poolSize=False,
            timeout=300,
            mute=False,
            progressBar=False,
            turnOffMP=self.debug,  # DISABLE MULTIPROCESSING IN DEBUG MODE
            windowHalf=self.windowHalf,
            iraf=iraf,
            sigmaLimit=sigmaLimit,
            brightest=brightest,
            exclude_border=exclude_border,
            iteration=iteration,
            returnAll=returnAll,
            debug=self.debug,
        )

        # AGGREGATE RESULTS FROM ALL STAMP MEASUREMENTS
        for rr in results:
            for k in predictedLines.keys():
                thisList = [r[k] for r in rr]
                predictedLines[k].append(thisList)
        print()

        # ADD MEASUREMENT RESULTS TO ORDER PIXEL TABLE
        for k, v in predictedLines.items():
            orderPixelTable[k] = v

        self.log.debug("completed the ``detect_pinhole_arc_lines`` method")
        return orderPixelTable

    def write_map_to_file(self, xcoeff, ycoeff, orderDeg, wavelengthDeg, slitDeg):
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
        self.log.debug("starting the ``write_map_to_file`` method")

        import pandas as pd
        from astropy.table import Table
        from astropy.io import fits
        from contextlib import suppress
        import copy
        import math

        arm = self.arm
        kw = self.kw

        # ORGANIZE X-AXIS POLYNOMIAL COEFFICIENTS INTO DICTIONARY
        # HANDLE BOTH SINGLE AND PER-AXIS POLYNOMIAL DEGREES
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
                        coeff_dict_x[f"c{i}{j}{k}"] = xcoeff[n_coeff]
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
                        coeff_dict_x[f"c{i}{j}{k}"] = xcoeff[n_coeff]
                        n_coeff += 1

        # ORGANIZE Y-AXIS POLYNOMIAL COEFFICIENTS INTO DICTIONARY
        # HANDLE BOTH SINGLE AND PER-AXIS POLYNOMIAL DEGREES
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
                        coeff_dict_y[f"c{i}{j}{k}"] = ycoeff[n_coeff]
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
                        coeff_dict_y[f"c{i}{j}{k}"] = ycoeff[n_coeff]
                        n_coeff += 1

        # GENERATE OUTPUT FILENAME FROM FRAME OR SOF NAME
        if not self.sofName:
            filename = filenamer(log=self.log, frame=self.pinholeFrame, settings=self.settings)
        else:
            filename = self.sofName + ".fits"

        # PREPARE HEADER BY COPYING FRAME HEADER AND REMOVING UNWANTED KEYWORDS
        header = copy.deepcopy(self.pinholeFrame.header)
        header.pop(kw("DPR_TECH"))
        header.pop(kw("DPR_CATG"))
        header.pop(kw("DPR_TYPE"))

        # REMOVE DETECTOR-SPECIFIC KEYWORDS (MAY NOT EXIST IN ALL HEADERS)
        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError, LookupError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        # SET PROCESSING TECHNIQUE BASED ON FRAME TYPE
        if slitDeg == 0 or slitDeg == [0, 0]:
            # SINGLE PINHOLE FRAME
            header[kw("PRO_TECH").upper()] = "ECHELLE,PINHOLE"
        else:
            # MULTI-PINHOLE FRAME
            header[kw("PRO_TECH").upper()] = "ECHELLE,MULTI-PINHOLE"
        filePath = f"{self.productDir}/{filename}"

        # CREATE BINARY TABLE WITH COEFFICIENT DICTIONARIES
        df = pd.DataFrame([coeff_dict_x, coeff_dict_y])
        t = Table.from_pandas(df)

        # SORT COLUMNS: METADATA FIRST, THEN COEFFICIENTS ALPHABETICALLY
        cols = list(t.columns)
        startList = ["axis", "order_deg", "wavelength_deg", "slit_deg"]
        cols = [c for c in cols if c not in startList]
        cols.sort()
        cols = startList + cols
        t = t[cols]

        # CONVERT TO FITS BINARY TABLE HDU
        BinTableHDU = fits.table_to_hdu(t)

        # UPDATE HEADER WITH STANDARD PROCESSING KEYWORDS
        header[kw("SEQ_ARM").upper()] = arm
        header[kw("PRO_TYPE").upper()] = "REDUCED"
        header[kw("PRO_CATG").upper()] = f"DISP_TAB_{arm}".upper()
        self.dispMapHeader = header

        # ADD QC METRICS TO HEADER
        for n, v, c, h in zip(
            self.qc["qc_name"].values,
            self.qc["qc_value"].values,
            self.qc["qc_comment"].values,
            self.qc["to_header"].values,
        ):
            if h:
                header[f"ESO QC {n}".upper()] = (v, c)

        # CREATE PRIMARY HDU AND WRITE FITS FILE
        priHDU = fits.PrimaryHDU(header=header)
        hduList = fits.HDUList([priHDU, BinTableHDU])
        hduList.writeto(filePath, checksum=True, overwrite=True)

        if False:
            # MAKE RELATIVE HOME PATH ABSOLUTE
            from os.path import expanduser

            home = expanduser("~")

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

        self.log.debug("completed the ``write_map_to_file`` method")
        return filePath

    def calculate_residuals(
        self, orderPixelTable, xcoeff, ycoeff, orderDeg, wavelengthDeg, slitDeg, writeQCs=False, pixelRange=False
    ):
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
        self.log.debug("starting the ``calculate_residuals`` method")

        import numpy as np
        import pandas as pd

        arm = self.arm

        utcnow = datetime.now(timezone.utc)
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
            log=self.log,
            orderDeg=orderDegx,
            wavelengthDeg=wavelengthDegx,
            slitDeg=slitDegx,
            exponentsIncluded=True,
            axis="x",
        ).poly
        polyy = chebyshev_order_wavelength_polynomials(
            log=self.log,
            orderDeg=orderDegy,
            wavelengthDeg=wavelengthDegy,
            slitDeg=slitDegy,
            exponentsIncluded=True,
            axis="y",
        ).poly

        # CALCULATE X & Y RESIDUALS BETWEEN OBSERVED LINE POSITIONS AND POLY
        # FITTED POSITIONS
        orderPixelTable["fit_x"] = polyx(orderPixelTable, *xcoeff)
        orderPixelTable["fit_y"] = polyy(orderPixelTable, *ycoeff)

        if pixelRange == True:
            polyx = chebyshev_order_wavelength_polynomials(
                log=self.log,
                orderDeg=orderDegx,
                wavelengthDeg=wavelengthDegx,
                slitDeg=slitDegx,
                exponentsIncluded=False,
                axis="x",
            ).poly
            polyy = chebyshev_order_wavelength_polynomials(
                log=self.log,
                orderDeg=orderDegy,
                wavelengthDeg=wavelengthDegy,
                slitDeg=slitDegy,
                exponentsIncluded=False,
                axis="y",
            ).poly
            # GET THE PIXEL SCALE
            orderPixelTableHigh = orderPixelTable.copy()
            nmRange = 4.0
            orderPixelTableHigh["wavelength"] = orderPixelTableHigh["wavelength"] + nmRange / 2.0
            orderPixelTableHigh["fit_x"] = polyx(orderPixelTableHigh, *xcoeff)
            orderPixelTableHigh["fit_y"] = polyy(orderPixelTableHigh, *ycoeff)

            orderPixelTableLow = orderPixelTable.copy()
            orderPixelTableLow["wavelength"] = orderPixelTableLow["wavelength"] - nmRange / 2.0
            orderPixelTableLow["fit_x"] = polyx(orderPixelTableLow, *xcoeff)
            orderPixelTableLow["fit_y"] = polyy(orderPixelTableLow, *ycoeff)

            orderPixelTable["fit_x_high"] = orderPixelTableHigh["fit_x"]
            orderPixelTable["fit_y_high"] = orderPixelTableHigh["fit_y"]
            orderPixelTable["fit_x_low"] = orderPixelTableLow["fit_x"]
            orderPixelTable["fit_y_low"] = orderPixelTableLow["fit_y"]

            orderPixelTable["pixelScaleNm"] = nmRange / np.power(
                np.power(orderPixelTable["fit_x_high"] - orderPixelTable["fit_x_low"], 2)
                + np.power(orderPixelTable["fit_y_high"] - orderPixelTable["fit_y_low"], 2),
                0.5,
            )
            orderPixelTable["delta_wavelength"] = orderPixelTable["pixelScaleNm"] * orderPixelTable["fwhm_pin_px"]
            orderPixelTable["R_pin"] = orderPixelTable["wavelength"] / orderPixelTable["delta_wavelength"]

            # REMOVE COLUMN FROM DATA FRAME
            orderPixelTable.drop(columns=["fit_x_high", "fit_y_high", "fit_x_low", "fit_y_low"], inplace=True)

        orderPixelTable["residuals_x"] = orderPixelTable["fit_x"] - orderPixelTable["observed_x"]
        orderPixelTable["residuals_y"] = orderPixelTable["fit_y"] - orderPixelTable["observed_y"]

        # CALCULATE COMBINED RESIDUALS AND STATS
        orderPixelTable["residuals_xy"] = np.sqrt(
            np.square(orderPixelTable["residuals_x"]) + np.square(orderPixelTable["residuals_y"])
        )
        combined_res_mean = np.mean(orderPixelTable["residuals_xy"])
        combined_res_std = np.std(orderPixelTable["residuals_xy"])
        combined_res_median = np.median(orderPixelTable["residuals_xy"])

        if False:
            import sqlite3 as sql

            # REGISTER SQL CONVERTERS
            sql.register_adapter(list, lambda arr: str(arr.tolist()))
            sql.register_adapter(np.array, lambda arr: str(arr.tolist()))
            sql.register_adapter(np.ndarray, lambda arr: str(arr.tolist()))
            sql.register_adapter(np.float64, lambda this: this.item())
            sql.register_adapter(np.ma.core.MaskedArray, lambda arr: str(arr.tolist()))

            # CONNECT TO THE DATABASE
            conn = sql.connect("pandas_export.db")
            # SEND TO DATABASE
            orderPixelTable.to_sql("my_export_table", con=conn, index=False, if_exists="replace")
            conn.commit()
            conn.close()
            sys.exit()

        if self.arcFrame:
            self.slitWidth = self.arcFrame.header[self.kw(f"SLIT_{arm}")].replace("SLIT", "").split("x")[0]
            orderPixelTable[["R_slit", "fwhm_slit_px"]] = orderPixelTable.apply(
                self._calculate_resolution_on_slit, axis=1
            )

        else:
            orderPixelTable["R_slit"] = 0.0
            orderPixelTable["fwhm_slit_px"] = 0.0
        if writeQCs:
            absx = abs(orderPixelTable["residuals_x"])
            absy = abs(orderPixelTable["residuals_y"])
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XRESMIN",
                            "qc_value": absx.min(),
                            "qc_comment": "[px] Minimum residual in dispersion solution fit along x-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XRESMAX",
                            "qc_value": absx.max(),
                            "qc_comment": "[px] Maximum residual in dispersion solution fit along x-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XRESRMS",
                            "qc_value": absx.std(),
                            "qc_comment": "[px] Std-dev of residual in dispersion solution fit along x-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "YRESMIN",
                            "qc_value": absy.min(),
                            "qc_comment": "[px] Minimum residual in dispersion solution fit along y-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "YRESMAX",
                            "qc_value": absy.max(),
                            "qc_comment": "[px] Maximum residual in dispersion solution fit along y-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "YRESRMS",
                            "qc_value": absy.std(),
                            "qc_comment": "[px] Std-dev of residual in dispersion solution fit along y-axis",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XYRESMIN",
                            "qc_value": orderPixelTable["residuals_xy"].min(),
                            "qc_comment": "[px] Minimum residual in dispersion solution fit (XY combined)",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XYRESMAX",
                            "qc_value": orderPixelTable["residuals_xy"].max(),
                            "qc_comment": "[px] Maximum residual in dispersion solution fit (XY combined)",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.Series(
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "XYRESRMS",
                            "qc_value": orderPixelTable["residuals_xy"].std(),
                            "qc_comment": "[px] Std-dev of residual in dispersion solution (XY combined)",
                            "qc_unit": "pixels",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    )
                    .to_frame()
                    .T,
                ],
                ignore_index=True,
            )

        self.log.debug("completed the ``calculate_residuals`` method")
        return combined_res_mean, combined_res_std, combined_res_median, orderPixelTable

    def _calculate_resolution_on_slit(self, row):
        import math
        from astropy.modeling import models, fitting
        import numpy as np
        import matplotlib.pyplot as plt
        import pandas as pd

        stdToFwhm = 2 * (2 * math.log(2)) ** 0.5

        # GOOD GUESS AT STD AND SLICE SIZE
        if self.arm == "NIR":
            guessStd = (float(self.slitWidth) / 0.25) / stdToFwhm
        else:
            guessStd = (float(self.slitWidth) / 0.29) / stdToFwhm
        sliceSize = int(guessStd * 5.0)

        if self.detectorParams["dispersion-axis"] == "y":
            detector_row = self.arcFrame.data[
                int(row["observed_y"]), int(row["observed_x"]) - sliceSize : int(row["observed_x"]) + sliceSize
            ]
            detector_row = detector_row - np.median(detector_row)
            detector_row = np.where(detector_row < 0, 0, detector_row)
        else:
            detector_row = self.arcFrame.data[
                int(row["observed_y"]) - sliceSize : int(row["observed_y"]) + sliceSize, int(row["observed_x"])
            ]
            detector_row = detector_row - np.median(detector_row)
            detector_row = np.where(detector_row < 0, 0, detector_row)

        # FITTING THE LINE WITH A GAUSSIAN PROFILE
        g_init = models.Gaussian1D(amplitude=1.0, mean=0, stddev=guessStd)
        fit_g = fitting.LevMarLSQFitter()

        try:
            g = fit_g(g_init, np.arange(0, len(detector_row)), detector_row)

            if g.stddev.value < guessStd / 2.0:
                # 1.0 SLIT CANNOT BE LESS THAN 1.0 px
                raise Exception

            # stddev_corrected = np.sqrt(g.stddev.value*g.stddev.value - np.abs(13*np.sin(row['tilt'])*np.sin(row['tilt'])))
            stddev_corrected = float(g.stddev.value)
            fwhm = stddev_corrected * stdToFwhm

            import random

            random_number = random.randint(1, 5)
            print("COME BACK HERE AND TRASH SOME LINES")
            if False and random_number == 2:
                gaussx = np.arange(0, len(detector_row), 0.05)
                plt.plot(np.arange(0, len(detector_row)), detector_row)
                stddev_corrected * stdToFwhm
                plt.plot(gaussx, g(gaussx), label=f"STD = {g.stddev.value:0.2f}, FWHM = {fwhm:0.2f}")
                plt.legend()
                plt.show()
        except Exception as e:
            return pd.Series([None, None])

        # stddev_corrected = np.sqrt(g.stddev.value*g.stddev.value - np.abs(13*np.sin(row['tilt'])*np.sin(row['tilt'])))
        stddev_corrected = float(g.stddev.value)

        # ENFORCING RESONABLE VALUES EXCLUDING OUTLIERS
        if stddev_corrected < guessStd / 2.0 or stddev_corrected > guessStd * 2.0:
            return pd.Series([None, None])

        delta_wavelength = row["pixelScaleNm"] * fwhm
        resolution_line = row["wavelength"] / delta_wavelength

        return pd.Series([resolution_line, fwhm])

    def fit_polynomials(self, orderPixelTable, wavelengthDeg, orderDeg, slitDeg, missingLines=False):
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
        self.log.debug("starting the ``fit_polynomials`` method")

        # FIRST REMOVE DROPPED LINES FILTERED ROWS FROM DATA FRAME
        allClippedLines = []
        mask = orderPixelTable["dropped"] == True
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
            log=self.log,
            orderDeg=orderDegx,
            wavelengthDeg=wavelengthDegx,
            slitDeg=slitDegx,
            exponentsIncluded=True,
            axis="x",
        ).poly
        polyy = chebyshev_order_wavelength_polynomials(
            log=self.log,
            orderDeg=orderDegy,
            wavelengthDeg=wavelengthDegy,
            slitDeg=slitDegy,
            exponentsIncluded=True,
            axis="y",
        ).poly

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
        mean_res = 100.0

        # orderPixelTable["observed_x"]=orderPixelTable["tilt_corrected_x"]
        # orderPixelTable["observed_y"]=orderPixelTable["tilt_corrected_y"]

        orderPixelTable["sigma_clipped"] = False
        while clippedCount > 0 and iteration < clippingIterationLimit:
            iteration += 1
            observed_x = orderPixelTable["observed_x"].to_numpy()
            observed_y = orderPixelTable["observed_y"].to_numpy()

            # IF mean_res > 10 WE WANT TO START FROM SCRATCH AGAIN SO NOT TO INFLUENCE THE FINAL RESULT
            if True or mean_res > 10:
                # FIND CACHED COEFF ELSE RETURN ARRAYS OF 1s
                xcoeff, ycoeff = get_cached_coeffs(
                    log=self.log,
                    arm=arm,
                    settings=self.settings,
                    recipeName=self.recipeName,
                    orderDeg=orderDeg,
                    wavelengthDeg=wavelengthDeg,
                    slitDeg=slitDeg,
                    reset=True,
                )

            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLYS
            # FIRST X
            self.log.info("""curvefit x""" % locals())

            try:
                xcoeff, pcov_x = curve_fit(polyx, xdata=orderPixelTable, ydata=observed_x, p0=xcoeff, maxfev=30000)
            except:
                return "xerror", None, None, None

            # NOW Y
            self.log.info("""curvefit y""" % locals())

            try:
                ycoeff, pcov_y = curve_fit(polyy, xdata=orderPixelTable, ydata=observed_y, p0=ycoeff, maxfev=30000)
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
                writeQCs=False,
            )

            # DO SOME CLIPPING ON THE PROFILES OF THE DETECTED LINES
            if iteration == -1:
                orderPixelTable = self._clip_on_measured_line_metrics(orderPixelTable)

            # COUNT THE CLIPPED LINES
            mask = orderPixelTable["sigma_clipped"] == True
            allClippedLines.append(orderPixelTable.loc[mask])
            mask = orderPixelTable["sigma_clipped"] == True
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
                lineGroups = orderPixelTable[columnsNoStrings].groupby(["wavelength", "order"]).mean()
                lineGroups = lineGroups.reset_index()

                # SIGMA-CLIP THE DATA ON SCATTER
                masked_residuals = sigma_clip(
                    lineGroups["residuals_x"].abs(),
                    sigma_lower=3000,
                    sigma_upper=clippingSigmaX,
                    maxiters=1,
                    cenfunc="mean",
                    stdfunc="std",
                )
                lineGroups["sigma_clipped_x"] = masked_residuals.mask
                masked_residuals = sigma_clip(
                    lineGroups["residuals_y"].abs(),
                    sigma_lower=3000,
                    sigma_upper=clippingSigmaY,
                    maxiters=1,
                    cenfunc="mean",
                    stdfunc="std",
                )
                lineGroups["sigma_clipped_y"] = masked_residuals.mask
                lineGroups.loc[
                    ((lineGroups["sigma_clipped_y"] == True) | (lineGroups["sigma_clipped_x"] == True)), "sigma_clipped"
                ] = True

                if True:
                    # CLIP ALSO ON COMBINED RESIDUALS
                    masked_residuals = sigma_clip(
                        lineGroups["residuals_xy"],
                        sigma_lower=5000,
                        sigma_upper=clippingSigma,
                        maxiters=1,
                        cenfunc="mean",
                        stdfunc="std",
                    )
                    lineGroups["sigma_clipped_xy"] = masked_residuals.mask
                    lineGroups.loc[
                        (
                            (lineGroups["sigma_clipped_y"] == True)
                            | (lineGroups["sigma_clipped_x"] == True)
                            | (lineGroups["sigma_clipped_xy"] == True)
                        ),
                        "sigma_clipped",
                    ] = True

                # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FLUX
                mask = lineGroups["sigma_clipped"] == True
                clippedGroups = lineGroups.loc[mask]
                clippedGroups = clippedGroups[["wavelength", "order"]]
                s = orderPixelTable[["wavelength", "order"]].merge(clippedGroups, indicator=True, how="left")
                s["clipped"] = False
                s.loc[(s["_merge"] == "both"), "clipped"] = True

                orderPixelTable["sigma_clipped"] = s["clipped"].values

                # CLIP THE MOST DEVIATE SINGLE PINHOLES
                if iteration == 1 and True:
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_x"].abs(),
                        sigma_lower=3000,
                        sigma_upper=clippingSigmaX * 3,
                        maxiters=1,
                        cenfunc="mean",
                        stdfunc="std",
                    )
                    orderPixelTable["sigma_clipped_x"] = masked_residuals.mask
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_y"].abs(),
                        sigma_lower=3000,
                        sigma_upper=clippingSigmaY * 3,
                        maxiters=1,
                        cenfunc="mean",
                        stdfunc="std",
                    )
                    orderPixelTable["sigma_clipped_y"] = masked_residuals.mask
                    masked_residuals = sigma_clip(
                        orderPixelTable["R_pin"],
                        sigma_lower=3000,
                        sigma_upper=10,
                        maxiters=1,
                        cenfunc="mean",
                        stdfunc="std",
                    )
                    orderPixelTable["sigma_clipped_R"] = masked_residuals.mask
                    orderPixelTable.loc[
                        (
                            (orderPixelTable["sigma_clipped_y"] == True)
                            | (orderPixelTable["sigma_clipped_x"] == True)
                            | (orderPixelTable["sigma_clipped_R"] == True)
                        ),
                        "sigma_clipped",
                    ] = True

            else:
                masked_residuals = sigma_clip(
                    orderPixelTable["residuals_x"].abs(),
                    sigma_lower=3000,
                    sigma_upper=clippingSigmaX * 1,
                    maxiters=3,
                    cenfunc="mean",
                    stdfunc="std",
                )
                orderPixelTable["sigma_clipped_x"] = masked_residuals.mask
                masked_residuals = sigma_clip(
                    orderPixelTable["residuals_y"].abs(),
                    sigma_lower=3000,
                    sigma_upper=clippingSigmaY * 1,
                    maxiters=3,
                    cenfunc="mean",
                    stdfunc="std",
                )
                orderPixelTable["sigma_clipped_y"] = masked_residuals.mask
                orderPixelTable.loc[
                    ((orderPixelTable["sigma_clipped_y"] == True) | (orderPixelTable["sigma_clipped_x"] == True)),
                    "sigma_clipped",
                ] = True

                if True:
                    # CLIP ALSO ON COMBINED RESIDUALS
                    masked_residuals = sigma_clip(
                        orderPixelTable["residuals_xy"],
                        sigma_lower=5000,
                        sigma_upper=clippingSigma,
                        maxiters=1,
                        cenfunc="mean",
                        stdfunc="std",
                    )
                    orderPixelTable["sigma_clipped_xy"] = masked_residuals.mask
                    try:
                        orderPixelTable.loc[
                            (
                                (orderPixelTable["sigma_clipped_y"] == True)
                                | (orderPixelTable["sigma_clipped_x"] == True)
                                | (orderPixelTable["sigma_clipped_xy"] == True)
                            ),
                            "sigma_clipped",
                        ] = True
                    except:
                        orderPixelTable.loc[(orderPixelTable["sigma_clipped_xy"] == True), "sigma_clipped"] = True

            # COUNT THE CLIPPED LINES
            mask = orderPixelTable["sigma_clipped"] == True
            allClippedLines.append(orderPixelTable.loc[mask])
            totalAllClippedLines = pd.concat(allClippedLines, ignore_index=True)

            # RETURN BREAKDOWN OF COLUMN VALUE COUNT
            valCounts = orderPixelTable["sigma_clipped"].value_counts(normalize=False)
            if True in valCounts:
                clippedCount = valCounts[True]
            else:
                clippedCount = 0

            if iteration > 1:
                # Cursor up one line and clear line
                sys.stdout.flush()
                sys.stdout.write("\x1b[1A\x1b[2K")

            try:
                self.log.print(
                    f"\tITERATION {iteration:02d}: {clippedCount} arc lines where clipped in this iteration of fitting a global dispersion map"
                )
            except:
                pass

            mask = orderPixelTable["sigma_clipped"] == True
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
            pixelRange=True,
        )

        self.log.debug("completed the ``fit_polynomials`` method")
        return xcoeff, ycoeff, orderPixelTable, allClippedLines

    def create_placeholder_images(self, order=False, plot=False, reverse=False):
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
        self.log.debug("starting the ``create_placeholder_images`` method")

        import numpy as np
        from astropy.nddata import CCDData

        kw = self.kw
        dp = self.detectorParams

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, extend=0.0, order=order
        )

        # CREATE THE IMAGE SAME SIZE AS DETECTOR - NAN INSIDE ORDERS, 0 OUTSIDE
        science_pixels = dp["science-pixels"]
        xlen = science_pixels["columns"]["end"] - science_pixels["columns"]["start"]
        ylen = science_pixels["rows"]["end"] - science_pixels["rows"]["start"]
        xlen, ylen

        if reverse:
            seedArray = np.empty((ylen, xlen))
            seedArray[:] = np.nan
        else:
            seedArray = np.zeros((ylen, xlen))
        wlMap = CCDData(seedArray, unit="adu")
        orderMap = wlMap.copy()
        uniqueOrders = orderPixelTable["order"].unique()
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
            axisBcoord = orderPixelTable.loc[(orderPixelTable["order"] == o)][f"{self.axisB}coord"]
            axisACoord_edgeup = (
                orderPixelTable.loc[(orderPixelTable["order"] == o)][f"{self.axisA}coord_edgeup"] + expandEdges
            )
            axisACoord_edgelow = (
                orderPixelTable.loc[(orderPixelTable["order"] == o)][f"{self.axisA}coord_edgelow"] - expandEdges
            )
            axisACoord_edgeup[axisACoord_edgeup > axisALen] = axisALen
            axisACoord_edgeup[axisACoord_edgeup < 0] = 0
            axisACoord_edgelow[axisACoord_edgelow > axisALen] = axisALen
            axisACoord_edgelow[axisACoord_edgelow < 0] = 0
            axisACoord_edgelow, axisACoord_edgeup, axisBcoord = zip(
                *[
                    (l, u, b)
                    for l, u, b in zip(axisACoord_edgelow, axisACoord_edgeup, axisBcoord)
                    if l >= 0 and l <= axisALen and u >= 0 and u <= axisALen and b >= 0 and b < axisBLen
                ]
            )
            if reverse:
                for b, u, l in zip(
                    axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)
                ):
                    if self.axisA == "x":
                        wlMap.data[b, l:u] = 0
                        orderMap.data[b, l:u] = o
                    else:
                        wlMap.data[l:u, b] = 0
                        orderMap.data[l:u, b] = o
            else:
                for b, u, l in zip(
                    axisBcoord, np.ceil(axisACoord_edgeup).astype(int), np.floor(axisACoord_edgelow).astype(int)
                ):
                    if self.axisA == "x":
                        wlMap.data[b, l:u] = np.nan
                        orderMap.data[b, l:u] = np.nan
                    else:
                        wlMap.data[l:u, b] = np.nan
                        orderMap.data[l:u, b] = np.nan

        # SLIT MAP PLACEHOLDER SAME AS WAVELENGTH MAP PLACEHOLDER
        slitMap = wlMap.copy()

        # PLOT CCDDATA OBJECT
        if self.debug and False:
            import matplotlib.pyplot as plt

            rotatedImg = np.rot90(slitMap.data, 1)
            rotatedImg = slitMap.data
            std = np.nanstd(slitMap.data)
            mean = np.nanmean(slitMap.data)
            vmax = mean + 3 * std
            vmin = mean - 3 * std
            plt.figure(figsize=(12, 5))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap="gray", alpha=1)
            plt.colorbar()
            plt.xlabel("y-axis", fontsize=10)
            plt.ylabel("x-axis", fontsize=10)
            plt.show()
            plt.close("all")

        self.log.debug("completed the ``create_placeholder_images`` method")
        return slitMap, wlMap, orderMap

    def map_to_image(self, dispersionMapPath, orders=False):
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
        self.log.debug("starting the ``map_to_image`` method")

        from soxspipe.commonutils.combiner import Combiner
        import numpy as np
        from astropy.io import fits
        import copy
        from fundamentals import fmultiprocess

        self.log.print("\n# CREATING 2D IMAGE MAP FROM DISPERSION SOLUTION\n\n")

        self.dispersionMapPath = dispersionMapPath
        kw = self.kw
        dp = self.detectorParams
        arm = self.arm

        self.map_to_image_displacement_threshold = self.recipeSettings["map_to_image_displacement_threshold"]
        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm
        )

        combinedSlitImage = False
        combinedWlImage = False

        if orders:
            theseOrders = orders
        else:
            theseOrders = orderNums

        # DEFINE AN INPUT ARRAY
        inputArray = [
            (order, minWl, maxWl)
            for order, minWl, maxWl in zip(orderNums, waveLengthMin, waveLengthMax)
            if order in theseOrders
        ]

        # NUMPY CAN BE TRICKY WITH MP
        numThreads = "1"
        os.environ["OPENBLAS_NUM_THREADS"] = numThreads
        os.environ["OMP_NUM_THREADS"] = numThreads
        os.environ["BLAS_NUM_THREADS"] = numThreads

        results = fmultiprocess(
            log=self.log,
            function=self.order_to_image,
            inputArray=inputArray,
            poolSize=6,
            timeout=3600,
            turnOffMP=self.debug,
        )
        del os.environ["OPENBLAS_NUM_THREADS"]
        del os.environ["OMP_NUM_THREADS"]
        del os.environ["BLAS_NUM_THREADS"]

        slitImages = [r[0] for r in results]
        wlImages = [r[1] for r in results]

        slitMap, wlMap, orderMap = self.create_placeholder_images(reverse=True)

        combinedSlitImage = Combiner(slitImages, dtype=np.float32)
        combinedSlitImage = combinedSlitImage.sum_combine()
        combinedWlImage = Combiner(wlImages, dtype=np.float32)
        combinedWlImage = combinedWlImage.sum_combine()

        combinedWlImage.data += wlMap.data
        combinedSlitImage.data += wlMap.data

        # GET THE EXTENSION (WITH DOT PREFIX)
        extension = os.path.splitext(dispersionMapPath)[1]
        filename = os.path.basename(dispersionMapPath).replace(extension, "_IMAGE.fits")

        dispersion_image_filePath = f"{self.productDir}/{filename}"
        # WRITE CCDDATA OBJECT TO FILE
        # from astropy.io import fits
        header = copy.deepcopy(self.dispMapHeader)
        header[kw("PRO_CATG")] = f"DISP_IMAGE_{arm}".upper()
        primary_hdu = fits.PrimaryHDU(combinedWlImage.data, header=header)
        primary_hdu.header["EXTNAME"] = "WAVELENGTH"
        image_hdu = fits.ImageHDU(combinedSlitImage.data)
        image_hdu.header["EXTNAME"] = "SLIT"
        image_hdu2 = fits.ImageHDU(orderMap.data)
        image_hdu2.header["EXTNAME"] = "ORDER"
        hdul = fits.HDUList([primary_hdu, image_hdu, image_hdu2])
        hdul.writeto(dispersion_image_filePath, output_verify="exception", overwrite=True, checksum=True)

        self.log.debug("completed the ``map_to_image`` method")
        return dispersion_image_filePath

    def order_to_image(self, orderInfo):
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
        self.log.debug("starting the ``order_to_image`` method")

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
        slitArray = np.arange(-halfGrid, halfGrid + grid_res_slit, grid_res_slit)
        wlArray = np.arange(minWl, maxWl, grid_res_wavelength)

        # ONE SINGLE-VALUE SLIT ARRAY FOR EVERY WAVELENGTH ARRAY
        bigSlitArray = np.concatenate([np.ones(wlArray.shape[0]) * slitArray[i] for i in range(0, slitArray.shape[0])])
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
                order=order,
                bigWlArray=bigWlArray,
                bigSlitArray=bigSlitArray,
                slitMap=slitMap,
                wlMap=wlMap,
                iteration=iteration,
                plots=False,
            )

            if remainingCount < 3:
                break

            orderPixelTable = orderPixelTable.drop_duplicates(subset=["pixel_x", "pixel_y", "order"])
            train_wlx = orderPixelTable["fit_x"].values
            train_wly = orderPixelTable["fit_y"].values
            train_wl = orderPixelTable["wavelength"].values
            train_sp = orderPixelTable["slit_position"].values

            targetX = orderPixelTable["pixel_x"].values
            targetY = orderPixelTable["pixel_y"].values

            if iteration == 1:
                # ADD MISSING PIXELS
                targetX = np.concatenate([targetX, allx])
                targetY = np.concatenate([targetY, ally])

            # USE CUBIC SPLINE NEAREST NEIGHBOUR TO SEED RESULTS
            bigWlArray = griddata((train_wlx, train_wly), train_wl, (targetX, targetY), method="cubic")
            bigSlitArray = griddata((train_wlx, train_wly), train_sp, (targetX, targetY), method="cubic")

        self.log.debug("completed the ``order_to_image`` method")
        return slitMap, wlMap

    def convert_and_fit(self, order, bigWlArray, bigSlitArray, slitMap, wlMap, iteration, plots=False):
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
        self.log.debug("starting the ``convert_and_fit`` method")

        import pandas as pd
        import numpy as np

        # CREATE PANDAS DATAFRAME WITH LARGE ARRAYS - ONE ROW PER
        # WAVELENGTH-SLIT GRID CELL
        myDict = {
            "order": np.ones(bigWlArray.shape[0]) * order,
            "wavelength": bigWlArray,
            "slit_position": bigSlitArray,
        }
        orderPixelTable = pd.DataFrame(myDict)

        # GET DETECTOR PIXEL POSITIONS FOR ALL WAVELENGTH-SLIT GRID CELLS
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log, dispersionMapPath=self.dispersionMapPath, orderPixelTable=orderPixelTable
        )

        # INTEGER PIXEL VALUES & FIT DISPLACEMENTS FROM PIXEL CENTRES
        orderPixelTable["pixel_x"] = np.round(orderPixelTable["fit_x"].values)
        orderPixelTable["pixel_y"] = np.round(orderPixelTable["fit_y"].values)
        orderPixelTable["residual_x"] = orderPixelTable["fit_x"] - orderPixelTable["pixel_x"]
        orderPixelTable["residual_y"] = orderPixelTable["fit_y"] - orderPixelTable["pixel_y"]
        orderPixelTable["residual_xy"] = np.sqrt(
            np.square(orderPixelTable["residual_x"]) + np.square(orderPixelTable["residual_y"])
        )

        # ADD A COUNT COLUMN FOR THE NUMBER OF SMALL SLIT/WL PIXELS FALLING IN
        # LARGE DETECTOR PIXELS
        count = orderPixelTable.groupby(["pixel_x", "pixel_y"]).size().reset_index(name="count")
        orderPixelTable = pd.merge(
            orderPixelTable, count, how="left", left_on=["pixel_x", "pixel_y"], right_on=["pixel_x", "pixel_y"]
        )
        orderPixelTable = orderPixelTable.sort_values(["order", "pixel_x", "pixel_y", "residual_xy"])

        # FILTER TO WL/SLIT POSITION CLOSE ENOUGH TO CENTRE OF PIXEL
        mask = orderPixelTable["residual_xy"] < self.map_to_image_displacement_threshold

        # KEEP ONLY VALUES CLOSEST TO CENTRE OF PIXEL
        newPixelValue = orderPixelTable.loc[mask].drop_duplicates(subset=["pixel_x", "pixel_y"], keep="first")

        # REMOVE PIXELS FOUND IN newPixelValue FROM orderPixelTable
        orderPixelTable = (
            newPixelValue[["pixel_x", "pixel_y"]]
            .merge(orderPixelTable, on=["pixel_x", "pixel_y"], how="right", indicator=True)
            .query('_merge == "right_only"')
            .drop(columns=["_merge"])
        )

        remainingCount = orderPixelTable.drop_duplicates(subset=["pixel_x", "pixel_y"], keep="first")

        # ADD FITTED PIXELS TO PLACE HOLDER IMAGES
        for xx, yy, wavelength, slit_position in zip(
            newPixelValue["pixel_x"].values.astype(int),
            newPixelValue["pixel_y"].values.astype(int),
            newPixelValue["wavelength"].values,
            newPixelValue["slit_position"].values,
        ):
            try:
                wlMap.data[yy, xx] = np.where(np.isnan(wlMap.data[yy, xx]), wavelength, wlMap.data[yy, xx])
                slitMap.data[yy, xx] = np.where(np.isnan(slitMap.data[yy, xx]), slit_position, slitMap.data[yy, xx])
            except IndexError:
                # PIXELS OUTSIDE OF DETECTOR EDGES - IGNORE
                pass

        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        percentageFound = (1 - (np.count_nonzero(np.isnan(wlMap.data)) / np.count_nonzero(wlMap.data))) * 100
        try:
            self.log.print(
                f"ORDER {order:02d}, iteration {iteration:02d}. {percentageFound:0.2f}% order pixels now fitted."
            )
        except:
            pass

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
            cmap.set_bad(color="#ADD8E6")
            vmax = np.nanmax(rotatedImg)
            vmin = np.nanmin(rotatedImg)
            plt.figure(figsize=(24, 10))
            plt.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap=cmap, alpha=1)
            if self.axisA == "x":
                plt.gca().invert_yaxis()
            plt.colorbar()
            plt.ylabel(f"{self.axisA}-axis", fontsize=12)
            plt.xlabel(f"{self.axisB}-axis", fontsize=12)
            plt.show()
            plt.close("all")

        remainingCount = len(remainingCount.index)

        self.log.debug("completed the ``convert_and_fit`` method")
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
        dispMapImage=False,
    ):
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
        self.log.debug("starting the ``create_dispersion_map_qc_plot`` method")

        import numpy as np
        from astropy.visualization import hist
        import matplotlib.pyplot as plt
        import pandas as pd
        from soxspipe.commonutils.toolkit import qc_settings_plot_tables
        from astropy.stats import sigma_clipped_stats

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
                slitPositions=self.uniqueSlitPos,
            )
        # DROP MISSING VALUES
        orderPixelTable.dropna(axis="index", how="any", subset=["residuals_x"], inplace=True)
        orderPixelTable["residuals_xy"] = np.sqrt(
            np.square(orderPixelTable["residuals_x"]) + np.square(orderPixelTable["residuals_y"])
        )
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
                aLen = rotatedImg.shape[0] - 1
                # OBSERVED VALUES
                orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]
                allClippedLines[f"observed_{self.axisA}"] = aLen - allClippedLines[f"observed_{self.axisA}"]
                # DETECTOR SHIFTED VALUES
                missingLines[f"detector_{self.axisA}_shifted"] = aLen - missingLines[f"detector_{self.axisA}_shifted"]
                allClippedLines[f"detector_{self.axisA}_shifted"] = (
                    aLen - allClippedLines[f"detector_{self.axisA}_shifted"]
                )
                orderPixelTable[f"detector_{self.axisA}_shifted"] = (
                    aLen - orderPixelTable[f"detector_{self.axisA}_shifted"]
                )
                # ORIGINAL DETECTOR VALUES
                allClippedLines[f"detector_{self.axisA}"] = aLen - allClippedLines[f"detector_{self.axisA}"]
                missingLines[f"detector_{self.axisA}"] = aLen - missingLines[f"detector_{self.axisA}"]
                orderPixelTable[f"detector_{self.axisA}"] = aLen - orderPixelTable[f"detector_{self.axisA}"]
                # FITTED VALUES
                orderPixelTable[f"fit_{self.axisA}"] = aLen - orderPixelTable[f"fit_{self.axisA}"]

                if not isinstance(gridLinePixelTable, bool):
                    gridLinePixelTable[f"fit_{self.axisA}"] = aLen - gridLinePixelTable[f"fit_{self.axisA}"]
                # orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]
                # orderPixelTable[f"observed_{self.axisA}"] = aLen - orderPixelTable[f"observed_{self.axisA}"]

        # a = plt.figure(figsize=(40, 15))

        if self.debug:
            fig = plt.figure(figsize=(6, 7), constrained_layout=True)
            gs = fig.add_gridspec(1, 1)
            toprow = fig.add_subplot(gs[:, :])

        elif rotatedImg.shape[0] / rotatedImg.shape[1] > 0.8:  # SOXS NIR
            fig = plt.figure(figsize=(6, 20), constrained_layout=True)
            # CREATE THE GRID OF AXES
            if self.arcFrame:
                gs = fig.add_gridspec(12, 4)
                sizeAx = fig.add_subplot(gs[8:9, :])
                gapAx = fig.add_subplot(gs[9:10, :])
                settingsAx = fig.add_subplot(gs[10:, 2:])
                qcAx = fig.add_subplot(gs[11:, 0:2])
            else:
                gs = fig.add_gridspec(10, 4)
                settingsAx = fig.add_subplot(gs[8:, 2:])
                qcAx = fig.add_subplot(gs[8:, 0:2])

            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:6, 0:2])
            bottomright = fig.add_subplot(gs[4:6, 2:])
            fwhmAx = fig.add_subplot(gs[6:7, :])
            resAx = fig.add_subplot(gs[7:8, :])

        else:  # SOXS VIS
            if self.firstGuessMap:
                fig = plt.figure(figsize=(6, 18), constrained_layout=True)
            else:
                fig = plt.figure(figsize=(6, 14), constrained_layout=True)
            # CREATE THE GRID OF AXES
            # CREATE THE GRID OF AXES
            if self.arcFrame:
                gs = fig.add_gridspec(12, 4)
                sizeAx = fig.add_subplot(gs[8:9, :])
                gapAx = fig.add_subplot(gs[9:10, :])
                settingsAx = fig.add_subplot(gs[10:, 2:])
                qcAx = fig.add_subplot(gs[11:, 0:2])
            else:
                gs = fig.add_gridspec(10, 4)
                settingsAx = fig.add_subplot(gs[8:, 2:])
                qcAx = fig.add_subplot(gs[8:, 0:2])

            toprow = fig.add_subplot(gs[0:2, :])
            midrow = fig.add_subplot(gs[2:4, :])
            bottomleft = fig.add_subplot(gs[4:6, 0:2])
            bottomright = fig.add_subplot(gs[4:6, 2:])
            fwhmAx = fig.add_subplot(gs[6:7, :])
            resAx = fig.add_subplot(gs[7:8, :])

        # TOP ROW - IMAGE WITH LINE POSITIONS OVERLAID
        vmax = self.meanFrameFlux + 25 * self.stdFrameFlux
        vmin = self.meanFrameFlux
        toprow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap="gray", alpha=0.5)
        toprow.set_title("observed arc-line positions (post-clipping)", fontsize=10)

        alphaBoost = 1.0
        if not self.firstGuessMap or self.debug:
            alphaBoost = 1.7

        if self.debug:
            # PLOT WHERE LINES GET SHIFTED TO
            toprow.scatter(
                missingLines[f"detector_{self.axisB}"],
                missingLines[f"detector_{self.axisA}"],
                marker="o",
                c="black",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="missed lines - original location",
            )
            toprow.scatter(
                allClippedLines[f"detector_{self.axisB}"],
                allClippedLines[f"detector_{self.axisA}"],
                marker="x",
                c="black",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="clipped lines - original location",
            )
            toprow.scatter(
                orderPixelTable[f"detector_{self.axisB}"],
                orderPixelTable[f"detector_{self.axisA}"],
                marker="o",
                c="black",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="detected lines - original location",
            )
            # PLOT WHERE LINES GET SHIFTED TO
            toprow.scatter(
                missingLines[f"detector_{self.axisB}_shifted"],
                missingLines[f"detector_{self.axisA}_shifted"],
                marker="o",
                c="red",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="missed lines - shifted location",
            )
            toprow.scatter(
                allClippedLines[f"detector_{self.axisB}_shifted"],
                allClippedLines[f"detector_{self.axisA}_shifted"],
                marker="x",
                c="red",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="clipped lines - shifted location",
            )
            toprow.scatter(
                orderPixelTable[f"detector_{self.axisB}_shifted"],
                orderPixelTable[f"detector_{self.axisA}_shifted"],
                marker="o",
                c="red",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="detected lines - shifted location",
            )
            # PLOT WHERE LINES GET SHIFTED TO
            toprow.scatter(
                allClippedLines[f"observed_{self.axisB}"],
                allClippedLines[f"observed_{self.axisA}"],
                marker="x",
                c="blue",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="clipped lines - observed location",
            )
            toprow.scatter(
                orderPixelTable[f"observed_{self.axisB}"],
                orderPixelTable[f"observed_{self.axisA}"],
                marker="o",
                c="blue",
                s=10,
                alpha=0.2 * alphaBoost,
                linewidths=0.5,
                label="detected lines - observed location",
            )

        else:

            if isinstance(missingLines, pd.core.frame.DataFrame):
                toprow.scatter(
                    missingLines[f"detector_{self.axisB}_shifted"],
                    missingLines[f"detector_{self.axisA}_shifted"],
                    marker="o",
                    c="black",
                    s=7,
                    alpha=0.1 * alphaBoost,
                    linewidths=0.5,
                    label="undetected line location",
                )

                # toprow.scatter(orderPixelTable[f"detector_{self.axisB}_shifted"], orderPixelTable[f"detector_{self.axisA}_shifted"], marker='o', c='black', s=20, alpha=0.4 * alphaBoost, linewidths=0.5, label="undetected line location")
                # toprow.scatter(allClippedLines[f"detector_{self.axisB}_shifted"], allClippedLines[f"detector_{self.axisA}_shifted"], marker='o', c='black', s=20, alpha=0.4 * alphaBoost, linewidths=0.5, label="undetected line location")
            if len(allClippedLines.index):
                pass
                mask = allClippedLines["dropped"] == True
                if self.firstGuessMap:
                    toprow.scatter(
                        allClippedLines.loc[mask][f"observed_{self.axisB}"],
                        allClippedLines.loc[mask][f"observed_{self.axisA}"],
                        marker="o",
                        c="blue",
                        s=5,
                        alpha=0.3 * alphaBoost,
                        linewidths=0.5,
                        label="dropped multi-pinhole set",
                    )
                else:
                    toprow.scatter(
                        allClippedLines.loc[mask][f"observed_{self.axisB}"],
                        allClippedLines.loc[mask][f"observed_{self.axisA}"],
                        marker="o",
                        c="blue",
                        s=5,
                        alpha=0.3 * alphaBoost,
                        linewidths=0.5,
                        label="dropped pinhole",
                    )

                toprow.scatter(
                    allClippedLines.loc[~mask][f"observed_{self.axisB}"],
                    allClippedLines.loc[~mask][f"observed_{self.axisA}"],
                    marker="o",
                    c="green",
                    s=5,
                    alpha=0.3 * alphaBoost,
                    linewidths=0.5 * alphaBoost,
                )

                toprow.scatter(
                    allClippedLines.loc[~mask][f"observed_{self.axisB}"],
                    allClippedLines.loc[~mask][f"observed_{self.axisA}"],
                    marker="x",
                    c="red",
                    s=5,
                    alpha=0.3 * alphaBoost,
                    linewidths=0.5,
                    label="clipped during dispersion solution fitting",
                )

            if len(orderPixelTable.index):
                toprow.scatter(
                    orderPixelTable[f"observed_{self.axisB}"],
                    orderPixelTable[f"observed_{self.axisA}"],
                    marker="o",
                    c="green",
                    s=5,
                    alpha=0.3 * alphaBoost,
                    linewidths=0.5,
                    label="detected line location",
                )
                # toprow.scatter(orderPixelTable[f"tilt_corrected_{self.axisB}"], orderPixelTable[f"tilt_corrected_{self.axisA}"], marker='o', c='purple', s=8, alpha=0.3 * alphaBoost, linewidths=0.5 * alphaBoost, )

        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        toprow.tick_params(axis="both", which="major", labelsize=9)
        if self.arm == "VIS":
            toprow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.2), fontsize=4)
        else:
            toprow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.15), fontsize=4)

        toprow.set_xlim([0, rotatedImg.shape[1]])
        if self.axisA == "x":
            toprow.invert_yaxis()
        toprow.set_ylim([0, rotatedImg.shape[0]])

        if self.debug:
            plt.show()

        midrow.imshow(rotatedImg, vmin=vmin, vmax=vmax, cmap="gray", alpha=0.5)
        midrow.set_title("global dispersion solution", fontsize=10)

        # ADD FULL DISPERSION SOLUTION GRID-LINES TO PLOT
        if not isinstance(gridLinePixelTable, bool):
            for l in range(int(gridLinePixelTable["line"].max())):
                mask = gridLinePixelTable["line"] == l
                if l == 1:
                    midrow.plot(
                        gridLinePixelTable.loc[mask][f"fit_{self.axisB}"],
                        gridLinePixelTable.loc[mask][f"fit_{self.axisA}"],
                        "w-",
                        linewidth=0.2,
                        alpha=0.5 * alphaBoost,
                        color="blue",
                        label="dispersion solution",
                    )
                    toprow.plot(
                        gridLinePixelTable.loc[mask][f"fit_{self.axisB}"],
                        gridLinePixelTable.loc[mask][f"fit_{self.axisA}"],
                        "w-",
                        linewidth=0.2,
                        alpha=0.5 * alphaBoost,
                        color="blue",
                        label="dispersion solution",
                    )
                else:
                    midrow.plot(
                        gridLinePixelTable.loc[mask][f"fit_{self.axisB}"],
                        gridLinePixelTable.loc[mask][f"fit_{self.axisA}"],
                        "w-",
                        linewidth=0.2,
                        alpha=0.5 * alphaBoost,
                        color="blue",
                    )
                    toprow.plot(
                        gridLinePixelTable.loc[mask][f"fit_{self.axisB}"],
                        gridLinePixelTable.loc[mask][f"fit_{self.axisA}"],
                        "w-",
                        linewidth=0.2,
                        alpha=0.5 * alphaBoost,
                        color="blue",
                    )
        else:
            midrow.scatter(
                orderPixelTable[f"fit_{self.axisB}"],
                orderPixelTable[f"fit_{self.axisA}"],
                marker="o",
                c="blue",
                s=orderPixelTable[f"residuals_xy"] * 30,
                alpha=0.1 * alphaBoost,
                label="fitted line (size proportional to line-fit residual)",
            )

        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis="both", which="major", labelsize=9)
        midrow.set_xlim([0, rotatedImg.shape[1]])
        if self.axisA == "x":
            midrow.invert_yaxis()
        midrow.set_ylim([0, rotatedImg.shape[0]])

        if self.arm == "VIS":
            midrow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.2), fontsize=4)
        else:
            midrow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.15), fontsize=4)

        # PLOT THE RESIDUALS
        plt.subplots_adjust(top=0.92)
        bottomleft.scatter(
            orderPixelTable[f"residuals_{self.axisA}"], orderPixelTable[f"residuals_{self.axisB}"], alpha=0.1
        )
        bottomleft.set_xlabel(f"{self.axisA} residual (mean: {mean_x_res:2.2f} pix)")
        bottomleft.set_ylabel(f"{self.axisB} residual (mean: {mean_y_res:2.2f} pix)")
        bottomleft.tick_params(axis="both", which="major", labelsize=9)
        hist(
            orderPixelTable["residuals_xy"],
            bins="scott",
            ax=bottomright,
            histtype="stepfilled",
            alpha=0.7,
            density=True,
        )
        bottomright.set_xlabel("xy residual")
        bottomright.tick_params(axis="both", which="major", labelsize=9)

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        if self.firstGuessMap:
            fig.suptitle(
                f"residuals of global dispersion solution fitting - {arm} multi-pinhole\n{subtitle}",
                fontsize=10,
                y=0.99,
            )
        else:
            fig.suptitle(
                f"residuals of global dispersion solution fitting - {arm} single pinhole\n{subtitle}",
                fontsize=10,
                y=0.99,
            )
        orderPixelTable_groups = orderPixelTable.groupby(["order"])

        if self.arcFrame:
            mean_fwhm, median, std_fwhm = sigma_clipped_stats(
                orderPixelTable["fwhm_slit_px"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
            )
            mean_r, median, std_r = sigma_clipped_stats(
                orderPixelTable["R_slit"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
            )
        else:
            mean_fwhm, median, std_fwhm = sigma_clipped_stats(
                orderPixelTable["fwhm_pin_px"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
            )
            mean_r, median, std_r = sigma_clipped_stats(
                orderPixelTable["R_pin"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
            )

        lowerFwhm = mean_fwhm - 3 * std_fwhm
        upperFwhm = mean_fwhm + 3 * std_fwhm
        lowerR = mean_r - 3 * std_r
        upperR = mean_r + 3 * std_r

        if self.arcFrame:

            # UNPACK ORDER TABLE TO GET SLIT HEIGHT IN PX
            orderPolyTable, orderGeoTable, orderMetaTable = unpack_order_table(
                log=self.log, orderTablePath=self.orderTable, extend=0.0
            )

            orderGeoTable[f"{self.axisA}_u"] = orderGeoTable[f"{self.axisA}coord_edgeup"].astype(int)
            orderGeoTable[f"{self.axisA}_l"] = orderGeoTable[f"{self.axisA}coord_edgelow"].astype(int)

            dispMapDF, interOrderMask = twoD_disp_map_image_to_dataframe(
                log=self.log, slit_length=14, twoDMapPath=dispMapImage, associatedFrame=self.arcFrame, kw=self.kw
            )

            # MERGE DATAFRAMES
            orderGeoTableLower = orderGeoTable.merge(
                dispMapDF, left_on=[f"{self.axisA}_l", f"{self.axisB}coord"], right_on=[self.axisA, self.axisB]
            )
            orderGeoTableLower.rename(columns={"order_x": "order"}, inplace=True)
            orderGeoTableLower = orderGeoTableLower[
                [
                    "order",
                    "x",
                    "y",
                    "wavelength",
                    "slit_position",
                    "flux",
                    "pixelScale",
                    f"{self.axisA}coord_centre",
                    f"{self.axisA}coord_edgelow",
                ]
            ]

            orderGeoTableUpper = orderGeoTable.merge(
                dispMapDF, left_on=[f"{self.axisA}_u", f"{self.axisB}coord"], right_on=[self.axisA, self.axisB]
            )
            orderGeoTableUpper.rename(columns={"order_x": "order"}, inplace=True)
            orderGeoTableUpper = orderGeoTableUpper[
                [
                    "order",
                    "x",
                    "y",
                    "wavelength",
                    "slit_position",
                    "flux",
                    "pixelScale",
                    f"{self.axisA}coord_centre",
                    f"{self.axisA}coord_edgeup",
                ]
            ]

            orderGeoTable = orderGeoTableLower.merge(
                orderGeoTableUpper,
                left_on=[f"{self.axisA}coord_centre", self.axisB],
                right_on=[f"{self.axisA}coord_centre", self.axisB],
                suffixes=("_l", "_u"),
            )
            orderGeoTable.rename(columns={"order_l": "order"}, inplace=True)

            orderGeoTable["wavelength"] = (orderGeoTable["wavelength_l"] + orderGeoTable["wavelength_u"]) / 2.0
            orderGeoTable["slitLengthPixelsInt"] = np.abs(
                orderGeoTable[f"{self.axisA}_u"] - orderGeoTable[f"{self.axisA}_l"]
            )
            orderGeoTable["slitLengthPixels"] = np.abs(
                orderGeoTable[f"{self.axisA}coord_edgeup"] - orderGeoTable[f"{self.axisA}coord_edgelow"]
            )
            orderGeoTable["slitLengthArcsec"] = np.abs(
                orderGeoTable[f"slit_position_u"] - orderGeoTable[f"slit_position_l"]
            )
            orderGeoTable["pixelScale"] = orderGeoTable["slitLengthArcsec"] / orderGeoTable["slitLengthPixelsInt"]
            orderGeoTable["slitLengthArcsec"] = orderGeoTable["slitLengthPixels"] * orderGeoTable["pixelScale"]

        for name, group in orderPixelTable_groups:

            thisOrder = group["order"].values[0]

            if self.arcFrame:

                slitWidth = self.arcFrame.header[kw(f"SLIT_{arm}")].replace("SLIT", "").split("x")[0]
                fwhmAx.set_title(f'Line FWHM measured via {slitWidth}" slit arc-lamp frame', fontsize=9)
                fwhmAx.scatter(group["wavelength"], group["fwhm_slit_px"], alpha=0.1)
                mean_fwhm, median, std_fwhm = sigma_clipped_stats(
                    group["fwhm_slit_px"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
                )
                fwhmAx.set_ylabel("FWHM (px)", fontsize=9)
                mean_wavelength = group["wavelength"].mean()
                fwhmAx.errorbar(mean_wavelength, mean_fwhm, yerr=std_fwhm, fmt="o", color="black", alpha=0.6)
                fwhmAx.tick_params(axis="both", which="major", labelsize=8)
                x_limits = fwhmAx.get_xlim()

                resAx.set_title(f'Resolution measured via {slitWidth}" slit arc-lamp frame', fontsize=9)
                resAx.scatter(group["wavelength"], group["R_slit"], alpha=0.1)
                mean_resol, median, std_resol = sigma_clipped_stats(
                    group["R_slit"].values, sigma=3.0, stdfunc="std", cenfunc="mean", maxiters=3
                )
                # ADD TO THE POINT THE ERROR BAR CONTAINED IN STD_RED
                resAx.errorbar(mean_wavelength, mean_resol, yerr=std_resol, fmt="o", color="black", alpha=0.6)

                if thisOrder not in [24]:

                    # FILTER DATA FRAME
                    # FIRST CREATE THE MASK
                    mask = orderGeoTable["order"] == thisOrder
                    filteredDf = orderGeoTable.loc[mask]

                    sizeAx.plot(filteredDf["wavelength"], filteredDf["slitLengthArcsec"])
                    sizeAx.set_xlabel("wavelength (nm)", fontsize=9)
                    sizeAx.set_ylabel("slit height\n(arcsec)", fontsize=9)
                    sizeAx.set_xlim(x_limits)
                    sizeAx.tick_params(axis="both", which="major", labelsize=8)
                    sizeAx.set_title(f"Slit height as measured between the lower and upper order edges", fontsize=9)

                    # PLOTTING THE ORDER GAP DISTANCES
                    if thisOrder not in [23]:
                        mask = orderGeoTable["order"] == thisOrder
                        order_l = orderGeoTable[mask]
                        mask = orderGeoTable["order"] == thisOrder + 1
                        order_u = orderGeoTable[mask]
                        if len(order_l.index) and len(order_u.index):
                            merged_ul = order_l.merge(order_u, on=[self.axisB], suffixes=("_l", "_u"))
                            merged_ul["order_distance"] = (
                                merged_ul[f"{self.axisA}coord_edgelow_u"] - merged_ul[f"{self.axisA}coord_edgeup_l"]
                            )
                            pixel_scale_mean = (merged_ul["pixelScale_l"] + merged_ul["pixelScale_u"]) / (2.0)
                            gapAx.plot(
                                merged_ul["wavelength_l"],
                                merged_ul["order_distance"] * pixel_scale_mean,
                                label=thisOrder,
                            )
                            gapAx.set_xlabel("wavelength (nm)", fontsize=9)
                            gapAx.set_ylabel("inter-order\ngap (arcsec)", fontsize=9)

                        gapAx.set_xlim(x_limits)
                        gapAx.tick_params(axis="both", which="major", labelsize=8)
                        gapAx.set_title(f"Inter-order gap measured between adjacent orders", fontsize=9)

            else:
                fwhmAx.set_title('Line FWHM measured via 0.5" pinhole  arc-lamp frame', fontsize=9)
                fwhmAx.scatter(group["wavelength"], group["fwhm_pin_px"], alpha=0.1)
                mean_fwhm = group["fwhm_pin_px"].mean()
                std_fwhm = group["fwhm_pin_px"].std()
                fwhmAx.set_ylabel("FWHM (px)", fontsize=9)
                mean_wavelength = group["wavelength"].mean()
                fwhmAx.errorbar(mean_wavelength, mean_fwhm, yerr=std_fwhm, fmt="o", color="black", alpha=0.6)
                fwhmAx.tick_params(axis="both", which="major", labelsize=8)
                x_limits = fwhmAx.get_xlim()

                resAx.set_title('Resolution as measured via 0.5" pinhole arc-lamp frame', fontsize=9)
                resAx.scatter(group["wavelength"], group["R_pin"], alpha=0.1)
                # CALCULATE THE MEAN AND STD DEV OF THE GROUP AND ADD TO THE PLOT
                mean_resol = group["R_pin"].mean()
                std_resol = group["R_pin"].std()
                # ADD TO THE POINT THE ERROR BAR CONTAINED IN STD_RED
                resAx.errorbar(mean_wavelength, mean_resol, yerr=std_resol, fmt="o", color="black", alpha=0.6)

        resAx.set_ylabel("R", fontsize=9)
        resAx.tick_params(axis="both", which="major", labelsize=8)

        # resAx.scatter(orderPixelTable["wavelength"], orderPixelTable["R"], alpha=0.1)
        resAx.set_xlabel("wavelength (nm)", fontsize=9)

        fwhmAx.set_ylim([lowerFwhm, upperFwhm])
        resAx.set_ylim([lowerR, upperR])
        resAx.set_xlim(x_limits)

        utcnow = datetime.now(timezone.utc)
        utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")

        # GET FILENAME FOR THE RESIDUAL PLOT
        if not self.sofName:
            res_plots = filenamer(log=self.log, frame=self.pinholeFrame, settings=self.settings)
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
        self.products = pd.concat(
            [
                self.products,
                pd.Series(
                    {
                        "soxspipe_recipe": self.recipeName,
                        "product_label": "DISP_MAP_RES",
                        "file_name": res_plots,
                        "file_type": "PDF",
                        "obs_date_utc": self.dateObs,
                        "reduction_date_utc": utcnow,
                        "product_desc": f"{self.arm} dispersion solution QC plots",
                        "file_path": filePath,
                        "label": "QC",
                    }
                )
                .to_frame()
                .T,
            ],
            ignore_index=True,
        )

        qc_settings_plot_tables(
            log=self.log,
            qc=self.qc,
            qcAx=qcAx,
            settings={**self.recipeSettings, **{"exptime": self.exptime}},
            settingsAx=settingsAx,
        )

        plt.tight_layout()
        if self.debug:
            plt.show()

        if not self.settings["tune-pipeline"]:
            plt.savefig(filePath, dpi=240)

        plt.close("all")

        if self.settings["tune-pipeline"]:
            import codecs

            filePath = f"residuals.txt"
            exists = os.path.exists(filePath)
            if not exists:
                with codecs.open(filePath, encoding="utf-8", mode="w") as writeFile:
                    writeFile.write(
                        f"polyOrders,mean_x_res,mean_y_res,mean_res,std_res,median_res,median_x_res,median_y_res,CLINE \n"
                    )
            with codecs.open(filePath, encoding="utf-8", mode="a") as writeFile:
                writeFile.write(
                    f"{polyOrders},{mean_x_res:0.4f},{mean_y_res:0.4f},{mean_res:2.4f},{std_res:2.4f},{median_res:2.4f},{median_x_res:2.4f},{median_y_res:2.4f},{self.CLINE}\n"
                )

        self.log.debug("completed the ``create_dispersion_map_qc_plot`` method")
        return res_plots

    def _clip_on_measured_line_metrics(self, orderPixelTable):
        """*clip lines & sets of lines based on measured line metrics (from daostarfinder etc)*

        **Key Arguments:**

        - `orderPixelTable` -- panda's data-frame containing measure line metrics

        **Return:**

        - `orderPixelTable` -- the data-frame with clipped lines indicated in the "sigma_clipped" column
        """
        self.log.debug("starting the ``_clip_on_measured_line_metrics`` method")

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

        toprow.set_title("Pre-fitting QC Clipping on Pinhole-Sets", fontsize=10)

        columnsNoStrings = list(orderPixelTable.columns)
        try:
            columnsNoStrings.remove("ion")
        except:
            pass

        if self.firstGuessMap:

            # GROUP BY ARC LINES (MPH SETS)

            lineGroups = (
                orderPixelTable.loc[(orderPixelTable["dropped"] == False)][columnsNoStrings]
                .groupby(["wavelength", "order"])
                .std()
            )
            lineGroups = lineGroups.reset_index()

            # SIGMA-CLIP THE DATA ON SCATTER
            # masked_residuals = sigma_clip(
            #     lineGroups["xy_diff"], sigma_lower=5000, sigma_upper=7, maxiters=5, cenfunc='median', stdfunc='mad_std')
            # lineGroups["sigma_clipped_scatter"] = masked_residuals.mask

            # SIGMA-CLIP THE DATA ON SCATTER
            lineGroups["sigma_clipped_scatter"] = False
            masked_residuals = sigma_clip(
                lineGroups["x_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc="median", stdfunc="mad_std"
            )
            lineGroups["sigma_clipped_x"] = masked_residuals.mask
            masked_residuals = sigma_clip(
                lineGroups["y_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc="median", stdfunc="mad_std"
            )
            lineGroups["sigma_clipped_y"] = masked_residuals.mask
            masked_residuals = sigma_clip(
                lineGroups["xy_diff"], sigma_lower=5000, sigma_upper=7, maxiters=3, cenfunc="median", stdfunc="mad_std"
            )
            lineGroups["sigma_clipped_xy"] = masked_residuals.mask
            lineGroups.loc[
                (
                    (lineGroups["sigma_clipped_y"] == True)
                    | (lineGroups["sigma_clipped_x"] == True)
                    | (lineGroups["sigma_clipped_xy"] == True)
                ),
                "sigma_clipped_scatter",
            ] = True

            bottomleft.scatter(
                x=lineGroups["x_diff"],  # numpy array of x-points
                y=lineGroups["y_diff"],  # numpy array of y-points
                # 1 number or array of areas for each datapoint (i.e. point size)
                s=1,
                c="black",  # color or sequence of color, optional, default
                marker="x",
                alpha=0.6,
                label="arc line shift from predicted position",
            )

            bottomleft.scatter(
                # numpy array of x-points
                x=lineGroups.loc[(lineGroups["sigma_clipped_scatter"] == True)]["x_diff"],
                # numpy array of y-points
                y=lineGroups.loc[(lineGroups["sigma_clipped_scatter"] == True)]["y_diff"],
                # 1 number or array of areas for each datapoint (i.e. point size)
                s=5,
                c="red",  # color or sequence of color, optional, default
                marker="x",
                alpha=0.6,
                label="clipped arc lines",
            )

            bottomleft.set_ylabel(f"y-shift rms (px)", fontsize=12)
            bottomleft.set_xlabel(f"x-shift rms (px)", fontsize=12)
            bottomleft.tick_params(axis="both", which="major", labelsize=9)
            bottomleft.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

            hist(
                lineGroups[(lineGroups["sigma_clipped_scatter"] == False)]["xy_diff"],
                bins="scott",
                ax=bottomright,
                histtype="stepfilled",
                alpha=0.7,
                density=True,
            )
            bottomright.set_xlabel("xy residual")
            bottomright.tick_params(axis="both", which="major", labelsize=9)

            # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FWHM
            mask = lineGroups["sigma_clipped_scatter"] == True
            dropGroups = lineGroups.loc[mask]
            setsToDrop = dropGroups[["wavelength", "order"]]
            s = orderPixelTable[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
            s["dropped"] = False
            s.loc[(s["_merge"] == "both"), "dropped"] = True
            orderPixelTable["droppedOnScatter"] = s["dropped"].values
            orderPixelTable.loc[(orderPixelTable["droppedOnScatter"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FWHM
        lineGroups = (
            orderPixelTable.loc[(orderPixelTable["dropped"] == False)][columnsNoStrings]
            .groupby(["wavelength", "order"])
            .mean()
        )
        lineGroups = lineGroups.reset_index()
        masked_residuals = sigma_clip(
            lineGroups["fwhm_pin_px"], sigma_lower=2.5, sigma_upper=5, maxiters=3, cenfunc="median", stdfunc="mad_std"
        )
        lineGroups["sigma_clipped_fwhm"] = masked_residuals.mask
        lineGroups["sigma_clipped"] = masked_residuals.mask

        toprow.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["fwhm_pin_px"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=1,
            c="black",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="measured pinhole lines",
        )

        toprow.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_fwhm"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_fwhm"] == True)]["fwhm_pin_px"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=5,
            c="red",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="clipped pinhole lines",
        )

        toprow.set_ylabel(f"fwhm (px)", fontsize=12)
        toprow.set_xlabel(f"wavelength (nm)", fontsize=12)
        toprow.tick_params(axis="both", which="major", labelsize=9)
        toprow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON FLUX
        mask = lineGroups["sigma_clipped_fwhm"] == True
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[["wavelength", "order"]]
        s = orderPixelTable[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnFWHM"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnFWHM"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FLUX
        lineGroups = lineGroups.loc[~mask]
        masked_residuals = sigma_clip(
            lineGroups["flux"], sigma_lower=5, sigma_upper=5, maxiters=5, cenfunc="mean", stdfunc="std"
        )
        lineGroups["sigma_clipped_flux"] = masked_residuals.mask

        midrow.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["flux"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=1,
            c="black",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="pinhole flux",
        )

        midrow.set_ylabel(f"pinhole flux", fontsize=12)
        midrow.set_xlabel(f"wavelength (nm)", fontsize=12)
        midrow.tick_params(axis="both", which="major", labelsize=9)
        midrow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

        midrow.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_flux"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_flux"] == True)]["flux"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=5,
            c="red",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="clipped pinhole lines",
        )

        midrow.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA BEFORE CLIPPING ON PEAK
        mask = lineGroups["sigma_clipped_flux"] == True
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[["wavelength", "order"]]
        s = orderPixelTable[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnFlux"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnFlux"] == True), "dropped"] = True

        # SIGMA-CLIP THE DATA ON FLUX
        lineGroups = lineGroups.loc[~mask]
        lineGroups["peak"] = lineGroups["peak"] / lineGroups["flux"]
        masked_residuals = sigma_clip(
            lineGroups["peak"], sigma_lower=5000, sigma_upper=7, maxiters=5, cenfunc="mean", stdfunc="std"
        )
        lineGroups["sigma_clipped_peak"] = masked_residuals.mask

        midrow2.scatter(
            x=lineGroups["wavelength"],  # numpy array of x-points
            y=lineGroups["peak"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=1,
            c="black",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="pinhole peak flux / mean flux",
        )

        midrow2.set_ylabel(f"pinhole peak flux", fontsize=12)
        midrow2.set_xlabel(f"wavelength (nm)", fontsize=12)
        midrow2.tick_params(axis="both", which="major", labelsize=9)
        midrow2.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

        midrow2.scatter(
            x=lineGroups[(lineGroups["sigma_clipped_peak"] == True)]["wavelength"],  # numpy array of x-points
            y=lineGroups[(lineGroups["sigma_clipped_peak"] == True)]["peak"],  # numpy array of y-points
            # 1 number or array of areas for each datapoint (i.e. point size)
            s=5,
            c="red",  # color or sequence of color, optional, default
            marker="x",
            alpha=0.6,
            label="clipped pinhole lines",
        )

        midrow2.legend(loc="upper right", bbox_to_anchor=(1.0, -0.05), fontsize=4)

        # REMOVE THE CLIPPED DATA
        mask = lineGroups["sigma_clipped_peak"] == True
        dropGroups = lineGroups.loc[mask]
        setsToDrop = dropGroups[["wavelength", "order"]]
        s = orderPixelTable[["wavelength", "order"]].merge(setsToDrop, indicator=True, how="left")
        s["dropped"] = False
        s.loc[(s["_merge"] == "both"), "dropped"] = True
        orderPixelTable["droppedOnPeak"] = s["dropped"].values
        orderPixelTable.loc[(orderPixelTable["droppedOnPeak"] == True), "dropped"] = True

        if self.debug:
            plt.show()
        plt.close("all")

        self.log.debug("completed the ``_clip_on_measured_line_metrics`` method")
        return orderPixelTable

    def update_static_line_list_detector_positions(self, originalOrderPixelTable, dispersionMapPath):
        """*using a first pass dispersion solution, update the original static line list*

        **Key Arguments:**

        - `originalOrderPixelTable` -- original order pixel table
        - `dispersionMapPath` -- path to the first pass dispersion solution

        **Return:**

        - `updatedLineList` -- updated static line list
        """
        self.log.debug("starting the ``update_static_line_list_detector_positions`` method")

        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        from soxspipe.commonutils.toolkit import read_spectral_format
        import pandas as pd
        from astropy.table import Table

        # GET UNIQUE VALUES OF order AND WAVELENGTH
        uniquecolNames = originalOrderPixelTable[["order", "wavelength"]].drop_duplicates()
        orders = uniquecolNames["order"]
        wavelengths = uniquecolNames["wavelength"]
        slit_positions = self.uniqueSlitPos
        slit_indexes = list(range(0, len(self.uniqueSlitPos), 1))

        dfCollection = []
        for si, sp in zip(slit_indexes, slit_positions):
            myDict = {
                "order": orders,
                "wavelength": wavelengths,
                "slit_position": [sp] * len(orders),
                "slit_index": [si] * len(orders),
            }

            orderPixelTable = pd.DataFrame(myDict)
            floats = ["wavelength", "slit_position"]
            ints = ["order", "slit_index"]
            for f in floats:
                orderPixelTable[f] = orderPixelTable[f].astype(float)
            for i in ints:
                orderPixelTable[i] = orderPixelTable[i].astype(int)

            orderPixelTable = dispersion_map_to_pixel_arrays(
                log=self.log,
                dispersionMapPath=dispersionMapPath,
                orderPixelTable=orderPixelTable,
                removeOffDetectorLocation=False,
            )

            orderPixelTable.rename(columns={"fit_x": "detector_x", "fit_y": "detector_y"}, inplace=True)
            orderPixelTable = orderPixelTable[
                ["wavelength", "order", "slit_index", "slit_position", "detector_x", "detector_y"]
            ]

            dfCollection.append(orderPixelTable)

        updatedLineList = pd.concat(dfCollection, ignore_index=True)

        self.log.debug("completed the ``update_static_line_list_detector_positions`` method")
        return updatedLineList

    def create_new_static_line_list(self, dispersionMapPath):
        """*using a first pass dispersion solution, use a line atlas to generate a more accurate and more complete static line list*

        **Key Arguments:**

        - `dispersionMapPath` -- path to the first pass dispersion solution

        **Return:**

        - `newPredictedLineList` -- a new predicted line list (to replace the static calibration line-list)
        """
        self.log.debug("starting the ``create_new_static_line_list`` method")

        from soxspipe.commonutils import dispersion_map_to_pixel_arrays
        from soxspipe.commonutils.toolkit import read_spectral_format
        import pandas as pd
        from astropy.table import Table

        dp = self.detectorParams

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax = read_spectral_format(
            log=self.log, settings=self.settings, arm=self.arm
        )

        # FIND THE LINE ATLAS
        calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)
        lineAtlas = calibrationRootPath + "/" + dp["line-atlas"]
        # LINE LIST TO PANDAS DATAFRAME
        lineAtlas = Table.read(lineAtlas, format="fits")
        lineAtlas = lineAtlas.to_pandas()

        # CLEAN UP LINE LIST DATA
        lineAtlas["ion"] = lineAtlas["ion"].str.decode("ascii")
        lineAtlas["source"] = lineAtlas["source"].str.decode("ascii")

        lineAtlas["wave"] = lineAtlas["wave"].round(5)

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK

        # GET UNIQUE VALUES IN COLUMN
        uniqueions = lineAtlas["ion"].unique()

        # # mask = (lineAtlas['ion'].isin(["XeI", 'UNK']))
        # mask = (lineAtlas['ion'].isin(["HgI", 'UNK']))
        # # mask = (lineAtlas['ion'].isin(["NeI", "NeII", 'UNK']))
        # # mask = (lineAtlas['ion'].isin(["Ar", "ArI", "ArII", "ArIII", 'UNK']))
        # lineAtlas = lineAtlas.loc[mask]

        # ORDER BY MOST INTENSE LINES
        lineAtlas.sort_values(["amplitude"], ascending=[False], inplace=True)
        lineAtlas = lineAtlas.head(500)

        # try:
        #     # FILTER DATAFRAME TO EXCLUDE source = "NIST" and amplitude > 1000
        #     mask = (lineAtlas['source'] != "NIST") & (lineAtlas['amplitude'] > 1000)
        #     lineAtlas = lineAtlas.loc[mask]
        # except:
        #     pass
        dfCollection = []
        for o, wmin, wmax in zip(orderNums, waveLengthMin, waveLengthMax):

            wrange = wmax - wmin
            wmin -= wrange / 5
            wmax += wrange / 5

            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            mask = lineAtlas["wave"].between(wmin, wmax)
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
                    "ion": ion,
                }

                orderPixelTable = pd.DataFrame(myDict)
                floats = ["wavelength", "slit_position"]
                ints = ["order", "slit_index"]
                for f in floats:
                    orderPixelTable[f] = orderPixelTable[f].astype(float)
                for i in ints:
                    orderPixelTable[i] = orderPixelTable[i].astype(int)

                orderPixelTable = dispersion_map_to_pixel_arrays(
                    log=self.log,
                    dispersionMapPath=dispersionMapPath,
                    orderPixelTable=orderPixelTable,
                    removeOffDetectorLocation=False,
                )

                orderPixelTable.rename(columns={"fit_x": "detector_x", "fit_y": "detector_y"}, inplace=True)
                orderPixelTable = orderPixelTable[
                    ["ion", "wavelength", "order", "slit_index", "slit_position", "detector_x", "detector_y"]
                ]

                dfCollection.append(orderPixelTable)

        newPredictedLineList = pd.concat(dfCollection, ignore_index=True)

        from astropy.table import Table

        t = Table.from_pandas(newPredictedLineList)
        # t.write("Xe.fits", overwrite=True)
        t.write("Hg.fits", overwrite=True)
        # t.write("Ne.fits", overwrite=True)
        # t.write("Ar.fits", overwrite=True)

        # sys.exit(0)

        self.log.debug("completed the ``create_new_static_line_list`` method")
        return newPredictedLineList

    # use the tab-trigger below for new method
    # xt-class-method


def measure_line_position(
    stampInfo,
    log,
    windowHalf,
    iraf,
    sigmaLimit,
    iteration,
    brightest=False,
    exclude_border=False,
    multipinhole=False,
    returnAll=True,
    debug=False,
):
    """*measure the line position on a given stamp*

    **Key Arguments:**

    - `stampInfo` -- stamp pixels, xlow, xup, ylow, yup
    - `log` -- logger
    - `windowHalf` -- half of the stamp side
    - `iraf` -- run IRAF source detection to get FWHM?
    - `sigmaLimit` -- stamp source detection minimum sigma
    - `brightest` -- find the brightest source?
    - `multipinhole` -- is this a multipinhole frame?

    **Returns:**

    - list of detected line positions with properties
    """
    log.debug("starting the ``measure_line_position`` function")

    import numpy as np
    from photutils import DAOStarFinder, IRAFStarFinder
    from astropy.stats import sigma_clipped_stats
    import logging

    # FIX ASTROPY LOGGING LEVEL RESET
    logging.getLogger().setLevel(logging.INFO + 5)

    # EXTRACT STAMP INFORMATION
    stamp, xlow, xup, ylow, yup = stampInfo[0], stampInfo[1], stampInfo[2], stampInfo[3], stampInfo[4]

    # CONVERT CCDDATA TO MASKED NUMPY ARRAY
    stamp = np.ma.array(stamp.data, mask=stamp.mask)

    # CALCULATE STAMP STATISTICS FOR DETECTION THRESHOLD
    mean, median, std = sigma_clipped_stats(stamp, sigma=3.0, stdfunc="mad_std", cenfunc="median")

    # PREPARE ITERATION TEXT FOR LOGGING
    iterationText = f", iter #{iteration}" if iteration is not False else ""

    # DETECT SOURCES USING DAO STARFINDER
    try:
        daofind = DAOStarFinder(
            fwhm=2.0,
            threshold=sigmaLimit * std,
            roundlo=-2.0,
            roundhi=2.0,
            sharplo=-0.5,
            sharphi=2.0,
            exclude_border=exclude_border,
        )
        # SUBTRACT MEDIAN FOR BETTER DETECTION IN LOW SIGNAL IMAGES
        sources = daofind(stamp.data - median, mask=stamp.mask)
    except Exception as e:
        sources = None

    # INITIALIZE DETECTION VARIABLES
    old_resid = windowHalf * 4
    selectedSource = None
    observed_x = np.nan
    observed_y = np.nan
    fwhm = np.nan
    old_detectionSigma = 0.0

    predictedLines = []
    keepValues = ["sharpness", "roundness1", "roundness2", "npix", "sky", "peak", "flux"]

    # PROCESS DETECTED SOURCES
    if sources:

        if returnAll:
            # RETURN ALL DETECTED SOURCES
            for source in sources:
                predictedLine = {
                    "observed_x": source["xcentroid"] + xlow,
                    "observed_y": source["ycentroid"] + ylow,
                    "stamp_x": source["xcentroid"],
                    "stamp_y": source["ycentroid"],
                    "detectionSigma": source["peak"] / std,
                    "fwhm_pin_px": np.nan,
                }
                # ADD PHOTOMETRIC PROPERTIES
                for k in keepValues:
                    predictedLine[k] = source[k]
                predictedLines.append(predictedLine)
        else:
            # SELECT BEST SINGLE SOURCE (CLOSEST TO CENTER OR BRIGHTEST)
            if len(sources) > 1:
                for source in sources:
                    tmp_x = source["xcentroid"]
                    tmp_y = source["ycentroid"]
                    new_resid = np.sqrt((windowHalf - tmp_x) ** 2 + (windowHalf - tmp_y) ** 2)
                    detectionSigma = source["peak"] / std

                    # SELECT BY BRIGHTNESS OR PROXIMITY TO CENTER
                    if (brightest and detectionSigma > old_detectionSigma) or (not brightest and new_resid < old_resid):
                        observed_x = tmp_x + xlow
                        observed_y = tmp_y + ylow
                        stamp_x = tmp_x
                        stamp_y = tmp_y
                        old_resid = new_resid
                        selectedSource = source
                        old_detectionSigma = detectionSigma
            else:
                # ONLY ONE SOURCE DETECTED
                observed_x = sources[0]["xcentroid"] + xlow
                observed_y = sources[0]["ycentroid"] + ylow
                stamp_x = sources[0]["xcentroid"]
                stamp_y = sources[0]["ycentroid"]
                selectedSource = sources[0]
                detectionSigma = sources[0]["peak"] / std

            # BUILD PREDICTED LINE DICTIONARY
            predictedLine = {
                "observed_x": observed_x,
                "observed_y": observed_y,
                "stamp_x": stamp_x,
                "stamp_y": stamp_y,
                "detectionSigma": detectionSigma,
            }

            # ADD PHOTOMETRIC PROPERTIES
            for k in keepValues:
                predictedLine[k] = selectedSource[k]

            # OPTIONALLY MEASURE FWHM USING IRAF STARFINDER
            if iraf:
                try:
                    iraf_find = IRAFStarFinder(
                        fwhm=3.0,
                        threshold=1 * std,
                        roundlo=-2.0,
                        roundhi=2.0,
                        sharplo=-1,
                        sharphi=3.0,
                        exclude_border=True,
                        xycoords=[(stamp_x, stamp_y)],
                    )
                    iraf_sources = iraf_find(stamp)
                    predictedLine["fwhm_pin_px"] = iraf_sources["fwhm"][0]
                except Exception:
                    predictedLine["fwhm_pin_px"] = np.nan
            else:
                predictedLine["fwhm_pin_px"] = np.nan
            predictedLines.append(predictedLine)

        # DEBUG VISUALIZATION (DISABLED)
        if False and debug and iteration == 0:
            import matplotlib.pyplot as plt

            plt.clf()
            plt.imshow(stamp)
            plt.scatter(0, 0, marker="x", s=30)
            plt.scatter(observed_x - xlow, observed_y - ylow, s=30)
            plt.text(
                windowHalf - 2,
                windowHalf - 2,
                f"{observed_x-xlow:0.2f},{observed_y - ylow:0.2f}",
                fontsize=16,
                c="black",
                verticalalignment="bottom",
            )
            plt.text(
                2,
                2,
                f"{sigmaLimit}$\\sigma$, {iterationText}\n,{detectionSigma:0.1f}",
                fontsize=16,
                c="white",
                verticalalignment="bottom",
            )
            plt.show()

    log.debug("completed the ``measure_line_position`` function")
    return predictedLines


def straighten_mph_sets(group):
    """
    *Straighten multi-pinhole sets by projecting points onto fitted line*

    **Key Arguments:**

    - ``group`` -- dataframe group containing observed_x and observed_y coordinates

    **Returns:**

    - ``group`` -- group with added tilt_corrected_x and tilt_corrected_y columns
    """
    import numpy as np

    # EXTRACT OBSERVED COORDINATES
    x = group["observed_x"].values
    y = group["observed_y"].values

    # FIT LINEAR MODEL: y = m*x + c
    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]

    # PROJECT EACH POINT ONTO THE FITTED LINE
    # For point (x0, y0), find closest point (xp, yp) on line y = m*x + c
    # Perpendicular projection formulas:
    #   xp = (m*(y0 - c) + x0) / (m**2 + 1)
    #   yp = m*xp + c
    xp = (m * (y - c) + x) / (m**2 + 1)
    yp = m * xp + c

    # ADD TILT-CORRECTED COORDINATES TO GROUP
    group = group.copy()
    group["tilt_corrected_x"] = xp
    group["tilt_corrected_y"] = yp

    # DEBUG VISUALIZATION (DISABLED)
    if False:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 4))
        plt.scatter(x, y, color="blue", label="Original points")
        plt.plot(x, m * x + c, color="green", label="Fitted line")
        plt.scatter(xp, yp, color="red", marker="x", label="Tilt-corrected points")

        if "detector_x_shifted" in group.columns and "detector_y_shifted" in group.columns:
            plt.scatter(
                group["detector_x_shifted"],
                group["detector_y_shifted"],
                color="orange",
                marker="^",
                label="Detector shifted",
            )

        # LABEL EACH POINT WITH SLIT INDEX
        for xi, yi, slit_idx in zip(x, y, group["slit_index"].values):
            plt.text(xi, yi, str(slit_idx), fontsize=8, color="black", ha="center", va="bottom")

        plt.xlabel("observed_x")
        plt.ylabel("observed_y")
        plt.legend()
        plt.title("Straighten Multi-Pinhole Set")
        plt.tight_layout()
        plt.text(
            0.95, 0.95, f"N={len(x)}", transform=plt.gca().transAxes, fontsize=10, verticalalignment="top", color="blue"
        )
        plt.show()

    return group


def _plot_slit_index_comparisons(df):
    """
    Plot slit_index vs xy_diff, x_diff, and y_diff, color-coded by slit_index.
    Replace legend labels with mean and std for each slit_index.
    The order number is shown as the figure title.
    Each panel title includes the global mean and std for that metric.
    Adds a fourth panel: scatter plot of x_diff vs y_diff, color-coded by slit_index.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    slit_indexes = np.sort(df["slit_index"].unique())
    try:
        order_num = df["order"].iloc[0] if "order" in df.columns else None
    except:
        return
    cmap = plt.get_cmap("tab10")
    colors = {idx: cmap(i % 10) for i, idx in enumerate(slit_indexes)}

    # Compute global mean and std for each metric
    global_xy_mean, global_xy_std = df["xy_diff"].mean(), df["xy_diff"].std()
    global_x_mean, global_x_std = df["x_diff"].mean(), df["x_diff"].std()
    global_y_mean, global_y_std = df["y_diff"].mean(), df["y_diff"].std()

    fig, axes = plt.subplots(1, 4, figsize=(24, 5), sharey=False)
    handles = [[], [], [], []]
    labels = [[], [], [], []]
    for idx in slit_indexes:
        mask = df["slit_index"] == idx
        xy_mean, xy_std = df.loc[mask, "xy_diff"].mean(), df.loc[mask, "xy_diff"].std()
        x_mean, x_std = df.loc[mask, "x_diff"].mean(), df.loc[mask, "x_diff"].std()
        y_mean, y_std = df.loc[mask, "y_diff"].mean(), df.loc[mask, "y_diff"].std()
        h0 = axes[0].scatter([idx] * np.sum(mask), df.loc[mask, "xy_diff"], color=colors[idx], alpha=0.7)
        h1 = axes[1].scatter([idx] * np.sum(mask), df.loc[mask, "x_diff"], color=colors[idx], alpha=0.7)
        h2 = axes[2].scatter([idx] * np.sum(mask), df.loc[mask, "y_diff"], color=colors[idx], alpha=0.7)
        h3 = axes[3].scatter(
            df.loc[mask, "x_diff"], df.loc[mask, "y_diff"], color=colors[idx], alpha=0.7, label=f"slit_index {idx}"
        )
        handles[0].append(h0)
        handles[1].append(h1)
        handles[2].append(h2)
        handles[3].append(h3)
        labels[0].append(f"slit_index {idx}: μ={xy_mean:.3f}, σ={xy_std:.3f}")
        labels[1].append(f"slit_index {idx}: μ={x_mean:.3f}, σ={x_std:.3f}")
        labels[2].append(f"slit_index {idx}: μ={y_mean:.3f}, σ={y_std:.3f}")
        labels[3].append(f"slit_index {idx}")

    axes[0].set_title(f"slit_index vs xy_diff\nGlobal μ={global_xy_mean:.3f}, σ={global_xy_std:.3f}")
    axes[1].set_title(f"slit_index vs x_diff\nGlobal μ={global_x_mean:.3f}, σ={global_x_std:.3f}")
    axes[2].set_title(f"slit_index vs y_diff\nGlobal μ={global_y_mean:.3f}, σ={global_y_std:.3f}")
    axes[3].set_title("x_diff vs y_diff (color: slit_index)")
    axes[3].set_xlabel("x_diff")
    axes[3].set_ylabel("y_diff")
    for i, ax in enumerate(axes):
        if i < 3:
            ax.set_xlabel("slit_index")
            ax.set_ylabel("difference")
            ax.legend(handles[i], labels[i], title="slit_index", bbox_to_anchor=(1.05, 1), loc="upper left")
        else:
            ax.legend(handles[i], labels[i], title="slit_index", bbox_to_anchor=(1.05, 1), loc="upper left")
    if order_num is not None:
        fig.suptitle(f"Order {order_num}. {len(df.index)} points", fontsize=16)
    plt.tight_layout()
    plt.show()
    return


def find_largest_cluster_center(x, y, eps=2.0, min_samples=5):
    """
    *Find the center of the largest cluster in (x, y) data using DBSCAN*

    **Key Arguments:**

    - ``x`` -- x-coordinates of points
    - ``y`` -- y-coordinates of points
    - ``eps`` -- maximum distance between points in a cluster
    - ``min_samples`` -- minimum samples required to form a cluster

    **Returns:**

    - ``center_x, center_y`` -- coordinates of largest cluster center, or (None, None) if no cluster found
    """
    import numpy as np
    from sklearn.cluster import DBSCAN

    # STACK COORDINATES AND RUN DBSCAN CLUSTERING
    points = np.column_stack((x, y))
    clustering = DBSCAN(eps=eps, min_samples=min_samples).fit(points)
    labels = clustering.labels_

    # FIND LARGEST CLUSTER (EXCLUDING NOISE LABEL -1)
    unique, counts = np.unique(labels[labels != -1], return_counts=True)

    if len(counts) == 0:
        # NO CLUSTERS FOUND - RETURN NONE
        return None, None

    # CALCULATE CENTER OF LARGEST CLUSTER
    largest_cluster = unique[np.argmax(counts)]
    mask = labels == largest_cluster
    center_x = np.mean(x[mask])
    center_y = np.mean(y[mask])

    return center_x, center_y
