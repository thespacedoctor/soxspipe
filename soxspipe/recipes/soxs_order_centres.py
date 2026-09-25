#!/usr/bin/env python
"""
*further constrain the first guess locations of the order centres derived in `soxs_disp_solution`*

Author
: David Young & Marco Landoni

Date Created
: September  8, 2020
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils import detect_continuum
from soxspipe.commonutils.toolkit import append_product, utcnow_string

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_order_centres(base_recipe):
    """
    *further constrain the first guess locations of the order centres derived in `soxs_disp_solution`*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``polyOrders`` -- the orders of the x-y polynomials used to fit the dispersion solution. Overrides parameters found in the yaml settings file. e.g 345400 is order_x=3, order_y=4 ,wavelength_x=5 ,wavelength_y=4. Default *False*.
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- debug mode. True or False. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.


    **Usage**

    ```python
    from soxspipe.recipes import soxs_order_centres
    productPath, qcTable = soxs_order_centres(
        log=log,
        settings=settings,
        inputFrames=a["inputFrames"]
    ).produce_product()
    ```
    """

    # INITIALISATION

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=[],
        verbose=False,
        overwrite=False,
        polyOrders=False,
        command=False,
        debug=False,
        turnOffMP=False,
    ):
        # INHERIT INITIALISATION FROM  base_recipe
        super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-order-centres",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_order_centres' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.polyOrders = polyOrders

        self._parse_poly_orders()
        self._collect_input_frames()
        self._verify_and_announce_input_frames()
        self._sort_and_report_input_frames()

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def _parse_poly_orders(self):
        """*validate the two-digit `polyOrders` override*

        Sets ``self.polyOrders``. A false value is left alone, so the recipe
        falls back to the degrees in the settings file.
        """
        if self.polyOrders is not False:
            value = self.polyOrders
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, str))
                or len(str(value)) != 2
                or not str(value).isascii()
                or not str(value).isdigit()
            ):
                raise TypeError("THE poly VALUE NEEDS TO BE A 2 DIGIT INTEGER")
            self.polyOrders = int(value)

    def _collect_input_frames(self):
        """*resolve the recipe's input into a ccdproc image collection*

        Sets ``self.inputFrames`` and ``self.supplementaryInput``.
        """
        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            ext=self.settings["data-extension"],
        )
        self.inputFrames, self.supplementaryInput = sof.get()

    def _verify_and_announce_input_frames(self):
        """*verify the collected frames and report the outcome to the terminal*

        Sets ``self.imageType``, through ``verify_input_frames``.
        """
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_order_centres - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

    def _sort_and_report_input_frames(self):
        """*sort the input frames by observation date and print the verbose summary*

        Sets no attribute; sorts ``self.inputFrames`` in place.
        """
        # SORT IMAGE COLLECTION
        self.inputFrames.sort(["MJD-OBS"])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary)
            self.log.print("\n")

    def verify_input_frames(self):
        """*verify input frames match those required by the soxs_order_centres recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            error = self._nir_input_frame_error(imageTypes, imageTech, imageCat, self.arm)
        else:
            error = self._uvb_vis_input_frame_error(imageTypes, imageCat, self.arm)

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]
        self.log.debug("completed the ``verify_input_frames`` method")
        return

    def _nir_input_frame_error(self, imageTypes, imageTech, imageCat, arm):
        """*report why a NIR frame set is not a single-pinhole flat-lamp set with its dispersion table*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageTech`` -- the image techniques the basic verification classified
        - ``imageCat`` -- the product categories the basic verification classified
        - ``arm`` -- the arm under reduction, which names the dispersion table required

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        # WANT ON AND OFF PINHOLE FRAMES
        # MIXED INPUT IMAGE TYPES ARE BAD
        if not error and len(imageTypes) > 1:
            joinedImageTypes = " and ".join(imageTypes)
            error = f"Input frames are a mix of {joinedImageTypes}"

        if not error:
            good = "FLAT" if self.inst == "SOXS" else "LAMP,ORDERDEF"
            if good not in imageTypes[0]:
                error = (
                    "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off "
                    "frames and a first-guess dispersion solution table for NIR"
                )

        if not error:
            for i in imageTech:
                if i not in ["ECHELLE,PINHOLE", "IMAGE"]:
                    error = (
                        "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off "
                        "frames a first-guess dispersion solution table for NIR"
                    )

        if not error:
            for i in [f"DISP_TAB_{arm}"]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe order_centres need to be single pinhole flat-lamp on and lamp off "
                        "frames a first-guess dispersion solution table for NIR"
                    )

        return error

    def _uvb_vis_input_frame_error(self, imageTypes, imageCat, arm):
        """*report why a UVB or VIS frame set is not a single-pinhole flat-lamp set with its calibrations*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageCat`` -- the product categories the basic verification classified
        - ``arm`` -- the arm under reduction, which names the master bias and dispersion table required

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        if not error:
            if self.inst == "SOXS":
                goodList = ["FLAT,LAMP", "LAMP,DFLAT", "LAMP,FLAT"]
            else:
                goodList = ["LAMP,ORDERDEF", "LAMP,DORDERDEF", "LAMP,QORDERDEF"]
            for i in imageTypes:
                if i not in goodList:
                    error = (
                        "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias "
                        "frame, a first-guess dispersion solution table and possibly a master dark for UVB/VIS. "
                        f"Found {i}"
                    )

        if not error:
            for i in [f"MASTER_BIAS_{arm}", f"DISP_TAB_{arm}"]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe order_centres need to be single pinhole flat-lamp, a master-bias "
                        "frame, a first-guess dispersion solution table and possibly a master dark for UVB/VIS."
                    )

        return error

    def produce_product(self):
        """*generate the order-table with polynomal fits of order-centres*

        **Return:**

        - ``productPath`` -- the path to the order-table
        - ``qcTable`` -- the reported quality-control table

        When the pipeline is tuning rather than reducing, the method returns a single
        ``None`` instead of the pair.

        **Usage:**

        ```python
        productPath, qcTable = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        arm = self.arm
        kw = self.kw

        productPath = None

        master_bias, dark = self._read_calibration_frames(kw, arm)
        orderDef_image = self._read_order_definition_frame(kw)

        add_filters = {kw("PRO_CATG"): f"DISP_TAB_{arm}".upper()}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            disp_map_table = i

        self._calibrate_order_frame(orderDef_image, master_bias, dark)
        binx, biny = self._read_trace_binning(kw, arm)

        if self.settings["tune-pipeline"]:
            from itertools import product

            digits = [2, 3, 4, 5, 6]
            perm = product(digits, repeat=2)
            try:
                os.remove("residuals.txt")
            except OSError as e:
                self.log.debug(f"produce_product: `os.remove('residuals.txt')` failed, continuing: {e}")

            self._tune_order_centre_parameters(perm, disp_map_table, binx, biny)
            return None

        if self.polyOrders:
            self.polyOrders = str(self.polyOrders)
            self.polyOrders = [int(digit) for digit in str(self.polyOrders)]
            self.recipeSettings["detect-continuum"]["order-deg"] = self.polyOrders[0]
            self.recipeSettings["detect-continuum"]["disp-axis-deg"] = self.polyOrders[1]

        productPath = self._fit_order_centres(kw, disp_map_table, binx, biny)

        qcTable = self.report_output()
        self.clean_up()

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _read_calibration_frames(self, kw, arm):
        """*read the master bias and the dark frame this reduction detrends with*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller
        - ``arm`` -- the arm under reduction, read by the caller

        **Return:**

        - ``master_bias`` -- the master bias frame, or False when the set carries none
        - ``dark`` -- the dark frame, or False when the set carries none. The lamp-off
          frame is read after the master dark, so it overrides one

        **Usage:**

        ```python
        master_bias, dark = self._read_calibration_frames(self.kw, self.arm)
        ```
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        master_bias = False
        dark = False

        add_filters = {kw("PRO_CATG"): "MASTER_BIAS_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_bias = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        # UVB/VIS DARK
        add_filters = {kw("PRO_CATG"): "MASTER_DARK_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        # NIR DARK
        if self.inst == "SOXS":
            add_filters = {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "IMAGE"}
        else:
            add_filters = {kw("DPR_TYPE"): "LAMP,ORDERDEF", kw("DPR_TECH"): "IMAGE"}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            dark = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        return master_bias, dark

    def _read_order_definition_frame(self, kw):
        """*read the single-pinhole flat-lamp frame whose order centres are traced*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller

        **Return:**

        - ``orderDef_image`` -- the last frame matching the instrument's filter list, or
          False when the set carries none

        **Usage:**

        ```python
        orderDef_image = self._read_order_definition_frame(self.kw)
        ```
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        orderDef_image = False

        if self.inst == "SOXS":
            filter_list = [
                {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                {kw("DPR_TYPE"): "LAMP,FLAT", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                {kw("DPR_TYPE"): "LAMP,DFLAT", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                # KEYWORD SCREW-UP DURING PAE MEANT WE HAD TO ADD BELOW WITH ECHELLE,SLIT ...
                # SHOULD REMOVE THIS EVENTUALLY
                {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "ECHELLE,SLIT"},
                {kw("DPR_TYPE"): "LAMP,FLAT", kw("DPR_TECH"): "ECHELLE,SLIT"},
                {kw("DPR_TYPE"): "LAMP,DFLAT", kw("DPR_TECH"): "ECHELLE,SLIT"},
            ]
        else:
            # UVB XSHOOTER - CHECK FOR D2 LAMP FIRST AND IF NOT FOUND USE THE QTH LAMP
            filter_list = [
                {kw("DPR_TYPE"): "LAMP,ORDERDEF", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                {kw("DPR_TYPE"): "LAMP,QORDERDEF", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                {kw("DPR_TYPE"): "LAMP,DORDERDEF", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
            ]

        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                orderDef_image = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        return orderDef_image

    def _calibrate_order_frame(self, orderDef_image, master_bias, dark):
        """*detrend the order-definition frame, stamp its keywords and optionally save it*

        **Key Arguments:**

        - ``orderDef_image`` -- the order-definition frame to calibrate
        - ``master_bias`` -- the master bias frame, or False
        - ``dark`` -- the dark frame, or False

        Sets ``self.orderFrame``.
        """
        self.orderFrame = self.detrend(inputFrame=orderDef_image, master_bias=master_bias, dark=dark)

        self.update_fits_keywords(frame=self.orderFrame)

        if self.settings["save-intermediate-products"]:
            fileDir = self.workspaceRootPath
            filepath = self._write(self.orderFrame, fileDir, filename=False, overwrite=True, product=False)
            self.log.print(f"\nCalibrated single pinhole frame frame saved to {filepath}\n")

    def _read_trace_binning(self, kw, arm):
        """*read the calibrated frame's binning, which selects the predicted line-list*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller
        - ``arm`` -- the arm under reduction, read by the caller

        **Return:**

        - ``binx`` -- the x binning; 1 on the NIR arm or when the header carries none
        - ``biny`` -- the y binning; 1 on the NIR arm or when the header carries none
        """
        # FIND THE APPROPRIATE PREDICTED LINE-LIST
        if arm != "NIR" and kw("WIN_BINX") in self.orderFrame.header:
            binx = int(self.orderFrame.header[kw("WIN_BINX")])
            biny = int(self.orderFrame.header[kw("WIN_BINY")])
        else:
            binx = 1
            biny = 1

        return binx, biny

    def _tune_order_centre_parameters(self, perm, disp_map_table, binx, biny):
        """*sample the trace once, then sweep the continuum polynomial degrees over the grid*

        **Key Arguments:**

        - ``perm`` -- the iterator of (order degree, dispersion-axis degree) pairs to try
        - ``disp_map_table`` -- the path to the first-guess dispersion table
        - ``binx`` -- the x binning of the calibrated frame
        - ``biny`` -- the y binning of the calibrated frame
        """
        # DETECT THE CONTINUUM OF ORDERE CENTRES - RETURN ORDER TABLE FILE PATH
        # self.log.print("\n# DETECTING ORDER CENTRE CONTINUUM\n")
        detector = detect_continuum(
            log=self.log,
            traceFrame=self.orderFrame,
            dispersion_map=disp_map_table,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName="soxs-order-centres",
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            binx=binx,
            biny=biny,
            startNightDate=self.startNightDate,
        )
        orderPixelTable, detectionPercentage = detector.sample_trace()

        print("\n\nTUNING SOXSPIPE\n")

        # DEFINE AN INPUT ARRAY
        from fundamentals import fmultiprocess

        # NOTE TO SELF: IF HAVING ISSUE WITH MULTIPROCESSING STALLING, TRY AND IMPORT REQUIRED MODULES INTO THE
        # METHOD/FUNCTION RUNNING THIS fmultiprocess FUNCTION INSTEAD OF AT THE MODULE LEVEL
        fmultiprocess(
            log=self.log,
            function=parameterTuning,
            inputArray=list(perm),
            poolSize=False,
            timeout=360000,
            recipeSettings=self.recipeSettings,
            settings=self.settings,
            orderFrame=self.orderFrame,
            disp_map_table=disp_map_table,
            orderPixelTable=orderPixelTable,
            qc=self.qc,
            products=self.products,
            sofName=self.sofName,
            binx=binx,
            biny=biny,
            turnOffMP=self.debug,
            mute=True,
            progressBar=True,
        )

    def _fit_order_centres(self, kw, disp_map_table, binx, biny):
        """*fit the order centre traces and record the order table as a product*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller
        - ``disp_map_table`` -- the path to the first-guess dispersion table
        - ``binx`` -- the x binning of the calibrated frame
        - ``biny`` -- the y binning of the calibrated frame

        **Return:**

        - ``productPath`` -- the path to the order table

        Sets ``self.products``, ``self.qc`` and ``self.dateObs``.

        A continuum fit that does not converge raises an ``ArithmeticError``, and no
        order table is recorded. The detector's quality-control rows are merged before
        the failure, so they survive it.
        """
        import pandas as pd

        # DETECT THE CONTINUUM OF ORDERE CENTRES - RETURN ORDER TABLE FILE PATH
        # self.log.print("\n# DETECTING ORDER CENTRE CONTINUUM\n")
        detector = detect_continuum(
            log=self.log,
            traceFrame=self.orderFrame,
            dispersion_map=disp_map_table,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            recipeName="soxs-order-centres",
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            binx=binx,
            biny=biny,
            startNightDate=self.startNightDate,
        )
        (
            productPath,
            qcTable,
            productsTable,
            orderPolyTable,
            orderPixelTable,
            orderMetaTable,
        ) = detector.get()

        self.products = pd.concat([self.products, productsTable])
        self.qc = pd.concat([self.qc, qcTable])

        # THE DETECTOR REPORTS A FAILED FIT AS A MISSING PRODUCT PATH. THE QUALITY-CONTROL
        # ROWS MERGE FIRST, SO THE EVIDENCE OF THE FAILED ATTEMPT SURVIVES THE FAILURE.
        if productPath is None:
            raise ArithmeticError(
                f"Could not converge on a good fit to the {self.arm} order-centre continuum. "
                "Please check the quality of your data or adjust your fitting parameters. "
                "No order table was produced."
            )

        filename = os.path.basename(productPath)

        utcnow = utcnow_string()
        self.dateObs = self.orderFrame.header[kw("DATE_OBS")]

        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel="ORDER_CENTRES",
            fileName=filename,
            filePath=productPath,
            productDesc=f"{self.arm} order centre traces",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        return productPath


# THE NAME IS PUBLIC: EACH RECIPE PASSES IT TO fmultiprocess AND tests/unit/test_parameter_tuning.py
# CALLS IT BY NAME, AS IN soxs_disp_solution AND soxs_spatial_solution
def parameterTuning(  # noqa: N802
    p,
    log,
    recipeSettings,
    settings,
    orderFrame,
    disp_map_table,
    orderPixelTable,
    qc,
    products,
    sofName,
    binx,
    biny,
):
    """*tuning the spatial solution*

    **Key Arguments:**

    - ``p`` -- one permutation of the order and dispersion-axis polynomial degrees
    - ``log`` -- logger
    - ``recipeSettings`` -- the recipe settings dictionary, rewritten in place with ``p``
    - ``settings`` -- the settings dictionary
    - ``orderFrame`` -- the calibrated order-definition frame to trace
    - ``disp_map_table`` -- the path to the first-guess dispersion table
    - ``orderPixelTable`` -- the sampled trace to reuse across permutations
    - ``qc`` -- the quality-control table to pass to the continuum detector
    - ``products`` -- the products table to pass to the continuum detector
    - ``sofName`` -- the name of the set-of-files this reduction came from
    - ``binx`` -- the x binning of the calibrated frame
    - ``biny`` -- the y binning of the calibrated frame

    The fit's own outputs are discarded.

    **Usage:**

    ```python
    parameterTuning(
        (3, 5),
        log=log,
        recipeSettings=recipeSettings,
        settings=settings,
        orderFrame=orderFrame,
        disp_map_table=disp_map_table,
        orderPixelTable=orderPixelTable,
        qc=qc,
        products=products,
        sofName=sofName,
        binx=binx,
        biny=biny,
    )
    ```
    """

    recipeSettings["detect-continuum"]["order-deg"] = p[0]
    recipeSettings["detect-continuum"]["disp-axis-deg"] = p[1]

    detector = detect_continuum(
        log=log,
        traceFrame=orderFrame,
        dispersion_map=disp_map_table,
        settings=settings,
        recipeSettings=recipeSettings,
        recipeName="soxs-order-centres",
        qcTable=qc,
        productsTable=products,
        sofName=sofName,
        binx=binx,
        biny=biny,
        orderPixelTable=orderPixelTable,
        startNightDate=self.startNightDate,
    )
    try:
        (
            productPath,
            qcTable,
            productsTable,
            orderPolyTable,
            orderPixelTable,
            orderMetaTable,
        ) = detector.get()
    except Exception as e:
        log.warning(f"parameterTuning: this tuning iteration failed and records nothing in the grid, continuing: {e}")

    return

    # USE THE TAB-TRIGGER BELOW FOR NEW METHOD
    # xt-class-method

    # OVERRIDE METHOD ATTRIBUTES
    # method-override-tmpx
