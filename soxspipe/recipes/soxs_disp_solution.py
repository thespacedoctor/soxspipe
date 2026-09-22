#!/usr/bin/env python
"""
*Recipe to generate a first approximation of the dispersion solution from single pinhole frames*

Author
: David Young & Marco Landoni

Date Created
: August 25, 2020
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils import create_dispersion_map
from soxspipe.commonutils.toolkit import append_product, utcnow_string

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_disp_solution(base_recipe):
    """
    *generate a first approximation of the dispersion solution from single pinhole frames*

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
    from soxspipe.recipes import soxs_disp_solution
    disp_map_path = soxs_disp_solution(
        log=log,
        settings=settings,
        inputFrames=sofPath
    ).produce_product()
    ```
    """

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
            recipeName="soxs-disp-solution",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_disp_solution' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.recipeName = "soxs-disp-solution"
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
        """*coerce the `polyOrders` override to an integer, or reject it*

        Sets ``self.polyOrders``. A false value is left alone, so the recipe
        falls back to the degrees in the settings file.
        """
        if self.polyOrders:
            try:
                self.polyOrders = int(self.polyOrders)
            except (ValueError, TypeError) as e:
                self.log.debug(f"__init__: `self.polyOrders = int(self.polyOrders)` failed, continuing: {e}")
            if not isinstance(self.polyOrders, int):
                raise TypeError("THE poly VALUE NEEDS TO BE A 4 DIGIT INTEGER")

    def _collect_input_frames(self):
        """*resolve the recipe's input into a ccdproc image collection*

        Sets ``self.inputFrames`` and ``self.supplementaryInput``.
        """
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
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_disp_solution - NO MORE, NO LESS.
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
        """*verify input frames match those required by the `soxs_disp_solution` recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if False:
            print(imageTypes)
            print(imageTech)
            print(imageCat)
            print(imageTypes[0])

        if self.arm == "NIR":
            error = self._nir_input_frame_error(imageTypes, imageTech)
        else:
            error = self._uvb_vis_input_frame_error(imageTypes, imageTech, imageCat, self.arm)

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

    def _nir_input_frame_error(self, imageTypes, imageTech):
        """*report why a NIR frame set is not a single-pinhole lamp-on and lamp-off pair*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageTech`` -- the image techniques the basic verification classified

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        # WANT ON AND OFF PINHOLE FRAMES
        # MIXED INPUT IMAGE TYPES ARE BAD
        if not error and len(imageTypes) > 1:
            imageTypes = " and ".join(imageTypes)
            imageTypes = " and ".join(imageTypes)
            error = (
                "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                "and lamp off frames for NIR"
            )

        # FIX ME!
        if not error and imageTypes[0] not in ["LAMP,FMTCHK", "LAMP,WAVE", "WAVE,LAMP"]:
            error = (
                "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                "and lamp off frames for NIR"
            )

        if not error:
            for i in imageTech:
                if i not in ["ECHELLE,PINHOLE", "IMAGE"]:
                    error = (
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                        "and lamp off frames for NIR"
                    )

        if not error:
            for i in ["ECHELLE,PINHOLE", "IMAGE"]:
                if i not in imageTech:
                    error = (
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                        "and lamp off frames for NIR"
                    )

        return error

    def _uvb_vis_input_frame_error(self, imageTypes, imageTech, imageCat, arm):
        """*report why a UVB or VIS frame set is not a single-pinhole set with its master bias*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageTech`` -- the image techniques the basic verification classified
        - ``imageCat`` -- the product categories the basic verification classified
        - ``arm`` -- the arm under reduction, which names the master bias required

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        if not error:
            for i in imageTypes:
                # FIX ME!
                if i not in ["LAMP,FMTCHK", "LAMP,WAVE", "WAVE,LAMP"]:
                    error = (
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                        "and a master-bias and possibly a master dark for UVB/VIS"
                    )

        if not error:
            for i in ["ECHELLE,PINHOLE"]:
                if i not in imageTech:
                    error = (
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                        "and a master-bias and possibly a master dark for UVB/VIS"
                    )

        if not error:
            for i in [f"MASTER_BIAS_{arm}"]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe disp_solution need to be single pinhole lamp on "
                        "and a master-bias and possibly a master dark for UVB/VIS"
                    )

        return error

    def produce_product(self):
        """*generate a fisrt guess of the dispersion solution*

        **Return:**

        - ``productPath`` -- the path to the first guess dispersion map, or None when
          the pipeline is tuning rather than reducing
        - ``qcTable`` -- the reported quality-control table

        **Usage:**

        ```python
        productPath, qcTable = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        arm = self.arm
        kw = self.kw

        # self.inputFrames.summary.pprint_all()

        master_bias, dark = self._read_calibration_frames(kw, arm)
        pinhole_image = self._read_pinhole_frame(kw)

        self._calibrate_pinhole_frame(master_bias, dark, pinhole_image)

        if self.settings["tune-pipeline"]:
            productPath, qcTable = self._tune_dispersion_parameters()
        else:
            productPath, qcTable = self._fit_dispersion_solution(kw)

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _read_calibration_frames(self, kw, arm):
        """*read the master bias and the dark frame this reduction detrends with*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller
        - ``arm`` -- the arm under reduction, read by the caller

        **Return:**

        - ``master_bias`` -- the master bias frame, or False when the set carries none
        - ``dark`` -- the dark frame, or False when the set carries none. The NIR
          lamp-off frame is read after the master dark, so it overrides one

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
        if self.inst.lower() == "soxs":
            add_filters = {kw("DPR_TYPE"): "WAVE,LAMP", kw("DPR_TECH"): "IMAGE"}
        else:
            add_filters = {kw("DPR_TYPE"): "LAMP,FMTCHK", kw("DPR_TECH"): "IMAGE"}
        # from tabulate import tabulate
        # self.log.print(tabulate(self.inputFrames.summary, headers='keys', tablefmt='github'))
        # self.log.print(self.inputFrames.files_filtered(include_path=True, **add_filters))
        # sys.exit(0)

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

    def _read_pinhole_frame(self, kw):
        """*read the single-pinhole frame this reduction fits a dispersion solution to*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller

        **Return:**

        - ``pinhole_image`` -- the last frame matching the arm's pinhole filters,
          or False when the set carries none

        **Usage:**

        ```python
        pinhole_image = self._read_pinhole_frame(self.kw)
        ```
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        pinhole_image = False

        if self.inst.lower() == "soxs":
            filter_list = [
                {kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
                {kw("DPR_TYPE"): "WAVE,LAMP", kw("DPR_TECH"): "ECHELLE,PINHOLE"},
            ]
        else:
            filter_list = [{kw("DPR_TYPE"): "LAMP,FMTCHK", kw("DPR_TECH"): "ECHELLE,PINHOLE"}]

        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                pinhole_image = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        return pinhole_image

    def _calibrate_pinhole_frame(self, master_bias, dark, pinhole_image):
        """*detrend the pinhole frame, update its keywords and, on request, write it out*

        **Key Arguments:**

        - ``master_bias`` -- the master bias frame, or False
        - ``dark`` -- the dark frame, or False
        - ``pinhole_image`` -- the raw single-pinhole frame

        Sets ``self.pinholeFrame``.

        **Usage:**

        ```python
        self._calibrate_pinhole_frame(master_bias, dark, pinhole_image)
        ```
        """
        self.pinholeFrame = self.detrend(inputFrame=pinhole_image, master_bias=master_bias, dark=dark)

        self.update_fits_keywords(frame=self.pinholeFrame)

        if self.settings["save-intermediate-products"]:
            outDir = self.workspaceRootPath
            filePath = self._write(
                frame=self.pinholeFrame,
                filedir=outDir,
                filename=False,
                overwrite=True,
                product=False,
            )
            self.log.print(f"\nCalibrated single pinhole frame: {filePath}\n")

    def _tune_dispersion_parameters(self):
        """*sweep the polynomial degree grid instead of producing a dispersion map*

        The line list is detected once and then reused across every permutation
        of four degrees drawn from 2 to 6. No product row is recorded and the
        output is not reported or cleaned up.

        **Return:**

        - ``productPath`` -- always None, because tuning writes no product
        - ``qcTable`` -- the quality-control table from the line-detection fit

        **Usage:**

        ```python
        productPath, qcTable = self._tune_dispersion_parameters()
        ```
        """
        from itertools import product

        digits = [2, 3, 4, 5, 6]
        perm = product(digits, repeat=4)
        try:
            os.remove("residuals.txt")
        except OSError as e:
            self.log.debug(f"produce_product: `os.remove('residuals.txt')` failed, continuing: {e}")

        # GET THE LINE DETECTION LIST BEFORE JUMPING TO PERMUTATIONS
        (
            mapPath,
            mapImagePath,
            res_plots,
            qcTable,
            productsTable,
            lineDetectionTable,
        ) = create_dispersion_map(
            log=self.log,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            pinholeFrame=self.pinholeFrame,
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            startNightDate=self.startNightDate,
        ).get()

        from fundamentals import fmultiprocess

        permList = list(perm)

        # DEFINE AN INPUT ARRAY

        print("TUNING SOXSPIPE\n")

        fmultiprocess(
            log=self.log,
            function=parameterTuning,
            inputArray=permList,
            poolSize=100,
            timeout=3600,
            recipeSettings=self.recipeSettings,
            settings=self.settings,
            pinholeFrame=self.pinholeFrame,
            qc=self.qc,
            products=self.products,
            sofName=self.sofName,
            lineDetectionTable=lineDetectionTable,
            turnOffMP=self.debug,
            mute=True,
            progressBar=True,
        )
        productPath = None

        return productPath, qcTable

    def _fit_dispersion_solution(self, kw):
        """*fit the dispersion solution and record its product and quality-control rows*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller

        Sets ``self.polyOrders``, ``self.products``, ``self.qc`` and ``self.dateObs``.

        **Return:**

        - ``productPath`` -- the path to the first guess dispersion map
        - ``qcTable`` -- the reported quality-control table

        **Usage:**

        ```python
        productPath, qcTable = self._fit_dispersion_solution(self.kw)
        ```
        """
        import pandas as pd

        if self.polyOrders:
            self.polyOrders = str(self.polyOrders)
            self.polyOrders = [int(digit) for digit in str(self.polyOrders)]
            self.recipeSettings["order-deg"] = self.polyOrders[:2]
            self.recipeSettings["wavelength-deg"] = self.polyOrders[2:4]

        (
            productPath,
            mapImagePath,
            res_plots,
            qcTable,
            productsTable,
            lineDetectionTable,
        ) = create_dispersion_map(
            log=self.log,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            pinholeFrame=self.pinholeFrame,
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            startNightDate=self.startNightDate,
            turnOffMP=self.turnOffMP,
        ).get()

        filename = os.path.basename(productPath)

        utcnow = utcnow_string()

        self.products = pd.concat([self.products, productsTable])
        self.qc = pd.concat([self.qc, qcTable])

        self.dateObs = self.pinholeFrame.header[kw("DATE_OBS")]

        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel="DISP_MAP",
            fileName=filename,
            filePath=productPath,
            productDesc=f"{self.arm} first pass dispersion solution",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS Table",
            label="PROD",
        )

        qcTable = self.report_output()
        self.clean_up()

        return productPath, qcTable


def parameterTuning(
    p,
    log,
    recipeSettings,
    settings,
    pinholeFrame,
    qc,
    products,
    sofName,
    lineDetectionTable,
):
    """*tuning the spatial solution*

    **Key Arguments:**

    - ``p`` -- one permutation of four polynomial degrees
    - ``log`` -- logger
    - ``recipeSettings`` -- the recipe settings dictionary, rewritten in place with ``p``
    - ``settings`` -- the settings dictionary
    - ``pinholeFrame`` -- the calibrated single-pinhole frame to fit
    - ``qc`` -- the quality-control table to pass to the dispersion map
    - ``products`` -- the products table to pass to the dispersion map
    - ``sofName`` -- the name of the set-of-files this reduction came from
    - ``lineDetectionTable`` -- the line detections to reuse across permutations

    The fit's own outputs are discarded. The residuals the dispersion map writes
    are the product of a tuning run.

    **Usage:**

    ```python
    parameterTuning(
        (3, 4, 5, 6),
        log=log,
        recipeSettings=recipeSettings,
        settings=settings,
        pinholeFrame=pinholeFrame,
        qc=qc,
        products=products,
        sofName=sofName,
        lineDetectionTable=lineDetectionTable,
    )
    ```
    """

    recipeSettings["order-deg"] = list(p[:2])
    recipeSettings["wavelength-deg"] = list(p[2:4])

    from soxspipe.commonutils import create_dispersion_map

    this = create_dispersion_map(
        log=log,
        settings=settings,
        recipeSettings=recipeSettings,
        pinholeFrame=pinholeFrame,
        qcTable=qc,
        productsTable=products,
        sofName=sofName,
        create2DMap=False,
        lineDetectionTable=lineDetectionTable,
        startNightDate=False,
    )
    (
        productPath,
        mapImagePath,
        res_plots,
        qcTable,
        productsTable,
        lineDetectionTable,
    ) = this.get()

    return
