#!/usr/bin/env python
"""
*enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit)*

Author
: David Young & Marco Landoni

Date Created
: March 17, 2021
"""

################# GLOBAL IMPORTS ####################
import os
import sys

from soxspipe.commonutils.toolkit import append_product, utcnow_string

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_spatial_solution(base_recipe):
    """
    *Enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit)*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``create2DMap`` -- create the 2D image map of wavelength, slit-position and order from disp solution.
    - ``polyOrders`` -- the orders of the x-y polynomials used to fit the dispersion solution. Overrides parameters found in the yaml settings file. e.g 345435 is order_x=3, order_y=4 ,wavelength_x=5 ,wavelength_y=4, slit_x=3 ,slit_y=5. Default *False*.
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- debug mode. True or False. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.


    See `produce_product` method for usage.

    """

    # INITIALISATION

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=[],
        verbose=False,
        overwrite=False,
        create2DMap=True,
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
            recipeName="soxs-spat-solution",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_spatial_solution' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose
        self.create2DMap = create2DMap
        self.polyOrders = polyOrders
        self.debug = debug

        self._parse_poly_orders()

        # xt-self-arg-tmpx

        self._collect_input_frames()
        self._verify_and_announce_input_frames()
        self._sort_and_report_input_frames()

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def _parse_poly_orders(self):
        """*validate the six-digit `polyOrders` override*

        Sets ``self.polyOrders``. A false value is left alone, so the recipe
        falls back to the degrees in the settings file.
        """
        if self.polyOrders is not False:
            value = self.polyOrders
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, str))
                or len(str(value)) != 6
                or not str(value).isascii()
                or not str(value).isdigit()
            ):
                raise TypeError("THE poly VALUE NEEDS TO BE A 6 DIGIT INTEGER")
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
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_spatial_solution - NO MORE, NO LESS.
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
        """*verify input frames match those required by the `soxs_spatial_solution` recipe*

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
        """*report why a NIR frame set is not a multi-pinhole arc set with its order and dispersion tables*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageTech`` -- the image techniques the basic verification classified
        - ``imageCat`` -- the product categories the basic verification classified
        - ``arm`` -- the arm under reduction, which names the tables required

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        # WANT ON AND OFF PINHOLE FRAMES
        if not error:
            for i in imageTypes:
                if i not in ["LAMP,WAVE", "LAMP,FLAT", "FLAT,LAMP", "WAVE,LAMP"]:
                    error = (
                        f"Found a {i} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE. Can "
                        "optionally supply a master-flat for NIR."
                    )

        if not error:
            for i in imageTech:
                if i not in [
                    "ECHELLE,MULTI-PINHOLE",
                    "IMAGE",
                    "ECHELLE,SLIT",
                    "ECHELLE,PINHOLE",
                ]:
                    error = (
                        f"Found a {i} file. Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on "
                        "and lamp off frames, a first-guess dispersion solution table and an order location table for "
                        "NIR. Can optionally supply a master-flat for NIR."
                    )

        if not error and "LAMP,WAVE" not in imageTypes and "WAVE,LAMP" not in imageTypes:
            error = (
                "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames, a "
                "first-guess dispersion solution table and an order location table for NIR. Can optionally supply a "
                "master-flat for NIR."
            )

        if not error and "ECHELLE,MULTI-PINHOLE" not in imageTech:
            error = (
                "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames, a "
                "first-guess dispersion solution table and an order location table for NIR. Can optionally supply a "
                "master-flat for NIR."
            )

        if not error:
            for i in [f"ORDER_TAB_{arm}", f"DISP_TAB_{arm}"]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe spatial_solution need to be LAMP,WAVE lamp on and lamp off frames, "
                        "a first-guess dispersion solution table and an order location table for NIR. Can optionally "
                        "supply a master-flat for NIR."
                    )

        return error

    def _uvb_vis_input_frame_error(self, imageTypes, imageCat, arm):
        """*report why a UVB or VIS frame set is not an arc set with its bias, order and dispersion tables*

        **Key Arguments:**

        - ``imageTypes`` -- the image types the basic verification classified
        - ``imageCat`` -- the product categories the basic verification classified
        - ``arm`` -- the arm under reduction, which names the calibrations required

        **Return:**

        - ``error`` -- the rejection message, or False when the set is acceptable
        """
        error = False

        if not error:
            for i in imageTypes:
                if i not in ["LAMP,WAVE", "LAMP,FLAT", "WAVE,LAMP"]:
                    error = (
                        f"Found a {i} frame. Input frames for soxspipe spatial_solution need to be LAMP,WAVE and a "
                        "master-bias, a first-guess dispersion solution table and an order location table. Can "
                        "optionally supply a master-flat and/or master-dark for UVB/VIS."
                    )

        if not error:
            for i in [
                f"MASTER_BIAS_{arm}",
                f"ORDER_TAB_{arm}",
                f"DISP_TAB_{arm}",
            ]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe spatial_solution need to be LAMP,WAVE, a master-bias, a "
                        "first-guess dispersion solution table and an order location table. Can optionally supply a "
                        "master-flat and/or master-dark for UVB/VIS."
                    )

        return error

    def produce_product(self):
        """*generate the 2D dispersion map*

        **Return:**

        - ``mapImagePath`` -- the path to the 2D detector map image, or None when no image was made
        - ``qcTable`` -- the reported quality-control table

        When the pipeline is tuning rather than reducing, the method returns
        ``(None, None, None)`` instead of the pair.

        **Usage**

        ```python
        from soxspipe.recipes import soxs_spatial_solution
        recipe = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=fileList
        )
        mapImagePath, qcTable = recipe.produce_product()
        ```
        """
        self.log.debug("starting the ``produce_product`` method")

        # TEMPORARY WARNING
        # if self.inst.upper() == "SOXS" and self.arm.upper() == "VIS":
        #    self.log.warning("The SOXS UVVIS Multi-Pinhole line-list is not yet ready. "
        #                     "It will be included in a future code release")
        #    return None, None, None

        arm = self.arm
        kw = self.kw

        master_bias, dark, master_flat = self._read_calibration_frames(kw, arm)
        multi_pinhole_image = self._read_multi_pinhole_frame(kw)

        self.dateObs = multi_pinhole_image.header[kw("DATE_OBS")]

        slit_arc = self._read_slit_arc_frame(kw)

        # FIND THE ORDER TABLE
        filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}
        order_table = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)[0]

        add_filters = {kw("PRO_CATG"): f"DISP_TAB_{arm}".upper()}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            disp_map_table = i

        self._calibrate_frames(multi_pinhole_image, slit_arc, master_bias, dark, master_flat, order_table)

        if self.settings["tune-pipeline"]:
            from itertools import product

            order = [2, 3, 4, 5]
            wavelength = [2, 3, 4, 5]
            slit = [1, 2, 3]
            # perm = product([self.recipeSettings["order-deg"][0]], [self.recipeSettings["order-deg"][1]],
            #     [self.recipeSettings["wavelength-deg"][0]], [self.recipeSettings["wavelength-deg"][1]], slit, slit)
            perm = product(order, order, wavelength, wavelength, slit, slit)
            try:
                os.remove("residuals.txt")
            except OSError as e:
                self.log.debug(f"produce_product: `os.remove('residuals.txt')` failed, continuing: {e}")

            self._tune_spatial_parameters(perm, disp_map_table, order_table)
            return None, None, None
        if self.polyOrders:
            self.polyOrders = str(self.polyOrders)
            self.polyOrders = [int(digit) for digit in str(self.polyOrders)]
            self.recipeSettings["order-deg"] = self.polyOrders[:2]
            self.recipeSettings["wavelength-deg"] = self.polyOrders[2:4]
            self.recipeSettings["slit-deg"] = self.polyOrders[4:]

        if self.debug:
            self.create2DMap = False
            self.slit_arc = False

        mapImagePath = self._fit_spatial_solution(disp_map_table, order_table)

        qcTable = self.report_output()
        self.clean_up()

        self.log.debug("completed the ``produce_product`` method")
        return mapImagePath, qcTable

    def _read_calibration_frames(self, kw, arm):
        """*read the master bias, the dark and the master flat this reduction detrends with*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller
        - ``arm`` -- the arm under reduction, read by the caller

        **Return:**

        - ``master_bias`` -- the master bias frame, or False when the set carries none
        - ``dark`` -- the dark frame, or False when the set carries none. The lamp-off
          frame is read after the master dark, so it overrides one
        - ``master_flat`` -- the master flat frame, or False when the set carries none

        **Usage:**

        ```python
        master_bias, dark, master_flat = self._read_calibration_frames(self.kw, self.arm)
        ```
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        master_bias = False
        dark = False
        master_flat = False

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
        if self.inst.upper() == "SOXS":
            add_filters = {kw("DPR_TYPE"): "WAVE,LAMP", kw("DPR_TECH"): "IMAGE"}
        else:
            add_filters = {kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "IMAGE"}
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

        # UVB/VIS/NIR FLAT
        add_filters = {kw("PRO_CATG"): "MASTER_FLAT_" + arm}
        for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
            master_flat = CCDData.read(
                i,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        return master_bias, dark, master_flat

    def _read_multi_pinhole_frame(self, kw):
        """*read the multi-pinhole arc frame the spatial solution is fitted to*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller

        **Return:**

        - ``multi_pinhole_image`` -- the last frame matching the instrument's filter list,
          or False when the set carries none
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        multi_pinhole_image = False

        # MULTIPINHOLE IMAGE
        if self.inst.upper() == "SOXS":
            filter_list = [
                {kw("DPR_TYPE"): "WAVE,LAMP", kw("DPR_TECH"): "ECHELLE,MULTI-PINHOLE"},
                {kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "ECHELLE,MULTI-PINHOLE"},
            ]
        else:
            filter_list = [{kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "ECHELLE,MULTI-PINHOLE"}]

        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                multi_pinhole_image = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        return multi_pinhole_image

    def _read_slit_arc_frame(self, kw):
        """*read the optional slit arc frame used for the resolution plot*

        **Key Arguments:**

        - ``kw`` -- the recipe's FITS keyword lookup, read by the caller

        **Return:**

        - ``slit_arc`` -- the last frame matching the instrument's filter list, or False
          when the set carries none

        Sets ``self.slit_arc`` to the raw frame, or None when there is none.
        """
        from astropy import units as u
        from astropy.nddata import CCDData

        # DO WE HAVE A SLIT ARC?
        slit_arc = False
        if self.inst.upper() == "SOXS":
            filter_list = [
                {kw("DPR_TYPE"): "WAVE,LAMP", kw("DPR_TECH"): "ECHELLE,SLIT"},
                {kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "ECHELLE,SLIT"},
            ]
        else:
            filter_list = [{kw("DPR_TYPE"): "LAMP,WAVE", kw("DPR_TECH"): "ECHELLE,SLIT"}]
        for add_filters in filter_list:
            for i in self.inputFrames.files_filtered(include_path=True, **add_filters):
                slit_arc = CCDData.read(
                    i,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        # CHECK IF SLIT ARC IS PRESENT
        if slit_arc:
            # print(f"The 1.0 arc slit arc-lamp frame is {i}. Resolution plot will be given for 1.0'' arcsec slit")
            self.slit_arc = slit_arc
        else:
            self.slit_arc = None

        return slit_arc

    def _calibrate_frames(self, multi_pinhole_image, slit_arc, master_bias, dark, master_flat, order_table):
        """*detrend the multi-pinhole frame and any slit arc, stamp keywords and optionally save*

        **Key Arguments:**

        - ``multi_pinhole_image`` -- the raw multi-pinhole frame
        - ``slit_arc`` -- the raw slit arc frame, or False
        - ``master_bias`` -- the master bias frame, or False
        - ``dark`` -- the dark frame, or False
        - ``master_flat`` -- the master flat frame, or False. Dropped when ``use_flat`` is off
        - ``order_table`` -- the path to the order table

        Sets ``self.multiPinholeFrame`` and ``self.slit_arc``.
        """
        from soxspipe.commonutils.toolkit import quicklook_image

        if not self.recipeSettings["use_flat"]:
            master_flat = False
        self.multiPinholeFrame = self.detrend(
            inputFrame=multi_pinhole_image,
            master_bias=master_bias,
            dark=dark,
            master_flat=master_flat,
            order_table=order_table,
        )

        # DETREND THE slit_arc IF PRESENT

        if slit_arc:
            # print('Detrending slit arc...')
            self.slit_arc = self.detrend(
                inputFrame=slit_arc,
                master_bias=master_bias,
                dark=dark,
                master_flat=master_flat,
                order_table=order_table,
            )
        else:
            self.slit_arc = None

        if False:
            quicklook_image(
                log=self.log,
                CCDObject=self.slit_arc,
                show=True,
                ext=False,
                stdWindow=1,
                title="Multi-pinhole Frame Overlaid with Dispersion Solution",
                settings=self.settings,
            )

        # INJECT KEYWORDS INTO HEADER
        self.update_fits_keywords(frame=self.multiPinholeFrame)

        if self.settings["save-intermediate-products"]:
            fileDir = self.workspaceRootPath
            filepath = self._write(
                self.multiPinholeFrame,
                fileDir,
                filename=False,
                overwrite=True,
                product=False,
            )
            self.log.print(f"\nCalibrated multi pinhole frame frame saved to {filepath}\n")

    def _tune_spatial_parameters(self, perm, disp_map_table, order_table):
        """*fit the line list once, then sweep the polynomial degrees over the grid*

        **Key Arguments:**

        - ``perm`` -- the iterator of six-degree permutations to try
        - ``disp_map_table`` -- the path to the first-guess dispersion table
        - ``order_table`` -- the path to the order table
        """
        from soxspipe.commonutils import create_dispersion_map

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
            pinholeFrame=self.multiPinholeFrame,
            firstGuessMap=disp_map_table,
            orderTable=order_table,
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            create2DMap=False,
            startNightDate=self.startNightDate,
            debug=self.debug,
        ).get()

        print("\n\nTUNING SOXSPIPE\n")

        from fundamentals import fmultiprocess

        # DEFINE AN INPUT ARRAY
        fmultiprocess(
            log=self.log,
            function=parameterTuning,
            inputArray=list(perm),
            poolSize=False,
            timeout=360000,
            recipeSettings=self.recipeSettings,
            settings=self.settings,
            multiPinholeFrame=self.multiPinholeFrame,
            disp_map_table=disp_map_table,
            order_table=order_table,
            qc=self.qc,
            products=self.products,
            sofName=self.sofName,
            lineDetectionTable=lineDetectionTable,
            turnOffMP=self.debug,
            mute=True,
            progressBar=True,
        )

    def _fit_spatial_solution(self, disp_map_table, order_table):
        """*fit the full dispersion-spatial solution, record its products and show the quick-look*

        **Key Arguments:**

        - ``disp_map_table`` -- the path to the first-guess dispersion table
        - ``order_table`` -- the path to the order table

        **Return:**

        - ``mapImagePath`` -- the path to the 2D map image, or None when none was made

        Sets ``self.products`` and ``self.qc``.
        """
        import pandas as pd

        from soxspipe.commonutils import create_dispersion_map
        from soxspipe.commonutils.toolkit import quicklook_image

        # GENERATE AN UPDATED DISPERSION MAP
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
            pinholeFrame=self.multiPinholeFrame,
            firstGuessMap=disp_map_table,
            orderTable=order_table,
            qcTable=self.qc,
            productsTable=self.products,
            sofName=self.sofName,
            create2DMap=self.create2DMap,
            startNightDate=self.startNightDate,
            arcFrame=self.slit_arc,
            debug=self.debug,
            turnOffMP=self.turnOffMP,
        ).get()

        filename = os.path.basename(mapPath)

        utcnow = utcnow_string()

        self.products = pd.concat([self.products, productsTable])
        self.qc = pd.concat([self.qc, qcTable])

        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel="SPAT_SOL",
            fileName=filename,
            filePath=mapPath,
            productDesc=f"{self.arm} full dispersion-spatial solution",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        if mapImagePath:
            filename = os.path.basename(mapImagePath)
            self.products = append_product(
                self.products,
                recipeName=self.recipeName,
                productLabel="2D_MAP",
                fileName=filename,
                filePath=mapImagePath,
                productDesc=f"{self.arm} 2D detector map of wavelength, slit position and order",
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow,
                fileType="FITS",
                label="PROD",
            )

        # INSPECT THE MAP AGAINST THE MULTIPINHOLE FRAME
        quicklook_image(
            log=self.log,
            CCDObject=self.multiPinholeFrame,
            show=False,
            ext=False,
            stdWindow=1,
            title="Multi-pinhole Frame Overlaid with Dispersion Solution",
            surfacePlot=True,
            dispMap=mapPath,
            dispMapImage=mapImagePath,
            settings=self.settings,
            skylines=False,
        )

        return mapImagePath


# THE NAME IS PUBLIC: EACH RECIPE PASSES IT TO fmultiprocess AND tests/unit/test_parameter_tuning.py
# CALLS IT BY NAME, AS IN soxs_disp_solution AND soxs_order_centres
def parameterTuning(  # noqa: N802
    p,
    log,
    recipeSettings,
    settings,
    multiPinholeFrame,
    disp_map_table,
    order_table,
    qc,
    products,
    sofName,
    lineDetectionTable,
):
    """*tuning the spatial solution*

    **Key Arguments:**

    - ``p`` -- one permutation of the six polynomial degrees: order, wavelength and slit pairs
    - ``log`` -- logger
    - ``recipeSettings`` -- the recipe settings dictionary, rewritten in place with ``p``
    - ``settings`` -- the settings dictionary
    - ``multiPinholeFrame`` -- the calibrated multi-pinhole frame to fit
    - ``disp_map_table`` -- the path to the first-guess dispersion table
    - ``order_table`` -- the path to the order table
    - ``qc`` -- the quality-control table to pass to the dispersion map
    - ``products`` -- the products table to pass to the dispersion map
    - ``sofName`` -- the name of the set-of-files this reduction came from
    - ``lineDetectionTable`` -- the line detections to reuse across permutations

    The fit's own outputs are discarded.

    **Usage:**

    ```python
    parameterTuning(
        (3, 4, 5, 4, 3, 5),
        log=log,
        recipeSettings=recipeSettings,
        settings=settings,
        multiPinholeFrame=multiPinholeFrame,
        disp_map_table=disp_map_table,
        order_table=order_table,
        qc=qc,
        products=products,
        sofName=sofName,
        lineDetectionTable=lineDetectionTable,
    )
    ```
    """

    recipeSettings["order-deg"] = list(p[:2])
    recipeSettings["wavelength-deg"] = list(p[2:4])
    recipeSettings["slit-deg"] = list(p[4:6])

    from soxspipe.commonutils import create_dispersion_map

    this = create_dispersion_map(
        log=log,
        settings=settings,
        recipeSettings=recipeSettings,
        pinholeFrame=multiPinholeFrame,
        firstGuessMap=disp_map_table,
        orderTable=order_table,
        qcTable=qc,
        productsTable=products,
        sofName=sofName,
        create2DMap=False,
        lineDetectionTable=lineDetectionTable,
        startNightDate=False,
        debug=self.debug,
    )
    try:
        (
            productPath,
            mapImagePath,
            res_plots,
            qcTable,
            productsTable,
            lineDetectionTable,
        ) = this.get()
    except Exception as e:
        log.warning(f"parameterTuning: this tuning iteration failed and records nothing in the grid, continuing: {e}")

    return
