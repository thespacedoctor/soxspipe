#!/usr/bin/env python
"""
*generate a single normalised master flat-field frame*

Author
: David Young & Marco Landoni

Date Created
: September 16, 2020
"""

#
import os
import sys
from os.path import expanduser

from soxspipe.commonutils import detect_order_edges, subtract_background
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.toolkit import (
    append_product,
    append_qc,
    generic_quality_checks,
    quicklook_image,
    spectroscopic_image_quality_checks,
    unpack_order_table,
    utcnow_string,
)

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"


class soxs_mflat(base_recipe):
    """
    *generate a single normalised master flat-field frame*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- generate debug plots. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.


    **Usage**

    ```python
    from soxspipe.recipes import soxs_mflat
    recipe = soxs_mflat(
        log=log,
        settings=settings,
        inputFrames=fileList
    )
    mflatFrame = recipe.produce_product()
    ```
    """

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=[],
        verbose=False,
        overwrite=False,
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
            recipeName="soxs-mflat",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        self.log = log
        log.debug("instantiating a new 'soxs_mflat' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose

        self._collect_input_frames()
        self._verify_and_announce_input_frames()
        self._sort_and_report_input_frames()

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def _collect_input_frames(self):
        """*convert the input files to a ccdproc image collection*

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

        return

    def _verify_and_announce_input_frames(self):
        """*verify the collected frames and report the result to the user*

        Sets ``self.imageType``, through ``verify_input_frames``.
        """
        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY soxs_mflat - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        return

    def _sort_and_report_input_frames(self):
        """*sort the image collection by observation date and, when verbose, print it*"""
        # SORT IMAGE COLLECTION BY MJD
        self.inputFrames.sort(["MJD-OBS"])
        if self.verbose:
            self.log.print("# RAW INPUT FRAMES - SUMMARY")
            self.log.print(self.inputFrames.summary)
            self.log.print("\n")

        return

    def verify_input_frames(self):
        """*verify the input frames match those required by the soxs_mflat recipe*

        If the fits files conform to required input for the recipe everything will pass silently, otherwise an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        kw = self.kw

        error = False

        import warnings

        from astropy.utils.exceptions import AstropyWarning

        warnings.simplefilter("ignore", AstropyWarning)

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        if self.arm == "NIR":
            error = self._nir_input_frame_error(imageTypes, imageTech)
        else:
            error = self._uvb_vis_input_frame_error(imageTypes, imageCat)

        # UV-VIS NEEDS BOTH D AND Q-LAMPS
        if not error and self.inst.upper() != "SOXS":
            if "LAMP,QFLAT" in imageTypes or "LAMP,DFLAT" in imageTypes:
                if "LAMP,QFLAT" in imageTypes and "LAMP,DFLAT" in imageTypes:
                    pass
                else:
                    for i in imageTypes:
                        if "LAMP" in i:
                            error = f'Only "{i}" image types found. Please include both D2 and QTH lamp flats'

        # LOOK FOR ORDER TABLE
        if not error:
            arm = self.arm
            if f"ORDER_TAB_{arm}" not in imageCat:
                error = f"Need an order centre for {arm} - none found with the input files"

        # CHECK EXPTIME
        if not error:
            lamps = ["LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT", "DOME,FLAT"]
            for l in lamps:
                filterDict = {kw("DPR_TYPE"): l, kw("DPR_TECH"): "ECHELLE,SLIT"}
                flatCollection = self.inputFrames.filter(**filterDict)
                if len(flatCollection.files):
                    exptime = flatCollection.values(keyword=kw("EXPTIME"), unique=True)
                    if len(exptime) > 1:
                        error = f"Input {l} frames for soxspipe mflat need to have a unique exptime"

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
        """*check the NIR input frame types and techniques against those the recipe accepts*

        **Key Arguments:**

        - ``imageTypes`` -- the image types of the input frames
        - ``imageTech`` -- the observing techniques of the input frames

        **Return:**

        - ``error`` -- the error message for the last failed check, or False when every check passes
        """
        error = False

        # WANT ON AND OFF PINHOLE FRAMES
        # MIXED INPUT IMAGE TYPES ARE BAD
        if not error:
            if len(imageTypes) > 1:
                # FIX ME
                if len(imageTypes) == 2 and ("DARK" in imageTypes):
                    pass
                else:
                    pass
                    # imageTypes = " and ".join(imageTypes)
                    # error = "Input frames are a mix of %(imageTypes)s" % locals()

        if not error:
            if "LAMP,FLAT" not in imageTypes and "FLAT,LAMP" not in imageTypes:
                error = (
                    "Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR" % locals()
                )

        if not error:
            for i in imageTech:
                if i not in ["ECHELLE,SLIT", "IMAGE"]:
                    error = (
                        f"Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. You have provided {i}"
                        % locals()
                    )

        if not error:
            for i in ["ECHELLE,SLIT", "IMAGE"]:
                if i not in imageTech:
                    error = (
                        f"Input frames for soxspipe mflat need to be flat-lamp on and lamp off frames for NIR. You have are missing TECH={i}"
                        % locals()
                    )

        return error

    def _uvb_vis_input_frame_error(self, imageTypes, imageCat):
        """*check the UVB/VIS input frame types and catalogues against those the recipe requires*

        **Key Arguments:**

        - ``imageTypes`` -- the image types of the input frames
        - ``imageCat`` -- the product categories of the input frames

        **Return:**

        - ``error`` -- the error message for the last failed check, or False when every check passes
        """
        error = False

        if not error:
            for i in imageTypes:
                if i not in [
                    "LAMP,FLAT",
                    "LAMP,QFLAT",
                    "LAMP,DFLAT",
                    "FLAT,LAMP",
                    "DOME,FLAT",
                ]:
                    error = (
                        "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS"
                        % locals()
                    )

        if not error:
            for i in [f"MASTER_BIAS_{self.arm}", f"ORDER_TAB_{self.arm}"]:
                if i not in imageCat:
                    error = (
                        "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS"
                        % locals()
                    )

        if not error:
            found = False
            for i in [
                "LAMP,FLAT",
                "LAMP,QFLAT",
                "LAMP,DFLAT",
                "FLAT,LAMP",
                "DOME,FLAT",
            ]:
                if i in imageTypes:
                    found = True
            if not found:
                error = (
                    "Input frames for soxspipe mflat need to be flat-lamp frames,a master-bias frame, an order-locations tables and possibly a master dark for UVB/VIS"
                    % locals()
                )

        return error

    def produce_product(self):
        """*generate the master flat frames updated order location table (with egde detection)*

        **Return:**

        - ``productPath`` -- the path to the master flat frame
        """
        self.log.debug("starting the ``produce_product`` method")

        import pandas as pd

        productPath = None
        arm = self.arm
        kw = self.kw

        home = expanduser("~")
        outDir = self.settings["workspace-root-dir"].replace("~", home)

        # CALIBRATE THE FRAMES BY SUBTRACTING BIAS AND/OR DARK
        calibratedFlats, dcalibratedFlats, qcalibratedFlats, domecalibratedFlats = self.calibrate_frame_set()

        allCalibratedFlats = calibratedFlats + dcalibratedFlats + qcalibratedFlats + domecalibratedFlats

        calibratedFlatSet = [
            calibratedFlats,
            dcalibratedFlats,
            qcalibratedFlats,
            domecalibratedFlats,
        ]
        flatKeywords = [
            "LAMP,ORDERDEF",
            "LAMP,DORDERDEF",
            "LAMP,QORDERDEF",
            "LAMP,QORDERDEF",
        ]
        lampTag = ["", "_DLAMP", "_QLAMP", "_DOME"]
        filelists = [
            self.calibratedFlatFiles,
            self.dFlatFiles,
            self.qFlatFiles,
            self.domeFlatFiles,
        ]
        normalisedFlatSet = []
        self.combinedNormalisedFlatSet = []
        self.masterFlatSet = []
        self.orderTableSet = []
        self.detectionCountSet = []
        medianOrderFluxDFExists = False

        productTable = self.products
        qcTable = self.qc

        for cf, fk, tag, files in zip(calibratedFlatSet, flatKeywords, lampTag, filelists):

            if len(cf) == 0:
                self.orderTableSet.append(None)
                normalisedFlatSet.append(None)
                self.combinedNormalisedFlatSet.append(None)
                self.masterFlatSet.append(None)
                continue

            if tag and self.inst.upper() != "SOXS":
                filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}", kw("OBJECT"): fk}
            else:
                filterDict = {kw("PRO_CATG"): f"ORDER_TAB_{arm}"}

            orderTablePaths = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)
            if len(orderTablePaths) > 0:
                orderTablePath = orderTablePaths[0]
                thisPath = orderTablePath
            else:
                self.orderTableSet.append(None)

            combined_normalised_flat = self._normalise_and_stack_lamp_flats(cf, orderTablePath, tag)

            self.combinedNormalisedFlatSet.append(combined_normalised_flat.copy())

            self.update_fits_keywords(frame=combined_normalised_flat, rawFrames=files)

            qcTable, orderTablePath = self._detect_lamp_order_edges(combined_normalised_flat, orderTablePath, tag)

            self.dateObs = combined_normalised_flat.header[self.kw("DATE_OBS")]

            # UNPACK THE ORDER TABLE
            orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
                log=self.log, orderTablePath=orderTablePath
            )

            if tag in ("_DLAMP", "_QLAMP"):
                writeQC = False
            else:
                writeQC = True

            if self.recipeSettings["subtract_background"]:
                backgroundFrame, combined_normalised_flat = self._subtract_lamp_background(
                    combined_normalised_flat, orderTablePath, tag, outDir
                )

            mflat, medianOrderFluxDF = self.mask_low_sens_pixels(
                frame=combined_normalised_flat,
                orderTablePath=orderTablePath,
                returnMedianOrderFlux=True,
                writeQC=writeQC,
            )

            self.masterFlatSet.append(mflat)

            productPath = self._write_lamp_master_flat(mflat, outDir, tag)

            if tag:
                medianOrderFluxDF.rename(columns={"medianFlux": tag}, inplace=True)

            if tag and not medianOrderFluxDFExists:
                medianOrderFluxDFExists = True
                medianOrderFluxDFFirst = medianOrderFluxDF.copy()
            elif tag:
                medianOrderFluxDF = pd.merge(medianOrderFluxDFFirst, medianOrderFluxDF)

        # UV-STITCHING
        if len(self.detectionCountSet) > 1:
            mflat = self.stitch_uv_mflats(medianOrderFluxDF, orderTablePath=thisPath)
        else:
            self.qc = pd.concat([self.qc, qcTable], ignore_index=True)

        quicklook_image(
            log=self.log,
            CCDObject=combined_normalised_flat,
            show=False,
            ext=None,
            surfacePlot=True,
            title="Final master flat frame",
        )

        self.update_fits_keywords(frame=mflat)

        # WRITE MFLAT TO FILE
        productPath = self._write(mflat, outDir, overwrite=True)

        utcnow = utcnow_string()
        basename = os.path.basename(productPath)
        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel="MFLAT",
            fileName=basename,
            filePath=productPath,
            productDesc=f"{self.arm} master spectroscopic flat frame",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        if 1 == 0:
            filename = filenamer(log=self.log, frame=mflat, settings=self.settings)
            filename = filename.replace(".fits", "_background.fits")
            filepath = self._write(backgroundFrame, outDir, filename=filename, overwrite=True)
            filepath = os.path.abspath(filepath)
            self.products = pd.concat(
                [
                    self.products,
                    pd.DataFrame([
                        {
                            "soxspipe_recipe": self.recipeName,
                            "product_label": "",
                            "file_name": filename,
                            "file_type": "FITS",
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "product_desc": "modelled scatter background light image (removed from master flat)",
                            "file_path": backgroundFrame,
                            "label": "PROD",
                        }
                    ]),
                ],
                ignore_index=True,
            )

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log,
            frame=mflat,
            settings=self.settings,
            recipeName=self.recipeName,
            qcTable=self.qc,
        )
        self.qc = spectroscopic_image_quality_checks(
            log=self.log,
            frame=mflat,
            settings=self.settings,
            recipeName=self.recipeName,
            qcTable=self.qc,
            orderTablePath=orderTablePath,
        )

        qcTable = self.report_output()
        self.clean_up()

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _normalise_and_stack_lamp_flats(self, cf, orderTablePath, tag):
        """*normalise and stack one lamp's flats, then renormalise against that first-pass stack and stack again*

        **Key Arguments:**

        - ``cf`` -- the calibrated flat frames of this lamp
        - ``orderTablePath`` -- path to the order table used to normalise the frames
        - ``tag`` -- the lamp tag, used in plot titles

        **Return:**

        - ``combined_normalised_flat`` -- the stack of the re-normalised flat frames
        """
        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL
        normalisedFlats = self.normalise_flats(cf, orderTablePath=orderTablePath, lamp=tag)

        quicklook_image(
            log=self.log,
            CCDObject=normalisedFlats[0],
            stdWindow=6,
            show=False,
            ext=None,
            surfacePlot=True,
            title=f"Single normalised flat frame {tag}",
        )
        # STACK THE NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats,
            recipe="soxs_mflat",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )
        quicklook_image(
            log=self.log,
            CCDObject=combined_normalised_flat,
            stdWindow=6,
            show=False,
            ext=None,
            surfacePlot=True,
            title=f"Combined normalised flat frames {tag}",
        )

        # DIVIDE THROUGH BY FIRST-PASS MASTER FRAME TO REMOVE CROSS-PLANE
        # ILLUMINATION VARIATIONS
        # DETERMINE THE MEDIAN EXPOSURE FOR EACH FLAT FRAME AND NORMALISE THE
        # FLUX TO THAT LEVEL (AGAIN!)
        self.log.print("\n# DIVIDING EACH ORIGINAL FLAT FRAME BY FIRST PASS MASTER FLAT")

        normalisedFlats = self.normalise_flats(
            cf,
            orderTablePath=orderTablePath,
            firstPassMasterFlat=combined_normalised_flat,
            lamp=tag,
        )

        quicklook_image(
            log=self.log,
            CCDObject=normalisedFlats[0],
            show=False,
            ext=None,
            surfacePlot=True,
            title=f"Single re-normalised flat frame {tag}",
        )

        # STACK THE RE-NORMALISED FLAT FRAMES
        combined_normalised_flat = self.clip_and_stack(
            frames=normalisedFlats,
            recipe="soxs_mflat",
            ignore_input_masks=False,
            post_stack_clipping=False,
        )

        quicklook_image(
            log=self.log,
            CCDObject=combined_normalised_flat,
            show=False,
            ext=None,
            surfacePlot=True,
            title=f"Recombined normalised flat frames {tag}",
        )

        return combined_normalised_flat

    def _detect_lamp_order_edges(self, combined_normalised_flat, orderTablePath, tag):
        """*detect the order edges on one lamp's combined flat and record the updated order locations table*

        Sets ``self.products`` and appends to ``self.detectionCountSet`` and ``self.orderTableSet``.

        **Key Arguments:**

        - ``combined_normalised_flat`` -- the combined normalised flat frame of this lamp
        - ``orderTablePath`` -- path to the order centre table
        - ``tag`` -- the lamp tag

        **Return:**

        - ``qcTable`` -- the QC table returned by the edge detection
        - ``orderTablePath`` -- path to the order locations table the edge detection wrote
        """
        # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
        edges = detect_order_edges(
            log=self.log,
            flatFrame=combined_normalised_flat,
            orderCentreTable=orderTablePath,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            qcTable=self.qc,
            productsTable=self.products,
            tag=tag,
            sofName=self.sofName,
            binx=self.binRatioX,
            biny=self.binRatioY,
            lampTag=tag,
            startNightDate=self.startNightDate,
        )
        self.products, qcTable, orderDetectionCounts = edges.get()

        if tag:
            # NEED TO TRY AND RENAME BOTH ORDER AND COUNT COLUMNS FOR PANDAS 1.X and 2.X
            orderDetectionCounts.rename(columns={"order": tag}, inplace=True)
            orderDetectionCounts.rename(columns={"count": tag}, inplace=True)
            orderDetectionCounts.index.names = ["order"]

        self.detectionCountSet.append(orderDetectionCounts)

        mask = self.products["product_label"] == f"ORDER_LOC{tag}"
        orderTablePath = self.products.loc[mask]["file_path"].values[0]

        self.orderTableSet.append(orderTablePath)

        return qcTable, orderTablePath

    def _subtract_lamp_background(self, combined_normalised_flat, orderTablePath, tag, outDir):
        """*remove the scattered background light from one lamp's combined flat and write the background image*

        Sets ``self.products``.

        **Key Arguments:**

        - ``combined_normalised_flat`` -- the combined normalised flat frame of this lamp
        - ``orderTablePath`` -- path to the order locations table
        - ``tag`` -- the lamp tag
        - ``outDir`` -- the directory the background image is written to

        **Return:**

        - ``backgroundFrame`` -- the modelled background light image
        - ``combined_normalised_flat`` -- the combined flat frame with the background removed
        """
        import copy

        background = subtract_background(
            log=self.log,
            frame=combined_normalised_flat,
            sofName=self.sofName,
            recipeName=self.recipeName,
            orderTable=orderTablePath,
            settings=self.settings,
            productsTable=self.products,
            qcTable=self.qc,
            lamp=tag,
            startNightDate=self.startNightDate,
        )
        backgroundFrame, combined_normalised_flat, self.products = background.subtract()

        quicklook_image(
            log=self.log,
            CCDObject=backgroundFrame,
            show=False,
            ext="data",
            stdWindow=3,
            title="Background Light",
            surfacePlot=True,
        )

        utcnow = utcnow_string()

        backgroundFrame.header = copy.deepcopy(combined_normalised_flat.header)
        backgroundQCImage = self.sofName + "_BKGROUND.fits"
        filepath = self._write(backgroundFrame, outDir, filename=backgroundQCImage, overwrite=True)
        # filepath = os.path.abspath(filepath)
        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel="BKGROUND",
            fileName=backgroundQCImage,
            filePath=filepath,
            productDesc="modelled scatter background light image (removed from master flat)",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="QC",
        )

        return backgroundFrame, combined_normalised_flat

    def _write_lamp_master_flat(self, mflat, outDir, tag):
        """*write one lamp's master flat to file and record it in the products table*

        Sets ``self.products``.

        **Key Arguments:**

        - ``mflat`` -- the master flat frame of this lamp
        - ``outDir`` -- the directory the frame is written to
        - ``tag`` -- the lamp tag

        **Return:**

        - ``productPath`` -- the path the master flat was written to
        """
        # WRITE MFLAT TO FILE
        productPath = self._write(mflat.copy(), outDir, filename=self.sofName + ".fits", overwrite=True)

        utcnow = utcnow_string()
        basename = os.path.basename(productPath)

        if len(tag):
            product_desc = f"{self.arm} master spectroscopic flat frame ({tag.replace('_', '')})"
        else:
            product_desc = f"{self.arm} master spectroscopic flat frame"

        self.products = append_product(
            self.products,
            recipeName=self.recipeName,
            productLabel=f"MFLAT{tag}",
            fileName=basename,
            filePath=productPath,
            productDesc=product_desc,
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        return productPath

    def calibrate_frame_set(self):
        """*given all of the input data calibrate the frames by subtracting bias and/or dark*

        **Return:**

        - ``calibratedFlats`` -- the calibrated frames
        """
        self.log.debug("starting the ``calibrate_frame_set`` method")

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        bias = self._find_master_bias(kw)

        darkCollection = self._find_dark_collection(kw)

        # FIND THE FLAT FRAMES
        if self.arm.upper() == "NIR" or (self.inst.upper() != "SOXS" and self.arm.upper() == "VIS"):
            filter_list = [
                {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "ECHELLE,SLIT"},
                {kw("DPR_TYPE"): "LAMP,FLAT", kw("DPR_TECH"): "ECHELLE,SLIT"},
            ]
            for add_filters in filter_list:
                flatCollection = self.inputFrames.filter(**add_filters)
                if len(flatCollection.files) > 0:
                    break
        else:
            filterDict = {kw("DPR_TYPE"): "JUNK_DO_NOT_MATCH"}
            flatCollection = self.inputFrames.filter(**filterDict)

        if self.inst.upper() == "SOXS":
            filterDict = {kw("LAMP2"): "Deut_Lamp", kw("DPR_TECH"): "ECHELLE,SLIT"}
            dflatCollection = self.inputFrames.filter(**filterDict)
            filterDict = {kw("LAMP1"): "Qth_Lamp", kw("DPR_TECH"): "ECHELLE,SLIT"}
            qflatCollection = self.inputFrames.filter(**filterDict)
            filterDict = {kw("DPR_TYPE"): "DOME,FLAT", kw("DPR_TECH"): "ECHELLE,SLIT"}
            domeflatCollection = self.inputFrames.filter(**filterDict)
        else:
            filterDict = {kw("DPR_TYPE"): "LAMP,DFLAT", kw("DPR_TECH"): "ECHELLE,SLIT"}
            dflatCollection = self.inputFrames.filter(**filterDict)
            filterDict = {kw("DPR_TYPE"): "LAMP,QFLAT", kw("DPR_TECH"): "ECHELLE,SLIT"}
            qflatCollection = self.inputFrames.filter(**filterDict)

        if (
            len(flatCollection.files) == 0
            and len(dflatCollection.files) == 0
            and len(qflatCollection.files) == 0
            and len(domeflatCollection.files) == 0
        ):
            raise FileNotFoundError("The mflat recipe needs flat-frames as input, none found")

        self.calibratedFlatFiles = flatCollection.files

        self.calibratedFlatFiles = []
        self.calibratedFlatFiles[:] = [os.path.basename(l).replace("_pre", "") for l in flatCollection.files]
        self.dFlatFiles = []
        self.dFlatFiles[:] = [os.path.basename(l).replace("_pre", "") for l in dflatCollection.files]
        self.qFlatFiles = []
        self.qFlatFiles[:] = [os.path.basename(l).replace("_pre", "") for l in qflatCollection.files]
        self.domeFlatFiles = []
        self.domeFlatFiles[:] = [os.path.basename(l).replace("_pre", "") for l in domeflatCollection.files]

        # LIST OF CCDDATA OBJECTS
        flats = [
            c
            for c in flatCollection.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]
        dflats = [
            c
            for c in dflatCollection.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]
        qflats = [
            c
            for c in qflatCollection.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]
        domeflats = [
            c
            for c in domeflatCollection.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        calibratedFlats, dcalibratedFlats, qcalibratedFlats, domecalibratedFlats = self._detrend_flat_sets(
            flats, dflats, qflats, domeflats, bias, darkCollection, kw
        )

        if 1 == 0:
            from os.path import expanduser

            home = expanduser("~")
            outDir = self.settings["workspace-root-dir"].replace("~", home)
            index = 1
            for frame in calibratedFlats:
                filePath = f"{outDir}/{index:02}_flat_{arm}_calibrated.fits"
                index += 1
                self._write(frame, filePath, overwrite=True)

        self.log.debug("completed the ``calibrate_frame_set`` method")
        return calibratedFlats, dcalibratedFlats, qcalibratedFlats, domecalibratedFlats

    def _find_master_bias(self, kw):
        """*find the master bias frame among the input frames*

        **Key Arguments:**

        - ``kw`` -- the keyword lookup

        **Return:**

        - ``bias`` -- the master bias frame, or None when there is none
        """
        # FIND THE BIAS FRAMES
        filterDict = {kw("PRO_CATG"): f"MASTER_BIAS_{self.arm.upper()}"}
        biasCollection = self.inputFrames.filter(**filterDict)
        # LIST OF CCDDATA OBJECTS
        biases = [
            c
            for c in biasCollection.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        if len(biasCollection.files) == 0:
            bias = None
            biasCollection = None
        else:
            bias = biases[0]

        return bias

    def _find_dark_collection(self, kw):
        """*find the frames to use as darks: master darks, then lamp-off frames, then raw darks*

        **Key Arguments:**

        - ``kw`` -- the keyword lookup

        **Return:**

        - ``darkCollection`` -- the collection of dark frames, or None when there are none
        """
        # FIND THE DARK FRAMES
        filterDict = {kw("PRO_CATG"): f"MASTER_DARK_{self.arm.upper()}"}
        darkCollection = self.inputFrames.filter(**filterDict)

        if len(darkCollection.files) == 0:
            if self.inst.upper() == "SOXS":
                filterDict = {kw("DPR_TYPE"): "FLAT,LAMP", kw("DPR_TECH"): "IMAGE"}
            else:
                filterDict = {kw("DPR_TYPE"): "LAMP,FLAT", kw("DPR_TECH"): "IMAGE"}
            darkCollection = self.inputFrames.filter(**filterDict)

        # FINAL ATTEMPT -- FIND RAW DARK
        if len(darkCollection.files) == 0:
            filterDict = {kw("DPR_TYPE"): "DARK", kw("DPR_TECH"): "IMAGE"}
            darkCollection = self.inputFrames.filter(**filterDict)
            if len(darkCollection.files) == 0:
                darkCollection = None

        return darkCollection

    def _detrend_flat_sets(self, flats, dflats, qflats, domeflats, bias, darkCollection, kw):
        """*subtract the master bias, or the dark nearest in time, from each set of flat frames*

        **Key Arguments:**

        - ``flats`` -- the flat frames
        - ``dflats`` -- the D-lamp flat frames
        - ``qflats`` -- the QTH-lamp flat frames
        - ``domeflats`` -- the dome flat frames
        - ``bias`` -- the master bias frame, or None
        - ``darkCollection`` -- the collection of dark frames, or None
        - ``kw`` -- the keyword lookup

        **Return:**

        - ``calibratedFlats`` -- the calibrated flat frames
        - ``dcalibratedFlats`` -- the calibrated D-lamp flat frames
        - ``qcalibratedFlats`` -- the calibrated QTH-lamp flat frames
        - ``domecalibratedFlats`` -- the calibrated dome flat frames
        """
        # IF NO DARK FRAMES EXIST - JUST A MASTER BIAS. SUBTRACT BIAS.
        calibratedFlats = []
        dcalibratedFlats = []
        qcalibratedFlats = []
        domecalibratedFlats = []
        if not darkCollection and bias:
            self.log.print("\n# SUBTRACTING MASTER BIAS FROM FRAMES")
            for flat in flats:
                calibratedFlats.append(self.detrend(inputFrame=flat, master_bias=bias, dark=None))
            for flat in dflats:
                dcalibratedFlats.append(self.detrend(inputFrame=flat, master_bias=bias, dark=None))
            for flat in qflats:
                qcalibratedFlats.append(self.detrend(inputFrame=flat, master_bias=bias, dark=None))
            for flat in domeflats:
                domecalibratedFlats.append(self.detrend(inputFrame=flat, master_bias=bias, dark=None))

        # IF DARKS EXIST - FIND CLOSEST IN TIME TO FLAT-FRAME. SUBTRACT BIAS
        # AND/OR DARK
        if darkCollection:
            darkMjds = [h[kw("MJDOBS")] for h in darkCollection.headers()]
            darks = [
                c
                for c in darkCollection.ccds(
                    ccd_kwargs={
                        "hdu_uncertainty": "ERRS",
                        "hdu_mask": "QUAL",
                        "hdu_flags": "FLAGS",
                        "key_uncertainty_type": "UTYPE",
                    }
                )
            ]
            self.log.print("\n# SUBTRACTING MASTER DARK/OFF-LAMP FROM FRAMES")
            for flat in flats:

                mjd = flat.header[kw("MJDOBS")]
                matchValue, matchIndex = nearest_neighbour(flat.header[kw("MJDOBS")], darkMjds)
                dark = darks[matchIndex]
                calibratedFlats.append(self.detrend(inputFrame=flat, master_bias=bias, dark=dark))

        return calibratedFlats, dcalibratedFlats, qcalibratedFlats, domecalibratedFlats

    def normalise_flats(self, inputFlats, orderTablePath, firstPassMasterFlat=False, lamp=""):
        """*determine the median exposure for each flat frame and normalise the flux to that level*

        **Key Arguments:**

        - ``inputFlats`` -- the input flat field frames
        - ``orderTablePath`` -- path to the order table
        - ``firstPassMasterFlat`` -- the first pass of the master flat. Default *False*
        - `lamp` -- a lamp tag for QL plots

        **Return:**

        - ``normalisedFrames`` -- the normalised flat-field frames (CCDData array)
        """
        self.log.debug("starting the ``normalise_flats`` method")

        kw = self.kw

        self._read_binning_ratios(inputFlats, orderTablePath, kw)

        window = int(self.recipeSettings["centre-order-window"] / 2)

        mask = self._order_centre_mask(inputFlats, orderTablePath, window)

        if self.debug:

            this = inputFlats[0].copy()
            this.mask = mask

            quicklook_image(
                log=self.log,
                CCDObject=this,
                stdWindow=6,
                show=False,
                ext=None,
                surfacePlot=True,
                title=f"Example input flat frame with order centre mask applied {lamp}",
            )

        if not firstPassMasterFlat:
            normalisedFrames = self._normalise_flats_first_pass(inputFlats, mask)
        else:
            normalisedFrames = self._normalise_flats_second_pass(inputFlats, mask, firstPassMasterFlat)

        # PLOT ONE OF THE NORMALISED FRAMES TO CHECK
        quicklook_image(
            log=self.log,
            CCDObject=normalisedFrames[0],
            show=False,
            ext=None,
            surfacePlot=False,
            title=f"Single normalised flat frame {lamp}",
        )

        self.log.debug("completed the ``normalise_flats`` method")
        return normalisedFrames

    def _read_binning_ratios(self, inputFlats, orderTablePath, kw):
        """*read the binning of the flat frames and of the order table, and the ratio between them*

        Sets ``self.binx``, ``self.biny``, ``self.binRatioX`` and ``self.binRatioY``.

        **Key Arguments:**

        - ``inputFlats`` -- the input flat field frames
        - ``orderTablePath`` -- path to the order table
        - ``kw`` -- the keyword lookup
        """
        from astropy.io import fits

        try:
            self.binx = inputFlats[0].header[kw("WIN_BINX")]
            self.biny = inputFlats[0].header[kw("WIN_BINY")]
        except KeyError as e:
            self.log.debug(f"normalise_flats: `self.binx = inputFlats[0].header[kw('WIN_BI...` failed, continuing: {e}")
            if self.arm.lower() == "nir":
                self.binx = 1
                self.biny = 1

        # GRAB HEADER FROM ORDER LOCATION .. WHAT IS BINNING?
        with fits.open(orderTablePath, memmap=True) as hdul:
            header = hdul[0].header
        try:
            dpBinx = header[kw("WIN_BINX")]
            dpBiny = header[kw("WIN_BINY")]
        except KeyError as e:
            self.log.debug(f"normalise_flats: `dpBinx = header[kw('WIN_BINX')]` failed, continuing: {e}")
            dpBinx = 1
            dpBiny = 1

        self.binRatioX = self.binx / dpBinx
        self.binRatioY = self.biny / dpBiny

        return

    def _order_centre_mask(self, inputFlats, orderTablePath, window):
        """*build a mask that leaves only a window around each order centre, combined with the bad-pixel mask*

        **Key Arguments:**

        - ``inputFlats`` -- the input flat field frames
        - ``orderTablePath`` -- path to the order table
        - ``window`` -- the half-width of the unmasked window around each order centre, in pixels

        **Return:**

        - ``mask`` -- boolean mask, True where a pixel is excluded
        """
        import numpy as np

        # UNPACK THE ORDER TABLE & CREATE ORDER CENTRE MASK
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=orderTablePath, binx=self.binx, biny=self.biny
        )
        mask = np.ones((inputFlats[0].data.shape[0], inputFlats[0].data.shape[1]), dtype=bool)
        axisAcoords = orderTablePixels[f"{self.axisA}coord_centre"].values
        axisBcoords = orderTablePixels[f"{self.axisB}coord"].values
        axisAcoords = axisAcoords.astype(int)
        # UPDATE THE MASK
        if self.axisA == "x":
            for x, y in zip(axisAcoords, axisBcoords):
                x_start = max(0, x - window)
                x_end = min(mask.shape[1], x + window)
                if 0 <= y < mask.shape[0] and x_start < x_end:
                    mask[y, x_start:x_end] = 0
        else:
            for y, x in zip(axisAcoords, axisBcoords):
                y_start = max(0, y - window)
                y_end = min(mask.shape[0], y + window)
                if 0 <= x < mask.shape[1] and y_start < y_end:
                    mask[y_start:y_end, x] = 0
        # COMBINE MASK WITH THE BAD PIXEL MASK
        mask = np.logical_or(mask, inputFlats[0].mask)

        return mask

    def _normalise_flats_first_pass(self, inputFlats, mask):
        """*normalise each flat frame to the sigma-clipped mean of its unmasked pixels, and record the ORDEXP QCs*

        Sets ``self.qc``.

        **Key Arguments:**

        - ``inputFlats`` -- the input flat field frames
        - ``mask`` -- boolean mask, True where a pixel is excluded

        **Return:**

        - ``normalisedFrames`` -- the normalised flat-field frames (CCDData array)
        """
        import numpy as np
        from astropy.stats import sigma_clipped_stats

        normalisedFrames = []
        self.log.print("\n# NORMALISING FLAT FRAMES TO THEIR MEAN EXPOSURE LEVEL - FIRST PASS")
        ORDEXP10list = []
        ORDEXP50list = []
        ORDEXP90list = []

        for i, frame in enumerate(inputFlats):
            nrows = frame.data.shape[0]
            chunk_size = 256  # tune to balance memory vs overhead
            sample_chunks = []
            rng = np.random.default_rng(seed=42)

            for row_start in range(0, nrows, chunk_size):
                row_end = min(row_start + chunk_size, nrows)
                chunk_data = frame.data[row_start:row_end].copy()
                chunk_mask = mask[row_start:row_end]
                chunk_data[chunk_mask] = np.nan
                valid = chunk_data.ravel()
                valid = valid[~np.isnan(valid)]
                if valid.size:
                    # # subsample to cap memory: keep at most 1000 values per chunk
                    if valid.size > 10000:
                        valid = rng.choice(valid, size=10000, replace=False)
                    sample_chunks.append(valid)
                del chunk_data, chunk_mask, valid

            all_valid = np.concatenate(sample_chunks)
            del sample_chunks

            ORDEXP10list.append(np.percentile(all_valid, 10))
            ORDEXP50list.append(np.percentile(all_valid, 50))
            ORDEXP90list.append(np.percentile(all_valid, 90))

            mean, median, std = sigma_clipped_stats(
                all_valid,
                sigma=25.0,
                stdfunc="mad_std",
                cenfunc="median",
                maxiters=3,
            )
            norm_level = mean
            del all_valid

            # Divide in-place to avoid allocating a full CCDData copy
            nframe = frame.copy()
            nframe.data /= norm_level
            if nframe.uncertainty is not None:
                nframe.uncertainty.array /= norm_level
            normalisedFrames.append(nframe)
        ORDEXP10 = np.median(ORDEXP10list)
        ORDEXP50 = np.median(ORDEXP50list)
        ORDEXP90 = np.median(ORDEXP90list)

        # if ORDEXP50 < 100:
        #     raise ValueError("FLUX IN THE INPUT FLAT FRAMES IS TOO LOW TO PROCEED. PLEASE CHECK THE RAW FRAMES")

        utcnow = utcnow_string()

        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="ORDEXP10",
            qcValue=f"{ORDEXP10:0.23f}",
            qcComment="[e-] 10th percentile inter-order flux",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            qcUnit="electrons",
            toHeader=True,
        )
        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="ORDEXP50",
            qcValue=f"{ORDEXP50:0.3f}",
            qcComment="[e-] 50th percentile inter-order flux",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            qcUnit="electrons",
            toHeader=True,
        )
        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="ORDEXP90",
            qcValue=f"{ORDEXP90:0.3f}",
            qcComment="[e-] 90th percentile inter-order flux",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow,
            qcUnit="electrons",
            toHeader=True,
        )

        return normalisedFrames

    def _normalise_flats_second_pass(self, inputFlats, mask, firstPassMasterFlat):
        """*normalise each flat frame, divided by the first-pass master flat, to its sigma-clipped unmasked mean*

        **Key Arguments:**

        - ``inputFlats`` -- the input flat field frames
        - ``mask`` -- boolean mask, True where a pixel is excluded
        - ``firstPassMasterFlat`` -- the first pass of the master flat

        **Return:**

        - ``normalisedFrames`` -- the normalised flat-field frames (CCDData array)
        """
        import numpy as np
        from astropy.stats import sigma_clipped_stats

        self.log.print("\n# NORMALISING FLAT FRAMES TO THEIR MEAN EXPOSURE LEVEL - SECOND PASS")

        # Process frames one-by-one to reduce peak memory usage
        normalisedFrames = []
        chunk_size = 256  # rows per chunk - tune to balance memory vs overhead

        for frame in inputFlats:

            nrows = frame.data.shape[0]
            # Compute median of (frame / firstPassMasterFlat) in chunks
            # to avoid allocating a full-size intermediate array
            rng = np.random.default_rng(seed=56)
            chunk_vals = []
            for row_start in range(0, nrows, chunk_size):
                row_end = min(row_start + chunk_size, nrows)
                chunk_data = frame.data[row_start:row_end] / firstPassMasterFlat.data[row_start:row_end]
                chunk_nan_mask = np.isnan(chunk_data)
                chunk_combined_mask = mask[row_start:row_end] | chunk_nan_mask
                valid = chunk_data[~chunk_combined_mask]
                if valid.size:
                    if valid.size > 10000:
                        valid = rng.choice(valid, size=10000, replace=False)
                    chunk_vals.append(valid)
                del chunk_data, chunk_nan_mask, chunk_combined_mask, valid

            if chunk_vals:
                all_valid = np.concatenate(chunk_vals)
                mean, median, std = sigma_clipped_stats(
                    all_valid,
                    sigma=25.0,
                    stdfunc="mad_std",
                    cenfunc="median",
                    maxiters=3,
                )
                norm_level = mean
                all_valid /= norm_level
                del all_valid, chunk_vals

            # Divide in-place to avoid allocating a full CCDData copy
            nframe = frame.copy()
            nframe.data /= norm_level
            if nframe.uncertainty is not None:
                nframe.uncertainty.array /= norm_level
            normalisedFrames.append(nframe)
            del norm_level

        return normalisedFrames

    def mask_low_sens_pixels(self, frame, orderTablePath, returnMedianOrderFlux=False, writeQC=True):
        """*add low-sensitivity pixels to bad-pixel mask*

        **Key Arguments:**

        - ``frame`` -- the frame to work on
        - ``orderTablePath`` -- path to the order table
        - ``returnMedianOrderFlux`` -- return a table of the median order fluxes. Default *False*.
        - ``writeQC`` -- add the QCs to the QC table?

        **Return:**

        - ``frame`` -- with BPM updated with low-sensitivity pixels
        - ``medianOrderFluxDF`` -- data-frame of the median order fluxes (if ``returnMedianOrderFlux`` is True)
        """
        self.log.debug("starting the ``mask_low_sens_pixels`` method")

        import numpy as np
        import numpy.ma as ma
        import pandas as pd
        from astropy.stats import sigma_clip

        self.log.print("\n# CLIPPING LOW-SENSITIVITY PIXELS AND SETTING INTER-ORDER AREA TO UNITY")

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log,
            orderTablePath=orderTablePath,
            binx=self.binx,
            biny=self.biny,
            prebinned=True,
        )

        # BAD PIXEL COUNT AT START
        originalBPM = np.copy(frame.mask)

        interOrderMask = np.ones_like(frame.data)
        orders = orderTablePixels["order"].values
        axisAcoords_up = orderTablePixels[f"{self.axisA}coord_edgeup"].values.round().astype(int)
        axisAcoords_low = orderTablePixels[f"{self.axisA}coord_edgelow"].values.round().astype(int)
        axisBcoords = orderTablePixels[f"{self.axisB}coord"].values
        uniqueOrders = orderTablePixels["order"].unique()

        # CALCULATE AND RETURN MEDIAN FLUXES FOR ORDERS
        if returnMedianOrderFlux:
            orderFluxes = {}
            bAxisMiddles = {}
            for o in orders:
                orderFluxes[o] = []
            for o in uniqueOrders:
                filteredDf = orderTablePixels.loc[(orderTablePixels["order"] == o)]
                bAxisMiddles[o] = int(filteredDf[f"{self.axisB}coord"].mean())
            medianFlux = []

        # UPDATE THE MASK TO ALLOW INTER-ORDER PIXELS
        for u, l, b, o in zip(axisAcoords_up, axisAcoords_low, axisBcoords, orders):
            if l < 0:
                l = 0
            if u < 0:
                u = 0
            if self.axisA == "x":
                interOrderMask[b, l:u] = 0
                if returnMedianOrderFlux and b > bAxisMiddles[o] - 3 and b < bAxisMiddles[o] + 3:
                    orderFluxes[o] = np.append(orderFluxes[o], frame.data[b, l:u])
            else:
                interOrderMask[l:u, b] = 0
                if returnMedianOrderFlux and b > bAxisMiddles[o] - 3 and b < bAxisMiddles[o] + 3:
                    orderFluxes[o] = np.append(orderFluxes[o], frame.data[b, l:u])

        # GET UNIQUE VALUES IN COLUMN
        if returnMedianOrderFlux:
            for o in uniqueOrders:
                medianFlux.append(np.median(orderFluxes[o]))

        # CONVERT TO BOOLEAN MASK AND MERGE WITH BPM
        interOrderMask = ma.make_mask(interOrderMask)
        frame.mask = (interOrderMask == 1) | (frame.mask == 1)

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(
            log=self.log,
            CCDObject=frame,
            show=False,
            ext=None,
            surfacePlot=True,
            title="Masking inter-order pixels",
        )

        beforeMask = np.copy(frame.mask)

        # SIGMA-CLIP THE LOW-SENSITIVITY PIXELS
        frameClipped = sigma_clip(
            frame,
            sigma_lower=self.recipeSettings["low-sensitivity-clipping-sigma"],
            sigma_upper=2000,
            maxiters=5,
            cenfunc="median",
            stdfunc="mad_std",
        )

        lowSensitivityPixelMask = (frameClipped.mask == 1) & (beforeMask != 1)
        lowSensPixelCount = lowSensitivityPixelMask.sum()

        if writeQC:
            utcnow = utcnow_string()
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName="N LOW SENS",
                qcValue=float(lowSensPixelCount),
                qcComment="Number of low-sensitivity pixels found in master flat",
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow,
                qcUnit="pixels",
            )
            self.log.print(f"        {lowSensPixelCount} low-sensitivity pixels added to bad-pixel mask")

        frame.mask = (lowSensitivityPixelMask == 1) | (originalBPM == 1)

        # SET INTRA-ORDER TO 1 OR ZERO
        if False:
            frame.data[interOrderMask] = 0
        elif False:
            frame.data[interOrderMask] = 1

        # PLOT MASKED FRAMES TO CHECK
        quicklook_image(
            log=self.log,
            CCDObject=frame,
            show=False,
            ext=None,
            surfacePlot=True,
            title="Low Sensitivity pixels masked and inter-order pixel set to 1",
        )

        self.log.debug("completed the ``mask_low_sens_pixels`` method")

        if returnMedianOrderFlux:
            medianOrderFluxDF = {"order": uniqueOrders, "medianFlux": medianFlux}
            medianOrderFluxDF = pd.DataFrame(medianOrderFluxDF)
            return frame, medianOrderFluxDF

        return frame

    def stitch_uv_mflats(self, medianOrderFluxDF, orderTablePath):
        """*return a master UV-VIS flat frame after slicing and stitch the UV-VIS D-Lamp and QTH-Lamp flat frames*

        **Key Arguments:**

        - ``medianOrderFluxDF`` -- data frame containing median order fluxes for D and QTH frames
        - ``orderTablePath`` -- the original order table paths from order-centre tracing

        **Return:**

        - ``stitchedFlat`` -- the stitch D and QTH-Lamp master flat frame

        **Usage:**

        ```python
        mflat = self.stitch_uv_mflats(medianOrderFluxDF)
        ```
        """
        self.log.debug("starting the ``stitch_uv_mflats`` method")

        import numpy as np

        kw = self.kw

        medianOrderFluxDF["_QLAMP_PREVIOUS"] = np.insert(medianOrderFluxDF["_QLAMP"].values[:-1], 0, 999)
        medianOrderFluxDF["scale"] = medianOrderFluxDF["_DLAMP"] / medianOrderFluxDF["_QLAMP_PREVIOUS"]
        medianOrderFluxDF["closest"] = abs(1 - medianOrderFluxDF["scale"])
        medianOrderFluxDF = medianOrderFluxDF.loc[(medianOrderFluxDF["closest"] == medianOrderFluxDF["closest"].min())]

        DQscale = medianOrderFluxDF["scale"].values[0]
        orderFlip = medianOrderFluxDF["order"].values[0]

        # SCALE D FRAME TO QTH FRAME
        if self.recipeSettings["scale-d2-to-qth"]:
            dmflatScaled = self.masterFlatSet[1].divide(DQscale)
            from soxspipe.commonutils.toolkit import quicklook_image

            quicklook_image(
                log=self.log,
                CCDObject=dmflatScaled,
                show=False,
                ext=False,
                stdWindow=9,
                surfacePlot=True,
                title="D Flat scaled to Q-Flat",
            )
        else:
            dmflatScaled = self.masterFlatSet[1]

        # UNPACK THE ORDER TABLE
        orderTableMeta, orderTablePixels, orderMetaTable = unpack_order_table(
            log=self.log,
            orderTablePath=self.orderTableSet[1],
            extend=3000,
            binx=self.binx,
            biny=self.biny,
        )

        # FIND THE LINE USED TO SLICE AND STITCH THE 2 FRAMES TOGETHER
        filteredDf = orderTablePixels.loc[(orderTablePixels["order"] == orderFlip)]
        axisAStitchCoords = filteredDf[f"{self.axisA}coord_edgeup"].values.astype(int) + 4
        axisBStitchCoords = filteredDf[f"{self.axisB}coord"].values
        stitchedFlat = self.masterFlatSet[2].copy()

        # STITCH FLAT FRAMES AND COMBINED NORMALISED FRAMES (NEEDED FOR BEST ORDER EDGE DETECTION) TOGETHER
        if self.axisA == "x":
            for x, y in zip(axisAStitchCoords, axisBStitchCoords):
                if y < stitchedFlat.data.shape[0] and x < stitchedFlat.data.shape[1]:
                    stitchedFlat.data[y, :x] = dmflatScaled.data[y, :x]
                    stitchedFlat.mask[y, :x] = dmflatScaled.mask[y, :x]
                    stitchedFlat.uncertainty.array[y, :x] = dmflatScaled.uncertainty.array[y, :x]
        else:
            for y, x in zip(axisAStitchCoords, axisBStitchCoords):
                stitchedFlat.data[y, x:] = dmflatScaled.data[y, x:]
                stitchedFlat.mask[y, x:] = dmflatScaled.mask[y, x:]
                stitchedFlat.uncertainty.array[y, x:] = dmflatScaled.uncertainty.array[y, x:]

        stitchedFlat.header[kw("DPR_TYPE")] = stitchedFlat.header[kw("DPR_TYPE")].replace(",D", ",").replace(",Q", ",")

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=stitchedFlat,
            show=False,
            ext=False,
            stdWindow=9,
            title=False,
            surfacePlot=True,
        )

        # DETECT THE ORDER EDGES AND UPDATE THE ORDER LOCATIONS TABLE
        edges = detect_order_edges(
            log=self.log,
            flatFrame=stitchedFlat,
            orderCentreTable=orderTablePath,
            settings=self.settings,
            recipeSettings=self.recipeSettings,
            qcTable=self.qc,
            productsTable=self.products,
            tag="",
            sofName=self.sofName,
            binx=self.binRatioX,
            biny=self.binRatioY,
            startNightDate=self.startNightDate,
        )
        self.products, self.qc, orderDetectionCounts = edges.get()
        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = self.products["product_label"] == "ORDER_LOC"
        orderTablePath = self.products.loc[mask]["file_path"].values[0]

        stitchedFlat = self.mask_low_sens_pixels(frame=stitchedFlat, orderTablePath=orderTablePath)

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=stitchedFlat,
            show=False,
            ext=False,
            stdWindow=3,
            title=False,
            surfacePlot=True,
        )

        self.log.debug("completed the ``stitch_uv_mflats`` method")
        return stitchedFlat

    def find_uvb_overlap_order_and_scale(self, dcalibratedFlats, qcalibratedFlats):
        """*find uvb order where both lamps produce a similar flux. This is the order at which the 2 lamp flats will be scaled and stitched together*

        **Key Arguments:**

        - ``qcalibratedFlats`` -- the QTH lamp calibration flats.
        - ``dcalibratedFlats`` -- D2 lamp calibration flats

        **Return:**

        - ``order`` -- the order number where the lamp fluxes are similar

        **Usage:**

        ```python
        overlapOrder = self.find_uvb_overlap_order_and_scale(dcalibratedFlats=dcalibratedFlats, qcalibratedFlats=qcalibratedFlats)
        ```
        """
        self.log.debug("starting the ``find_uvb_overlap_order_and_scale`` method")

        import pandas as pd

        # USE THIS METHOD TO FIND THE MEAN FLUX PER ORDER FOR BOTH LAMPS
        filterDict = {
            self.kw("PRO_CATG"): f"ORDER_TAB_{self.arm}",
            self.kw("OBJECT"): "LAMP,DORDERDEF",
        }
        orderTablePaths = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)
        if len(orderTablePaths) == 1:
            orderTablePath = orderTablePaths[0]
            thisPath = orderTablePath
        normalisedFlats, DorderMeanFluxes = self.normalise_flats(dcalibratedFlats, orderTablePath=orderTablePath)
        DorderMeanFluxes.rename(columns={"90_perc": "D2"}, inplace=True)

        filterDict = {
            self.kw("PRO_CATG"): f"ORDER_TAB_{self.arm}",
            self.kw("OBJECT"): "LAMP,QORDERDEF",
        }
        orderTablePaths = self.inputFrames.filter(**filterDict).files_filtered(include_path=True)
        if len(orderTablePaths) == 1:
            orderTablePath = orderTablePaths[0]
            thisPath = orderTablePath
        normalisedFlats, QorderMeanFluxes = self.normalise_flats(qcalibratedFlats, orderTablePath=orderTablePath)
        QorderMeanFluxes.rename(columns={"90_perc": "QTH"}, inplace=True)

        # MERGE MEAN ORDER FLUX DATAFRAMES FOR BOTH LAMPS
        bothOrderMeanFluxes = pd.merge(DorderMeanFluxes, QorderMeanFluxes, on=["order"])

        # NOW FIND THE ORDER FOR WHICH THE FLUXES ARE MOST SIMILAR IN BOTH LAMPS
        mask = bothOrderMeanFluxes["QTH"] == bothOrderMeanFluxes["D2"]
        bothOrderMeanFluxes["scale"] = bothOrderMeanFluxes["D2"] / bothOrderMeanFluxes["QTH"]
        bothOrderMeanFluxes["best_frame"] = bothOrderMeanFluxes.idxmax(axis=1)

        from tabulate import tabulate

        print(tabulate(bothOrderMeanFluxes, headers="keys", tablefmt="pretty"))

        mask = bothOrderMeanFluxes["best_frame"] == "QTH"
        orderFlip = bothOrderMeanFluxes.loc[mask]
        # THIS IS THE ORDER THAT WE USE TO RESCALE ONE FLAT TO ANOTHER
        orderFlip = orderFlip["order"].max() + 1

        print(f"THE D2 and QTH FRAMES ARE FOUND TO OVERLAP AT ORDER {orderFlip}")

        self.log.debug("completed the ``find_uvb_overlap_order_and_scale`` method")
        return orderFlip

    # use the tab-trigger below for new method
    # xt-class-method


def nearest_neighbour(singleValue, listOfValues):
    import numpy as np

    arrayOfValues = np.asarray(listOfValues)
    dist = np.square(arrayOfValues - singleValue)
    minDist = np.amin(dist)
    minIndex = np.where(dist == minDist)[0][0]
    matchValue = listOfValues[minIndex]
    return matchValue, minIndex


def print_memory_usage(pprint=False, message=""):
    if pprint:
        import humanize
        import psutil

        process = psutil.Process()
        print(humanize.naturalsize(process.memory_info().rss), message)
