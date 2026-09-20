#!/usr/bin/env python
"""
*The base recipe class which all other recipes inherit*

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""

################# GLOBAL IMPORTS ####################

import os
import sqlite3
import sys

from soxspipe.commonutils import detector_lookup, filenamer, keyword_lookup, subtract_background

# VALIDATES A COLUMN OR TABLE NAME AT THE SQL TRUST BOUNDARY BEFORE IT IS
# INTERPOLATED. A LEAF MODULE WITH NO IMPORTS OF ITS OWN, SO IT CANNOT
# DISTURB THE PACKAGE'S LOAD-BEARING IMPORT ORDER.
from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

# THE "CALLER DID NOT PASS THIS KEYWORD" SENTINEL USED BY THE add_qc AND
# add_product DELEGATORS BELOW. IMPORTED FROM `toolkit` RATHER THAN REDEFINED
# HERE SO THE DELEGATORS FORWARD THE SAME OBJECT `append_qc` AND
# `append_product` TEST FOR, AND SO A CALLER MAY PASS `OMITTED` EXPLICITLY.
from soxspipe.commonutils.toolkit import OMITTED

os.environ["TERM"] = "vt100"


# THE FORMER `imstats` CLOSURE INSIDE `qc_ron`. A MODULE-LEVEL FUNCTION RATHER
# THAN A METHOD, SO THAT THE READ-OUT-NOISE CALCULATION STAYS AS UNOVERRIDEABLE
# AS THE CLOSURE IT REPLACES.
def _image_stats(dat):
    """*report the minimum, maximum, mean and standard deviation of an image*

    **Key Arguments:**

    - ``dat`` -- the image data to report on. Masked array or numpy array.

    **Return:**

    - ``stats`` -- the minimum, maximum, mean and standard deviation of ``dat``

    **Usage:**

    ```python
    dmin, dmax, dmean, dstd = _image_stats(maskedFrameData)
    ```
    """
    return (dat.min(), dat.max(), dat.mean(), dat.std())



class base_recipe:
    """
    The base recipe class which all other recipes inherit

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- debug mode. True or False. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*


    **Usage**

    To use this base recipe to create a new `soxspipe` recipe, have a look at the code for one of the simpler recipes (e.g. `soxs_mbias`) - copy and modify the code.
    """

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=False,
        verbose=False,
        overwrite=False,
        recipeName=False,
        command=False,
        debug=False,
        turnOffMP=False,
    ):
        import random

        import matplotlib
        import pandas as pd

        log.debug("instantiating a new '__init__' object")
        self.recipeName = recipeName
        self.settings = settings
        self.debug = debug
        self.workspaceRootPath = self._absolute_path(settings["workspace-root-dir"])
        self.turnOffMP = turnOffMP

        if self.debug:
            matplotlib.use("TkAgg")
        else:
            matplotlib.use("Agg")

        self.darkDetrendWarningIssued1 = False
        self.darkDetrendWarningIssued2 = False

        if isinstance(inputFrames, str) and "_STD_" in inputFrames:

            self.recipeName = self.recipeName.replace("soxs-nod", "soxs-nod-std")
            self.recipeName = self.recipeName.replace("soxs-stare", "soxs-stare-std")
            self.recipeName = self.recipeName.replace("soxs-offset", "soxs-offset-std")

        self._resolve_product_path(log, inputFrames, verbose, overwrite)

        if command:
            self.log.print(f"\nRecipe Command: {command}")

        from soxspipe.commonutils.toolkit import get_calibrations_path

        self.calibrationRootPath = get_calibrations_path(log=self.log, settings=self.settings)

        self.verbose = verbose
        # SET LATER WHEN VERIFYING FRAMES
        self.arm = None
        self.detectorParams = None
        self.dateObs = None

        # `/tmp/` HERE IS A SUBDIRECTORY OF THE USER'S OWN WORKSPACE, NOT THE
        # SYSTEM TEMPORARY DIRECTORY, AND THE RANDOM NAME ONLY HAS TO DIFFER
        # BETWEEN CONCURRENT RECIPES, NOT RESIST AN ATTACKER. THE UNSEEDED
        # DRAW ITSELF IS DY-49.
        self.outDir = self.workspaceRootPath + "/tmp/" + str(random.randint(100000, 999999))  # noqa: S108, S311

        # FIND THE CURRENT SESSION
        from os.path import expanduser

        home = expanduser("~")
        from soxspipe.commonutils import data_organiser

        do = data_organiser(
            log=self.log,
            rootDir=self.settings["workspace-root-dir"].replace("~", home),
            dbConnect=False,
        )
        self.currentSession, allSessions = do.session_list(silent=True)
        do.close()

        # INITIATE A DB CONNECTION
        self.conn = None
        self.status = None
        if not self.turnOffMP:
            self._open_session_database(home)

        # MERGE ADVANCED SETTINGS AND USER SETTINGS (USER SETTINGS OVERRIDE)
        self.settings = {**self._advanced_settings(), **self.settings}

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        self.qc, self.products = self._empty_qc_and_product_tables(pd)

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=self.recipeName,
            startNightDate=self.startNightDate,
        )

        self.generateReponseCurve = False

        return

    def _resolve_product_path(self, log, inputFrames, verbose, overwrite):
        """*predict where this recipe's product will be written, and refuse to run over one that is already there*

        **Key Arguments:**

        - ``log`` -- the logger the recipe was handed
        - ``inputFrames`` -- the input frames the recipe was handed
        - ``verbose`` -- also print the refusal to the terminal
        - ``overwrite`` -- run again over an existing product or error log

        Sets ``self.sofName``, ``self.productPath`` and ``self.log``. A product
        or an error log already in place raises ``FileExistsError`` unless
        ``overwrite`` is set.
        """
        from soxspipe.commonutils import toolkit

        # CHECK IF PRODUCT ALREADY EXISTS
        if inputFrames and not isinstance(inputFrames, list) and inputFrames.split(".")[-1].lower() == "sof":
            self.sofName = os.path.basename(inputFrames).replace(".sof", "")
            self.productPath, self.startNightDate = toolkit.predict_product_path(inputFrames, self.recipeName)
            errorLog = os.path.splitext(self.productPath)[0] + "_ERROR.log"
            basename = os.path.basename(self.productPath)

            if os.path.exists(errorLog) and not overwrite:
                if verbose:
                    print(
                        f"This recipe previously failed (see `{basename}`). To rerun the recipe, run the recipe command with the overwrite flag (-x)."
                    )
                raise FileExistsError(
                    f"This recipe previously failed (see `{basename}`). To rerun the recipe, run the recipe command with the overwrite flag (-x)."
                )
            if os.path.exists(self.productPath) and not overwrite:
                basename = os.path.basename(self.productPath)
                if verbose:
                    print(
                        f"The product of this recipe already exists: `{basename}`. To overwrite this product, rerun the pipeline command with the overwrite flag (-x)."
                    )
                raise FileExistsError(
                    f"The product of this recipe already exists: `{basename}`. To overwrite this product, rerun the pipeline command with the overwrite flag (-x)."
                )

            self.log = toolkit.add_recipe_logger(log, self.productPath)
        else:
            self.sofName = False
            self.productPath = False
            self.log = log

        return

    def _open_session_database(self, home):
        """*open the workspace database and mark this recipe as failed until it completes*

        **Key Arguments:**

        - ``home`` -- the user's home directory, used to expand the workspace path

        Sets ``self.conn`` and ``self.status``. The stored status is set to
        'fail' up front, so a recipe that crashes leaves a failure behind
        rather than its previous result.
        """
        import sqlite3 as sql

        if self.currentSession and self.sofName:
            self.sessionDb = self.settings["workspace-root-dir"].replace("~", home) + "/soxspipe.db"

            def dict_factory(cursor, row):
                d = {}
                for idx, col in enumerate(cursor.description):
                    d[col[0]] = row[idx]
                return d

            self.conn = sql.connect(
                self.sessionDb,
                check_same_thread=False,
                timeout=300,
                autocommit=True,
            )
            c = self.conn.cursor()
            c.execute("PRAGMA busy_timeout = 100000")
            c.execute("PRAGMA synchronous = OFF")

            self.conn.row_factory = dict_factory

        # SET RECIPE TO 'FAIL' AND SWITCH TO 'PASS' ONLY IF RECIPE COMPLETES
        if self.conn:
            c = self.conn.cursor()
            # THE SESSION STATUS COLUMN NAME IS A COLUMN NAME, WHICH SQLITE
            # CANNOT PARAMETERISE, AND IT COMES FROM THE WORKSPACE DATABASE
            # RATHER THAN FROM USER INPUT. THE SESSION NAME IS NEVER USED AS A
            # BARE IDENTIFIER -- IT IS VALIDATED ONLY AFTER BEING COMPOSED WITH
            # THE LITERAL `status_` PREFIX, SO A DIGIT-LEADING SESSION NAME
            # (THE DEFAULT `%Y%m%dt%H%M%S` SHAPE) STILL PASSES. THE SOF NAME IS
            # BOUND RATHER THAN INTERPOLATED.
            statusColumn = validate_sql_identifier(f"status_{self.currentSession}", "session status column")
            sqlQuery = f"select {statusColumn} as status from product_frames where sof = ?"  # noqa: S608
            c.execute(sqlQuery, (f"{self.sofName}.sof",))
            try:
                self.status = c.fetchone()["status"]
                sqlQuery = f"update product_frames set {statusColumn} = 'fail' where sof = ?"  # noqa: S608
                c.execute(sqlQuery, (f"{self.sofName}.sof",))
            except (sqlite3.Error, TypeError) as e:
                self.log.warning(f"__init__: `self.status = c.fetchone()['status']` failed, continuing: {e}")
                self.status = None

            c.close()

        return

    def _advanced_settings(self):
        """*read the advanced settings shipped with the package*

        **Return:**

        - ``advs`` -- the advanced settings, or an empty dictionary if the file cannot be found
        """
        import yaml

        # COLLECT ADVANCED SETTINGS IF AVAILABLE
        parentDirectory = os.path.dirname(__file__)
        advs = parentDirectory + "/advanced_settings.yaml"
        level = 0
        exists = False
        count = 1
        while not exists and len(advs) and count < 10:
            count += 1
            level -= 1
            exists = os.path.exists(advs)
            if not exists:
                advs = "/".join(parentDirectory.split("/")[:level]) + "/advanced_settings.yaml"
        if not exists:
            advs = {}
        else:
            with open(advs) as stream:
                advs = yaml.safe_load(stream)

        return advs

    def _empty_qc_and_product_tables(self, pd):
        """*build the empty tables this recipe collects its QCs and products in*

        **Key Arguments:**

        - ``pd`` -- the pandas module, imported by the caller so that the import happens where it always did

        **Return:**

        - ``qc`` -- the empty quality-control table
        - ``products`` -- the empty product table

        The column order set here is the order every row appended later takes.
        """
        qc = pd.DataFrame(
            {
                "soxspipe_recipe": [],
                "qc_name": [],
                "qc_value": [],
                "qc_unit": [],
                "qc_order": [],
                "qc_comment": [],
                "obs_date_utc": [],
                "reduction_date_utc": [],
                "to_header": [],
            }
        )
        products = pd.DataFrame(
            {
                "soxspipe_recipe": [],
                "product_label": [],
                "file_name": [],
                "file_type": [],
                "obs_date_utc": [],
                "reduction_date_utc": [],
                "file_path": [],
                "label": [],
            }
        )

        return qc, products

    def _prepare_single_frame(self, frame, save=False):
        """*prepare a single raw frame by converting pixel data from ADU to electrons and adding mask and uncertainty extensions*

        **Key Arguments:**

        - ``frame`` -- the path to the frame to prepare, of a CCDData object
        - ``save`` -- save the prepared frame to disk. Default: False

        **Return:**

        - ``frame`` -- the prepared frame with mask and uncertainty extensions (CCDData object)

        :::{todo}
            - write a command-line tool for this method
        :::
        """
        self.log.debug("starting the ``_prepare_single_frame`` method")

        import logging
        import warnings

        import ccdproc
        from astropy import units as u
        from astropy.nddata import CCDData

        from soxspipe.commonutils import toolkit

        warnings.filterwarnings(action="ignore")
        logging.captureWarnings(True)

        kw = self.kw
        dp = self.detectorParams

        # STORE FILEPATH FOR LATER USE
        filepath = frame

        # CONVERT FILEPATH TO CCDDATA OBJECT
        if isinstance(frame, str):
            # CONVERT RELATIVE TO ABSOLUTE PATHS
            frame = self._absolute_path(frame)
            # OPEN THE RAW FRAME - MASK AND UNCERT TO BE POPULATED LATER
            try:
                frame = CCDData.read(
                    frame,
                    hdu=0,
                    unit=u.adu,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )
            except TypeError as e:
                if "buffer is too small" in str(e):
                    self.log.warning(
                        f"Buffer is too small for frame {filepath}. The frame is likely corrupted and will not be used in the reduction.\n"
                    )
                    return None
                self.log.info(f"{filepath} is a FITS Binary Table")
                return filepath

        # CHECK THE NUMBER OF EXTENSIONS IS ONLY 1 AND "SXSPRE" DOES NOT
        # EXIST. i.e. THIS IS A RAW UNTOUCHED FRAME
        if len(frame.to_hdu()) > 1 or "SXSPRE" in frame.header:
            return filepath

        # MANIPULATE XSH DATA
        frame = self.xsh2soxs(frame)
        frame = self._trim_frame(frame)

        # CORRECT FOR GAIN - CONVERT DATA FROM ADU TO ELECTRONS
        frame = ccdproc.gain_correct(frame, dp["gain"], add_keyword=None)
        toolkit.frame_to_32(frame)

        frame = self._add_uncertainty_map(frame)

        frame.mask = self._bad_pixel_mask(frame)

        # THE NESTED `if` IS KEPT SO THE TWO CONDITIONS STAY SEPARATELY
        # COMMENTED, AND BECAUSE COLLAPSING IT IS COSMETIC WORK BELONGING TO
        # DY-88.
        if self.recipeName in ["soxs-nod-std", "soxs-stare-std", "soxs-offset-std"] and self.recipeSettings["use_flat"]:  # noqa: SIM102
            # OBJECT/STANDARD FRAMES
            if frame.meta[kw("DPR_TYPE")] == "STD,FLUX" or "STD_stare" in frame.meta[kw("OBS_NAME")]:
                # ASSUMING WE HAVE ONLY STANDARD A-B CYCLES AND NOT JITTER.
                self.generateReponseCurve = True

        frame = self._clean_cosmic_rays(frame)

        filePath = self._write_prepared_frame(frame, filepath, save)

        ## KEYWORDS FOR LATER QCs
        if self.inst == "SOXS":
            self.cptemp = frame.header[kw("CP_TEMP_C")]
            if self.arm == "NIR":
                self.detectorTemp = frame.header[kw("NIR_TEMP_K")]
            elif self.arm == "VIS":
                self.detectorTemp = frame.header[kw("VIS_TEMP_C")]
            else:
                self.detectorTemp = None

        self.log.debug("completed the ``_prepare_single_frame`` method")
        return filePath

    def _add_uncertainty_map(self, frame):
        """*add the uncertainty extension to a gain-corrected frame*

        **Key Arguments:**

        - ``frame`` -- the gain-corrected frame (CCDData object)

        **Return:**

        - ``frame`` -- the frame with its uncertainty extension populated

        A bias frame's uncertainty is the readnoise alone. Every other frame
        also carries the photon noise of its own counts.
        """
        import ccdproc
        import numpy as np

        from soxspipe.commonutils import toolkit

        kw = self.kw
        dp = self.detectorParams

        # GENERATE UNCERTAINTY MAP AS EXTENSION
        if frame.header[kw("DPR_TYPE")] == "BIAS":
            # ERROR IS ONLY FROM READNOISE FOR BIAS FRAMES
            errorMap = np.ones_like(frame.data).astype(np.float32) * dp["ron"].astype(np.float32)
            # errorMap = StdDevUncertainty(errorMap)
            frame.uncertainty = errorMap.astype(np.float32)
        else:
            # GENERATE UNCERTAINTY MAP AS EXTENSION
            frame = ccdproc.create_deviation(frame, readnoise=dp["ron"], disregard_nan=True, add_keyword=None)
        toolkit.frame_to_32(frame)

        return frame

    def _bad_pixel_mask(self, frame):
        """*read the bad-pixel bitmap matching the frame's binning and flatten it to a boolean mask*

        **Key Arguments:**

        - ``frame`` -- the frame the mask is built for (CCDData object)

        **Return:**

        - ``boolMask`` -- False where the pixel is good, True where it is bad

        A missing bitmap is written out as an all-good map before the frame is failed.
        """
        import numpy as np
        from astropy import units as u
        from astropy.nddata import CCDData

        from soxspipe.commonutils import toolkit

        kw = self.kw
        dp = self.detectorParams

        # FIND THE APPROPRIATE BAD-PIXEL BITMAP AND APPEND AS 'FLAG' EXTENSION
        # NOTE FLAGS NOT YET SUPPORTED BY CCDPROC THIS THIS WON'T GET SAVED OUT
        # AS AN EXTENSION
        arm = self.arm
        if arm != "NIR" and kw("WIN_BINX") in frame.header:
            binx = int(frame.header[kw("WIN_BINX")])
            biny = int(frame.header[kw("WIN_BINY")])
        else:
            binx = 1
            biny = 1

        bitMapPath = self.calibrationRootPath + "/" + dp["bad-pixel map"][f"{binx}x{biny}"]

        if not os.path.exists(bitMapPath):
            message = f"the path to the bitMapPath {bitMapPath} does not exist on this machine"

            if True:
                # CREATE A DUMMY BAD-PIXEL MAP
                import numpy as np
                from astropy.nddata import CCDData

                frame = CCDData(np.full_like(frame.data, 0), unit="adu")
                toolkit.frame_to_32(frame)
                # WRITE CCDDATA OBJECT TO FILE
                HDUList = frame.to_hdu()
                HDUList.verify("fix")
                HDUList.writeto(bitMapPath, output_verify="exception", overwrite=True, checksum=True)

            self.log.critical(message)
            raise OSError(message)

        bitMap = CCDData.read(bitMapPath, hdu=0, unit=u.dimensionless_unscaled)

        # frame.flags = bitMap.data

        # FLATTEN BAD-PIXEL BITMAP TO BOOLEAN FALSE (GOOD) OR TRUE (BAD) AND
        # APPEND AS 'UNCERT' EXTENSION
        boolMask = bitMap.data.astype(bool).data
        try:
            # FAILS IN PYTHON 2.7 AS BOOLMASK IS A BUFFER - NEED TO CONVERT TO
            # 2D ARRAY
            boolMask.shape  # noqa: B018

        except AttributeError as e:
            self.log.debug(f"_prepare_single_frame: `boolMask.shape` failed, continuing: {e}")
            arr = np.frombuffer(boolMask, dtype=np.uint8)
            arr.shape = frame.data.shape
            boolMask = arr

        return boolMask

    def _clean_cosmic_rays(self, frame):
        """*flag cosmic rays in the frame, if the recipe settings ask for it*

        **Key Arguments:**

        - ``frame`` -- the frame to clean (CCDData object)

        **Return:**

        - ``frame`` -- the cleaned frame, or the frame unchanged when cleaning is switched off
        """
        from soxspipe.commonutils import toolkit

        if (
            "use_lacosmic" in self.recipeSettings
            and self.recipeSettings["use_lacosmic"]
            and "lacosmic-sigma" in self.recipeSettings[self.arm.lower()]
            and self.recipeSettings[self.arm.lower()]["lacosmic-sigma"]
        ):
            # oldCount = frame.mask.sum()
            # oldMask = frame.mask.copy()

            from ccdproc import cosmicray_lacosmic

            frame = cosmicray_lacosmic(
                frame,
                sigclip=self.recipeSettings[self.arm.lower()]["lacosmic-sigma"],
                sigfrac=0.3,
                objlim=5.0,
                gain_apply=False,
                niter=2,
                verbose=False,
                cleantype="meanmask",
            )
            toolkit.frame_to_32(frame)
            # newCount = self.skySubtractedFrame.mask.sum()
            # self.laCosmicClippedCount = newCount - oldCount
            # frame.mask = oldMask | frame.mask
            from soxspipe.commonutils.toolkit import quicklook_image

            quicklook_image(
                log=self.log,
                CCDObject=frame,
                show=self.debug,
                ext=False,
                stdWindow=3,
                title="L.A.Cosmic cleaned (red = masked)",
                surfacePlot=False,
                settings=self.settings,
                skylines=False,
            )

        return frame

    def _write_prepared_frame(self, frame, filepath, save):
        """*stamp the prepared frame with its preparation time and write it to disk*

        **Key Arguments:**

        - ``frame`` -- the prepared frame (CCDData object)
        - ``filepath`` -- the path of the raw frame the prepared frame came from
        - ``save`` -- write to the workspace root instead of the recipe's scratch directory

        **Return:**

        - ``filePath`` -- the path the prepared frame was written to
        """
        from soxspipe.commonutils import toolkit

        if save:
            outDir = self.workspaceRootPath
        else:
            outDir = self.outDir

        # INJECT THE PRE KEYWORD
        frame.header["SXSPRE"] = (
            toolkit.utcnow_string(microseconds=True),
            "UTC timestamp",
        )

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outDir):
            try:
                os.makedirs(outDir)
            except OSError as e:
                self.log.debug(f"_prepare_single_frame: `os.makedirs(outDir)` failed, continuing: {e}")
        # CONVERT CCDData TO FITS HDU (INCLUDING HEADER) AND SAVE WITH PRE TAG
        # PREPENDED TO FILENAME
        basename = os.path.basename(filepath)
        filenameNoExtension = os.path.splitext(basename)[0]
        extension = os.path.splitext(basename)[1]
        filePath = outDir + "/" + filenameNoExtension + "_pre" + extension

        # SAVE TO DISK
        self._write(
            frame=frame,
            filedir=outDir,
            filename=filenameNoExtension + "_pre" + extension,
            overwrite=True,
            product=False,
        )

        return filePath

    def _absolute_path(self, path):
        """*convert paths from home directories to absolute paths*

        **Key Arguments:**

        - ``path`` -- path possibly relative to home directory

        **Return:**

        - ``absolutePath`` -- absolute path

        **Usage**

        ```python
        myPath = self._absolute_path(myPath)
        ```
        """

        from os.path import expanduser

        home = expanduser("~")
        if path[0] == "~":
            path = home + "/" + path[1:]

        return path.replace("//", "/")

    def prepare_frames(self, save=False):
        """*prepare raw frames by converting pixel data from ADU to electrons and adding mask and uncertainty extensions*

        **Key Arguments:**

        - ``save`` -- save out the prepared frame to the intermediate products directory. Default False.

        **Return:**

        - ``preframes`` -- the new image collection containing the prepared frames

        **Usage**

        Usually called within a recipe class once the input frames have been selected and verified (see `soxs_mbias` code for example):

        ```python
        self.inputFrames = self.prepare_frames(
            save=self.settings["save-intermediate-products"])
        ```
        """
        self.log.debug("starting the ``prepare_frames`` method")

        import numpy as np

        from soxspipe.commonutils.set_of_files import set_of_files

        kw = self.kw

        filepaths = self.inputFrames.files_filtered(include_path=True)

        frameCount = len(filepaths)

        if "use_lacosmic" in self.recipeSettings and self.recipeSettings["use_lacosmic"]:
            laC = ", RUNNING L.A.COSMIC"
        else:
            laC = ""

        self.log.print(
            f"\n# PREPARING {frameCount} RAW FRAMES - TRIMMING OVERSCAN, CONVERTING TO ELECTRON COUNTS, GENERATING UNCERTAINTY MAPS{laC} AND APPENDING DEFAULT BAD-PIXEL MASK"
        )
        preframes = []
        preframes[:] = [self._prepare_single_frame(frame=frame, save=save) for frame in filepaths]
        preframes = [f for f in preframes if f is not None]

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=preframes,
            recipeName=self.recipeName,
            verbose=self.verbose,
        )
        preframes, supplementaryInput = sof.get()
        preframes.sort([kw("MJDOBS")])

        self.log.print("# PREPARED FRAMES - SUMMARY")

        slitname = kw(f"SLIT_{self.arm}".upper())
        try:
            preframes.summary["SLIT"] = preframes.summary[slitname]
        except KeyError as e:
            self.log.debug(f"prepare_frames: `preframes.summary['SLIT'] = preframes.summar...` failed, continuing: {e}")

        preframes.summary["LAMP"] = "------------"
        columns = preframes.summary.colnames

        for i in range(7):
            thisLamp = kw(f"LAMP{i+1}")
            # FIRST FIND THE NAME OF THE LAMP
            newLamp = preframes.summary[thisLamp][np.where(preframes.summary[thisLamp].filled(999) != 999)]
            if len(newLamp):
                newLamp = newLamp[0]
                newLamp = newLamp.replace("_Lamp", "")
                newLamp = (
                    newLamp.replace("Argo", "Ar").replace("Neon", "Ne").replace("Merc", "Hg").replace("Xeno", "Xe")
                )

                updatedList = list(
                    preframes.summary["LAMP"][np.where(preframes.summary[thisLamp].filled(999) != 999)].data
                )
                updatedList[:] = [u.replace("-", "") + newLamp for u in updatedList]
                preframes.summary["LAMP"][np.where(preframes.summary[thisLamp].filled(999) != 999)] = updatedList
            columns.remove(thisLamp)

        preframes.summary["LAMP"][np.where(preframes.summary["LAMP"] == "------------")] = "--"

        try:
            columns.remove(kw("SLIT_NIR"))
        except ValueError as e:
            self.log.debug(f"prepare_frames: `columns.remove(kw('SLIT_NIR'))` failed, continuing: {e}")
        try:
            columns.remove(kw("SLIT_VIS"))
        except ValueError as e:
            self.log.debug(f"prepare_frames: `columns.remove(kw('SLIT_VIS'))` failed, continuing: {e}")
        try:
            columns.remove(kw("SLIT_UVB"))
        except ValueError as e:
            self.log.debug(f"prepare_frames: `columns.remove(kw('SLIT_UVB'))` failed, continuing: {e}")

        if "filename" in columns:
            # columns.remove("file")
            columns.remove("filename")
            columns = ["filename"] + columns

        # MAKE A COPY TO RENAME COLUMNS
        this = preframes.summary.copy()

        newColumns = []
        for c in columns:
            newC = c
            if "ESO SEQ " in c:
                newC = c.replace("ESO SEQ ", "")
                this.rename_column(c, newC)
            if "ESO DPR " in c:
                newC = c.replace("ESO DPR ", "")
                this.rename_column(c, newC)
            newColumns.append(newC)

        cleanColumns = newColumns.copy()
        cleanColumns.remove("file")

        self.log.print(this[cleanColumns])
        self.log.print("\n")

        self.rawFrames = this[newColumns]

        # SORT RECIPE AND ARM SETTINGS
        self.recipeSettings = self.get_recipe_settings()

        self.log.debug("completed the ``prepare_frames`` method")
        return preframes

    def _verify_input_frames_basics(self):
        """*the basic verifications that needs done for all recipes*

        **Return:**

        - None

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``_verify_input_frames_basics`` method")

        kw = self.kw

        # CHECK WE ACTUALLY HAVE IMAGES
        if not len(self.inputFrames.files_filtered(include_path=True)):
            self._report_verification_error()
            raise FileNotFoundError("No image frames where passed to the recipe")

        arm = self.inputFrames.values(keyword=kw("SEQ_ARM"), unique=True)

        # SORT RECIPE AND ARM SETTINGS
        self.recipeSettings = self.get_recipe_settings()

        self._verify_single_instrument_and_arm(arm)

        # CREATE DETECTOR LOOKUP DICTIONARY - SOME VALUES CAN BE OVERWRITTEN
        # WITH WHAT IS FOUND HERE IN FITS HEADERS
        self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        # BE CAREFUL WHAT TO MATCH WHEN INSPECTING BINNING ... DON'T BE TO STRICT
        matches = (
            (self.inputFrames.summary[kw("PRO_CATG")] != f"DISP_TAB_{self.arm}".upper())
            & (self.inputFrames.summary[kw("PRO_CATG")] != f"ORDER_TAB_{self.arm}".upper())
            & (self.inputFrames.summary[kw("PRO_CATG")] != f"DISP_IMAGE_{self.arm}".upper())
            & (self.inputFrames.summary[kw("PRO_CATG")] != f"RESP_TAB_{self.arm}".upper())
        )
        binningMatch = self.inputFrames.summary[matches]

        self._verify_single_binning(binningMatch)
        self._verify_single_read_speed(binningMatch)
        self._verify_single_gain()
        self._verify_single_slit_width()
        self._verify_single_readnoise()

        imageTypes, imageTech, imageCat = self._collect_frame_classifications()

        self.log.debug("completed the ``_verify_input_frames_basics`` method")
        return imageTypes, imageTech, imageCat

    def _report_verification_error(self, showSummary=False, trailingNewlines=False):
        """*print the frame-verification error banner, and optionally the frame summary*

        **Key Arguments:**

        - ``showSummary`` -- also print the input-frame summary table. Default *False*
        - ``trailingNewlines`` -- print two blank lines after the summary. Default *False*

        """
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - **ERROR**\n")
        if showSummary:
            self.log.print(self.inputFrames.summary)
        if trailingNewlines:
            self.log.print("\n\n")

        return

    def _verify_single_instrument_and_arm(self, arm):
        """*record the instrument and the arm the input frames were taken with*

        **Key Arguments:**

        - ``arm`` -- the unique arm values read from the input frames

        A mix of arms raises an exception; otherwise ``self.inst`` and ``self.arm`` are set.
        """
        from contextlib import suppress

        kw = self.kw

        inst = self.inputFrames.values(kw("INSTRUME"), unique=True)
        with suppress(ValueError):
            inst.remove(None)
        self.inst = inst[0]

        # MIXED INPUT ARMS ARE BAD
        if None in arm:
            arm.remove(None)
        if len(arm) > 1:
            arms = " and ".join(arm)
            self._report_verification_error(showSummary=True)
            # THIS INTERPOLATION IS A DEFECT, NOT A STYLE CHOICE: `imageTypes`
            # DOES NOT EXIST YET, SO IT RAISES `KeyError`. REWRITING IT IS
            # DY-89, AND THE CHARACTERIZATION TEST PINS TODAY'S BEHAVIOUR.
            raise TypeError("Input frames are a mix of %(imageTypes)s" % locals())  # noqa: UP031
        self.arm = arm[0]

        return

    def _verify_single_binning(self, binningMatch):
        """*record the detector binning the input frames were read out with*

        **Key Arguments:**

        - ``binningMatch`` -- the input-frame summary with the table products filtered out

        A mix of binnings raises an exception; otherwise ``self.detectorParams["binning"]`` is set.
        """
        import numpy as np

        kw = self.kw

        # MIXED BINNING IS BAD
        if self.arm == "NIR":
            # NIR ARRAY NEVER BINNED
            cdelt1 = [1]
            cdelt2 = [1]
        else:
            cdelt1 = np.unique(binningMatch[kw("CDELT1")].data)
            cdelt2 = np.unique(binningMatch[kw("CDELT2")].data)
            try:
                cdelt1.remove(None)
                cdelt2.remove(None)
            except (AttributeError, ValueError) as e:
                self.log.debug(f"_verify_input_frames_basics: `cdelt1.remove(None)` failed, continuing: {e}")

        if len(cdelt1) > 1 or len(cdelt2) > 1:
            self._report_verification_error()
            raise TypeError("Input frames are a mix of binnings" % locals())  # noqa: F507, UP031

        if cdelt1[0] and cdelt2[0]:
            self.detectorParams["binning"] = [int(cdelt2[0]), int(cdelt1[0])]

        return

    def _verify_single_read_speed(self, binningMatch):
        """*check the input frames were all read out at the same speed*

        **Key Arguments:**

        - ``binningMatch`` -- the input-frame summary with the table products filtered out

        A mix of readout speeds raises an exception.
        """
        import numpy as np

        kw = self.kw

        # MIXED READOUT SPEEDS IS BAD
        readSpeed = np.unique(binningMatch[kw("DET_READ_SPEED")].data)

        if len(readSpeed) > 1:
            self._report_verification_error(showSummary=True, trailingNewlines=True)
            raise TypeError(f"Input frames are a mix of readout speeds. {readSpeed}" % locals())

        return

    def _verify_single_gain(self):
        """*record the detector gain the input frames were read out with*

        A mix of gains raises an exception; otherwise ``self.detectorParams["gain"]`` is set.
        """
        from contextlib import suppress

        import numpy as np
        from astropy import units as u

        kw = self.kw

        # MIXED GAIN SPEEDS IS BAD
        # HIERARCH ESO DET OUT1 CONAD - Electrons/ADU
        # CONAD IS REALLY GAIN AND HAS UNIT OF Electrons/ADU
        if self.settings["instrument"] == "xsh":
            gain = self.inputFrames.values(keyword=kw("CONAD"), unique=True)
        else:
            mask = self.inputFrames.summary[kw("PRO_CATG")] != f"RESP_TAB_{self.arm}"
            checkFrames = self.inputFrames.summary[mask]
            # REMOVE MASKED/EMPTY VALUES (WORKS FOR ASTROPY TABLE COLUMNS TOO)
            conad = np.ma.asarray(checkFrames[kw("CONAD")])
            gainVals = np.ma.asarray(checkFrames[kw("GAIN")])
            valid = (~np.ma.getmaskarray(conad)) & (~np.ma.getmaskarray(gainVals))
            a = conad[valid].tolist()
            b = gainVals[valid].tolist()
            gain = [max(x, y) for x, y in zip(a, b) if x is not None and y is not None]
            gain = list(set(gain))

        with suppress(ValueError):
            gain.remove(None)

        if len(gain) > 1:
            self._report_verification_error(showSummary=True)
            # gain = np.unique(gain)
            raise TypeError(f"Input frames are a mix of gain {gain}" % locals())
        if len(gain) and gain[0]:
            # UVB & VIS
            self.detectorParams["gain"] = gain[0] * u.electron / u.adu
        else:
            # NIR
            self.log.print("\n\tGain is being read from the detector parameter file (not the FITS header)\n")
            self.detectorParams["gain"] = self.detectorParams["gain"] * u.electron / u.adu

        return

    def _verify_single_slit_width(self):
        """*check the input science and flat frames all used the same slit*

        A mix of slit widths raises an exception, except for the SOXS NIR
        nodding, staring and offset recipes, which are allowed the 5.0 and 1.5
        arcsecond pair.
        """
        from contextlib import suppress

        kw = self.kw

        # CONVERT TO DATAFRAME AND FILTER TO CHECK SLIT WIDTHS
        filteredDf = self.inputFrames.summary.to_pandas()
        matchList = [
            "LAMP,FLAT",
            "LAMP,QFLAT",
            "LAMP,DFLAT",
            "OBJECT",
            "STD,FLUX",
            f"MASTER_FLAT_{self.arm}",
        ]
        mask = ((filteredDf[kw("DPR_TYPE")].isin(matchList)) | (filteredDf[kw("PRO_CATG")].isin(matchList))) & (
            ~filteredDf[kw("PRO_CATG")].isin([f"RESP_TAB_{self.arm}"])
        )
        filteredDf = filteredDf.loc[mask]

        # MIXED SLIT-WIDTH IS BAD
        slitWidth = list(filteredDf[kw(f"SLIT_{self.arm}")].unique())
        with suppress(ValueError):
            slitWidth.remove(None)

        # ITEMS TO REMOVE
        removeItems = ["pin", "blind"]
        for i in slitWidth:
            for j in removeItems:
                if j in str(i).lower():
                    slitWidth.remove(i)

        slitWidth = [str(s).replace("JH", "") for s in slitWidth]
        slitWidth = list(set(slitWidth))

        if len(slitWidth) > 1:
            if (
                self.inst == "SOXS"
                and self.arm == "NIR"
                and self.recipeName.replace("-std", "").replace("-obj", "") in ["soxs-nod", "soxs-stare", "soxs-offset"]
                and "SLIT5.0" in slitWidth
                and "SLIT1.5" in slitWidth
                and len(slitWidth) == 2
            ):
                pass
            else:
                self._report_verification_error(showSummary=True)
                raise TypeError(f"Input frames are a mix of slit-width ({slitWidth})" % locals())

        return

    def _verify_single_readnoise(self):
        """*record the detector readnoise the input frames were read out with*

        A mix of readnoise values raises an exception; otherwise ``self.detectorParams["ron"]`` is set.
        """
        from contextlib import suppress

        from astropy import units as u

        kw = self.kw

        # HIERARCH ESO DET OUT1 RON - Readout noise in electrons
        ron = self.inputFrames.values(keyword=kw("RON"), unique=True)
        with suppress(ValueError):
            ron.remove(None)

        # MIXED NOISE
        if len(ron) > 1:
            self._report_verification_error(showSummary=True)
            raise TypeError(f"Input frames are a mix of readnoise. {ron}" % locals())
        if len(ron) and ron[0]:
            # UVB & VIS
            self.detectorParams["ron"] = ron[0] * u.electron
        else:
            # NIR
            self.detectorParams["ron"] = self.detectorParams["ron"] * u.electron

        return

    def _collect_frame_classifications(self):
        """*collect the type, technology and category of the input frames*

        **Return:**

        - ``imageTypes`` -- the unique frame types, raw and reduced
        - ``imageTech`` -- the unique frame technologies, raw and reduced
        - ``imageCat`` -- the unique frame categories, raw and reduced
        """
        kw = self.kw

        imageTypes = self.inputFrames.values(keyword=kw("DPR_TYPE"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_TYPE"), unique=True
        )
        imageTech = self.inputFrames.values(keyword=kw("DPR_TECH"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_TECH"), unique=True
        )
        imageCat = self.inputFrames.values(keyword=kw("DPR_CATG"), unique=True) + self.inputFrames.values(
            keyword=kw("PRO_CATG"), unique=True
        )

        def clean_list(myList):
            myList = list(set(myList))
            try:
                myList.remove(None)
            except ValueError as e:
                self.log.debug(f"clean_list: `myList.remove(None)` failed, continuing: {e}")
            try:
                myList.remove("REDUCED")
            except ValueError as e:
                self.log.debug(f"clean_list: `myList.remove('REDUCED')` failed, continuing: {e}")

            return myList

        return clean_list(imageTypes), clean_list(imageTech), clean_list(imageCat)

    def clean_up(self, forceFail=False):
        """*update product status in DB and remove intermediate files once recipe is complete*

        **Key Arguments:**
        - ``forceFail`` -- force the recipe to be marked as failed in the DB

        **Usage**

        ```python
        recipe.clean_up(forceFail=True)
        ```
        """
        self.log.debug("starting the ``clean_up`` method")

        import shutil

        # FILTER QC TABLE ON qc_flag
        mask = self.qc["qc_flag"] == "fail"
        passToFail = False
        failedQcs = self.qc.loc[mask]
        failedQcs = failedQcs["qc_name"].tolist()
        if len(failedQcs):
            forceFail = True
            if self.status == "pass":
                passToFail = True

        # SET RECIPE PRODUCTS TO 'PASS'
        if self.conn:

            if not passToFail and not forceFail:
                c = self.conn.cursor()
                # THE SESSION STATUS COLUMN NAME CANNOT BE A BOUND PARAMETER. THE
                # SESSION NAME IS NEVER VALIDATED BARE -- IT IS COMPOSED WITH THE
                # LITERAL `status_` PREFIX FIRST, SO A DIGIT-LEADING SESSION NAME
                # (THE DEFAULT `%Y%m%dt%H%M%S` SHAPE) STILL PASSES THE GRAMMAR.
                # THE SOF NAME IS BOUND.
                statusColumn = validate_sql_identifier(f"status_{self.currentSession}", "session status column")
                sqlQuery = f"update product_frames set {statusColumn} = 'pass' where sof = ?"  # noqa: S608
                c.execute(sqlQuery, (f"{self.sofName}.sof",))
                c.close()

            # PREVIOUSLY FAILED RECIPE THAT HAS NOW PASSED
            if (self.status == "fail" and passToFail) or (forceFail and self.status == "pass"):
                from soxspipe.commonutils import data_organiser

                if passToFail or forceFail:
                    failure = True
                else:
                    failure = False
                do = data_organiser(log=self.log, rootDir=self.workspaceRootPath)
                do.session_refresh(failure=failure)
                do.close()

            if forceFail and isinstance(forceFail, str):
                c = self.conn.cursor()
                # BOTH THE FAILURE MESSAGE AND THE SOF NAME ARE BOUND PARAMETERS
                # RATHER THAN INTERPOLATED INTO THE SQL TEXT.
                sqlQuery = "update product_frames set error_message = ? where sof = ?"
                c.execute(sqlQuery, (forceFail, f"{self.sofName}.sof"))
                c.close()

            self.conn.close()

        del self.conn

        try:
            shutil.rmtree(self.outDir)
        except OSError as e:
            self.log.debug(f"clean_up: `shutil.rmtree(self.outDir)` failed, continuing: {e}")

        if forceFail and isinstance(forceFail, str):
            self.log.error(f"\nRecipe marked as failed in the database. {forceFail}")


        elif forceFail:
            self.log.error(
                f"\nRecipe marked as failed in the database as the following QC values are outside of the acceptable limits: {', '.join(failedQcs)}."
            )

        self.log.debug("completed the ``clean_up`` method")
        return

    def xsh2soxs(self, frame):
        """*perform some massaging of the xshooter data so it more closely resembles soxs data -  this function can be removed once code is production ready*

        **Key Arguments:**

        - ``frame`` -- the CCDDate frame to manipulate

        **Return:**

        - ``frame`` -- the manipulated soxspipe-ready frame

        **Usage:**

        ```python
        frame = self.xsh2soxs(frame)
        ```
        """
        self.log.debug("starting the ``xsh2soxs`` method")
        import numpy as np

        dp = self.detectorParams

        # NP ROTATION OF ARRAYS IS IN COUNTER-CLOCKWISE DIRECTION
        rotationIndex = int(dp["clockwise-rotation"] / 90.0)

        if rotationIndex > 0:
            frame.data = np.rot90(frame.data, rotationIndex)

        self.log.debug("completed the ``xsh2soxs`` method")
        return frame

    def _trim_frame(self, frame):
        """*return frame with pre-scan and overscan regions removed*

        **Key Arguments:**

        - ``frame`` -- the CCDData frame to be trimmed

        **Return:**

        - ``trimmed_frame`` -- the frame with its pre-scan and overscan regions removed (CCDData object)
        """
        self.log.debug("starting the ``_trim_frame`` method")

        import ccdproc

        dp = self.detectorParams

        rs, re, cs, ce = (
            dp["science-pixels"]["rows"]["start"],
            dp["science-pixels"]["rows"]["end"],
            dp["science-pixels"]["columns"]["start"],
            dp["science-pixels"]["columns"]["end"],
        )

        binning = dp["binning"]
        if binning[0] > 1:
            rs = int(rs / binning[0])
            re = int(re / binning[0])
        if binning[1] > 1:
            cs = int(cs / binning[1])
            ce = int(ce / binning[1])

        trimmed_frame = ccdproc.trim_image(frame[rs:re, cs:ce], add_keyword=None)

        self.log.debug("completed the ``_trim_frame`` method")
        return trimmed_frame

    def _write(
        self,
        frame,
        filedir,
        filename=False,
        overwrite=True,
        product=True,
        maskToZero=False,
    ):
        """*write frame to disk at the specified location*

        **Key Arguments:**

        - ``frame`` -- the frame to save to disk (CCDData object)
        - ``filedir`` -- the location to save the frame
        - ``filename`` -- the filename to save the file as. Default: **False** (standardised filename generated in code)
        - ``overwrite`` -- if a file exists at the filepath then choose to overwrite the file. Default: True
        - ``product`` -- is this a recipe product?
        - ``maskToZero`` -- set masked pixels to zero before writing to file?

        **Return:**

        - ``filepath`` -- the absolute path of the file written to disk

        **Usage:**

        Use within a recipe like so:

        ```python
        self._write(frame, filePath)
        ```
        """
        self.log.debug("starting the ``write`` method")

        from soxspipe.commonutils.phase3 import basic_header_scrubbing, sort_keywords

        # WRITE QCs TO HEADERS
        for n, v, c, h in zip(
            self.qc["qc_name"].values,
            self.qc["qc_value"].values,
            self.qc["qc_comment"].values,
            self.qc["to_header"].values,
        ):
            if h:
                frame.header[f"ESO QC {n}".upper()] = (v, c)

        if product:
            frame.header = basic_header_scrubbing(log=self.log, settings=self.settings, header=frame.header)
        frame.header = sort_keywords(log=self.log, header=frame.header)

        if not filename and self.sofName:
            filename = self.sofName + ".fits"
        if not filename:
            filename = filenamer(log=self.log, frame=frame, settings=self.settings)

        if product:
            filedir += f"/reduced/{self.startNightDate}/{self.recipeName}/"
            filedir = filedir.replace("//", "/")
            # Recursively create missing directories
            if not os.path.exists(filedir):
                try:
                    os.makedirs(filedir)
                except OSError as e:
                    self.log.debug(f"_write: `os.makedirs(filedir)` failed, continuing: {e}")

        filepath = filedir + "/" + filename

        # SET BAD-PIXELS TO 0 IN DATA FRAME
        if maskToZero:
            self.log.print(f"\nSetting {frame.mask.sum()} bad-pixels to a value of 0 while saving '{filename}'.")
            frame.data[frame.mask] = 1

        HDUList = frame.to_hdu(hdu_mask="QUAL", hdu_uncertainty="ERRS", hdu_flags=None)
        HDUList[0].name = "FLUX"
        HDUList.verify("fix")
        if product:
            HDUList.writeto(filepath, output_verify="fix+warn", overwrite=overwrite, checksum=True)
        else:
            HDUList.writeto(filepath, overwrite=overwrite, checksum=False)

        filepath = os.path.abspath(filepath)

        self.log.debug("completed the ``write`` method")
        return filepath

    def _ccds_for_stacking(self, frames):
        """*convert the frames to stack into a list of 32-bit CCDData objects*

        **Key Arguments:**

        - ``frames`` -- an ImageFileCollection of the frames to stack or a list of CCDData objects

        **Return:**

        - ``ccds`` -- the frames as a list of CCDData objects

        **Usage:**

        ```python
        ccds = self._ccds_for_stacking(frames)
        ```
        """
        from astropy import units as u

        from soxspipe.commonutils import toolkit

        # LIST OF CCDDATA OBJECTS NEEDED BY COMBINER OBJECT
        if not isinstance(frames, list):
            return [
                toolkit.frame_to_32(c)
                # c
                for c in frames.ccds(
                    ccd_kwargs={
                        "hdu_uncertainty": "ERRS",
                        "hdu_mask": "QUAL",
                        "hdu_flags": "FLAGS",
                        "key_uncertainty_type": "UTYPE",
                        "unit": u.electron,
                    }
                )
            ]
        return [toolkit.frame_to_32(c) for c in frames]

    def _clip_individual_frames(
        self,
        combiner,
        stacked_clipping_sigma,
        stacked_clipping_iterations,
        totalPixels,
    ):
        """*sigma-clip the individual input frames, replacing the combiner's masks*

        **Key Arguments:**

        - ``combiner`` -- the Combiner object holding the frames to clip. Its mask is replaced in place.
        - ``stacked_clipping_sigma`` -- the sigma level to clip the individual frames at
        - ``stacked_clipping_iterations`` -- the number of clipping iterations to run
        - ``totalPixels`` -- the pixel count of one frame, used to report the clipped fraction

        **Usage:**

        ```python
        self._clip_individual_frames(
            combiner,
            stacked_clipping_sigma=3.0,
            stacked_clipping_iterations=2,
            totalPixels=np.size(combinedMask),
        )
        ```
        """
        import numpy as np
        from astropy.stats import sigma_clip

        # GENERATE A MASK FOR EACH OF THE INDIVIDUAL INPUT FRAMES - USING
        # MEDIAN WITH MEDIAN ABSOLUTE DEVIATION (MAD) AS THE DEVIATION FUNCTION
        # THIS IS THE SUM OF BAD-PIXELS IN ALL INDIVIDUAL FRAME MASKS
        old_n_masked = combiner.data_arr.mask.sum()

        ## Reduce memory by avoiding an extra full-array copy during clipping
        combiner.data_arr.mask = sigma_clip(
            np.asarray(combiner.data_arr.data, dtype=np.float32),
            sigma_lower=stacked_clipping_sigma,
            sigma_upper=stacked_clipping_sigma,
            axis=0,
            copy=False,
            maxiters=stacked_clipping_iterations,
            cenfunc="median",
            stdfunc="mad_std",
            masked=True,
        ).mask

        # RECOUNT BAD-PIXELS NOW CLIPPING HAS RUN
        new_n_masked = combiner.data_arr.mask.sum()
        diff = new_n_masked - old_n_masked
        if self.verbose:
            percent = 100 * combiner.data_arr.mask[0].sum() / totalPixels
            self.log.print(
                f"\tClipping found {diff} more rogue pixels in the set of all input frames (~{percent:0.2}% per-frame)"
            )
        return

    def _mean_combine_frames(self, combiner, ccds, combinedMask):
        """*mean combine the clipped frames, and combine their error maps the same way*

        **Key Arguments:**

        - ``combiner`` -- the Combiner object holding the clipped frames. Its data is overwritten with the error maps.
        - ``ccds`` -- the input frames as a list of CCDData objects
        - ``combinedMask`` -- the bad-pixel mask shared by every input frame

        **Return:**

        - ``combined_frame`` -- the combined frame, with its uncertainty map attached

        **Usage:**

        ```python
        combined_frame = self._mean_combine_frames(combiner, ccds, combinedMask)
        ```
        """
        import numpy as np

        from soxspipe.commonutils import toolkit

        # GENERATE THE COMBINED MEAN
        # self.log.print("\n# MEAN COMBINING FRAMES - WITH UPDATED BAD-PIXEL MASKS")
        combined_frame = combiner.average_combine()

        # RECOMBINE THE COMBINED MASK FROM ABOVE
        combined_frame.mask = combined_frame.mask | combinedMask

        # INDIVIDUAL UPDATED MASKS (POST CLIPPING)
        new_individual_masks = combiner.data_arr.mask
        masked_values = new_individual_masks.sum(axis=0)

        # A HACK TO THE COMBINER OBJECT TO COMBINE ERROR MAPS EXACTLY AS DATA WAS COMBINED
        for i, ccd in enumerate(ccds):
            combiner.data_arr.data[i] = ccd.uncertainty.array
        combined_uncertainty = combiner.average_combine()
        combined_frame.uncertainty = combined_uncertainty.data / (np.sqrt(len(new_individual_masks) - masked_values))
        toolkit.frame_to_32(combined_frame)

        return combined_frame

    def clip_and_stack(self, frames, recipe, ignore_input_masks=False, post_stack_clipping=True):
        """*mean combine input frames after sigma-clipping outlying pixels using a median value with median absolute deviation (mad) as the deviation function*

        **Key Arguments:**

        - ``frames`` -- an ImageFileCollection of the frames to stack or a list of CCDData objects
        - ``recipe`` -- the name of recipe needed to read the correct settings from the yaml files
        - ``ignore_input_masks`` -- ignore the input masks during clip and stacking?
        - ``post_stack_clipping`` -- allow cross-plane clipping on combined frame. Clipping settings in setting file. Default *True*.

        **Return:**

        - ``combined_frame`` -- the combined master frame (with updated bad-pixel and uncertainty maps)

        **Usage:**

        This snippet can be used within the recipe code to combine individual (using bias frames as an example):

        ```python
        combined_bias_mean = self.clip_and_stack(
            frames=self.inputFrames, recipe="soxs_mbias", ignore_input_masks=False, post_stack_clipping=True)
        ```
        """
        self.log.debug("starting the ``clip_and_stack`` method")

        import numpy as np
        from astropy.stats import sigma_clip

        from soxspipe.commonutils.combiner import Combiner

        if len(frames) == 1:
            self.log.info(
                "Only 1 frame was sent to the clip and stack method. Returning the frame with no further processing."
            )
            return frames[0]
        if len(frames) == 0:
            self.log.critical("No frames were sent to the clip and stack method. Cannot proceed.")
            raise ValueError("No frames were sent to the clip and stack method.")

        arm = self.arm
        kw = self.kw
        imageType = self.imageType

        # ALLOW FOR UNDERSCORE AND HYPHENS
        recipe = recipe.replace("soxs_", "soxs-")

        # UNPACK SETTINGS
        stacked_clipping_sigma = self.recipeSettings["stacked-clipping-sigma"]
        stacked_clipping_iterations = self.recipeSettings["stacked-clipping-iterations"]
        if post_stack_clipping:
            frame_clipping_sigma = self.recipeSettings["frame-clipping-sigma"]
            frame_clipping_iterations = self.recipeSettings["frame-clipping-iterations"]
        else:
            frame_clipping_sigma = None
            frame_clipping_iterations = None

        ccds = self._ccds_for_stacking(frames)

        imageType = ccds[0].header[kw("DPR_TYPE")].replace(",", "-")
        imageTech = ccds[0].header[kw("DPR_TECH")].replace(",", "-")
        imageCat = ccds[0].header[kw("DPR_CATG")].replace(",", "-")

        self.log.print(f"\n# MEAN COMBINING {len(ccds)} {arm} {imageCat} {imageTech} {imageType} FRAMES")

        # COMBINE MASKS AND THEN RESET
        combinedMask = ccds[0].mask
        for c in ccds:
            combinedMask = c.mask & combinedMask
            if ignore_input_masks:
                c.mask[:, :] = False

        # COMBINER OBJECT WILL FIRST GENERATE MASKS FOR INDIVIDUAL IMAGES VIA
        # CLIPPING AND THEN COMBINE THE IMAGES WITH THE METHOD SELECTED. PIXEL
        # MASKED IN ALL INDIVIDUAL IMAGES ARE MASK IN THE FINAL COMBINED IMAGE
        combiner = Combiner(ccds, dtype=np.float32)

        # self.log.print(f"\n# SIGMA-CLIPPING PIXEL WITH OUTLYING VALUES IN INDIVIDUAL {imageType} FRAMES")
        # PRINT SOME INFO FOR USER
        badCount = combinedMask.sum()
        totalPixels = np.size(combinedMask)
        percent = (float(badCount) / float(totalPixels)) * 100.0
        if imageType != "BIAS":
            self.log.print(
                f"\tThe basic bad-pixel mask for the {arm} detector {imageType} frames contains {badCount} pixels ({percent:0.2}% of all pixels)"
            )

        self._clip_individual_frames(
            combiner,
            stacked_clipping_sigma=stacked_clipping_sigma,
            stacked_clipping_iterations=stacked_clipping_iterations,
            totalPixels=np.size(combinedMask),
        )

        combined_frame = self._mean_combine_frames(combiner, ccds, combinedMask)

        # MASSIVE FUDGE - NEED TO CORRECTLY WRITE THE HEADER FOR COMBINED
        # IMAGES
        combined_frame.header = ccds[0].header
        try:
            combined_frame.wcs = ccds[0].wcs
        except (AttributeError, IndexError) as e:
            self.log.debug(f"clip_and_stack: `combined_frame.wcs = ccds[0].wcs` failed, continuing: {e}")

        if post_stack_clipping:
            maskedFrame = sigma_clip(
                combined_frame.data.astype(np.float32),
                sigma=frame_clipping_sigma,
                maxiters=frame_clipping_iterations,
                cenfunc="mean",
                stdfunc="std",
            )
            # RECOMBINE THE COMBINED MASK FROM ABOVE
            combined_frame.mask = combined_frame.mask | maskedFrame.mask

        # CALCULATE NEW PIXELS ADDED TO MASK
        if imageType != "BIAS":
            newBadCount = combined_frame.mask.sum()
            diff = newBadCount - badCount
            totalPixels = np.size(combinedMask)
            percent = (float(newBadCount) / float(totalPixels)) * 100.0
            self.log.print(
                f"\t{diff} new pixels made it into the combined bad-pixel map (bad pixels now account for {percent:0.2f}% of all pixels)"  # noqa: E501
            )

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=combined_frame,
            show=self.debug,
            ext=False,
            stdWindow=3,
            title=False,
            surfacePlot=True,
        )
        self.log.debug("completed the ``clip_and_stack`` method")
        return combined_frame

    def _subtract_dark_frame(self, processedFrame, dark, inputFrame):
        """*subtract a dark frame, scaling it first when its exposure time does not match*

        **Key Arguments:**

        - ``processedFrame`` -- the frame to have the dark subtracted. CCDData object.
        - ``dark`` -- the dark frame to be subtracted. CCDData object.
        - ``inputFrame`` -- the frame as it was handed to `detrend`, reported on when the dark is
          scaled. CCDData object.

        **Return:**

        - ``processedFrame`` -- the frame with the dark subtracted. CCDData object.

        **Usage:**

        ```python
        processedFrame = self._subtract_dark_frame(processedFrame, dark, inputFrame)
        ```
        """
        import ccdproc
        from astropy import units as u

        kw = self.kw

        # DARK WITH MATCHING EXPOSURE TIME
        tolerence = 0.5
        if (int(dark.header[kw("EXPTIME")]) < int(processedFrame.header[kw("EXPTIME")]) + tolerence) and (
            int(dark.header[kw("EXPTIME")]) > int(processedFrame.header[kw("EXPTIME")]) - tolerence
        ):
            return ccdproc.subtract_dark(
                processedFrame,
                dark,
                exposure_time=kw("EXPTIME"),
                exposure_unit=u.second,
                add_keyword=None,
            )

        # THE `and False` DISABLES THIS BRANCH DELIBERATELY. REMOVING IT IS A
        # BEHAVIOUR CHANGE, NOT A LINT FIX, SO IT IS DY-88'S TO DECIDE.
        if self.inst == "SOXS" and False:  # noqa: SIM223
            if not self.darkDetrendWarningIssued2:
                self.log.warning(
                    "Dark and science/calibration frame have differing exposure-times. SOXS dark noise does not scale linearly with time. Skipping dark subtraction."  # noqa: E501
                )
                self.darkDetrendWarningIssued2 = True
            return processedFrame

        if not self.darkDetrendWarningIssued2:
            self.log.warning(
                "Dark and science/calibration frame have differing exposure-times. Scaling dark to match science/calibration frame."  # noqa: E501
            )
            self.darkDetrendWarningIssued2 = True
            self.log.print(f"Scaling the dark to the exposure time of {inputFrame.header[kw('EXPTIME')]}s")
        processedFrame = ccdproc.subtract_dark(
            processedFrame,
            dark,
            exposure_time=kw("EXPTIME"),
            exposure_unit=u.second,
            scale=True,
            add_keyword=None,
        )
        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(log=self.log, CCDObject=dark, show=True, ext="mask", stdWindow=3, title=False, surfacePlot=True)
        quicklook_image(
            log=self.log,
            CCDObject=processedFrame,
            show=True,
            ext="mask",
            stdWindow=3,
            title=False,
            surfacePlot=True,
        )
        return processedFrame

    def _subtract_scattered_light(self, processedFrame, order_table):
        """*subtract the scattered-light background from a frame, using the order edges*

        Sets `self.products` to the product table the background subtraction returns.

        **Key Arguments:**

        - ``processedFrame`` -- the frame to have the background subtracted. CCDData object.
        - ``order_table`` -- order table with order edges defined.

        **Return:**

        - ``processedFrame`` -- the background-subtracted frame. CCDData object.

        **Usage:**

        ```python
        processedFrame = self._subtract_scattered_light(processedFrame, order_table)
        ```
        """
        background = subtract_background(
            log=self.log,
            frame=processedFrame,
            sofName=self.sofName,
            recipeName=self.recipeName,
            orderTable=order_table,
            settings=self.settings,
            productsTable=self.products,
            qcTable=self.qc,
            startNightDate=self.startNightDate,
        )
        backgroundFrame, processedFrame, self.products = background.subtract()

        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=backgroundFrame,
            show=False,
            ext="data",
            stdWindow=3,
            title="Background Light",
            surfacePlot=True,
        )

        return processedFrame

    def detrend(
        self,
        inputFrame,
        master_bias=False,
        dark=False,
        master_flat=False,
        order_table=False,
    ):
        """*subtract calibration frames from an input frame*

        **Key Arguments:**

        - ``inputFrame`` -- the input frame to have calibrations subtracted. CCDData object.
        - ``master_bias`` -- the master bias frame to be subtracted. CCDData object. Default *False*.
        - ``dark`` -- a dark frame to be subtracted. CCDData object. Default *False*.
        - ``master_flat`` -- divided input frame by this master flat frame. CCDData object. Default *False*.
        - ``order_table`` -- order table with order edges defined. Used to subtract scattered light background from frames. Default *False*.

        **Return:**

        - ``calibration_subtracted_frame`` -- the input frame with the calibration frame(s) subtracted. CCDData object.

        **Usage:**

        Within a soxspipe recipe use `detrend` like so:

        ```python
        myCalibratedFrame = self.detrend(
            inputFrame=inputFrameCCDObject, master_bias=masterBiasCCDObject, dark=darkCCDObject)
        ```
        """
        self.log.debug("starting the ``detrend`` method")

        from datetime import datetime

        import ccdproc

        from soxspipe.commonutils import toolkit

        # `arm` AND `dp` ARE UNUSED, TWO OF THE MODULE'S `F841` FINDINGS. DELETING
        # THEM IS DY-88'S.
        arm = self.arm  # noqa: F841
        kw = self.kw
        dp = self.detectorParams  # noqa: F841

        if master_bias == None:
            master_bias = False
        if dark == None:
            dark = False

        # VERIFY DATA IS IN ORDER
        # EACH OF THESE IS EITHER `False` OR A CCDData FRAME, WHOSE TRUTH VALUE IS
        # AMBIGUOUS, SO NONE OF THE COMPARISONS CAN BECOME A TRUTH CHECK.
        if master_bias == False and dark == False and master_flat == False:  # noqa: E712
            raise TypeError("detrend method needs at least a master-bias frame, a dark frame or a master flat frame")
        if master_bias == False and dark != False and dark.header[kw("EXPTIME")] != inputFrame.header[kw("EXPTIME")]:
            if not self.darkDetrendWarningIssued1:
                self.log.warning("Dark and science/calibration frame have differing exposure-times.")
                self.darkDetrendWarningIssued1 = True

        processedFrame = inputFrame

        if master_bias != False:  # noqa: E712
            processedFrame = ccdproc.subtract_bias(processedFrame, master_bias, add_keyword=None)
            toolkit.frame_to_32(processedFrame)

        # `dark` IS EITHER `False` OR A CCDData FRAME, WHOSE TRUTH VALUE IS
        # AMBIGUOUS, SO THE COMPARISON TO `False` CANNOT BECOME A TRUTH CHECK.
        if dark != False:  # noqa: E712
            processedFrame = self._subtract_dark_frame(processedFrame, dark, inputFrame)
        toolkit.frame_to_32(processedFrame)

        doSubtraction = True
        if "subtract_background" in self.recipeSettings and not self.recipeSettings["subtract_background"]:
            doSubtraction = False

        # `order_table` IS EITHER `False` OR A TABLE PATH.
        if order_table != False and doSubtraction:  # noqa: E712
            processedFrame = self._subtract_scattered_light(processedFrame, order_table)

            # ASSIGNED AND NEVER READ. IT IS THE MODULE'S LAST DUPLICATION HIT AND
            # ONE OF ITS `F841` FINDINGS, BOTH OF WHICH ARE DY-88'S.
            utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")  # noqa: F841

        if master_flat != False:  # noqa: E712
            processedFrame = ccdproc.flat_correct(processedFrame, master_flat, norm_value=1.0, add_keyword=None)
            toolkit.frame_to_32(processedFrame)

        self.log.debug("completed the ``detrend`` method")
        return processedFrame

    def _qc_report_columns(self):
        """*select the QC columns to report to the terminal and to the database*

        **Return:**

        - ``columns`` -- the QC columns to print to the terminal
        - ``dbColumns`` -- the QC columns to send to the database

        **Usage:**

        ```python
        columns, dbColumns = self._qc_report_columns()
        ```
        """
        columns = list(self.qc.columns)
        columns.remove("to_header")
        columns.remove("obs_date_utc")
        columns.remove("qc_order")
        columns.remove("reduction_date_utc")
        columns.remove("soxspipe_recipe")
        columns.remove("sof_name")
        try:
            columns.remove("qc_value_min")
            columns.remove("qc_value_max")
        except ValueError as e:
            self.log.debug(f"report_output: `columns.remove('qc_value_min')` failed, continuing: {e}")
        dbColumns = list(self.qc.columns)
        dbColumns.remove("to_header")

        return columns, dbColumns

    @staticmethod
    def _format_qc_for_display(qcRows):
        """*render the numeric values of a QC table to four decimal places*

        **Key Arguments:**

        - ``qcRows`` -- the QC rows to render. Pandas DataFrame.

        **Return:**

        - ``qc_display`` -- a copy of ``qcRows`` with its numeric values rendered as strings

        **Usage:**

        ```python
        qc_display = self._format_qc_for_display(self.qc.loc[mask][columns])
        ```
        """
        # Format float values to 3 decimal places
        qc_display = qcRows.copy()
        for col in qc_display.columns:
            qc_display[col] = qc_display[col].apply(
                lambda x: (
                    f"{float(x):.4f}"
                    if isinstance(x, (int, float))
                    or (isinstance(x, str) and x.replace(".", "", 1).replace("-", "", 1).isdigit())
                    else x
                )
            )
        return qc_display

    def report_output(self, rformat="stdout"):
        """*a method to report QC values alongside intermediate and final products*

        **Key Arguments:**

        - ``rformat`` -- the format to outout reports as. Default *stdout*. [stdout|....]

        **Return:**

        - ``qc`` -- the QC dataframe, with the columns that are written to the database

        **Usage:**

        ```python
        self.report_output(rformat="stdout")
        ```
        """
        self.log.debug("starting the ``report_output`` method")

        from tabulate import tabulate

        if not self.verbose:
            # REMOVE COLUMN FROM DATA FRAME
            self.products.drop(columns=["file_path"], inplace=True)

        self.flag_poor_data()

        # REMOVE DUPLICATE ENTRIES IN COLUMN 'qc_name' AND KEEP THE LAST ENTRY
        self.qc = self.qc.drop_duplicates(subset=["qc_name", "qc_order"], keep="last")

        # SEND TO DATABASE
        self.qc["sof_name"] = self.sofName + ".sof"
        self.qc["obs_date_utc"] = self.dateObs

        # SORT BY COLUMN NAME
        self.qc.sort_values(["qc_name"], inplace=True, kind="stable")
        columns, dbColumns = self._qc_report_columns()

        # SORT BY COLUMN NAME
        self.products.sort_values(["label"], ascending=[True], inplace=True, kind="stable")

        self.products.drop_duplicates(inplace=True)

        columns2 = list(self.products.columns)
        columns2.remove("reduction_date_utc")
        columns2.remove("soxspipe_recipe")

        try:
            soxspipe_recipe = self.qc["soxspipe_recipe"].values[0].upper()
        except (IndexError, KeyError, AttributeError) as e:
            self.log.debug(f"report_output: `soxspipe_recipe = self.qc['soxspipe_recipe']....` failed, continuing: {e}")
            soxspipe_recipe = self.recipeName.upper()

        if rformat == "stdout":
            self.log.print(f"\n# {soxspipe_recipe} QC METRICS")

            mask = self.qc["qc_order"] == "-1"

            qc_display = self._format_qc_for_display(self.qc.loc[mask][columns])
            self.log.print(
                tabulate(
                    qc_display,
                    headers="keys",
                    tablefmt="pretty",
                    showindex=False,
                    stralign="right",
                )
            )
            self.log.print(f"\n# {soxspipe_recipe} RECIPE PRODUCTS & QC OUTPUTS")
            self.log.print(
                tabulate(
                    self.products[columns2],
                    headers="keys",
                    tablefmt="pretty",
                    showindex=False,
                    stralign="right",
                )
            )
            if self.conn:
                sofNames = self.qc[dbColumns]["sof_name"].values.tolist()
                # A FALSE POSITIVE: THE F-STRING INTERPOLATES ONLY `?` PLACEHOLDERS,
                # AND EVERY SOF NAME IS PASSED TO `execute` AS A BOUND PARAMETER.
                sqlQuery = f"delete from quality_control where sof_name in ({', '.join(['?']*len(sofNames))})"  # noqa: S608
                c = self.conn.cursor()
                c.execute(sqlQuery, sofNames)
                c.close()
                self._dataframe_to_sqlite(self.qc[dbColumns], "quality_control", replace=False)

        self.log.debug("completed the ``report_output`` method")
        return self.qc[dbColumns]

    def flag_poor_data(self):
        """*a method to flag data as 'poor' quality based on QC values exceeding thresholds defined in the settings file*

        **Usage:**

        ```python
        self.flag_poor_data()
        ```
        """
        self.log.debug("starting the ``flag_poor_data`` method")

        from soxspipe.commonutils import toolkit

        # ONE TIMESTAMP COVERS BOTH TEMPERATURE ROWS, AS IT DID INLINE.
        utcnow = toolkit.utcnow_string()

        if self.inst.upper() == "SOXS":
            self.add_qc(
                qcName="DETECTOR TEMP",
                qcValue=self.detectorTemp,
                qcComment="[K] temp of detector",
                qcUnit="kelvin",
                reductionDateUtc=utcnow,
                toHeader=False,
            )
            self.add_qc(
                qcName="CPATH TEMP",
                qcValue=self.cptemp,
                qcComment="[C] temp of common path",
                qcUnit="celsius",
                reductionDateUtc=utcnow,
                toHeader=False,
            )

        # FILTER DATA FRAME
        mask = self.qc["qc_order"].isna()
        self.qc.loc[mask, "qc_order"] = "-1"
        self.qc["qc_flag"] = "pass"

        if "qc-acceptable-ranges" not in self.recipeSettings:
            self.log.debug("No acceptable ranges defined in settings file. Skipping the ``flag_poor_data`` method.")
            return

        for k, v in self.recipeSettings["qc-acceptable-ranges"].items():
            matchName = k.lower().replace("-", " ").replace("_", " ")
            mask = (
                (self.qc["qc_name"].str.lower() == matchName)
                & (self.qc["qc_order"] == "-1")
                & ((self.qc["qc_value"].astype(float) <= v[0]) | (self.qc["qc_value"].astype(float) >= v[1]))
            )
            self.qc.loc[mask, "qc_flag"] = "fail"
            # ADD GUARDRAIL VALUES TO QC TABLE FOR REPORTING PURPOSES
            self.qc.loc[(self.qc["qc_name"].str.lower() == matchName), "qc_value_min"] = v[0]
            self.qc.loc[(self.qc["qc_name"].str.lower() == matchName), "qc_value_max"] = v[1]

        self.log.debug("completed the ``flag_poor_data`` method")
        return

    def _measure_raw_frame_ron(self):
        """*measure the read-out noise in a single raw frame, from the first two input frames*

        The noise is measured on the difference of the first two input frames,
        so the mask returned alongside it is the mask of that difference and is
        needed by the caller to measure the master frame against the same pixels.

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons
        - ``combinedMask`` -- the bad-pixel mask of the sigma-clipped frame difference

        **Usage:**

        ```python
        rawRon, combinedMask = self._measure_raw_frame_ron()
        ```
        """
        import math

        import numpy as np
        from astropy.stats import sigma_clip

        from soxspipe.commonutils import toolkit

        # LIST OF RAW CCDDATA OBJECTS
        # THE COMPREHENSION IS THE BASE REVISION'S. COLLAPSING IT TO `list()` IS
        # DY-88'S LINT SWEEP, NOT THIS COMMIT'S.
        ccds = [  # noqa: C416
            c
            for c in self.inputFrames.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        # SINGLE FRAME RON
        raw_one = ccds[0]
        raw_two = ccds[1]
        raw_diff = raw_one.subtract(raw_two)
        toolkit.frame_to_32(raw_diff)

        # SIGMA-CLIP THE DATA (AT HIGH LEVEL)
        masked_diff = sigma_clip(
            raw_diff.data.astype(np.float32),
            sigma_lower=10,
            sigma_upper=10,
            maxiters=2,
            cenfunc="median",
            stdfunc="mad_std",
        )
        combinedMask = raw_diff.mask | masked_diff.mask

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        raw_diff = np.ma.array(raw_diff.data, mask=combinedMask)

        dmin, dmax, dmean, dstd = _image_stats(raw_diff)

        if dstd == 0:
            message = "The raw input frames appear to be corrupted. Cannot calculate the read-out noise. Please check the raw frames."  # noqa: E501
            raise ValueError(message)

        # ACCOUNT FOR EXTRA NOISE ADDED FROM SUBTRACTING FRAMES
        return dstd / math.sqrt(2), combinedMask

    def qc_ron(
        self,
        frameType=False,
        frameName=False,
        masterFrame=False,
        rawRon=False,
        masterRon=False,
    ):
        """*calculate the read-out-noise from bias/dark frames*

        **Key Arguments:**

        - ``frameType`` -- the type of the frame for reporting QC values. Default *False*
        - ``frameName`` -- the name of the frame in human readable words. Default *False*
        - ``masterFrame`` -- the master frame (only makes sense to measure RON on master bias). Default *False*
        - ``rawRon`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*
        - ``masterRon`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons
        - ``masterRon`` -- combined read-out-noise in mbias

        **Usage:**

        ```python
        rawRon, mbiasRon = self.qc_ron(
            frameType="MBIAS",
            frameName="master bias",
            masterFrame=masterFrame
        )
        ```
        """
        self.log.debug("starting the ``qc_bias_ron`` method")

        import numpy as np

        from soxspipe.commonutils import toolkit

        # ONE TIMESTAMP COVERS BOTH THE RAW AND MASTER RON ROWS, AS IT DID INLINE.
        utcnow = toolkit.utcnow_string()

        if not rawRon and len(self.inputFrames.files) > 1:
            rawRon, combinedMask = self._measure_raw_frame_ron()

        if rawRon:
            singleFrameType = frameType
            if frameType[0] == "M":
                singleFrameType = frameType[1:]

            self.add_qc(
                qcName="RAW RON",
                qcValue=rawRon,
                qcComment=f"[e-] RON in single {singleFrameType}",
                qcUnit="electrons",
                reductionDateUtc=utcnow,
                toHeader=True,
            )

        if masterFrame and not masterRon:
            # PREDICTED MASTER NOISE
            # predictedMasterRon = rawRon / math.sqrt(len(ccds))

            # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
            tmp = np.ma.array(masterFrame.data, mask=combinedMask)

            dmin, dmax, dmean, dstd = _image_stats(tmp)
            masterRon = float(dstd)

        elif masterRon:
            self.add_qc(
                qcName="MASTER RON",
                qcValue=float(masterRon),
                qcComment=f"[e-] Combined RON in {frameType}",
                qcUnit="electrons",
                reductionDateUtc=utcnow,
                toHeader=True,
            )
        else:
            masterRon = None

        self.log.debug("completed the ``qc_bias_ron`` method")
        return rawRon, masterRon

    def qc_median_flux_level(self, frame, frameType="MBIAS", frameName="master bias", medianFlux=False):
        """*calculate the median flux level in the frame, excluding masked pixels*

        **Key Arguments:**

        - ``frame`` -- the frame (CCDData object) to determine the median level.
        - ``frameType`` -- the type of the frame for reporting QC values Default "MBIAS"
        - ``frameName`` -- the name of the frame in human readable words. Default "master bias"
        - ``medianFlux`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*

        **Return:**

        - ``medianFlux`` -- median flux level in electrons

        **Usage:**

        ```python
        medianFlux = self.qc_median_flux_level(
            frame=myFrame,
            frameType="MBIAS",
            frameName="master bias")
        ```
        """
        self.log.debug("starting the ``qc_median_flux_level`` method")

        import numpy as np

        from soxspipe.commonutils import toolkit

        if not medianFlux:
            # DETERMINE MEDIAN BIAS LEVEL
            maskedDataArray = np.ma.array(frame.data, mask=frame.mask)
            medianFlux = np.ma.median(maskedDataArray)

        # THE VALUE IS UNREAD, BUT THE EXPTIME LOOKUP RAISES KeyError FOR A FRAME
        # WITHOUT THAT KEYWORD. DELETING THE ASSIGNMENT WOULD REMOVE THAT FAILURE,
        # WHICH IS PINNED BY test_base_recipe_lint_characterization.py.
        fluxRange = (np.nanpercentile(frame.data, 95) - np.nanpercentile(frame.data, 5)) / frame.header[  # noqa: F841
            self.kw("EXPTIME")
        ]

        utcnow = toolkit.utcnow_string()

        self.add_qc(
            qcName=f"{frameType} MEDIAN".upper(),
            qcValue=medianFlux,
            qcComment=f"[e-] Median flux level of {frameName}",
            qcUnit="electrons",
            reductionDateUtc=utcnow,
            toHeader=True,
        )

        self.log.debug("completed the ``qc_median_flux_level`` method")
        return medianFlux

    def subtract_mean_flux_level(self, rawFrame):
        """*iteratively median sigma-clip raw bias data frames before calculating and removing the mean bias level*

        **Key Arguments:**

        - ``rawFrame`` -- the raw bias frame

        **Return:**

        - `meanFluxLevel` -- the frame mean bias level
        - `fluxStd` -- the standard deviation of the flux distribution (RON)
        - `noiseFrame` -- the raw bias frame with mean bias level removed

        **Usage:**

        ```python
        meanFluxLevel, fluxStd, noiseFrame = self.subtract_mean_flux_level(rawFrame)
        ```
        """
        self.log.debug("starting the ``subtract_mean_flux_level`` method")

        import numpy as np
        from astropy.stats import sigma_clip

        from soxspipe.commonutils import toolkit

        # UNPACK SETTINGS
        clipping_lower_sigma = self.recipeSettings["frame-clipping-sigma"]
        clipping_iteration_count = self.recipeSettings["frame-clipping-iterations"]

        maskedFrame = sigma_clip(
            rawFrame.data.astype(np.float32),
            sigma=clipping_lower_sigma,
            maxiters=clipping_iteration_count,
            cenfunc="mean",
            stdfunc="std",
        )

        rawFrame.mask = maskedFrame.mask
        # DETERMINE MEDIAN BIAS LEVEL
        maskedDataArray = np.ma.array(maskedFrame.data, mask=maskedFrame.mask).astype(np.float32)

        meanFluxLevel = np.ma.mean(maskedDataArray)
        fluxStd = np.ma.std(maskedDataArray)
        rawFrame.data -= meanFluxLevel
        toolkit.frame_to_32(rawFrame)

        self.log.debug("completed the ``subtract_mean_flux_level`` method")
        return (meanFluxLevel, fluxStd, rawFrame)

    def _stamp_raw_frame_records(self, frame, tableData, rawFrames=False):
        """*record the raw frames used by this recipe in the product header*

        A raw frame is a row of ``tableData`` with no PRO TYPE value.

        **Key Arguments:**

        - ``frame`` -- the frame whose header is updated
        - ``tableData`` -- the recipe's raw-frame summary, as a pandas DataFrame
        - ``rawFrames`` -- limit the raw frames listed in the header to only these frames (list)

        **Usage:**

        ```python
        self._stamp_raw_frame_records(frame, tableData, rawFrames)
        ```
        """
        import math

        kw = self.kw

        iterator = 1
        for f, t, z in zip(
            tableData["filename"].values,
            tableData["tag"].values,
            tableData[kw("PRO_TYPE")].values,
        ):
            if rawFrames and f not in rawFrames:
                continue
            if isinstance(z, float) and math.isnan(z):
                valueLen = 80 - len(f"ESO PRO REC1 RAW{iterator} NAME" + "HIERARCH  = '")
                if len(f) > valueLen:
                    self.log.warning(f"The filename {f} has been trucated to {f[:valueLen]} in the FITS header")
                frame.header[f"ESO PRO REC1 RAW{iterator} NAME"] = f[:valueLen]
                frame.header[f"ESO PRO REC1 RAW{iterator} CATG"] = t
                iterator += 1
        return

    def _stamp_calibration_frame_records(self, frame, tableData):
        """*record the calibration frames used by this recipe in the product header*

        A calibration frame is a row of ``tableData`` that carries a PRO TYPE value.
        Each record carries the MD5 hash of the calibration file on disk.

        **Key Arguments:**

        - ``frame`` -- the frame whose header is updated
        - ``tableData`` -- the recipe's raw-frame summary, as a pandas DataFrame

        **Usage:**

        ```python
        self._stamp_calibration_frame_records(frame, tableData)
        ```
        """
        import math

        from astropy.utils.data import compute_hash

        kw = self.kw

        iterator = 1

        for f, c, z, p in zip(
            tableData["filename"].values,
            tableData[kw("PRO_CATG")].values,
            tableData[kw("PRO_TYPE")].values,
            tableData["file"].values,
        ):
            if not isinstance(z, float) or not math.isnan(z):
                valueLen = 80 - len(f"ESO PRO REC1 CAL{iterator} NAME" + "HIERARCH  = '")
                # if len(f) > valueLen:
                #     self.log.warning(f"The filename {f} has been trucated to {f[:valueLen]} in the FITS header")
                frame.header[f"ESO PRO REC1 CAL{iterator} NAME"] = f[:valueLen]
                frame.header[f"ESO PRO REC1 CAL{iterator} CATG"] = c
                frame.header[f"ESO PRO REC1 CAL{iterator} DATAMD5"] = compute_hash(p)
                iterator += 1
        return

    def _stamp_recipe_parameter_records(self, frame):
        """*record the recipe's own settings in the product header*

        Arm-specific settings and the QC ranges are skipped; a nested settings
        block is flattened into one record per inner key.

        **Key Arguments:**

        - ``frame`` -- the frame whose header is updated

        **Usage:**

        ```python
        self._stamp_recipe_parameter_records(frame)
        ```
        """
        iterator = 1
        recipeSettings = self.get_recipe_settings()
        for k, v in recipeSettings.items():
            if k.lower() in ["uvb", "vis", "nir", "qc-acceptable-ranges"]:
                continue

            if not isinstance(v, dict):
                if isinstance(v, list):
                    v = ", ".join(map(str, v)).strip()

                frame.header[f"ESO PRO REC1 PARAM{iterator} NAME"] = k[:40]
                frame.header[f"ESO PRO REC1 PARAM{iterator} VALUE"] = v
                iterator += 1
            else:
                for k2, v2 in v.items():
                    frame.header[f"ESO PRO REC1 PARAM{iterator} NAME"] = k2[:40]
                    frame.header[f"ESO PRO REC1 PARAM{iterator} VALUE"] = v2
                    iterator += 1
        return

    def update_fits_keywords(self, frame, rawFrames=False):
        """*update fits keywords to comply with ESO Phase 3 standards*

        **Key Arguments:**

        - ``frame`` -- the frame to update
        - ``rawFrames`` -- limit the raw frames to be listed in the fits header to only these frames (list)

        **Usage:**

        ```python
        usage code
        ```

        :::{todo}
            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        :::
        """
        self.log.debug("starting the ``update_fits_keywords`` method")

        import soxspipe.__version__ as version

        arm = self.arm
        kw = self.kw
        # UNUSED, AND ONE OF THE MODULE'S ELEVEN `F841` FINDINGS. DELETING IT IS DY-88'S.
        dp = self.detectorParams  # noqa: F841
        imageType = self.imageType
        if "FLAT" in imageType:
            imageType = "FLAT"

        frame.header[kw("SEQ_ARM").upper()] = arm
        frame.header[kw("PRO_TYPE").upper()] = "REDUCED"

        # PROD CATG
        if imageType in ["BIAS", "DARK", "FLAT"]:
            frame.header[kw("PRO_CATG")] = f"MASTER_{imageType}_{arm}".replace("QLAMP", "LAMP").replace("DLAMP", "LAMP")
            frame.header[kw("PRO_TECH")] = "IMAGE"

        # SPEC FORMAT TO PANDAS DATAFRAME
        tableData = self.rawFrames.to_pandas()

        tableData["filename"] = tableData["filename"].str.replace("_pre", "")
        tableData["tag"] = tableData["TYPE"] + "_" + tableData["ARM"]

        self._stamp_raw_frame_records(frame, tableData, rawFrames)
        self._stamp_calibration_frame_records(frame, tableData)
        self._stamp_recipe_parameter_records(frame)

        # SOXSPIPE VERSION
        frame.header["ESO PRO REC1 PIPE ID"] = f"soxspipe/v{version}"

        # RECIPE
        if self.recipeName:
            frame.header["ESO PRO REC1 ID"] = self.recipeName

        # from tabulate import tabulate
        # print(tabulate(tableData, headers='keys', tablefmt='github'))

        self.log.debug("completed the ``update_fits_keywords`` method")
        return

    def get_recipe_settings(self):
        """*get the recipe and arm specific settings*

        **Return:**

        - ``recipeSettings`` -- the recipe specific settings

        **Usage:**

        ```python
        usage code
        ```
        """
        self.log.debug("starting the ``get_recipe_settings`` method")

        recipeSettings = False
        if self.recipeName:
            recipeSettings = self.settings[self.recipeName]
        if recipeSettings and self.arm and self.arm.lower() in recipeSettings:
            for k, v in recipeSettings[self.arm.lower()].items():
                recipeSettings[k] = v

        self.log.debug("completed the ``get_recipe_settings`` method")
        return recipeSettings

    def _dataframe_to_sqlite(self, dataframe, table_name, replace=False):
        """*write a dataframe to a database table, retrying the insert up to seven times*

        **Key Arguments:**

        - ``dataframe`` -- the dataframe containing the rows to insert
        - ``table_name`` -- the name of the database table to insert into
        - ``replace`` -- if True, delete the table's existing rows first; otherwise append. Default: False

        **Raises:**

        - Exception if the insert fails after seven attempts.
        - ``UnsafeSqlIdentifierError`` if ``table_name`` fails the safe-identifier grammar.
        """
        import time

        # A TABLE NAME CANNOT BE A BOUND PARAMETER, SO IT IS VALIDATED AGAINST
        # THE SAFE-IDENTIFIER GRAMMAR ONCE HERE, BEFORE EITHER SQL SITE BELOW
        # INTERPOLATES IT.
        table_name = validate_sql_identifier(table_name, "table name")

        if replace:
            c = self.conn.cursor()
            sqlQuery = f"delete from {table_name};"  # noqa: S608
            try:
                c.execute(sqlQuery)
            except sqlite3.OperationalError as e:
                self.log.debug(f"_dataframe_to_sqlite: `c.execute(sqlQuery)` failed, continuing: {e}")
            c.close()

        keepTrying = 0
        keepTryingMax = 7
        while keepTrying < keepTryingMax:
            try:
                dataframe.replace(["--"], None).to_sql(
                    table_name,
                    con=self.conn,
                    index=False,
                    if_exists="append",
                    method="multi",
                )
                keepTrying = keepTryingMax
            except Exception as e:
                if keepTrying > keepTryingMax - 1:
                    raise Exception(e)
                time.sleep(1)
                keepTrying += 1

    def add_qc(
        self,
        qcName,
        qcValue,
        qcComment,
        reductionDateUtc,
        qcUnit=OMITTED,
        toHeader=OMITTED,
        qcOrder=OMITTED,
        recipeName=None,
        obsDateUtc=None,
    ):
        """*append a QC row to `self.qc`, saving recipe call sites from repeating `self.qc`,
        `self.recipeName` and `self.dateObs`*

        This is a thin delegator around
        `soxspipe.commonutils.toolkit.append_qc`. It builds no row and makes
        no decision about one: the only thing it does beyond forwarding is
        fall back to `self.recipeName` and `self.dateObs` when ``recipeName``
        and ``obsDateUtc`` are not supplied. Both can be overridden, because
        some call sites hardcode a different recipe name, or need to record
        a different observation date. The optional
        ``qcUnit``/``toHeader``/``qcOrder`` keywords are only forwarded when
        the caller passes them, so an unpassed optional column stays absent
        from the appended row exactly as `append_qc` would leave it.

        **Key Arguments:**

        - ``qcName`` -- the QC metric name
        - ``qcValue`` -- the QC metric value
        - ``qcComment`` -- the QC metric comment
        - ``reductionDateUtc`` -- the reduction date (UTC) to record against the QC row. Shared
          across several rows by the caller, never generated here
        - ``qcUnit`` -- the QC metric unit. Omit to leave the column absent for this row
        - ``toHeader`` -- whether the QC metric should be written to the FITS header. Omit to
          leave the column absent for this row
        - ``qcOrder`` -- the echelle order the QC metric applies to. Omit to leave the column absent for this row
        - ``recipeName`` -- override for the recipe name recorded against the QC row. Defaults to `self.recipeName`
        - ``obsDateUtc`` -- override for the observation date recorded against the QC row. Defaults
          to `self.dateObs`

        **Usage:**

        ```python
        self.add_qc(
            qcName="RON",
            qcValue=1.2,
            qcComment="[e-] RON in single BIAS",
            reductionDateUtc=utcnow,
            qcUnit="electron",
            toHeader=True,
        )
        ```
        """
        # DELIBERATELY A FUNCTION-LOCAL IMPORT, MATCHING EVERY OTHER `toolkit`
        # IMPORT IN THIS FILE. ONLY `OMITTED` IS IMPORTED AT MODULE LEVEL, AND
        # ONLY BECAUSE A DEFAULT ARGUMENT VALUE IS EVALUATED AT MODULE LOAD AND
        # SO CANNOT COME FROM A LOCAL IMPORT. THE ASYMMETRY IS INTENTIONAL: DO
        # NOT HOIST THESE WITHOUT RE-TESTING THE COLD-IMPORT PATHS, SINCE THIS
        # PACKAGE'S `__init__` IMPORT ORDER IS LOAD-BEARING.
        from soxspipe.commonutils.toolkit import append_qc

        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName if recipeName is None else recipeName,
            qcName=qcName,
            qcValue=qcValue,
            qcComment=qcComment,
            obsDateUtc=self.dateObs if obsDateUtc is None else obsDateUtc,
            reductionDateUtc=reductionDateUtc,
            qcUnit=qcUnit,
            toHeader=toHeader,
            qcOrder=qcOrder,
        )
        return

    def add_product(
        self,
        productLabel,
        fileName,
        filePath,
        productDesc,
        reductionDateUtc,
        fileType=OMITTED,
        label=OMITTED,
        recipeName=None,
        obsDateUtc=None,
    ):
        """*append a product row to `self.products`, saving recipe call sites from repeating
        `self.products`, `self.recipeName` and `self.dateObs`*

        This is a thin delegator around
        `soxspipe.commonutils.toolkit.append_product`. It builds no row and
        makes no decision about one: the only thing it does beyond
        forwarding is fall back to `self.recipeName` and `self.dateObs` when
        ``recipeName`` and ``obsDateUtc`` are not supplied. Both can
        be overridden: several product rows hardcode a literal recipe name
        regardless of the running recipe, and `base_recipe` itself rewrites
        `self.recipeName` to a `-std` variant for standard-star input. The
        optional ``fileType``/``label`` keywords are only forwarded when the
        caller passes them, so an unpassed optional column stays absent
        from the appended row exactly as `append_product` would leave it.

        **Key Arguments:**

        - ``productLabel`` -- the product label
        - ``fileName`` -- the product file name
        - ``filePath`` -- the product file path
        - ``productDesc`` -- the product description
        - ``reductionDateUtc`` -- the reduction date (UTC) to record against the product row. Shared
          across several rows by the caller, never generated here
        - ``fileType`` -- the product file type. Omit to leave the column absent for this row
        - ``label`` -- the product label category (e.g. `PROD`). Omit to leave the column absent for this row
        - ``recipeName`` -- override for the recipe name recorded against the product row. Defaults to `self.recipeName`
        - ``obsDateUtc`` -- override for the observation date recorded against the product row. Defaults
          to `self.dateObs`

        **Usage:**

        ```python
        self.add_product(
            productLabel="MBIAS",
            fileName=filename,
            filePath=productPath,
            productDesc=f"{self.arm} Master bias frame",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )
        ```
        """
        # DELIBERATELY A FUNCTION-LOCAL IMPORT, MATCHING EVERY OTHER `toolkit`
        # IMPORT IN THIS FILE. ONLY `OMITTED` IS IMPORTED AT MODULE LEVEL, AND
        # ONLY BECAUSE A DEFAULT ARGUMENT VALUE IS EVALUATED AT MODULE LOAD AND
        # SO CANNOT COME FROM A LOCAL IMPORT. THE ASYMMETRY IS INTENTIONAL: DO
        # NOT HOIST THESE WITHOUT RE-TESTING THE COLD-IMPORT PATHS, SINCE THIS
        # PACKAGE'S `__init__` IMPORT ORDER IS LOAD-BEARING.
        from soxspipe.commonutils.toolkit import append_product

        self.products = append_product(
            self.products,
            recipeName=self.recipeName if recipeName is None else recipeName,
            productLabel=productLabel,
            fileName=fileName,
            filePath=filePath,
            productDesc=productDesc,
            obsDateUtc=self.dateObs if obsDateUtc is None else obsDateUtc,
            reductionDateUtc=reductionDateUtc,
            fileType=fileType,
            label=label,
        )
        return

    # use the tab-trigger below for new method
    # xt-class-method
