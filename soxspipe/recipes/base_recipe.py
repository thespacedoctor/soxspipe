#!/usr/bin/env python
"""
*The base recipe class which all other recipes inherit*

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""

from __future__ import annotations

################# GLOBAL IMPORTS ####################
import os
import sqlite3
import sys
from typing import Any

from soxspipe.commonutils import detector_lookup, filenamer, keyword_lookup, subtract_background
from soxspipe.commonutils.scratch import (
    absolute_path,
    register_scratch_directory,
    remove_scratch_tree,
    scratch_root,
    unregister_scratch_directory,
)

# VALIDATES A COLUMN OR TABLE NAME AT THE SQL TRUST BOUNDARY BEFORE IT IS
# INTERPOLATED. A LEAF MODULE WITH NO IMPORTS OF ITS OWN, SO IT CANNOT
# DISTURB THE PACKAGE'S LOAD-BEARING IMPORT ORDER.
from soxspipe.commonutils.sql_identifiers import validate_sql_identifier
from soxspipe.commonutils.sqlite_dataframe_writer import write_dataframe_to_sqlite

# THE `_STD_` RULE THAT NAMES A STANDARD-STAR SOF. ANOTHER LEAF MODULE, WITH NO IMPORTS AT MODULE LEVEL
from soxspipe.commonutils.stare_mode import is_standard_input

# THE "CALLER DID NOT PASS THIS KEYWORD" SENTINEL USED BY THE add_qc AND
# add_product DELEGATORS BELOW. IMPORTED FROM `toolkit` RATHER THAN REDEFINED
# HERE SO THE DELEGATORS FORWARD THE SAME OBJECT `append_qc` AND
# `append_product` TEST FOR, AND SO A CALLER MAY PASS `OMITTED` EXPLICITLY.
from soxspipe.commonutils.toolkit import OMITTED, PREPARED_FRAME_KEYWORD
from soxspipe.commonutils.unusable_flat_pixels import (
    flag_unusable_pixels,
    neutralise_unusable_flat,
    unusable_flat_pixel_mask,
)
from soxspipe.commonutils.workspace_paths import expand_workspace_root

os.environ["TERM"] = "vt100"

# DARK AND FRAME EXPOSURE TIMES CLOSER THAN THIS (IN SECONDS) COUNT AS EQUAL. IT ONLY
# ABSORBS FLOAT NOISE IN THE HEADER VALUES; ANY REAL DIFFERENCE IS A MISMATCH.
EXPTIME_MATCH_TOLERANCE_SECONDS = 1e-3


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


class RecipeProductExistsError(FileExistsError):
    """A recipe refused to run because its product, or the error log of an earlier failed run, is already in place."""


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
    - ``productNameMarker`` -- a marker the product file names carry after the SOF name. Default "" (no marker)


    **Usage**

    To use this base recipe to create a new `soxspipe` recipe, have a look at the code for one of the simpler recipes (e.g. `soxs_mbias`) - copy and modify the code.
    """

    # THE READ NOISE RESOLVED FROM THE MASTER BIAS IN THE SOF (A ``ReadNoise``); NONE WITHOUT A MASTER BIAS (DY-1458)
    readNoise = None
    # THE NOMINAL READ NOISE (E-) IN THE RAW-FRAME HEADERS (NONE IF THEY CARRY NONE) AND THE DETECTOR DEFAULT (E-)
    rawHeaderRon = None
    defaultRon = None
    # THE COSMIC-RAY MASK PRODUCT WRITTEN BY FRAME PREPARATION, WAITING TO BE REGISTERED BY THE FIRST `report_output`
    cosmicRayProduct = None
    # A MARKER THE PRODUCT FILE NAMES CARRY AFTER THE SOF NAME, SO TWO REDUCTIONS OF ONE SOF NEVER SHARE A FILE NAME
    # (DY-2194). DATABASE AND QC KEYS KEEP THE BARE SOF NAME. A CLASS ATTRIBUTE, SO A RECIPE OBJECT BUILT WITHOUT
    # `__init__` HAS ONE TOO
    productNameMarker = ""
    # THE ARM OF THE INPUT FRAMES, SET WHEN THEY ARE VERIFIED. A CLASS ATTRIBUTE, SO A RECIPE OBJECT BUILT WITHOUT
    # `__init__` READS AS NO ARM (DY-2322)
    arm = None
    # WHETHER A NIR RECIPE SUBTRACTS THE MASTER DARK IN THE SOF. A SOXS RECIPE THEN ADDS THE DARK'S NOISE TO ITS SKY
    # FIT (THE READ NOISE ITSELF COMES FROM THE FRAME'S NDSAMPLES, DY-2341); AN X-SHOOTER RECIPE TAKES ITS READ NOISE
    # FROM THE DARK. A NODDING RECIPE SETS IT FALSE, AS A-B REMOVES THE DARK (DY-2335)
    usesMasterDarkReadNoise: bool = True

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
        productNameMarker="",
    ):
        import tempfile

        import matplotlib
        import pandas as pd

        from soxspipe.commonutils.toolkit import ignore_expected_library_warnings

        log.debug("instantiating a new '__init__' object")
        self.recipeName = recipeName
        self.productNameMarker = productNameMarker
        # Merge before any constructor work reads settings.  Keep the merge
        # non-mutating so callers can safely reuse their input dictionary.
        settings = {**self._advanced_settings(), **settings}
        self.settings = settings
        self.debug = debug
        self.workspaceRootPath = self._absolute_path(self.settings["workspace-root-dir"])
        self.turnOffMP = turnOffMP

        if self.debug:
            matplotlib.use("TkAgg")
        else:
            matplotlib.use("Agg")

        self.darkDetrendWarningIssued2 = False

        # EVERY REDUCTION PATH BUILDS A RECIPE, INCLUDING SPAWNED WORKERS, WHICH DO NOT INHERIT THE PARENT'S FILTERS
        ignore_expected_library_warnings()

        if is_standard_input(inputFrames):
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

        # `tmp` HERE IS A SUBDIRECTORY OF THE USER'S OWN WORKSPACE, NOT THE
        # SYSTEM TEMPORARY DIRECTORY. `mkdtemp` CREATES A FRESH DIRECTORY, SO
        # CONCURRENT RECIPES NEVER SHARE ONE AND ONE RECIPE'S `clean_up` CANNOT
        # DELETE ANOTHER'S FRAMES (DY-49).
        self.scratchRoot = scratch_root(self.settings)
        os.makedirs(self.scratchRoot, exist_ok=True)
        self.outDir = tempfile.mkdtemp(dir=self.scratchRoot)
        register_scratch_directory(self.outDir)

        # FIND THE CURRENT SESSION
        from soxspipe.commonutils import data_organiser

        do = data_organiser(
            log=self.log,
            rootDir=expand_workspace_root(self.settings["workspace-root-dir"]),
            dbConnect=False,
        )
        self.currentSession, allSessions = do.session_list(silent=True)
        do.close()

        # INITIATE A DB CONNECTION
        self.conn = None
        self.status = None
        if not self.turnOffMP:
            self._open_session_database()

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        self.qc, self.products = self._empty_qc_and_product_tables(pd)

        # THE COSMIC-RAY PIXELS L.A.COSMIC FLAGGED IN EACH RAW FRAME, WRITTEN OUT ONCE FRAME PREPARATION ENDS
        self.cosmicRayMasks = []
        self.cosmicRayProduct = None

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=self.settings,
            recipeName=self.recipeName,
            startNightDate=self.startNightDate,
        )

        self.generateReponseCurve = False

        return

    @property
    def productStem(self) -> str:  # noqa: N802
        """*the stem of every product file name: the SOF name followed by the product-name marker*

        **Return:**

        - ``productStem`` -- the SOF name and the marker
        """
        return f"{self.sofName}{self.productNameMarker}"

    @property
    def runSerially(self) -> bool:  # noqa: N802
        """*True when work must run in this process: multiprocessing is turned off or debug mode is on*"""
        return bool(self.turnOffMP or self.debug)

    def _complete_product_exists(self, log):
        """*test whether a complete product is already in place at* ``self.productPath``

        **Key Arguments:**

        - ``log`` -- the logger to warn on

        **Return:**

        - ``exists`` -- *True* only when the product is a whole FITS file. A file that is cut short by an
          interrupted write is warned about by name and counted as absent, so the rerun rebuilds it.
        """
        from soxspipe.commonutils.atomic_fits import is_complete_fits

        if not os.path.exists(self.productPath):
            return False
        if is_complete_fits(self.productPath):
            return True
        log.warning(
            f"The existing product `{os.path.basename(self.productPath)}` is not a complete FITS file "
            "(an earlier run was probably interrupted). It is being replaced."
        )
        return False

    def _resolve_product_path(self, log, inputFrames, verbose, overwrite):
        """*predict where this recipe's product will be written, and refuse to run over one that is already there*

        **Key Arguments:**

        - ``log`` -- the logger the recipe was handed
        - ``inputFrames`` -- the input frames the recipe was handed
        - ``verbose`` -- also print the refusal to the terminal
        - ``overwrite`` -- run again over an existing product or error log

        Sets ``self.sofName``, ``self.productPath``, ``self.startNightDate`` and
        ``self.log``. A product or an error log already in place raises
        ``RecipeProductExistsError`` (a ``FileExistsError``) unless ``overwrite`` is set. That refusal applies to
        a ``.sof`` input only.

        A list of frame paths or a directory has no ``.sof`` name to read, so
        ``startNightDate`` is the start-of-night date of the earliest
        ``MJD-OBS`` among the frames (``toolkit.start_night_date``) and
        ``sofName`` is a synthetic name, ``<recipeName>_<ARM>_<startNightDate>``
        (for example ``soxs-mbias_VIS_2024-01-01``). The arm is that of the
        earliest frame, so two arms reduced on one night do not delete each
        other's QC rows. The QC table and the session database key on that name
        with a ``.sof`` suffix, exactly as for a real set-of-files.
        ``productPath`` is *False*: the product path is not predicted. Input
        with no frames, or with no readable ``MJD-OBS``, raises ``ValueError``.
        """
        from soxspipe.commonutils import toolkit

        # CHECK IF PRODUCT ALREADY EXISTS
        if inputFrames and not isinstance(inputFrames, list) and inputFrames.split(".")[-1].lower() == "sof":
            self.sofName = os.path.basename(inputFrames).replace(".sof", "")
            self.productPath, self.startNightDate = toolkit.predict_product_path(
                inputFrames,
                self.recipeName,
                self.productNameMarker,
                outputDir=self.settings["workspace-root-dir"],
            )
            errorLog = os.path.splitext(self.productPath)[0] + "_ERROR.log"
            basename = os.path.basename(self.productPath)

            if os.path.exists(errorLog) and not overwrite:
                if verbose:
                    print(
                        f"This recipe previously failed (see `{basename}`). To rerun the recipe, run the recipe command with the overwrite flag (-x)."
                    )
                raise RecipeProductExistsError(
                    f"This recipe previously failed (see `{basename}`). To rerun the recipe, run the recipe command with the overwrite flag (-x)."
                )
            if not overwrite and self._complete_product_exists(log):
                basename = os.path.basename(self.productPath)
                if verbose:
                    print(
                        f"The product of this recipe already exists: `{basename}`. To overwrite this product, rerun the pipeline command with the overwrite flag (-x)."
                    )
                raise RecipeProductExistsError(
                    f"The product of this recipe already exists: `{basename}`. To overwrite this product, rerun the pipeline command with the overwrite flag (-x)."
                )

            self.log = toolkit.add_recipe_logger(log, self.productPath)
        else:
            self.startNightDate, arm = toolkit.start_night_date_from_frames(
                log=log, settings=self.settings, inputFrames=inputFrames
            )
            self.sofName = f"{self.recipeName}_{arm}_{self.startNightDate}"
            self.productPath = False
            self.log = log

        return

    def _open_session_database(self):
        """*open the workspace database and mark this recipe as failed until it completes*

        Sets ``self.conn`` and ``self.status``. The stored status is set to
        'fail' up front, so a recipe that crashes leaves a failure behind
        rather than its previous result. If that write errors, the connection
        is closed and a ``RuntimeError`` stops the recipe (DY-2263).
        """
        from soxspipe.commonutils.workspace_database import connect_workspace_database

        if self.currentSession and self.sofName:
            self.sessionDb = expand_workspace_root(self.settings["workspace-root-dir"]) + "/soxspipe.db"

            def dict_factory(cursor, row):
                d = {}
                for idx, col in enumerate(cursor.description):
                    d[col[0]] = row[idx]
                return d

            self.conn = connect_workspace_database(self.sessionDb, rowFactory=dict_factory)

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
            row = c.fetchone()
            if row is None:
                # NO ROW FOR THIS SOF, SO NO EARLIER STATUS CAN GO STALE
                self.log.warning(
                    f"__init__: `self.status = c.fetchone()['status']` failed, continuing: "
                    f"no product_frames row for `{self.sofName}.sof`"
                )
                self.status = None
            else:
                self.status = row["status"]
                self._write_start_of_run_fail_status(c, statusColumn)

            c.close()

        return

    def _write_start_of_run_fail_status(self, cursor: sqlite3.Cursor, statusColumn: str) -> None:
        """*write* ``fail`` *to this recipe's product row, and stop the recipe if the write errors*

        **Key Arguments:**

        - ``cursor`` -- a cursor on the open workspace database
        - ``statusColumn`` -- the validated name of the session status column

        A recipe that ran on without the write could crash before ``clean_up`` and leave the ``pass`` of an
        earlier run in place, so downstream recipes would trust a product that the last run did not finish. The
        connection is closed before the error is raised, and the write is not retried: the connection
        already waits for the busy timeout (DY-2263).
        """
        sofFile = f"{self.sofName}.sof"
        sqlQuery = f"update product_frames set {statusColumn} = 'fail' where sof = ?"  # noqa: S608
        try:
            cursor.execute(sqlQuery, (sofFile,))
        except sqlite3.Error as e:
            from contextlib import suppress

            # BEST-EFFORT CLEANUP, SO A FAILING CLOSE CANNOT HIDE THE WRITE ERROR
            with suppress(sqlite3.Error):
                cursor.close()
            with suppress(sqlite3.Error):
                self.conn.close()
            raise RuntimeError(
                f"The recipe could not write the `fail` status for `{sofFile}` to the workspace database "
                f"and was stopped, because a crash would leave the status of an earlier run in place: {e}"
            ) from e

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
                "soxspipe_recipe": pd.Series([], dtype="object"),
                "qc_name": pd.Series([], dtype="object"),
                "qc_value": pd.Series([], dtype="object"),
                "qc_unit": pd.Series([], dtype="object"),
                "qc_order": pd.Series([], dtype="object"),
                "qc_comment": pd.Series([], dtype="object"),
                "obs_date_utc": pd.Series([], dtype="object"),
                "reduction_date_utc": pd.Series([], dtype="object"),
                "to_header": pd.Series([], dtype="object"),
            }
        )
        products = pd.DataFrame(
            {
                "soxspipe_recipe": pd.Series([], dtype="object"),
                "product_label": pd.Series([], dtype="object"),
                "file_name": pd.Series([], dtype="object"),
                "file_type": pd.Series([], dtype="object"),
                "obs_date_utc": pd.Series([], dtype="object"),
                "reduction_date_utc": pd.Series([], dtype="object"),
                "file_path": pd.Series([], dtype="object"),
                "label": pd.Series([], dtype="object"),
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

        import ccdproc
        from astropy import units as u

        from soxspipe.commonutils import toolkit
        from soxspipe.commonutils.fits_headers import read_ccd, set_header_card
        from soxspipe.commonutils.read_noise import SCIENCE_RON_CARD

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
                frame = read_ccd(
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

        # CHECK THE NUMBER OF EXTENSIONS IS ONLY 1 AND THE PREPARED-FRAME MARKER DOES NOT
        # EXIST. i.e. THIS IS A RAW UNTOUCHED FRAME
        if len(frame.to_hdu()) > 1 or PREPARED_FRAME_KEYWORD in frame.header:
            return filepath

        # MANIPULATE XSH DATA
        frame = self.xsh2soxs(frame)
        # THE PRESCAN IS TRIMMED AWAY, SO ITS NOISE IS MEASURED FIRST
        prescanRon = self._measure_prescan_read_noise(frame, filepath)
        frame = self._trim_frame(frame)

        # CORRECT FOR GAIN - CONVERT DATA FROM ADU TO ELECTRONS
        frame = ccdproc.gain_correct(frame, dp["gain"], add_keyword=None)
        toolkit.frame_to_32(frame)

        frame = self._add_uncertainty_map(frame, readNoise=prescanRon)
        if prescanRon is not None:
            set_header_card(
                frame.header, SCIENCE_RON_CARD, prescanRon, "[e-] RON of the science frame, from its prescan"
            )

        frame.mask = self._bad_pixel_mask(frame)

        # THE NESTED `if` IS KEPT SO THE TWO CONDITIONS STAY SEPARATELY
        # COMMENTED, AND BECAUSE COLLAPSING IT IS COSMETIC WORK BELONGING TO
        # DY-88.
        if self.recipeName in ["soxs-nod-std", "soxs-stare-std", "soxs-offset-std"] and self.recipeSettings["use_flat"]:  # noqa: SIM102
            # OBJECT/STANDARD FRAMES
            if frame.meta[kw("DPR_TYPE")] == "STD,FLUX" or "STD_stare" in frame.meta[kw("OBS_NAME")]:
                # ASSUMING WE HAVE ONLY STANDARD A-B CYCLES AND NOT JITTER.
                self.generateReponseCurve = True

        frame = self._clean_cosmic_rays(frame, filepath)

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

    def _add_uncertainty_map(self, frame, readNoise=None):
        """*add the uncertainty extension to a gain-corrected frame*

        **Key Arguments:**

        - ``frame`` -- the gain-corrected frame (CCDData object)
        - ``readNoise`` -- the read noise (e-) of this frame, measured from its prescan. Default *None*, for the
          read noise of the reduction (``self.detectorParams["ron"]``)

        **Return:**

        - ``frame`` -- the frame with its uncertainty extension populated

        A bias frame's uncertainty is the readnoise alone. Every other frame
        also carries the photon noise of its own counts. A frame that is the
        mean of NDIT integrations has its variance divided by NDIT (DY-2317).
        """
        import ccdproc
        import numpy as np
        from astropy import units as u

        from soxspipe.commonutils import toolkit
        from soxspipe.commonutils.ndit import read_ndit

        kw = self.kw
        dp = self.detectorParams
        ron = dp["ron"] if readNoise is None else readNoise * u.electron

        # GENERATE UNCERTAINTY MAP AS EXTENSION
        if frame.header[kw("DPR_TYPE")] == "BIAS":
            # ERROR IS ONLY FROM READNOISE FOR BIAS FRAMES
            errorMap = np.ones_like(frame.data).astype(np.float32) * ron.astype(np.float32)
            # errorMap = StdDevUncertainty(errorMap)
            frame.uncertainty = errorMap.astype(np.float32)
        else:
            # GENERATE UNCERTAINTY MAP AS EXTENSION
            frame = ccdproc.create_deviation(frame, readnoise=ron, disregard_nan=True, add_keyword=None)

        # A RAW NIR FILE IS THE MEAN OF ITS NDIT INTEGRATIONS, SO BOTH THE POISSON AND THE READ-NOISE TERMS SCALE
        ndit = read_ndit(frame.header, kw, self.log)
        if ndit > 1:
            frame.uncertainty.array = frame.uncertainty.array / np.float32(np.sqrt(ndit))
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

    def _clean_cosmic_rays(self, frame, filepath):
        """*flag cosmic rays in the frame, if the recipe settings ask for it, and keep the pixels newly flagged*

        **Key Arguments:**

        - ``frame`` -- the frame to clean (CCDData object)
        - ``filepath`` -- the path of the raw frame it was read from

        **Return:**

        - ``frame`` -- the cleaned frame, or the frame unchanged when cleaning is switched off

        The pixels L.A.Cosmic newly flagged, without the pixels already in the bad-pixel mask, are appended to
        ``self.cosmicRayMasks`` for `_write_cosmic_ray_mask_product`.
        """
        from soxspipe.commonutils import toolkit

        if (
            "use_lacosmic" in self.recipeSettings
            and self.recipeSettings["use_lacosmic"]
            and "lacosmic-sigma" in self.recipeSettings[self.arm.lower()]
            and self.recipeSettings[self.arm.lower()]["lacosmic-sigma"]
        ):
            maskBefore = frame.mask.copy()

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
            self.cosmicRayMasks.append(self._cosmic_ray_entry(frame, maskBefore, filepath))

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

    def _cosmic_ray_entry(self, frame, maskBefore, filepath):
        """*describe the cosmic rays L.A.Cosmic flagged in one frame*

        **Key Arguments:**

        - ``frame`` -- the cleaned frame (CCDData object)
        - ``maskBefore`` -- the mask of the frame before cleaning
        - ``filepath`` -- the path of the raw frame it was read from

        **Return:**

        - ``entry`` -- dictionary with the keys ``fileName``, ``mjdObs`` (None when the header has no finite
          number), ``dateObs`` and ``mask`` (uint8, 1 where cleaning newly flagged a pixel)
        """
        import math

        from soxspipe.commonutils.cosmic_ray_mask import new_cosmic_ray_mask

        kw = self.kw

        try:
            mjdObs = float(frame.header.get(kw("MJDOBS")))
        except (TypeError, ValueError):
            mjdObs = None
        if mjdObs is not None and not math.isfinite(mjdObs):
            mjdObs = None

        isPath = isinstance(filepath, (str, os.PathLike))
        return {
            "fileName": os.path.basename(os.fspath(filepath)) if isPath else "unknown",
            "mjdObs": mjdObs,
            "dateObs": frame.header.get(kw("DATE_OBS")),
            "mask": new_cosmic_ray_mask(maskBefore, frame.mask),
        }

    def _write_cosmic_ray_mask_product(self):
        """*write the cosmic-ray masks of the cleaned frames as one QC product*

        Writes nothing when no frame was cleaned. Otherwise writes ``<sof name>_LACOSMIC_MASK.fits`` to the QC
        directory, one image extension per cleaned frame in observation order. The product row and the QC row
        ``COSMIC RAY PIX NUM`` (the total number of pixels flagged, over all frames) are added by the first
        `report_output`, so they carry the name of the SOF whatever pass of the recipe reports first. The number
        in each frame is in the header of its extension.
        """
        self.cosmicRayProduct = None
        if not self.cosmicRayMasks:
            return

        from soxspipe.commonutils.cosmic_ray_mask import cosmic_ray_entry_order, cosmic_ray_mask_hdulist

        # OBSERVATION ORDER, NOT THE ORDER THE FRAMES WERE LISTED IN, SO THE EXTENSIONS ARE ALWAYS IN THE SAME ORDER
        entries = sorted(self.cosmicRayMasks, key=cosmic_ray_entry_order)
        filename = f"{self.productStem}_LACOSMIC_MASK.fits"
        filePath = os.path.join(self.qcDir, filename)
        cosmic_ray_mask_hdulist(entries).writeto(filePath, overwrite=True)

        self.cosmicRayProduct = {
            "sofName": self.sofName,
            "fileName": filename,
            "filePath": filePath,
            "dateObs": entries[0]["dateObs"],
            "count": sum(int(entry["mask"].sum()) for entry in entries),
        }

    def _register_cosmic_ray_product(self):
        """*add the cosmic-ray mask product and its QC to the tables, once*

        **Return:**

        - ``sofName`` -- the name of the SOF the mask was written for, or None when there is nothing to register
        """
        product = self.cosmicRayProduct
        if product is None:
            return None
        self.cosmicRayProduct = None

        from soxspipe.commonutils.toolkit import utcnow_string

        utcnow = utcnow_string()
        self.add_product(
            productLabel="LACOSMIC_MASK",
            fileName=product["fileName"],
            filePath=product["filePath"],
            productDesc="Pixels L.A.Cosmic flagged as cosmic rays in each raw frame, one image extension per frame",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="QC",
            obsDateUtc=product["dateObs"],
        )
        self.add_qc(
            qcName="COSMIC RAY PIX NUM",
            qcValue=product["count"],
            qcComment="Number of pixels L.A.Cosmic flagged as cosmic rays, over all raw frames",
            reductionDateUtc=utcnow,
            qcUnit="pixel",
            toHeader=False,
        )
        return product["sofName"]

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
        from soxspipe.commonutils.fits_frame_names import (
            FITS_SUFFIX,
            fits_frame_stem,
            is_fits_frame,
        )
        from soxspipe.commonutils.fits_headers import set_header_card

        if save:
            outDir = self.workspaceRootPath
        else:
            outDir = self.outDir

        # INJECT THE PRE KEYWORD
        set_header_card(frame.header, PREPARED_FRAME_KEYWORD, toolkit.utcnow_string(microseconds=True), "UTC timestamp")

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(outDir):
            try:
                os.makedirs(outDir)
            except OSError as e:
                self.log.debug(f"_prepare_single_frame: `os.makedirs(outDir)` failed, continuing: {e}")
        # CONVERT CCDData TO FITS HDU (INCLUDING HEADER) AND SAVE WITH PRE TAG
        # PREPENDED TO FILENAME
        # A .fits.Z INPUT IS WRITTEN UNCOMPRESSED, SO ITS NAME DROPS THE .Z
        basename = os.path.basename(filepath)
        if is_fits_frame(basename):
            filenameNoExtension = fits_frame_stem(basename)
            extension = FITS_SUFFIX
        else:
            filenameNoExtension, extension = os.path.splitext(basename)
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

        return absolute_path(path)

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

        # A SECOND CALL STARTS AGAIN, SO THE MASKS OF THE FIRST ARE NOT WRITTEN OR COUNTED TWICE
        self.cosmicRayMasks = []
        self.cosmicRayProduct = None

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
        self._write_cosmic_ray_mask_product()

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
            thisLamp = kw(f"LAMP{i + 1}")
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
        self._apply_measured_read_noise()

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
            raise TypeError(f"Input frames are a mix of arms: {arms}")
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
            raise TypeError("Input frames are a mix of binnings")

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
            raise TypeError(f"Input frames are a mix of readout speeds. {readSpeed}")

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
            raise TypeError(f"Input frames are a mix of gain {gain}")
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
                raise TypeError(f"Input frames are a mix of slit-width ({slitWidth})")

        return

    def _verify_single_readnoise(self):
        """*record the detector readnoise the input frames were read out with*

        Only raw frames (no ``PRO_CATG``) are read: a calibration product's read noise describes the frames it
        was built from, which may have used another readout mode or night (DY-1457).

        A mix of raw-frame readnoise values raises an exception; otherwise ``self.detectorParams["ron"]`` is set.
        """
        import numpy as np
        from astropy import units as u

        from soxspipe.commonutils.read_noise import in_electrons

        kw = self.kw

        def is_missing(value: object) -> bool:
            return value is None or value is np.ma.masked

        # HIERARCH ESO DET OUT1 RON - READOUT NOISE IN ELECTRONS
        summary = self.inputFrames.summary
        rawRon = [
            value
            for catg, value in zip(summary[kw("PRO_CATG")], summary[kw("RON")], strict=True)
            if is_missing(catg) and not is_missing(value)
        ]
        # UNIQUE VALUES IN ORDER OF FIRST APPEARANCE
        ron = list(dict.fromkeys(rawRon))
        self.defaultRon = in_electrons(self.detectorParams["ron"])
        self.rawHeaderRon = float(ron[0]) if len(ron) and ron[0] else None

        # MIXED NOISE
        if len(ron) > 1:
            self._report_verification_error(showSummary=True)
            raise TypeError(f"Input frames are a mix of readnoise. {ron}")
        if len(ron) and ron[0]:
            # UVB & VIS
            self.detectorParams["ron"] = ron[0] * u.electron
        else:
            # NIR
            self.log.debug(
                f"no raw frame carries a usable {kw('RON')}; using the detector-parameter read noise "
                f"({self.detectorParams['ron']} e-)"
            )
            self.detectorParams["ron"] = self.detectorParams["ron"] * u.electron

        return

    def _apply_measured_read_noise(self):
        """*set the read noise from the sample count of the SOXS NIR frames, or the master frame the SOF names*

        A SOXS NIR recipe looks the read noise up from the ``ESO DET NDSAMPLES`` of its own raw frames, with or
        without a master dark in the SOF, through ``_apply_ndsamples_read_noise`` (DY-2341). Every other case takes
        it from a master frame: UVB and VIS frames from the master bias (``MASTER_BIAS_<arm>``), X-Shooter NIR
        frames from the master dark (``MASTER_DARK_NIR``). A UVB/VIS SOF never takes its read noise from a master
        dark, and neither does an X-Shooter NIR recipe with ``usesMasterDarkReadNoise`` False (the nodding recipe,
        DY-2335). Without that master frame in the SOF, the nominal read noise stays. Otherwise
        ``self.detectorParams["ron"]`` is set to the single-frame read noise, which the prepared uncertainty and the
        bias-pedestal correction use. The master frame's noise reaches a detrended frame through
        ``ccdproc.subtract_bias`` or ``ccdproc.subtract_dark``, so it is not added here. The sky-subtraction noise
        model has no such propagation; ``self.readNoise.skyFit`` holds the value it needs (DY-1458, DY-1581).
        """
        from astropy import units as u
        from astropy.io import fits

        from soxspipe.commonutils.read_noise import SOURCE_MASTER_BIAS, SOURCE_MASTER_DARK, resolve_read_noise

        arm = self.arm.upper()
        if arm == "NIR" and self.settings.get("instrument", "soxs") == "soxs":
            self._apply_ndsamples_read_noise()
            return
        if arm in ("UVB", "VIS"):
            masterCatg, masterLabel, masterPlural = f"MASTER_BIAS_{arm}", SOURCE_MASTER_BIAS, "master biases"
        elif arm == "NIR" and self.usesMasterDarkReadNoise:
            masterCatg, masterLabel, masterPlural = "MASTER_DARK_NIR", SOURCE_MASTER_DARK, "master darks"
        else:
            return

        masterPath = self._select_master_frame(masterCatg, masterPlural)
        if masterPath is None:
            return

        try:
            masterHeader = fits.getheader(masterPath)
        except OSError as e:
            self.log.warning(f"The {masterLabel} header could not be read ({e}).")
            masterHeader = None

        self.readNoise = resolve_read_noise(
            masterHeader, self.rawHeaderRon, self.defaultRon, self.kw, self.log, masterLabel=masterLabel
        )
        # THE MASTER DARK'S QC RAW RON IS TAKEN AS THE NOISE OF ONE INTEGRATION. A MASTER BUILT FROM NDIT > 1 DARKS
        # (X-SHOOTER NIR ONLY; EVERY SOXS DARK IS NDIT 1) WOULD NEED IT SCALED BY SQRT(NDIT) OF THE MASTER
        # (DY-2317, DY-2322)
        self.detectorParams["ron"] = self.readNoise.singleFrame * u.electron

        return

    def _extraction_read_noise(self, frame: Any = None) -> float:
        """*the read noise of one integration (e-) that the Horne extraction uses for its sky-line signal-to-noise*

        **Key Arguments:**

        - ``frame`` -- the frame the extraction works on. Default *None*

        **Return:**

        - ``readNoise`` -- the prescan read noise of the frame when it has one (DY-2437). Otherwise the ``ron`` of
          the recipe's detector parameters, which the read-noise lookup (the ``ESO DET NDSAMPLES`` table for SOXS NIR)
          has set, as a plain float (DY-2323)
        """
        from soxspipe.commonutils.read_noise import in_electrons

        if frame is not None:
            return self._frame_read_noise(frame)
        return in_electrons(self.detectorParams["ron"])

    def _select_master_frame(self, masterCatg, masterPlural):
        """*pick the master frame of a product category from the SOF, warning when it names several*

        **Key Arguments:**

        - ``masterCatg`` -- the ``PRO_CATG`` of the master frame, e.g. ``MASTER_DARK_NIR``
        - ``masterPlural`` -- what several of them are called in the warning

        **Return:**

        - ``masterPath`` -- the path of the last master frame by name, so the choice does not depend on directory
          order, or *None* when the SOF names none
        """
        masterPaths = self.inputFrames.files_filtered(include_path=True, **{self.kw("PRO_CATG"): masterCatg})
        if not len(masterPaths):
            return None

        masterPath = sorted(masterPaths)[-1]
        if len(masterPaths) > 1:
            self.log.warning(
                f"The SOF names {len(masterPaths)} {masterPlural}; the read noise is taken from "
                f"{os.path.basename(masterPath)}."
            )
        return masterPath

    def _apply_ndsamples_read_noise(self):
        """*set the SOXS NIR read noise from the up-the-ramp sample count of the raw frames (DY-2341)*

        ``self.detectorParams["ron"]`` is the table value at the ``ESO DET NDSAMPLES`` of the raw frames, whether or
        not the SOF names a master dark; the dark's ``QC RAW RON`` sets nothing. ``self.readNoise.skyFit`` adds the
        noise of the master dark that the recipe subtracts. A nodding recipe (``usesMasterDarkReadNoise`` False)
        subtracts none, as A-B removes it, so its sky fit is the frame's own read noise. A SOF with no raw frame
        keeps the detector read noise.
        """
        import numpy as np
        from astropy import units as u

        from soxspipe.commonutils.read_noise import resolve_nir_read_noise

        kw = self.kw
        summary = self.inputFrames.summary
        # A SETTINGS FILE WITH ITS OWN summary-keys HAS NO SAMPLE-COUNT COLUMN: EVERY FRAME THEN HAS NO CARD
        counts = summary[kw("DET_NDSAMPLES")] if kw("DET_NDSAMPLES") in summary.colnames else [None] * len(summary)
        rawFrames = [
            (os.path.basename(str(name)), None if count is None or count is np.ma.masked else count)
            for name, catg, count in zip(summary["file"], summary[kw("PRO_CATG")], counts, strict=True)
            if catg is None or catg is np.ma.masked
        ]
        masterPath = None
        if self.usesMasterDarkReadNoise:
            masterPath = self._select_master_frame("MASTER_DARK_NIR", "master darks")

        self.readNoise = resolve_nir_read_noise(
            rawFrames,
            masterPath,
            self.detectorParams["read-noise-by-ndsamples"],
            self.defaultRon,
            self.workspaceRootPath,
            kw,
            self.log,
        )
        if self.readNoise is not None:
            self.detectorParams["ron"] = self.readNoise.singleFrame * u.electron

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

        A recipe whose stored status was not ``pass`` before this run (it failed or was interrupted) and which has
        now passed refreshes the session's SOF files once, so its dependents can run in the same ``reduce``. A
        recipe that passed before, or never ran, does not refresh. A recipe that passed before and now fails
        refreshes as a failure.

        **Key Arguments:**
        - ``forceFail`` -- force the recipe to be marked as failed in the DB

        **Usage**

        ```python
        recipe.clean_up(forceFail=True)
        ```
        """
        self.log.debug("starting the ``clean_up`` method")

        # FILTER QC TABLE ON qc_flag
        mask = self.qc["qc_flag"] == "fail"
        passToFail = False
        failedQcs = self.qc.loc[mask]
        failedQcs = failedQcs["qc_name"].tolist()
        if len(failedQcs):
            # KEEP AN EXPLICIT FAILURE MESSAGE SO IT STILL REACHES THE DATABASE
            if not isinstance(forceFail, str):
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

            # SELF.STATUS IS THE STATUS STORED BEFORE THIS RUN, READ BEFORE THE START-OF-RUN WRITE. A RECIPE THAT
            # WAS NOT 'PASS' (FAILED OR INTERRUPTED) AND HAS NOW PASSED UNBLOCKS ITS DEPENDENTS, SO THE SOF FILES
            # ARE REBUILT. A FIRST RUN (NO STORED STATUS) OR A REPEAT PASS CHANGES NOTHING DOWNSTREAM, SO IT IS
            # SKIPPED: THE REBUILD IS EXPENSIVE.
            isRecovery = self.status not in (None, "pass") and not forceFail
            isRegression = bool(forceFail) and self.status == "pass"
            if isRecovery or isRegression:
                from soxspipe.commonutils import data_organiser

                do = data_organiser(log=self.log, rootDir=self.workspaceRootPath)
                do.session_refresh(failure=isRegression)
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

        scratchPath = getattr(self, "scratchRoot", None) or scratch_root(self.settings)
        remove_scratch_tree(self.log, self.outDir, requiredParent=scratchPath)
        unregister_scratch_directory(self.outDir)

        if forceFail and isinstance(forceFail, str):
            self.log.error(f"\nRecipe marked as failed in the database. {forceFail}")

        elif forceFail:
            self.log.warning(
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

    def _measure_prescan_read_noise(self, frame: Any, filepath: Any) -> float | None:
        """*measure the read noise of a raw science frame from its prescan (DY-2437)*

        Only a UVB or VIS science frame of a science recipe is measured (see ``uses_prescan_read_noise``). Calibration
        frames in the SOF keep the read noise of the reduction. The frame must still be untrimmed. A frame whose prescan
        cannot be measured keeps the read noise of the reduction too, with a warning.

        **Key Arguments:**

        - ``frame`` -- the raw, untrimmed frame in ADU (CCDData object)
        - ``filepath`` -- the path of the raw frame, or the frame itself, for the log

        **Return:**

        - ``ron`` -- the read noise (e-) of the frame, or *None* when the frame is not measured or cannot be
        """
        from soxspipe.commonutils.prescan_read_noise import (
            PRESCAN_KEY,
            measure_prescan_read_noise,
            uses_prescan_read_noise,
        )
        from soxspipe.commonutils.read_noise import in_electrons

        dp = self.detectorParams
        dprType = frame.header.get(self.kw("DPR_TYPE"))
        if not uses_prescan_read_noise(self.arm, self.recipeName, dp, dprType):
            return None

        name = os.path.basename(filepath) if isinstance(filepath, str) else "a frame"
        baselineRon = in_electrons(dp["ron"])
        try:
            ron = measure_prescan_read_noise(frame.data, dp[PRESCAN_KEY], dp["binning"], dp["gain"])
        except ValueError as e:
            self.log.warning(
                f"The read noise of {name} cannot be measured from its prescan: {e}. "
                f"The frame takes the read noise of the reduction, {baselineRon:0.2f} e-."
            )
            return None

        self.log.info(f"The prescan read noise of {name} is {ron:0.2f} e- (the reduction's is {baselineRon:0.2f} e-).")
        return ron

    @staticmethod
    def _prescan_read_noise_of(frame: Any) -> float | None:
        """*the prescan read noise (e-) in the header of a frame or stack, or None without a usable one*

        The card holds anything ``float()`` accepts, except a boolean, that is a finite number above zero.
        The card is only ever written for a frame the science gate passed, so it needs no second check. The gate
        reads the DPR type, which the first product write scrubs from the frame header.
        """
        import math

        from soxspipe.commonutils.read_noise import SCIENCE_RON_CARD

        header = getattr(frame, "header", None)
        value = None if header is None else header.get(SCIENCE_RON_CARD)
        if value is None or isinstance(value, bool):
            return None
        try:
            ron = float(value)
        except (TypeError, ValueError):
            return None
        return ron if math.isfinite(ron) and ron > 0 else None

    def _frame_read_noise(self, frame: Any) -> float:
        """*the read noise (e-) of a frame: its prescan value in the header, else the read noise of the reduction*

        **Key Arguments:**

        - ``frame`` -- a prepared frame, or a stack of them (CCDData object)

        **Return:**

        - ``ron`` -- the ``QC SCI RON`` of the frame header, or ``self.detectorParams["ron"]`` without one
        """
        from soxspipe.commonutils.read_noise import in_electrons

        ron = self._prescan_read_noise_of(frame)
        if ron is not None:
            return ron
        return in_electrons(self.detectorParams["ron"])

    def _sky_fit_read_noise(self, frame: Any) -> float | None:
        """*the read noise (e-) the sky-fit noise model uses for a frame with the master bias subtracted*

        A frame with a prescan read noise adds the noise of the master bias in quadrature (``science_sky_fit_ron``).
        Any other frame uses what the master bias gave the reduction (``self.readNoise.skyFit``).

        **Key Arguments:**

        - ``frame`` -- the detrended frame (CCDData object)

        **Return:**

        - ``ron`` -- the read noise (e-), or *None* when the reduction has none, so that the frame header gives it
        """
        from soxspipe.commonutils.read_noise import science_sky_fit_ron

        ron = self._prescan_read_noise_of(frame)
        if ron is not None:
            return science_sky_fit_ron(ron, self.readNoise)
        return None if self.readNoise is None else self.readNoise.skyFit

    def _record_science_read_noise_qc(self, frame: Any) -> None:
        """*record the prescan read noise of a frame as the ``SCI RON`` QC, written to the product headers*

        A pass records the value of the frame (or stack, or difference) its products are built from, before it writes
        them, and replaces the row of an earlier pass. A frame with no prescan read noise records nothing.

        **Key Arguments:**

        - ``frame`` -- the frame the pass reduces (CCDData object)
        """
        ron = self._prescan_read_noise_of(frame)
        if ron is not None:
            self._set_science_read_noise_qc(ron)

    def _set_science_read_noise_qc(self, ron: float) -> None:
        """*record a read noise (e-) as the ``SCI RON`` QC, replacing the row of an earlier pass*

        **Key Arguments:**

        - ``ron`` -- the read noise (e-)
        """
        from soxspipe.commonutils import toolkit

        self.qc = self.qc[self.qc["qc_name"] != "SCI RON"]
        self.add_qc(
            qcName="SCI RON",
            qcValue=ron,
            qcComment="[e-] RON of the science frame, from its prescan",
            qcUnit="electrons",
            reductionDateUtc=toolkit.utcnow_string(),
            toHeader=True,
        )

    def _trim_frame(self, frame):
        """*return frame with pre-scan and overscan regions removed*

        **Key Arguments:**

        - ``frame`` -- the CCDData frame to be trimmed

        **Return:**

        - ``trimmed_frame`` -- the frame with its pre-scan and overscan regions removed (CCDData object)
        """
        self.log.debug("starting the ``_trim_frame`` method")

        import ccdproc

        from soxspipe.commonutils.prescan_read_noise import binned_window

        dp = self.detectorParams

        rs, re, cs, ce = binned_window(dp["science-pixels"], dp["binning"])

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

        from soxspipe.commonutils.atomic_fits import write_hdulist_atomically
        from soxspipe.commonutils.phase3 import add_qc_to_header, basic_header_scrubbing, sort_keywords

        # WRITE QCs TO HEADERS
        add_qc_to_header(frame.header, self.qc)

        if product:
            frame.header = basic_header_scrubbing(log=self.log, settings=self.settings, header=frame.header)
        frame.header = sort_keywords(log=self.log, header=frame.header)

        if not filename and self.sofName:
            filename = self.productStem + ".fits"
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
            frame.data[frame.mask] = 0

        HDUList = frame.to_hdu(hdu_mask="QUAL", hdu_uncertainty="ERRS", hdu_flags=None)
        HDUList[0].name = "FLUX"
        HDUList.verify("fix")
        # WRITE VIA A TEMPORARY NAME, SO A KILLED RUN NEVER LEAVES A PARTIAL FILE AT THE FINAL PATH (DY-1697)
        if product:
            write_hdulist_atomically(HDUList, filepath, output_verify="fix+warn", overwrite=overwrite, checksum=True)
        else:
            write_hdulist_atomically(HDUList, filepath, overwrite=overwrite, checksum=False)

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
        """*sigma-clip the individual input frames, adding the clipped pixels to each frame's mask*

        **Key Arguments:**

        - ``combiner`` -- the Combiner object holding the frames to clip. Its mask is replaced in place by the
          union of each frame's existing mask and its clip mask.
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

        # PIXELS ALREADY MASKED IN A FRAME TAKE NO PART IN THE CLIPPING STATISTICS AND STAY
        # MASKED IN THAT FRAME ONLY. THE CLIP MASK RETURNED INCLUDES THEM
        ## Reduce memory by avoiding an extra full-array copy during clipping
        combiner.data_arr.mask = sigma_clip(
            np.ma.array(
                np.asarray(combiner.data_arr.data, dtype=np.float32),
                mask=np.ma.getmaskarray(combiner.data_arr),
            ),
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
            percent = 100 * diff / (len(combiner.data_arr) * totalPixels)
            self.log.print(
                f"\tClipping found {diff} more rogue pixels in the set of all input frames (~{percent:0.2}% per-frame)"
            )
        return

    @staticmethod
    def _frame_mask(ccd):
        """*return the bad-pixel mask of a frame as a new full-shape boolean array*

        **Key Arguments:**

        - ``ccd`` -- the frame, a CCDData object. A frame with no mask has no flagged pixels.

        **Return:**

        - ``mask`` -- a boolean array the shape of the frame data, True where the pixel is flagged bad

        **Usage:**

        ```python
        combinedMask = self._frame_mask(ccds[0])
        ```
        """
        import numpy as np

        if ccd.mask is None:
            return np.zeros(ccd.data.shape, dtype=bool)
        return np.array(ccd.mask, dtype=bool)

    def _set_input_masks(self, combiner, combinedMask, ignore_input_masks):
        """*set the per-frame masks the clipping starts from*

        **Key Arguments:**

        - ``combiner`` -- the Combiner object holding the frames. Its per-frame masks are replaced in place;
          the input frames themselves are not changed.
        - ``combinedMask`` -- the pixels flagged bad in every input frame
        - ``ignore_input_masks`` -- drop the input masks, so that only the clip mask applies?

        **Usage:**

        ```python
        self._set_input_masks(combiner, combinedMask, ignore_input_masks=False)
        ```
        """
        import numpy as np

        if ignore_input_masks:
            combiner.data_arr.mask = np.zeros(combiner.data_arr.shape, dtype=bool)
            return

        # A PIXEL FLAGGED IN EVERY FRAME HAS NO FRAME TO STACK. IT KEEPS ITS DATA AND STAYS
        # FLAGGED IN THE COMBINED MASK, SO IT IS NOT LEFT WITH A ZERO FRAME COUNT. NON-FINITE
        # DATA IN THE FRAMES STILL GIVES A NAN MEAN THERE
        combiner.data_arr.mask = np.ma.getmaskarray(combiner.data_arr) & ~combinedMask

    def _stacked_header(self, ccds: list[Any]) -> Any:
        """*the header of a stack: that of the first frame, with the read noise of the whole stack (DY-2437)*

        A stack of science frames carries the root mean square of their prescan read noises (``QC SCI RON``), the
        read noise of the mean. A frame with no prescan value counts at the read noise of the reduction. Without a
        prescan value in any frame the first header is returned as it is, as before.

        **Key Arguments:**

        - ``ccds`` -- the frames that were stacked (CCDData objects)

        **Return:**

        - ``header`` -- the header of the stack. A copy when it carries a new read noise, so frame 1 keeps its own
        """
        from soxspipe.commonutils.fits_headers import set_header_card
        from soxspipe.commonutils.read_noise import SCIENCE_RON_CARD, rms_read_noise

        header = ccds[0].header
        if all(self._prescan_read_noise_of(ccd) is None for ccd in ccds):
            return header
        header = header.copy()
        set_header_card(
            header,
            SCIENCE_RON_CARD,
            rms_read_noise([self._frame_read_noise(ccd) for ccd in ccds]),
            "[e-] RON of the science frame, from its prescan",
        )
        return header

    def _mean_combine_frames(self, combiner, ccds, combinedMask):
        """*mean combine the clipped frames, and combine their error maps the same way*

        **Key Arguments:**

        - ``combiner`` -- the Combiner object holding the clipped frames. Its data is overwritten with the error maps.
        - ``ccds`` -- the input frames as a list of CCDData objects
        - ``combinedMask`` -- the pixels flagged bad in every input frame

        **Return:**

        - ``combined_frame`` -- the combined frame, with its uncertainty map attached

        The error map is the mean of the frames' errors over the square root of the number of frames that contributed
        at the pixel. For frames with different errors this is a little below the quadrature sum (by 5 % for errors
        of 4 and 8), which would be ``sqrt(sum(errors**2)) / N``. The mean is kept: it is what every stacked
        calibration product has used, and a quadrature sum would move them all.

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
        - ``ignore_input_masks`` -- ignore the input masks during clip and stacking? If *False*, a pixel flagged
          bad in a frame is left out of that frame only (its input mask is joined to its clip mask), and the
          uncertainty uses the number of frames that contributed at the pixel. If *True*, only the clip mask
          applies. Either way a pixel flagged in every frame stays flagged in the combined frame. The input
          frames' masks are never changed.
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

        # PIXELS FLAGGED BAD IN EVERY INPUT FRAME. THE CALLER'S MASKS ARE NEVER CHANGED
        combinedMask = self._frame_mask(ccds[0])
        for c in ccds[1:]:
            combinedMask &= self._frame_mask(c)

        # COMBINER OBJECT WILL FIRST GENERATE MASKS FOR INDIVIDUAL IMAGES VIA
        # CLIPPING AND THEN COMBINE THE IMAGES WITH THE METHOD SELECTED. PIXEL
        # MASKED IN ALL INDIVIDUAL IMAGES ARE MASK IN THE FINAL COMBINED IMAGE
        combiner = Combiner(ccds, dtype=np.float32)
        self._set_input_masks(combiner, combinedMask, ignore_input_masks)

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
        combined_frame.header = self._stacked_header(ccds)
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
        """*subtract a dark frame, scaling an X-Shooter dark first when its exposure time does not match*

        SOXS darks are never scaled by exposure time, because the SOXS NIR dark level does not
        change with exposure time. Darks of every other instrument (X-Shooter) are scaled whenever the two
        exposure times differ.

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

        # SOXS DARKS ARE NEVER SCALED BY EXPOSURE TIME: THE NIR DARK LEVEL IS FLAT ACROSS
        # EXPOSURE TIMES (DY-1263 RULING). X-SHOOTER DARKS ARE SCALED WHENEVER THE TIMES DIFFER.
        darkExptime = float(dark.header[kw("EXPTIME")])
        frameExptime = float(processedFrame.header[kw("EXPTIME")])
        exptimesMatch = abs(darkExptime - frameExptime) <= EXPTIME_MATCH_TOLERANCE_SECONDS
        if self.inst == "SOXS" or exptimesMatch:
            return ccdproc.subtract_dark(
                processedFrame,
                dark,
                exposure_time=kw("EXPTIME"),
                exposure_unit=u.second,
                add_keyword=None,
            )

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
            sofName=self.productStem,
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

    def _remove_bias_pedestal_from_uncertainty(self, frame, master_bias):
        """*re-derive the photon noise of a prepared frame from its bias-subtracted counts*

        **Key Arguments:**

        - ``frame`` -- the prepared frame, whose uncertainty was derived from its raw counts. CCDData object.
        - ``master_bias`` -- the master bias that is about to be subtracted from the frame. CCDData object.

        **Return:**

        - ``correctedFrame`` -- a copy of the frame whose uncertainty carries the photon noise of the bias-subtracted
          counts. CCDData object.

        The uncertainty map is built at preparation, before any master bias is available, so its photon term
        counts the electrons of the bias pedestal (about 650 e- in the VIS arm) as signal. This rescales each
        pixel's variance by the ratio of ``max(data - bias, 0) + ron^2`` to ``max(data, 0) + ron^2``. A ratio
        keeps the scaling right for a mean-combined stack, whose variance is the single-frame variance over the
        number of frames. ``ccdproc.subtract_bias`` then adds the variance of the master bias itself. Pixels
        that cosmic-ray cleaning changed after preparation are masked, so a ratio from their cleaned counts
        does not matter. A master dark still carries the pedestal variance of its own raw frames.

        **Usage:**

        ```python
        frame = self._remove_bias_pedestal_from_uncertainty(frame, master_bias)
        ```
        """
        import numpy as np
        from astropy.nddata import StdDevUncertainty

        if frame.uncertainty is None:
            return frame

        ronValue = self._frame_read_noise(frame)
        rawCounts = np.asarray(frame.data, dtype=np.float64)
        biasCounts = np.asarray(master_bias.data, dtype=np.float64)
        with np.errstate(invalid="ignore", divide="ignore"):
            varianceScale = (np.clip(rawCounts - biasCounts, 0, None) + ronValue**2) / (
                np.clip(rawCounts, 0, None) + ronValue**2
            )
        # LEAVE PIXELS WITH NO USABLE COUNTS UNCHANGED
        varianceScale = np.where(np.isfinite(varianceScale), varianceScale, 1.0)

        correctedFrame = frame.copy()
        correctedError = np.asarray(frame.uncertainty.array, dtype=np.float64) * np.sqrt(varianceScale)
        correctedFrame.uncertainty = StdDevUncertainty(correctedError.astype(np.float32), unit=frame.uncertainty.unit)
        return correctedFrame

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
        # UVB/VIS MASTER DARKS KEEP THE BIAS PEDESTAL (SOXS_MDARK SUBTRACTS NO BIAS), AS DO THE LAMP-OFF AND RAW DARK
        # FRAMES SOXS_MFLAT FALLS BACK TO, SO A BIAS PLUS A DARK WOULD SUBTRACT THE PEDESTAL TWICE (DY-1580)
        if self.arm.upper() in ("UVB", "VIS") and master_bias != False and dark != False:  # noqa: E712
            raise ValueError(
                f"{self.arm} frames cannot be detrended with both a master bias and a master dark: "
                f"{self.arm} master darks, lamp-off frames and raw darks still carry the bias pedestal, "
                "so this combination is not supported"
            )
        processedFrame = inputFrame

        if master_bias != False:  # noqa: E712
            processedFrame = self._remove_bias_pedestal_from_uncertainty(processedFrame, master_bias)
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
            # A FLAT PIXEL THAT IS ZERO, NEGATIVE OR NOT FINITE CANNOT FLAT-FIELD THE FRAME (DY-2279)
            unusableFlat = unusable_flat_pixel_mask(master_flat.data)
            processedFrame = ccdproc.flat_correct(
                processedFrame,
                neutralise_unusable_flat(master_flat, unusableFlat),
                norm_value=1.0,
                add_keyword=None,
            )
            processedFrame = flag_unusable_pixels(processedFrame, unusableFlat)
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

        cosmicRaySofName = self._register_cosmic_ray_product()

        if not self.verbose:
            # REMOVE COLUMN FROM DATA FRAME
            self.products.drop(columns=["file_path"], inplace=True)

        self.flag_poor_data()

        # REMOVE DUPLICATE ENTRIES IN COLUMN 'qc_name' AND KEEP THE LAST ENTRY
        self.qc = self.qc.drop_duplicates(subset=["qc_name", "qc_order"], keep="last")

        # SEND TO DATABASE
        self.qc["sof_name"] = self.sofName + ".sof"
        if cosmicRaySofName is not None:
            # THE MASK BELONGS TO THE WHOLE SOF, EVEN WHEN A PASS NAMED `<SOF>_F<N>` REPORTS FIRST
            self.qc.loc[self.qc["qc_name"] == "COSMIC RAY PIX NUM", "sof_name"] = cosmicRaySofName + ".sof"
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
                # THE OLD ROWS FOR THESE SOFS ARE DELETED IN THE SAME TRANSACTION AS THE INSERT, SO A READER NEVER
                # SEES THEM MISSING AND A FAILED INSERT KEEPS THEM
                self._dataframe_to_sqlite(self.qc[dbColumns], "quality_control", delete_where_in=("sof_name", sofNames))

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

        import pandas as pd

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

        # TEXT QCS (E.G. UG SHIFT SOURCE) BECOME NaN HERE, SO THEY ARE NEVER RANGE-FLAGGED
        numericValues = pd.to_numeric(self.qc["qc_value"], errors="coerce")
        for k, v in self.recipeSettings["qc-acceptable-ranges"].items():
            matchName = k.lower().replace("-", " ").replace("_", " ")
            mask = (
                (self.qc["qc_name"].str.lower() == matchName)
                & (self.qc["qc_order"] == "-1")
                & ((numericValues <= v[0]) | (numericValues >= v[1]))
            )
            self.qc.loc[mask, "qc_flag"] = "fail"
            # ADD GUARDRAIL VALUES TO QC TABLE FOR REPORTING PURPOSES
            self.qc.loc[(self.qc["qc_name"].str.lower() == matchName), "qc_value_min"] = v[0]
            self.qc.loc[(self.qc["qc_name"].str.lower() == matchName), "qc_value_max"] = v[1]

        self.log.debug("completed the ``flag_poor_data`` method")
        return

    def _measure_raw_frame_ron(self):
        """*measure the read-out noise in a single raw frame*

        The mask returned alongside the noise is the mask of the difference of the first two input frames and is
        needed by the caller to measure the master frame against the same pixels. At least two input frames are
        required.

        UVB and VIS (bias) noise is the clipped standard deviation of the difference of the first two input
        frames, so any fixed pattern the frames share cancels. NIR (dark) noise is measured by
        ``_measure_nir_dark_ron``, because the first two NIR darks carry a heavy tail and a first frame that is
        still settling.

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons
        - ``combinedMask`` -- the bad-pixel mask of the sigma-clipped frame difference

        **Usage:**

        ```python
        rawRon, combinedMask = self._measure_raw_frame_ron()
        ```
        """
        from itertools import islice

        ccds = self.inputFrames.ccds(
            ccd_kwargs={
                "hdu_uncertainty": "ERRS",
                "hdu_mask": "QUAL",
                "hdu_flags": "FLAGS",
                "key_uncertainty_type": "UTYPE",
            }
        )

        if self.arm == "NIR":
            return self._measure_nir_dark_ron(ccds)

        # ONLY THE FIRST TWO RAW FRAMES ARE NEEDED, SO NO MORE ARE LOADED
        raw_one, raw_two = islice(ccds, 2)
        return self._first_pair_ron(raw_one, raw_two)

    def _cross_check_nir_dark_ron(self, rawRon, header):
        """*warn when a SOXS NIR master dark's measured read noise is far from the table value at its sample count*

        The measured ``QC RAW RON`` is a QC value only; it never sets the read noise (DY-2341). It is compared with
        the table value at the dark's own ``ESO DET NDSAMPLES``, and a warning is logged when they differ by more than
        ``DARK_RON_TABLE_TOLERANCE`` (20%) of the table value. X-Shooter and the UVB and VIS arms are not checked, and
        neither is a dark whose sample count the table cannot serve.

        **Key Arguments:**

        - ``rawRon`` -- the read noise (e-) measured from the raw darks
        - ``header`` -- the header of the master dark

        **Usage:**

        ```python
        self._cross_check_nir_dark_ron(rawRon, frame.header)
        ```
        """
        from soxspipe.commonutils.read_noise import DARK_RON_TABLE_TOLERANCE, table_read_noise

        if self.arm != "NIR" or self.settings.get("instrument", "soxs") != "soxs":
            return

        ndsamples = header.get(self.kw("DET_NDSAMPLES"))
        tableRon = table_read_noise(ndsamples, self.detectorParams["read-noise-by-ndsamples"])
        if tableRon is None:
            self.log.info(
                f"The master dark has no usable {self.kw('DET_NDSAMPLES')}; its QC RAW RON is not cross-checked."
            )
            return

        deviation = abs(rawRon - tableRon) / tableRon
        if deviation > DARK_RON_TABLE_TOLERANCE:
            self.log.warning(
                f"The master dark QC RAW RON of {rawRon:0.2f} e- is {deviation:.0%} from the table value of "
                f"{tableRon:0.2f} e- at {int(ndsamples)} samples (limit {DARK_RON_TABLE_TOLERANCE:.0%}). "
                "The read noise of the reduction does not use it; check the raw darks."
            )

    def _first_pair_ron(self, raw_one, raw_two):
        """*measure the read-out noise and the bad-pixel mask from the difference of two raw frames*

        **Key Arguments:**

        - ``raw_one`` -- the first raw frame (CCDData object)
        - ``raw_two`` -- the second raw frame (CCDData object)

        **Return:**

        - ``rawRon`` -- the clipped standard deviation of the difference, divided by the square root of two (e-)
        - ``combinedMask`` -- the bad-pixel mask of the sigma-clipped frame difference

        **Usage:**

        ```python
        rawRon, combinedMask = self._first_pair_ron(raw_one, raw_two)
        ```
        """
        import math

        import numpy as np
        from astropy.stats import sigma_clip

        from soxspipe.commonutils import toolkit

        # SINGLE FRAME RON
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

    def _measure_nir_dark_ron(self, ccds):
        """*measure the NIR read-out noise as a robust median over consecutive raw-dark pairs*

        The noise of a pair is the MAD sigma of the difference of two consecutive raw darks, divided by the square
        root of two. The result is the median of the pair noises. The first dark is still settling, so with three
        or more darks it is left out of the pairs; with two darks the one pair is used. A pair with no spread (a
        fully masked difference or identical darks) is left out of the median. The darks are read one at a time, so
        only the previous dark is held in memory. The returned mask is that of the difference of the first two
        darks, as for the bias, so the master noise is measured against the same pixels.

        A ValueError is raised when no pair gives a usable noise, for example with a single dark.

        **Key Arguments:**

        - ``ccds`` -- the raw dark frames (CCDData objects), in order, as an iterable that is read once

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons
        - ``combinedMask`` -- the bad-pixel mask of the sigma-clipped difference of the first two darks

        **Usage:**

        ```python
        rawRon, combinedMask = self._measure_nir_dark_ron(ccds)
        ```
        """
        import math

        import numpy as np
        from astropy.stats import mad_std

        combinedMask = None
        pairRons = []
        previous = None
        for frame in ccds:
            if previous is not None:
                if combinedMask is None:
                    _firstPairRon, combinedMask = self._first_pair_ron(previous, frame)
                # MASKED ARITHMETIC TAKES THE UNION OF THE TWO MASKS
                difference = np.ma.array(frame.data, mask=frame.mask, dtype=np.float64) - np.ma.array(
                    previous.data, mask=previous.mask, dtype=np.float64
                )
                # A FULLY MASKED OR NaN DIFFERENCE GIVES 0.0, OR A MASKED CONSTANT ON SOME ASTROPY
                # VERSIONS (FILLED AS NaN HERE); BOTH ARE DROPPED BELOW
                pairSigma = float(np.ma.filled(mad_std(difference, ignore_nan=True), np.nan))
                pairRons.append(pairSigma / math.sqrt(2))
            previous = frame

        # THE FIRST DARK IS STILL SETTLING, SO PAIR 1-2 IS USED ONLY WHEN IT IS THE ONLY PAIR
        usedRons = [ron for ron in (pairRons[1:] or pairRons) if np.isfinite(ron) and ron > 0]
        if not usedRons:
            message = "The raw input frames appear to be corrupted. Cannot calculate the read-out noise. Please check the raw frames."  # noqa: E501
            raise ValueError(message)
        return float(np.median(usedRons)), combinedMask

    def qc_ron(
        self,
        frameType=False,
        frameName=False,
        masterFrame=False,
        masterRon=False,
    ):
        """*calculate the read-out-noise from bias/dark frames*

        For UVB and VIS bias frames the raw read-out-noise is measured on the difference
        of the first two raw frames, so fixed pattern noise cancels. For NIR dark frames
        it is the median, over consecutive raw-frame pairs, of the MAD sigma of the pair
        difference divided by the square root of two, with the first frame skipped when
        there are at least three. The ``MASTER RON`` mask is that of the frame 1-2
        difference in both cases. With fewer than two raw frames the raw noise cannot be
        measured, and no ``RAW RON`` QC is recorded.

        **Key Arguments:**

        - ``frameType`` -- the type of the frame for reporting QC values. Default *False*
        - ``frameName`` -- the name of the frame in human readable words. Default *False*
        - ``masterFrame`` -- the master frame (only makes sense to measure RON on master bias). Default *False*
        - ``masterRon`` -- if serendipitously calculated elsewhere don't recalculate. Default *False*

        **Return:**

        - ``rawRon`` -- raw read-out-noise in electrons, or *None* if fewer than two raw frames were supplied
        - ``masterRon`` -- combined read-out-noise in mbias, or *None* if there is no master noise

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

        rawRon = None
        combinedMask = None
        if len(self.inputFrames.files) > 1:
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
            # WITHOUT A RAW PAIR THERE IS NO RAW MASK, SO USE THE MASTER FRAME'S OWN
            if combinedMask is None:
                combinedMask = masterFrame.mask

            # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
            tmp = np.ma.array(masterFrame.data, mask=combinedMask)

            dmin, dmax, dmean, dstd = _image_stats(tmp)
            masterRon = float(dstd)

        # RECORD THE MASTER RON ROW ONCE, WHETHER MEASURED HERE OR SUPPLIED
        if masterRon:
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

    def _fill_masked_uncertainty(self, frame, fillValue):
        """*set the uncertainty of the masked pixels of a frame to a single value*

        **Key Arguments:**

        - ``frame`` -- the frame (CCDData object) whose uncertainty is updated in place
        - ``fillValue`` -- the uncertainty (e-) given to every masked pixel

        **Usage:**

        ```python
        self._fill_masked_uncertainty(frame, fillValue=rawRon)
        ```
        """
        import numpy as np

        frame.uncertainty = np.ma.array(
            frame.uncertainty.array,
            mask=frame.mask,
            fill_value=fillValue,
        ).filled()

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

        # THE EXPTIME LOOKUP RAISES KeyError FOR A FRAME WITHOUT THAT KEYWORD.
        # THE LINT CHARACTERIZATION TEST PINS THAT FAILURE.
        _ = frame.header[self.kw("EXPTIME")]

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
        import re

        from soxspipe.commonutils.fits_headers import set_header_card

        kw = self.kw

        # A HEADER STAMPED BEFORE KEEPS ITS RECORDS, AND A SHORTER LIST WOULD LEAVE THE LAST ONES OF THE OLD LIST
        for keyword in [key for key in frame.header if re.fullmatch(r"ESO PRO REC1 RAW\d+ (NAME|CATG)", key)]:
            del frame.header[keyword]

        if rawFrames:
            unmatched = [name for name in rawFrames if name not in set(tableData["filename"].values)]
            if unmatched:
                self.log.warning(
                    f"No raw frame of this recipe has the name {', '.join(unmatched)}, so the product header does "
                    "not list it in ESO PRO REC1 RAWi."
                )

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
                set_header_card(frame.header, f"ESO PRO REC1 RAW{iterator} NAME", f[:valueLen])
                set_header_card(frame.header, f"ESO PRO REC1 RAW{iterator} CATG", t)
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

        from soxspipe.commonutils.fits_headers import set_header_card

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
                set_header_card(frame.header, f"ESO PRO REC1 CAL{iterator} NAME", f[:valueLen])
                set_header_card(frame.header, f"ESO PRO REC1 CAL{iterator} CATG", c)
                set_header_card(frame.header, f"ESO PRO REC1 CAL{iterator} DATAMD5", compute_hash(p))
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
        from soxspipe.commonutils.fits_headers import set_header_card

        iterator = 1
        for k, v in self.recipeSettings.items():
            if k.lower() in ["uvb", "vis", "nir", "qc-acceptable-ranges"]:
                continue

            if not isinstance(v, dict):
                if isinstance(v, list):
                    v = ", ".join(map(str, v)).strip()

                set_header_card(frame.header, f"ESO PRO REC1 PARAM{iterator} NAME", k[:40])
                set_header_card(frame.header, f"ESO PRO REC1 PARAM{iterator} VALUE", v)
                iterator += 1
            else:
                for k2, v2 in v.items():
                    set_header_card(frame.header, f"ESO PRO REC1 PARAM{iterator} NAME", k2[:40])
                    set_header_card(frame.header, f"ESO PRO REC1 PARAM{iterator} VALUE", v2)
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

        from soxspipe.commonutils.fits_headers import set_header_card
        from soxspipe.commonutils.phase3 import pipe_id

        arm = self.arm
        kw = self.kw
        # UNUSED, AND ONE OF THE MODULE'S ELEVEN `F841` FINDINGS. DELETING IT IS DY-88'S.
        dp = self.detectorParams  # noqa: F841
        imageType = self.imageType
        if "FLAT" in imageType:
            imageType = "FLAT"

        set_header_card(frame.header, kw("SEQ_ARM").upper(), arm)
        set_header_card(frame.header, kw("PRO_TYPE").upper(), "REDUCED")

        # PROD CATG
        if imageType in ["BIAS", "DARK", "FLAT"]:
            productCategory = f"MASTER_{imageType}_{arm}".replace("QLAMP", "LAMP").replace("DLAMP", "LAMP")
            set_header_card(frame.header, kw("PRO_CATG"), productCategory)
            set_header_card(frame.header, kw("PRO_TECH"), "IMAGE")

        # SPEC FORMAT TO PANDAS DATAFRAME
        tableData = self.rawFrames.to_pandas()

        tableData["filename"] = tableData["filename"].str.replace("_pre", "")
        tableData["tag"] = tableData["TYPE"] + "_" + tableData["ARM"]

        self._stamp_raw_frame_records(frame, tableData, rawFrames)
        self._stamp_calibration_frame_records(frame, tableData)
        self._stamp_recipe_parameter_records(frame)

        # SOXSPIPE VERSION
        set_header_card(frame.header, "ESO PRO REC1 PIPE ID", pipe_id())

        # RECIPE
        if self.recipeName:
            set_header_card(frame.header, "ESO PRO REC1 ID", self.recipeName)

        # from tabulate import tabulate
        # print(tabulate(tableData, headers='keys', tablefmt='github'))

        self.log.debug("completed the ``update_fits_keywords`` method")
        return

    @staticmethod
    def _concat_qc_tables(qcTables):
        """*join QC tables into one, skipping empty tables so their all-NA columns do not decide the dtypes*

        **Key Arguments:**

        - ``qcTables`` -- the QC tables to join

        **Return:**

        - ``qcTable`` -- the joined QC table. Empty when every table is empty, or when there are no tables
        """
        import pandas as pd

        if not qcTables:
            return pd.DataFrame()
        return pd.concat([table for table in qcTables if len(table.index)] or qcTables, ignore_index=True)

    def _stacking_override_requested(self):
        """*check whether the user asked on the command line to extract each frame or cycle on its own (``--no-stack``)*

        A ``stack_before_extraction`` value in the settings file does not count. Only the command-line option sets the
        ``no-stack`` setting, and only then may a standard star be extracted frame by frame or cycle by cycle.

        **Return:**

        - ``isRequested`` -- True when ``--no-stack`` was given
        """
        return bool(self.settings.get("no-stack", False))

    def get_recipe_settings(self):
        """*get the recipe and arm specific settings*

        The settings are a deep copy of the recipe's block, with the arm's values merged in. Changing them does not
        change ``self.settings``, which the reducer shares between the recipes it runs in one process.

        **Return:**

        - ``recipeSettings`` -- the recipe specific settings, or False when the recipe has no settings block

        **Usage:**

        ```python
        usage code
        ```
        """
        import copy

        self.log.debug("starting the ``get_recipe_settings`` method")

        recipeSettings = False
        if self.recipeName:
            recipeSettings = copy.deepcopy(self.settings[self.recipeName])
        if recipeSettings and self.arm and self.arm.lower() in recipeSettings:
            for k, v in recipeSettings[self.arm.lower()].items():
                recipeSettings[k] = v

        self.log.debug("completed the ``get_recipe_settings`` method")
        return recipeSettings

    def _dataframe_to_sqlite(self, dataframe, table_name, replace=False, delete_where_in=None):
        """*write a dataframe to a database table in one transaction, retrying a failed write*

        The number of attempts is `sqlite_dataframe_writer.MAX_ATTEMPTS`.

        **Key Arguments:**

        - ``dataframe`` -- the dataframe containing the rows to insert
        - ``table_name`` -- the name of the database table to insert into
        - ``replace`` -- if True, delete the table's existing rows in the same transaction; otherwise append. Default:
          False
        - ``delete_where_in`` -- a ``(column, values)`` pair: delete the rows whose column value is in ``values`` in
          the same transaction. Default: None

        **Raises:**

        - Exception if every attempt fails.
        - ``SqliteTransactionOpenError`` if the connection is already inside a transaction.
        - ``UnsafeSqlIdentifierError`` if ``table_name`` or the delete column fails the safe-identifier grammar.
        """
        write_dataframe_to_sqlite(
            self.conn,
            self.log,
            dataframe,
            table_name,
            replace=replace,
            method="multi",
            deleteWhereIn=delete_where_in,
        )

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
