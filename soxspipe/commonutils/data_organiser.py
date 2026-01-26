#!/usr/bin/env python
# encoding: utf-8
"""
*The SOXSPIPE Data Organiser*

Author
: David Young

Date Created
: March  9, 2023
"""
from line_profiler import profile
from fundamentals import tools
from builtins import object
import sys
import os
from soxspipe.commonutils import uncompress
from soxspipe.commonutils.toolkit import get_calibrations_path

os.environ["TERM"] = "vt100"


class data_organiser(object):
    """
    *The `soxspipe` Data Organiser*

    **Key Arguments:**

    - ``log`` -- logger
    - ``rootDir`` -- the root directory of the data to process
    - ``vlt`` -- prepare the workspace using the standard vlt /data directory

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/master/initialisation.html).

    To initiate a data_organiser object, use the following:

    ```python
    from soxspipe.commonutils import data_organiser
    do = data_organiser(
        log=log,
        rootDir="/path/to/workspace/root/"
    )
    do.prepare()
    ```
    """

    def __init__(self, log, rootDir, vlt=False):
        from os.path import expanduser
        import codecs
        from fundamentals.logs import emptyLogger
        import warnings
        from astropy.utils.exceptions import AstropyWarning
        import sqlite3 as sql

        warnings.simplefilter("ignore", AstropyWarning)

        self.vlt = vlt
        self.PAE = False
        log.debug("instantiating a new 'data_organiser' object")
        self.log = log

        # MAKE RELATIVE HOME PATH ABSOLUTE
        if rootDir[0] == "~":
            home = expanduser("~")
            directory = directory.replace("~", home)

        self.rootDir = rootDir
        self.rawDir = rootDir + "/raw"
        self.miscDir = rootDir + "/misc"
        self.sessionsDir = rootDir + "/sessions"

        if self.vlt:
            self.vltReduced = self.use_vlt_environment_folders()

        # SESSION ID PLACEHOLDER FILE
        self.sessionIdFile = self.sessionsDir + "/.sessionid"
        exists = os.path.exists(self.sessionIdFile)
        if exists:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                sessionId = readFile.read()
                self.sessionPath = self.sessionsDir + "/" + sessionId
                self.sessionId = sessionId

        # DATABASE FILE
        self.rootDbPath = rootDir + "/soxspipe.db"

        # RETURN HERE: add these to yaml file
        # A LIST OF FITS HEADER KEYWORDS LOOKUP KEYS. THESE KEYWORDS WILL BE LIFTED FROM ALL FITS FILES
        self.keyword_lookups = [
            "MJDOBS",
            "DATE_OBS",
            "SEQ_ARM",
            "DPR_CATG",
            "DPR_TECH",
            "DPR_TYPE",
            "PRO_CATG",
            "PRO_TECH",
            "PRO_TYPE",
            "EXPTIME",
            "WIN_BINX",
            "WIN_BINY",
            "DET_READ_SPEED",
            "SLIT_UVB",
            "SLIT_VIS",
            "SLIT_NIR",
            "LAMP1",
            "LAMP2",
            "LAMP3",
            "LAMP4",
            "LAMP5",
            "LAMP6",
            "LAMP7",
            "DET_READ_TYPE",
            "CONAD",
            "RON",
            "OBS_ID",
            "OBS_NAME",
            "NAXIS",
            "OBJECT",
            "TPL_ID",
            "INSTRUME",
            "ABSROT",
            "EXPTIME2",
            "TPL_NAME",
            "TPL_NEXP",
            "TPL_EXPNO",
            "ACFW_ID",
            "RA",
            "DEC",
            "DET",
            "SWSIM_NIR",
            "SWSIM_VIS",
            "SWSIM_UVB",
        ]

        # THE MINIMUM SET OF KEYWORD WE EVER WANT RETURNED
        # USED IN GROUPING RAW FRAMES
        # RETURN HERE: add these to yaml file
        self.rawFrameGroupKeywords = [
            "file",
            "eso seq arm",
            "eso dpr catg",
            "eso dpr type",
            "eso dpr tech",
            "eso pro catg",
            "eso pro tech",
            "eso pro type",
            "eso obs id",
            "eso obs name",
            "exptime",
            "binning",
            "rospeed",
            "slit",
            "slitmask",
            "lamp",
            "night start date",
            "night start mjd",
            "mjd-date",
            "mjd-obs",
            "date-obs",
            "object",
            "template",
            "instrume",
            "absrot",
            "eso tpl name",
            "eso tpl nexp",
            "eso tpl expno",
            "filter",
            "ra",
            "dec",
            "simulation",
        ]

        # THIS TYPE MAP WILL BE USED TO GROUP SET OF FILES TOGETHER
        self.typeMapXSH = {
            "bias": [
                {"tech": None, "slitmask": None, "recipe": "mbias"}
            ],  # XSH/SOXS BIAS CAN BE DEFINED WITH JUST DPR TYPE
            "dark": [
                {"tech": None, "slitmask": None, "recipe": "mdark"}
            ],  # XSH/SOXS DARK CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,fmtchk": [
                {"tech": None, "slitmask": None, "recipe": "disp_sol"}
            ],  # XSH disp_sol CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,orderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,dorderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,qorderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,flat": [
                {"tech": None, "slitmask": None, "recipe": "mflat"}
            ],  # XSH flats CAN BE DEFINED WITH JUST DPR TYPE
            "flat,lamp": [
                {"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"},
            ],
            "flat,dome": [{"tech": None, "slitmask": None, "recipe": "mflat"}],
            "lamp,dflat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],
            "lamp,qflat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],
            "lamp,wave": [
                {"tech": ["echelle,multi-pinhole", "image"], "slitmask": None, "recipe": "spat_sol"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": None, "recipe": "disp_sol"},
            ],
            "wave,lamp": [
                {"tech": ["echelle,multi-pinhole", "image"], "slitmask": ["MPH"], "recipe": "spat_sol"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "disp_sol"},
            ],
            "object": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std,flux": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std,telluric": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std,sky": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
        }

        # THIS TYPE MAP WILL BE USED TO GROUP SET OF FILES TOGETHER
        self.typeMapSOXS = {
            "bias": [
                {"tech": None, "slitmask": None, "recipe": "mbias"}
            ],  # XSH/SOXS BIAS CAN BE DEFINED WITH JUST DPR TYPE
            "dark": [
                {"tech": None, "slitmask": None, "recipe": "mdark"}
            ],  # XSH/SOXS DARK CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,fmtchk": [
                {"tech": None, "slitmask": None, "recipe": "disp_sol"}
            ],  # XSH disp_sol CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,orderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,dorderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,qorderdef": [
                {"tech": None, "slitmask": None, "recipe": "order_centres"}
            ],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,flat": [
                {"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"},
            ],
            "flat,lamp": [
                {"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"},
            ],
            "dome,flat": [{"tech": None, "slitmask": ["SLIT"], "recipe": "mflat"}],
            "lamp,dflat": [
                {"tech": None, "slitmask": ["SLIT"], "recipe": "mflat"},
                {"tech": None, "slitmask": ["PH"], "recipe": "order_centres"},
            ],
            "lamp,qflat": [
                {"tech": None, "slitmask": ["SLIT"], "recipe": "mflat"},
                {"tech": None, "slitmask": ["PH"], "recipe": "order_centres"},
            ],
            "lamp,wave": [
                {"tech": ["echelle,multi-pinhole", "image"], "slitmask": ["MPH"], "recipe": "spat_sol"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "disp_sol"},
            ],
            "wave,lamp": [
                {"tech": ["echelle,multi-pinhole", "image"], "slitmask": ["MPH"], "recipe": "spat_sol"},
                {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "disp_sol"},
            ],
            "object": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "object,async": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std,flux": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
            "std,telluric": [
                {"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"},
                {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"},
                {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"},
            ],
        }

        # THIS PRODUCT MAP IS USED TO PREDICT THE PRODUCTS THAT WILL RESULT FROM REDUCING EACH SOF
        # LIST ORDER: ['pro type kw', 'pro tech kw', 'pro catg kw', "pixels/table", "find in sof", "replace/suffix in sof", "recipe"]
        self.productMap = {
            "mbias": {"ALL": [["REDUCED", "IMAGE", "MASTER_BIAS", "PIXELS", None, None, "soxs-mbias"]]},
            "mdark": {"ALL": [["REDUCED", "IMAGE", "MASTER_DARK", "PIXELS", None, None, "soxs-mdark"]]},
            "disp_sol": {
                "ALL": [["REDUCED", "ECHELLE,PINHOLE", "DISP_TAB", "TABLE", None, None, "soxs-disp-solution"]]
            },
            "order_centres": {
                "ALL": [["REDUCED", "ECHELLE,SLIT", "ORDER_TAB", "TABLE", None, None, "soxs-order-centres"]]
            },
            "mflat": {
                "ALL": [
                    ["REDUCED", "ECHELLE,SLIT", "MASTER_FLAT", "PIXELS", None, None, "soxs-mflat"],
                    ["REDUCED", "ECHELLE,SLIT", "ORDER_TAB", "TABLE", "MFLAT", "OLOC", "soxs-mflat"],
                ]
            },
            "spat_sol": {
                "ALL": [
                    ["REDUCED", "ECHELLE,PINHOLE", "DISP_TAB", "TABLE", None, None, "soxs-spatial-solution"],
                    [
                        "REDUCED",
                        "ECHELLE,PINHOLE",
                        "DISP_IMAGE",
                        "PIXELS",
                        ".fits",
                        "_IMAGE.fits",
                        "soxs-spatial-solution",
                    ],
                ]
            },
            "stare": {
                "ALL": [["REDUCED", "ECHELLE,SLIT,STARE", "OBJECT_TAB", "TABLE", None, None, "soxs-stare"]],
                "STD,FLUX": [
                    ["REDUCED", "ECHELLE,SLIT,STARE", "RESP_TAB", "TABLE", ".fits", "_RESP.fits", "soxs-stare"]
                ],
            },
            "nod": {
                "ALL": [["REDUCED", "ECHELLE,SLIT,NODDING", "OBJECT_TAB", "TABLE", None, None, "soxs-nod"]],
                "STD,FLUX": [
                    ["REDUCED", "ECHELLE,SLIT,NODDING", "RESP_TAB", "TABLE", ".fits", "_RESP.fits", "soxs-nod"]
                ],
            },
        }

        self.proKeywords = ["eso pro type", "eso pro tech", "eso pro catg"]

        # THESE ARE KEYS WE NEED TO FILTER ON, AND SO NEED TO CREATE ASTROPY TABLE INDEXES
        self.filterKeywords = [
            "eso seq arm",
            "eso dpr catg",
            "eso dpr tech",
            "eso dpr type",
            "eso pro catg",
            "eso pro tech",
            "eso pro type",
            "exptime",
            "rospeed",
            "slit",
            "slitmask",
            "binning",
            "mjd-date",
            "night start mjd",
            "night start date",
            "instrume",
            "lamp",
            "template",
            "eso obs name",
            "eso obs id",
            "eso tpl name",
            "filter",
            "object",
        ]

        self.productFilterKeywords = [
            "eso seq arm",
            "eso pro catg",
            "eso pro tech",
            "eso pro type",
        ]

        # THIS IS THE ORDER TO PROCESS THE FRAME TYPES
        self.reductionOrder = [
            "BIAS",
            "DARK",
            "LAMP,FMTCHK",
            "LAMP,ORDERDEF",
            "LAMP,DORDERDEF",
            "LAMP,QORDERDEF",
            "LAMP,FLAT",
            "FLAT,LAMP",
            "DOME,FLAT",
            "LAMP,DFLAT",
            "LAMP,QFLAT",
            "WAVE,LAMP",
            "LAMP,WAVE",
            "STD,FLUX",
            "STD",
            "STD,TELLURIC",
            "OBJECT",
            "OBJECT,ASYNC",
        ]

        # THIS IS THE ORDER THE RECIPES NEED TO BE RUN IN (MAKE SURE THE REDUCTION SCRIPT HAS RECIPES IN THE CORRECT ORDER)
        self.recipeOrder = [
            "mbias",
            "mdark",
            "disp_sol",
            "order_centres",
            "mflat",
            "spat_sol",
            "stare",
            "nod",
            "offset",
        ]

        # DECOMPRESS .Z FILES
        from soxspipe.commonutils import uncompress

        uncompress(log=self.log, directory=self.rootDir)

        exists = os.path.exists(self.rootDbPath)
        if exists:
            self.conn, reset = self._get_or_create_db_connection()
        else:
            self.conn = None

        return None

    def prepare(self, refresh=False, report=True):
        """*Prepare the workspace for data reduction by generating all SOF files and reduction scripts.*

        **Key Arguments:**
        - ``refresh`` -- trigger a complete refresh the workspace during preparation (delete database and do a complete prepare)
        """
        self.log.debug("starting the ``prepare`` method")
        import codecs
        import glob

        if refresh:
            # DELETE THE SQLITE DATABASE IF IT EXISTS
            exists = os.path.exists(self.rootDbPath)
            self.conn = None
            if exists:
                os.remove(self.rootDbPath)
                print("The existing database has been removed to allow a complete refresh of the workspace.")
                try:
                    os.remove(self.rootDbPath + "-shm")
                except:
                    pass
                try:
                    os.remove(self.rootDbPath + "-wal")
                except:
                    pass
            # DELETE ALL ERROR LOG AND SOF FILES
            for root, dirs, files in os.walk(os.path.abspath(self.rootDir)):
                for file in files:
                    if file.endswith("ERROR.log") or file.endswith(".sof"):
                        os.remove(os.path.join(root, file))

        self._select_instrument()

        # TEST FITS FILES OR raw_frames DIRECT EXISTS
        fitsExist = self._fits_files_exist()
        # EXIST IF NO FITS FILES EXIST - SOME PROTECTION AGAINST MOVING USER FILES IF THEY MAKE A MISTAKE PREPARE A WORKSPACE IN THE WRONG LOCATION
        if fitsExist == False:
            print("There are no FITS files in this directory. Please add your data before running `soxspipe prep`")
            sys.exit()
            return None

        # MK RAW FRAME DIRECTORY
        if not os.path.exists(self.rawDir):
            os.makedirs(self.rawDir)

        # TEST FOR SQLITE DATABASE - ADD IF MISSING
        self.conn, reset = self._get_or_create_db_connection()

        # MK SESSION DIRECTORY
        if not os.path.exists(self.sessionsDir):
            os.makedirs(self.sessionsDir)

        basename = os.path.basename(self.rootDir)
        print(f"PREPARING THE `{basename}` WORKSPACE FOR DATA-REDUCTION")

        self._sync_raw_frames()
        self._move_misc_files()

        # IF SESSION ID FILE DOES NOT EXIST, CREATE A NEW SESSION
        # OTHERWISE USE CURRENT SESSION
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            sessionId = self.session_create(sessionId="base")
            self.sessionId = sessionId
        else:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                sessionId = readFile.read()
                self.sessionPath = self.sessionsDir + "/" + sessionId

        # GET SETTINGS
        settingsPath = self.sessionPath + "/soxspipe.yaml"
        su = tools(
            arguments={"<workspaceDirectory>": self.sessionPath, "settingsFile": settingsPath},
            docString=False,
            logLevel="WARNING",
            options_first=False,
            projectName="soxspipe",
            createLogger=False,
        )
        arguments, self.settings, replacedLog, dbConn = su.setup()

        self._flag_files_to_ignore()
        self.build_sof_files()

        if report:

            print(f"\nTHE `{basename}` WORKSPACE FOR HAS BEEN PREPARED FOR DATA-REDUCTION\n")
            print(f"In this workspace you will find:\n")
            print(f"   - `misc/`: a lost-and-found archive of non-fits files")
            print(f"   - `qc/`: nested folders, ordered by date, containing quality-control plots and tables.")
            print(f"   - `{self.rawDir}/`: nested folders, ordered by date, containing raw-frames.")
            print(f"   - `sessions/`: directory of data-reduction sessions")
            print(f"   - `sof/`: the set-of-files (sof) files required for each reduction step")
            print(f"   - `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete")
            print(f"   - `reductions/`: nested folders, ordered by date, containing reduced data.\n")

            incompleteSets = self.get_incomplete_raw_frames_set()
            if len(incompleteSets.index):
                from tabulate import tabulate

                print("SOME CALIBRATION FRAMES ARE NOT PRESENT FOR THE FOLLOWING DATA SETS AND THEY CANNOT BE REDUCED:")
                print(tabulate(incompleteSets, headers="keys", tablefmt="psql", showindex=False))

            self.conn.close()

        self.log.debug("completed the ``prepare`` method")
        return None

    def list_obs(self):
        """*list all observation names and IDs in the current workspace*"""
        import pandas as pd

        self.log.debug("starting the ``list_obs`` method")

        import pandas as pd

        query = 'select distinct "eso seq arm", "night start date", "night start mjd", "eso obs name","eso obs id" from raw_frame_sets where complete = 1 and recipe in ("nod-obj","stare-obj","offset-obj") and `eso dpr type` like "%OBJECT%" order by "mjd-obs" asc'
        obsDf = pd.read_sql(query, con=self.conn)

        if len(obsDf.index):
            from tabulate import tabulate

            print(f"THE CURRENT WORKSPACE CONTAINS {len(obsDf.index)} SCIENCE OBSERVATION BLOCKS:")
            print(tabulate(obsDf, headers="keys", tablefmt="psql", showindex=False))

        self.log.debug("completed the ``list_obs`` method")
        return obsDf

    def list_sofs(self):
        """*list all science object SOF files in the current workspace*"""
        import pandas as pd

        self.log.debug("starting the ``list_sofs`` method")

        query = 'select distinct "eso seq arm", "night start date", round("mjd-obs",1) as "mjd-obs", "eso obs name","eso obs id", "slit", "sof" from raw_frame_sets where complete = 1 and recipe in ("nod-obj","stare-obj","offset-obj") and `eso dpr type` like "%OBJECT%" order by "mjd-obs" asc'
        sofDf = pd.read_sql(query, con=self.conn)

        if len(sofDf.index):
            from tabulate import tabulate

            print(f"THE CURRENT WORKSPACE CONTAINS {len(sofDf.index)} SCIENCE SOF FILES:")
            print(tabulate(sofDf, headers="keys", tablefmt="psql", showindex=False))

        self.log.debug("completed the ``list_sofs`` method")
        return sofDf

    def _sync_raw_frames(self, skipSqlSync=False):
        """*sync the raw frames between the project folder and the database*

        **Key Arguments:**

        - ``skipSqlSync`` -- skip the SQL db sync (used only in secondary clean-up scan)

        **Return:**

        - None

        **Usage:**

        ```python
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            rootDir="/path/to/root/folder/"
        )
        do._sync_raw_frames()
        ```
        """
        self.log.debug("starting the ``_sync_raw_frames`` method")

        import shutil
        import pandas as pd
        import re

        remainingFiles = 1
        firstPass = True

        while remainingFiles > 0:
            # GENERATE AN ASTROPY TABLE OF FITS FRAMES WITH ALL INDEXES NEEDED
            rawFrames, fitsPaths, remainingFiles = self._create_directory_table(
                pathToDirectory=self.rootDir, filterKeys=self.filterKeywords
            )
            if not remainingFiles:
                remainingFiles = 0
            elif remainingFiles > 0:
                if firstPass:
                    firstPass = False
                else:
                    sys.stdout.flush()
                    sys.stdout.write("\x1b[1A\x1b[2K")
                print(remainingFiles, "FITS files remaining to be indexed")

            if fitsPaths:

                conn = self.conn
                knownRawFrames = pd.read_sql("SELECT * FROM raw_frames", con=conn)

                # SPLIT INTO RAW, REDUCED PIXELS, REDUCED TABLES
                rawFrames = self._populate_raw_frames_extra_columns(rawFrames)
                # FIND AND REMOVE DUPLICATE FILES
                if len(rawFrames.index):
                    rawFrames["filepath"] = f"./raw/" + rawFrames["mjd-date"] + "/" + rawFrames["file"]
                    # FIND AND REMOVE DUPLICATE FILES
                    matchedFiles = pd.merge(rawFrames, knownRawFrames, on=["file", "eso dpr tech"], how="inner")
                    if len(matchedFiles.index):
                        for file in matchedFiles["file"]:
                            try:
                                os.remove(file)
                            except:
                                pass
                        # FIND RECORDS IN THE FILE SYSTEM NOT YET IN THE DATABASE
                        rawFrames = rawFrames[
                            ~rawFrames.set_index(["file", "eso dpr tech"]).index.isin(
                                knownRawFrames.set_index(["file", "eso dpr tech"]).index
                            )
                        ]

                # ADD THE NEWLY FOUND FRAMES TO THE DATABASE
                databaseDeletes = []
                if len(rawFrames.index):

                    # NOW MAKE FILEPATHS RELATIVE TO THE rawDir
                    rawFrames["filepath"] = f"./raw/" + rawFrames["mjd-date"] + "/" + rawFrames["file"]
                    rawFrames.replace(["--", -99.99], None).to_sql(
                        "raw_frames", con=self.conn, index=False, if_exists="append"
                    )

                    # NOW MAKE THE DESTINATION FILEPATHS ABSOLUTE
                    rawFrames["filepath"] = f"{self.rawDir}/" + rawFrames["mjd-date"] + "/" + rawFrames["file"]

                    # MOVE THE FILES TO THE CORRECT LOCATION
                    filepaths = rawFrames["filepath"]
                    filenames = rawFrames["file"]
                    for p, n in zip(filepaths, filenames):
                        parentDirectory = os.path.dirname(p)
                        if not os.path.exists(parentDirectory):
                            # Recursively create missing directories
                            if not os.path.exists(parentDirectory):
                                os.makedirs(parentDirectory)
                        if os.path.exists(self.rootDir + "/" + n):
                            realSource = os.path.realpath(self.rootDir + "/" + n)
                            realDest = os.path.realpath(p)

                            matchObject = re.match(r".*?(raw\/\d{4}-\d{2}-\d{2}.*)", realSource)

                            if matchObject:
                                # FILE NOT WHERE THEY SHOULD BE - DELETE FROM DATABASE
                                databaseDeletes.append(realDest)
                            elif realSource != realDest:
                                shutil.move(realSource, realDest)
                            if os.path.islink(self.rootDir + "/" + n):
                                os.remove(self.rootDir + "/" + n)

                if len(databaseDeletes):
                    databaseDeletes = (", ").join(databaseDeletes)
                    c = self.conn.cursor()
                    sqlQuery = f'delete from raw_frames where filepath in ("{databaseDeletes}");'
                    c.execute(sqlQuery)
                    c.close()

        if not skipSqlSync:
            self._sync_sql_table_to_directory(self.rawDir, "raw_frames", recursive=False)

        self.log.debug("completed the ``_sync_raw_frames`` method")
        return None

    def _create_directory_table(self, pathToDirectory, filterKeys, limit=10000):
        """*create an astropy table based on the contents of a directory*

        **Key Arguments:**

        - `log` -- logger
        - `pathToDirectory` -- path to the directory containing the FITS frames
        - `filterKeys` -- these are the keywords we want to filter on later
        - `limit` -- maximum number of files to process in one go (to avoid memory issues)

        **Return**

        - `rawFrames` -- the primary dataframe table listing all FITS files in the directory (including indexes on `filterKeys` columns)
        - `fitsPaths` -- a simple list of all FITS file paths
        - `remainingFiles` -- number of remaining files not included in the table due to the limit

        **Usage:**

        ```python
        # GENERATE AN ASTROPY TABLES OF FITS FRAMES WITH ALL INDEXES NEEDED
        rawFrames, fitsPaths, remainingFiles = _create_directory_table(
            log=log,
            pathToDirectory="/my/directory/path",
            keys=["file","mjd-obs", "exptime","cdelt1", "cdelt2"],
            filterKeys=["mjd-obs","exptime"]
        )
        ```
        """
        self.log.debug("starting the ``_create_directory_table`` function")

        from ccdproc import ImageFileCollection
        from astropy.time import Time
        import numpy as np
        import pandas as pd
        from soxspipe.commonutils import keyword_lookup

        # GENERATE A LIST OF FITS FILE PATHS
        fitsPaths = []
        fitsNames = []
        for entry in os.scandir(pathToDirectory):
            if (
                not entry.name.startswith(".")
                and entry.is_file()
                and (os.path.splitext(entry.name)[1] == ".fits" or ".fits.Z" in entry.name)
            ):
                fitsPaths.append(entry.path)
                fitsNames.append(entry.name)

        remainingFiles = max(0, len(fitsPaths) - limit)
        fitsPaths = fitsPaths[:limit]
        fitsNames = fitsNames[:limit]

        recursive = False

        if len(fitsPaths) == 0:
            return None, None, None

        # INSTRUMENT CHECK
        if recursive:
            allFrames = ImageFileCollection(filenames=fitsPaths[:3], keywords=["instrume"])
        else:
            allFrames = ImageFileCollection(location=pathToDirectory, filenames=fitsNames[:3], keywords=["instrume"])

        tmpTable = allFrames.summary
        tmpTable["instrume"].fill_value = "--"
        instrument = tmpTable["instrume"].filled()
        instrument = list(set(instrument))
        if "--" in instrument:
            instrument.remove("--")

        if len(instrument) == 2 and "SHOOT" in instrument and "XSHOOTER" in instrument:
            instrument = ["XSH"]

        if len(instrument) > 1:
            self.log.error(
                f"The directory contains data from a mix of instruments. Please only provide data from either SOXS or XSH"
            )
            raise AssertionError
        else:
            self.instrument = instrument[0]

        self._select_instrument(inst=self.instrument)

        if remainingFiles < 1:
            print(f"The instrument has been set to '{self.instrument}'")

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, instrument=self.instrument).get
        self.keywords = ["file"]
        for k in self.keyword_lookups:
            try:
                self.keywords.append(self.kw(k).lower())
            except Exception as e:
                self.log.warning(f"Keyword '{k}' not found in lookup table.")

        # TOP-LEVEL COLLECTION
        if recursive:
            allFrames = ImageFileCollection(filenames=fitsPaths, keywords=self.keywords)
            rawFrames = allFrames.summary
        else:
            # Split fitsNames into batches of 100
            batch_size = 1000
            batches = [fitsNames[i : i + batch_size] for i in range(0, len(fitsNames), batch_size)]

            from fundamentals import fmultiprocess

            results = fmultiprocess(
                log=self.log,
                function=_harvest_fits_headers,
                inputArray=batches,
                poolSize=False,
                timeout=300,
                pathToDirectory=pathToDirectory,
                keywords=self.keywords,
                filterKeys=filterKeys,
                instrument=self.instrument,
                kw=self.kw,
                turnOffMP=False,
                progressBar=True,
            )
            rawFrames = pd.concat(results)

        # FIX BOUNDARY GROUP FILES -- MOVE TO NEXT DAY SO THEY GET COMBINED WITH THE REST OF THEIR GROUP
        # E.G. UVB BIASES TAKEN ACROSS THE BOUNDARY BETWEEN 2 NIGHTS
        # FIRST FIND END OF NIGHT DATA - AND PUSH TO THE NEXT DAY
        mask = rawFrames["boundary"] > 0.96
        filteredDf = rawFrames.loc[mask].copy()
        filteredDf["night start mjd"] = filteredDf["night start mjd"] + 1
        mask = filteredDf["eso dpr type"].isin(["LAMP,DFLAT", "LAMP,QFLAT"])
        filteredDf.loc[mask, "eso dpr type"] = "LAMP,FLAT"
        # NOW FIND START OF NIGHT DATA
        mask = rawFrames["boundary"] < 0.04
        filteredDf2 = rawFrames.loc[mask].copy()
        mask = filteredDf2["eso dpr type"].isin(["LAMP,DFLAT", "LAMP,QFLAT"])
        filteredDf2.loc[mask, "eso dpr type"] = "LAMP,FLAT"

        # NOW FIND MATCHES BETWEEN 2 DATASETS
        theseKeys = [
            "eso seq arm",
            "eso dpr catg",
            "eso dpr tech",
            "eso dpr type",
            "eso pro catg",
            "eso pro tech",
            "eso pro type",
            "night start mjd",
        ]
        matched = pd.merge(filteredDf, filteredDf2, on=theseKeys)
        boundaryFiles = np.unique(matched["file_x"].values)
        mask = rawFrames["file"].isin(boundaryFiles)
        rawFrames.loc[mask, "night start mjd"] += 1

        rawFrames["night start date"] = Time(rawFrames["night start mjd"], format="mjd").to_value("iso", subfmt="date")

        self.log.debug("completed the ``_create_directory_table`` function")
        return rawFrames, fitsPaths, remainingFiles

    def _sync_sql_table_to_directory(self, directory, tableName, recursive=False):
        """*sync sql table content to files in a directory (add and delete from table as appropriate)*

        **Key Arguments:**

        - ``directory`` -- the directory of fits file to inspect.
        - ``tableName`` -- the sqlite table to sync.
        - ``recursive`` -- recursively dig into the directory to find FITS files? Default *False*.

        **Return:**

        - None

        **Usage:**

        ```python
        do._sync_sql_table_to_directory('/raw/directory/', 'raw_frames', recursive=False)
        ```
        """
        self.log.debug("starting the ``_sync_sql_table_to_directory`` method")

        import sqlite3 as sql
        import shutil

        # GENERATE A LIST OF FITS FILE PATHS IN RAW DIR
        from fundamentals.files import recursive_directory_listing

        fitsPaths = recursive_directory_listing(
            log=self.log, baseFolderPath=directory, whatToList="files"  # all | files | dirs
        )

        c = self.conn.cursor()

        sqlQuery = f"select filepath from {tableName};"
        c.execute(sqlQuery)

        # MAKE PATHS ABSOLUTE
        dbFiles = [r[0].replace("//", "/") for r in c.fetchall()]

        # DELETED FILES
        filesNotInDB = list(set(fitsPaths) - set(dbFiles))
        filesNotInFS = list(set(dbFiles) - set(fitsPaths))
        # MAKE PATHS RELATIVE TO rawDir
        filesNotInFS = [f.replace(self.rawDir + "/", "./raw/").replace("//", "/") for f in filesNotInFS]

        if len(filesNotInFS):
            filesNotInFS = ("','").join(filesNotInFS)
            sqlQuery = f"delete from {tableName} where filepath in ('{filesNotInFS}');"
            c.execute(sqlQuery)
            sqlQuery = f"delete from sof_map_{self.sessionId} where sof in (select sof from sof_map_{self.sessionId} where filepath in ('{filesNotInFS}'));"
            c.execute(sqlQuery)

        if len(filesNotInDB):

            for f in filesNotInDB:
                # GET THE EXTENSION (WITH DOT PREFIX)
                basename = os.path.basename(f)
                extension = os.path.splitext(basename)[1]
                if extension.lower() != ".fits":
                    pass
                elif self.rootDir in f:
                    exists = os.path.exists(os.path.abspath(self.rootDir) + "/" + basename)
                    if not exists:
                        os.symlink(os.path.realpath(f), os.path.abspath(self.rootDir) + "/" + basename)
                else:
                    shutil.move(self.rootDir + "/" + f, self.rootDir)
            self._sync_raw_frames(skipSqlSync=True)

        c.close()

        self.log.debug("completed the ``_sync_sql_table_to_directory`` method")
        return None

    def _populate_raw_frames_extra_columns(self, filteredFrames, verbose=False):
        """*populate extra columns for raw frames to later filter on*

        **Key Arguments:**

        - ``filteredFrames`` -- the dataframe from which to split frames into categorise.
        - ``verbose`` -- print results to stdout.

        **Return:**

        - ``rawFrames`` -- dataframe of raw frames only

        **Usage:**

        ```python
        rawFrames = self._populate_raw_frames_extra_columns(filteredFrames)
        ```
        """
        self.log.debug("starting the ``catagorise_frames`` method")

        from astropy.table import Table, unique
        import numpy as np
        import pandas as pd
        from tabulate import tabulate

        # SPLIT INTO RAW, REDUCED PIXELS, REDUCED TABLES
        rawFrameGroupKeywords = self.rawFrameGroupKeywords[:]
        filterKeywordsRaw = self.filterKeywords[:]

        filteredFrames["slit"] = "--"
        filteredFrames["slitmask"] = "--"
        filteredFrames["lamp"] = "--"
        filteredFrames["simulation"] = "--"

        # ADD SLIT FOR SPECTROSCOPIC DATA
        filteredFrames.loc[(filteredFrames["eso seq arm"] == "NIR"), "slit"] = filteredFrames.loc[
            (filteredFrames["eso seq arm"] == "NIR"), self.kw("SLIT_NIR").lower()
        ]
        filteredFrames.loc[(filteredFrames["eso seq arm"] == "VIS"), "slit"] = filteredFrames.loc[
            (filteredFrames["eso seq arm"] == "VIS"), self.kw("SLIT_VIS").lower()
        ]
        filteredFrames.loc[(filteredFrames["eso seq arm"] == "UVB"), "slit"] = filteredFrames.loc[
            (filteredFrames["eso seq arm"] == "UVB"), self.kw("SLIT_UVB").lower()
        ]

        # ADD SIMULATION FLAG FOR SPECTROSCOPIC DATA
        if self.instrument.lower() == "soxs":
            filteredFrames.loc[(filteredFrames["eso seq arm"] == "NIR"), "simulation"] = filteredFrames.loc[
                (filteredFrames["eso seq arm"] == "NIR"), self.kw("SWSIM_NIR").lower()
            ]
            filteredFrames.loc[(filteredFrames["eso seq arm"] == "VIS"), "simulation"] = filteredFrames.loc[
                (filteredFrames["eso seq arm"] == "VIS"), self.kw("SWSIM_VIS").lower()
            ]
            filteredFrames.loc[(filteredFrames["eso seq arm"] == "UVB"), "simulation"] = filteredFrames.loc[
                (filteredFrames["eso seq arm"] == "UVB"), self.kw("SWSIM_UVB").lower()
            ]
            filteredFrames.loc[(filteredFrames["simulation"] == "--"), "simulation"] = 0
            filteredFrames.loc[(filteredFrames["simulation"] == -99.99), "simulation"] = 0
            filteredFrames.loc[(filteredFrames["simulation"] == "T"), "simulation"] = 1
        else:
            filteredFrames["simulation"] = 0

        filteredFrames.loc[
            ((filteredFrames["slit"].str.contains("MULT")) & (filteredFrames["slitmask"] == "--")), "slitmask"
        ] = "MPH"
        filteredFrames.loc[
            ((filteredFrames["slit"].str.contains("PINHOLE")) & (filteredFrames["slitmask"] == "--")), "slitmask"
        ] = "PH"
        filteredFrames.loc[
            ((filteredFrames["slit"].str.contains("SLIT")) & (filteredFrames["slitmask"] == "--")), "slitmask"
        ] = "SLIT"

        lampLong = ["argo", "merc", "neon", "xeno", "qth", "deut", "thar"]
        lampEle = ["Ar", "Hg", "Ne", "Xe", "QTH", "D", "ThAr"]

        for i in [1, 2, 3, 4, 5, 6, 7]:
            lamp = self.kw(f"LAMP{i}").lower()

            if self.instrument.lower() == "soxs":
                for l, e in zip(lampLong, lampEle):
                    if l in lamp:
                        lamp = e
                filteredFrames.loc[
                    ((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] != "--")), "lamp"
                ] += lamp
                filteredFrames.loc[
                    ((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] == "--")), "lamp"
                ] = lamp
                filteredFrames.loc[
                    ((filteredFrames[self.kw("DPR_TYPE").lower()] == "DOME,FLAT") & (filteredFrames["lamp"] == "--")),
                    "lamp",
                ] = "DOME"

            else:
                filteredFrames.loc[((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99)), "lamp"] = (
                    filteredFrames.loc[
                        ((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99)), self.kw(f"LAMP{i}").lower()
                    ]
                )
        mask = []
        for i in self.proKeywords:
            rawFrameGroupKeywords.remove(i)
            filterKeywordsRaw.remove(i)
            if not len(mask):
                mask = filteredFrames[i] == "--"
            else:
                mask = np.logical_and(mask, (filteredFrames[i] == "--"))

        filteredFrames["lamp"] = filteredFrames["lamp"].str.replace("_lamp", "")
        filteredFrames["lamp"] = filteredFrames["lamp"].str.replace("_Lamp", "")

        rawFrames = filteredFrames.loc[mask]

        # MATCH OFF FRAMES TO ADD THE MISSING LAMPS
        mask = rawFrames["eso obs name"] == "Maintenance"
        rawFrames.loc[mask, "eso obs name"] = rawFrames.loc[mask, "eso obs name"] + rawFrames.loc[mask, "eso dpr type"]
        if self.instrument.lower() == "soxs":
            groupBy = "eso obs name"
        else:
            groupBy = "template"
        rawFrames.loc[(rawFrames["lamp"] == "--"), "lamp"] = np.nan
        rawFrames.loc[(rawFrames["eso seq arm"].str.lower() == "nir"), "lamp"] = rawFrames.loc[
            (rawFrames["eso seq arm"].str.lower() == "nir"), "lamp"
        ].fillna(
            rawFrames.loc[(rawFrames["eso seq arm"].str.lower() == "nir")].groupby(groupBy)["lamp"].transform("first")
        )
        rawFrames.loc[(rawFrames["lamp"].isnull()), "lamp"] = "--"
        rawFrames.loc[(rawFrames["absrot"] == -99.99), "absrot"] = 0.0

        rawFrames["exptime"] = rawFrames["exptime"].apply(lambda x: round(x, 2))
        rawGroups = rawFrames.groupby(filterKeywordsRaw).size().reset_index(name="counts")
        rawGroups.style.hide(axis="index")
        pd.options.mode.chained_assignment = None

        if verbose:
            print("\n# CONTENT FILE INDEX\n")
        if verbose and len(rawGroups.index):
            print("\n## ALL RAW FRAMES\n")
            print(
                tabulate(
                    rawFrames[rawFrameGroupKeywords],
                    headers="keys",
                    tablefmt="github",
                    showindex=False,
                    stralign="right",
                    floatfmt=".3f",
                )
            )

        self.log.debug("completed the ``catagorise_frames`` method")
        return rawFrames[rawFrameGroupKeywords].replace(["--"], None)

    def _populate_product_frames_db_table(self):
        """*scan the raw frame table to generate the listing of products that are expected to be created*

        **Key Arguments:**
            # -

        **Return:**

        - None

        **Usage:**

        ```python
        usage code
        ```
        """
        self.log.debug("starting the ``_populate_product_frames_db_table`` method")

        import pandas as pd
        import sqlite3 as sql
        import time
        from collections import defaultdict
        from itertools import product

        conn = self.conn
        rawFrames = pd.read_sql("SELECT * FROM raw_frames_valid", con=conn)

        c = conn.cursor()
        sqlQuery = "select count(*) from raw_frame_sets where complete = 1"
        c.execute(sqlQuery)
        completeCountStart = c.fetchall()[0][0]
        c.close()

        rawFrames = rawFrames.astype({"exptime": float, "ra": float, "dec": float})
        rawFrames.fillna({"exptime": -99.99, "ra": -99.99, "dec": -99.99}, inplace=True)
        rawFrames.fillna("--", inplace=True)
        filterKeywordsRaw = self.filterKeywords[:]

        for i in self.proKeywords:
            filterKeywordsRaw.remove(i)

        # rawFrames.replace("LAMP,DFLAT", "LAMP,FLAT", inplace=True)
        # rawFrames.replace("LAMP,QFLAT", "LAMP,FLAT", inplace=True)

        # HIDE OFF FRAMES FROM GROUPS
        mask = (
            (rawFrames["eso dpr tech"] == "IMAGE")
            & (rawFrames["eso seq arm"] == "NIR")
            & (rawFrames["eso dpr type"] != "DARK")
        )
        rawFramesNoOffFrames = rawFrames.loc[~mask]
        rawGroups = rawFramesNoOffFrames.groupby(filterKeywordsRaw)
        mjds = rawGroups.mean(numeric_only=True)["mjd-obs"].values

        rawGroups = rawGroups.size().reset_index(name="counts")
        rawGroups["mjd-obs"] = mjds

        # REMOVE GROUPED STARE - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = rawGroups["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"])
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD SCIENCE FRAMES AS ONE ENTRY PER EXPOSURE
        rawScienceFrames = pd.read_sql(
            "SELECT * FROM raw_frames_valid where `eso dpr tech` in ('ECHELLE,SLIT,STARE')", con=conn
        )

        # rawScienceFrames.fillna("--", inplace=True)
        rawScienceFrames = rawScienceFrames.groupby(filterKeywordsRaw + ["mjd-obs"])
        rawScienceFrames = rawScienceFrames.size().reset_index(name="counts")

        # MERGE DATAFRAMES
        if len(rawScienceFrames.index):
            rawGroups = pd.concat([rawGroups, rawScienceFrames], ignore_index=True)

        # REMOVE GROUPED SINGLE PINHOLE ARCS - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = rawGroups["eso dpr tech"].isin(["ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE"])
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD PINHOLE FRAMES AS ONE ENTRY PER EXPOSURE
        if self.instrument.upper() == "SOXS":
            rawPinholeFrames = pd.read_sql(
                "SELECT * FROM raw_frames_valid where `eso dpr tech` in ('ECHELLE,PINHOLE','ECHELLE,MULTI-PINHOLE') and ('eso seq arm' = 'NIR' or ('lamp' not in ('Xe', 'Ar', 'Hg', 'Ne', 'ArNeHgXe' )))",
                con=conn,
            )
        else:
            rawPinholeFrames = pd.read_sql(
                "SELECT * FROM raw_frames_valid where `eso dpr tech` in ('ECHELLE,PINHOLE','ECHELLE,MULTI-PINHOLE')",
                con=conn,
            )

        with pd.option_context("future.no_silent_downcasting", True):
            rawPinholeFrames = rawPinholeFrames.fillna({"exptime": -99.99, "ra": -99.99, "dec": -99.99}).infer_objects(
                copy=False
            )
        rawPinholeFrames.fillna("--", inplace=True)
        rawPinholeFrames = rawPinholeFrames.groupby(filterKeywordsRaw + ["mjd-obs"])
        rawPinholeFrames = rawPinholeFrames.size().reset_index(name="counts")

        # MERGE DATAFRAMES
        if len(rawPinholeFrames.index):
            rawGroups = pd.concat([rawGroups, rawPinholeFrames], ignore_index=True)

        rawGroups["recipe"] = None
        rawGroups["sof"] = None

        calibrationFrames = pd.read_sql(
            f"SELECT * FROM product_frames where `eso pro catg` not like '%_TAB_%' and (status_{self.sessionId} != 'fail' or status_{self.sessionId} is null)",
            con=conn,
        )
        calibrationFrames.fillna("--", inplace=True)

        calibrationTables = pd.read_sql(
            f"SELECT * FROM product_frames where `eso pro catg` like '%_TAB_%' and (status_{self.sessionId} != 'fail' or status_{self.sessionId} is null)",
            con=conn,
        )
        calibrationTables.fillna("--", inplace=True)

        # SEND TO DATABASE
        c = self.conn.cursor()
        sqlQuery = f"delete from sof_map_{self.sessionId};"
        c.execute(sqlQuery)
        c.close()

        if "complete" not in rawGroups.columns:
            rawGroups["complete"] = None
        if "recipe_order" not in rawGroups.columns:
            rawGroups["recipe_order"] = None

        # BULK ADD TAG COLUMN
        rawFrames["tag"] = rawFrames["eso dpr type"].replace(",", "_") + "_" + rawFrames["eso seq arm"]

        for o in self.reductionOrder:

            if "FLAT" in o:
                rawMask = rawFrames["eso dpr type"].str.contains("FLAT")
            else:
                rawMask = rawFrames["eso dpr type"].str.contains(o, case=False)
            rawFramesFiltered = rawFrames.loc[rawMask]

            keywords = ["slit", "eso seq arm", "binning", "rospeed", "eso dpr type", "lamp", "exptime", "slit"]
            rawFramesSubset = rawFramesFiltered[keywords].drop_duplicates().copy()
            # EXTRACT VALUES FOR EACH KEYWORD WITHOUT APPLYING UNIQUE FILTERING
            arms = rawFramesSubset["eso seq arm"]
            binnings = rawFramesSubset["binning"]
            rospeeds = rawFramesSubset["rospeed"]
            types = rawFramesSubset["eso dpr type"]
            lamps = rawFramesSubset["lamp"]
            exptimes = rawFramesSubset["exptime"]
            slits = rawFramesSubset["slit"]

            ## ALSO FILTER CALIBRATIONS

            # OPTIMISE: 95%
            if not "FLAsT" in o:
                # GET UNIQUE VALUES IN COLUMN
                slits = [False] * len(slits)

            if o.lower() not in self.typeMap:
                continue

            for row in self.typeMap[o.lower()]:
                slitmasks = False
                techs = False
                if row["slitmask"]:
                    slitmasks = row["slitmask"]
                if row["tech"]:
                    techs = [item.upper() for item in row["tech"]]
                recipe = row["recipe"]
                # Iterate over all permutations
                # for slit, arm, binning, rospeed, type_, lamp, exptime in product(
                #     slits, arms, binnings, rospeeds, types, lamps, exptimes
                # ):
                #     print(slit, arm, binning, rospeed, type_, lamp, exptime)
                #     pass

                # ITERATE OVER ALL PERMUTATIONS
                for slit, arm, binning, rospeed in zip(slits, arms, binnings, rospeeds):

                    rawGroups = self.generate_sof_names_and_maps(
                        reductionOrder=o,
                        rawFrames=rawFramesFiltered,
                        calibrationFrames=calibrationFrames,
                        calibrationTables=calibrationTables,
                        rawGroups=rawGroups,
                        recipe=recipe,
                        slit=slit,
                        slitmasks=slitmasks,
                        techs=techs,
                        arm=arm,
                        binning=binning,
                        rospeed=rospeed,
                    )
                    rawGroups = self.populate_products_table(reductionOrder=o, rawGroups=rawGroups)

        mask = rawGroups["complete"] == 1
        completeCountEnd = len(rawGroups.loc[mask].index)

        if completeCountStart == completeCountEnd:
            return completeCountStart

        # SEND TO DATABASE
        c = self.conn.cursor()
        sqlQuery = f"delete from raw_frame_sets;"
        c.execute(sqlQuery)
        c.close()

        keepTrying = 0
        while keepTrying < 7:
            try:
                rawGroups.replace(["--"], None).to_sql("raw_frame_sets", con=self.conn, index=False, if_exists="append")
                keepTrying = 10
            except Exception as e:
                if keepTrying > 5:
                    raise Exception(e)
                time.sleep(1)
                keepTrying += 1

        self.log.debug("completed the ``_populate_product_frames_db_table`` method")
        return None

    @profile
    def _generate_sof_name_and_map(self, series, rawFrames, calibrationFrames, calibrationTables, recipe):
        """*add a recipe name and SOF filename to all rows in the raw_frame_sets DB table*

        **Key Arguments:**

        - ``series`` -- the dataframe row/series to apply work on
        """

        import pandas as pd
        import astropy
        import numpy as np

        incomplete = False

        sofName = []
        matchDict = {}
        sofName.append(self.sofName)
        matchDict["eso seq arm"] = series["eso seq arm"].upper()
        filteredFrames = rawFrames

        seriesRecipe = recipe

        # print(series["eso dpr type"], series["eso seq arm"], series["eso dpr tech"], series["slitmask"])

        # GENERATE SOF FILENAME AND MATCH DICTIONARY TO FILTER ON
        if series["rospeed"] != "--":
            matchDict["rospeed"] = series["rospeed"]
        if series["binning"] != "--":
            matchDict["binning"] = series["binning"]

        if series["eso dpr type"].lower() in self.typeMap:
            for i in self.typeMap[series["eso dpr type"].lower()]:
                if i["recipe"] == seriesRecipe:
                    sofName.append(i["recipe"].replace("_centres", "_locations"))

            if "DORDER" in series["eso dpr type"].upper():
                sofName.append("dlamp")
            if "QORDER" in series["eso dpr type"].upper():
                sofName.append("qlamp")
        # if True and ("PINHOLE" in series["eso dpr tech"].upper() or (series["instrume"] == "SOXS" and "FLAT" in series["eso dpr type"].upper())):
        if True and (series["instrume"] == "SOXS" and "FLAT" in series["eso dpr type"].upper()):
            if series["lamp"] != "--":
                matchDict["lamp"] = series["lamp"]
                sofName.append(series["lamp"])

        if series["instrume"] == "SOXS":
            sofName.append(series["slit"])
        if series["exptime"] and (
            series["eso seq arm"].lower() == "nir"
            or (
                series["eso seq arm"].lower() == "vis"
                and ("FLAT" in series["eso dpr type"].upper() or "DARK" in series["eso dpr type"].upper())
            )
        ):
            matchDict["exptime"] = float(series["exptime"])
            sofName.append(str(series["exptime"]) + "S")
        elif (
            series["exptime"]
            and "BIAS" not in series["eso dpr type"].upper()
            and not ("FLAT" in series["eso dpr type"].upper() and series["eso seq arm"].upper() == "UVB")
        ):
            sofName.append(str(series["exptime"]) + "S")

        if series["eso obs name"] != "--":
            matchDict["eso obs name"] = series["eso obs name"]

        if "std" in series["eso dpr type"].lower() or "object" in series["eso dpr type"].lower():
            matchDict["slit"] = series["slit"]

        sofName.append(str(series["instrume"]))

        for k, v in matchDict.items():

            if "type" in k.lower() and "lamp" in v.lower() and "flat" in v.lower():
                mask = filteredFrames[k].isin(["LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT", "FLAT,LAMP", "DOME,FLAT"])
            elif k not in ["rospeed", "binning", "eso seq arm"]:
                print(k)
                mask = filteredFrames[k].isin([v])
            else:
                continue
            filteredFrames = filteredFrames.loc[mask]

        # INITIAL CALIBRATIONS FILTERING
        if series["eso seq arm"].upper() in ["UVB", "VIS"]:
            for k, v in matchDict.items():
                if k in ["exptime", "lamp", "eso obs name", "eso obs id"]:
                    continue
                if "type" in k.lower():
                    mask = (calibrationFrames["eso pro catg"].str.contains("MASTER_")) | (
                        calibrationFrames["eso pro catg"].str.contains("DISP_IMAGE")
                    )
                elif "slit" in k.lower():
                    mask = ~(
                        ~calibrationFrames[k].isin([v])
                        & (calibrationFrames["eso pro catg"].str.contains("MASTER_MFLAT"))
                    )
                elif "rospeed" in k.lower() or "binning" in k.lower():
                    mask = calibrationFrames[k].isin([v]) | (
                        calibrationFrames["eso pro catg"].str.contains("DISP_IMAGE")
                    )
                else:
                    mask = calibrationFrames[k].isin([v])
                calibrationFrames = calibrationFrames.loc[mask]

        elif series["eso seq arm"].upper() in ["NIR"]:
            for k, v in matchDict.items():
                if k in ["binning", "rospeed", "exptime", "lamp", "eso obs name", "eso obs id"]:
                    continue
                if "type" in k.lower():
                    mask = calibrationFrames["eso pro catg"].str.contains("MASTER_") | (
                        calibrationFrames["eso pro catg"].str.contains("DISP_IMAGE")
                    )
                elif "slit" in k.lower():
                    mask = ~(
                        ~calibrationFrames[k].isin([v])
                        & (calibrationFrames["eso pro catg"].str.contains("MASTER_MFLAT"))
                    )
                else:
                    mask = calibrationFrames[k].isin([v])
                calibrationFrames = calibrationFrames.loc[mask]

        # EXTRA CALIBRATION TABLES
        for k, v in matchDict.items():
            if k in ["rospeed", "exptime", "lamp", "eso obs name", "slit", "eso obs id"]:
                continue
            if k in ["binning"] and seriesRecipe in ["mflat"]:
                continue
            if k in ["binning"]:
                mask = calibrationTables[k].isin([v]) | calibrationTables["eso pro catg"].str.contains("DISP_TAB_")
            elif "type" in k.lower() and series["eso seq arm"] in ["UVB", "VIS", "NIR"]:
                mask = calibrationTables["eso pro catg"].str.contains("_TAB_")
            else:
                try:
                    mask = calibrationTables[k].isin([v])
                except Exception as e:
                    print(series)
                    print(k, v)
                    print(e)
                    sys.exit(0)

            calibrationTables = calibrationTables.loc[mask]

        # NIGHT START
        # YYYY.MM.DDThh.mm.xxx
        offFrameCount = 0
        if series["night start mjd"]:

            if "PINHOLE" in series["eso dpr tech"].upper():
                if "NIR" not in series["eso seq arm"].upper():
                    mask = filteredFrames["mjd-obs"] == series["mjd-obs"]
                    filteredFrames = filteredFrames.loc[mask]
                else:
                    # FURTHER FILTERING OF NIR ON/OFF FRAMES
                    filteredFrames["obs-delta"] = filteredFrames["mjd-obs"] - series["mjd-obs"]
                    filteredFrames["obs-delta"] = filteredFrames["obs-delta"].abs()
                    filteredFrames.sort_values(["obs-delta"], inplace=True)
                    mask = filteredFrames["eso dpr tech"].isin(["IMAGE"])
                    offFrame = filteredFrames.loc[mask].head(1)
                    onFrame = filteredFrames.loc[~mask].head(1)
                    filteredFrames = pd.concat([onFrame, offFrame], ignore_index=True)
                    offFrameCount = len(offFrame.index)

            if "FLAT" in series["eso dpr tech"].upper() and "NIR" not in series["eso seq arm"].upper():
                mask = filteredFrames["eso dpr tech"].isin(["IMAGE"])
                offFrame = filteredFrames.loc[mask]

            if series["eso dpr tech"] in ["ECHELLE,SLIT,STARE"]:
                mask = filteredFrames["mjd-obs"] == series["mjd-obs"]
                filteredFrames = filteredFrames.loc[mask]
            else:
                mask = filteredFrames["night start mjd"] == int(series["night start mjd"])
                filteredFrames = filteredFrames.loc[mask]

            if series["eso dpr tech"] in ["ECHELLE,SLIT,NODDING"]:
                mask = filteredFrames["exptime"] == series["exptime"]
                filteredFrames = filteredFrames.loc[mask]

            try:
                frameMjd = filteredFrames["mjd-obs"].values[0]
            except:
                print(series)

            calibrationFrames = calibrationFrames.copy()
            calibrationTables = calibrationTables.copy()
            calibrationFrames["obs-delta"] = calibrationFrames["mjd-obs"] - frameMjd
            calibrationTables["obs-delta"] = calibrationTables["mjd-obs"] - frameMjd
            # dispImages["obs-delta"] = dispImages['mjd-obs'] - frameMjd
            calibrationFrames["obs-delta"] = calibrationFrames["obs-delta"].abs()
            calibrationTables["obs-delta"] = calibrationTables["obs-delta"].abs()
            # dispImages["obs-delta"] = dispImages["obs-delta"].abs()
            if isinstance(filteredFrames, astropy.table.row.Row):
                from astropy.table import Table

                filteredFrames = Table(filteredFrames)

            if seriesRecipe not in ["mbias", "mdark"]:
                mask = filteredFrames["eso dpr tech"].isin(["IMAGE"])
            else:
                mask = filteredFrames["eso dpr tech"].isin(["NONSENSE"])

            firstDate = filteredFrames.loc[~mask]["date-obs"].values[0].split(".")[0]
            firstDate = firstDate.replace("-", "").replace(":", "")
            sofName.insert(0, firstDate)

        # NEED SOME FINAL FILTERING ON UVB FLATS
        if (
            "lamp" in series["eso dpr type"].lower()
            and "flat" in series["eso dpr type"].lower()
            and "uvb" in series["eso seq arm"].lower()
        ):
            mask = filteredFrames["eso dpr type"].isin(["LAMP,DFLAT"])
            dFrames = filteredFrames.loc[mask]
            dexptime = dFrames["exptime"].max()
            mask = filteredFrames["eso dpr type"].isin(["LAMP,QFLAT"])
            qFrames = filteredFrames.loc[mask]
            qexptime = qFrames["exptime"].max()
            mask = ((filteredFrames["eso dpr type"] == "LAMP,DFLAT") & (filteredFrames["exptime"] == dexptime)) | (
                (filteredFrames["eso dpr type"] == "LAMP,QFLAT") & (filteredFrames["exptime"] == qexptime)
            )
            filteredFrames = filteredFrames.loc[mask]

        # LAMP ON AND OFF TAGS
        if series["eso seq arm"].lower() == "nir":
            if seriesRecipe in ["disp_sol", "order_centres", "mflat", "spat_sol"]:
                mask = filteredFrames["eso dpr tech"].isin(["IMAGE"])
                filteredFrames.loc[mask, "tag"] += ",OFF"
                filteredFrames.loc[~mask, "tag"] += ",ON"

        files = filteredFrames["file"].values
        filepaths = filteredFrames["filepath"].values
        tags = filteredFrames["tag"].values

        if seriesRecipe in ("stare", "nod", "offset"):
            objectt = filteredFrames["object"].values[0].replace(" ", "_")[:8]
            sofName.append(objectt)

        # COMBINE AND SHORTEN SOF NAME
        sofName = ("_").join(sofName).replace("-", "").replace(",", "_").upper() + ".sof"
        sofName = sofName.replace("XSHOOTER", "XSH")
        sofName = sofName.replace("FAST", "F")
        sofName = sofName.replace("SLOW", "S")
        sofName = sofName.replace("TELLURIC", "TELL")
        sofName = sofName.replace("ORDER_LOCATIONS", "OLOC")
        sofName = sofName.replace("DISP_SOL", "DSOL")
        sofName = sofName.replace("SPAT_SOL", "SSOL")

        series["sof"] = sofName
        series["recipe"] = seriesRecipe
        series["recipe_order"] = self.recipeOrder.index(seriesRecipe) + 1

        # CALIBRATIONS NEEDED?
        # BIAS FRAMES
        if series["recipe"] in ["disp_sol", "order_centres", "mflat", "spat_sol", "stare", "nod", "offset"]:

            mask = calibrationFrames["eso pro catg"].isin([f"MASTER_BIAS_{series['eso seq arm'].upper()}"])
            df = calibrationFrames.loc[mask]

            if len(df.index):
                df = df.sort_values(by=["obs-delta"])
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])
            elif series["eso seq arm"].upper() != "NIR":
                incomplete = True

        # DISP SOLS
        if series["recipe"] in ["order_centres", "spat_sol", "stare", "nod", "offset"]:
            mask = calibrationTables["eso pro catg"].str.contains("DISP_TAB")
            df = calibrationTables.loc[mask]

            if len(df.index):
                df = df.sort_values(by=["obs-delta"])
                if series["recipe"] in ["stare", "nod", "offset"]:
                    mask = df["recipe"] == "spat_sol"
                    if len(df.loc[mask, "file"].index):
                        files = np.append(files, df.loc[mask, "file"].values[0])
                        tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                        filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                    else:
                        incomplete = True
                else:
                    mask = df["recipe"] == "disp_sol"
                    if len(df.loc[mask, "file"].index):
                        files = np.append(files, df.loc[mask, "file"].values[0])
                        tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                        filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                    else:
                        incomplete = True
            else:
                incomplete = True

        # DISP SOLS IMAGE
        if series["recipe"] in ["stare", "nod", "offset"]:
            mask = calibrationFrames["eso pro catg"].str.contains("DISP_IMAGE")
            df = calibrationFrames.loc[mask]
            if len(df.index):
                df = df.sort_values(by=["obs-delta"])
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])
            if "object" in series["eso dpr type"].lower():
                mask = calibrationTables["eso pro catg"].str.contains("RESP_TAB") & calibrationTables["slit"].isin(
                    [series["slit"]]
                )
                df = calibrationTables.loc[mask]
                if len(df.index):
                    df = df.sort_values(by=["obs-delta"])
                    files = np.append(files, df["file"].values[0])
                    tags = np.append(tags, df["eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df["filepath"].values[0])

        # ORDER TAB
        if series["recipe"] in ["mflat", "spat_sol", "stare", "nod", "offset"]:
            mask = calibrationTables["eso pro catg"].str.contains("ORDER_TAB")
            df = calibrationTables.loc[mask]
            if series["recipe"] in ["mflat"]:
                mask = df["recipe"] == "order_centres"
            elif series["recipe"] in ["spat_sol"]:
                mask = df["recipe"] == "mflat"
            elif series["recipe"] in ["stare", "nod", "offset"]:
                mask = (df["recipe"] == "mflat") & (df["slit"] == series["slit"])
                if (
                    len(df.loc[mask].index) == 0
                    and series["eso seq arm"].lower() == "nir"
                    and self.instrument.upper() == "SOXS"
                    and series["slit"] == "SLIT5.0"
                ):
                    mask = (df["recipe"] == "mflat") & (df["slit"] == "SLIT1.5")
            df = df.loc[mask]

            if series["recipe"] in ["spat_sol"]:
                # REMOVE BLOCKING FILTERS
                mask = df["slit"].str.contains("JH")
                df = df.loc[~mask]

            if len(df.index):
                df = df.sort_values(by=["obs-delta"])
                if series["eso seq arm"].upper() in ["UVB"] and series["recipe"] == "mflat":
                    lamps = ["QLAMP", "DLAMP"]
                    for lamp in lamps:
                        mask = df["file"].str.contains(lamp)
                        if len(df.loc[mask, "file"].index):
                            files = np.append(files, df.loc[mask, "file"].values[0])
                            tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                            filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                elif series["recipe"] == "mflat":
                    files = np.append(files, df["file"].values[0])
                    tags = np.append(tags, df["eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df["filepath"].values[0])
                else:
                    mask = df["filepath"].str.contains("soxs-mflat/")
                    files = np.append(files, df.loc[mask, "file"].values[0])
                    tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
            elif series["recipe"] in ["spat_sol", "mflat", "stare", "nod", "offset"]:
                incomplete = True

        # FLAT FRAMES
        if series["recipe"] in ["spat_sol", "stare", "nod", "offset"]:
            mask = calibrationFrames["eso pro catg"].str.contains("MASTER_FLAT")
            df = calibrationFrames.loc[mask]

            if series["recipe"] in ["spat_sol"]:
                # REMOVE BLOCKING FILTERS
                mask = df["slit"].str.contains("JH")
                df = df.loc[~mask]

            if series["recipe"] in ["stare", "nod", "offset"]:

                if len(filteredFrames["slit"].values):
                    if self.PAE and self.instrument.upper() == "SOXS":
                        pass
                    else:
                        mask = df["slit"] == filteredFrames["slit"].values[0]
                        if (len(df.loc[mask].index) == 0) and (
                            series["eso seq arm"].lower() == "nir"
                            and self.instrument.upper() == "SOXS"
                            and series["slit"] == "SLIT5.0"
                        ):
                            mask = df["slit"] == "SLIT1.5"
                        df = df.loc[mask]

            # FAVOUR DOME FLATS IF AVAILABLE
            if False:
                mask = df["lamp"] == "DOME"
                if len(df.loc[mask].index):
                    df = df.loc[mask]

            if len(df.index):
                df = df.sort_values(by=["obs-delta"])
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # DARK FRAMES
        if series["eso seq arm"].lower() == "nir" and series["recipe"] in [
            "stare",
            "disp_sol",
            "order_centres",
            "mflat",
            "spat_sol",
        ]:
            moveOn = False
            if series["recipe"] == "mflat":
                if "FLAT" in series["eso dpr type"].upper() and "NIR" in series["eso seq arm"].upper():
                    mask = filteredFrames["eso dpr tech"].isin(["IMAGE"])
                    offFrame = filteredFrames.loc[mask]
                    offFrameCount = len(offFrame.index)
            if series["recipe"] in ["disp_sol", "order_centres", "mflat", "spat_sol"] and offFrameCount > 0:
                moveOn = True
            elif series["recipe"] in ["stare", "disp_sol", "order_centres", "mflat", "spat_sol"] and offFrameCount == 0:
                mask = (calibrationFrames["exptime"] == filteredFrames["exptime"].max()) & (
                    calibrationFrames["eso pro catg"].str.contains("MASTER_DARK")
                )
                df = calibrationFrames.loc[mask]
                if len(df.index) == 0:
                    mask = calibrationFrames["eso pro catg"].str.contains("MASTER_DARK")

            else:
                mask = calibrationFrames["eso pro catg"].str.contains("MASTER_DARK")

            if not moveOn:
                df = calibrationFrames.loc[mask]
                if len(df.index):
                    df = df.sort_values(by=["obs-delta"])
                    files = np.append(files, df["file"].values[0])
                    tags = np.append(tags, df["eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df["filepath"].values[0])

        # SLIT ARC-LAMP FRAMES
        if series["recipe"] in ["spat_sol"]:
            if True:
                mask = (
                    (rawFrames["binning"] == series["binning"])
                    & (rawFrames["night start mjd"] == series["night start mjd"])
                    & (rawFrames["eso seq arm"] == series["eso seq arm"])
                    & (rawFrames["rospeed"] == series["rospeed"])
                    & (rawFrames["eso dpr type"].isin(["LAMP,WAVE", "WAVE,LAMP"]))
                    & (rawFrames["eso dpr tech"].isin(["ECHELLE,SLIT"]))
                    & (rawFrames["slit"].isin(["SLIT1.0"]))
                )
            else:
                # THIS IS FOR TESTING XSHOOTER RESOLUTION MEASUREMENTS ONLY
                mask = (
                    (rawFrames["eso seq arm"] == series["eso seq arm"])
                    & (rawFrames["rospeed"] == series["rospeed"])
                    & (rawFrames["eso dpr type"].isin(["LAMP,WAVE", "WAVE,LAMP"]))
                    & (rawFrames["eso dpr tech"].isin(["ECHELLE,SLIT"]))
                )
            df = rawFrames.loc[mask].copy()
            if len(df.index):
                df["obs-delta"] = df["mjd-obs"] - series["mjd-obs"]
                df = df.sort_values(by=["obs-delta"])
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, "SLIT_ARC")
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {
            "file": files,
            "filepath": filepaths,
            "tag": tags,
            "sof": [series["sof"]] * len(tags),
        }

        if incomplete:
            series["complete"] = 0
        else:
            series["complete"] = 1

        if incomplete:
            myDict["complete"] = [0] * len(tags)
        else:
            myDict["complete"] = [1] * len(tags)

        self.sofMaps.append(myDict)

        return series

    def _populate_products_table_by_row(self, series, reductionOrder):
        """*determine what the products should be for a given recipe and SOF file and populate the products table*

        **Key Arguments:**

        - ``recipeName`` -- the name of the recipe.
        - ``sofName`` -- the name of the sof file.

        **Return:**

        - None

        """
        self.log.debug("starting the ``_populate_products_table_by_row`` method")

        import pandas as pd
        import time

        products = series.to_dict()

        if products["recipe"] and products["recipe"].lower() in self.productMap:

            for k, v in self.productMap[products["recipe"].lower()].items():
                if k == "ALL" or k == reductionOrder:
                    for i in v:
                        products["eso pro type"] = i[0]
                        products["eso pro tech"] = i[1]
                        products["eso pro catg"] = i[2] + f"_{products['eso seq arm'].upper()}"
                        products["file"] = products["sof"].replace(".sof", ".fits")
                        if i[4] and i[5]:
                            products["file"] = products["file"].replace(i[4], i[5]) + ".fits"
                            products["file"] = products["file"].replace(".fits.fits", ".fits")
                        products["filepath"] = (
                            f"./reduced/{products['night start date']}/" + i[6] + "/" + products["file"]
                        )
                        myDict = {k: [v] for k, v in products.items()}
                        self.products.append(myDict)

        self.log.debug("completed the ``_populate_products_table_by_row`` method")
        return series

    def _move_misc_files(self):
        """*move extra/miscellaneous files to a misc directory*"""
        self.log.debug("starting the ``_move_misc_files`` method")

        import shutil

        if not os.path.exists(self.miscDir):
            os.makedirs(self.miscDir)

        # GENERATE A LIST OF FILE PATHS
        pathToDirectory = self.rootDir
        allowlistExtensions = [".db", ".yaml", ".log", ".sh"]
        for d in os.listdir(pathToDirectory):
            filepath = os.path.join(pathToDirectory, d)
            if os.path.splitext(filepath)[1] in allowlistExtensions:
                continue
            if os.path.isfile(filepath) and os.path.splitext(filepath)[1] != ".db" and "readme." not in d.lower():
                shutil.move(filepath, self.miscDir + "/" + d)

        self.log.debug("completed the ``_move_misc_files`` method")
        return None

    def _write_sof_files(self):
        """*Write out all possible SOF files from the sof_map database table*

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug("starting the ``_write_sof_files`` method")

        import pandas as pd
        from tabulate import tabulate

        conn, reset = self._get_or_create_db_connection()

        # RECURSIVELY CREATE MISSING DIRECTORIES
        self.sofDir = self.sessionPath + "/sof"
        if not os.path.exists(self.sofDir):
            os.makedirs(self.sofDir)

        df = pd.read_sql_query(f"select * from sof_map_{self.sessionId} where complete = 1;", conn)

        # GROUP RESULTS
        for name, group in df.groupby("sof"):
            sofPath = self.sofDir + "/" + name
            if os.path.exists(sofPath):
                continue
            myFile = open(sofPath, "w")
            content = tabulate(group[["filepath", "tag"]], tablefmt="plain", showindex=False)
            myFile.write(content)
            myFile.close()

        self.log.debug("completed the ``_write_sof_files`` method")
        return None

    def _write_reduction_shell_scripts(self):
        """*write the reduction shell scripts*

        _reduce_all.sh HAS BEEN REPLACED WITH THE `SOXSPIPE REDUCE` COMMAND

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug("starting the ``_write_reduction_shell_scripts`` method")

        import pandas as pd
        import sqlite3 as sql

        conn, reset = self._get_or_create_db_connection()

        rawGroups = pd.read_sql(
            "SELECT * FROM raw_frame_sets where recipe_order is not null order by recipe_order", con=conn
        )

        rawGroups["command"] = "soxspipe " + rawGroups["recipe"] + " sof/" + rawGroups["sof"]

        # WRITE FULL REDUCTION SCRIPT
        myFile = open(self.sessionPath + "/_reduce_all.sh", "w")
        myFile.write(("\n").join(pd.unique(rawGroups["command"])))
        myFile.close()

        self.log.debug("completed the ``_write_reduction_shell_scripts`` method")
        return None

    def session_create(self, sessionId=False):
        """*create a data-reduction session with accompanying settings file and required directories*

        **Key Arguments:**

        - ``sessionId`` -- optionally provide a sessionId (A-Z, a-z 0-9 and/or _- allowed, 16 character limit)

        **Return:**

        - ``sessionId`` -- the unique ID of the data-reduction session

        **Usage:**

        ```python
        do = data_organiser(
            log=log,
            rootDir="/path/to/workspace/root/"
        )
        sessionId = do.session_create(sessionId="my_supernova")
        ```
        """
        self.log.debug("starting the ``session_create`` method")

        import re
        import shutil
        import sqlite3 as sql

        rootDbExists = os.path.exists(self.rootDbPath)
        if rootDbExists:
            # CREATE THE DATABASE CONNECTION
            self.conn, reset = self._get_or_create_db_connection()
            c = self.conn.cursor()
            sqlQuery = "select distinct instrume from raw_frames"
            c.execute(sqlQuery)
            inst = c.fetchall()[0][0].lower()
            c.close()

        # TEST SESSION DIRECTORY EXISTS
        exists = os.path.exists(self.sessionsDir)
        if not exists:
            print("Please prepare your workspace using the `soxspipe prep` command before creating a new session.")
            sys.exit(0)

        if sessionId:
            if len(sessionId) > 16:
                print("Session ID must be 16 characters long or shorter, consisting of A-Z, a-z, 0-9 and/or _-")
            matchObjectList = re.findall(r"[^0-9a-zA-Z\-\_]+", sessionId)
            if matchObjectList:
                print("Session ID must be 16 characters long or shorter, consisting of A-Z, a-z, 0-9 and/or _-")
        else:
            # CREATE SESSION ID FROM TIME STAMP
            from datetime import datetime, date, time

            now = datetime.now()
            sessionId = now.strftime("%Y%m%dt%H%M%S")

        self.sessionId = sessionId

        # MAKE THE SESSION DIRECTORY
        self.sessionPath = self.sessionsDir + "/" + sessionId
        if not os.path.exists(self.sessionPath):
            os.makedirs(self.sessionPath)

        # SETUP SESSION SETTINGS AND LOGGING
        testPath = self.sessionPath + "/soxspipe.yaml"
        exists = os.path.exists(testPath)

        if not exists:
            if "shoo" in inst:
                inst = "xsh"
            su = tools(
                arguments={"<workspaceDirectory>": self.sessionPath, "init": True, "settingsFile": None},
                docString=False,
                logLevel="WARNING",
                options_first=False,
                projectName="soxspipe",
                defaultSettingsFile=f"{inst}_default_settings.yaml",
                createLogger=False,
            )
            arguments, settings, replacedLog, dbConn = su.setup()

        # MAKE ASSET PLACEHOLDERS
        if self.vlt:
            dest = self.sessionPath + "/reduced"
            try:
                os.symlink(self.vltReduced, dest)
            except:
                pass

        folders = ["sof", "qc", "reduced"]
        for f in folders:
            if not os.path.exists(self.sessionPath + f"/{f}"):
                os.makedirs(self.sessionPath + f"/{f}")

        # ADD A NEW STATUS COLUMN IN product_frames FOR THIS SESSION
        import sqlite3 as sql

        conn, reset = self._get_or_create_db_connection()
        c = conn.cursor()
        sqlQuery = f"ALTER TABLE product_frames ADD status_{sessionId} TEXT;"
        try:
            c.execute(sqlQuery)
        except:
            pass

        # DUPLICATE TEH SOF_MAP TABLE
        sqlQuery = "SELECT sql FROM sqlite_master WHERE type='table' AND name='z_sof_map'"
        c.execute(sqlQuery)
        sqlQuery = c.fetchall()[0][0]
        sqlQuery = sqlQuery.replace("z_sof_map", f"sof_map_{sessionId}")
        try:
            c.execute(sqlQuery)
        except:
            pass

        sqlQueries = [f"DROP VIEW IF EXISTS sof_map;", f"CREATE VIEW sof_map as select * from sof_map_{sessionId};"]
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)

        c.close()

        self._write_sof_files()

        # _reduce_all.sh HAS BEEN REPLACED WITH THE `SOXSPIPE REDUCE` COMMAND
        # self._write_reduction_shell_scripts()

        self._symlink_session_assets_to_workspace_root()

        # WRITE THE SESSION ID FILE
        import codecs

        with codecs.open(self.sessionIdFile, encoding="utf-8", mode="w") as writeFile:
            writeFile.write(sessionId)

        message = f"A new data-reduction session has been created with sessionId '{sessionId}'"
        try:
            self.log.print(message)
        except:
            print(message)
        self.log.debug("completed the ``session_create`` method")

        return sessionId

    def session_list(self, silent=False):
        """*list the sessions available to the user*

        **Key Arguments:**

        - ``silent`` -- don't print listings if True

        **Return:**

        - ``currentSession`` -- the single ID of the currently used session
        - ``allSessions`` -- the IDs of the other sessions

        **Usage:**

        ```python
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            rootDir="."
        )
        currentSession, allSessions = do.session_list()
        ```
        """
        self.log.debug("starting the ``session_list`` method")

        import codecs

        # IF SESSION ID FILE DOES NOT EXIST, REPORT
        self.sessionIdFile = self.sessionsDir + "/.sessionid"
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        else:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                currentSession = readFile.read()

        # LIST ALL SESSIONS
        allSessions = [d for d in os.listdir(self.sessionsDir) if os.path.isdir(os.path.join(self.sessionsDir, d))]
        allSessions.sort()

        if not silent:
            for s in allSessions:
                if s == currentSession.strip():
                    print(f"\033[0;32m*{s}*\u001b[38;5;15m")
                else:
                    print(s)

        self.log.debug("completed the ``session_list`` method")
        return currentSession, allSessions

    def session_switch(self, sessionId):
        """*switch to an existing workspace data-reduction session*

        **Key Arguments:**

        - ``sessionId`` -- the sessionId to switch to

        **Usage:**

        ```python
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            rootDir="."
        )
        do.session_switch(mySessionId)
        ```
        """
        self.log.debug("starting the ``session_switch`` method")
        import codecs

        currentSession, allSessions = self.session_list(silent=True)

        if sessionId == currentSession:
            print(f"Session '{sessionId}' is already in use.")
            return None
        elif sessionId in allSessions:
            # WRITE THE SESSION ID FILE
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="w") as writeFile:
                writeFile.write(sessionId)
        else:
            print(f"There is no session with the ID '{sessionId}'. List existing sessions with `soxspipe session ls`.")
            return None

        self.sessionPath = self.sessionsDir + "/" + sessionId
        self._symlink_session_assets_to_workspace_root()
        print(f"Session successfully switched to '{sessionId}'.")

        self.log.debug("completed the ``session_switch`` method")
        return None

    def _symlink_session_assets_to_workspace_root(self):
        """*symlink session QC, product, SOF directories, database and scripts to workspace root*

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug("starting the ``_symlink_session_assets_to_workspace_root`` method")

        import shutil
        import os

        # SYMLINK FILES AND FOLDERS
        toLink = ["reduced", "qc", "soxspipe.yaml", "sof", "soxspipe.log"]
        for l in toLink:
            dest = self.rootDir + f"/{l}"
            src = self.sessionPath + f"/{l}"
            try:
                os.symlink(src, dest)
            except:
                os.unlink(dest)
                os.symlink(src, dest)

        # REDUCTION SCRIPTS
        for d in os.listdir(self.sessionPath):
            filepath = os.path.join(self.sessionPath, d)
            if os.path.isfile(filepath) and os.path.splitext(filepath)[1] == ".sh":
                dest = self.rootDir + f"/{d}"
                src = filepath
                try:
                    os.symlink(src, dest)
                except:
                    os.unlink(dest)
                    os.symlink(src, dest)

        self.log.debug("completed the ``_symlink_session_assets_to_workspace_root`` method")
        return None

    def session_refresh(self, silent=False, failure=True):
        """*refresh a session's SOF files (needed if a recipe fails)*

        **Usage:**

        ```python
        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            rootDir="."
        )
        do.session_refresh()
        ```
        """
        self.log.debug("starting the ``session_refresh`` method")

        import sys
        import os
        import pandas as pd

        if failure is True:
            self.log.print("\nRefeshing SOF files due to recipe failure\n")
        elif failure is False:
            self.log.print("\nRefreshing SOF files, as a previously failed recipe is now passing.\n")
        import codecs

        # IF SESSION ID FILE DOES NOT EXIST, REPORT
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        else:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                sessionId = readFile.read()
        self.sessionPath = self.sessionsDir + "/" + sessionId
        self.sessionPath = self.sessionsDir + "/" + sessionId
        self.sessionId = sessionId

        self.conn, reset = self._get_or_create_db_connection()

        # SELECT INSTR
        self._select_instrument()

        self.build_sof_files()

        if failure in [True, False]:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print("SOF file refresh complete")

        self.log.debug("completed the ``session_refresh`` method")
        return reset

    def close(self):
        """*close the database connection*

        **Usage:**

        ```python
        do.close()
        ```
        """
        self.log.debug("starting the ``session_refresh`` method")

        try:
            self.conn.close()
        except:
            pass

        self.log.debug("completed the ``session_refresh`` method")
        return

    def use_vlt_environment_folders(self):
        """*use vlt environment folders*

        **Key Arguments:**
            # -

        **Return:**
            - None

        **Usage:**

        ```python
        usage code
        ```

        ---

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        ```
        """
        self.log.debug("starting the ``use_vlt_environment_folders`` method")

        import yaml

        # COLLECT ADVANCED SETTINGS
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
            with open(advs, "r") as stream:
                advs = yaml.safe_load(stream)

        vltRaw = advs["vlt-data-raw"]
        vltReduced = advs["vlt-data-reduced"]

        # TEST THE VLT FOLDERS EXIST
        if not os.path.exists(vltRaw) or not os.path.exists(vltReduced):
            print(
                "The VLT data structure does not seem to exist on this machine. Are you sure you need to use the --vlt flag?"
            )
            sys.exit(0)

        try:
            os.symlink(vltRaw, self.rawDir)
        except:
            os.unlink(self.rawDir)
            os.symlink(vltRaw, self.rawDir)

        self.log.debug("completed the ``use_vlt_environment_folders`` method")
        return vltReduced

    def generate_sof_names_and_maps(
        self,
        reductionOrder,
        rawFrames,
        calibrationFrames,
        calibrationTables,
        rawGroups,
        recipe,
        slit=False,
        slitmasks=False,
        techs=False,
        arm=False,
        binning=False,
        rospeed=False,
    ):
        """*generate sof names and maps*

        **Key Arguments:**
            - ``reductionOrder`` -- what files are we generating SOF sets for?
            - ``rawFrames`` -- dataframe of raw frames
            - ``calibrationFrames`` -- dataframe of calibration frames
            - ``calibrationTables`` -- dataframe of calibration tables
            - ``rawGroups`` -- raw frames grouped together
            - ``recipe`` -- recipe name
            - ``slit`` -- optional slit to filter by
            - ``slitmasks`` -- optional slit masks to filter by
            - ``techs`` -- optional techs to filter by
            - ``arm`` -- optional arm to filter by
            - ``binning`` -- optional binning to filter by
            - ``rospeed`` -- optional readout speed to filter by

        **Return:**
            - ``rawGroups`` -- updated rawGroups

        **Usage:**

        ```python
        rawGroups = self.generate_sof_names_and_maps(reductionOrder=reductionOrder, rawFrames=rawFrames, calibrationFrames=calibrationFrames, calibrationTables=calibrationTables, rawGroups=rawGroups)
        ```
        """
        self.log.debug("starting the ``generate_sof_names_and_maps`` method")

        from collections import defaultdict
        import pandas as pd

        self.sofMaps = []
        sofName = []
        armKey = self.kw("SEQ_ARM").lower()

        # INITIAL FILTER MASKS
        groupMask = rawGroups[armKey] == arm  # ALL TRUE
        rawMask = rawFrames[armKey] == arm  # ALL TRUE
        calibrationFrameMask = calibrationFrames[armKey] == arm  # ALL TRUE
        calibrationTableMask = calibrationTables[armKey] == arm  # ALL TRUE
        # REMOVE THIS WHEN PROCESSING ACQUISITION FRAMES
        groupMask = groupMask & (~rawGroups["eso seq arm"].str.contains("ACQ|--", case=False, regex=True))

        ## REDUCTION ORDER
        groupMask = groupMask & rawGroups["eso dpr type"].str.contains(reductionOrder, case=False)
        if "FLAT" in reductionOrder.upper():
            rawMask = rawMask & rawFrames["eso dpr type"].str.contains("FLAT")
        else:
            rawMask = rawMask & rawFrames["eso dpr type"].str.contains(reductionOrder, case=False)

        filteringMap = {
            "eso seq arm": {
                "applyMask": ["sofName", "rawFrames", "calibrationFrames", "calibrationTables"],
                "value": arm,
            },
            "binning": {"applyMask": ["sofName", "rawFrames", "calibrationFrames"], "value": binning},
            "rospeed": {"applyMask": ["sofName", "rawFrames", "calibrationFrames"], "value": rospeed},
            "slit": {"applyMask": ["sofName", "rawFrames"], "value": slit},
        }

        for key, value in filteringMap.items():
            applyMasks = value["applyMask"]
            value = value["value"]
            if value:
                print(value)
                groupMask = groupMask & (rawGroups[key] == value)

                if "rawFrames" in applyMasks:
                    rawMask = rawMask & (rawFrames[key] == value)
                if "calibrationFrames" in applyMasks:
                    calibrationFrameMask = calibrationFrameMask & (calibrationFrames[key] == value)
                if "calibrationTables" in applyMasks:
                    calibrationTableMask = calibrationTableMask & (calibrationTables[key] == value)
                if "sofName" in applyMasks and value != "--":
                    if isinstance(value, str):
                        sofName.append(value.upper())
                    else:
                        sofName.append(value)

        ## LIST FILTERS
        listFilters = {
            "slitmask": slitmasks,
            "eso dpr tech": techs,
        }
        for key, value in listFilters.items():
            if value:
                groupMask = groupMask & (rawGroups[key].isin(value))

        print(reductionOrder, arm, binning, rospeed, slit, slitmasks, techs)

        rawGroupsFiltered = rawGroups.loc[groupMask]
        calibrationFramesFiltered = calibrationFrames.loc[calibrationFrameMask]
        calibrationTablesFiltered = calibrationTables.loc[calibrationTableMask]

        # FILTER RAW FRAMES

        if slit:
            rawMask = rawMask & (rawFrames["slit"] == slit)
        # RETURN HERE: broken
        rawFramesFiltered = rawFrames.loc[rawMask]

        # GET UNIQUE VALUES IN COLUMN
        uniqueTemplateNames = rawGroupsFiltered["template"].unique()
        if len(uniqueTemplateNames):
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            mask2 = rawFramesFiltered["template"].isin(uniqueTemplateNames)
            rawFramesFiltered = rawFramesFiltered.loc[mask2]

        if not len(rawGroupsFiltered.index):
            self.log.debug("completed the ``generate_sof_names_and_maps`` method")
            return rawGroups

        # sofName = ("_").join(sofName).replace("-", "").replace(",", "_").upper() + ".sof"
        sofName = ("_").join(sofName).replace("-", "").replace(",", "_").upper()
        sofName = sofName.replace("XSHOOTER", "XSH")
        sofName = sofName.replace("FAST", "F")
        sofName = sofName.replace("SLOW", "S")
        sofName = sofName.replace("TELLURIC", "TELL")
        sofName = sofName.replace("ORDER_LOCATIONS", "OLOC")
        sofName = sofName.replace("DISP_SOL", "DSOL")
        sofName = sofName.replace("SPAT_SOL", "SSOL")
        self.sofName = sofName

        # OPTIMISE: 95%
        rawGroupsFiltered = rawGroupsFiltered.apply(
            self._generate_sof_name_and_map,
            axis=1,
            rawFrames=rawFramesFiltered,
            calibrationFrames=calibrationFramesFiltered,
            calibrationTables=calibrationTablesFiltered,
            recipe=recipe,
        )

        rawGroups.loc[groupMask] = rawGroupsFiltered

        if len(self.sofMaps):
            # FLATTENING THE LIST OF DICTIONARIES
            flattened_dict = defaultdict(list)
            for record in self.sofMaps:
                for key in record:
                    flattened_dict[key].extend(record[key])
            # CONVERT DEFAULTDICT BACK TO A REGULAR DICT (OPTIONAL)
            self.sofMaps = pd.DataFrame(dict(flattened_dict))
            keepTrying = 0
            while keepTrying < 7:
                try:
                    self.sofMaps.to_sql(f"sof_map_{self.sessionId}", con=self.conn, index=False, if_exists="append")
                    keepTrying = 10
                except Exception as e:
                    if keepTrying > 5:
                        raise Exception(e)
                    keepTrying += 1

        self.log.debug("completed the ``generate_sof_names_and_maps`` method")
        return rawGroups

    def populate_products_table(self, reductionOrder, rawGroups):
        """*populate products table*

        **Key Arguments:**
            - ``reductionOrder`` -- what files are we generating SOF sets for?
            - ``rawGroups`` -- raw frames grouped together

        **Return:**
            - ``rawGroups`` -- updated rawGroups

        **Usage:**

        ```python
        rawGroups = self.populate_products_table(reductionOrder=reductionOrder, rawGroups=rawGroups)
        ```
        """
        self.log.debug("starting the ``populate_products_table`` method")

        from collections import defaultdict
        import pandas as pd
        import time

        self.products = []

        # FIRST CREATE THE MASK
        mask = (
            (rawGroups["eso dpr type"].str.contains(reductionOrder, case=False))
            & ~(rawGroups["eso seq arm"].str.contains("ACQ|--", case=False, regex=True))
            & (rawGroups["complete"] == 1)
        )
        rawGroupsFiltered = rawGroups.loc[mask]

        # RETURN HERE: _populate_products_table_by_row does not return anything new (just self.products)
        rawGroupsFiltered = rawGroupsFiltered.apply(
            self._populate_products_table_by_row, axis=1, reductionOrder=reductionOrder
        )
        rawGroups.loc[mask] = rawGroupsFiltered

        if len(self.products):
            # FLATTENING THE LIST OF DICTIONARIES
            flattened_dict = defaultdict(list)
            for record in self.products:
                for key in record:

                    flattened_dict[key].extend(record[key])
            # CONVERT DEFAULTDICT BACK TO A REGULAR DICT (OPTIONAL)
            self.products = pd.DataFrame(dict(flattened_dict))
            # REMOVE COLUMN FROM DATA FRAME
            self.products.drop(
                columns=["eso dpr catg", "eso dpr tech", "eso dpr type", "counts", "complete"], inplace=True
            )

            # RETURN HERE: does this still need done?
            try:
                self.products.drop(columns=["product", "command"], inplace=True)
            except:
                pass

            keepTrying = 0
            while keepTrying < 7:
                try:
                    self.products.replace(["--"], None).to_sql(
                        "product_frames", con=self.conn, index=False, if_exists="append"
                    )
                    keepTrying = 10
                except Exception as e:

                    if keepTrying > 5:
                        raise Exception(e)
                    time.sleep(1)
                    keepTrying += 1

        self.log.debug("completed the ``populate_products_table`` method")
        return rawGroups

    def _fits_files_exist(self):
        """Check if any FITS files exist in rawDir or rootDir."""
        fitsExist = False
        exists = os.path.exists(self.rawDir)
        if exists:
            from fundamentals.files import recursive_directory_listing

            theseFiles = recursive_directory_listing(
                log=self.log, baseFolderPath=self.rawDir, whatToList="files"  # all | files | dirs
            )
            for f in theseFiles:
                if os.path.splitext(f)[1] == ".fits" or ".fits.gz" in os.path.splitext(f):
                    fitsExist = True
                    break
        if not fitsExist:
            for d in os.listdir(self.rootDir):
                filepath = os.path.join(self.rootDir, d)
                if os.path.isfile(filepath) and (
                    os.path.splitext(filepath)[1] == ".fits" or ".fits.gz" in os.path.splitext(filepath)
                ):
                    fitsExist = True
                    break
        return fitsExist

    def _get_or_create_db_connection(self):
        """Private method to get or create the SQLite database connection, copying the template if missing."""
        import shutil
        import sqlite3 as sql
        import time

        reset = False

        conn = None
        i = 0

        while i < 4:
            if not conn:
                try:
                    if self.conn:
                        conn = self.conn
                except:
                    pass

            if not conn:
                try:
                    with open(self.rootDbPath):
                        pass
                    self.freshRun = False
                except IOError:
                    # try:
                    #     os.remove(self.sessionIdFile)
                    # except Exception:
                    #     pass
                    self.freshRun = True
                    emptyDb = os.path.dirname(os.path.dirname(__file__)) + "/resources/soxspipe.db"
                    shutil.copyfile(emptyDb, self.rootDbPath)
                conn = sql.connect(self.rootDbPath, timeout=30, autocommit=True, check_same_thread=False)
            c = conn.cursor()

            try:
                c.execute("PRAGMA integrity_check;")
                this = c.fetchall()
                i = 10
            except Exception as e:
                # DATABASE IS BROKEN, REPLACE WITH EMPTY ONE
                i += 1
                c.close()
                conn.close()
                try:
                    del conn
                    del self.conn
                except:
                    pass

                time.sleep(1)

                if i > 3:
                    self.prepare(refresh=True)
                    if not reset:
                        reset = True
                conn = sql.connect(self.rootDbPath, timeout=30, autocommit=True)
                c = conn.cursor()

        c.execute("PRAGMA busy_timeout=10")
        c.close()

        return conn, reset

    # METHOD TO RETURN ALL THE RAW FRAME SETS THAT ARE NOT COMPLETE FROM THE DATABASE AS A PANDAS TABLE
    def get_incomplete_raw_frames_set(self):
        import pandas as pd

        query = 'select distinct "eso seq arm", round("mjd-obs",1) as "mjd-obs", "eso dpr tech","eso dpr type","slit","eso obs name","eso obs id"  from raw_frame_sets where complete = 0 and recipe in ("nod-obj","stare-obj","offset-obj")'
        return pd.read_sql(query, con=self.conn)

    def _select_instrument(self, inst=False):
        """Select the instrument and set related attributes."""
        from soxspipe.commonutils import keyword_lookup
        import yaml

        if inst:
            self.instrument = inst
        else:
            try:
                c = self.conn.cursor()
                sqlQuery = "select instrume from raw_frames where instrume is not null limit 1"
                c.execute(sqlQuery)
                self.instrument = c.fetchall()[0][0]
                c.close()
            except:
                return

        if "SOXS" not in self.instrument.upper():
            self.instrument = "XSH"

        # SETUP THE KEYWORD LOOKUP FUNCTION
        self.kw = keyword_lookup(log=self.log, instrument=self.instrument).get

        # SETUP SOF MAP
        yamlFilePath = (
            os.path.dirname(os.path.dirname(__file__)) + "/resources/" + self.instrument.lower() + "_sof_map.yaml"
        )

        # YAML CONTENT TO DICTIONARY
        with open(yamlFilePath, "r") as stream:
            self.sofMapLookup = yaml.safe_load(stream)

        # RETURN HERE .. delete
        if "SOXS" in self.instrument.upper():
            self.typeMap = self.typeMapSOXS
        else:
            self.typeMap = self.typeMapXSH
            if "EXPTIME2" in self.keyword_lookups:
                self.keyword_lookups.remove("EXPTIME2")

    def _flag_files_to_ignore(self):
        """*Flag files to ignore based on settings and reduction order*"""
        import sqlite3 as sql

        # FLAG FILES TO IGNORE BASED ON REDUCTION ORDER
        c = self.conn.cursor()
        theseKeywords = "','".join(self.reductionOrder)
        sqlQuery = f"update raw_frames set ignore = 1 WHERE `eso dpr type` not in ('{theseKeywords}')"
        c.execute(sqlQuery)
        self.conn.commit()

        # FLAG SIMULATION FILES TO IGNORE
        sqlQuery = f"update raw_frames set ignore = 1 WHERE `eso dpr type` not like '%OBJECT%' and `eso dpr type` not like '%STD%' and simulation = 1"
        c.execute(sqlQuery)
        self.conn.commit()

        # FLAG DFLATS TO IGNORE IF SPECIFIED IN SETTINGS
        if "ignore-dflats" in self.settings and self.settings["ignore-dflats"]:
            sqlQuery = "update raw_frames set ignore = 1 WHERE `eso dpr type` like '%DFLAT%'"
            c.execute(sqlQuery)
            self.conn.commit()

        # FLAG DOME FLATS TO IGNORE IF SPECIFIED IN SETTINGS
        if "ignore-dome-flats" in self.settings and self.settings["ignore-dome-flats"]:
            sqlQuery = "update raw_frames set ignore = 1 WHERE `eso dpr type` like '%DOME,FLAT%'"
            c.execute(sqlQuery)
            self.conn.commit()

        sqlQueries = ["update raw_frames set ignore = 1 WHERE `slit` = 'UNDEFINED'"]
        sqlQueries.append("update raw_frames set ignore = 1 WHERE `eso dpr type` like '%FLAT%' and `slit` = 'BLANK'")
        sqlQueries.append("update raw_frames set ignore = 1 WHERE `eso seq arm` = 'ACQ'")
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)
            self.conn.commit()
        c.close()

    def build_sof_files(self):
        """*scan the raw frame table to generate the listing of products that are expected to be created and then write out all of the needed SOF files*

        **Usage:**

        ```python
        self.build_sof_files()
        ```
        """
        self.log.debug("starting the ``_populate_product_frames_db_table`` method")

        import pandas as pd

        c = self.conn.cursor()
        sqlQuery = f"update product_frames set status = status_{self.sessionId};"
        c.execute(sqlQuery)
        sqlQuery = f"update raw_frames set processed = 0 where processed < 0;"
        c.execute(sqlQuery)

        # CLEAN UP FAILED FILES
        # DELETE FROM
        count = 0
        oldCount = -1
        while count != oldCount:
            oldCount = count
            c = self.conn.cursor()
            sqlQuery = f"select distinct sof from sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1));"
            compromisedSofs = pd.read_sql(sqlQuery, con=self.conn)["sof"].tolist()
            count = len(compromisedSofs)

            sqlQuery = f"update product_frames set complete = 0 where (status != 'fail' or status is null) and sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1)));"
            c.execute(sqlQuery)

        sqlQueries = [
            f"update raw_frames set processed = 0 where file in (select file from sof_map where sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1))));",
            f"update raw_frame_sets set complete = 0 where sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1)));",
            f"delete from sof_map_{self.sessionId} where sof in (  select s.sof from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1));",
            f"update raw_frames set processed = -1 where file in (select distinct s.file from sof_map s, product_frames p where p.sof=s.sof and p.status = 'fail');",
        ]
        for sqlQuery in sqlQueries:
            c = self.conn.cursor()
            c.execute(sqlQuery)
            c.close()

        # DELETE COMPROMISED SOF FILES
        for sof in compromisedSofs:
            sofPath = self.sessionPath + "/sof/" + sof
            try:
                os.remove(sofPath)
            except:
                pass

        # RESET ALL PRODUCTS TO INCOMPLETE
        c = self.conn.cursor()

        allRawGroups = []
        for recipeOrder, (name, filters) in enumerate(self.sofMapLookup.items()):

            # READ FILTER VARIABLES
            ttypes = filters["eso dpr type"]
            tech = filters["eso dpr tech"]
            productTypes = filters["products"] or []
            recipe = filters["recipe"]
            if "calibrations" in filters:
                calibrationTypes = filters["calibrations"]
            else:
                calibrationTypes = []

            for ttype in ttypes:

                rawFrames, rawGroups = self.get_raw_frames_and_groups(
                    ttype=ttype,
                    tech=tech,
                    recipe=recipe,
                    recipeOrder=recipeOrder,
                    filterName=name,
                    unprocessedOnly=True,
                )

                # ADD PREDICTED PRODUCT TO PRODUCT TABLE - DETERMINE IF COMPLETE LATER
                incompleteProducts = self.predict_product_frames(productTypes, rawGroups, recipe)

                if not incompleteProducts:
                    continue

                allRawGroups.append(rawGroups)

                if not len(calibrationTypes):
                    # MBIAS AND MDARK -- ALWAYS COMPLETE (NO PRIOR CALIBRATION REQUIRED)
                    sqlQuery = f"select sof from product_frames where recipe = '{recipe}' and complete = 0;"
                    containerSofs = pd.read_sql(sqlQuery, con=self.conn)["sof"].tolist()
                    self.raw_frames_to_sof_map(rawGroups=rawGroups, containerSofs=containerSofs)
                    sqlQuery = f"update product_frames set complete = 1 where recipe = '{recipe}' and complete = 0;"
                    c.execute(sqlQuery)

                else:

                    if isinstance(calibrationTypes, dict):
                        for arm, calType in calibrationTypes.items():
                            ffrom = ", ".join([f"cal_{ct}" for ct in calType])
                            where1 = " and ".join([f"product_frames.sof=cal_{ct}.sof" for ct in calType])
                            where2 = " and ".join(
                                [
                                    f"(cal_{ct}.upstream_status = 'pass' or cal_{ct}.upstream_status is null)"
                                    for ct in calType
                                ]
                            )

                            sqlQuery = f"""select sof from product_frames where uuid in 
                            (select product_frames.uuid from product_frames, {ffrom} where product_frames.recipe = '{recipe}' and product_frames.'eso seq arm' = '{arm}' and {where1} and {where2}) and complete < 1;"""

                            containerSofs = pd.read_sql(sqlQuery, con=self.conn)["sof"].tolist()

                            self.raw_frames_to_sof_map(rawGroups=rawGroups, containerSofs=containerSofs)

                            sqlQuery = f"""update product_frames set complete = -1 where uuid in 
                            (select product_frames.uuid from product_frames, {ffrom} where product_frames.recipe = '{recipe}' and product_frames.'eso seq arm' = '{arm}' and {where1} and {where2}) and complete < 1;
                            """
                            c.execute(sqlQuery)

                            # FOR COMPLETE PRODUCTS, ADD CALIBRATION FILES TO SOF MAP
                            # NEED TO ALSO ADD THE RAW FILES TOO ... ADD RAW FRAMES, SET COMPLETE = 1 WHERE PRODUCT FRAMES COMPLETE = 1
                            for ct in calType:

                                sqlQuery = f"""select cal_{ct}.file, cal_{ct}.upstream_tag as tag, product_frames.sof, cal_{ct}.filepath, product_frames.complete from product_frames, cal_{ct} where product_frames.complete = -1 and product_frames.sof=cal_{ct}.sof;"""
                                newSof = pd.read_sql(sqlQuery, con=self.conn)

                                if len(newSof):
                                    self._dataframe_to_sqlite(newSof, f"sof_map_{self.sessionId}", replace=False)

            sqlQuery = f"""update product_frames set complete = 1 where complete = -1;"""
            c.execute(sqlQuery)

        self.conn.commit()
        c.close()

        # CONCAT ALL RAW GROUPS
        if len(allRawGroups):
            rawGroups = pd.concat(allRawGroups, ignore_index=True)
            if len(rawGroups.index):
                rawGroups = rawGroups.drop(columns=["filepaths"])
                self._dataframe_to_sqlite(rawGroups, "raw_frame_sets", replace=False)

        c = self.conn.cursor()
        sqlQueries = [
            f"update sof_map_{self.sessionId} set complete = 1 where complete = -1;",
            "update raw_frame_sets set complete = 1 where sof in (select r.sof from raw_frame_sets r, product_frames p where p.sof = r.sof and p.complete = 1);",
            "update raw_frame_sets set complete = 0 where sof in (select r.sof from raw_frame_sets r, product_frames p where p.sof = r.sof and p.complete = 0);",
            "update raw_frames set processed = 1 where processed = 0 and filepath in (select filepath from sof_map);",
        ]
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)
        self.conn.commit()
        c.close()

        self._write_sof_files()

        return

    def get_raw_frames_and_groups(
        self, ttype=None, arm=None, tech=None, recipe=None, recipeOrder=None, filterName=None, unprocessedOnly=False
    ):
        """*Process raw frames to group and calculate mean MJD values.*

        **Key Arguments:**
            - ``ttype`` -- optional data product `eso dpr type` to filter by
            - ``arm`` -- optional instrument `eso seq arm` to filter by
            - ``tech`` -- optional list of `eso dpr tech` to filter by
            - ``recipe`` -- recipe name to assign to groups
            - ``recipeOrder`` -- recipe reduction order
            - ``filterName`` -- optional name to filter the groups
            - ``unprocessedOnly`` -- if True, only return unprocessed raw frames

        **Return:**
            - `rawFrames` -- processed raw frames dataframe
            - `rawGroups` -- grouped raw frames with calculated MJD values

        **Usage:**

        ```python
        rawFrames, rawGroups = self.get_raw_frames_and_groups()
        ```
        """
        import pandas as pd

        # CLEAN DATABASE TABLE
        c = self.conn.cursor()
        sqlQueries = [
            "update raw_frames set lamp = null, slit = null, slitmask = null where `eso dpr type` in ('BIAS','DARK');",
            "update raw_frames set rospeed = null where rospeed = -1;",
        ]
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)
        c.close()

        # IF NONE, SET TO EMPTY STRING
        ttype, arm, tech = ttype or "", arm or "", tech or ""

        if ttype or arm:
            where = "where"
        else:
            where = ""
        if ttype:
            ttype = "and `eso dpr type` = '" + ttype + "'"
        if arm:
            arm = "and `eso seq arm` = '" + arm + "'"
        if tech:
            # JOIN ITEMS IN TECH LIST TO A COMMA-SEPARATED STRING
            tech = "and `eso dpr tech` in (" + ",".join(["'" + t + "'" for t in tech]) + ")"

        if unprocessedOnly:
            where = where + " and processed = 0"

        # READ IN RAW FRAMES TABLE
        conn = self.conn
        rawFrames = pd.read_sql(
            f"SELECT * FROM raw_frames_valid {where} {ttype} {arm} {tech} order by `mjd-obs` asc".replace(
                "where and", "where"
            ),
            con=conn,
        )

        rawFrames = rawFrames.astype({"exptime": float, "ra": float, "dec": float})
        rawFrames.fillna({"exptime": -99.99, "ra": -99.99, "dec": -99.99}, inplace=True)
        rawFrames.fillna("--", inplace=True)
        filterKeywordsRaw = self.filterKeywords[:]
        for i in self.proKeywords:
            filterKeywordsRaw.remove(i)

        # HIDE OFF FRAMES FROM GROUPS
        mask = (
            (rawFrames["eso dpr tech"] == "IMAGE")
            & (rawFrames["eso seq arm"] == "NIR")
            & (rawFrames["eso dpr type"] != "DARK")
        )
        rawFramesNoOffFrames = rawFrames.loc[~mask]
        rawGroups = rawFramesNoOffFrames.groupby(filterKeywordsRaw)

        if not len(rawFramesNoOffFrames.index):
            return pd.DataFrame(), pd.DataFrame()

        mjds = rawGroups.mean(numeric_only=True)["mjd-obs"].values
        # MANIPULATE ALL STRINGS IN startTime TO BE OF FORMAT YYYY-MM-DDTHH:MM:SS.SSS
        startTime = rawGroups.min()["date-obs"].values
        filepaths = [list(group["filepath"].values) for name, group in rawGroups]

        # PRINT EACH GROUP IN THE GROUPED DATAFRAME

        startTime = [str(s).split(".")[0].replace("-", "").replace(":", "") for s in startTime]
        rawGroups = rawGroups.size().reset_index(name="counts")
        rawGroups["mjd-obs"] = mjds
        rawGroups["date-obs"] = startTime
        rawGroups["filepaths"] = filepaths

        # REMOVE GROUPED STARE - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = rawGroups["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"])
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD SCIENCE FRAMES AS ONE ENTRY PER EXPOSURE
        rawScienceFrames = rawFrames.loc[rawFrames["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"])]
        if len(rawScienceFrames.index):
            rawScienceFrames = rawScienceFrames.groupby(filterKeywordsRaw + ["mjd-obs", "date-obs"])
            filepaths = [list(group["filepath"].values) for name, group in rawScienceFrames]
            rawScienceFrames = rawScienceFrames.size().reset_index(name="counts")
            rawScienceFrames["filepaths"] = filepaths
            startTime = rawScienceFrames["date-obs"].values
            startTime = [str(s).split(".")[0].replace("-", "").replace(":", "") for s in startTime]
            rawScienceFrames["date-obs"] = startTime
            # MERGE DATAFRAMES
            rawGroups = pd.concat([rawGroups, rawScienceFrames], ignore_index=True)

        # REMOVE GROUPED SINGLE PINHOLE ARCS - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = rawGroups["eso dpr tech"].isin(["ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE"])
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD PINHOLE FRAMES AS ONE ENTRY PER EXPOSURE
        if self.instrument.upper() == "SOXS":
            rawPinholeFrames = rawFrames.loc[
                (rawFrames["eso dpr tech"].isin(["ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE"]))
                & (
                    (rawFrames["eso seq arm"] == "NIR")
                    | (~rawFrames["lamp"].isin(["Xe", "Ar", "Hg", "Ne", "ArNeHgXe"]))
                )
            ]
        else:
            rawPinholeFrames = rawFrames.loc[
                rawFrames["eso dpr tech"].isin(["ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE"])
            ]
        if len(rawPinholeFrames.index):
            rawPinholeFrames = rawPinholeFrames.groupby(filterKeywordsRaw + ["mjd-obs", "date-obs"])
            filepaths = [list(group["filepath"].values) for name, group in rawPinholeFrames]
            rawPinholeFrames = rawPinholeFrames.size().reset_index(name="counts")
            rawPinholeFrames["filepaths"] = filepaths
            startTime = rawPinholeFrames["date-obs"].values
            startTime = [str(s).split(".")[0].replace("-", "").replace(":", "") for s in startTime]
            rawPinholeFrames["date-obs"] = startTime
            # MERGE DATAFRAMES
            rawGroups = pd.concat([rawGroups, rawPinholeFrames], ignore_index=True)

        rawGroups["recipe"] = recipe
        if "STD,FLUX" in ttype:
            recipe += "_std_flux"
        rawGroups["sof"] = (
            rawGroups["date-obs"].astype(str)
            + "_"
            + rawGroups["eso seq arm"].astype(str)
            + "_"
            + rawGroups["binning"].astype(str)
            + "_"
            + rawGroups["rospeed"].astype(str)
            + "_"
            + recipe
            + "_"
            + rawGroups["lamp"].astype(str)
            + "_"
            + rawGroups["slit"].astype(str)
            + "_"
            + rawGroups["exptime"].astype(str)
            + "s"
            + "_"
            + rawGroups["instrume"].astype(str)
        )
        rawGroups["sof"] = rawGroups["sof"].str.upper()
        if "science" in filterName.lower():
            rawGroups["sof"] += "_" + rawGroups["object"].astype(str).replace("-", "_")
        for _ in range(5):
            rawGroups["sof"] = rawGroups["sof"].str.replace("_--_", "_", regex=False)
            rawGroups["sof"] = rawGroups["sof"].str.replace(" ", "", regex=False)
            rawGroups["sof"] = rawGroups["sof"].str.replace(".", "_", regex=False)
            rawGroups["sof"] = rawGroups["sof"].str.replace("__", "_", regex=False)
        rawGroups["sof"] += ".sof"

        rawGroups["sof"] = rawGroups["sof"].str.replace("MBIAS_0_0S_", "MBIAS_")
        rawGroups["sof"] = rawGroups["sof"].str.replace("DISP_SOLUTION.*?_PINHOLE", "DSOL_PINHOLE", regex=True)
        rawGroups["sof"] = rawGroups["sof"].str.replace("SPAT_SOLUTION.*?_MULTPIN", "SSOL_MULTPIN", regex=True)
        rawGroups["sof"] = rawGroups["sof"].str.replace("ORDER_CENTRES", "OLOC", regex=True)
        if recipe in ["mbias", "mdark"]:
            rawGroups["complete"] = 1
        else:
            rawGroups["complete"] = 0
        rawGroups["recipe_order"] = recipeOrder

        return rawFrames, rawGroups

    def get_calibration_frames(self, recipes=None):
        """
        Fetch calibration frames and tables from the database.

        **Key Arguments:**
        - `recipes` -- Optional list of recipe types to filter calibration frames.

        **Return:**
        - `calibrationFrames` -- DataFrame containing calibration frames.
        - `calibrationTables` -- DataFrame containing calibration tables.
        """
        import pandas as pd

        if recipes and len(recipes):
            recipes = f"and recipe in (\'{'\',\''.join(recipes)}\')"
        else:
            recipes = ""

        conn = self.conn

        # RETURN HERE
        calibrationFrames = pd.read_sql(
            f"SELECT * FROM product_frames WHERE `eso pro catg` NOT LIKE '%_TAB_%' AND (status_{self.sessionId} != 'fail' OR status_{self.sessionId} IS NULL ) {recipes}",
            con=conn,
        )
        calibrationFrames = pd.read_sql(
            f"SELECT * FROM product_frames WHERE (status_{self.sessionId} != 'fail' OR status_{self.sessionId} IS NULL ) {recipes}",
            con=conn,
        )

        calibrationFrames.fillna("--", inplace=True)
        self.productFilterKeywords = [
            "eso seq arm",
            "eso pro catg",
            "eso pro tech",
            "eso pro type",
        ]
        calibrationFrames = calibrationFrames.groupby(self.productFilterKeywords).size().reset_index(name="counts")

        calibrationTables = pd.read_sql(
            f"SELECT * FROM product_frames WHERE `eso pro catg` LIKE '%_TAB_%' AND (status_{self.sessionId} != 'fail' OR status_{self.sessionId} IS NULL) {recipes}",
            con=conn,
        )
        calibrationTables.fillna("--", inplace=True)

        return calibrationFrames, calibrationTables

    def predict_product_frames(self, productTypes, rawGroups, recipe):
        """
        Process product frames for a given set of product types and raw groups.

        **Key Arguments:**
        - `productTypes` -- List of product types to process.
        - `rawGroups` -- DataFrame containing raw groups.
        - `recipe` -- Recipe name.

        **Return:**
        - `incompleteProducts` -- Number of incomplete products.
        """
        import pandas as pd

        if not len(rawGroups.index):
            sqlQuery = f"select count(*) from product_frames where recipe = '{recipe}' and complete< 1;"
            c = self.conn.cursor()
            c.execute(sqlQuery)
            incompleteProducts = c.fetchall()[0][0]
            c.close()
            return incompleteProducts

        for product in productTypes:
            proKeys = list(product.values())[0]
            product = list(product.keys())[0]
            productFrames = rawGroups.copy()

            productFrames.drop(
                columns=[
                    "eso dpr catg",
                    "eso dpr tech",
                    "eso dpr type",
                    "counts",
                    "complete",
                    "date-obs",
                    "filepaths",
                ],
                inplace=True,
            )

            productFrames["eso pro type"] = proKeys["eso pro type"]
            productFrames["eso pro tech"] = proKeys["eso pro tech"]
            productFrames["eso pro catg"] = proKeys["eso pro catg"]
            productFrames["eso pro catg"] = (
                productFrames["eso pro catg"].astype(str) + "_" + productFrames["eso seq arm"].astype(str).str.upper()
            )
            if product in ["fits image", "fits table"]:
                productFrames["file"] = productFrames["sof"].str.replace(".sof", f".fits")
                if "replace" in proKeys:
                    for item in proKeys["replace"]:
                        productFrames["file"] = productFrames["file"].str.replace(item["from"], item["to"])
            else:
                productFrames["file"] = "XXXX"

            productFrames["filepath"] = (
                "./reduced/"
                + productFrames["mjd-date"].astype(str)
                + "/soxs-"
                + recipe.replace("_", "-")
                + "/"
                + productFrames["file"].astype(str)
            )

            self._dataframe_to_sqlite(productFrames, "product_frames")

        return 1

    def raw_frames_to_sof_map(self, rawGroups, containerSofs):
        """
        Generate the SOF map from raw groups and complete product SOFs.

        **Key Arguments:**
        - `rawGroups` -- DataFrame containing raw frame groups.
        - `containerSofs` -- array of complete product SOFs.

        **Return:**
        - `sofMapDF` -- DataFrame containing the generated SOF map.
        """
        import pandas as pd

        if not len(rawGroups.index):
            return

        # BUILD SOF MAP TABLE
        mask = rawGroups["sof"].isin(containerSofs)
        sofMapDF = rawGroups.loc[mask]
        sofMapDF = sofMapDF.explode("filepaths")
        sofMapDF["file"] = sofMapDF["filepaths"].apply(lambda x: os.path.basename(x) if pd.notnull(x) else x)
        sofMapDF = sofMapDF.rename(columns={"filepaths": "filepath"})
        sofMapDF["tag"] = sofMapDF["eso dpr type"].replace(",", "_") + "_" + sofMapDF["eso seq arm"]
        sofMapDF = sofMapDF[["file", "tag", "sof", "filepath", "complete"]]
        sofMapDF["complete"] = 1
        # ADD TO DATABASE
        self._dataframe_to_sqlite(sofMapDF, f"sof_map_{self.sessionId}", replace=False)

        # UPDATE RAW FRAMES AS PROCESSED
        processedRawFiles = sofMapDF["file"].unique().tolist()
        if len(processedRawFiles):
            c = self.conn.cursor()
            placeholders = ",".join(["?"] * len(processedRawFiles))
            sqlQuery = f"update raw_frames set processed=1 where file in ({placeholders});"
            c.execute(sqlQuery, processedRawFiles)
            self.conn.commit()
            c.close()

        return

    def _dataframe_to_sqlite(self, dataframe, table_name, replace=False):
        """
        Retry inserting into the database with a maximum of 7 attempts.

        **Key Arguments:**
        - `dataframe` -- DataFrame containing rows to insert.
        - `table_name` -- Name of the database table to insert into.
        - `replace` -- If True, replace existing entries; otherwise, append.

        **Raises:**
        - Exception if the insertion fails after 7 attempts.
        """
        import time

        if replace:
            c = self.conn.cursor()
            sqlQuery = f"delete from {table_name};"
            try:
                c.execute(sqlQuery)
            except:
                pass
            c.close()

        keepTrying = 0
        while keepTrying < 7:
            try:
                dataframe.replace(["--"], None).to_sql(table_name, con=self.conn, index=False, if_exists="append")
                keepTrying = 10
            except Exception as e:
                if keepTrying > 5:
                    raise Exception(e)
                time.sleep(1)
                keepTrying += 1


def _harvest_fits_headers(batch, log, pathToDirectory, keywords, filterKeys, instrument, kw):
    from ccdproc import ImageFileCollection
    import numpy as np
    from astropy.time import Time, TimeDelta

    masterTable = ImageFileCollection(location=pathToDirectory, filenames=batch, keywords=keywords)
    masterTable = masterTable.summary
    # ADD FILLED VALUES FOR MISSING CELLS
    for fil in keywords:
        if fil in filterKeys and fil not in ["exptime"]:

            try:
                masterTable[fil].fill_value = "--"
            except:
                masterTable.replace_column(fil, masterTable[fil].astype(str))
                masterTable[fil].fill_value = "--"
        # elif fil in ["exptime"]:
        #     masterTable[fil].fill_value = "--"
        else:
            try:
                masterTable[fil].fill_value = -99.99
            except:
                masterTable[fil].fill_value = "--"
    masterTable = masterTable.filled()

    # FIX ACQ CAM EXPTIME & ARM & FILTER
    if "SOXS" in instrument.upper():
        matches = (masterTable["exptime"] == -99.99) & (masterTable[kw("EXPTIME2").lower()] != -99.99)
        masterTable["exptime"][matches] = masterTable[kw("EXPTIME2").lower()][matches]
        matches = (masterTable["eso seq arm"] == "--") & (masterTable[kw("DET").lower()] == "ACQ")
        masterTable["eso seq arm"][matches] = "ACQ"
        matches = (masterTable["eso seq arm"] != "ACQ") & (masterTable[kw("ACFW_ID").lower()] != "--")
        masterTable[kw("ACFW_ID").lower()][matches] = "--"
        matches = (masterTable["eso seq arm"] == "ACQ") & (
            (masterTable["eso dpr type"] == "BIAS") | (masterTable["eso dpr type"] == "DARK")
        )
        masterTable[kw("ACFW_ID").lower()][matches] = "--"

    # FILTER OUT FRAMES WITH NO MJD
    matches = (
        (masterTable["mjd-obs"] == -99.99)
        | (masterTable["eso dpr catg"] == "--")
        | (masterTable["eso dpr tech"] == "--")
        | (masterTable["eso dpr type"] == "--")
        | (masterTable["exptime"] == -99.99)
    )
    missingMJDFiles = masterTable["file"][matches]
    if len(missingMJDFiles):
        print("\nThe following FITS files are missing DPR keywords and will be ignored:\n\n")
        print(missingMJDFiles)
        masterTable = masterTable[~matches]

    # SETUP A NEW COLUMN GIVING THE INT MJD THE CHILEAN NIGHT BEGAN ON
    # 12:00 NOON IN CHILE IS TYPICALLY AT 16:00 UTC (CHILE = UTC - 4)
    # SO COUNT CHILEAN OBSERVING NIGHTS AS 15:00 UTC-15:00 UTC (11am-11am)
    if "mjd-obs" in masterTable.colnames:
        chile_offset = TimeDelta(4.0 * 60 * 60, format="sec")
        night_start_offset = TimeDelta(15.0 * 60 * 60, format="sec")
        mjd_ofset = TimeDelta(12.0 * 60 * 60, format="sec")
        masterTable["mjd-obs"] = masterTable["mjd-obs"].astype(float)
        chileTimes = Time(masterTable["mjd-obs"], format="mjd", scale="utc") - chile_offset
        startNightDate = Time(masterTable["mjd-obs"], format="mjd", scale="utc") - night_start_offset
        utcDate = Time(masterTable["mjd-obs"], format="mjd", scale="utc")
        # masterTable["utc-4hrs"] = (masterTable["mjd-obs"] - 2 / 3).astype(int)
        mjdDate = Time(masterTable["mjd-obs"], format="mjd", scale="utc")
        masterTable["mjd-date"] = mjdDate.strftime("%Y-%m-%d")
        masterTable["utc-4hrs"] = chileTimes.strftime("%Y-%m-%dt%H:%M:%S")
        masterTable["night start date"] = startNightDate.strftime("%Y-%m-%d")
        masterTable["night start mjd"] = startNightDate.mjd.astype(int)
        masterTable["boundary"] = startNightDate.mjd - startNightDate.mjd.astype(int)
        masterTable.add_index("night start date")
        masterTable.add_index("night start mjd")
        masterTable.add_index("mjd-date")

    if instrument.upper() != "SOXS":
        if kw("DET_READ_SPEED").lower() in masterTable.colnames:
            masterTable["rospeed"] = np.copy(masterTable[kw("DET_READ_SPEED").lower()])
            try:
                masterTable["rospeed"][masterTable["rospeed"] == -99.99] = "--"
            except:
                masterTable["rospeed"] = masterTable["rospeed"].astype(str)
                masterTable["rospeed"][masterTable["rospeed"] == -99.99] = "--"
            masterTable["rospeed"][masterTable["rospeed"] == "1pt/400k/lg"] = "fast"
            masterTable["rospeed"][masterTable["rospeed"] == "1pt/400k/lg/AFC"] = "fast"
            masterTable["rospeed"][masterTable["rospeed"] == "1pt/100k/hg"] = "slow"
            masterTable["rospeed"][masterTable["rospeed"] == "1pt/100k/hg/AFC"] = "slow"
            masterTable.add_index("rospeed")
    else:
        if kw("DET_READ_SPEED").lower() in masterTable.colnames:
            masterTable["rospeed"] = np.copy(masterTable[kw("DET_READ_SPEED").lower()])

            try:
                masterTable["rospeed"][masterTable["rospeed"] == -99.99] = -1
            except:
                masterTable["rospeed"] = masterTable["rospeed"].astype(str)
                masterTable["rospeed"][masterTable["rospeed"] == -99.99] = -1

            masterTable.add_index("rospeed")

    if kw("TPL_ID").lower() in masterTable.colnames:
        masterTable["template"] = np.copy(masterTable[kw("TPL_ID").lower()])

    if kw("ACFW_ID").lower() in masterTable.colnames:
        masterTable["filter"] = np.copy(masterTable[kw("ACFW_ID").lower()])

    if "naxis" in masterTable.colnames:
        masterTable["table"] = np.copy(masterTable["naxis"]).astype(str)
        masterTable["table"][masterTable["table"] == "0"] = "T"
        masterTable["table"][masterTable["table"] != "T"] = "F"

    if kw("WIN_BINX").lower() in masterTable.colnames:
        masterTable["binning"] = np.core.defchararray.add(
            masterTable[kw("WIN_BINX").lower()].astype("int").astype("str"), "x"
        )
        masterTable["binning"] = np.core.defchararray.add(
            masterTable["binning"], masterTable[kw("WIN_BINY").lower()].astype("int").astype("str")
        )
        masterTable["binning"][masterTable["binning"] == "-99x-99"] = "--"
        masterTable["binning"][masterTable["binning"] == "1x-99"] = "--"
        masterTable.add_index("binning")

    if kw("ABSROT").lower() in masterTable.colnames:
        masterTable["absrot"] = masterTable[kw("ABSROT").lower()].astype(float)
        masterTable.add_index("absrot")

    # ADD INDEXES ON ALL KEYS
    for k in keywords:
        try:
            masterTable.add_index(k)
        except:
            pass

    # SORT IMAGE COLLECTION
    masterTable.sort(
        [
            "eso pro type",
            "eso seq arm",
            "eso dpr catg",
            "eso dpr tech",
            "eso dpr type",
            "eso pro catg",
            "eso pro tech",
            "mjd-obs",
        ]
    )
    return masterTable.to_pandas(index=False)
