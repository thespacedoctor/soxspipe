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

    def __init__(self, log, rootDir, vlt=False, dbConnect=True):
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
            "GAIN",
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
            "PAR_ANG_END",
            "PAR_ANG_START",
            "AZ_ANG",
            "ALT_ANG",
            "SEEING_END",
            "SEEING_START",
            "AIRMASS_END",
            "AIRMASS_START",
            "TARG_NAME",
            "OB_TPL_NO",
            "OB_NTPL",
            "OB_START",
            "TPL_START",
            "NIR_TEMP_K",
            "VIS_TEMP_C",
            "CP_TEMP_C",
            "AFC1_POS1",
            "AFC1_POS2",
            "AFC2_POS1",
            "AFC2_POS2",
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
            "gain",
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
            "filepath",
            "eso tel parang end",
            "eso tel parang start",
            "eso tel az",
            "eso tel alt",
            "eso tel ambi fwhm end",
            "eso tel ambi fwhm start",
            "eso tel airm end",
            "eso tel airm start",
            "eso obs targ name",
            "eso obs tplno",
            "eso obs ntpl",
            "eso tpl start",
            "eso obs start",
            "nir temp k",
            "vis temp c",
            "cp temp c",
            "afc1 pos1",
            "afc1 pos2",
            "afc2 pos1",
            "afc2 pos2",
        ]

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
            "gain",
            "binning",
            "night start mjd",
            "night start date",
            "instrume",
            "lamp",
            "template",
            "eso obs name",
            "eso obs id",
            "eso tpl name",
            "eso tpl nexp",
            "filter",
            "object",
        ]

        self.filterKeywordsExtras = [
            "mjd-obs",
            "nir temp k",
            "vis temp c",
            "cp temp c",
            "afc1 pos1",
            "afc1 pos2",
            "afc2 pos1",
            "afc2 pos2",
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
        if exists and dbConnect:
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
            if False:
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

        if True:
            c = self.conn.cursor()
            sqlQueries = [
                f'update product_frames set status_{self.sessionId} = "pass" where status_{self.sessionId} = "fail" and sof in (select sof_name from quality_control);',
                f"update quality_control set qc_value_min = null, qc_value_max = null, qc_flag = 'pass';",
            ]

            for k, v in self.settings.items():
                if k[:5] == "soxs-":
                    recipe = k
                    for a in ["acq", "vis", "nir"]:
                        if a in v and "qc-acceptable-ranges" in v[a]:
                            arm = a.upper()
                            for kk, vv in v[a]["qc-acceptable-ranges"].items():
                                qc_name = kk.upper().replace("-", " ")
                                qc_min = vv[0]
                                qc_max = vv[1]
                                sqlQueries.append(
                                    f'update quality_control set qc_value_min = {qc_min}, qc_value_max = {qc_max} where soxspipe_recipe = "{recipe}" and sof_name like "%{arm}%" and qc_name = "{qc_name}"  and qc_order = "-1";'
                                )
                    if "qc-acceptable-ranges" in v:
                        for kk, vv in v["qc-acceptable-ranges"].items():
                            qc_name = kk.upper().replace("-", " ")
                            qc_min = vv[0]
                            qc_max = vv[1]
                            sqlQueries.append(
                                f'update quality_control set qc_value_min = {qc_min}, qc_value_max = {qc_max} where soxspipe_recipe = "{recipe}" and qc_name = "{qc_name}"  and qc_order = "-1";'
                            )

            sqlQueries += [
                'update quality_control set qc_flag = "pass" where CAST(qc_value as float) < qc_value_max and CAST(qc_value as float) > qc_value_min and qc_flag != "pass" and qc_order = "-1";',
                'update quality_control set qc_flag = "fail" where (CAST(qc_value as float) > qc_value_max or CAST(qc_value as float) < qc_value_min) and qc_flag != "fail" and qc_order = "-1";',
                'update product_frames set status_base = "fail" where sof in (select sof_name from quality_control where qc_flag = "fail");',
            ]

            for sqlQuery in sqlQueries:
                c.execute(sqlQuery)
            self.conn.commit()
            c.close()

        self._flag_files_to_ignore()
        self.build_sof_files()

        if report:

            rawDirStr = self.rawDir.replace("./", "")

            print(f"\nTHE `{basename}` WORKSPACE FOR HAS BEEN PREPARED FOR DATA-REDUCTION\n")
            print(f"In this workspace you will find:\n")
            print(f"   - `misc/`: a lost-and-found archive of non-fits files")
            print(f"   - `qc/`: nested folders, ordered by date, containing quality-control plots and tables.")
            print(f"   - `{rawDirStr}/`: nested folders, ordered by date, containing raw-frames.")
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

    def list_raw(self, sofFile):
        """*list the all the raw frames associated with a given science object SOF file*"""
        import pandas as pd

        self.log.debug("starting the ``list_raw`` method")

        sqlQuery = f"select sof from product_frames where sof = '{sofFile}' and complete = 1"

        for _ in range(4):  # Recursively query up to 5 times
            sqlQuery = f"SELECT distinct sof FROM product_frames WHERE file IN (SELECT file FROM sof_map_base WHERE sof in ({sqlQuery})) or sof in ({sqlQuery})"

        sqlQuery = f"SELECT filepath from sof_map WHERE sof in ({sqlQuery}) and filepath like '%./raw/%' order by sof"

        filepaths = pd.read_sql(sqlQuery, con=self.conn)["filepath"].tolist()

        self.log.debug("completed the ``list_raw`` method")
        return list(set(filepaths))

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
                    mask = rawFrames["filepath"].isnull()
                    rawFrames.loc[mask, "filepath"] = (
                        "./raw/" + rawFrames.loc[mask, "mjd-date"] + "/" + rawFrames.loc[mask, "file"]
                    )

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
                    rawFrames.replace(["--", -99.99], None).to_sql(
                        "raw_frames", con=self.conn, index=False, if_exists="append"
                    )

                    # MOVE THE FILES TO THE CORRECT LOCATION
                    filepaths = rawFrames["filepath"]
                    filenames = rawFrames["file"]
                    for p, n in zip(filepaths, filenames):
                        parentDirectory = os.path.dirname(p)
                        if not os.path.exists(parentDirectory):
                            # Recursively create missing directories
                            os.makedirs(parentDirectory)
                        if os.path.exists(self.rootDir + "/" + n):
                            realSource = os.path.realpath(self.rootDir + "/" + n)
                            realDest = os.path.realpath(p)

                            matchObject = re.match(r".*?(raw\/\d{4}-\d{2}-\d{2}.*)", realSource)

                            if matchObject and realSource != realDest:
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
        fitsPathsRel = []
        fitsNames = []
        for entry in os.scandir(pathToDirectory):
            if (
                not entry.name.startswith(".")
                and entry.is_file()
                and (os.path.splitext(entry.name)[1] == ".fits" or ".fits.Z" in entry.name)
            ):
                # fitsPaths.append(entry.path)
                if os.path.islink(entry.path):
                    fp = "./" + os.path.relpath(os.path.realpath(entry.path), pathToDirectory)

                else:
                    fp = os.path.relpath(entry.path, pathToDirectory)
                fitsPaths.append(fp)
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
            batches = [fitsPaths[i : i + batch_size] for i in range(0, len(fitsPaths), batch_size)]

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

        # # FIX BOUNDARY GROUP FILES -- MOVE TO NEXT DAY SO THEY GET COMBINED WITH THE REST OF THEIR GROUP
        # # E.G. UVB BIASES TAKEN ACROSS THE BOUNDARY BETWEEN 2 NIGHTS
        # # FIRST FIND END OF NIGHT DATA - AND PUSH TO THE NEXT DAY
        # mask = rawFrames["boundary"] > 0.96
        # filteredDf = rawFrames.loc[mask].copy()
        # filteredDf["night start mjd"] = filteredDf["night start mjd"] + 1
        # mask = filteredDf["eso dpr type"].isin(["LAMP,DFLAT", "LAMP,QFLAT"])
        # filteredDf.loc[mask, "eso dpr type"] = "LAMP,FLAT"
        # # NOW FIND START OF NIGHT DATA
        # mask = rawFrames["boundary"] < 0.04
        # filteredDf2 = rawFrames.loc[mask].copy()
        # mask = filteredDf2["eso dpr type"].isin(["LAMP,DFLAT", "LAMP,QFLAT"])
        # filteredDf2.loc[mask, "eso dpr type"] = "LAMP,FLAT"

        # # NOW FIND MATCHES BETWEEN 2 DATASETS
        # theseKeys = [
        #     "eso seq arm",
        #     "eso dpr catg",
        #     "eso dpr tech",
        #     "eso dpr type",
        #     "eso pro catg",
        #     "eso pro tech",
        #     "eso pro type",
        #     "night start mjd",
        # ]
        # matched = pd.merge(filteredDf, filteredDf2, on=theseKeys)
        # boundaryFiles = np.unique(matched["file_x"].values)
        # mask = rawFrames["file"].isin(boundaryFiles)
        # rawFrames.loc[mask, "night start mjd"] += 1

        # rawFrames["night start date"] = Time(rawFrames["night start mjd"], format="mjd").to_value("iso", subfmt="date")

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
        filteredFrames["gain"] = -99.99

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

        # CHECK GAIN AND CONAD ARE CORRECTLY POPULATED
        filteredFrames["gain"] = filteredFrames[self.kw("CONAD").lower()]
        mask = filteredFrames[self.kw("GAIN").lower()] > filteredFrames[self.kw("CONAD").lower()]
        filteredFrames.loc[mask, "gain"] = filteredFrames.loc[mask, self.kw("GAIN").lower()]

        # ADD SIMULATION FLAG FOR SPECTROSCOPIC DATA (AND MORE)
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
            filteredFrames = filteredFrames.rename(
                columns={
                    "eso ins temp217 val": "nir temp k",
                    "eso ins temp104 val": "vis temp c",
                    "eso ins temp301 val": "cp temp c",
                    "eso ins afc1 pos1": "afc1 pos1",
                    "eso ins afc1 pos2": "afc1 pos2",
                    "eso ins afc2 pos1": "afc2 pos1",
                    "eso ins afc2 pos2": "afc2 pos2",
                }
            )
        else:
            filteredFrames["simulation"] = 0
            filteredFrames["nir temp k"] = 0
            filteredFrames["vis temp c"] = 0
            filteredFrames["cp temp c"] = 0
            filteredFrames["afc1 pos1"] = 0
            filteredFrames["afc1 pos2"] = 0
            filteredFrames["afc2 pos1"] = 0
            filteredFrames["afc2 pos2"] = 0

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

        rawFrames["exptime"] = rawFrames["exptime"].apply(lambda x: round(x, 2))

        rawGroups = self._group_raw_frames(rawFrames, filterKeywordsRaw, addFilepaths=False, addStartDate=False)

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
            if (
                os.path.isfile(filepath)
                and os.path.splitext(filepath)[1] != ".db"
                and "readme." not in d.lower()
                and "soxspipe" not in d.lower()
            ):
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

        tries = 50

        while i < tries:
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
                    self.freshRun = True
                    emptyDb = os.path.dirname(os.path.dirname(__file__)) + "/resources/soxspipe.db"
                    shutil.copyfile(emptyDb, self.rootDbPath)
                conn = sql.connect(self.rootDbPath, timeout=300, autocommit=True, check_same_thread=False)
            c = conn.cursor()

            try:
                c.execute("PRAGMA integrity_check;")
                c.execute("PRAGMA busy_timeout = 100000")
                c.execute("PRAGMA synchronous = OFF")

                this = c.fetchall()
                i = tries + 1
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

                if i > tries - 1:
                    self.prepare(refresh=True)
                    if not reset:
                        reset = True
                conn = sql.connect(self.rootDbPath, timeout=300, autocommit=True, check_same_thread=False)
                c = conn.cursor()
        c.execute("PRAGMA busy_timeout = 100000")
        c.execute("PRAGMA synchronous = OFF")

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
            "update raw_frames set lamp = null, slit = null, slitmask = null where `eso dpr type` in ('BIAS','DARK');",
            "update raw_frames set rospeed = null where rospeed = -1;",
            """WITH s AS (
                    SELECT
                        uuid,
                        "file",
                        "eso tpl start",
                        "eso tpl expno" AS expno,
                        CASE
                            WHEN LAG("eso tpl expno") OVER (ORDER BY "eso tpl start", "eso tpl expno") IS NULL THEN 1
                            WHEN "eso tpl expno" <= LAG("eso tpl expno") OVER (ORDER BY "eso tpl start", "eso tpl expno") THEN 1
                            WHEN "eso tpl name" != LAG("eso tpl name") OVER (ORDER BY "eso tpl start", "eso tpl expno") THEN 1
                            ELSE 0
                        END AS new_set
                    FROM raw_frames
                ),
                g AS (
                    SELECT
                        uuid,
                        "file",
                        "eso tpl start", 
                        expno,
                        SUM(new_set) OVER (ORDER BY "eso tpl start", expno ROWS UNBOUNDED PRECEDING) AS grp
                    FROM s
                ),
                -- Combine labeling and sizing in one pass
                set_info AS (
                    SELECT
                        grp,
                        FIRST_VALUE("file") OVER (PARTITION BY grp ORDER BY "eso tpl start", expno) AS set_file,
                        COUNT(*) OVER (PARTITION BY grp) AS set_size
                    FROM g
                ),
                -- Create final mapping
                final_mapping AS (
                    SELECT DISTINCT
                        g.uuid,
                        si.set_file,
                        si.set_size
                    FROM g
                    JOIN set_info si ON g.grp = si.grp
                )
                UPDATE raw_frames
                SET
                    set_first_file = fm.set_file,
                    set_size = fm.set_size
                FROM final_mapping fm
                WHERE raw_frames.uuid = fm.uuid;
         """,
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

                            exists = " AND ".join(
                                [
                                    f"""EXISTS (
                                SELECT 1 FROM cal_{ct} 
                                WHERE cal_{ct}.sof = product_frames.sof 
                                AND (cal_{ct}.upstream_status = 'pass' OR cal_{ct}.upstream_status IS NULL)
                            )"""
                                    for ct in calType
                                ]
                            )

                            sqlQuery = f"""select sof from product_frames
                                WHERE complete < 1 
                                AND recipe = '{recipe}'
                                AND "eso seq arm" = '{arm}'
                                AND {exists};"""

                            containerSofs = pd.read_sql(sqlQuery, con=self.conn)["sof"].tolist()

                            self.raw_frames_to_sof_map(rawGroups=rawGroups, containerSofs=containerSofs)

                            sqlQuery = f"""UPDATE product_frames 
                                SET complete = -1 
                                WHERE complete < 1 
                                AND recipe = '{recipe}'
                                AND "eso seq arm" = '{arm}'
                                AND {exists};"""

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
            f"UPDATE sof_map_{self.sessionId} SET complete = 1 WHERE complete = -1;",
            "UPDATE raw_frame_sets SET complete = 1 WHERE sof IN (SELECT r.sof FROM raw_frame_sets r JOIN product_frames p ON p.sof = r.sof WHERE p.complete = 1);",
            "UPDATE raw_frame_sets SET complete = 0 WHERE sof IN (SELECT r.sof FROM raw_frame_sets r JOIN product_frames p ON p.sof = r.sof WHERE p.complete = 0);",
            "UPDATE raw_frames SET processed = 1 WHERE processed = 0 AND filepath IN (SELECT filepath FROM sof_map);",
            "UPDATE product_frames SET set_first_file = (SELECT s.file FROM sof_map s WHERE s.sof = product_frames.sof AND s.file LIKE 'SOXS%' LIMIT 1) WHERE set_first_file IS NULL;",
            "UPDATE raw_frame_sets SET set_first_file = (SELECT s.file FROM sof_map s WHERE s.sof = raw_frame_sets.sof AND s.file LIKE 'SOXS%' LIMIT 1) WHERE set_first_file IS NULL;",
            "UPDATE product_frames SET absrot=(SELECT r.absrot FROM raw_frames r WHERE r.file=product_frames.set_first_file);",
            "UPDATE raw_frame_sets SET absrot=(SELECT r.absrot FROM raw_frames r WHERE r.file=raw_frame_sets.set_first_file);",
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
        import numpy as np

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

        rawFrames = rawFrames.astype(
            {
                "exptime": float,
                "gain": float,
                "ra": float,
                "dec": float,
                "eso tel parang end": float,
                "eso tel parang start": float,
                "eso tel az": float,
                "eso tel alt": float,
                "eso tel ambi fwhm end": float,
                "eso tel ambi fwhm start": float,
                "eso tel airm end": float,
                "eso tel airm start": float,
                "absrot": float,
                "nir temp k": float,
                "vis temp c": float,
                "cp temp c": float,
                "afc1 pos1": float,
                "afc1 pos2": float,
                "afc2 pos1": float,
                "afc2 pos2": float,
            }
        )
        rawFrames.fillna(
            {
                "exptime": -99.99,
                "gain": -99.99,
                "ra": -99.99,
                "dec": -99.99,
                "eso tel parang end": -99.99,
                "eso tel parang start": -99.99,
                "eso tel az": -99.99,
                "eso tel alt": -99.99,
                "eso tel ambi fwhm end": -99.99,
                "eso tel ambi fwhm start": -99.99,
                "eso tel airm end": -99.99,
                "eso tel airm start": -99.99,
                "absrot": -99.99,
                "nir temp k": -99.99,
                "vis temp c": -99.99,
                "cp temp c": -99.99,
                "afc1 pos1": -99.99,
                "afc1 pos2": -99.99,
                "afc2 pos1": -99.99,
                "afc2 pos2": -99.99,
            },
            inplace=True,
        )
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

        if not len(rawFramesNoOffFrames.index):
            return pd.DataFrame(), pd.DataFrame()

        rawGroups = self._group_raw_frames(rawFramesNoOffFrames, filterKeywordsRaw + ["set_first_file"])

        # REMOVE GROUPED STARE - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = rawGroups["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"])
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD SCIENCE FRAMES AS ONE ENTRY PER EXPOSURE
        rawScienceFrames = rawFrames.loc[rawFrames["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"])]
        if len(rawScienceFrames.index):
            rawScienceFrames = self._group_raw_frames(rawScienceFrames, filterKeywordsRaw + ["mjd-obs"])
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
            rawPinholeFrames = self._group_raw_frames(rawPinholeFrames, filterKeywordsRaw + ["mjd-obs"])
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

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        mask = (rawGroups["recipe"].isin(("mbias", "mdark"))) & (rawGroups["counts"] < rawGroups["eso tpl nexp"])
        mask = mask | ((rawGroups["recipe"] == "mflat") & (rawGroups["counts"] < 5))
        rawGroups = rawGroups.loc[~mask]

        return rawFrames, rawGroups

    def _group_raw_frames(self, rawFrames, filterKeywordsRaw, addFilepaths=True, addStartDate=True):
        """Group raw frames and return grouped rows with aggregation metadata."""
        import pandas as pd

        # Create aggregation dictionary
        agg_dict = {col: "mean" for col in self.filterKeywordsExtras if col not in filterKeywordsRaw}
        agg_dict["file"] = "size"  # for counting rows

        rawGroups = rawFrames.groupby(filterKeywordsRaw)

        if addFilepaths:
            filepaths = [list(group["filepath"].values) for name, group in rawGroups]
        if addStartDate:
            startTime = rawGroups.min()["date-obs"].values

        # Group and aggregate
        rawGroups = rawFrames.groupby(filterKeywordsRaw).agg(agg_dict).rename(columns={"file": "counts"}).reset_index()
        rawGroups.style.hide(axis="index")
        pd.options.mode.chained_assignment = None

        # Normalise timestamps to compact YYYYMMDDTHHMMSS format for SOF naming.
        if addFilepaths:
            rawGroups["filepaths"] = filepaths
        if addStartDate:
            startTime = [str(s).split(".")[0].replace("-", "").replace(":", "") for s in startTime]
            rawGroups["date-obs"] = startTime

        return rawGroups

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
                + productFrames["night start date"].astype(str)
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
                dataframe.replace(["--", -99.99], None).to_sql(
                    table_name, con=self.conn, index=False, if_exists="append"
                )
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

    masterTable = ImageFileCollection(filenames=batch, keywords=keywords)
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
        # masterTable["utc-4hrs"] = (masterTable["mjd-obs"] - 2 / 3).astype(int)
        mjdDate = Time(masterTable["mjd-obs"], format="mjd", scale="utc") - mjd_ofset
        masterTable["mjd-date"] = mjdDate.strftime("%Y-%m-%d")
        masterTable["utc-4hrs"] = chileTimes.strftime("%Y-%m-%dt%H:%M:%S")
        masterTable["night start date"] = startNightDate.strftime("%Y-%m-%d")
        masterTable["night start mjd"] = startNightDate.mjd.astype(int)
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

    rawFrames = masterTable.to_pandas(index=False)

    # ADD FILEPATHS IF IN ./raw/ FOLDER
    rawFrames["filepath"] = "--"
    rawFrames["file"] = (
        rawFrames["file"].astype(str).str.replace(r"^.*?(raw/\d{4}-\d{2}-\d{2}.*)$", r"./\1", regex=True)
    )
    mask = rawFrames["file"].str.contains(r"\.\/raw\/\d{4}\-\d{2}\-\d{2}.*$", regex=True, na=False)
    rawFrames.loc[mask, "filepath"] = rawFrames.loc[mask, "file"]

    # MAKE FILE NAME ONLY THE BASENAME IF IN ./raw/ FOLDER
    rawFrames.loc[mask, "file"] = rawFrames.loc[mask, "file"].apply(lambda x: os.path.basename(x))

    return rawFrames
