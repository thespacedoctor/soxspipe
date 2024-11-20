#!/usr/bin/env python
# encoding: utf-8
"""
*The SOXSPIPE Data Organiser*

Author
: David Young

Date Created
: March  9, 2023
"""
from fundamentals import tools
from builtins import object
import sys
import os
from soxspipe.commonutils import uncompress
from soxspipe.commonutils.toolkit import get_calibrations_path
os.environ['TERM'] = 'vt100'


class data_organiser(object):
    """
    *The `soxspipe` Data Organiser*

    **Key Arguments:**

    - ``log`` -- logger
    - ``rootDir`` -- the root directory of the data to process

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

    def __init__(
            self,
            log,
            rootDir
    ):
        from os.path import expanduser
        import codecs
        from fundamentals.logs import emptyLogger
        import warnings
        from astropy.utils.exceptions import AstropyWarning
        warnings.simplefilter('ignore', AstropyWarning)

        self.PAE = True

        log.debug("instantiating a new 'data_organiser' object")

        self.log = log

        # if rootDir == ".":
        #     rootDir = os.getcwd()
        # MAKE RELATIVE HOME PATH ABSOLUTE
        if rootDir[0] == "~":
            home = expanduser("~")
            directory = directory.replace("~", home)

        self.rootDir = rootDir
        self.rawDir = rootDir + "/raw"
        self.miscDir = rootDir + "/misc"
        self.sessionsDir = rootDir + "/sessions"

        # SESSION ID PLACEHOLDER FILE
        self.sessionIdFile = self.sessionsDir + "/.sessionid"

        self.rootDbPath = rootDir + "/soxspipe.db"

        exists = os.path.exists(self.sessionIdFile)
        if exists:
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='r') as readFile:
                sessionId = readFile.read()
                self.sessionPath = self.sessionsDir + "/" + sessionId
                self.sessionId = sessionId

        self.keyword_lookups = [
            'MJDOBS',
            'DATE_OBS',
            'SEQ_ARM',
            'DPR_CATG',
            'DPR_TECH',
            'DPR_TYPE',
            'PRO_CATG',
            'PRO_TECH',
            'PRO_TYPE',
            'EXPTIME',
            'WIN_BINX',
            'WIN_BINY',
            'DET_READ_SPEED',
            'SLIT_UVB',
            'SLIT_VIS',
            'SLIT_NIR',
            'LAMP1',
            'LAMP2',
            'LAMP3',
            'LAMP4',
            'LAMP5',
            'LAMP6',
            'LAMP7',
            'DET_READ_TYPE',
            'CONAD',
            'RON',
            'OBS_ID',
            'OBS_NAME',
            "NAXIS",
            "OBJECT",
            "TPL_ID",
            "INSTRUME",
            "ABSROT"
        ]

        # THE MINIMUM SET OF KEYWORD WE EVER WANT RETURNED
        self.keywordsTerse = [
            'file',
            'eso seq arm',
            'eso dpr catg',
            'eso dpr type',
            'eso dpr tech',
            'eso pro catg',
            'eso pro tech',
            'eso pro type',
            'eso obs id',
            'eso obs name',
            'exptime',
            'binning',
            'rospeed',
            'slit',
            'slitmask',
            'lamp',
            'night start date',
            'night start mjd',
            'mjd-obs',
            'date-obs',
            'object',
            "template",
            "instrume",
            "absrot"
        ]

        # THIS TYPE MAP WILL BE USED TO GROUP SET OF FILES TOGETHER
        self.typeMapXSH = {
            "bias": [{"tech": None, "slitmask": None, "recipe": "mbias"}],  # XSH/SOXS BIAS CAN BE DEFINED WITH JUST DPR TYPE
            "dark": [{"tech": None, "slitmask": None, "recipe": "mdark"}],  # XSH/SOXS DARK CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,fmtchk": [{"tech": None, "slitmask": None, "recipe": "disp_sol"}],  # XSH disp_sol CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,orderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,dorderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,qorderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,flat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],  # XSH flats CAN BE DEFINED WITH JUST DPR TYPE
            "flat,lamp": [{"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"}],
            "lamp,dflat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],
            "lamp,qflat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],
            "lamp,wave": [{"tech": ["echelle,multi-pinhole", "image"], "slitmask": None, "recipe": "spat_sol"}, {"tech": ["echelle,pinhole", "image"], "slitmask": None, "recipe": "disp_sol"}],
            "wave,lamp": [{"tech": ["echelle,multi-pinhole", "image"], "slitmask": ["MPH"], "recipe": "spat_sol"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "disp_sol"}],
            "object": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std,flux": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std,telluric": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}]
        }

        # THIS TYPE MAP WILL BE USED TO GROUP SET OF FILES TOGETHER
        self.typeMapSOXS = {
            "bias": [{"tech": None, "slitmask": None, "recipe": "mbias"}],  # XSH/SOXS BIAS CAN BE DEFINED WITH JUST DPR TYPE
            "dark": [{"tech": None, "slitmask": None, "recipe": "mdark"}],  # XSH/SOXS DARK CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,fmtchk": [{"tech": None, "slitmask": None, "recipe": "disp_sol"}],  # XSH disp_sol CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,orderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,dorderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,qorderdef": [{"tech": None, "slitmask": None, "recipe": "order_centres"}],  # XSH order_centres CAN BE DEFINED WITH JUST DPR TYPE
            # "lamp,flat": [{"tech": None, "slitmask": None, "recipe": "mflat"}],  # XSH flats CAN BE DEFINED WITH JUST DPR TYPE
            "lamp,flat": [{"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"}],
            "flat,lamp": [{"tech": ["echelle,slit", "image"], "slitmask": ["SLIT"], "recipe": "mflat"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"}],
            "lamp,dflat": [{"tech": None, "slitmask": ["SLIT"], "recipe": "mflat"}, {"tech": None, "slitmask": ["PH"], "recipe": "order_centres"}],
            "lamp,qflat": [{"tech": None, "slitmask": ["SLIT"], "recipe": "mflat"}, {"tech": None, "slitmask": ["PH"], "recipe": "order_centres"}],
            "lamp,wave": [{"tech": ["echelle,multi-pinhole", "image"], "slitmask": None, "recipe": "spat_sol"}, {"tech": ["echelle,pinhole", "image"], "slitmask": None, "recipe": "disp_sol"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "order_centres"}],
            "wave,lamp": [{"tech": ["echelle,multi-pinhole", "image"], "slitmask": ["MPH"], "recipe": "spat_sol"}, {"tech": ["echelle,pinhole", "image"], "slitmask": ["PH"], "recipe": "disp_sol"}],
            "object": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "object,async": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std,flux": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}],
            "std,telluric": [{"tech": ["echelle,slit,stare"], "slitmask": None, "recipe": "stare"}, {"tech": ["echelle,slit,nodding"], "slitmask": None, "recipe": "nod"}, {"tech": ["echelle,slit,offset"], "slitmask": None, "recipe": "offset"}]
        }

        # THIS PRODUCT MAP IS USED TO PREDICT THE PRODUCTS THAT WILL RESULTS FROM REDUCING EACH SOFs
        # LIST ORDER: ['pro type kw', 'pro tech kw', 'pro catg kw', "pixels/table", "find in sof", "replace/suffix in sof", "recipe"]
        self.productMap = {
            "mbias": [
                ["REDUCED", "IMAGE", "MASTER_BIAS", "PIXELS", None, None, "soxs-mbias"]
            ],
            "mdark": [
                ["REDUCED", "IMAGE", "MASTER_DARK", "PIXELS", None, None, "soxs-mdark"]
            ],
            "disp_sol": [
                ["REDUCED", "ECHELLE,PINHOLE", "DISP_TAB", "TABLE", None, None, "soxs-disp-solution"]
            ],
            "order_centres": [
                ["REDUCED", "ECHELLE,SLIT", "ORDER_TAB", "TABLE", None, None, "soxs-order-centre"]
            ],
            "mflat": [
                ["REDUCED", "ECHELLE,SLIT", "MASTER_FLAT", "PIXELS", None, None, "soxs-mflat"],
                ["REDUCED", "ECHELLE,SLIT", "ORDER_TAB", "TABLE", "MFLAT", "ORDER_LOCATIONS", "soxs-mflat"]
            ],
            "spat_sol": [
                ["REDUCED", "ECHELLE,PINHOLE", "DISP_TAB", "TABLE", None, None, "soxs-spatial-solution"],
                ["REDUCED", "ECHELLE,PINHOLE", "DISP_IMAGE", "PIXELS", ".fits", "_IMAGE.fits", "soxs-spatial-solution"]
            ],
            "stare": [
                ["REDUCED", "ECHELLE,SLIT,STARE", "OBJECT_TAB", "TABLE", None, None, "soxs-stare"]
            ],
        }

        self.proKeywords = ['eso pro type', 'eso pro tech', 'eso pro catg']

        # THESE ARE KEYS WE NEED TO FILTER ON, AND SO NEED TO CREATE ASTROPY TABLE
        # INDEXES
        self.filterKeywords = ['eso seq arm', 'eso dpr catg',
                               'eso dpr tech', 'eso dpr type', 'eso pro catg', 'eso pro tech', 'eso pro type', 'exptime', 'rospeed', 'slit', 'slitmask', 'binning', 'night start mjd', 'night start date', 'instrume', "lamp", 'template', 'eso obs name']

        # THIS IS THE ORDER TO PROCESS THE FRAME TYPES
        self.reductionOrder = ["BIAS", "DARK", "LAMP,FMTCHK", "LAMP,ORDERDEF", "LAMP,DORDERDEF", "LAMP,QORDERDEF", "LAMP,FLAT", "FLAT,LAMP", "LAMP,DFLAT", "LAMP,QFLAT", "WAVE,LAMP", "LAMP,WAVE", "STD,FLUX", "STD", "STD,TELLURIC", "OBJECT", "OBJECT,ASYNC"]

        # THIS IS THE ORDER THE RECIPES NEED TO BE RUN IN (MAKE SURE THE REDUCTION SCRIPT HAS RECIPES IN THE CORRECT ORDER)
        self.recipeOrder = ["mbias", "mdark", "disp_sol", "order_centres", "mflat", "spat_sol", "stare", "nod", "offset"]

        # DECOMPRESS .Z FILES
        from soxspipe.commonutils import uncompress
        uncompress(
            log=self.log,
            directory=self.rootDir
        )

        return None

    def prepare(
            self):
        """*Prepare the workspace for data reduction by generating all SOF files and reduction scripts.*
        """
        self.log.debug('starting the ``prepare`` method')
        import codecs
        import shutil
        import sqlite3 as sql

        # TEST FITS FILES OR raw_frames DIRECT EXISTS
        fitsExist = False
        exists = os.path.exists(self.rawDir)
        if exists:
            from fundamentals.files import recursive_directory_listing
            theseFiles = recursive_directory_listing(
                log=self.log,
                baseFolderPath=self.rawDir,
                whatToList="files"  # all | files | dirs
            )
            for f in theseFiles:
                if os.path.splitext(f)[1] == ".fits" or ".fits.gz" in os.path.splitext(f):
                    fitsExist = True
                    break
        if not fitsExist:
            for d in os.listdir(self.rootDir):
                filepath = os.path.join(self.rootDir, d)
                if os.path.isfile(filepath) and (os.path.splitext(filepath)[1] == ".fits" or ".fits.gz" in os.path.splitext(filepath)):
                    fitsExist = True
                    break

        # EXIST IF NO FITS FILES EXIST - SOME PROTECT AGAINST MOVING USER FILES IF THEY MAKE A MISTAKE PREPARE A WORKSPACE IN THE WRONG LOCATION
        if fitsExist == False:
            print("There are no FITS files in this directory. Please add your data before running `soxspipe prep`")
            return None

        # MK RAW FRAME DIRECTORY
        if not os.path.exists(self.rawDir):
            os.makedirs(self.rawDir)

        # TEST FOR SQLITE DATABASE - ADD IF MISSING
        try:
            with open(self.rootDbPath):
                pass
            self.freshRun = False
        except IOError:
            try:
                os.remove(self.sessionIdFile)
            except:
                pass
            self.freshRun = True
            emptyDb = os.path.dirname(os.path.dirname(__file__)) + "/resources/soxspipe.db"
            shutil.copyfile(emptyDb, self.rootDbPath)

        def dict_factory(cursor, row):
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        # CREATE THE DATABASE CONNECTION
        self.conn = sql.connect(
            self.rootDbPath)
        # self.conn.row_factory = dict_factory

        # SELECT INSTR
        try:
            c = self.conn.cursor()
            sqlQuery = "select instrume from raw_frames where instrume is not null limit 1"
            c.execute(sqlQuery)
            self.instrument = c.fetchall()[0][0]
            c.close()
            if "SOXS" in self.instrument.upper():
                self.typeMap = self.typeMapSOXS
            else:
                self.typeMap = self.typeMapXSH
        except:
            pass

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
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='r') as readFile:
                sessionId = readFile.read()
                self.sessionPath = self.sessionsDir + "/" + sessionId

        i = 0
        while i < 6:
            self._populate_product_frames_db_table()
            i += 1
        self._write_sof_files()

        print(f"\nTHE `{basename}` WORKSPACE FOR HAS BEEN PREPARED FOR DATA-REDUCTION\n")
        print(f"In this workspace you will find:\n")
        print(f"   - `misc/`: a lost-and-found archive of non-fits files")
        print(f"   - `{self.rawDir}/`: all raw-frames to be reduced")
        print(f"   - `sessions/`: directory of data-reduction sessions")
        print(f"   - `sof/`: the set-of-files (sof) files required for each reduction step")
        print(f"   - `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete\n")

        self.conn.close()

        self.log.debug('completed the ``prepare`` method')
        return None

    def _sync_raw_frames(
            self,
            skipSqlSync=False):
        """*sync the raw frames between the project folder and the database *

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
        self.log.debug('starting the ``_sync_raw_frames`` method')

        import sqlite3 as sql
        import shutil
        import pandas as pd

        # GENERATE AN ASTROPY TABLES OF FITS FRAMES WITH ALL INDEXES NEEDED
        filteredFrames, fitsPaths, fitsNames = self._create_directory_table(pathToDirectory=self.rootDir, filterKeys=self.filterKeywords)

        if fitsPaths:

            conn = self.conn
            knownRawFrames = pd.read_sql('SELECT * FROM raw_frames', con=conn)

            # SPLIT INTO RAW, REDUCED PIXELS, REDUCED TABLES
            rawFrames, reducedFramesPixels, reducedFramesTables = self._categorise_frames(filteredFrames)

            # FILTER DATA FRAME
            if self.PAE and self.instrument.upper() == "SOXS":
                mask = ((rawFrames["eso dpr tech"] == "ECHELLE,PINHOLE") & (rawFrames["eso dpr type"] == "LAMP,FLAT"))
                filteredDf = rawFrames.loc[mask]
                filteredDf["eso dpr tech"] = "ECHELLE,SLIT,STARE"
                filteredDf["eso dpr type"] = "OBJECT"

                rawFrames = pd.concat([rawFrames, filteredDf], ignore_index=True)

            # xpd-update-filter-dataframe-column-values

            # FIND AND REMOVE DUPLICATE FILES
            if len(rawFrames.index):
                rawFrames["filepath"] = f"{self.rawDir}/" + rawFrames['night start date'] + "/" + rawFrames['file']
                # FIND AND REMOVE DUPLICATE FILES
                matchedFiles = pd.merge(rawFrames, knownRawFrames, on=['file', 'eso dpr tech'], how='inner')
                if len(matchedFiles.index):
                    for file in matchedFiles["file"]:
                        try:
                            os.remove(file)
                        except:
                            pass
                    # FIND RECORDS IN rawFrames NOT IN rawFrames
                    rawFrames = rawFrames[~rawFrames.set_index(['file', 'eso dpr tech']).index.isin(knownRawFrames.set_index(['file', 'eso dpr tech']).index)]

            if len(rawFrames.index):
                rawFrames["filepath"] = f"{self.rawDir}/" + rawFrames['night start date'] + "/" + rawFrames['file']

                rawFrames.to_sql('raw_frames', con=self.conn,
                                 index=False, if_exists='append')
                filepaths = rawFrames['filepath']
                filenames = rawFrames['file']
                for p, n in zip(filepaths, filenames):
                    parentDirectory = os.path.dirname(p)
                    if not os.path.exists(parentDirectory):
                        # Recursively create missing directories
                        if not os.path.exists(parentDirectory):
                            os.makedirs(parentDirectory)
                    if os.path.exists(self.rootDir + "/" + n):
                        shutil.move(self.rootDir + "/" + n, p)

        if not skipSqlSync:
            self._sync_sql_table_to_directory(self.rawDir, 'raw_frames', recursive=False)

        self.log.debug('completed the ``_sync_raw_frames`` method')
        return None

    def _create_directory_table(
            self,
            pathToDirectory,
            filterKeys):
        """*create an astropy table based on the contents of a directory*

        **Key Arguments:**

        - `log` -- logger
        - `pathToDirectory` -- path to the directory containing the FITS frames
        - `filterKeys` -- these are the keywords we want to filter on later

        **Return**

        - `masterTable` -- the primary dataframe table listing all FITS files in the directory (including indexes on `filterKeys` columns)
        - `fitsPaths` -- a simple list of all FITS file paths
        - `fitsNames` -- a simple list of all FITS file name

        **Usage:**

        ```python
        # GENERATE AN ASTROPY TABLES OF FITS FRAMES WITH ALL INDEXES NEEDED
        masterTable, fitsPaths, fitsNames = _create_directory_table(
            log=log,
            pathToDirectory="/my/directory/path",
            keys=["file","mjd-obs", "exptime","cdelt1", "cdelt2"],
            filterKeys=["mjd-obs","exptime"]
        )
        ```
        """
        self.log.debug('starting the ``_create_directory_table`` function')

        from ccdproc import ImageFileCollection
        from astropy.time import Time, TimeDelta
        import numpy as np
        from fundamentals.files import recursive_directory_listing
        import pandas as pd
        from soxspipe.commonutils import keyword_lookup

        # GENERATE A LIST OF FITS FILE PATHS
        fitsPaths = []
        fitsNames = []
        for d in os.listdir(pathToDirectory):
            filepath = os.path.join(pathToDirectory, d)
            if d[0] != "." and os.path.isfile(filepath) and (os.path.splitext(filepath)[1] == ".fits" or ".fits.Z" in filepath):
                fitsPaths.append(filepath)
                fitsNames.append(d)

        recursive = False
        # for d in os.listdir(pathToDirectory):
        #     if os.path.isdir(os.path.join(pathToDirectory, d)) and d in ('raw_frames', 'product'):
        #         recursive = True
        #         theseFiles = recursive_directory_listing(
        #             log=self.log,
        #             baseFolderPath=os.path.join(pathToDirectory, d),
        #             whatToList="files"  # all | files | dirs
        #         )
        #         newFitsPaths = [n for n in theseFiles if ".fits" in n]
        #         newFitsNames = [os.path.basename(n) for n in theseFiles if ".fits" in n]
        #         fitsPaths += newFitsPaths
        #         fitsNames += newFitsNames

        if len(fitsPaths) == 0:
            # print(f"No fits files found in directory `{pathToDirectory}`")
            return None, None, None

        # INSTRUMENT CHECK
        if recursive:
            allFrames = ImageFileCollection(filenames=fitsPaths, keywords=["instrume"])
        else:
            allFrames = ImageFileCollection(
                location=pathToDirectory, filenames=fitsNames, keywords=["instrume"])

        tmpTable = allFrames.summary
        tmpTable['instrume'].fill_value = "--"
        instrument = tmpTable['instrume'].filled()
        instrument = list(set(instrument))
        if "--" in instrument:
            instrument.remove("--")

        self.instrument = None
        if len(instrument) > 1:
            self.log.error(f'The directory contains data from a mix of instruments. Please only provide data from either SOXS or XSH')
            raise AssertionError
        else:
            self.instrument = instrument[0]

        if "XSH" in self.instrument.upper() or "SHOOT" in self.instrument.upper():
            self.instrument = "XSH"
        print(f"The instrument has been set to '{self.instrument}'")

        if "SOXS" in self.instrument.upper():
            self.typeMap = self.typeMapSOXS
        else:
            self.typeMap = self.typeMapXSH

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            instrument=self.instrument
        ).get
        self.keywords = ['file']
        for k in self.keyword_lookups:
            self.keywords.append(self.kw(k).lower())

        # TOP-LEVEL COLLECTIONi
        if recursive:
            allFrames = ImageFileCollection(filenames=fitsPaths, keywords=self.keywords)
        else:
            allFrames = ImageFileCollection(
                location=pathToDirectory, filenames=fitsNames, keywords=self.keywords)

        masterTable = allFrames.summary

        # ADD FILLED VALUES FOR MISSING CELLS
        for fil in self.keywords:
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

        # FILTER OUT FRAMES WITH NO MJD
        matches = ((masterTable["mjd-obs"] == -99.99) | (masterTable["eso dpr catg"] == "--") | (masterTable["eso dpr tech"] == "--") | (masterTable["eso dpr type"] == "--") | (masterTable["exptime"] == -99.99))
        missingMJDFiles = masterTable['file'][matches]
        if len(missingMJDFiles):
            print("\nThe following FITS files are missing DPR keywords and will be ignored:\n\n")
            print(missingMJDFiles)
            masterTable = masterTable[~matches]

        # SETUP A NEW COLUMN GIVING THE INT MJD THE CHILEAN NIGHT BEGAN ON
        # 12:00 NOON IN CHILE IS TYPICALLY AT 16:00 UTC (CHILE = UTC - 4)
        # SO COUNT CHILEAN OBSERVING NIGHTS AS 15:00 UTC-15:00 UTC (11am-11am)
        if "mjd-obs" in masterTable.colnames:
            chile_offset = TimeDelta(4.0 * 60 * 60, format='sec')
            night_start_offset = TimeDelta(15.0 * 60 * 60, format='sec')
            masterTable["mjd-obs"] = masterTable["mjd-obs"].astype(float)
            chileTimes = Time(masterTable["mjd-obs"],
                              format='mjd', scale='utc') - chile_offset
            startNightDate = Time(masterTable["mjd-obs"],
                                  format='mjd', scale='utc') - night_start_offset
            # masterTable["utc-4hrs"] = (masterTable["mjd-obs"] - 2 / 3).astype(int)
            masterTable["utc-4hrs"] = chileTimes.strftime("%Y-%m-%dt%H:%M:%S")
            masterTable["night start date"] = startNightDate.strftime("%Y-%m-%d")
            masterTable["night start mjd"] = startNightDate.mjd.astype(int)
            masterTable["boundary"] = startNightDate.mjd - startNightDate.mjd.astype(int)
            masterTable.add_index("night start date")
            masterTable.add_index("night start mjd")

        if self.instrument.upper() != "SOXS":
            if self.kw("DET_READ_SPEED").lower() in masterTable.colnames:
                masterTable["rospeed"] = np.copy(masterTable[self.kw("DET_READ_SPEED").lower()])
                try:
                    masterTable["rospeed"][masterTable[
                        "rospeed"] == -99.99] = '--'
                except:
                    masterTable["rospeed"] = masterTable["rospeed"].astype(str)
                    masterTable["rospeed"][masterTable[
                        "rospeed"] == -99.99] = '--'
                masterTable["rospeed"][masterTable[
                    "rospeed"] == '1pt/400k/lg'] = 'fast'
                masterTable["rospeed"][masterTable[
                    "rospeed"] == '1pt/400k/lg/AFC'] = 'fast'
                masterTable["rospeed"][masterTable[
                    "rospeed"] == '1pt/100k/hg'] = 'slow'
                masterTable["rospeed"][masterTable[
                    "rospeed"] == '1pt/100k/hg/AFC'] = 'slow'
                masterTable.add_index("rospeed")
        else:
            if self.kw("DET_READ_SPEED").lower() in masterTable.colnames:
                masterTable["rospeed"] = np.copy(masterTable[self.kw("DET_READ_SPEED").lower()])

                try:
                    masterTable["rospeed"][masterTable[
                        "rospeed"] == -99.99] = -1
                except:
                    masterTable["rospeed"] = masterTable["rospeed"].astype(str)
                    masterTable["rospeed"][masterTable[
                        "rospeed"] == -99.99] = -1

                masterTable.add_index("rospeed")

        if self.kw("TPL_ID").lower() in masterTable.colnames:
            masterTable["template"] = np.copy(masterTable[self.kw("TPL_ID").lower()])

        if "naxis" in masterTable.colnames:
            masterTable["table"] = np.copy(masterTable["naxis"]).astype(str)
            masterTable["table"][masterTable[
                "table"] == '0'] = 'T'
            masterTable["table"][masterTable[
                "table"] != 'T'] = 'F'

        if self.kw("WIN_BINX").lower() in masterTable.colnames:
            masterTable["binning"] = np.core.defchararray.add(
                masterTable[self.kw("WIN_BINX").lower()].astype('int').astype('str'), "x")
            masterTable["binning"] = np.core.defchararray.add(masterTable["binning"],
                                                              masterTable[self.kw("WIN_BINY").lower()].astype('int').astype('str'))
            masterTable["binning"][masterTable[
                "binning"] == '-99x-99'] = '--'
            masterTable["binning"][masterTable[
                "binning"] == '1x-99'] = '--'
            masterTable.add_index("binning")

        if self.kw("ABSROT").lower() in masterTable.colnames:
            masterTable["absrot"] = masterTable[self.kw("ABSROT").lower()].astype(float)
            masterTable.add_index("absrot")

        # ADD INDEXES ON ALL KEYS
        for k in self.keywords:
            try:
                masterTable.add_index(k)
            except:
                pass

        # SORT IMAGE COLLECTION
        masterTable.sort(['eso pro type', 'eso seq arm', 'eso dpr catg', 'eso dpr tech', 'eso dpr type', 'eso pro catg', 'eso pro tech', 'mjd-obs'])

        # FIX BOUNDARY GROUP FILES -- MOVE TO NEXT DAY SO THEY GET COMBINED WITH THE REST OF THEIR GROUP
        # E.G. UVB BIASES TAKEN ACROSS THE BOUNDARY BETWEEN 2 NIGHTS
        masterTable = masterTable.to_pandas(index=False)
        # FIRST FIND END OF NIGHT DATA - AND PUSH TO THE NEXT DAY
        mask = (masterTable['boundary'] > 0.96)
        filteredDf = masterTable.loc[mask].copy()
        filteredDf['night start mjd'] = filteredDf['night start mjd'] + 1
        mask = (filteredDf["eso dpr type"].isin(['LAMP,DFLAT', 'LAMP,QFLAT']))
        filteredDf.loc[mask, "eso dpr type"] = 'LAMP,FLAT'
        # NOW FIND START OF NIGHT DATA
        mask = (masterTable['boundary'] < 0.04)
        filteredDf2 = masterTable.loc[mask].copy()
        mask = (filteredDf2["eso dpr type"].isin(['LAMP,DFLAT', 'LAMP,QFLAT']))
        filteredDf2.loc[mask, "eso dpr type"] = 'LAMP,FLAT'

        # NOW FIND MATCHES BETWEEN 2 DATASETS
        theseKeys = ['eso seq arm', 'eso dpr catg', 'eso dpr tech', 'eso dpr type', 'eso pro catg', 'eso pro tech', 'eso pro type', 'night start mjd']
        matched = pd.merge(filteredDf, filteredDf2, on=theseKeys)
        boundaryFiles = np.unique(matched['file_x'].values)
        mask = (masterTable['file'].isin(boundaryFiles))
        masterTable.loc[mask, 'night start mjd'] += 1

        masterTable['night start date'] = Time(masterTable['night start mjd'], format='mjd').to_value('iso', subfmt='date')

        self.log.debug('completed the ``_create_directory_table`` function')
        return masterTable, fitsPaths, fitsNames

    def _sync_sql_table_to_directory(
            self,
            directory,
            tableName,
            recursive=False):
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
        self.log.debug('starting the ``_sync_sql_table_to_directory`` method')

        import sqlite3 as sql
        import shutil
        import time

        # GENERATE A LIST OF FITS FILE PATHS IN RAW DIR
        from fundamentals.files import recursive_directory_listing
        fitsPaths = recursive_directory_listing(
            log=self.log,
            baseFolderPath=directory,
            whatToList="files"  # all | files | dirs
        )

        c = self.conn.cursor()

        sqlQuery = f"select filepath from {tableName};"
        c.execute(sqlQuery)

        dbFiles = [r[0] for r in c.fetchall()]

        # DELETED FILES
        filesNotInDB = set(fitsPaths) - set(dbFiles)
        filesNotInFS = set(dbFiles) - set(fitsPaths)

        if len(filesNotInFS):
            filesNotInFS = ("','").join(filesNotInFS)
            sqlQuery = f"delete from {tableName} where filepath in ('{filesNotInFS}');"
            c.execute(sqlQuery)

        if len(filesNotInDB):
            for f in filesNotInDB:
                # GET THE EXTENSION (WITH DOT PREFIX)
                basename = os.path.basename(f)
                extension = os.path.splitext(basename)[1]
                if extension.lower() != ".fits":
                    pass
                elif self.rootDir in f:
                    shutil.move(f, self.rootDir)
                else:
                    shutil.move(self.rootDir + "/" + f, self.rootDir)
            self._sync_raw_frames(skipSqlSync=True)

        c.close()

        self.log.debug('completed the ``_sync_sql_table_to_directory`` method')
        return None

    def _categorise_frames(
            self,
            filteredFrames,
            verbose=False):
        """*given a dataframe of frame, categorise frames into raw, reduced pixels, reduced tables*

        **Key Arguments:**

        - ``filteredFrames`` -- the dataframe from which to split frames into categorise.
        - ``verbose`` -- print results to stdout.

        **Return:**

        - ``rawFrames`` -- dataframe of raw frames only
        - ``reducedFramesPixels`` -- dataframe of reduced images only
        - ``reducedFramesTables`` -- dataframe of reduced tables only

        **Usage:**

        ```python
        rawFrames, reducedFramesPixels, reducedFramesTables = self._categorise_frames(filteredFrames)
        ```
        """
        self.log.debug('starting the ``catagorise_frames`` method')

        from astropy.table import Table, unique
        import numpy as np
        import pandas as pd
        from tabulate import tabulate

        # SPLIT INTO RAW, REDUCED PIXELS, REDUCED TABLES
        keywordsTerseRaw = self.keywordsTerse[:]
        keywordsTerseReduced = self.keywordsTerse[:]
        filterKeywordsRaw = self.filterKeywords[:]
        filterKeywordsReduced = self.filterKeywords[:]

        filteredFrames['slit'] = "--"
        filteredFrames['slitmask'] = "--"
        filteredFrames['lamp'] = "--"

        # ADD SLIT FOR SPECTROSCOPIC DATA
        filteredFrames.loc[(filteredFrames['eso seq arm'] == "NIR"), "slit"] = filteredFrames.loc[(filteredFrames['eso seq arm'] == "NIR"), self.kw("SLIT_NIR").lower()]
        filteredFrames.loc[(filteredFrames['eso seq arm'] == "VIS"), "slit"] = filteredFrames.loc[(filteredFrames['eso seq arm'] == "VIS"), self.kw("SLIT_VIS").lower()]
        filteredFrames.loc[(filteredFrames['eso seq arm'] == "UVB"), "slit"] = filteredFrames.loc[(filteredFrames['eso seq arm'] == "UVB"), self.kw("SLIT_UVB").lower()]

        filteredFrames["slit"] = filteredFrames["slit"].str.upper()

        filteredFrames.loc[((filteredFrames['slit'].str.contains("MULT")) & (filteredFrames['slitmask'] == "--")), "slitmask"] = "MPH"
        filteredFrames.loc[((filteredFrames['slit'].str.contains("PINHOLE")) & (filteredFrames['slitmask'] == "--")), "slitmask"] = "PH"
        filteredFrames.loc[((filteredFrames['slit'].str.contains("SLIT")) & (filteredFrames['slitmask'] == "--")), "slitmask"] = "SLIT"

        lampLong = ["argo", "merc", "neon", "xeno", "qth", "deut", "thar"]
        lampEle = ["Ar", "Hg", "Ne", "Xe", "QTH", "D", "ThAr"]

        for i in [1, 2, 3, 4, 5, 6, 7]:
            lamp = self.kw(f"LAMP{i}").lower()

            if self.instrument.lower() == "soxs":
                for l, e in zip(lampLong, lampEle):
                    if l in lamp:
                        lamp = e
                filteredFrames.loc[((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] != "--")), "lamp"] += lamp
                filteredFrames.loc[((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] == "--")), "lamp"] = lamp
            else:
                filteredFrames.loc[((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99)), "lamp"] = filteredFrames.loc[((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99)), self.kw(f"LAMP{i}").lower()]
        mask = []
        for i in self.proKeywords:
            keywordsTerseRaw.remove(i)
            filterKeywordsRaw.remove(i)
            if not len(mask):
                mask = (filteredFrames[i] == "--")
            else:
                mask = np.logical_and(mask, (filteredFrames[i] == "--"))

        filteredFrames["lamp"] = filteredFrames["lamp"].str.replace("_lamp", "")
        filteredFrames["lamp"] = filteredFrames["lamp"].str.replace("_Lamp", "")

        rawFrames = filteredFrames.loc[mask]

        # MATCH OFF FRAMES TO ADD THE MISSING LAMPS
        mask = (rawFrames["eso obs name"] == "Maintenance")
        rawFrames.loc[mask, "eso obs name"] = rawFrames.loc[mask, "eso obs name"] + rawFrames.loc[mask, "eso dpr type"]
        if self.instrument.lower() == "soxs":
            groupBy = 'eso obs name'
        else:
            groupBy = 'template'
        rawFrames.loc[(rawFrames['lamp'] == "--"), 'lamp'] = np.nan
        rawFrames.loc[(rawFrames['eso seq arm'].str.lower() == "nir"), "lamp"] = rawFrames.loc[(rawFrames['eso seq arm'].str.lower() == "nir"), 'lamp'].fillna(rawFrames.loc[(rawFrames['eso seq arm'].str.lower() == "nir")].groupby(groupBy)['lamp'].transform('first'))
        rawFrames.loc[(rawFrames['lamp'].isnull()), 'lamp'] = "--"
        rawFrames.loc[(rawFrames['absrot'] == -99.99), 'absrot'] = 0.0

        rawFrames["exptime"] = rawFrames["exptime"].apply(lambda x: round(x, 2))

        reducedFrames = filteredFrames.loc[~mask]
        pd.options.display.float_format = '{:,.4f}'.format

        mask = (reducedFrames["naxis"] == 0)
        reducedFramesTables = reducedFrames.loc[mask]
        reducedFramesPixels = reducedFrames.loc[~mask]

        rawGroups = rawFrames.groupby(filterKeywordsRaw).size().reset_index(name='counts')
        rawGroups.style.hide(axis='index')
        pd.options.mode.chained_assignment = None

        dprKeywords = ['eso dpr type', 'eso dpr tech', 'eso dpr catg']
        for i in dprKeywords:
            keywordsTerseReduced.remove(i)
            filterKeywordsReduced.remove(i)
        filterKeywordsReducedTable = filterKeywordsReduced[:]
        filterKeywordsReducedTable.remove("binning")
        keywordsTerseReducedTable = keywordsTerseReduced[:]
        keywordsTerseReducedTable.remove("binning")

        reducedPixelsGroups = reducedFramesPixels.groupby(filterKeywordsReduced).size().reset_index(name='counts')
        reducedPixelsGroups.style.hide(axis='index')

        # SORT BY COLUMN NAME
        reducedPixelsGroups.sort_values(by=['eso pro type', 'eso seq arm', 'eso pro catg', 'eso pro tech'], inplace=True)

        reducedTablesGroups = reducedFramesTables.groupby(filterKeywordsReducedTable).size().reset_index(name='counts')
        reducedTablesGroups.style.hide(axis='index')
        reducedTablesGroups.sort_values(by=['eso pro type', 'eso seq arm', 'eso pro catg', 'eso pro tech'], inplace=True)

        if verbose and (len(reducedPixelsGroups.index) or len(reducedTablesGroups.index)):
            print("\n# CONTENT SETS INDEX\n")

        if verbose and len(reducedPixelsGroups.index):
            print("\n## REDUCED PIXEL-FRAME-SET SUMMARY\n")
            print(tabulate(reducedPixelsGroups, headers='keys', tablefmt='github', showindex=False, stralign="right"))

        if verbose and len(reducedTablesGroups.index):
            print("\n## REDUCED TABLES-SET SUMMARY\n")
            print(tabulate(reducedTablesGroups, headers='keys', tablefmt='github', showindex=False, stralign="right"))

        if verbose:
            print("\n# CONTENT FILE INDEX\n")
        if verbose and len(rawGroups.index):
            print("\n## ALL RAW FRAMES\n")
            print(tabulate(rawFrames[keywordsTerseRaw], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

        if verbose and len(reducedPixelsGroups.index):
            print("\n## ALL REDUCED PIXEL-FRAMES\n")
            print(tabulate(reducedFramesPixels[keywordsTerseReduced], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

        if verbose and len(reducedTablesGroups.index):
            print("\n## ALL REDUCED TABLES\n")
            print(tabulate(reducedFramesTables[keywordsTerseReducedTable], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

        self.log.debug('completed the ``catagorise_frames`` method')
        return rawFrames[keywordsTerseRaw].replace(['--'], None), reducedFramesPixels[keywordsTerseReduced], reducedFramesTables[keywordsTerseReducedTable]

    def _populate_product_frames_db_table(
            self):
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
        self.log.debug('starting the ``_populate_product_frames_db_table`` method')

        import pandas as pd
        import sqlite3 as sql
        import time

        conn = self.conn
        rawFrames = pd.read_sql('SELECT * FROM raw_frames', con=conn)

        rawFrames.fillna("--", inplace=True)
        filterKeywordsRaw = self.filterKeywords[:]

        for i in self.proKeywords:
            filterKeywordsRaw.remove(i)

        # rawFrames.replace("LAMP,DFLAT", "LAMP,FLAT", inplace=True)
        # rawFrames.replace("LAMP,QFLAT", "LAMP,FLAT", inplace=True)

        # HIDE OFF FRAMES FROM GROUPS
        mask = ((rawFrames["eso dpr tech"] == "IMAGE") & (rawFrames['eso seq arm'] == "NIR") & (rawFrames['eso dpr type'] != "DARK"))
        rawFramesNoOffFrames = rawFrames.loc[~mask]
        rawGroups = rawFramesNoOffFrames.groupby(filterKeywordsRaw)
        mjds = rawGroups.mean(numeric_only=True)["mjd-obs"].values

        rawGroups = rawGroups.size().reset_index(name='counts')
        rawGroups['mjd-obs'] = mjds

        # REMOVE GROUPED STARE - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = (rawGroups["eso dpr tech"].isin(["ECHELLE,SLIT,STARE"]))
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD SCIENCE FRAMES AS ONE ENTRY PER EXPOSURE
        rawScienceFrames = pd.read_sql(
            'SELECT * FROM raw_frames where "eso dpr tech" in ("ECHELLE,SLIT,STARE")', con=conn)

        rawScienceFrames.fillna("--", inplace=True)
        rawScienceFrames = rawScienceFrames.groupby(filterKeywordsRaw + ["mjd-obs"])
        rawScienceFrames = rawScienceFrames.size().reset_index(name='counts')

        # MERGE DATAFRAMES
        if len(rawScienceFrames.index):
            rawGroups = pd.concat([rawGroups, rawScienceFrames], ignore_index=True)

        # REMOVE GROUPED SINGLE PINHOLE ARCS - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = (rawGroups["eso dpr tech"].isin(["ECHELLE,PINHOLE", "ECHELLE,MULTI-PINHOLE"]))
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD PINHOLE FRAMES AS ONE ENTRY PER EXPOSURE
        if self.instrument.upper() == "SOXS":
            rawPinholeFrames = pd.read_sql(
                'SELECT * FROM raw_frames where "eso dpr tech" in ("ECHELLE,PINHOLE","ECHELLE,MULTI-PINHOLE") and ("eso seq arm" = "NIR" or ("lamp" not in ("Xe", "Ar", "Hg", "Ne", "ArNeHgXe" )))', con=conn)
        else:
            rawPinholeFrames = pd.read_sql(
                'SELECT * FROM raw_frames where "eso dpr tech" in ("ECHELLE,PINHOLE","ECHELLE,MULTI-PINHOLE")', con=conn)
        rawPinholeFrames.fillna("--", inplace=True)
        rawPinholeFrames = rawPinholeFrames.groupby(filterKeywordsRaw + ["mjd-obs"])
        rawPinholeFrames = rawPinholeFrames.size().reset_index(name='counts')

        # MERGE DATAFRAMES
        if len(rawPinholeFrames.index):
            rawGroups = pd.concat([rawGroups, rawPinholeFrames], ignore_index=True)

        rawGroups['recipe'] = None
        rawGroups['sof'] = None

        calibrationFrames = pd.read_sql(f'SELECT * FROM product_frames where `eso pro catg` not like "%_TAB_%" and (status_{self.sessionId} != "fail" or status_{self.sessionId} is null)', con=conn)
        calibrationFrames.fillna("--", inplace=True)

        calibrationTables = pd.read_sql(f'SELECT * FROM product_frames where `eso pro catg` like "%_TAB_%" and (status_{self.sessionId} != "fail" or status_{self.sessionId} is null)', con=conn)
        calibrationTables.fillna("--", inplace=True)

        # _generate_sof_and_product_names SHOULD TAKE ROW OF DF AS INPUT

        # SEND TO DATABASE
        c = self.conn.cursor()
        sqlQuery = f"delete from sof_map_{self.sessionId};"
        c.execute(sqlQuery)
        c.close()

        repeat = 7
        while repeat:
            for o in self.reductionOrder:
                rawGroups = rawGroups.apply(self._generate_sof_and_product_names, axis=1, reductionOrder=o, rawFrames=rawFrames, calibrationFrames=calibrationFrames, calibrationTables=calibrationTables)
                rawGroups = rawGroups.apply(self._populate_products_table, axis=1, reductionOrder=o)
            repeat -= 1

        # xpd-update-filter-dataframe-column-values

        # SEND TO DATABASE
        c = self.conn.cursor()
        sqlQuery = f"delete from raw_frame_sets;"
        c.execute(sqlQuery)
        c.close()

        keepTrying = 0
        while keepTrying < 6:
            try:
                rawGroups.replace(['--'], None).to_sql('raw_frame_sets', con=self.conn,
                                                       index=False, if_exists='append')
                keepTrying = 10
            except Exception as e:
                if keepTrying > 5:
                    raise Exception(e)
                time.sleep(1)
                keepTrying += 1

        self.log.debug('completed the ``_populate_product_frames_db_table`` method')
        return None

    def _generate_sof_and_product_names(
            self,
            series,
            reductionOrder,
            rawFrames,
            calibrationFrames,
            calibrationTables):
        """*add a recipe name and SOF filename to all rows in the raw_frame_sets DB table*

        **Key Arguments:**

        - ``series`` -- the dataframe row/series to apply work on
        """

        import pandas as pd
        import astropy
        import numpy as np
        import time

        incomplete = False

        sofName = []
        matchDict = {}
        sofName.append(series['eso seq arm'].upper())
        matchDict['eso seq arm'] = series['eso seq arm'].upper()
        filteredFrames = rawFrames.copy()

        if series["eso dpr type"].lower() != reductionOrder.lower():
            return series

        # FILTER BY TEMPLATE NAME
        mask = (filteredFrames["template"].isin([series["template"]]))
        filteredFrames = filteredFrames.loc[mask]

        # FILTER BY TYPE FIRST
        if "FLAT" in series["eso dpr type"].upper():
            mask = ((filteredFrames["eso dpr type"].str.contains("FLAT")) & (filteredFrames["slit"] == series["slit"]))
        else:
            mask = (filteredFrames["eso dpr type"].isin([series["eso dpr type"].upper()]))
        filteredFrames = filteredFrames.loc[mask]

        seriesRecipe = None

        # CHECK SLIT
        if self.typeMap[series["eso dpr type"].lower()][0]["slitmask"]:
            match = False
            for row in self.typeMap[series["eso dpr type"].lower()]:
                rowSlit = [item.upper() for item in row["slitmask"]]
                if not match and series["slitmask"] in rowSlit:
                    match = True

                    mask = (filteredFrames["slitmask"].isin(rowSlit))
                    filteredFrames = filteredFrames.loc[mask]
                    seriesRecipe = row["recipe"]
            if not match:
                return series

        # CHECK TECH
        if self.typeMap[series["eso dpr type"].lower()][0]["tech"] and not seriesRecipe:
            match = False
            for row in self.typeMap[series["eso dpr type"].lower()]:
                rowTech = [item.upper() for item in row["tech"]]
                if not match and series["eso dpr tech"].upper() in rowTech:
                    match = True
                    mask = (filteredFrames["eso dpr tech"].isin(rowTech))
                    filteredFrames = filteredFrames.loc[mask]
                    seriesRecipe = row["recipe"]
            if not match:
                return series
        elif not seriesRecipe:
            seriesRecipe = self.typeMap[series["eso dpr type"].lower()][0]["recipe"]

        # GENEREATE SOF FILENAME AND MATCH DICTIONARY TO FILTER ON
        if series["binning"] != "--":
            matchDict['binning'] = series["binning"]
            sofName.append(series['binning'].upper())
        if series["rospeed"] != "--":
            matchDict['rospeed'] = series["rospeed"]
            sofName.append(series["rospeed"])
        if series["eso dpr type"].lower() in self.typeMap:
            matchDict['eso dpr type'] = series["eso dpr type"]
            for i in self.typeMap[series["eso dpr type"].lower()]:
                if i["recipe"] == seriesRecipe:
                    sofName.append(i["recipe"].replace("_centres", "_locations"))

            if "DORDER" in series["eso dpr type"].upper():
                sofName.append("dlamp")
            if "QORDER" in series["eso dpr type"].upper():
                sofName.append("qlamp")
        if True and ("PINHOLE" in series["eso dpr tech"].upper() or (series["instrume"] == "SOXS" and "FLAT" in series["eso dpr type"].upper())):
            if series["lamp"] != "--":
                matchDict['lamp'] = series["lamp"]
                sofName.append(series["lamp"])

        if series["instrume"] == "SOXS":
            sofName.append(series["slit"])
        if series["exptime"] and (series["eso seq arm"].lower() == "nir" or (series["eso seq arm"].lower() == "vis" and ("FLAT" in series["eso dpr type"].upper() or "DARK" in series["eso dpr type"].upper()))):
            matchDict['exptime'] = float(series["exptime"])
            sofName.append(str(series["exptime"]) + "S")
        elif series["exptime"] and "BIAS" not in series["eso dpr type"].upper() and not ("FLAT" in series["eso dpr type"].upper() and series['eso seq arm'].upper() == "UVB"):
            sofName.append(str(series["exptime"]) + "S")

        if series["eso obs name"] != "--":
            matchDict["eso obs name"] = series["eso obs name"]

        if "std" in series["eso dpr type"].lower() or "object" in series["eso dpr type"].lower():
            matchDict["slit"] = series["slit"]

        sofName.append(str(series["instrume"]))

        for k, v in matchDict.items():

            if "type" in k.lower() and "lamp" in v.lower() and "flat" in v.lower():
                mask = (filteredFrames[k].isin(["LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT", "FLAT,LAMP"]))
            else:
                mask = (filteredFrames[k].isin([v]))
            filteredFrames = filteredFrames.loc[mask]

        # INITIAL CALIBRATIONS FILTERING
        if series['eso seq arm'].upper() in ["UVB", "VIS"]:
            for k, v in matchDict.items():
                if k in ["exptime", "lamp", "eso obs name"]:
                    continue
                if "type" in k.lower():
                    mask = (calibrationFrames['eso pro catg'].str.contains("MASTER_")) | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE"))
                elif "slit" in k.lower():
                    mask = ~(~calibrationFrames[k].isin([v]) & (calibrationFrames['eso pro catg'].str.contains("MASTER_MFLAT")))
                elif "rospeed" in k.lower() or "binning" in k.lower():
                    mask = (calibrationFrames[k].isin([v]) | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")))
                else:
                    mask = (calibrationFrames[k].isin([v]))
                calibrationFrames = calibrationFrames.loc[mask]

        elif series['eso seq arm'].upper() in ["NIR"]:
            for k, v in matchDict.items():
                if k in ["binning", "rospeed", "exptime", "lamp", "eso obs name"]:
                    continue
                if "type" in k.lower():
                    mask = (calibrationFrames['eso pro catg'].str.contains("MASTER_") | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")))
                elif "slit" in k.lower():
                    mask = ~(~calibrationFrames[k].isin([v]) & (calibrationFrames['eso pro catg'].str.contains("MASTER_MFLAT")))
                else:
                    mask = (calibrationFrames[k].isin([v]))
                calibrationFrames = calibrationFrames.loc[mask]

        # EXTRA CALIBRATION TABLES
        for k, v in matchDict.items():
            if k in ["rospeed", "exptime", "lamp", "eso obs name", "slit"]:
                continue
            if k in ["binning"] and seriesRecipe in ["mflat"]:
                continue
            if k in ["binning"]:
                mask = (calibrationTables[k].isin([v]) | calibrationTables['eso pro catg'].str.contains("DISP_TAB_"))
            elif "type" in k.lower() and series['eso seq arm'] in ["UVB", "VIS", "NIR"]:
                mask = (calibrationTables['eso pro catg'].str.contains("_TAB_"))
            else:
                try:
                    mask = (calibrationTables[k].isin([v]))
                except:
                    print(k, v)
                    sys.exit(0)
            calibrationTables = calibrationTables.loc[mask]

        # NIGHT START
        # YYYY.MM.DDThh.mm.xxx
        offFrameCount = 0
        if series["night start mjd"]:

            if "PINHOLE" in series["eso dpr tech"].upper():
                if "NIR" not in series['eso seq arm'].upper():
                    mask = (filteredFrames['mjd-obs'] == series["mjd-obs"])
                    filteredFrames = filteredFrames.loc[mask]
                else:
                    # FURTHER FILTERING OF NIR ON/OFF FRAMES
                    filteredFrames["obs-delta"] = filteredFrames['mjd-obs'] - series["mjd-obs"]
                    filteredFrames["obs-delta"] = filteredFrames["obs-delta"].abs()
                    filteredFrames.sort_values(['obs-delta'], inplace=True)
                    mask = (filteredFrames['eso dpr tech'].isin(["IMAGE"]))
                    offFrame = filteredFrames.loc[mask].head(1)
                    onFrame = filteredFrames.loc[~mask].head(1)
                    filteredFrames = pd.concat([onFrame, offFrame], ignore_index=True)
                    offFrameCount = len(offFrame.index)

            if "FLAT" in series["eso dpr tech"].upper() and "NIR" not in series['eso seq arm'].upper():
                mask = (filteredFrames['eso dpr tech'].isin(["IMAGE"]))
                offFrame = filteredFrames.loc[mask]
                from tabulate import tabulate
                print(tabulate(offFrame, headers='keys', tablefmt='psql'))
                sys.exit(0)

            if series["eso dpr tech"] in ["ECHELLE,SLIT,STARE"]:
                mask = (filteredFrames['mjd-obs'] == series["mjd-obs"])
                filteredFrames = filteredFrames.loc[mask]
            else:
                mask = (filteredFrames['night start mjd'] == int(series["night start mjd"]))
                filteredFrames = filteredFrames.loc[mask]
            try:
                frameMjd = filteredFrames["mjd-obs"].values[0]
            except:
                print(series)

            calibrationFrames["obs-delta"] = calibrationFrames['mjd-obs'] - frameMjd
            calibrationTables["obs-delta"] = calibrationTables['mjd-obs'] - frameMjd
            # dispImages["obs-delta"] = dispImages['mjd-obs'] - frameMjd
            calibrationFrames["obs-delta"] = calibrationFrames["obs-delta"].abs()
            calibrationTables["obs-delta"] = calibrationTables["obs-delta"].abs()
            # dispImages["obs-delta"] = dispImages["obs-delta"].abs()
            if isinstance(filteredFrames, astropy.table.row.Row):
                filteredFrames = Table(filteredFrames)

            if seriesRecipe not in ["mbias", "mdark"]:
                mask = (filteredFrames['eso dpr tech'].isin(["IMAGE"]))
            else:
                mask = (filteredFrames['eso dpr tech'].isin(["NONSENSE"]))
            firstDate = filteredFrames.loc[~mask]['date-obs'].values[0].replace("-", ".").replace(":", ".")
            sofName.insert(0, firstDate)

        # NEED SOME FINAL FILTERING ON UVB FLATS
        if "lamp" in series["eso dpr type"].lower() and "flat" in series["eso dpr type"].lower() and "uvb" in series["eso seq arm"].lower():
            mask = (filteredFrames['eso dpr type'].isin(["LAMP,DFLAT"]))
            dFrames = filteredFrames.loc[mask]
            dexptime = dFrames["exptime"].max()
            mask = (filteredFrames['eso dpr type'].isin(["LAMP,QFLAT"]))
            qFrames = filteredFrames.loc[mask]
            qexptime = qFrames["exptime"].max()
            mask = ((filteredFrames['eso dpr type'] == "LAMP,DFLAT") & (filteredFrames['exptime'] == dexptime)) | ((filteredFrames['eso dpr type'] == "LAMP,QFLAT") & (filteredFrames['exptime'] == qexptime))
            filteredFrames = filteredFrames.loc[mask]

        # ADDING TAG FOR SOF FILES
        filteredFrames["tag"] = filteredFrames["eso dpr type"].replace(",", "_") + "_" + filteredFrames["eso seq arm"]
        # LAMP ON AND OFF TAGS
        if series["eso seq arm"].lower() == "nir":
            if seriesRecipe in ["disp_sol", "order_centres", "mflat" "spat_sol"]:
                mask = (filteredFrames['eso dpr tech'].isin(["IMAGE"]))
                filteredFrames.loc[mask, "tag"] += ",OFF"
                filteredFrames.loc[~mask, "tag"] += ",ON"

        # SORT BY COLUMN NAME
        filteredFrames.sort_values(['tag', 'filepath', 'file'],
                                   ascending=[True, True, True], inplace=True)

        files = filteredFrames["file"].values
        filepaths = filteredFrames["filepath"].values
        tags = filteredFrames["tag"].values

        # if series["eso seq arm"].lower() in ("uvb", "vis") and seriesRecipe in ["disp_sol", "order_centres", "spat_sol", "stare"]:
        #     if series['eso dpr type'].lower() == "std,flux":
        #         files = [files[0]]
        #         tags = [tags[0]]
        #         filepaths = [filepaths[0]]

        # if series["eso seq arm"].lower() == "nir" and seriesRecipe in ["disp_sol", "order_centres", "spat_sol", "stare"]:
        #     if series["eso seq arm"].lower() == "std,flux":
        #         files = [files[0]]
        #         tags = [tags[0]]
        #         filepaths = [filepaths[0]]
        #     else:
        #         if len(files) > 1:
        #             files = files[:2]
        #             tags = tags[:2]
        #             filepaths = filepaths[:2]
        #         else:
        #             return series

        if seriesRecipe in ("stare", "nod"):
            object = filteredFrames['object'].values[0].replace(" ", "_")
            sofName.append(object)

        series['sof'] = ("_").join(sofName).replace("-", "").replace(",", "_").upper() + ".sof"
        series["recipe"] = seriesRecipe
        series["recipe_order"] = self.recipeOrder.index(seriesRecipe) + 1

        # CALIBRATIONS NEEDED?
        # BIAS FRAMES
        if series["recipe"] in ["disp_sol", "order_centres", "mflat", "spat_sol", "stare", "nod", "offset"]:

            mask = calibrationFrames['eso pro catg'].isin([f"MASTER_BIAS_{series['eso seq arm'].upper()}"])
            df = calibrationFrames.loc[mask]

            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])
            elif series["eso seq arm"].upper() != "NIR":
                incomplete = True

        # DISP SOLS
        if series["recipe"] in ["order_centres", "spat_sol", "stare", "nod", "offset"]:
            mask = calibrationTables['eso pro catg'].str.contains("DISP_TAB")
            df = calibrationTables.loc[mask]

            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                if series["recipe"] in ["stare", "nod", "offset"]:
                    mask = (df['recipe'] == "spat_sol")
                    if len(df.loc[mask, "file"].index):
                        files = np.append(files, df.loc[mask, "file"].values[0])
                        tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                        filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                    else:
                        incomplete = True
                else:
                    mask = (df['recipe'] == "disp_sol")
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
            mask = calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")
            df = calibrationFrames.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # ORDER TAB
        if series["recipe"] in ["mflat", "spat_sol", "stare", "nod", "offset"]:
            mask = calibrationTables['eso pro catg'].str.contains('ORDER_TAB')
            df = calibrationTables.loc[mask]
            if series["recipe"] in ["mflat"]:
                mask = (df['recipe'] == "order_centres")
            else:
                mask = (df['recipe'] == "mflat")
            df = df.loc[mask]

            if series["recipe"] in ["spat_sol"]:
                # REMOVE BLOCKING FILTERS
                mask = df['slit'].str.contains('JH')
                df = df.loc[~mask]

            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                if series["eso seq arm"].upper() in ["UVB"] and series["recipe"] == "mflat":
                    lamps = ["QLAMP", "DLAMP"]
                    for lamp in lamps:
                        mask = df['file'].str.contains(lamp)
                        if len(df.loc[mask, "file"].index):
                            files = np.append(files, df.loc[mask, "file"].values[0])
                            tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                            filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                elif series["recipe"] == "mflat":
                    files = np.append(files, df["file"].values[0])
                    tags = np.append(tags, df["eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df["filepath"].values[0])
                else:
                    mask = df['filepath'].str.contains("soxs-mflat/")
                    files = np.append(files, df.loc[mask, "file"].values[0])
                    tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
            elif series["recipe"] in ["spat_sol", "mflat"]:
                incomplete = True

        # FLAT FRAMES
        if series["recipe"] in ["spat_sol", "stare", "nod"]:
            mask = calibrationFrames['eso pro catg'].str.contains('MASTER_FLAT')
            df = calibrationFrames.loc[mask]

            if series["recipe"] in ["spat_sol"]:
                # REMOVE BLOCKING FILTERS
                mask = df['slit'].str.contains('JH')
                df = df.loc[~mask]

            if series["recipe"] in ["stare", "nod"]:
                from tabulate import tabulate

                if len(filteredFrames["slit"].values):
                    if self.PAE and self.instrument.upper() == "SOXS":
                        pass
                    else:
                        df = df.loc[(df["slit"] == filteredFrames["slit"].values[0])]

            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # DARK FRAMES
        if series["eso seq arm"].lower() == "nir" and series["recipe"] in ["stare", "disp_sol", "order_centres", "mflat", "spat_sol"]:
            moveOn = False
            if series["recipe"] == "mflat":
                if "FLAT" in series["eso dpr type"].upper() and "NIR" in series['eso seq arm'].upper():
                    mask = (filteredFrames['eso dpr tech'].isin(["IMAGE"]))
                    offFrame = filteredFrames.loc[mask]
                    offFrameCount = len(offFrame.index)
            if series["recipe"] in ["disp_sol", "order_centres", "mflat", "spat_sol"] and offFrameCount > 0:
                moveOn = True
            elif series["recipe"] in ["stare", "disp_sol", "order_centres", "mflat", "spat_sol"] and offFrameCount == 0:
                mask = ((calibrationFrames['exptime'] == filteredFrames['exptime'].max()) & (calibrationFrames['eso pro catg'].str.contains('MASTER_DARK')))
                df = calibrationFrames.loc[mask]
                if len(df.index) == 0:
                    mask = calibrationFrames['eso pro catg'].str.contains('MASTER_DARK')
            else:
                mask = calibrationFrames['eso pro catg'].str.contains('MASTER_DARK')

            if not moveOn:
                df = calibrationFrames.loc[mask]
                if len(df.index):
                    df.sort_values(by=['obs-delta'], inplace=True)
                    files = np.append(files, df["file"].values[0])
                    tags = np.append(tags, df["eso pro catg"].values[0])
                    filepaths = np.append(filepaths, df["filepath"].values[0])

        # CREATE DATA FRAME FROM A DICTIONARY OF LISTS
        myDict = {
            "file": files,
            "filepath": filepaths,
            "tag": tags,
            "sof": [series['sof']] * len(tags),
        }

        if incomplete:
            series["complete"] = 0
        else:
            series["complete"] = 1

        if incomplete:
            myDict["complete"] = [0] * len(tags)
        else:
            myDict["complete"] = [1] * len(tags)

        sofMap = pd.DataFrame(myDict)

        keepTrying = 0
        while keepTrying < 6:
            try:
                sofMap.to_sql(f'sof_map_{self.sessionId}', con=self.conn,
                              index=False, if_exists='append')
                keepTrying = 10
            except Exception as e:
                if keepTrying > 5:
                    raise Exception(e)
                keepTrying += 1

        return series

    def _populate_products_table(
            self,
            series,
            reductionOrder):
        """*determine what the products should be for a given recipe and SOF file and populate the products table*

        **Key Arguments:**

        - ``recipeName`` -- the name of the recipe.
        - ``sofName`` -- the name of the sof file.

        **Return:**

        - None

        """
        self.log.debug('starting the ``_populate_products_table`` method')

        import pandas as pd
        import time

        if series["eso dpr type"].lower() != reductionOrder.lower():
            return series

        if series["complete"] == 0:
            return series

        template = series.copy()
        removeColumns = ["counts", "product", "command", "complete"]

        if template["recipe"] and template["recipe"].lower() in self.productMap:
            for i in removeColumns:
                if i in template:
                    template.pop(i)
            template["eso pro catg"] = template.pop("eso dpr catg")
            template["eso pro tech"] = template.pop("eso dpr tech")
            template["eso pro type"] = template.pop("eso dpr type")

            for i in self.productMap[template["recipe"].lower()]:
                # if template["recipe"].lower() == "mflat" and template["binning"] in ["1x2", "2x2"] and i[2] == "ORDER_TAB":
                #     continue
                products = template.copy()
                products["eso pro type"] = i[0]
                products["eso pro tech"] = i[1]
                products["eso pro catg"] = i[2] + f"_{products['eso seq arm'].upper()}"
                products["file"] = products["sof"].replace(".sof", ".fits")
                if i[4] and i[5]:
                    products["file"] = products["file"].split(i[4])[0] + i[5] + ".fits"
                    products["file"] = products["file"].replace(".fits.fits", ".fits")
                products["filepath"] = "./product/" + i[6] + "/" + products["file"]
                myDict = {k: [v] for k, v in products.items()}
                products = pd.DataFrame(myDict)

                keepTrying = 0
                while keepTrying < 6:
                    try:
                        products.to_sql('product_frames', con=self.conn,
                                        index=False, if_exists='append')
                        keepTrying = 10
                    except Exception as e:

                        if keepTrying > 5:
                            raise Exception(e)
                        time.sleep(1)
                        keepTrying += 1

        self.log.debug('completed the ``_populate_products_table`` method')
        return series

    def _move_misc_files(
            self):
        """*move extra/miscellaneous files to a misc directory*
        """
        self.log.debug('starting the ``_move_misc_files`` method')

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

        self.log.debug('completed the ``_move_misc_files`` method')
        return None

    def _write_sof_files(
            self):
        """*Write out all possible SOF files from the sof_map database table*

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug('starting the ``_write_sof_files`` method')

        import pandas as pd
        from tabulate import tabulate
        import sqlite3 as sql

        conn = sql.connect(
            self.rootDbPath)

        # RECURSIVELY CREATE MISSING DIRECTORIES
        self.sofDir = self.sessionPath + "/sof"
        if not os.path.exists(self.sofDir):
            os.makedirs(self.sofDir)

        df = pd.read_sql_query(f"select * from sof_map_{self.sessionId} where complete = 1;", conn)

        # GROUP RESULTS
        for name, group in df.groupby('sof'):
            myFile = open(self.sofDir + "/" + name, 'w')
            content = tabulate(group[["filepath", "tag"]], tablefmt='plain', showindex=False)
            myFile.write(content)
            myFile.close()

        self.log.debug('completed the ``_write_sof_files`` method')
        return None

    def _write_reduction_shell_scripts(
            self):
        """*write the reduction shell scripts*

        _reduce_all.sh HAS BEEN REPLACED WITH THE `SOXSPIPE REDUCE` COMMAND

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug('starting the ``_write_reduction_shell_scripts`` method')

        import pandas as pd
        import sqlite3 as sql

        conn = sql.connect(
            self.rootDbPath)

        rawGroups = pd.read_sql(
            'SELECT * FROM raw_frame_sets where recipe_order is not null order by recipe_order', con=conn)

        rawGroups["command"] = "soxspipe " + rawGroups["recipe"] + " sof/" + rawGroups["sof"]

        # WRITE FULL REDUCTION SCRIPT
        myFile = open(self.sessionPath + "/_reduce_all.sh", 'w')
        myFile.write(("\n").join(pd.unique(rawGroups["command"])))
        myFile.close()

        self.log.debug('completed the ``_write_reduction_shell_scripts`` method')
        return None

    def session_create(
            self,
            sessionId=False):
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
        self.log.debug('starting the ``session_create`` method')

        import re
        import shutil
        import sqlite3 as sql

        rootDbExists = os.path.exists(self.rootDbPath)
        if rootDbExists:
            # CREATE THE DATABASE CONNECTION
            self.conn = sql.connect(
                self.rootDbPath)
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
            matchObjectList = re.findall(r'[^0-9a-zA-Z\-\_]+', sessionId)
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
                createLogger=False
            )
            arguments, settings, replacedLog, dbConn = su.setup()

        # MAKE ASSET PLACEHOLDERS
        folders = ["sof", "qc", "product"]
        for f in folders:
            if not os.path.exists(self.sessionPath + f"/{f}"):
                os.makedirs(self.sessionPath + f"/{f}")

        # ADD A NEW STATUS COLUMN IN product_frames FOR THIS SESSION
        import sqlite3 as sql
        conn = sql.connect(self.rootDbPath)
        c = conn.cursor()
        sqlQuery = f'ALTER TABLE product_frames ADD status_{sessionId} TEXT;'
        c.execute(sqlQuery)

        # DUPLICATE TEH SOF_MAP TABLE
        sqlQuery = "SELECT sql FROM sqlite_master WHERE type='table' AND name='sof_map'"
        c.execute(sqlQuery)
        sqlQuery = c.fetchall()[0][0]
        sqlQuery = sqlQuery.replace('sof_map', f'sof_map_{sessionId}')
        c.execute(sqlQuery)

        c.close()

        self._write_sof_files()

        # _reduce_all.sh HAS BEEN REPLACED WITH THE `SOXSPIPE REDUCE` COMMAND
        # self._write_reduction_shell_scripts()

        self._symlink_session_assets_to_workspace_root()

        # WRITE THE SESSION ID FILE
        import codecs
        with codecs.open(self.sessionIdFile, encoding='utf-8', mode='w') as writeFile:
            writeFile.write(sessionId)

        message = f"A new data-reduction session has been created with sessionId '{sessionId}'"
        try:
            self.log.print(message)
        except:
            print(message)
        self.log.debug('completed the ``session_create`` method')

        return sessionId

    def session_list(
            self,
            silent=False):
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
        self.log.debug('starting the ``session_list`` method')

        import codecs

        # IF SESSION ID FILE DOES NOT EXIST, REPORT
        self.sessionIdFile = self.sessionsDir + "/.sessionid"
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        else:
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='r') as readFile:
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

        self.log.debug('completed the ``session_list`` method')
        return currentSession, allSessions

    def session_switch(
            self,
            sessionId):
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
        self.log.debug('starting the ``session_switch`` method')
        import codecs

        currentSession, allSessions = self.session_list(silent=True)

        if sessionId == currentSession:
            print(f"Session '{sessionId}' is already in use.")
            return None
        elif sessionId in allSessions:
            # WRITE THE SESSION ID FILE
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='w') as writeFile:
                writeFile.write(sessionId)
        else:
            print(f"There is no session with the ID '{sessionId}'. List existing sessions with `soxspipe session ls`.")
            return None

        self.sessionPath = self.sessionsDir + "/" + sessionId
        self._symlink_session_assets_to_workspace_root()
        print(f"Session successfully switched to '{sessionId}'.")

        self.log.debug('completed the ``session_switch`` method')
        return None

    def _symlink_session_assets_to_workspace_root(
            self):
        """*symlink session QC, product, SOF directories, database and scripts to workspace root*

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug('starting the ``_symlink_session_assets_to_workspace_root`` method')

        import shutil
        import os

        # UNLINK SYMLINK IN ROOT
        for d in os.listdir(self.rootDir):
            filepath = os.path.join(self.rootDir, d)
            if os.path.islink(filepath):
                os.unlink(filepath)

        # SYMLINK FILES AND FOLDERS
        toLink = ["product", "qc", "soxspipe.yaml", "sof", "soxspipe.log"]
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

        self.log.debug('completed the ``_symlink_session_assets_to_workspace_root`` method')
        return None

    def session_refresh(
            self):
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
        self.log.debug('starting the ``session_refresh`` method')

        self.log.print("Refeshing SOF file due to recipe failure\n")

        import codecs
        import sqlite3 as sql

        # IF SESSION ID FILE DOES NOT EXIST, REPORT
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        else:
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='r') as readFile:
                sessionId = readFile.read()
        self.sessionPath = self.sessionsDir + "/" + sessionId
        self.sessionPath = self.sessionsDir + "/" + sessionId
        self.sessionId = sessionId

        self.conn = sql.connect(
            self.rootDbPath)

        # SELECT INSTR
        c = self.conn.cursor()
        sqlQuery = "select instrume from raw_frames where instrume is not null limit 1"
        c.execute(sqlQuery)
        self.instrument = c.fetchall()[0][0]
        c.close()

        if "SOXS" in self.instrument.upper():
            self.typeMap = self.typeMapSOXS
        else:
            self.typeMap = self.typeMapXSH

        # CLEAN UP FAILED FILES
        c = self.conn.cursor()
        sqlQuery = f'delete from sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and p.status_{sessionId} = "fail");'
        c.execute(sqlQuery)
        c.close()

        self._populate_product_frames_db_table()
        self._write_sof_files()

        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("SOF file refresh complete")

        self.log.debug('completed the ``session_refresh`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
