#!/usr/bin/env python
# encoding: utf-8
"""
*The SOXSPIPE Data Organiser*

:Author:
    David Young

:Date Created:
    March  9, 2023
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
    *The worker class for the data_organiser module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``rootDir`` -- the root directory of the data to process

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_).

    To initiate a data_organiser object, use the following:

    ```eval_rst
    .. todo::

        - create cl-util for this class
        - add a tutorial about ``data_organiser`` to documentation
        - create a blog post about what ``data_organiser`` does
    ```

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

        import shutil
        import sqlite3 as sql
        from os.path import expanduser
        import codecs
        from fundamentals.logs import emptyLogger

        log.debug("instantiating a new 'data_organiser' object")
        self.log = emptyLogger()

        # if rootDir == ".":
        #     rootDir = os.getcwd()
        # MAKE RELATIVE HOME PATH ABSOLUTE
        if rootDir[0] == "~":
            home = expanduser("~")
            directory = directory.replace("~", home)

        self.rootDir = rootDir
        self.rawDir = rootDir + "/raw_frames"
        self.miscDir = rootDir + "/misc"
        self.sofDir = rootDir + "/sof"
        self.sessionsDir = rootDir + "/sessions"

        # MK RAW FRAME DIRECTORY
        if not os.path.exists(self.rawDir):
            os.makedirs(self.rawDir)

        # MK RAW FRAME DIRECTORY
        if not os.path.exists(self.sessionsDir):
            os.makedirs(self.sessionsDir)

        # SESSION ID PLACEHOLDER FILE
        self.sessionIdFile = self.sessionsDir + "/.sessionid"

        # TEST FOR SQLITE DATABASE
        self.dbPath = rootDir + "/soxspipe.db"
        try:
            with open(self.dbPath):
                pass
            self.freshRun = False
        except IOError:
            self.freshRun = True
            emptyDb = os.path.dirname(os.path.dirname(__file__)) + "/resources/soxspipe.db"
            shutil.copyfile(emptyDb, self.dbPath)

        def dict_factory(cursor, row):
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        # CREATE THE DATABASE CONNECTION
        self.conn = sql.connect(
            self.dbPath)
        self.conn.row_factory = dict_factory

        # HERE ARE THE KEYS WE WANT PRESENTED IN THE SUMMARY OUTPUT TABLES
        self.keywords = [
            'file',
            'mjd-obs',
            'date-obs',
            'eso seq arm',
            'eso dpr catg',
            'eso dpr tech',
            'eso dpr type',
            'eso pro catg',
            'eso pro tech',
            'eso pro type',
            'exptime',
            'cdelt1',
            'cdelt2',
            'eso det read speed',
            'eso ins opti3 name',
            'eso ins opti4 name',
            'eso ins opti5 name',
            'eso det ncorrs name',
            'eso det out1 conad',
            'eso det out1 ron',
            'eso obs id',
            'eso obs name',
            "naxis",
            "object"
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
            'night start date',
            'night start mjd',
            'mjd-obs',
            'date-obs',
            'object'
        ]

        # THIS TYPE MAP WILL BE USED TO GROUP SET OF FILES TOGETHER
        self.typeMap = {
            "bias": [{"tech": None, "catg": None, "recipe": "mbias"}],
            "dark": [{"tech": None, "catg": None, "recipe": "mdark"}],
            "lamp,fmtchk": [{"tech": None, "catg": None, "recipe": "disp_sol"}],
            "lamp,orderdef": [{"tech": None, "catg": None, "recipe": "order_centres"}],
            "lamp,dorderdef": [{"tech": None, "catg": None, "recipe": "order_centres"}],
            "lamp,qorderdef": [{"tech": None, "catg": None, "recipe": "order_centres"}],
            "lamp,flat": [{"tech": None, "catg": None, "recipe": "mflat"}],
            "lamp,dflat": [{"tech": None, "catg": None, "recipe": "mflat"}],
            "lamp,qflat": [{"tech": None, "catg": None, "recipe": "mflat"}],
            "lamp,wave": [{"tech": ["echelle,multi-pinhole", "image"], "catg": None, "recipe": "spat_sol"}],
            "object": [{"tech": ["echelle,slit,stare"], "catg": None, "recipe": "stare"}],
            "std,flux": [{"tech": ["echelle,slit,stare"], "catg": None, "recipe": "stare"}],
            "std,telluric": [{"tech": ["echelle,slit,stare"], "catg": None, "recipe": "stare"}],
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
                ["REDUCED", "ECHELLE,SLIT", "OBJECT_TAB", "TABLE", None, None, "soxs-stare"]
            ],
        }

        self.proKeywords = ['eso pro type', 'eso pro tech', 'eso pro catg']

        # THESE ARE KEYS WE NEED TO FILTER ON, AND SO NEED TO CREATE ASTROPY TABLE
        # INDEXES
        self.filterKeywords = ['eso seq arm', 'eso dpr catg',
                               'eso dpr tech', 'eso dpr type', 'eso pro catg', 'eso pro tech', 'eso pro type', 'exptime', 'rospeed', 'binning', 'night start mjd', 'night start date']

        # THIS IS THE ORDER TO PROCESS THE FRAME TYPES
        self.reductionOrder = ["BIAS", "DARK", "LAMP,FMTCHK", "LAMP,ORDERDEF", "LAMP,DORDERDEF", "LAMP,QORDERDEF", "LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT", "LAMP,WAVE", "STD,FLUX", "STD,TELLURIC", "OBJECT"]

        # THIS IS THE ORDER THE RECIPES NEED TO BE RUN IN (MAKE SURE THE REDUCTION SCRIPT HAS RECIPES IN THE CORRECT ORDER)
        self.recipeOrder = ["mbias", "mdark", "disp_sol", "order_centres", "mflat", "spat_sol", "stare"]

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

        # IF SESSION ID FILE DOES NOT EXIST, CREATE A NEW SESSION
        # OTHERWISE USE CURRENT SESSION
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            sessionId = self.session_create(sessionId="base")
        else:
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='r') as readFile:
                sessionId = readFile.read()

        basename = os.path.basename(self.rootDir)
        print(f"PREPARING THE `{basename}` WORKSPACE FOR DATA-REDUCTION")
        self._sync_raw_frames()
        self._move_misc_files()
        self._populate_product_frames_db_table()
        self._populate_product_frames_db_table()
        self._write_sof_files()
        self._write_reduction_shell_scripts()

        print(f"\nTHE `{basename}` WORKSPACE FOR HAS BEEN PREPARED FOR DATA-REDUCTION\n")
        print(f"In this workspace you will find:\n")
        print(f"   - `misc/`: a lost-and-found archive of non-fits files")
        print(f"   - `raw_frames/`: all raw-frames to be reduced")
        print(f"   - `sessions/`: directory of data-reduction sessions")
        print(f"   - `sof/`: the set-of-files (sof) files required for each reduction step")
        print(f"   - `_reduce_all.sh`: a single script to reduce all the data in the workspace")
        print(f"   - `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete\n")

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

        # GENERATE AN ASTROPY TABLES OF FITS FRAMES WITH ALL INDEXES NEEDED
        filteredFrames, fitsPaths, fitsNames = self.create_directory_table(pathToDirectory=self.rootDir, keys=self.keywords, filterKeys=self.filterKeywords)

        if fitsPaths:
            # SPLIT INTO RAW, REDUCED PIXELS, REDUCED TABLES
            rawFrames, reducedFramesPixels, reducedFramesTables = self.categorise_frames(filteredFrames)

            if len(rawFrames.index):
                rawFrames["filepath"] = "./raw_frames/" + rawFrames['file']
                rawFrames.to_sql('raw_frames', con=self.conn,
                                 index=False, if_exists='append')
                filepaths = rawFrames['file'].values
                for f in filepaths:
                    shutil.move(self.rootDir + "/" + f, self.rawDir)

        if not skipSqlSync:
            self.sync_sql_table_to_directory(self.rawDir, 'raw_frames', recursive=False)

        self.log.debug('completed the ``_sync_raw_frames`` method')
        return None

    def create_directory_table(
            self,
            pathToDirectory,
            keys,
            filterKeys):
        """*create an astropy table based on the contents of a directory*

        **Key Arguments:**

        - `log` -- logger
        - `pathToDirectory` -- path to the directory containing the FITS frames
        - `keys` -- the keys needed to be returned for the imageFileCollection
        - `filterKeys` -- these are the keywords we want to filter on later

        **Return**

        - `masterTable` -- the primary dataframe table listing all FITS files in the directory (including indexes on `filterKeys` columns)
        - `fitsPaths` -- a simple list of all FITS file paths
        - `fitsNames` -- a simple list of all FITS file name

        **Usage:**

        ```python
        # GENERATE AN ASTROPY TABLES OF FITS FRAMES WITH ALL INDEXES NEEDED
        masterTable, fitsPaths, fitsNames = create_directory_table(
            log=log,
            pathToDirectory="/my/directory/path",
            keys=["file","mjd-obs", "exptime","cdelt1", "cdelt2"],
            filterKeys=["mjd-obs","exptime"]
        )
        ```
        """
        self.log.debug('starting the ``create_directory_table`` function')

        from ccdproc import ImageFileCollection
        from astropy.time import Time, TimeDelta
        import numpy as np
        from fundamentals.files import recursive_directory_listing
        import pandas as pd

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

        # TOP-LEVEL COLLECTION
        if recursive:
            allFrames = ImageFileCollection(filenames=fitsPaths, keywords=keys)
        else:
            allFrames = ImageFileCollection(
                location=pathToDirectory, filenames=fitsNames, keywords=keys)
        masterTable = allFrames.summary

        # ADD FILLED VALUES FOR MISSING CELLS
        for fil in keys:
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

        # SETUP A NEW COLUMN GIVING THE INT MJD THE CHILEAN NIGHT BEGAN ON
        # 12:00 NOON IN CHILE IS TYPICALLY AT 16:00 UTC (CHILE = UTC - 4)
        # SO COUNT CHILEAN OBSERVING NIGHTS AS 15:00 UTC-15:00 UTC (11am-11am)
        if "mjd-obs" in masterTable.colnames:
            chile_offset = TimeDelta(4.0 * 60 * 60, format='sec')
            night_start_offset = TimeDelta(15.0 * 60 * 60, format='sec')
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

        if "eso det read speed" in masterTable.colnames:
            masterTable["rospeed"] = np.copy(masterTable["eso det read speed"])
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

        if "naxis" in masterTable.colnames:
            masterTable["table"] = np.copy(masterTable["naxis"]).astype(str)
            masterTable["table"][masterTable[
                "table"] == '0'] = 'T'
            masterTable["table"][masterTable[
                "table"] != 'T'] = 'F'

        if "cdelt2" in masterTable.colnames:
            masterTable["binning"] = np.core.defchararray.add(
                masterTable["cdelt1"].astype('int').astype('str'), "x")
            masterTable["binning"] = np.core.defchararray.add(masterTable["binning"],
                                                              masterTable["cdelt2"].astype('int').astype('str'))
            masterTable["binning"][masterTable[
                "binning"] == '-99x-99'] = '--'
            masterTable["binning"][masterTable[
                "binning"] == '1x-99'] = '--'
            masterTable.add_index("binning")

        # ADD INDEXES ON ALL KEYS
        for k in keys:
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

        self.log.debug('completed the ``create_directory_table`` function')
        return masterTable, fitsPaths, fitsNames

    def sync_sql_table_to_directory(
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
        do.sync_sql_table_to_directory('/raw/directory/', 'raw_frames', recursive=False)
        ```
        """
        self.log.debug('starting the ``sync_sql_table_to_directory`` method')

        import sqlite3 as sql
        import shutil
        import time

        # GENERATE A LIST OF FITS FILE PATHS IN RAW DIR
        fitsPaths = []
        fitsNames = []

        folderName = os.path.basename(directory)
        for d in os.listdir(directory):
            absfilepath = os.path.join(directory, d)
            filepath = f"./{folderName}/{d}"
            if os.path.isfile(absfilepath) and (os.path.splitext(absfilepath)[1] == ".fits" or ".fits.Z" in absfilepath):
                fitsPaths.append(filepath)
                fitsNames.append(d)

        c = self.conn.cursor()

        sqlQuery = f"select filepath from {tableName};"
        c.execute(sqlQuery)
        dbFiles = [r["filepath"] for r in c.fetchall()]

        # DELETED FILES
        filesNotInDB = set(fitsPaths) - set(dbFiles)
        filesNotInFS = set(dbFiles) - set(fitsPaths)
        if len(filesNotInFS):
            filesNotInFS = ("','").join(filesNotInFS)
            sqlQuery = f"delete from {tableName} where filepath in ('{filesNotInFS}');"
            c.execute(sqlQuery)

        if len(filesNotInDB):
            for f in filesNotInDB:
                shutil.move(self.rootDir + "/" + f, self.rootDir)
            self._sync_raw_frames(skipSqlSync=True)

        c.close()

        self.log.debug('completed the ``sync_sql_table_to_directory`` method')
        return None

    def categorise_frames(
            self,
            filteredFrames,
            verbose=False):
        """*given a dataframe of frame, categorise frames into raw, reduced pixels, reduced tables*

        **Key Arguments:**
            - ``filteredFrames`` -- the dataframe from which to split frames into categorise.
            - ``verbose`` -- print restuls to stdout.

        **Return:**
            - ``rawFrames`` -- dataframe of raw frames only
            - ``reducedFramesPixels`` -- dataframe of reduced images only
            - ``reducedFramesTables`` -- dataframe of reduced tables only

        **Usage:**

        ```python
        rawFrames, reducedFramesPixels, reducedFramesTables = self.categorise_frames(filteredFrames)
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

        mask = []
        for i in self.proKeywords:
            keywordsTerseRaw.remove(i)
            filterKeywordsRaw.remove(i)
            if not len(mask):
                mask = (filteredFrames[i] == "--")
            else:
                mask = np.logical_and(mask, (filteredFrames[i] == "--"))

        rawFrames = filteredFrames.loc[mask]

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
            self.log.print("\n# CONTENT SETS INDEX\n")

        if verbose and len(reducedPixelsGroups.index):
            self.log.print("\n## REDUCED PIXEL-FRAME-SET SUMMARY\n")
            self.log.print(tabulate(reducedPixelsGroups, headers='keys', tablefmt='github', showindex=False, stralign="right"))

        if verbose and len(reducedTablesGroups.index):
            self.log.print("\n## REDUCED TABLES-SET SUMMARY\n")
            self.log.print(tabulate(reducedTablesGroups, headers='keys', tablefmt='github', showindex=False, stralign="right"))

        if verbose:
            self.log.print("\n# CONTENT FILE INDEX\n")
        if verbose and len(rawGroups.index):
            self.log.print("\n## ALL RAW FRAMES\n")
            self.log.print(tabulate(rawFrames[keywordsTerseRaw], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

        if verbose and len(reducedPixelsGroups.index):
            self.log.print("\n## ALL REDUCED PIXEL-FRAMES\n")
            self.log.print(tabulate(reducedFramesPixels[keywordsTerseReduced], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

        if verbose and len(reducedTablesGroups.index):
            self.log.print("\n## ALL REDUCED TABLES\n")
            self.log.print(tabulate(reducedFramesTables[keywordsTerseReducedTable], headers='keys', tablefmt='github', showindex=False, stralign="right", floatfmt=".3f"))

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

        ---

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - write a command-line tool for this method
            - update package tutorial with command-line tool info if needed
        ```
        """
        self.log.debug('starting the ``_populate_product_frames_db_table`` method')

        import pandas as pd
        import sqlite3 as sql

        conn = sql.connect(
            self.dbPath)

        rawFrames = pd.read_sql('SELECT * FROM raw_frames', con=conn)

        rawFrames.fillna("--", inplace=True)
        filterKeywordsRaw = self.filterKeywords[:]

        for i in self.proKeywords:
            filterKeywordsRaw.remove(i)

        # rawFrames.replace("LAMP,DFLAT", "LAMP,FLAT", inplace=True)
        # rawFrames.replace("LAMP,QFLAT", "LAMP,FLAT", inplace=True)
        rawGroups = rawFrames.groupby(filterKeywordsRaw)

        mjds = rawGroups.mean(numeric_only=True)["mjd-obs"].values

        rawGroups = rawGroups.size().reset_index(name='counts')
        rawGroups['mjd-obs'] = mjds

        # REMOVE GROUPED STARE AND NODDING - NEED TO ADD INDIVIDUAL FRAMES TO GROUP
        mask = (rawGroups["eso dpr tech"].isin(["ECHELLE,SLIT,STARE", "ECHELLE,SLIT,NODDING"]))
        rawGroups = rawGroups.loc[~mask]
        # NOW ADD SCIENCE FRAMES AS ONE ENTRY PER EXPOSURE
        rawScienceFrames = pd.read_sql(
            'SELECT * FROM raw_frames where "eso dpr tech" in ("ECHELLE,SLIT,STARE","ECHELLE,SLIT,NODDING")', con=conn)

        rawScienceFrames.fillna("--", inplace=True)
        rawScienceFrames = rawScienceFrames.groupby(filterKeywordsRaw + ["mjd-obs"])
        rawScienceFrames = rawScienceFrames.size().reset_index(name='counts')

        # MERGE DATAFRAMES
        rawGroups = pd.concat([rawGroups, rawScienceFrames], ignore_index=True)
        rawGroups['recipe'] = None
        rawGroups['sof'] = None

        calibrationFrames = pd.read_sql('SELECT * FROM product_frames where `eso pro catg` not like "%_TAB_%"', con=conn)
        calibrationFrames.fillna("--", inplace=True)

        calibrationTables = pd.read_sql('SELECT * FROM product_frames where `eso pro catg` like "%_TAB_%"', con=conn)
        calibrationTables.fillna("--", inplace=True)

        # generate_sof_and_product_names SHOULD TAKE ROW OF DF AS INPUT
        for o in self.reductionOrder:
            rawGroups = rawGroups.apply(self.generate_sof_and_product_names, axis=1, reductionOrder=o, rawFrames=rawFrames, calibrationFrames=calibrationFrames, calibrationTables=calibrationTables)
            rawGroups = rawGroups.apply(self.populate_products_table, axis=1, reductionOrder=o)

        # SEND TO DATABASE
        c = self.conn.cursor()
        sqlQuery = f"delete from raw_frame_sets;"
        c.execute(sqlQuery)
        c.close()
        rawGroups.replace(['--'], None).to_sql('raw_frame_sets', con=self.conn,
                                               index=False, if_exists='append')

        self.log.debug('completed the ``_populate_product_frames_db_table`` method')
        return None

    def generate_sof_and_product_names(
            self,
            series,
            reductionOrder,
            rawFrames,
            calibrationFrames,
            calibrationTables):
        """*add a recipe name and SOF filename to all rows in the raw_frame_sets DB table*

        **Key Arguments:**

        - ``series`` -- the dataframe row/series to apply work on

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
        ```

        ```python
        usage code
        ```
        """

        import pandas as pd
        import astropy
        import numpy as np

        sofName = []
        matchDict = {}
        sofName.append(series['eso seq arm'].upper())
        matchDict['eso seq arm'] = series['eso seq arm'].upper()
        filteredFrames = rawFrames.copy()

        if series["eso dpr type"].lower() != reductionOrder.lower():
            return series

        # FILTER BY TYPE FIRST
        if "FLAT" in series["eso dpr type"].upper():
            mask = (filteredFrames["eso dpr type"].str.contains("FLAT"))
        else:
            mask = (filteredFrames["eso dpr type"].isin([series["eso dpr type"].upper()]))
        filteredFrames = filteredFrames.loc[mask]

        # CHECK TECH
        if self.typeMap[series["eso dpr type"].lower()][0]["tech"]:
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
        else:
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
            sofName.append(self.typeMap[series["eso dpr type"].lower()][0]["recipe"].replace("_centres", "_locations"))
            if "DORDER" in series["eso dpr type"].upper():
                sofName.append("dlamp")
            if "QORDER" in series["eso dpr type"].upper():
                sofName.append("qlamp")
        if series["exptime"] and (series["eso seq arm"].lower() == "nir" or (series["eso seq arm"].lower() == "vis" and "FLAT" in series["eso dpr type"].upper())):
            matchDict['exptime'] = float(series["exptime"])
            sofName.append(str(series["exptime"]).replace(".", "pt"))

        for k, v in matchDict.items():
            if "type" in k.lower() and "lamp" in v.lower() and "flat" in v.lower():
                mask = (filteredFrames[k].isin(["LAMP,FLAT", "LAMP,DFLAT", "LAMP,QFLAT"]))
            else:
                mask = (filteredFrames[k].isin([v]))
            filteredFrames = filteredFrames.loc[mask]

        # INITIAL CALIBRATIONS FILTERING
        if series['eso seq arm'].upper() in ["UVB", "VIS"]:
            for k, v in matchDict.items():
                if k in ["exptime"]:
                    continue
                if "type" in k.lower():
                    mask = (calibrationFrames['eso pro catg'].str.contains("MASTER_")) | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE"))
                elif "rospeed" in k.lower() or "binning" in k.lower():
                    mask = (calibrationFrames[k].isin([v]) | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")))
                else:
                    mask = (calibrationFrames[k].isin([v]))
                calibrationFrames = calibrationFrames.loc[mask]

        elif series['eso seq arm'].upper() in ["NIR"]:
            for k, v in matchDict.items():
                if k in ["binning", "rospeed", "exptime"]:
                    continue
                if "type" in k.lower():
                    mask = (calibrationFrames['eso pro catg'].str.contains("MASTER_") | (calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")))
                else:
                    mask = (calibrationFrames[k].isin([v]))
                calibrationFrames = calibrationFrames.loc[mask]

        # EXTRA CALIBRATION TABLES
        for k, v in matchDict.items():
            if k in ["rospeed", "exptime"]:
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
        if series["night start mjd"]:
            if series["eso dpr tech"] in ["ECHELLE,SLIT,STARE", "ECHELLE,SLIT,NODDING"]:
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

            # SELECT FIRST DATE
            # SORT BY COLUMN NAME
            filteredFrames.sort_values(['date-obs'], inplace=True)
            firstDate = filteredFrames['date-obs'].values[0].replace("-", ".").replace(":", ".")
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

        if seriesRecipe == "stare":
            object = filteredFrames['object'].values[0].replace(" ", "_")
            sofName.append(object)

        series['sof'] = ("_").join(sofName).replace("-", "").replace(",", "_").upper() + ".sof"
        series["recipe"] = seriesRecipe
        series["recipe_order"] = self.recipeOrder.index(seriesRecipe) + 1

        # CALIBRATIONS NEEDED?
        # BIAS FRAMES
        if series["recipe"] in ["disp_sol", "order_centres", "mflat", "spat_sol", "stare"]:
            mask = calibrationFrames['eso pro catg'].isin([f"MASTER_BIAS_{series['eso seq arm'].upper()}"])
            df = calibrationFrames.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # DISP SOLS
        if series["recipe"] in ["order_centres", "spat_sol", "stare"]:
            mask = calibrationTables['eso pro catg'].str.contains("DISP_TAB")
            df = calibrationTables.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                if series["recipe"] in ["stare"]:
                    mask = (df['recipe'] == "spat_sol")
                    if len(df.loc[mask, "file"].index):
                        files = np.append(files, df.loc[mask, "file"].values[0])
                        tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                        filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])
                else:
                    mask = (df['recipe'] == "disp_sol")
                    if len(df.loc[mask, "file"].index):
                        files = np.append(files, df.loc[mask, "file"].values[0])
                        tags = np.append(tags, df.loc[mask, "eso pro catg"].values[0])
                        filepaths = np.append(filepaths, df.loc[mask, "filepath"].values[0])

        # DISP SOLS IMAGE
        if series["recipe"] in ["stare"]:
            mask = calibrationFrames['eso pro catg'].str.contains("DISP_IMAGE")
            df = calibrationFrames.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # ORDER TAB
        if series["recipe"] in ["mflat", "spat_sol", "stare"]:
            mask = calibrationTables['eso pro catg'].str.contains('ORDER_TAB')
            df = calibrationTables.loc[mask]
            if series["recipe"] in ["mflat"]:
                mask = (df['recipe'] == "order_centres")
            else:
                mask = (df['recipe'] == "mflat")
            df = df.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                if series["eso seq arm"].upper() in ["UVB"] and series["recipe"] == "mflat":
                    lamps = ["QLAMP", "DLAMP"]
                    for lamp in lamps:
                        mask = df['file'].str.contains(lamp)
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

        # FLAT FRAMES
        if series["recipe"] in ["spat_sol", "stare"]:
            mask = calibrationFrames['eso pro catg'].str.contains('MASTER_FLAT')
            df = calibrationFrames.loc[mask]
            if len(df.index):
                df.sort_values(by=['obs-delta'], inplace=True)
                files = np.append(files, df["file"].values[0])
                tags = np.append(tags, df["eso pro catg"].values[0])
                filepaths = np.append(filepaths, df["filepath"].values[0])

        # DARK FRAMES
        if series["eso seq arm"].lower() == "nir" and series["recipe"] in ["stare"]:
            mask = calibrationFrames['eso pro catg'].str.contains('MASTER_DARK')
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
            "sof": [series['sof']] * len(tags)
        }

        sofMap = pd.DataFrame(myDict)
        sofMap.to_sql('sof_map', con=self.conn,
                      index=False, if_exists='append')

        return series

    def populate_products_table(
            self,
            series,
            reductionOrder):
        """*determine what the products should be for a given recipe and SOF file and ppulate the products table*

        **Key Arguments:**
            - ``recipeName`` -- the name of the recipe.
            - ``sofName`` -- the name of the sof file.

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
        self.log.debug('starting the ``populate_products_table`` method')

        import pandas as pd

        if series["eso dpr type"].lower() != reductionOrder.lower():
            return series

        template = series.copy()
        removeColumns = ["counts", "product", "command"]

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
                products.to_sql('product_frames', con=self.conn,
                                index=False, if_exists='append')

        self.log.debug('completed the ``populate_products_table`` method')
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
            if os.path.isfile(filepath) and os.path.splitext(filepath)[1] != ".db" and "readme" not in d.lower():
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
        self.log.debug('starting the ``_write_sof_files`` method')

        import pandas as pd
        from tabulate import tabulate
        import sqlite3 as sql

        conn = sql.connect(
            self.dbPath)

        # RECURSIVELY CREATE MISSING DIRECTORIES
        if not os.path.exists(self.sofDir):
            os.makedirs(self.sofDir)

        df = pd.read_sql_query("select * from sof_map;", conn)

        # GROUP RESULTS
        for name, group in df.groupby('sof'):
            myFile = open(self.sofDir + "/" + name, 'w')
            content = tabulate(group[["filepath", "tag"]], tablefmt='plain', showindex=False)
            myFile.write(content)
            myFile.close()

        self.log.debug('completed the ``_write_sof_files`` method')
        return None

    # use the tab-trigger below for new method
    def _write_reduction_shell_scripts(
            self):
        """*write the reduction shell scripts*

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
        self.log.debug('starting the ``_write_reduction_shell_scripts`` method')

        import pandas as pd
        import sqlite3 as sql

        conn = sql.connect(
            self.dbPath)

        rawGroups = pd.read_sql(
            'SELECT * FROM raw_frame_sets where recipe_order is not null order by recipe_order', con=conn)

        rawGroups["command"] = "soxspipe " + rawGroups["recipe"] + " sof/" + rawGroups["sof"]

        # WRITE FULL REDUCTION SCRIPT
        myFile = open(self.rootDir + "/_reduce_all.sh", 'w')
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

        # MAKE THE SESSION DIRECTORY
        sessionPath = self.sessionsDir + "/" + sessionId
        if not os.path.exists(sessionPath):
            os.makedirs(sessionPath)

        # SETUP SESSION SETTINGS AND LOGGING
        testPath = sessionPath + "/soxspipe.yaml"
        exists = os.path.exists(testPath)
        if not exists:
            su = tools(
                arguments={"<workspaceDirectory>": sessionPath, "init": True, "settingsFile": None},
                docString=False,
                logLevel="WARNING",
                options_first=False,
                projectName="soxspipe",
                defaultSettingsFile=True
            )
            arguments, settings, replacedLog, dbConn = su.setup()

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
        elif sessionId in allSessions:
            # WRITE THE SESSION ID FILE
            with codecs.open(self.sessionIdFile, encoding='utf-8', mode='w') as writeFile:
                writeFile.write(sessionId)
            print(f"Session successfully switched to '{sessionId}'.")
        else:
            print(f"There is no session with the ID '{sessionId}'. List existing sessions with `soxspipe session ls`.")

        self.log.debug('completed the ``session_switch`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
