#!/usr/bin/env python
"""
*The SOXSPIPE Data Organiser*

Author
: David Young

Date Created
: March  9, 2023
"""

import os
import shutil
import sqlite3
import sys
import traceback
from contextlib import closing
from pathlib import Path

from fundamentals import tools

from soxspipe.commonutils.late_frame_rejoin import (
    delete_stale_files,
    find_rejoined_sofs,
    format_rejoin_summary,
    pin_rejoined_sof_names,
    reset_rejoined_sofs,
)
from soxspipe.commonutils.missing_calibrations import (
    CALIBRATION_RECIPES,
    FAILED_QC,
    FAILED_RUN,
    NOT_YET_REDUCED,
    QC_FAILURE_MESSAGE_PREFIX,
    calibration_match_clause,
    describe_missing,
    find_missing_calibrations,
    missing_reasons,
)
from soxspipe.commonutils.session_status_restore import (
    classify_session_sofs,
    format_session_summary,
    read_frame_sets,
    read_session_snapshots,
)
from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

os.environ["TERM"] = "vt100"


def print_incomplete_sets_report(incompleteSets, blockingSets=None):
    """*print the "cannot be reduced" report shared by `data_organiser.prepare` and `reducer.reduce`*

    **Key Arguments:**

    - ``incompleteSets`` -- dataframe from `data_organiser.get_incomplete_raw_frames_set`
    - ``blockingSets`` -- optional dataframe from `data_organiser.get_blocking_calibration_sets`. Printed as a
      second table when given and non-empty.
    """
    if incompleteSets is None or not len(incompleteSets.index):
        return

    from tabulate import tabulate

    print(
        "\nSOME CALIBRATION FRAMES ARE NOT PRESENT (OR FAILED TO BE BUILT) FOR THE FOLLOWING DATA SETS "
        "AND THEY CANNOT BE REDUCED:"
    )
    print(
        tabulate(
            incompleteSets,
            headers="keys",
            tablefmt="pretty",
            showindex=False,
        )
    )

    if blockingSets is not None and len(blockingSets.index):
        print("\nCALIBRATIONS BLOCKING THESE DATA SETS:")
        print(
            tabulate(
                blockingSets,
                headers="keys",
                tablefmt="pretty",
                showindex=False,
            )
        )


class _UnsafePathError(ValueError):
    """Report an unsafe path or workspace identifier at a trust boundary."""


class DatabasePreservationError(sqlite3.DatabaseError):
    """The workspace database could not be copied into `backups/`, so it was left in place and not rebuilt."""


class _UnreadableDatabaseError(sqlite3.DatabaseError):
    """A database that opens but whose SQLite snapshot fails its `quick_check`."""


# LABEL IN THE NAME OF A RAW, BYTE-FOR-BYTE BACKUP OF A DATABASE THAT FAILED TO OPEN
_RAW_BACKUP_LABEL = "corrupt"
# SQLITE PRIMARY RESULT CODES THAT MEAN THE FILE ITSELF IS DAMAGED (NOT LOCKED OR BUSY)
_UNREADABLE_DATABASE_CODES = {sqlite3.SQLITE_NOTADB, sqlite3.SQLITE_CORRUPT}
# SQLITE PRIMARY RESULT CODES THAT MEAN ANOTHER CONNECTION HOLDS THE DATABASE
_LOCKED_DATABASE_CODES = {sqlite3.SQLITE_BUSY, sqlite3.SQLITE_LOCKED}
_SQLITE_PRIMARY_CODE_MASK = 0xFF


def _is_unreadable_database_error(error):
    """Return True when a SQLite error says the file is not a database or is malformed."""
    if isinstance(error, _UnreadableDatabaseError):
        return True
    errorCode = getattr(error, "sqlite_errorcode", None)
    return errorCode is not None and (errorCode & _SQLITE_PRIMARY_CODE_MASK) in _UNREADABLE_DATABASE_CODES


def _is_locked_database_error(error):
    """Return True when a SQLite error says the database is locked or busy."""
    errorCode = getattr(error, "sqlite_errorcode", None)
    return errorCode is not None and (errorCode & _SQLITE_PRIMARY_CODE_MASK) in _LOCKED_DATABASE_CODES


def _fsync_path(path):
    """Flush a file or directory to disk."""
    descriptor = os.open(path, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def _validate_owned_path(path, owner, label):
    """Return a path only when its resolved target remains below its owner."""
    candidatePath = Path(path)
    ownerPath = Path(owner)
    try:
        resolvedOwnerPath = ownerPath.resolve()
        resolvedCandidatePath = candidatePath.resolve()
        resolvedCandidatePath.relative_to(resolvedOwnerPath)
    except (OSError, RuntimeError, ValueError) as error:
        raise _UnsafePathError(f"Unsafe {label}: path resolves outside {ownerPath}") from error
    return resolvedCandidatePath


def _validate_session_id(sessionId):
    """Enforce the documented grammar for workspace session identifiers."""
    import re

    if not isinstance(sessionId, str) or re.fullmatch(r"[0-9A-Za-z_]{1,16}", sessionId) is None:
        raise _UnsafePathError("Session ID must be 16 characters long or shorter, consisting of A-Z, a-z, 0-9 and/or _")
    return sessionId


class data_organiser:
    """
    *The `soxspipe` Data Organiser*

    **Key Arguments:**

    - ``log`` -- logger
    - ``rootDir`` -- the root directory of the data to process
    - ``vlt`` -- prepare the workspace using the standard vlt /data directory

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/main/initialisation.html).

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

    # ATTEMPTS TO OPEN AND CHECK THE DATABASE BEFORE IT IS REBUILT, AND SECONDS EACH ATTEMPT WAITS ON A LOCK
    _DB_OPEN_ATTEMPTS = 50
    _DB_BUSY_TIMEOUT_SECONDS = 300
    # MOST RESTORE-AND-REBUILD PASSES TO RUN WHILE RESTORED FAILURES KEEP RE-MATCHING DOWNSTREAM CALIBRATIONS
    _STATUS_RESTORE_MAX_PASSES = 10

    def __init__(self, log, rootDir, vlt=False, dbConnect=True):
        import codecs
        import warnings
        from os.path import expanduser

        from astropy.utils.exceptions import AstropyWarning

        warnings.simplefilter("ignore", AstropyWarning)

        self.vlt = vlt
        self.PAE = False
        log.debug("instantiating a new 'data_organiser' object")
        self.log = log
        # NAMES PINNED FOR REJOINED SETS LIVE ONLY ON THIS INSTANCE; THEY ARE LOST WHEN THE PROCESS EXITS
        self.pinnedSofNames = {}

        # MAKE RELATIVE HOME PATH ABSOLUTE
        if rootDir[0] == "~":
            home = expanduser("~")
            directory = directory.replace("~", home)

        self.rootDir = rootDir
        self.rawDir = str(_validate_owned_path(Path(rootDir) / "raw", rootDir, "raw directory"))
        self.miscDir = str(_validate_owned_path(Path(rootDir) / "misc", rootDir, "misc directory"))
        self.sessionsDir = str(_validate_owned_path(Path(rootDir) / "sessions", rootDir, "sessions directory"))

        if self.vlt:
            self.vltReduced = self.use_vlt_environment_folders()

        # SESSION ID PLACEHOLDER FILE
        self.sessionIdFile = str(
            _validate_owned_path(
                Path(self.sessionsDir) / ".sessionid",
                self.sessionsDir,
                "session ID path",
            )
        )
        exists = os.path.exists(self.sessionIdFile)
        if exists:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                sessionId = _validate_session_id(readFile.read())
                self.sessionPath = str(
                    _validate_owned_path(
                        Path(self.sessionsDir) / sessionId,
                        self.sessionsDir,
                        "session path",
                    )
                )
                self.sessionId = sessionId

        # DATABASE FILE
        self.rootDbPath = str(_validate_owned_path(Path(rootDir) / "soxspipe.db", rootDir, "database path"))
        # COPIES OF THE DATABASE KEPT WHENEVER IT IS REBUILT (NEVER PRUNED)
        self.dbBackupsDir = str(_validate_owned_path(Path(rootDir) / "backups", rootDir, "database backups directory"))
        # PATH OF THE BACKUP MADE BY THE MOST RECENT DATABASE REBUILD, IF ANY
        self.databaseBackupPath = None

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

        return

    def prepare(self, refresh=False, report=True, *, _failedToOpen=False):
        """*Prepare the workspace for data reduction by generating all SOF files and reduction scripts.*

        **Key Arguments:**
        - ``refresh`` -- trigger a complete refresh the workspace during preparation (rebuild the database and do a
          complete prepare). The old database is first copied into `backups/` in the workspace root, its
          `quality_control` rows are restored into the rebuilt database, and so are the pass/fail statuses of every
          SOF that still holds the same frames (see `_restore_session_statuses`).
        - ``_failedToOpen`` -- private, set only by `_rebuild_database_that_failed_to_open`: keep a raw copy of the
          database instead of a SQLite snapshot.

        **Raises:**

        - `DatabasePreservationError` when ``refresh`` is set and the database cannot be copied into `backups/`.
          The database is then left in place and nothing is rebuilt.
        """
        self.log.debug("starting the ``prepare`` method")

        backupPath = None
        if refresh:
            # PRESERVE THE DATABASE BEFORE REMOVING IT; RAISES, LEAVING IT IN PLACE, IF IT CANNOT BE PRESERVED
            backupPath = self._preserve_database(failedToOpen=_failedToOpen)
            # THE AUTOMATIC-REBUILD ROUTE HAS ALREADY DELETED `self.conn`
            if getattr(self, "conn", None):
                self.conn.close()
            self.conn = None
            self._remove_database_files()
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
            return

        # MK RAW FRAME DIRECTORY
        if not os.path.exists(self.rawDir):
            os.makedirs(self.rawDir)

        # TEST FOR SQLITE DATABASE - ADD IF MISSING
        self.conn, reset = self._get_or_create_db_connection()

        # MK SESSION DIRECTORY
        if not os.path.exists(self.sessionsDir):
            os.makedirs(self.sessionsDir)

        # A newly copied database contains only the template's base-session
        # objects.  Restore the empty schema objects for every session before
        # synchronising frames, since that synchronisation can consult the
        # current session's SOF map.
        if self.freshRun:
            self._restore_session_database_objects()

        # QUALITY-CONTROL ROWS EXIST ONLY IN THE DATABASE, SO BRING THEM BACK FROM THE BACKUP
        # BEFORE THE QC GUARDRAIL QUERIES BELOW RECOMPUTE PRODUCT STATUS FROM THEM
        if backupPath:
            self._restore_quality_control_history(backupPath)
            # THE PER-SESSION STATUSES ARE RESTORED AFTER THE SOF BUILDS BELOW (SEE `_restore_session_statuses`)

        basename = os.path.basename(self.rootDir)
        print(f"PREPARING THE `{basename}` WORKSPACE FOR DATA-REDUCTION")

        self._sync_raw_frames()
        self._move_misc_files()

        # IF SESSION ID FILE DOES NOT EXIST, CREATE A NEW SESSION
        # OTHERWISE USE CURRENT SESSION
        self._select_session()

        # GET SETTINGS
        settingsPath = self.sessionPath + "/soxspipe.yaml"
        su = tools(
            arguments={
                "<workspaceDirectory>": self.sessionPath,
                "settingsFile": settingsPath,
            },
            docString=False,
            logLevel="WARNING",
            options_first=False,
            projectName="soxspipe",
            createLogger=False,
        )
        arguments, self.settings, replacedLog, dbConn = su.setup()

        self._apply_qc_acceptable_ranges()

        self._flag_files_to_ignore()
        self.build_sof_files(rejoinLateFrames=True)
        self.build_sof_files(rejoinLateFrames=True)

        if backupPath:
            self._restore_session_statuses(backupPath)

        if report:
            self._print_prepare_report(basename)

        self.log.debug("completed the ``prepare`` method")
        return

    def _select_session(self):
        import codecs

        if not os.path.exists(self.sessionIdFile):
            self.sessionId = self.session_create(sessionId="base")
        else:
            with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
                sessionId = _validate_session_id(readFile.read())
                self.sessionPath = str(
                    _validate_owned_path(
                        Path(self.sessionsDir) / sessionId,
                        self.sessionsDir,
                        "session path",
                    )
                )
        return

    def _print_prepare_report(self, basename):
        rawDirStr = self.rawDir.replace("./", "")

        print(f"\nTHE `{basename}` WORKSPACE FOR HAS BEEN PREPARED FOR DATA-REDUCTION\n")
        print("In this workspace you will find:\n")
        print("   - `backups/`: copies of `soxspipe.db` kept each time the database is rebuilt")
        print("   - `misc/`: a lost-and-found archive of non-fits files")
        print("   - `qc/`: nested folders, ordered by date, containing quality-control plots and tables.")
        print(f"   - `{rawDirStr}/`: nested folders, ordered by date, containing raw-frames.")
        print("   - `sessions/`: directory of data-reduction sessions")
        print("   - `sof/`: the set-of-files (sof) files required for each reduction step")
        print("   - `soxspipe.db`: a sqlite database needed by the data-organiser, please do not delete")
        print("   - `reduced/`: nested folders, ordered by date, containing reduced data.\n")

        incompleteSets, blockingSets = self.get_incomplete_sets_report()
        print_incomplete_sets_report(incompleteSets, blockingSets)

        self.conn.close()
        return

    def _apply_qc_acceptable_ranges(self):
        """*apply this session's QC-acceptable-range settings to the QC rows and product statuses*"""
        with closing(self.conn.cursor()) as c:
            for sqlQuery, sqlParams in self._qc_acceptable_range_queries():
                c.execute(sqlQuery, sqlParams)
            self.conn.commit()

    def _qc_acceptable_range_queries(self):
        """*build the `(sql, params)` pairs that apply this session's QC-acceptable-range settings*

        Recipe names and QC-range keys come from the workspace YAML settings,
        so they are bound as `?` parameters rather than interpolated into the
        SQL text.

        **Return:**

        - a list of `(sqlQuery, sqlParams)` tuples, ready for `cursor.execute(sqlQuery, sqlParams)`
        """
        statusColumn = validate_sql_identifier(f"status_{self.sessionId}", "status column")
        queries = [
            (
                f'update product_frames set {statusColumn} = "pass" '  # noqa: S608
                f'where {statusColumn} = "fail" and sof in (select sof_name from quality_control);',
                (),
            ),
            (
                "update quality_control set qc_value_min = null, qc_value_max = null, qc_flag = 'pass';",
                (),
            ),
        ]

        for k, v in self.settings.items():
            if k[:5] != "soxs-":
                continue
            recipe = k
            for a in ["acq", "vis", "nir"]:
                if a in v and "qc-acceptable-ranges" in v[a]:
                    arm = a.upper()
                    for kk, vv in v[a]["qc-acceptable-ranges"].items():
                        qc_name = kk.upper().replace("-", " ")
                        qc_min = vv[0]
                        qc_max = vv[1]
                        queries.append(
                            (
                                "update quality_control set qc_value_min = ?, qc_value_max = ? "
                                'where soxspipe_recipe = ? and sof_name like ? and qc_name = ?  and qc_order = "-1";',
                                (qc_min, qc_max, recipe, f"%{arm}%", qc_name),
                            )
                        )
            if "qc-acceptable-ranges" in v:
                for kk, vv in v["qc-acceptable-ranges"].items():
                    qc_name = kk.upper().replace("-", " ")
                    qc_min = vv[0]
                    qc_max = vv[1]
                    queries.append(
                        (
                            "update quality_control set qc_value_min = ?, qc_value_max = ? "
                            'where soxspipe_recipe = ? and qc_name = ?  and qc_order = "-1";',
                            (qc_min, qc_max, recipe, qc_name),
                        )
                    )

        queries += [
            (
                'update quality_control set qc_flag = "pass" '
                "where CAST(qc_value as float) < qc_value_max and CAST(qc_value as float) > qc_value_min "
                'and qc_flag != "pass" and qc_order = "-1";',
                (),
            ),
            (
                'update quality_control set qc_flag = "fail" '
                "where (CAST(qc_value as float) > qc_value_max or CAST(qc_value as float) < qc_value_min) "
                'and qc_flag != "fail" and qc_order = "-1";',
                (),
            ),
            (
                'update product_frames set status_base = "fail" '
                'where sof in (select sof_name from quality_control where qc_flag = "fail");',
                (),
            ),
        ]
        return queries

    def list_obs(self):
        """*list all observation names and IDs in the current workspace*"""
        import pandas as pd

        self.log.debug("starting the ``list_obs`` method")

        import pandas as pd

        query = 'select distinct "eso seq arm", "night start date", "night start mjd", "eso obs name","eso obs id" from raw_frame_sets where complete = 1 and recipe in ("nod_obj","stare_obj","offset_obj") and `eso dpr type` like "%OBJECT%" order by "mjd-obs" asc'
        obsDf = pd.read_sql(query, con=self.conn)

        if len(obsDf.index):
            from tabulate import tabulate

            print(f"THE CURRENT WORKSPACE CONTAINS {len(obsDf.index)} SCIENCE OBSERVATION BLOCKS:")
            print(tabulate(obsDf, headers="keys", tablefmt="pretty", showindex=False))

        self.log.debug("completed the ``list_obs`` method")
        return obsDf

    def list_sofs(self):
        """*list all science object SOF files in the current workspace*"""
        import pandas as pd

        self.log.debug("starting the ``list_sofs`` method")

        query = 'select distinct "eso seq arm", upper(replace("recipe", "_obj", "")) as "recipe", "night start date", round("mjd-obs",1) as "mjd-obs", "eso obs name","eso obs id", "slit", "sof" from raw_frame_sets where complete = 1 and recipe in ("nod_obj","stare_obj","offset_obj") and `eso dpr type` like "%OBJECT%" order by "mjd-obs" asc'
        sofDf = pd.read_sql(query, con=self.conn)

        if len(sofDf.index):
            from tabulate import tabulate

            print(f"# THE CURRENT WORKSPACE CONTAINS {len(sofDf.index)} SCIENCE SOF FILES:\n")
            print(tabulate(sofDf, headers="keys", tablefmt="pretty", showindex=False))
            print()

        self.log.debug("completed the ``list_sofs`` method")
        return sofDf

    def list_raw(self, sofFile):
        """*list the all the raw frames associated with a given science object SOF file*"""
        import pandas as pd

        self.log.debug("starting the ``list_raw`` method")

        sqlQuery = "select sof from product_frames where sof = :sofFile and complete = 1"

        # `sofFile` IS ALREADY BOUND ABOVE AS A NAMED PARAMETER. THE LOOP AND
        # THE LINE BELOW ONLY RE-WRAP THE `sqlQuery` TEXT AROUND ITSELF --
        # NO NEW EXTERNAL VALUE EVER RE-ENTERS THE STRING -- SO NEITHER LINE
        # IS AN INJECTION VECTOR, DESPITE THE F-STRING/SQL-KEYWORD SHAPE.
        for _ in range(4):  # Recursively query up to 5 times
            sqlQuery = f"SELECT distinct sof FROM product_frames WHERE file IN (SELECT file FROM sof_map_base WHERE sof in ({sqlQuery})) or sof in ({sqlQuery})"  # noqa: S608

        sqlQuery = (
            f"SELECT * from sof_map WHERE sof in ({sqlQuery}) and filepath "  # noqa: S608
            "like '%./raw/%' order by sof"
        )

        table = pd.read_sql(sqlQuery, con=self.conn, params={"sofFile": sofFile})

        filepaths = table["filepath"].tolist()

        self.log.debug("completed the ``list_raw`` method")
        return list(set(filepaths)), table

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

        import re
        import shutil

        import pandas as pd

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
                    matchedFiles = pd.merge(
                        rawFrames,
                        knownRawFrames,
                        on=["file", "eso dpr tech"],
                        how="inner",
                    )
                    if len(matchedFiles.index):
                        for file in matchedFiles["file"]:
                            try:
                                os.remove(file)
                            except OSError as e:
                                self.log.debug(f"_sync_raw_frames: `os.remove(file)` failed, continuing: {e}")
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
                    c = self.conn.cursor()
                    placeholders = ", ".join("?" for _ in databaseDeletes)
                    sqlQuery = f"delete from raw_frames where filepath in ({placeholders});"  # noqa: S608
                    c.execute(sqlQuery, databaseDeletes)
                    c.close()

        if not skipSqlSync:
            self._sync_sql_table_to_directory(self.rawDir, "raw_frames", recursive=False)

        self.log.debug("completed the ``_sync_raw_frames`` method")
        return

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

        import pandas as pd
        from ccdproc import ImageFileCollection

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
                "The directory contains data from a mix of instruments. Please only provide data from either SOXS or XSH"
            )
            raise AssertionError
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
            except Exception:
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

        import shutil

        if tableName != "raw_frames":
            raise ValueError("Only the raw_frames table can be synchronized")
        # THE EQUALITY CHECK ABOVE ALREADY RESTRICTS `tableName` TO THE SINGLE
        # LITERAL "raw_frames"; THIS ADDS THE SAME GRAMMAR CHECK USED FOR EVERY
        # OTHER INTERPOLATED IDENTIFIER IN THIS FILE, FOR CONSISTENCY.
        tableName = validate_sql_identifier(tableName, "table name")

        # GENERATE A LIST OF FITS FILE PATHS IN RAW DIR
        from fundamentals.files import recursive_directory_listing

        fitsPaths = recursive_directory_listing(
            log=self.log,
            baseFolderPath=directory,
            whatToList="files",  # all | files | dirs
        )

        c = self.conn.cursor()

        sqlQuery = f"select filepath from {tableName};"  # noqa: S608
        c.execute(sqlQuery)

        # NORMALIZE DATABASE PATHS BEFORE COMPARING THEM WITH THE ABSOLUTE LISTING.
        dbFiles = [r[0].replace("//", "/") for r in c.fetchall()]
        normalizedDbFiles = {
            filePath: str(
                _validate_owned_path(
                    Path(filePath) if os.path.isabs(filePath) else Path(self.rootDir) / filePath,
                    self.rootDir,
                    "database filepath",
                )
            )
            for filePath in dbFiles
        }
        absoluteDbFiles = set(normalizedDbFiles.values())

        # DELETED FILES
        filesNotInDB = list(set(fitsPaths) - absoluteDbFiles)
        filesNotInFS = [
            filePath for filePath, normalizedPath in normalizedDbFiles.items() if normalizedPath not in fitsPaths
        ]
        if len(filesNotInFS):
            placeholders = ", ".join("?" for _ in filesNotInFS)
            sqlQuery = f"delete from {tableName} where filepath in ({placeholders});"  # noqa: S608
            c.execute(sqlQuery, filesNotInFS)
            sessionId = _validate_session_id(self.sessionId)
            sofTableName = validate_sql_identifier(f"sof_map_{sessionId}", "sof map table name")
            sqlQuery = f"delete from {sofTableName} where sof in (select sof from {sofTableName} where filepath in ({placeholders}));"  # noqa: S608
            c.execute(sqlQuery, filesNotInFS)

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
                        os.symlink(
                            os.path.realpath(f),
                            os.path.abspath(self.rootDir) + "/" + basename,
                        )
                else:
                    shutil.move(self.rootDir + "/" + f, self.rootDir)
            self._sync_raw_frames(skipSqlSync=True)

        c.close()

        self.log.debug("completed the ``_sync_sql_table_to_directory`` method")
        return

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

        import numpy as np
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
            ((filteredFrames["slit"].str.contains("MULT")) & (filteredFrames["slitmask"] == "--")),
            "slitmask",
        ] = "MPH"
        filteredFrames.loc[
            ((filteredFrames["slit"].str.contains("PINHOLE")) & (filteredFrames["slitmask"] == "--")),
            "slitmask",
        ] = "PH"
        filteredFrames.loc[
            ((filteredFrames["slit"].str.contains("SLIT")) & (filteredFrames["slitmask"] == "--")),
            "slitmask",
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
                    ((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] != "--")),
                    "lamp",
                ] += lamp
                filteredFrames.loc[
                    ((filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99) & (filteredFrames["lamp"] == "--")),
                    "lamp",
                ] = lamp
                filteredFrames.loc[
                    ((filteredFrames[self.kw("DPR_TYPE").lower()] == "DOME,FLAT") & (filteredFrames["lamp"] == "--")),
                    "lamp",
                ] = "DOME"

            else:
                filteredFrames.loc[(filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99), "lamp"] = (
                    filteredFrames.loc[
                        (filteredFrames[self.kw(f"LAMP{i}").lower()] != -99.99),
                        self.kw(f"LAMP{i}").lower(),
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
                    tablefmt="pretty",
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
        allowlistExtensions = [".db", ".yaml", ".log", ".sh", ".py"]
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
        return

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
        self.sessionPath = str(_validate_owned_path(self.sessionPath, self.sessionsDir, "session path"))
        self.sofDir = str(_validate_owned_path(self.sofDir, self.sessionPath, "SOF directory"))
        if not os.path.exists(self.sofDir):
            os.makedirs(self.sofDir)

        sofMapTableName = validate_sql_identifier(f"sof_map_{self.sessionId}", "sof map table name")
        df = pd.read_sql_query(f"select * from {sofMapTableName} where complete = 1;", conn)  # noqa: S608

        # GROUP RESULTS
        for name, group in df.groupby("sof"):
            if not isinstance(name, str) or Path(name).name != name:
                raise _UnsafePathError("SOF filename must be a filename without directory components")
            sofPath = str(
                _validate_owned_path(
                    Path(self.sofDir) / name,
                    self.sofDir,
                    "SOF path",
                )
            )
            if os.path.exists(sofPath):
                continue
            myFile = open(sofPath, "w")
            content = tabulate(group[["filepath", "tag"]], tablefmt="plain", showindex=False)

            myFile.write(content)
            myFile.close()

        self.log.debug("completed the ``_write_sof_files`` method")
        return

    def session_create(self, sessionId=False):
        """*create a data-reduction session with accompanying settings file and required directories*

        **Key Arguments:**

        - ``sessionId`` -- optionally provide a sessionId (A-Z, a-z, 0-9 and/or _ allowed, 16 character limit)

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

        if sessionId:
            sessionId = _validate_session_id(sessionId)

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

        if not sessionId:
            # CREATE SESSION ID FROM TIME STAMP
            from datetime import datetime

            now = datetime.now()
            sessionId = now.strftime("%Y%m%dt%H%M%S")
            sessionId = _validate_session_id(sessionId)

        self.sessionId = sessionId

        # MAKE THE SESSION DIRECTORY
        self.sessionPath = str(
            _validate_owned_path(
                Path(self.sessionsDir) / sessionId,
                self.sessionsDir,
                "session path",
            )
        )
        if not os.path.exists(self.sessionPath):
            os.makedirs(self.sessionPath)

        # SETUP SESSION SETTINGS AND LOGGING
        testPath = self.sessionPath + "/soxspipe.yaml"
        exists = os.path.exists(testPath)

        if not exists:
            if "shoo" in inst:
                inst = "xsh"
            su = tools(
                arguments={
                    "<workspaceDirectory>": self.sessionPath,
                    "init": True,
                    "settingsFile": None,
                },
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
            except OSError as e:
                self.log.debug(f"session_create: `os.symlink(self.vltReduced, dest)` failed, continuing: {e}")

        folders = ["sof", "qc", "reduced"]
        for f in folders:
            if not os.path.exists(self.sessionPath + f"/{f}"):
                os.makedirs(self.sessionPath + f"/{f}")

        self._ensure_session_database_objects(sessionId)

        self._write_sof_files()

        # _reduce_all.sh HAS BEEN REPLACED WITH THE `SOXSPIPE REDUCE` COMMAND
        # self._write_reduction_shell_scripts()

        self._symlink_session_assets_to_workspace_root()

        # WRITE THE SESSION ID FILE
        import codecs

        sessionIdFile = _validate_owned_path(self.sessionIdFile, self.sessionsDir, "session ID path")
        with codecs.open(sessionIdFile, encoding="utf-8", mode="w") as writeFile:
            writeFile.write(sessionId)

        message = f"A new data-reduction session has been created with sessionId '{sessionId}'"
        try:
            self.log.print(message)
        except (AttributeError, OSError, ValueError):
            print(message)
        self.log.debug("completed the ``session_create`` method")

        return sessionId

    def _ensure_session_database_objects(self, sessionId):
        """Create the empty database objects required by one session."""
        sessionId = _validate_session_id(sessionId)
        statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
        sofMapTableName = validate_sql_identifier(f"sof_map_{sessionId}", "sof map table name")

        conn, reset = self._get_or_create_db_connection()
        c = conn.cursor()
        sqlQuery = f"ALTER TABLE product_frames ADD {statusColumn} TEXT;"  # noqa: S608
        try:
            c.execute(sqlQuery)
        except sqlite3.OperationalError as e:
            self.log.debug(f"_ensure_session_database_objects: status column already exists: {e}")

        # DUPLICATE TEH SOF_MAP TABLE
        sqlQuery = "SELECT sql FROM sqlite_master WHERE type='table' AND name='z_sof_map'"
        c.execute(sqlQuery)
        sqlQuery = c.fetchall()[0][0]
        sqlQuery = sqlQuery.replace("z_sof_map", sofMapTableName)
        try:
            c.execute(sqlQuery)
        except sqlite3.OperationalError as e:
            self.log.debug(f"_ensure_session_database_objects: SOF map table already exists: {e}")

        sqlQueries = [
            "DROP VIEW IF EXISTS sof_map;",
            f"CREATE VIEW sof_map as select * from {sofMapTableName};",  # noqa: S608
        ]
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)

        c.close()

    def _session_ids_on_disk(self):
        """*list the valid session directories in the sessions directory*

        **Return:**

        - ``sessionIds`` -- the session IDs, in directory order; entries that are not valid sessions are skipped
        """
        sessionIds = []
        for entry in Path(self.sessionsDir).iterdir():
            if not entry.is_dir():
                self.log.debug(f"Skipping non-session entry while rebuilding database: {entry.name}")
                continue
            try:
                sessionIds.append(_validate_session_id(entry.name))
            except _UnsafePathError:
                self.log.debug(f"Skipping invalid session directory while rebuilding database: {entry.name}")
        return sessionIds

    def _restore_session_database_objects(self):
        """Restore schema objects for session directories after a database rebuild."""
        sessionIds = self._session_ids_on_disk()

        currentSession = None
        if Path(self.sessionIdFile).is_file():
            currentSession = _validate_session_id(Path(self.sessionIdFile).read_text(encoding="utf-8"))

        # The helper updates the shared view, so process the current session
        # last to leave the view pointing at the correct table.
        for sessionId in sorted(sessionIds, key=lambda value: value == currentSession):
            self._ensure_session_database_objects(sessionId)

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
        self.sessionIdFile = str(
            _validate_owned_path(
                Path(self.sessionsDir) / ".sessionid",
                self.sessionsDir,
                "session ID path",
            )
        )
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
            currentSession = _validate_session_id(readFile.read())

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

        sessionId = _validate_session_id(sessionId)

        currentSession, allSessions = self.session_list(silent=True)

        if sessionId == currentSession:
            print(f"Session '{sessionId}' is already in use.")
            return
        if sessionId in allSessions:
            sessionPath = _validate_owned_path(
                Path(self.sessionsDir) / sessionId,
                self.sessionsDir,
                "session path",
            )
            # POINT THE SHARED SOF_MAP VIEW AT THIS SESSION BEFORE MAKING IT ACTIVE
            self._ensure_session_database_objects(sessionId)
            # WRITE THE SESSION ID FILE
            sessionIdFile = _validate_owned_path(self.sessionIdFile, self.sessionsDir, "session ID path")
            with codecs.open(sessionIdFile, encoding="utf-8", mode="w") as writeFile:
                writeFile.write(sessionId)
        else:
            print(f"There is no session with the ID '{sessionId}'. List existing sessions with `soxspipe session ls`.")
            return

        self.sessionPath = str(sessionPath)
        self._symlink_session_assets_to_workspace_root()
        print(f"Session successfully switched to '{sessionId}'.")

        self.log.debug("completed the ``session_switch`` method")
        return

    def _symlink_session_assets_to_workspace_root(self):
        """*symlink session QC, product, SOF directories, database and scripts to workspace root*

        **Key Arguments:**
            # -

        **Return:**

        - None

        """
        self.log.debug("starting the ``_symlink_session_assets_to_workspace_root`` method")

        import os

        # SYMLINK FILES AND FOLDERS
        toLink = ["reduced", "qc", "soxspipe.yaml", "sof", "soxspipe.log"]
        for l in toLink:
            dest = self.rootDir + f"/{l}"
            src = self.sessionPath + f"/{l}"
            try:
                os.symlink(src, dest)
            except OSError as e:
                self.log.debug(f"_symlink_session_assets_to_workspace_root: `os.symlink(sr...` failed, continuing: {e}")
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
                except OSError as e:
                    self.log.debug(
                        f"_symlink_session_assets_to_workspace_root: `os.symlink(src, dest)` failed, continuing: {e}"
                    )
                    os.unlink(dest)
                    os.symlink(src, dest)

        self.log.debug("completed the ``_symlink_session_assets_to_workspace_root`` method")
        return

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

        import os
        import sys

        if failure is True:
            self.log.print("\nRefeshing SOF files due to recipe failure\n")
        elif failure is False:
            self.log.print("\nRefreshing SOF files, as a previously failed recipe is now passing.\n")
        import codecs

        # IF SESSION ID FILE DOES NOT EXIST, REPORT
        self.sessionIdFile = str(_validate_owned_path(self.sessionIdFile, self.sessionsDir, "session ID path"))
        exists = os.path.exists(self.sessionIdFile)
        if not exists:
            if not silent:
                print("No reduction sessions exist in this workspace yet.")
            return None, None
        with codecs.open(self.sessionIdFile, encoding="utf-8", mode="r") as readFile:
            sessionId = _validate_session_id(readFile.read())
        self.sessionPath = str(
            _validate_owned_path(
                Path(self.sessionsDir) / sessionId,
                self.sessionsDir,
                "session path",
            )
        )
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
        except (AttributeError, sqlite3.ProgrammingError) as e:
            self.log.debug(f"close: `self.conn.close()` failed, continuing: {e}")

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
            with open(advs) as stream:
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
        except OSError as e:
            self.log.debug(f"use_vlt_environment_folders: `os.symlink(vltRaw, self.rawDir)` failed, continuing: {e}")
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
                log=self.log,
                baseFolderPath=self.rawDir,
                whatToList="files",  # all | files | dirs
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

    def _preserve_database(self, failedToOpen=False):
        """*keep a complete, durable copy of the workspace database in the backups directory before it is rebuilt*

        The `quality_control` table exists nowhere but the database, so the database is never deleted until a
        verified copy of it is safe on disk in `backups/` inside the workspace root. Backups are never pruned. The
        caller deletes the original (see `_remove_database_files`) only after this method returns.

        A readable database is copied as a single-file SQLite snapshot named `<stem>_refresh_<UTC stamp>.db`. A
        database that failed to open, or that SQLite reports is not a database or is malformed, is copied byte
        for byte (with any `-wal`/`-shm` sidecar files) as `<stem>_corrupt_<UTC stamp>.db`. Any other failure,
        such as a locked or busy database, refuses the rebuild.

        **Key Arguments:**

        - ``failedToOpen`` -- the database failed to open, so skip the snapshot and keep a raw copy.

        **Return:**

        - ``backupPath`` -- path of the preserved database file, or None when there is no database to preserve.

        **Raises:**

        - `DatabasePreservationError` when the database cannot be preserved. Nothing in the workspace is changed.
        """
        from datetime import UTC, datetime

        if not os.path.exists(self.rootDbPath):
            return None

        dbPath = Path(self.rootDbPath)
        stamp = datetime.now(UTC).strftime("%Y%m%dT%H%M%S%fZ")
        backupPath = None
        try:
            os.makedirs(self.dbBackupsDir, exist_ok=True)
            if not failedToOpen:
                backupPath = self._backup_path(dbPath, "refresh", stamp)
                try:
                    self._snapshot_database(backupPath)
                except sqlite3.DatabaseError as error:
                    if not _is_unreadable_database_error(error):
                        raise
                    self.log.warning(f"cannot snapshot `{self.rootDbPath}` ({error}); keeping a raw copy instead")
                    failedToOpen = True
            if failedToOpen:
                backupPath = self._backup_path(dbPath, _RAW_BACKUP_LABEL, stamp)
                self._copy_database_files(backupPath)
            _fsync_path(self.dbBackupsDir)
        except (OSError, sqlite3.Error, _UnsafePathError) as error:
            message = (
                f"Could not preserve the database `{self.rootDbPath}` in `{self.dbBackupsDir}` ({error}). "
                "It has been left in place and the workspace has not been rebuilt. Fix the cause (for example "
                "free disk space, or stop any other soxspipe process using the workspace), then run "
                "`soxspipe prep --refresh` again."
            )
            self.log.error(message)
            raise DatabasePreservationError(message) from error

        self.databaseBackupPath = str(backupPath)
        return self.databaseBackupPath

    def _backup_path(self, dbPath, label, stamp):
        """*return the validated path of a new database backup file inside the backups directory*"""
        return _validate_owned_path(
            Path(self.dbBackupsDir) / f"{dbPath.stem}_{label}_{stamp}{dbPath.suffix}",
            self.dbBackupsDir,
            "database backup path",
        )

    def _snapshot_database(self, backupPath):
        """*write a consistent, single-file SQLite snapshot of the workspace database to `backupPath`*

        Uses the SQLite online-backup API from a read-only source, so the copy is transactionally consistent
        even while other connections hold the database open, and includes any content still in the `-wal` file.
        The snapshot is switched to `DELETE` journal mode so it is one self-contained file, checked with
        `quick_check`, and flushed to disk. A partial or unverifiable snapshot is removed and the error re-raised.
        """
        if os.path.exists(backupPath):
            raise FileExistsError(f"database backup `{backupPath}` already exists")
        try:
            # READ-ONLY, SO A FAILED SNAPSHOT CANNOT ALTER THE DATABASE OR ITS SIDECARS BEFORE THE RAW COPY
            source = sqlite3.connect(Path(self.rootDbPath).resolve().as_uri() + "?mode=ro", uri=True)
            try:
                target = sqlite3.connect(backupPath)
                try:
                    source.backup(target)
                    target.execute("PRAGMA journal_mode=DELETE;")
                    check = target.execute("PRAGMA quick_check;").fetchall()
                finally:
                    target.close()
            finally:
                source.close()
            if check != [("ok",)]:
                raise _UnreadableDatabaseError(f"database snapshot `{backupPath}` failed its quick check: {check}")
            _fsync_path(backupPath)
        except (OSError, sqlite3.Error):
            self._remove_partial_backup(backupPath)
            raise

    def _copy_database_files(self, backupPath):
        """*copy the raw database file and any `-wal`/`-shm` sidecar files to `backupPath`, byte for byte*

        Used for a database that will not open, so no SQLite call is made on it. Each copy must match its
        source's size and is flushed to disk. Nothing is removed from the workspace here; on a partial failure
        the copies already made are removed and the error re-raised.
        """
        if os.path.exists(backupPath):
            raise FileExistsError(f"database backup `{backupPath}` already exists")
        try:
            for suffix in ("", "-wal", "-shm"):
                source = self.rootDbPath + suffix
                if not os.path.exists(source):
                    continue
                target = f"{backupPath}{suffix}"
                shutil.copy2(source, target)
                sourceSize = os.path.getsize(source)
                targetSize = os.path.getsize(target)
                if sourceSize != targetSize:
                    raise OSError(f"copy of `{source}` has size {targetSize}, expected {sourceSize}")
                _fsync_path(target)
        except OSError:
            self._remove_partial_backup(backupPath)
            raise

    def _remove_partial_backup(self, backupPath):
        """*remove what a failed backup attempt left behind in the backups directory*"""
        for suffix in ("", "-wal", "-shm", "-journal"):
            try:
                os.remove(f"{backupPath}{suffix}")
            except FileNotFoundError:
                continue
            except OSError as e:
                self.log.debug(f"_remove_partial_backup: removing `{backupPath}{suffix}` failed, continuing: {e}")

    def _remove_database_files(self):
        """*delete the workspace database and any `-shm`/`-wal` sidecar files*

        Sidecar files are removed even when the main file is missing, so a stale `-wal` can never attach itself
        to the fresh database copied from the template.
        """
        if os.path.exists(self.rootDbPath):
            os.remove(self.rootDbPath)
            print("The existing database has been removed to allow a complete refresh of the workspace.")
        for suffix in ("-shm", "-wal"):
            try:
                os.remove(self.rootDbPath + suffix)
            except FileNotFoundError:
                continue
            except OSError as e:
                self.log.debug(f"_remove_database_files: removing `{self.rootDbPath + suffix}` failed, continuing: {e}")

    def _read_preserved_quality_control(self, backupPath):
        """*read the `quality_control` rows from a preserved copy of the database*

        The copy is opened read-only, so reading it never changes the preserved file, and it must pass
        `PRAGMA quick_check` before any row is read.

        **Key Arguments:**

        - ``backupPath`` -- path of the preserved database file

        **Return:**

        - ``qcRows`` -- pandas DataFrame of every `quality_control` row in the preserved database

        **Raises:**

        - `sqlite3.DatabaseError` when the preserved file fails its `quick_check`
        """
        import pandas as pd

        source = sqlite3.connect(Path(backupPath).resolve().as_uri() + "?mode=ro", uri=True)
        try:
            check = source.execute("PRAGMA quick_check;").fetchall()
            if check != [("ok",)]:
                raise sqlite3.DatabaseError(f"the preserved database failed its quick check: {check}")
            return pd.read_sql_query("select * from quality_control;", source)
        finally:
            source.close()

    def _write_quality_control_rows(self, qcRows):
        """*append preserved `quality_control` rows to the current database, skipping rows already present*

        A row is skipped when a row with the same `(qc_name, sof_name, qc_order)` key already exists, with
        NULL matching NULL (the table's `UNIQUE ... ON CONFLICT IGNORE` constraint treats NULLs as distinct,
        so it alone would duplicate rows whose `qc_order` is NULL). Only the columns the current schema knows
        are written, so a backup made by an older schema still restores; any column dropped this way is named
        in a log warning.

        The values pass unchanged through `_dataframe_to_sqlite`: its `"--"`/`-99.99` to NULL substitution
        cannot match, because every QC writer already applied it before the rows were stored, and the text
        columns read back as strings.

        **Key Arguments:**

        - ``qcRows`` -- pandas DataFrame from `_read_preserved_quality_control`

        **Return:**

        - ``insertedCount`` -- the number of rows actually added to the table
        """
        if qcRows.empty:
            return 0
        currentColumns = {row[1] for row in self.conn.execute("PRAGMA table_info(quality_control);")}
        keptColumns = [column for column in qcRows.columns if column in currentColumns]
        droppedColumns = sorted(set(qcRows.columns) - currentColumns)
        if droppedColumns:
            self.log.warning(f"quality-control columns not in the current schema were not restored: {droppedColumns}")

        keyColumns = ["qc_name", "sof_name", "qc_order"]
        existingKeys = set(self.conn.execute("select qc_name, sof_name, qc_order from quality_control;"))
        rowKeys = qcRows.reindex(columns=keyColumns).astype(object)
        rowKeys = rowKeys.where(rowKeys.notna(), None)
        isNew = [key not in existingKeys for key in rowKeys.itertuples(index=False, name=None)]
        newRows = qcRows.loc[isNew, keptColumns]
        if newRows.empty:
            return 0

        countQuery = "select count(*) from quality_control;"
        before = self.conn.execute(countQuery).fetchone()[0]
        self._dataframe_to_sqlite(newRows, "quality_control", replace=False)
        return self.conn.execute(countQuery).fetchone()[0] - before

    @staticmethod
    def _describe_preserved_database(backupPath):
        """*return the phrase that names a preserved database in messages to the user*

        **Key Arguments:**

        - ``backupPath`` -- path of the preserved database file

        **Return:**

        - ``source`` -- the phrase, marking a raw copy of a database that failed to open
        """
        source = f"the preserved database `{backupPath}`"
        if f"_{_RAW_BACKUP_LABEL}_" in Path(backupPath).name:
            source += " (a raw copy of a database that failed to open)"
        return source

    def _restore_quality_control_history(self, backupPath):
        """*restore `quality_control` rows from a preserved database into the rebuilt one, best effort*

        Never raises for a source it cannot read. When the preserved file cannot be read or fails its
        `quick_check` (for example, a raw copy of a database that failed to open), a warning naming the
        preserved path is printed and the rebuild continues without QC history.

        **Key Arguments:**

        - ``backupPath`` -- path of the preserved database file

        **Return:**

        - ``restoredCount`` -- the number of rows added to the rebuilt database (0 on failure)
        """
        import pandas as pd

        source = self._describe_preserved_database(backupPath)
        try:
            qcRows = self._read_preserved_quality_control(backupPath)
            restoredCount = self._write_quality_control_rows(qcRows)
        except (sqlite3.Error, pd.errors.DatabaseError, OSError, ValueError, TypeError, KeyError) as error:
            self.log.warning(f"could not restore quality-control history from `{backupPath}`: {error}")
            self.log.debug(traceback.format_exc())
            print(
                f"WARNING: could not restore quality-control history from {source} ({error}). "
                "The preserved file has been kept."
            )
            return 0

        print(f"Restored {restoredCount} of {len(qcRows)} quality-control rows from {source}.")
        return restoredCount

    def _restore_session_statuses(self, backupPath):
        """*restore the per-session `status_<id>` columns from a preserved database and print a summary per session*

        **Key Arguments:**

        - ``backupPath`` -- path of the preserved database file

        **Return:**

        - ``counts`` -- session ID to the `(restored, changed, dropped)` counts, or an empty dict when the preserved
          file cannot be read. Never raises for the preserved file: a warning naming it is printed and the rebuild
          continues with empty statuses.
        """
        sessionIds = self._session_ids_on_disk()
        try:
            snapshots = read_session_snapshots(backupPath, sessionIds, self.log)
        except (sqlite3.Error, OSError, ValueError) as error:
            self.log.warning(f"could not restore session statuses from `{backupPath}`: {error}")
            self.log.debug(traceback.format_exc())
            print(
                "WARNING: could not restore session statuses from "
                f"{self._describe_preserved_database(backupPath)} ({error}). "
                "Status columns are left empty; the preserved file has been kept."
            )
            return {}
        counts = self._converge_session_statuses(snapshots)
        for sessionId in sorted(counts):
            print(format_session_summary(sessionId, *counts[sessionId], backupPath))
        return counts

    def _converge_session_statuses(self, snapshots):
        """*restore statuses until no restored failure changes what the remaining SOFs are built from*

        The SOF maps are built with every status empty, so a SOF downstream of a failed calibration may hold a
        different calibration than it did before the rebuild. Each pass writes the statuses of the SOFs that now
        hold the frames they held before, re-applies the QC ranges, and, when that adds a failure, rebuilds the SOF
        maps so downstream SOFs are re-matched away from the failed product. SOFs still differing afterwards are
        counted as changed.

        **Key Arguments:**

        - ``snapshots`` -- session ID to `session_status_snapshot`

        **Return:**

        - ``counts`` -- session ID to the `(restored, changed, dropped)` counts
        """
        candidates = {sessionId: frozenset(snapshot.statuses) for sessionId, snapshot in snapshots.items()}
        counts = dict.fromkeys(snapshots, (0, 0, 0))
        for _ in range(self._STATUS_RESTORE_MAX_PASSES):
            failuresBefore = self._current_session_failures()
            results = self._restore_status_pass(snapshots, candidates)
            for sessionId, result in results.items():
                self._write_session_statuses(sessionId, result.toRestore)
            counts = self._add_pass_counts(counts, results)
            candidates = {sessionId: result.pending for sessionId, result in results.items()}
            # A RESTORED FAIL CAN BE OVERTURNED BY QC ROWS THAT NOW PASS, AND QC FAILURES MARK THEIR PRODUCTS AGAIN
            self._apply_qc_acceptable_ranges()
            hasNewFailures = not self._current_session_failures() <= failuresBefore
            if hasNewFailures:
                self.build_sof_files()
            if not hasNewFailures or not any(candidates.values()):
                break
        else:
            self.log.debug("_converge_session_statuses: stopped at the pass limit with SOFs still differing")
        return {
            sessionId: (restored, changed + len(candidates[sessionId]), dropped)
            for sessionId, (restored, changed, dropped) in counts.items()
        }

    @staticmethod
    def _add_pass_counts(counts, results):
        """*add one pass's results to the running `(restored, changed, dropped)` counts*

        **Key Arguments:**

        - ``counts`` -- session ID to the counts so far
        - ``results`` -- session ID to the `status_pass_result` of the latest pass

        **Return:**

        - ``counts`` -- new session ID to counts mapping; the inputs are not changed
        """
        return {
            sessionId: (
                counts[sessionId][0] + len(result.toRestore),
                counts[sessionId][1] + len(result.changed),
                counts[sessionId][2] + len(result.dropped),
            )
            for sessionId, result in results.items()
        }

    def _current_session_failures(self):
        """*list the SOFs the current session has marked `fail`*

        **Return:**

        - ``failures`` -- frozenset of SOF names
        """
        statusColumn = validate_sql_identifier(f"status_{self.sessionId}", "status column")
        # COLUMN NAME CANNOT BE BOUND; COMPOSED FROM A VALIDATED SESSION ID AND CHECKED BY validate_sql_identifier
        sqlQuery = f"SELECT DISTINCT sof FROM product_frames WHERE {statusColumn} = 'fail';"  # noqa: S608
        return frozenset(row[0] for row in self.conn.execute(sqlQuery))

    def _restore_status_pass(self, snapshots, candidates):
        """*classify each session's candidate SOFs against the rebuilt database*

        The frames are compared with the current session's `sof_map_<id>` table, not the shared `sof_map` view.

        **Key Arguments:**

        - ``snapshots`` -- session ID to `session_status_snapshot`
        - ``candidates`` -- session ID to the SOF names still to classify

        **Return:**

        - ``results`` -- session ID to `status_pass_result`
        """
        rebuiltSofs = frozenset(row[0] for row in self.conn.execute("SELECT DISTINCT sof FROM product_frames;"))
        sofMapTable = validate_sql_identifier(f"sof_map_{self.sessionId}", "sof map table name")
        rebuiltFrameSets = read_frame_sets(self.conn, sofMapTable)
        return {
            sessionId: classify_session_sofs(snapshot, candidates[sessionId], rebuiltSofs, rebuiltFrameSets)
            for sessionId, snapshot in snapshots.items()
        }

    def _write_session_statuses(self, sessionId, sofStatuses):
        """*write SOF statuses into one session's `status_<id>` column*

        **Key Arguments:**

        - ``sessionId`` -- the session whose column is written
        - ``sofStatuses`` -- SOF name to ``pass`` or ``fail``
        """
        statusColumn = validate_sql_identifier(f"status_{_validate_session_id(sessionId)}", "status column")
        # COLUMN NAME CANNOT BE BOUND; COMPOSED FROM A VALIDATED SESSION ID AND CHECKED BY validate_sql_identifier;
        # SOF VALUES ARE BOUND
        sqlQuery = f"UPDATE product_frames SET {statusColumn} = ? WHERE sof = ?;"  # noqa: S608
        self.conn.executemany(sqlQuery, [(status, sof) for sof, status in sofStatuses.items()])
        self.conn.commit()

    def _rebuild_database_that_failed_to_open(self, error):
        """*rebuild the workspace database after it repeatedly failed to open or pass its integrity check*

        A database that SQLite reports is not a database or is malformed is kept as a raw copy in `backups/`.
        Any other failure goes through the normal snapshot-first preservation. A locked or busy database is
        never rebuilt, because another process may still be writing to it.

        **Key Arguments:**

        - ``error`` -- the exception from the last failed attempt to open and check the database

        **Raises:**

        - `DatabasePreservationError` when the database is locked or busy, or could not be copied into
          `backups/`. It is then left in place, never deleted, and there is no usable database to continue with.
        """
        if _is_locked_database_error(error):
            message = (
                f"The database `{self.rootDbPath}` is locked or busy ({error}), so it has not been rebuilt. "
                "Stop any other soxspipe process using this workspace, then try again."
            )
            self.log.error(message)
            raise DatabasePreservationError(message) from error
        self.prepare(refresh=True, _failedToOpen=_is_unreadable_database_error(error))

    def _get_or_create_db_connection(self):
        """Private method to get or create the SQLite database connection, copying the template if missing."""
        import shutil
        import sqlite3 as sql
        import time

        reset = False

        self.rootDbPath = str(_validate_owned_path(self.rootDbPath, self.rootDir, "database path"))

        conn = None
        i = 0

        tries = self._DB_OPEN_ATTEMPTS

        while i < tries:
            if not conn:
                try:
                    if self.conn:
                        conn = self.conn
                except AttributeError as e:
                    self.log.debug(f"_get_or_create_db_connection: `if self.conn: conn = s...` failed, continuing: {e}")

            if not conn:
                try:
                    with open(self.rootDbPath):
                        pass
                    self.freshRun = False
                except OSError:
                    self.freshRun = True
                    emptyDb = os.path.dirname(os.path.dirname(__file__)) + "/resources/soxspipe.db"
                    shutil.copyfile(emptyDb, self.rootDbPath)
                conn = sql.connect(
                    self.rootDbPath,
                    timeout=self._DB_BUSY_TIMEOUT_SECONDS,
                    autocommit=True,
                    check_same_thread=False,
                )
            c = conn.cursor()

            try:
                c.execute("PRAGMA integrity_check;")
                integrityCheckRows = c.fetchall()
                if integrityCheckRows != [("ok",)]:
                    raise sql.DatabaseError(f"database integrity check failed: {integrityCheckRows}")
                c.execute("PRAGMA busy_timeout = 100000")
                c.execute("PRAGMA synchronous = OFF")

                i = tries + 1
            except Exception as error:
                # DATABASE IS BROKEN, REPLACE WITH EMPTY ONE
                lastError = error
                i += 1
                c.close()
                conn.close()
                try:
                    del conn
                    del self.conn
                except (AttributeError, NameError, UnboundLocalError) as e:
                    self.log.debug(f"_get_or_create_db_connection: `del conn` failed, continuing: {e}")

                # ANOTHER CONNECTION HOLDS THE DATABASE, SO REFUSE AT ONCE INSTEAD OF RETRYING
                if _is_locked_database_error(lastError):
                    self._rebuild_database_that_failed_to_open(lastError)

                time.sleep(1)

                if i > tries - 1:
                    self._rebuild_database_that_failed_to_open(lastError)
                    if not reset:
                        reset = True
                conn = sql.connect(
                    self.rootDbPath,
                    timeout=self._DB_BUSY_TIMEOUT_SECONDS,
                    autocommit=True,
                    check_same_thread=False,
                )
                c = conn.cursor()
        c.execute("PRAGMA busy_timeout = 100000")
        c.execute("PRAGMA synchronous = OFF")

        c.close()

        return conn, reset

    # DISPLAY COLUMNS SHARED BY get_incomplete_raw_frames_set AND get_blocking_calibration_sets
    _INCOMPLETE_SET_DISPLAY_COLUMNS = [
        "eso seq arm",
        "mjd-obs",
        "eso dpr tech",
        "eso dpr type",
        "slit",
        "eso obs name",
        "eso obs id",
    ]

    # METHOD TO GROUP INCOMPLETE RAW-FRAME SETS AND WORK OUT WHICH CALIBRATIONS EACH GROUP LACKS, AND WHY
    def _missing_calibrations_by_group(self):
        """*group incomplete science sets and find which calibrations each group lacks, and why*

        Called once per report by `get_incomplete_sets_report`, which both `get_incomplete_raw_frames_set`
        and `get_blocking_calibration_sets` build on, so a report never re-runs this grouping-plus-reason
        pass twice.

        **Return:**

        - ``groups`` -- list of `(displayValues, mergedPositions, reasons, matchContext)` tuples, one per
          distinct combination of `_INCOMPLETE_SET_DISPLAY_COLUMNS`. `mergedPositions` is `{calibration
          type: position in the sof map}`, unioned across every sof sharing that display row. `reasons` is
          `{calibration type: reason string}` from `missing_reasons`, empty when nothing is missing.
          `matchContext` is the group's own `(binning, rospeed, gain)`, taken from its first sof, passed
          to `missing_reasons` so `mbias`/`mflat` are classified against the same mode the group is in.
        """
        import pandas as pd

        query = (
            'select "eso seq arm", round("mjd-obs",1) as "mjd-obs", "eso dpr tech","eso dpr type","slit",'
            '"eso obs name","eso obs id", sof, recipe, binning, rospeed, gain from raw_frame_sets '
            'where complete = 0 and recipe in ("nod_obj","stare_obj","offset_obj")'
        )
        rawSets = pd.read_sql(query, con=self.conn)

        # SOF MAP MAY NOT BE LOADED YET IF THIS data_organiser WAS BUILT WITHOUT CALLING prepare() FIRST
        if not hasattr(self, "sofMapLookup"):
            self._select_instrument()

        # SOF MAP IS STILL ABSENT IF _select_instrument RETURNED EARLY (NO raw_frames TABLE, NO INSTRUMENT YET)
        missingBySof = find_missing_calibrations(self.conn, rawSets, getattr(self, "sofMapLookup", None))

        armIndex = self._INCOMPLETE_SET_DISPLAY_COLUMNS.index("eso seq arm")
        groups = []
        for displayValues, group in rawSets.groupby(self._INCOMPLETE_SET_DISPLAY_COLUMNS, sort=False, dropna=False):
            mergedPositions = {}
            for sof in group["sof"]:
                for calibrationType, position in missingBySof.get(sof, {}).items():
                    mergedPositions[calibrationType] = min(position, mergedPositions.get(calibrationType, position))
            # A DISPLAY GROUP IS ONE OBSERVATION, SO EVERY SOF IN IT SHARES ONE MODE -- THE FIRST ROW STANDS FOR ALL
            matchContext = tuple(
                None if pd.isna(value) else value for value in group[["binning", "rospeed", "gain"]].iloc[0]
            )
            reasons = (
                missing_reasons(self.conn, mergedPositions, displayValues[armIndex], *matchContext)
                if mergedPositions
                else {}
            )
            groups.append((displayValues, mergedPositions, reasons, matchContext))

        return groups

    # METHOD TO BUILD THE INCOMPLETE-SETS AND BLOCKING-SETS TABLES FROM ONE SHARED GROUPING PASS
    def get_incomplete_sets_report(self):
        """*return the incomplete-sets and blocking-sets tables together, from one grouping pass*

        `get_incomplete_raw_frames_set` and `get_blocking_calibration_sets` each call
        `_missing_calibrations_by_group` on their own for standalone use; this method computes it once
        and builds both tables from that one result, which is what `prepare` and `reducer.reduce` use so
        printing the report never re-runs the grouping-plus-reason SQL pass twice.

        **Return:**

        - ``(incompleteSets, blockingSets)`` -- the same two dataframes `get_incomplete_raw_frames_set`
          and `get_blocking_calibration_sets` return
        """
        groups = self._missing_calibrations_by_group()
        return self._incomplete_sets_frame(groups), self._blocking_sets_frame(groups)

    # METHOD TO RETURN ALL THE RAW FRAME SETS THAT ARE NOT COMPLETE FROM THE DATABASE AS A PANDAS TABLE
    def get_incomplete_raw_frames_set(self):
        """*return the science raw-frame sets that cannot be reduced yet, and the calibrations each one lacks*

        **Return:**

        - ``incompleteSets`` -- dataframe with one row per distinct combination of the display columns
          (`eso seq arm`, `mjd-obs`, `eso dpr tech`, `eso dpr type`, `slit`, `eso obs name`, `eso obs id`)
          plus a `missing calibrations` column. That column lists, in sof-map order, the calibrations with no
          passing row in their `cal_<type>` view, each with why it is missing (`failed QC`, `failed run`,
          `not yet reduced`, `not observed` or `no match`), or `unknown` if none can be found. When several
          sofs share one display row, the column lists every calibration that any one of them lacks. The
          frame is empty when every science set is complete.
        """
        return self._incomplete_sets_frame(self._missing_calibrations_by_group())

    def _incomplete_sets_frame(self, groups):
        """*build the `get_incomplete_raw_frames_set` dataframe from an already-computed `groups` list*"""
        import pandas as pd

        rows = [
            [*displayValues, describe_missing(mergedPositions, reasons)]
            for displayValues, mergedPositions, reasons, _ in groups
        ]
        return pd.DataFrame(rows, columns=[*self._INCOMPLETE_SET_DISPLAY_COLUMNS, "missing calibrations"])

    # METHOD TO NAME THE RAW SOFS RESPONSIBLE FOR EACH MISSING CALIBRATION, WHERE THAT IS KNOWN
    def get_blocking_calibration_sets(self):
        """*name the raw sofs blocking each missing calibration currently reported as failed QC, a failed
        run, or not yet reduced*

        For `mbias`/`mflat` this is scoped to the group's own binning/readout speed/gain (see
        `missing_calibrations.calibration_match_clause`); every other type is arm-only, since its
        `cal_<type>` view does not hard-gate on mode either.

        **Return:**

        - ``blockingSets`` -- dataframe with columns `eso seq arm`, `calibration`, `sof`, `recipe`, `reason`,
          `detail` (the first line of the failure message, or `None` when the reason is `not yet reduced`).
          Empty when nothing is currently missing for any of these three reasons.
        """
        return self._blocking_sets_frame(self._missing_calibrations_by_group())

    def _blocking_sets_frame(self, groups):
        """*build the `get_blocking_calibration_sets` dataframe from an already-computed `groups` list*"""
        import pandas as pd

        columns = ["eso seq arm", "calibration", "sof", "recipe", "reason", "detail"]

        pairs = set()
        for displayValues, _, reasons, matchContext in groups:
            arm = displayValues[self._INCOMPLETE_SET_DISPLAY_COLUMNS.index("eso seq arm")]
            for calibrationType, reason in reasons.items():
                if reason in (FAILED_QC, FAILED_RUN, NOT_YET_REDUCED):
                    pairs.add((arm, calibrationType, reason, matchContext))

        rows = []
        for arm, calibrationType, reason, matchContext in sorted(pairs, key=lambda pair: pair[:3]):
            recipes = CALIBRATION_RECIPES.get(calibrationType, (calibrationType,))
            placeholders = ", ".join("?" for _ in recipes)
            matchClause, matchParams = calibration_match_clause(calibrationType, arm, *matchContext)
            if reason in (FAILED_QC, FAILED_RUN):
                # GROUPED: A FAILED RECIPE RUN CAN LEAVE SEVERAL PRODUCT FILES (E.G. A TABLE AND A RESPONSE
                # CURVE) UNDER ONE sof, EACH CARRYING THE SAME error_message -- ONE ROW PER sof, NOT PER FILE
                sqlQuery = (
                    f"SELECT sof, recipe, MIN(error_message) FROM failed_products "  # noqa: S608
                    f'WHERE recipe IN ({placeholders}) AND "eso seq arm" = ?{matchClause} GROUP BY sof, recipe'
                )
                # RECIPE NAMES COME FROM CALIBRATION_RECIPES (A MODULE CONSTANT), NOT USER INPUT; ALL VALUES ARE BOUND
                wantQc = reason == FAILED_QC
                candidates = self.conn.execute(sqlQuery, (*recipes, arm, *matchParams))
                matches = [
                    (sof, recipe, errorMessage)
                    for sof, recipe, errorMessage in candidates
                    if bool((errorMessage or "").startswith(QC_FAILURE_MESSAGE_PREFIX)) == wantQc
                ]
            else:
                sqlQuery = (
                    f"SELECT sof, recipe, NULL FROM raw_frame_sets "  # noqa: S608
                    f'WHERE recipe IN ({placeholders}) AND "eso seq arm" = ? '
                    f"AND complete = 0{matchClause} GROUP BY sof, recipe"
                )
                matches = list(self.conn.execute(sqlQuery, (*recipes, arm, *matchParams)))
            for sof, recipe, errorMessage in matches:
                detail = errorMessage.splitlines()[0] if errorMessage else None
                rows.append([arm, calibrationType, sof, recipe, reason, detail])

        return pd.DataFrame(rows, columns=columns)

    def _select_instrument(self, inst=False):
        """Select the instrument and set related attributes."""
        import yaml

        from soxspipe.commonutils import keyword_lookup

        if inst:
            self.instrument = inst
        else:
            try:
                c = self.conn.cursor()
                sqlQuery = "select instrume from raw_frames where instrume is not null limit 1"
                c.execute(sqlQuery)
                self.instrument = c.fetchall()[0][0]
                c.close()
            except (AttributeError, IndexError, sqlite3.OperationalError) as e:
                self.log.debug(f"_select_instrument: `c = self.conn.cursor()` failed, continuing: {e}")
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
        with open(yamlFilePath) as stream:
            self.sofMapLookup = yaml.safe_load(stream)

    def _flag_files_to_ignore(self):
        """*Flag files to ignore based on settings and reduction order*"""

        # FLAG FILES TO IGNORE BASED ON REDUCTION ORDER
        c = self.conn.cursor()
        placeholders = ", ".join("?" for _ in self.reductionOrder)
        sqlQuery = f"update raw_frames set ignore = 1 WHERE `eso dpr type` not in ({placeholders})"  # noqa: S608
        c.execute(sqlQuery, self.reductionOrder)
        self.conn.commit()

        # FLAG STANDARDS NOT IN STATIC LIBRARY TO IGNORE
        sqlQuery = """UPDATE raw_frames
        SET IGNORE = 1
        WHERE rowid IN (
            SELECT rowid
            FROM (
                SELECT
                    rowid,
                    replace(replace(replace(IFNULL("eso obs name","") || IFNULL("eso obs targ name",""), "_", ""), " ", ""), "-", "") AS matchStr,
                    "eso dpr type" AS dprType
                FROM raw_frames
            )
            WHERE
                dprType LIKE "%STD,FLUX%"
                AND NOT (
                    matchStr LIKE "%GD71%" OR
                    matchStr LIKE "%LTT3218%" OR
                    matchStr LIKE "%GD153%" OR
                    matchStr LIKE "%EG274%" OR
                    matchStr LIKE "%LTT7987%" OR
                    matchStr LIKE "%FEIGE110%" OR
                    matchStr LIKE "%EG21%" OR
                    matchStr LIKE "%CD3017706%" OR
                    matchStr LIKE "%CD3810980%" OR
                    matchStr LIKE "%CD325613%" OR
                    matchStr LIKE "%CPD69177%"
            )
        );"""
        c.execute(sqlQuery)
        self.conn.commit()

        # FLAG SIMULATION FILES TO IGNORE
        sqlQuery = "update raw_frames set ignore = 1 WHERE `eso dpr type` not like '%OBJECT%' and `eso dpr type` not like '%STD%' and simulation = 1"
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

    def _calibration_completeness_select_query(self, recipe, arm, ttype, calType):
        """*build the `(sql, params)` pair that finds product SOFs whose calibration prerequisites are satisfied*

        `recipe`, `arm` and `ttype` come from the instrument's `sof_map.yaml`
        resource, so they are bound as `?` parameters rather than
        interpolated into the SQL text. `calType` entries are calibration
        type names used as *table-name* identifiers (`cal_<type>`), which
        SQLite cannot bind, so each is checked against the safe-identifier
        grammar instead.

        **Key Arguments:**

        - ``recipe`` -- the recipe name to filter product SOFs by
        - ``arm`` -- the instrument arm to filter product SOFs by
        - ``ttype`` -- the raw `eso dpr type` this group was built from (adds an extra filter for `STD` types)
        - ``calType`` -- the list of calibration type names this recipe depends on. This is only ever called for a recipe that has at least one, so an empty list is not a supported input and produces invalid SQL (`AND ;`), matching the pre-refactor behaviour it replaces.

        **Return:**

        - ``sqlQuery`` -- the parameterized select
        - ``sqlParams`` -- the parameters to bind against ``sqlQuery``

        **Raises:**

        - ``UnsafeSqlIdentifierError`` -- if any entry in ``calType`` fails the safe-identifier grammar
        """
        params = [recipe, arm]

        extraType = ""
        if "STD" in ttype:
            extraType = 'AND "eso dpr type" = ?'
            params.append(ttype)

        calTables = [validate_sql_identifier(f"cal_{ct}", "calibration table name") for ct in calType]
        exists = " AND ".join(
            f"EXISTS (SELECT 1 FROM {calTable} WHERE {calTable}.sof = p.sof "  # noqa: S608
            f"AND ({calTable}.upstream_status = 'pass' OR {calTable}.upstream_status IS NULL))"
            for calTable in calTables
        )

        sqlQuery = f"""select sof from product_frames_plus as p
            WHERE complete < 1
            AND recipe = ?
            AND "eso seq arm" = ?
            {extraType}
            AND {exists};"""  # noqa: S608

        return sqlQuery, tuple(params)

    def _calibration_completeness_update_query(self, containerSofs):
        """*build the `(sql, params)` pair that marks product SOFs calibration-complete*

        `containerSofs` are sof filenames read back from a prior query, so
        they are bound as `?` parameters rather than joined into the SQL
        text. An empty list matches no row, since `sof in ()` is invalid
        SQLite syntax.

        **Key Arguments:**

        - ``containerSofs`` -- the list of sof filenames to mark complete

        **Return:**

        - ``sqlQuery`` -- the parameterized update
        - ``sqlParams`` -- the parameters to bind against ``sqlQuery``
        """
        if not containerSofs:
            return (
                "UPDATE product_frames as p SET complete = -1 WHERE complete < 1 And sof in (NULL)",
                (),
            )

        placeholders = ", ".join("?" for _ in containerSofs)
        sqlQuery = (
            f"UPDATE product_frames as p SET complete = -1 WHERE complete < 1 "  # noqa: S608
            f"And sof in ({placeholders})"
        )
        return sqlQuery, tuple(containerSofs)

    def _calibration_raw_frames_query(self, calibrationType):
        """*build the select that pulls raw calibration frames newly marked complete for one calibration type*

        `calibrationType` is used as a *table-name* identifier (`cal_<type>`),
        which SQLite cannot bind, so it is checked against the safe-identifier
        grammar instead.

        **Key Arguments:**

        - ``calibrationType`` -- the calibration type name (e.g. `bias`, `dark`)

        **Return:**

        - ``sqlQuery`` -- the select, safe to run with no bound parameters

        **Raises:**

        - ``UnsafeSqlIdentifierError`` -- if ``calibrationType`` fails the safe-identifier grammar
        """
        calTable = validate_sql_identifier(f"cal_{calibrationType}", "calibration table name")
        return (
            f"select {calTable}.file, {calTable}.upstream_tag as tag, product_frames.sof, "  # noqa: S608
            f"{calTable}.filepath, product_frames.complete from product_frames, {calTable} "
            f"where product_frames.complete = -1 and product_frames.sof={calTable}.sof;"
        )

    def build_sof_files(self, rejoinLateFrames=False):  # noqa: PLR0915
        """*scan the raw frame table to generate the listing of products that are expected to be created and then write out all of the needed SOF files*

        **Key Arguments:**

        - ``rejoinLateFrames`` -- when True, a raw frame added after its set was grouped rejoins that set, which is
          reset to be grouped and reduced again. Only `prepare` asks for this; `session_refresh` does not

        **Usage:**

        ```python
        self.build_sof_files()
        ```
        """
        self.log.debug("starting the ``_populate_product_frames_db_table`` method")

        import pandas as pd

        statusColumn = validate_sql_identifier(f"status_{self.sessionId}", "status column")

        c = self.conn.cursor()
        sqlQuery = f"update product_frames set status = {statusColumn};"  # noqa: S608
        c.execute(sqlQuery)
        sqlQuery = "update raw_frames set processed = 0 where processed < 0;"
        c.execute(sqlQuery)

        # CLEAN UP FAILED FILES
        # DELETE FROM
        count = 0
        oldCount = -1
        while count != oldCount:
            oldCount = count
            c = self.conn.cursor()
            sqlQuery = "select distinct sof from sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1));"
            compromisedSofs = pd.read_sql(sqlQuery, con=self.conn)["sof"].tolist()
            count = len(compromisedSofs)

            sqlQuery = "update product_frames set complete = 0 where (status != 'fail' or status is null) and sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1)));"
            c.execute(sqlQuery)

        sofMapTableName = validate_sql_identifier(f"sof_map_{self.sessionId}", "sof map table name")
        sqlQueries = [
            "update raw_frames set processed = 0 where file in (select file from sof_map where sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1))));",
            "update raw_frame_sets set complete = 0 where sof in (select distinct sof from  sof_map where filepath in (  select p.filepath from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1)));",
            f"delete from {sofMapTableName} where sof in (  select s.sof from sof_map s, product_frames p where p.filepath=s.filepath and (p.status = 'fail' or p.complete < 1));",  # noqa: S608
            "update raw_frames set processed = -1 where file in (select distinct s.file from sof_map s, product_frames p where p.sof=s.sof and p.status = 'fail');",
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

        # A FRAME ADDED AFTER ITS SET WAS GROUPED REJOINS THAT SET, WHICH IS THEN GROUPED AND REDUCED AGAIN
        rejoined = self._rejoin_late_frames(sofMapTableName) if rejoinLateFrames else {}
        compromisedSofs = compromisedSofs + sorted(rejoined)

        # DELETE COMPROMISED SOF FILES
        sofDir = Path(self.sessionPath) / "sof"
        for sof in compromisedSofs:
            try:
                os.remove(_validate_owned_path(sofDir / sof, sofDir, "SOF path"))
            except _UnsafePathError as e:
                self.log.warning(f"build_sof_files: the SOF name `{sof}` is not a safe path, not deleted: {e}")
            except OSError as e:
                self.log.debug(f"build_sof_files: the SOF file `{sof}` could not be deleted, continuing: {e}")

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
                rawGroups = pin_rejoined_sof_names(rawGroups, self.pinnedSofNames)

                # ADD PREDICTED PRODUCT TO PRODUCT TABLE - DETERMINE IF COMPLETE LATER
                incompleteProducts = self.predict_product_frames(productTypes, rawGroups, recipe)

                if not incompleteProducts:
                    continue

                allRawGroups.append(rawGroups)

                if not len(calibrationTypes):
                    # MBIAS AND MDARK -- ALWAYS COMPLETE (NO PRIOR CALIBRATION REQUIRED)
                    sqlQuery = "select sof from product_frames where recipe = ? and complete = 0;"
                    containerSofs = pd.read_sql(sqlQuery, con=self.conn, params=(recipe,))["sof"].tolist()
                    self.raw_frames_to_sof_map(rawGroups=rawGroups, containerSofs=containerSofs)
                    sqlQuery = "update product_frames set complete = 1 where recipe = ? and complete = 0;"
                    c.execute(sqlQuery, (recipe,))

                else:

                    if isinstance(calibrationTypes, dict):
                        for arm, calType in calibrationTypes.items():

                            sqlQuery, sqlParams = self._calibration_completeness_select_query(
                                recipe, arm, ttype, calType
                            )

                            containerSofs = pd.read_sql(sqlQuery, con=self.conn, params=sqlParams)["sof"].tolist()

                            self.raw_frames_to_sof_map(rawGroups=rawGroups, containerSofs=containerSofs)

                            sqlQuery, sqlParams = self._calibration_completeness_update_query(containerSofs)
                            c.execute(sqlQuery, sqlParams)

                            # FOR COMPLETE PRODUCTS, ADD CALIBRATION FILES TO SOF MAP
                            # NEED TO ALSO ADD THE RAW FILES TOO ... ADD RAW FRAMES, SET COMPLETE = 1 WHERE PRODUCT FRAMES COMPLETE = 1
                            for ct in calType:

                                sqlQuery = self._calibration_raw_frames_query(ct)
                                newSof = pd.read_sql(sqlQuery, con=self.conn)

                                if len(newSof):
                                    self._dataframe_to_sqlite(
                                        newSof,
                                        sofMapTableName,
                                        replace=False,
                                    )

            sqlQuery = """update product_frames set complete = 1 where complete = -1;"""
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
            f"UPDATE {sofMapTableName} SET complete = 1 WHERE complete = -1;",  # noqa: S608
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

    def _rejoin_late_frames(self, sofMapTableName):
        """*reset every set that a late-arriving raw frame belongs to, so it is grouped and reduced again*

        **Key Arguments:**

        - ``sofMapTableName`` -- the current session's validated `sof_map_<id>` table name

        **Return:**

        - ``rejoined`` -- SOF name to (recipe, frozenset of raw frame filepaths) for the sets that gained a frame in
          this call, empty when none did. Their names stay pinned in `self.pinnedSofNames` for the rest of this
          instance's life

        The pin lives only on this instance. A rejoined set that is still unprocessed after both `prepare()` passes
        (for example, while a calibration is missing) loses its pin when the process exits, and a later `prepare()`
        then names the set after its earliest frame.
        """
        groupingKeys = [key for key in self.filterKeywords if key not in self.proKeywords]
        rejoined = find_rejoined_sofs(self.conn, sofMapTableName, groupingKeys)
        if not rejoined:
            return rejoined
        self.pinnedSofNames = {**self.pinnedSofNames, **rejoined}

        stalePaths = reset_rejoined_sofs(
            self.conn, sofMapTableName, rejoined, self.rootDir, _validate_owned_path, self.log
        )
        self.conn.commit()
        # DELETE ONLY AFTER THE COMMIT, SO A FAILED COMMIT NEVER LEAVES A PRODUCT THE DATABASE STILL EXPECTS DELETED
        delete_stale_files(stalePaths, self.log)
        self.log.debug(f"_rejoin_late_frames: reset {len(rejoined)} SOF(s) that gained late frames")
        print(format_rejoin_summary(rejoined))
        return rejoined

    def get_raw_frames_and_groups(
        self,
        ttype=None,
        arm=None,
        tech=None,
        recipe=None,
        recipeOrder=None,
        filterName=None,
        unprocessedOnly=False,
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

        # IF NONE, SET TO EMPTY STRING
        ttype, arm, tech = ttype or "", arm or "", tech or ""

        whereClauses = []
        params: list = []
        if ttype:
            whereClauses.append("`eso dpr type` = ?")
            params.append(ttype)
        if arm:
            whereClauses.append("`eso seq arm` = ?")
            params.append(arm)
        if tech:
            placeholders = ", ".join("?" for _ in tech)
            whereClauses.append(f"`eso dpr tech` in ({placeholders})")
            params.extend(tech)
        if unprocessedOnly:
            whereClauses.append("processed = 0")

        where = ""
        if whereClauses:
            where = "where " + " and ".join(whereClauses)

        # READ IN RAW FRAMES TABLE
        conn = self.conn
        rawFrames = pd.read_sql(
            f"SELECT * FROM raw_frames_valid {where} order by `mjd-obs` asc",  # noqa: S608
            con=conn,
            params=params,
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

            recipe = recipe.replace("_std", "") + "_std_flux"

        if "STD,TELLURIC" in ttype:

            recipe = recipe.replace("_std", "") + "_std_tell"

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

        if "set_first_file" in filterKeywordsRaw:
            filterKeywordsRaw.remove("set_first_file")

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

        if not len(rawGroups.index):
            sqlQuery = "select count(*) from product_frames where recipe = ? and complete< 1;"
            c = self.conn.cursor()
            c.execute(sqlQuery, (recipe,))
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
                productFrames["file"] = productFrames["sof"].str.replace(".sof", ".fits")
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
            # `placeholders` IS ALWAYS LITERAL `?` MARKS -- THE ACTUAL VALUES
            # ARE BOUND BELOW VIA `processedRawFiles`, NEVER INTERPOLATED.
            sqlQuery = f"update raw_frames set processed=1 where file in ({placeholders});"  # noqa: S608
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
        - `UnsafeSqlIdentifierError` if `table_name` fails the safe-identifier grammar.
        """
        import time

        # A TABLE NAME CANNOT BE A BOUND PARAMETER, SO IT IS VALIDATED AGAINST
        # THE SAFE-IDENTIFIER GRAMMAR BEFORE IT IS INTERPOLATED BELOW.
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
        while keepTrying < 7:
            try:
                dataframe.replace(["--", -99.99], None).to_sql(
                    table_name, con=self.conn, index=False, if_exists="append"
                )
                keepTrying = 10
            except Exception:
                if keepTrying > 5:
                    raise
                time.sleep(1)
                keepTrying += 1


def _harvest_fits_headers(batch, log, pathToDirectory, keywords, filterKeys, instrument, kw):
    import numpy as np
    from astropy.time import Time, TimeDelta
    from ccdproc import ImageFileCollection

    masterTable = ImageFileCollection(filenames=batch, keywords=keywords)
    masterTable = masterTable.summary
    # ADD FILLED VALUES FOR MISSING CELLS
    for fil in keywords:
        if fil in filterKeys and fil not in ["exptime"]:

            try:
                masterTable[fil].fill_value = "--"
            except (TypeError, ValueError) as e:
                log.debug(f"_harvest_fits_headers: `masterTable[fil].fill_value = '--'` failed, continuing: {e}")
                masterTable.replace_column(fil, masterTable[fil].astype(str))
                masterTable[fil].fill_value = "--"
        # elif fil in ["exptime"]:
        #     masterTable[fil].fill_value = "--"
        else:
            try:
                masterTable[fil].fill_value = -99.99
            except (TypeError, ValueError) as e:
                log.debug(f"_harvest_fits_headers: `masterTable[fil].fill_value = -99.99` failed, continuing: {e}")
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
            except (TypeError, ValueError) as e:
                log.debug(f"_harvest_fits_headers: `masterTable['rospeed'][masterTable['ro...` failed, continuing: {e}")
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
            except (TypeError, ValueError) as e:
                log.debug(f"_harvest_fits_headers: `masterTable['rospeed'][masterTable['ro...` failed, continuing: {e}")
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
            masterTable["binning"],
            masterTable[kw("WIN_BINY").lower()].astype("int").astype("str"),
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
        except (TypeError, ValueError, KeyError) as e:
            log.debug(f"_harvest_fits_headers: `masterTable.add_index(k)` failed, continuing: {e}")

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
