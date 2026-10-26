#!/usr/bin/env python
"""
*reduce all the data in a workspace, or target specific obs and files for reduction*

Author
: David Young

Date Created
: January 17, 2024
"""

# from memory_profiler import profile
import contextlib
import logging
import os
import sys
import time
from collections.abc import Callable
from dataclasses import dataclass

os.environ["TERM"] = "vt100"


def _sql_in_placeholders(values):
    """*build a bound `?` placeholder clause and its parameter list for a SQL `in (...)` test*

    An empty `values` produces a clause that matches no row, since `sof in
    ()` is invalid SQLite syntax.

    **Key Arguments:**

    - ``values`` -- the list of values the `in (...)` clause should match

    **Return:**

    - ``clause`` -- the `in (...)` clause body, either `?, ?, ...` or `NULL`
    - ``params`` -- the parameters to bind against `clause`, in order
    """
    if not values:
        return "NULL", []
    return ", ".join("?" for _ in values), list(values)


class reducer:
    """
        *reduce all the data in a workspace, or target specific obs and files for reduction*

        **Key Arguments:**

        - ``log`` -- logger
        - ``workspaceDirectory`` -- path to the root of the workspace
        - ``reductionTarget`` -- target for reduction: "all", "sof", "ob" (default: "all")
        - ``settings`` -- the settings dictionary
        - ``pathToSettings`` -- path to the settings file.
        - ``quitOnFail`` -- quit the pipeline on any recipe failure
        - ``overwrite`` -- overwrite existing reductions. Default *False*.
        - ``daemon`` -- run in daemon mode (no terminal output). Default *False*.
        - ``verbose`` -- print verbose output to terminal. Default *False*.
        - ``refreshWorkspace`` -- refresh the workspace before reducing to collect new files. Default *False*.
    `

        **Usage:**

        ```python
        from soxspipe.commonutils import reducer
        collection = reducer(
            log=log,
            workspaceDirectory="/path/to/workspace/root/",
            reductionTarget="all",
            settings=settings,
            pathToSettings="/path/to/settings.yaml"
        )
        collection.reduce()
        ```

    """

    def __init__(
        self,
        log,
        workspaceDirectory,
        reductionTarget="all",
        settings=False,
        pathToSettings=False,
        quitOnFail=False,
        overwrite=False,
        daemon=False,
        verbose=False,
        refreshWorkspace=False,
    ):
        self.log = log
        log.debug("instantiating a new 'reducer' object")
        self.settings = settings
        self.workspaceDirectory = workspaceDirectory
        self.reductionTarget = reductionTarget
        self.overwrite = overwrite
        self.pathToSettings = pathToSettings
        self.quitOnFail = quitOnFail
        self.daemon = daemon
        self.verbose = verbose

        # REQUEST THE WORKSPACE PARAMETERS FROM THE DATA-ORGANISER
        from soxspipe.commonutils import data_organiser

        do = data_organiser(
            log=log,
            rootDir=workspaceDirectory,
            dbConnect=False,
        )
        self.sessionId, allSessions = do.session_list(silent=True)
        do.close()

        if self.sessionId is None:
            return

        self.recipeList = [
            "mbias",
            "mdark",
            "disp_solution",
            "order_centres",
            "mflat",
            "spat_solution",
            "nod_std",
            "stare_std",
            "offset_std",
            "nod_obj",
            "stare_obj",
            "offset_obj",
        ]

        self.sessionPath = workspaceDirectory + "/sessions/" + self.sessionId
        self.sessionDB = workspaceDirectory + "/soxspipe.db"

        if refreshWorkspace:
            self.log.print("Refreshing the workspace before reduction...")
            do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
            do.prepare(refresh=False, report=False)
            do.close()

        return

    def reduce(self, batch=False, multiprocess=False, poolSize=None):  # noqa: PLR0915
        """
        *reduce the selected data*

        **Key Arguments:**

        - ``batch`` -- the maximum number of recipe batches to run. Default *False* (no limit).
        - ``multiprocess`` -- run each recipe's reductions in a worker pool. Default *False*.
        - ``poolSize`` -- the number of worker processes in multiprocessing mode. Default *None*, which keeps
          the built-in choice of ``run_recipe_bulk``.
        """
        self.log.debug("starting the ``reduce`` method")

        if self.sessionId is None:
            print("Please prepare this workspace using `soxspipe prep` before attempting to reduce the data.")
            return

        import traceback

        from fundamentals import times

        from soxspipe.commonutils import data_organiser

        do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
        do.session_refresh(failure=None)
        do.close()

        if self.reductionTarget != "all":
            self.recipeList = [False]

        if not batch:
            batch = 100000000

        batchCount = 0
        for rootRecipe in self.recipeList:
            while True and batchCount < batch:
                # rawGroups WILL CONTAIN ONE RECIPE COMMAND PER ENTRY
                rawGroups = self.select_sof_files_to_process(
                    recipe=rootRecipe, reductionTarget=self.reductionTarget, arm=False
                )

                if rawGroups.empty:
                    break

                if multiprocess:
                    import sqlite3 as sql

                    from soxspipe.commonutils.data_organiser import (
                        _validate_owned_path,
                    )

                    databasePath = _validate_owned_path(
                        self.sessionDB,
                        self.workspaceDirectory,
                        "database path",
                    )
                    conn = sql.connect(
                        databasePath,
                        timeout=300,
                        autocommit=True,
                        check_same_thread=False,
                    )
                    c = conn.cursor()
                    c.execute("PRAGMA busy_timeout = 100000")
                    c.execute("PRAGMA synchronous = OFF")

                    sofList = rawGroups["sof"].tolist()
                    sofList = [self.sessionPath + "/sof/" + sof for sof in sofList]
                    # time the bulk run and report when finished
                    startTime = times.get_now_sql_datetime()
                    run_recipe_bulk(
                        log=self.log,
                        recipe=rootRecipe,
                        sofList=sofList,
                        commandList=rawGroups["command"].tolist(),
                        settings=self.settings,
                        overwrite=self.overwrite,
                        workspaceDirectory=self.workspaceDirectory,
                        conn=conn,
                        sessionId=self.sessionId,
                        poolSize=poolSize,
                    )
                    endTime = times.get_now_sql_datetime()
                    runningTime = times.calculate_time_difference(startTime, endTime)
                    self.log.print(f"Multiprocess for {rootRecipe} recipe completed for {len(sofList)} files.")
                    self.log.print(f"Multiprocess Recipe Run Time: {runningTime}\n\n")
                    break

                fail = False
                for index, row in rawGroups.iterrows():
                    if batchCount >= batch:
                        self.log.print(f"Batch limit of {batch} reached, pausing reductions.")
                        break

                    recipe = row["recipe"].replace("_obj", "")
                    sof = row["sof"]
                    startTime = times.get_now_sql_datetime()
                    sof = self.sessionPath + "/sof/" + sof

                    try:
                        run_recipe(
                            self.log,
                            recipe,
                            sof,
                            settings=self.settings,
                            overwrite=self.overwrite,
                            command=row["command"],
                            verbose=self.verbose,
                        )
                        batchCount += 1
                    except FileExistsError:
                        continue
                    except Exception:
                        # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
                        self.log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
                        self.log.error(
                            f"\nRecipe Command: {row['command'].replace('-obj ', ' ').replace('-std ', ' ')}\n\n"
                        )
                        fail = True

                        if self.quitOnFail:
                            sys.exit(1)

                        if self.reductionTarget != "all":
                            self.overwrite = False

                        if not self.daemon:
                            print(f"{'=' * 70}\n")

                    ## FINISH LOGGING ##
                    endTime = times.get_now_sql_datetime()
                    runningTime = times.calculate_time_difference(startTime, endTime)
                    sys.argv[0] = os.path.basename(sys.argv[0])

                    self.log.print(f"\nRecipe Command: {row['command'].replace('_obj ', ' ')} ")
                    self.log.print(f"Recipe Run Time: {runningTime}\n\n")
                    if not self.daemon:
                        print(f"{'=' * 70}\n")

                if fail:
                    do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
                    reset = do.session_refresh()
                    do.close()
                    if reset:
                        print(f"BACK TO THE START! {rootRecipe}\n\n")
                        break
                break

        if self.reductionTarget == "all":
            from soxspipe.commonutils.data_organiser import print_incomplete_sets_report

            do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
            incompleteSets, blockingSets = do.get_incomplete_sets_report()
            do.close()
            print_incomplete_sets_report(incompleteSets, blockingSets)

        do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
        do.session_refresh(failure=None)
        do.close()

        self.log.debug("completed the ``reduce`` method")
        return

    def _passed_this_session_clause(self):
        """*return the SQL condition that leaves out the SOFs that already passed this session*

        The product file cannot mark a SOF as done, as a stare SOF reduced frame by frame never writes the predicted
        product. The session status can. ``product_frames`` holds one row per product file, so a SOF is left out only
        when it has a passing row and no row that is not passing. The match is on ``sof`` with ``=``, so a row with
        no ``sof`` matches nothing.

        **Return:**

        - ``clause`` -- text starting with ``and``, to append to a query on ``raw_frame_sets``
        """
        from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

        statusColumn = validate_sql_identifier(f"status_{self.sessionId}", "status column")
        rowsOfThisSof = "select 1 from product_frames p where p.sof = raw_frame_sets.sof"
        return (
            f" and not (exists ({rowsOfThisSof} and p.{statusColumn} = 'pass') "  # noqa: S608
            f"and not exists ({rowsOfThisSof} and p.{statusColumn} is not 'pass')) "  # noqa: S608
        )

    def select_sof_files_to_process(self, recipe=False, reductionTarget=False, batch=False, arm=False):
        """*select all of the SOF files still requiring processing*

        With the target ``all``, a SOF whose ``status_<session>`` is ``pass`` is left out unless the reducer was
        given ``overwrite``. A SOF that failed or has not run is selected again. A named SOF is always selected.

        **Key Arguments:**
            - ``recipe`` -- the name of the recipe to filter by (optional)
            - ``reductionTarget`` -- target for reduction: "all", "sof", "ob" (default: False)
            - ``batch`` -- number of SOF files to return (default: False, all)
            - ``arm`` -- filter by arm (default: False, all)

        **Return:**

        - `rawGroups` -- a dataframe of the containing a list of recipes and sof file paths

        **Usage:**

        ```python
        rawGroups = reducer.select_sof_files_to_process()
        ```
        """
        self.log.debug("starting the ``select_sof_files_to_process`` method")

        import sqlite3 as sql

        import pandas as pd

        from soxspipe.commonutils.data_organiser import _validate_owned_path

        databasePath = _validate_owned_path(self.sessionDB, self.workspaceDirectory, "database path")
        conn = sql.connect(databasePath, timeout=300, autocommit=True, check_same_thread=False)
        c = conn.cursor()
        c.execute("PRAGMA busy_timeout = 100000")
        c.execute("PRAGMA synchronous = OFF")

        allBranchParams = {}

        if batch:
            limitText = " LIMIT :batch "
            allBranchParams["batch"] = batch
        else:
            limitText = ""

        if arm:
            armText = " and `eso seq arm` = :arm "
            allBranchParams["arm"] = arm
        else:
            armText = ""

        if reductionTarget == "all":
            # GET THE GROUPS OF FILES NEEDING REDUCED, ASSIGN THE CORRECT COMMAND TO EXECUTE THE RECIPE
            if not recipe:
                recipeText = "is not null"
            else:
                recipeText = "= :recipe"
                allBranchParams["recipe"] = recipe

            # `recipeText`, `armText` AND `limitText` ARE ALWAYS EITHER STATIC
            # TEXT OR A `:name` PLACEHOLDER -- NEVER A RAW VALUE -- SO THE
            # ACTUAL `recipe`/`arm`/`batch` VALUES REACH SQLITE ONLY THROUGH
            # `allBranchParams` BELOW.
            passedText = "" if self.overwrite else self._passed_this_session_clause()

            rawGroups = pd.read_sql(
                f"SELECT * FROM raw_frame_sets where recipe_order is not null "  # noqa: S608
                f"and complete = 1 and recipe {recipeText} {armText} {passedText} "  # noqa: S608
                f"order by recipe_order, sof {limitText}",
                con=conn,
                params=allBranchParams,
            )

        elif reductionTarget.split(".")[-1].lower() == "sof":
            # `reductionTarget` IS THE ONLY EXTERNAL VALUE HERE. IT IS BOUND
            # ONCE, AS A NAMED PARAMETER, AND REUSED ACROSS EVERY NESTED
            # SUBQUERY BELOW -- SQLITE ALLOWS THE SAME NAMED PLACEHOLDER TO
            # REPEAT WITHIN ONE STATEMENT.
            sqlQuery = "select sof from product_frames where sof = :reductionTarget and complete = 1"

            for _ in range(4):  # Recursively query up to 5 times
                sqlQuery = f"SELECT distinct sof FROM product_frames WHERE file IN (SELECT file FROM sof_map_base WHERE sof in ({sqlQuery})) or sof in ({sqlQuery})"  # noqa: S608

            sqlQuery = (
                f"SELECT distinct sof, recipe from raw_frame_sets WHERE sof in ({sqlQuery}) "  # noqa: S608
                "order by recipe_order, sof"
            )

            rawGroups = pd.read_sql(sqlQuery, con=conn, params={"reductionTarget": reductionTarget})

        if not len(rawGroups.index):
            if reductionTarget != "all":
                self.log.warning("The SOF file selected for processing is either missing or incomplete.")
            else:
                self.log.info("No SOF files require processing.")
            conn.close()
            return pd.DataFrame(columns=["recipe", "sof", "command"])

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        if recipe:
            mask = rawGroups["recipe"] == recipe
            rawGroups = rawGroups.loc[mask]

        rawGroups["command"] = (
            "soxspipe "
            + rawGroups["recipe"].str.replace("-obj", "").str.replace("_obj", "")
            + " sof/"
            + rawGroups["sof"]
        )
        if self.pathToSettings:
            rawGroups["command"] += f" -s {self.pathToSettings}"
        conn.close()

        self.log.debug("completed the ``select_sof_files_to_process`` method")
        return rawGroups[["recipe", "sof", "command"]].drop_duplicates()


def run_recipe(log, recipe, sof, settings, overwrite, command=False, verbose=False, turnOffMP=False):
    """*execute a pipeline recipe*

    **Key Arguments:**

    - ``recipe`` -- the name of the recipe tp execute
    - ``sof`` -- path to the sof file containing the files the recipe requires
    - ``command`` -- the command used to run the recipe
    - ``settings`` -- the settings dictionary
    - ``overwrite`` -- overwrite existing reductions. Default *False*.
    - ``verbose`` -- print verbose output to terminal. Default *False*.
    - ``turnOffMP`` -- turn off multiprocessing mode. Default *False*.

    **Usage:**

    ```python
    reducer.run_recipe("mbias", "/path/to/sofs/my_bias_files.sof")
    ```
    """
    log.debug("starting the ``run_recipe`` method")

    if recipe == "mbias":
        from soxspipe.recipes import soxs_mbias

        soxs_recipe = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if recipe == "mdark":
        from soxspipe.recipes import soxs_mdark

        soxs_recipe = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if recipe == "disp_solution":
        from soxspipe.recipes import soxs_disp_solution

        soxs_recipe = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if recipe == "order_centres":
        from soxspipe.recipes import soxs_order_centres

        soxs_recipe = soxs_order_centres(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if recipe == "mflat":
        from soxspipe.recipes import soxs_mflat

        soxs_recipe = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if recipe == "spat_solution":
        from soxspipe.recipes import soxs_spatial_solution

        soxs_recipe = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if "stare" in recipe:
        from soxspipe.recipes import soxs_stare

        soxs_recipe = soxs_stare(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if "nod" in recipe:
        from soxspipe.recipes import soxs_nod

        soxs_recipe = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    if "offset" in recipe:
        from soxspipe.recipes import soxs_offset

        soxs_recipe = soxs_offset(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )

    productPath, qcTable = soxs_recipe.produce_product()
    del soxs_recipe

    log.debug("completed the ``run_recipe`` method")
    return productPath, qcTable


def _worker_log(log, settings, logQueue):
    """*return a logger that can write from inside a spawned worker process*

    A `logging.Logger` pickles by name, so a spawned worker receives a bare logger with no
    handlers and no ``print`` level. This adds the ``print`` level, sets the root level from the
    ``logging settings`` block and sends every record to the parent process through ``logQueue``,
    where `_worker_log_queue` writes it. A worker owns no file handler, so only the parent writes
    and rotates the shared log. A logger that already works (it has handlers and a ``print``
    method, as in a forked worker, a reused spawned worker, or a test double) is returned
    unchanged.

    The changes are process-wide: the ``print`` level, the root logger's level and queue handler, and
    ``logging.captureWarnings(True)``. They are safe in a spawned worker, which exists only to run reductions,
    but do not call this in a process that has logging of its own.

    **Key Arguments:**

    - ``log`` -- the logger the worker was given
    - ``settings`` -- the settings dictionary, holding the ``logging settings`` block
    - ``logQueue`` -- the queue that carries the worker's log records to the parent process
    """
    from logging.handlers import QueueHandler

    from fundamentals import logs

    isWorkingLogger = not isinstance(log, logging.Logger) or (hasattr(log, "print") and log.hasHandlers())
    if isWorkingLogger or "logging settings" not in settings:
        return log

    logs.addLoggingLevel("PRINT", logging.INFO + 5)
    rootLog = logging.getLogger()
    rootLog.setLevel(settings["logging settings"].get("root", {}).get("level", logging.INFO))
    # NO LEVEL IS SET ON THE QUEUE HANDLER: THE PARENT'S FILE HANDLERS NEED THE INFO AND PRINT RECORDS
    rootLog.addHandler(QueueHandler(logQueue))
    logging.captureWarnings(True)
    return logging.getLogger(logs.__name__)


class _parent_log_forwarder(logging.Handler):
    """*write the records that workers queued to the handlers the parent process already has*

    Each record goes to the handlers of the logger it was logged on and of that logger's parents, as the parent's
    own logging would send it. The walk is needed because the shipped setup puts the file handler on the root
    logger and the terminal handler on the ``fundamentals.logs`` logger, so neither logger alone reaches both.
    Each handler keeps its own level.

    Bulk workers are muted: a handler that is not a ``logging.FileHandler`` is treated as a terminal handler,
    and gets only warnings and above. This is the rule the workers' own handlers used before the queue.

    A `QueueHandler` formats a record before it queues it, so the traceback of an ``exc_info`` is folded into the
    message text and ``exc_info`` and ``stack_info`` are dropped. The line the parent's file handler writes is
    the same as before.
    """

    def emit(self, record):
        logger = logging.getLogger(record.name)
        while logger is not None:
            for handler in logger.handlers:
                isShown = isinstance(handler, logging.FileHandler) or record.levelno >= logging.WARNING
                if isShown and record.levelno >= handler.level:
                    handler.handle(record)
            logger = logger.parent if logger.propagate else None


def _report_log_failure(log, action, function):
    """*run ``function``, and log an error on ``log`` instead of raising if it fails*

    The log plumbing runs after the reductions finish, and a failure in it must not hide the results or
    replace an exception that is already propagating.
    """
    try:
        function()
    except Exception as error:
        log.error(f"Could not {action}, so some worker log records may be missing: {error}")


@contextlib.contextmanager
def _worker_log_queue(log):
    """*run a listener in this process that writes the log records spawned workers queue, and yield the queue*

    The queue comes from a manager process, so it pickles into any spawned worker and a worker killed
    mid-record cannot leave a lock held that blocks the others. The listener stops, after it has written
    every queued record, however the block ends. If the listener or the manager fails to stop (for example
    because the manager process was killed), the failure is logged on ``log`` and not raised, so it can never
    replace an exception from the block.

    Records cross the queue as formatted text: see `_parent_log_forwarder`.

    **Key Arguments:**

    - ``log`` -- the parent logger, used to report a failure to stop the listener or the manager
    """
    import multiprocessing
    from logging.handlers import QueueListener

    manager = multiprocessing.get_context("spawn").Manager()
    try:
        logQueue = manager.Queue()
        listener = QueueListener(logQueue, _parent_log_forwarder())
        listener.start()
        try:
            yield logQueue
        finally:
            _report_log_failure(log, "stop the worker log listener", listener.stop)
    finally:
        _report_log_failure(log, "shut down the worker log manager", manager.shutdown)


POOL_TIMEOUT_SECONDS = 36000
WORKER_CRASH_MESSAGE = (
    "The reduction worker process crashed (killed by a signal) while this SOF was in progress. "
    "Check the system crash reports, then run this SOF again on its own."
)
WORKER_START_MESSAGE = "The reduction worker processes could not start, so this SOF was not run."
WORKER_TIMEOUT_MESSAGE = (
    "The batch deadline passed before this reduction finished. Its worker process was stopped, or it never started."
)
# THE TIME A WORKER GETS TO EXIT AFTER SIGTERM BEFORE IT IS KILLED
WORKER_STOP_SECONDS = 10


def _bulk_worker(
    inputDict,
    log,
    recipe,
    settings,
    overwrite,
    workspaceDirectory,
    logQueue,
    wrapperTurnOffMP=True,
):
    """*run one recipe for one SOF inside a bulk worker process and describe the outcome*

    **Key Arguments:**

    - ``inputDict`` -- a dictionary holding the ``sof`` path and the ``command`` that runs it
    - ``log`` -- logger
    - ``recipe`` -- the name of the recipe to execute
    - ``settings`` -- the settings dictionary
    - ``overwrite`` -- overwrite existing reductions
    - ``workspaceDirectory`` -- path to the root of the workspace
    - ``logQueue`` -- the queue that carries this worker's log records to the parent process
    - ``wrapperTurnOffMP`` -- turn off the multiprocessing inside the recipe. Default *True*.

    **Return:**

    - ``returnDict`` -- the ``status``, ``sof``, ``productPath``, ``qcTable`` and ``error_message`` of the run.
      The ``error_message`` is a string, so the result always pickles back to the parent process
    """
    import traceback

    log = _worker_log(log, settings, logQueue)

    returnDict = {
        "status": None,
        "sof": inputDict["sof"],
        "productPath": None,
        "qcTable": None,
        "error_message": None,
    }

    try:
        productPath, qcTable = run_recipe(
            log=log,
            recipe=recipe,
            sof=inputDict["sof"],
            settings=settings,
            overwrite=overwrite,
            command=inputDict["command"],
            turnOffMP=wrapperTurnOffMP,
        )
        returnDict["status"] = "pass"
        returnDict["productPath"] = productPath
        returnDict["qcTable"] = qcTable
        mask = qcTable["qc_flag"] == "fail"
        failedQcs = qcTable[mask]
        if len(failedQcs):
            returnDict["status"] = "fail"
            failedQcs = failedQcs["qc_name"].tolist()
            returnDict["error_message"] = (
                f"The following QC values are outside of acceptable limits: {', '.join(failedQcs)}."
            )
    except FileExistsError as e:
        if "previously failed" in str(e):
            returnDict["status"] = "previous-fail"
        elif "product of this recipe" in str(e):
            returnDict["status"] = "previous-pass"
        else:
            returnDict["status"] = "fail"
            returnDict["error_message"] = str(e)
    except Exception as e:
        # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
        log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
        log.error(f"\nRecipe Command: {inputDict['command'].replace('-obj ', ' ').replace('_obj ', ' ')}\n\n")
        returnDict["status"] = "fail"
        returnDict["error_message"] = str(e)

    return returnDict


def _failed_result(inputDict, message):
    """*the result of a SOF that failed outside the recipe, in the shape ``_bulk_worker`` returns*"""
    return {
        "status": "fail",
        "sof": inputDict["sof"],
        "productPath": None,
        "qcTable": None,
        "error_message": message,
    }


def _mute_worker_stdout():
    """*silence the printed output of a bulk worker process, as ``fmultiprocess(mute=True)`` did*"""
    sys.stdout = open(os.devnull, "w")  # noqa: SIM115


def _run_marked(markerPath, function, inputDict, log, kwargs):
    """*touch a marker file to show the task started, then run it*

    A pool that breaks cannot say which tasks were running, because every unfinished future fails the same way.
    The marker file tells the parent which SOFs had started.
    """
    open(markerPath, "w").close()
    return function(inputDict, log=log, **kwargs)


@dataclass(frozen=True)
class _pool_task:
    """The fixed parts of a pool run: what to call, on what, and where the start markers go."""

    function: Callable[..., dict]
    inputDicts: list[dict]
    log: logging.Logger
    kwargs: dict
    markerDirectory: str

    def marker_path(self, index):
        """Return the path of the file that shows the task for ``inputDicts[index]`` has started."""
        return os.path.join(self.markerDirectory, f"{index}.started")


@dataclass(frozen=True)
class _round_outcome:
    """What one pool round settled.

    ``resolved`` maps an input index to its result. ``unstarted`` lists the inputs whose function body was never
    entered, so they are safe to run again. ``inFlight`` lists the inputs that had started when the pool broke.
    """

    resolved: dict[int, dict]
    unstarted: list[int]
    inFlight: list[int]


def _run_in_spawn_pool(function, inputDicts, poolSize, log, **kwargs):
    """*run ``function`` for each input in a pool of spawned worker processes, one task per input*

    Each task is its own future, so a worker killed by a signal (a segfault, an out-of-memory kill) does not lose the
    whole batch. Workers are spawned, not forked, because forked workers crash in Apple Accelerate (GCD is not
    fork-safe).

    A pool that breaks fails every unfinished future the same way, so each task first touches a marker file. A
    marker means the function body was entered. No marker means the function never ran, so the input is safe to run
    again. After a break:

    - inputs with no marker are run again in a fresh pool
    - a single started input is the one that crashed, and fails with ``WORKER_CRASH_MESSAGE``
    - several started inputs are each run once more, alone, in a pool of one worker. The input that kills its
      worker there fails with ``WORKER_CRASH_MESSAGE``, and the others keep their results

    Every round settles at least one input or ends the batch, and an input is run in isolation at most once, so the
    loop ends. When no input can start (for example a worker cannot unpickle ``function``), the remaining inputs fail
    with ``WORKER_START_MESSAGE``. When ``POOL_TIMEOUT_SECONDS`` passes, the workers are stopped and the unfinished
    inputs fail with ``WORKER_TIMEOUT_MESSAGE``. An interrupt stops the workers and is raised again.

    A SOF that is run again after a crash may meet partial products written by its first, killed run. The
    function must classify those itself, as ``_bulk_worker`` does with its ``FileExistsError`` handling.

    **Key Arguments:**

    - ``function`` -- a module-level function, called as ``function(inputDict, log=log, **kwargs)``. It must return
      a dictionary and may not raise. The function, ``log``, ``kwargs`` and the inputs are sent to the spawned
      workers, so all of them must be picklable
    - ``inputDicts`` -- the list of input dictionaries, each with a ``sof`` key
    - ``poolSize`` -- the number of worker processes
    - ``log`` -- logger
    - ``kwargs`` -- further keyword arguments for ``function``

    **Return:**

    - ``results`` -- one result per input, in the order of ``inputDicts``
    """
    import tempfile

    from tqdm import tqdm

    results = [None] * len(inputDicts)
    pending = list(range(len(inputDicts)))
    deadline = time.monotonic() + POOL_TIMEOUT_SECONDS
    with tempfile.TemporaryDirectory() as markerDirectory, tqdm(total=len(inputDicts)) as progressBar:
        task = _pool_task(function, inputDicts, log, kwargs, markerDirectory)
        while pending:
            if deadline - time.monotonic() <= 0:
                resolved = {index: _timed_out_result(task, index) for index in pending}
                unstarted = []
            else:
                outcome = _run_pool_round(task, pending, poolSize, deadline - time.monotonic())
                resolved = {**outcome.resolved, **_settle_in_flight(task, outcome.inFlight, deadline)}
                unstarted = outcome.unstarted
                if not resolved and len(unstarted) == len(pending):
                    # NO TASK STARTED, SO A NEW POOL WOULD BREAK THE SAME WAY
                    log.error(WORKER_START_MESSAGE)
                    resolved = {index: _failed_result(inputDicts[index], WORKER_START_MESSAGE) for index in pending}
                    unstarted = []
            for index, result in resolved.items():
                results[index] = result
            progressBar.update(len(resolved))
            pending = unstarted
    return results


def _crashed_result(task, index):
    """*log and return the failure of an input whose worker process was killed*"""
    sof = task.inputDicts[index]["sof"]
    task.log.error(f"{WORKER_CRASH_MESSAGE} SOF: {sof}")
    return _failed_result(task.inputDicts[index], WORKER_CRASH_MESSAGE)


def _timed_out_result(task, index):
    """*log and return the failure of an input that the batch deadline caught*"""
    task.log.error(f"{WORKER_TIMEOUT_MESSAGE} SOF: {task.inputDicts[index]['sof']}")
    return _failed_result(task.inputDicts[index], WORKER_TIMEOUT_MESSAGE)


def _settle_in_flight(task, inFlight, deadline):
    """*decide which started inputs crashed the pool, running each of several once more on its own*

    A lone started input is taken to be the one that crashed, and is not run again. That is wrong only if a
    different worker died before it touched its marker. Several started inputs are each run alone in a pool of one.
    Once the deadline has passed, the input and all that follow fail as timed out and no pool is started.
    """
    if len(inFlight) == 1:
        return {inFlight[0]: _crashed_result(task, inFlight[0])}
    settled = {}
    for index in inFlight:
        if deadline - time.monotonic() <= 0:
            settled[index] = _timed_out_result(task, index)
            continue
        outcome = _run_pool_round(task, [index], 1, deadline - time.monotonic())
        if index in outcome.resolved:
            settled[index] = outcome.resolved[index]
        elif outcome.inFlight:
            settled[index] = _crashed_result(task, index)
        else:
            settled[index] = _failed_result(task.inputDicts[index], WORKER_START_MESSAGE)
    return settled


def _run_pool_round(task, indexes, poolSize, timeoutSeconds):
    """*run the inputs at ``indexes`` in one fresh pool and return what the round settled*"""
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor

    for index in indexes:
        # A MARKER LEFT BY AN EARLIER ROUND MUST NOT MAKE A NEW ROUND LOOK STARTED
        if os.path.exists(task.marker_path(index)):
            os.remove(task.marker_path(index))
    # FORKED WORKERS CRASH IN APPLE ACCELERATE (GCD IS NOT FORK-SAFE); SPAWN GIVES EACH WORKER A CLEAN PROCESS
    executor = ProcessPoolExecutor(
        max_workers=max(1, min(poolSize, len(indexes))),
        mp_context=multiprocessing.get_context("spawn"),
        initializer=_mute_worker_stdout,
    )
    try:
        futures, unsubmitted = _submit_all(executor, task, indexes)
        outcome = _collect_round(executor, task, futures, unsubmitted, timeoutSeconds)
        executor.shutdown()
    except BaseException:
        # AN INTERRUPT, EVEN DURING THE SHUTDOWN, MUST NOT LEAVE WORKERS RUNNING AFTER THE MARKER DIRECTORY IS DELETED
        _stop_workers(executor)
        raise
    return outcome


def _submit_all(executor, task, indexes):
    """*submit one future per input, and return the futures with the indexes the broken pool refused*"""
    from concurrent.futures.process import BrokenProcessPool

    futures = {}
    for index in indexes:
        try:
            future = executor.submit(
                _run_marked, task.marker_path(index), task.function, task.inputDicts[index], task.log, task.kwargs
            )
        except BrokenProcessPool:
            break
        futures[future] = index
    submitted = set(futures.values())
    return futures, [index for index in indexes if index not in submitted]


def _collect_round(executor, task, futures, unsubmitted, timeoutSeconds):
    """*gather the futures of one round into a ``_round_outcome``, stopping the workers on a timeout*"""
    from concurrent.futures import as_completed
    from concurrent.futures.process import BrokenProcessPool

    resolved, inFlight, unstarted = {}, [], list(unsubmitted)
    try:
        for future in as_completed(futures, timeout=timeoutSeconds):
            index = futures[future]
            try:
                resolved[index] = future.result()
            except BrokenProcessPool:
                (inFlight if os.path.exists(task.marker_path(index)) else unstarted).append(index)
            except Exception as error:
                task.log.error(f"The reduction of {task.inputDicts[index]['sof']} returned no result: {error}")
                resolved[index] = _failed_result(task.inputDicts[index], str(error))
    except TimeoutError:
        _stop_workers(executor)
        settled = set(resolved) | set(inFlight) | set(unstarted)
        for index in futures.values():
            if index not in settled:
                resolved[index] = _timed_out_result(task, index)
    return _round_outcome(resolved, sorted(unstarted), sorted(inFlight))


def _stop_workers(executor):
    """*terminate and join every worker process of an executor, then shut the executor down*"""
    # THE EXECUTOR HAS NO PUBLIC WAY TO STOP RUNNING WORKERS, SO THIS READS ITS PRIVATE _processes MAPPING
    processes = list((getattr(executor, "_processes", None) or {}).values())
    for process in processes:
        process.terminate()
    for process in processes:
        process.join(WORKER_STOP_SECONDS)
        if process.is_alive():
            process.kill()
            process.join()
    executor.shutdown(wait=True, cancel_futures=True)


def run_recipe_bulk(
    log,
    recipe,
    sofList,
    commandList,
    settings,
    overwrite,
    workspaceDirectory,
    conn,
    sessionId,
    poolSize=None,
):
    """*execute a pipeline recipe in multiprocessing mode*

    **Key Arguments:**

    - ``log`` -- logger
    - ``recipe`` -- the name of the recipe tp execute
    - ``sofList`` -- a list of paths to the sof files containing the files the recipe requires
    - ``commandList`` -- a list of the commands used to run the recipe for each
    - ``settings`` -- the settings dictionary
    - ``overwrite`` -- overwrite existing reductions. Default *False*.
    - ``workspaceDirectory`` -- path to the root of the workspace
    - ``conn`` -- a connection to the workspace database
    - ``sessionId`` -- the session ID of the workspace
    - ``poolSize`` -- the number of worker processes. Default *None*, which uses 3 for the memory-heavy recipes
      (mflat, stare, nod and offset) and one per CPU for every other recipe.

    The reductions run in a pool of spawned worker processes, one task per SOF (see ``_run_in_spawn_pool``).
    A SOF whose worker process crashes is recorded as a failure and does not stop the rest of the batch.
    The workers send their log records to this process, which writes them through its own log handlers (see
    ``_worker_log_queue``), so only one process writes and rotates the log file.
    A ``ValueError`` is raised if ``settings`` has no ``logging settings`` block, because spawned workers
    take their log level from it.
    """
    log.debug("starting the ``run_recipe_bulk`` method")

    if "logging settings" not in settings:
        raise ValueError("The settings have no `logging settings` block, which spawned reduction workers need to log.")

    import shutil

    import pandas as pd

    from soxspipe.commonutils import data_organiser
    from soxspipe.commonutils.sql_identifiers import validate_sql_identifier

    inputDicts = [{"sof": sof, "command": command} for sof, command in zip(sofList, commandList)]

    wrapperTurnOffMP = True

    runningMessage = f"Running {len(inputDicts)} reductions for the {recipe.upper()} recipe in multiprocessing mode"
    if poolSize is not None:
        effectivePoolSize = poolSize
        log.print(f"{runningMessage} with a pool size of {effectivePoolSize}...")
    elif "mflat" in recipe or "stare" in recipe or "nod" in recipe or "offset" in recipe:
        effectivePoolSize = 3
        log.print(f"{runningMessage} with a pool size of {effectivePoolSize} to avoid memory issues...")
    else:
        effectivePoolSize = os.cpu_count() or 1
        log.print(f"{runningMessage} with a pool size of {os.cpu_count()}...")

    with _worker_log_queue(log) as logQueue:
        results = _run_in_spawn_pool(
            function=_bulk_worker,
            inputDicts=inputDicts,
            poolSize=effectivePoolSize,
            log=log,
            settings=settings,
            overwrite=overwrite,
            recipe=recipe,
            workspaceDirectory=workspaceDirectory,
            logQueue=logQueue,
            wrapperTurnOffMP=wrapperTurnOffMP,
        )

    passing = []
    failing = []
    skipped = []
    for result in results:
        sof = os.path.basename(result["sof"])
        sofList.append(sof)
        if result["status"] == "pass":
            passing.append(sof)
        elif result["status"] == "fail":
            failing.append(sof)
        else:
            skipped.append(sof)

    print(
        f"Number of successful {recipe} reductions: {len(passing)}. Number of failed {recipe} reductions: {len(failing)}. Number of pre-existing {recipe} reductions: {len(skipped)}.\n"
    )

    ## COLLECT TOGETHER THE RESULTS AND UPDATE THE DATABASE
    passing = []
    failing = []
    skipped = []
    qcTables = []
    sofList = []
    errorMessages = []
    errorSOF = []
    for result in results:
        sof = os.path.basename(result["sof"])
        sofList.append(sof)
        if result["status"] == "pass":
            passing.append(sof)
            qcTables.append(result["qcTable"])
        elif result["status"] == "previous-pass":
            passing.append(sof)
        elif result["status"] == "previous-fail":
            failing.append(sof)
        elif result["status"] == "fail":
            failing.append(sof)
            errorSOF.append(sof)
            errorMessages.append(result["error_message"])
        if result["qcTable"] is not None:
            qcTables.append(result["qcTable"])

    c = conn.cursor()
    statusColumn = validate_sql_identifier(f"status_{sessionId}", "status column")
    passingClause, passingParams = _sql_in_placeholders(passing)
    failingClause, failingParams = _sql_in_placeholders(failing)
    sqlQuery = (
        "select count(*) from product_frames where "  # noqa: S608
        f"(sof in ({passingClause}) and {statusColumn} = 'fail') "
        f"or (sof in ({failingClause}) and ({statusColumn} != 'fail' or {statusColumn} is null))"
    )
    c.execute(sqlQuery, passingParams + failingParams)
    count = c.fetchone()[0]
    if len(passing) or len(failing):
        sqlQueries = [
            (
                f"update product_frames set {statusColumn} = 'pass' where sof in ({passingClause}) "  # noqa: S608
                f"and ({statusColumn} != 'pass' or {statusColumn} is null)",
                passingParams,
            ),
            (
                f"update product_frames set {statusColumn} = 'fail' where sof in ({failingClause}) "  # noqa: S608
                f"and ({statusColumn} != 'fail' or {statusColumn} is null)",
                failingParams,
            ),
        ]
        for sqlQuery, sqlParams in sqlQueries:
            c.execute(sqlQuery, sqlParams)
    if len(errorSOF):
        error_rows = [(str(error), sof) for sof, error in zip(errorSOF, errorMessages)]
        c.executemany(
            "update product_frames set error_message = ? where sof = ?",
            error_rows,
        )
    c.close()
    conn.close()

    if count:
        do = data_organiser(log=log, rootDir=workspaceDirectory)
        reset = do.session_refresh(failure=True)
        do.close()

    if len(qcTables):
        qcTables = pd.concat(qcTables, ignore_index=True)
        do = data_organiser(log=log, rootDir=workspaceDirectory)
        do._dataframe_to_sqlite(qcTables, "quality_control", replace=False)
        do.close()

    exists = os.path.exists(workspaceDirectory + "/tmp/")
    if exists:
        shutil.rmtree(workspaceDirectory + "/tmp/")

    log.debug("completed the ``run_recipe_bulk`` method")
    return
