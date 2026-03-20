#!/usr/bin/env python
# encoding: utf-8
"""
*reduce all the data in a workspace, or target specific obs and files for reduction*

Author
: David Young

Date Created
: January 17, 2024
"""
# from memory_profiler import profile
from fundamentals import tools
from builtins import object
import sys
import os

os.environ["TERM"] = "vt100"


class reducer(object):
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
            return None

        self.recipeList = [
            "mbias",
            "mdark",
            "disp_solution",
            "order_centres",
            "mflat",
            "spat_solution",
            "nod",
            "stare",
            "offset",
            "nod-obj",
            "stare-obj",
            "offset-obj",
        ]

        self.sessionPath = workspaceDirectory + "/sessions/" + self.sessionId
        self.sessionDB = workspaceDirectory + "/soxspipe.db"

        if refreshWorkspace:
            self.log.print("Refreshing the workspace before reduction...")
            do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
            do.prepare(refresh=False, report=False)
            do.close()

        return None

    def reduce(self, batch=False, multiprocess=False):
        """
        *reduce the selected data*
        """
        self.log.debug("starting the ``reduce`` method")

        if self.sessionId is None:
            print("Please prepare this workspace using `soxspipe prep` before attempting to reduce the data.")
            return None

        from fundamentals import times
        import traceback
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

                    conn = sql.connect(self.sessionDB, timeout=300, autocommit=True, check_same_thread=False)
                    c = conn.cursor()
                    c.execute("PRAGMA busy_timeout = 100000")
                    c.execute("PRAGMA synchronous = OFF")

                    sofList = rawGroups["sof"].tolist()
                    sofList = [self.sessionPath + "/sof/" + sof for sof in sofList]
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
                    )
                    break
                else:

                    fail = False
                    for index, row in rawGroups.iterrows():
                        if batchCount >= batch:
                            self.log.print(f"Batch limit of {batch} reached, pausing reductions.")
                            break

                        recipe = row["recipe"].replace("-std", "").replace("-obj", "")
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
                        except FileExistsError as e:
                            continue
                        except Exception as e:
                            # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
                            self.log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
                            self.log.error(f'\nRecipe Command: {row["command"].replace("-obj ", " ")}\n\n')
                            fail = True

                            if self.quitOnFail:
                                sys.exit(0)

                            if self.reductionTarget != "all":
                                self.overwrite = False

                            if not self.daemon:
                                print(f"{'='*70}\n")

                        ## FINISH LOGGING ##
                        endTime = times.get_now_sql_datetime()
                        runningTime = times.calculate_time_difference(startTime, endTime)
                        sys.argv[0] = os.path.basename(sys.argv[0])

                        self.log.print(f'\nRecipe Command: {row["command"].replace("-obj ", " ")} ')
                        self.log.print(f"Recipe Run Time: {runningTime}\n\n")
                        if not self.daemon:
                            print(f"{'='*70}\n")

                    if fail:
                        do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
                        reset = do.session_refresh()
                        do.close()
                        if reset:
                            print(f"BACK TO THE START! {rootRecipe}\n\n")
                            break
                    break

        if self.reductionTarget == "all":
            do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
            incompleteSets = do.get_incomplete_raw_frames_set()
            do.close()
            if len(incompleteSets.index):
                from tabulate import tabulate

                print(
                    "\nSOME CALIBRATION FRAMES ARE NOT PRESENT (OR FAILED TO BE BUILT) FOR THE FOLLOWING DATA SETS AND THEY CANNOT BE REDUCED:"
                )
                print(tabulate(incompleteSets, headers="keys", tablefmt="psql", showindex=False))

        self.log.debug("completed the ``reduce`` method")
        return None

    def select_sof_files_to_process(self, recipe=False, reductionTarget=False, batch=False, arm=False):
        """*select all of the SOF files still requiring processing*

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

        import pandas as pd
        import sqlite3 as sql

        conn = sql.connect(self.sessionDB, timeout=300, autocommit=True, check_same_thread=False)
        c = conn.cursor()
        c.execute("PRAGMA busy_timeout = 100000")
        c.execute("PRAGMA synchronous = OFF")

        if batch:
            limitText = f" LIMIT {batch} "
        else:
            limitText = ""

        if arm:
            armText = f" and `eso seq arm` = '{arm}' "
        else:
            armText = ""

        if reductionTarget == "all":
            # GET THE GROUPS OF FILES NEEDING REDUCED, ASSIGN THE CORRECT COMMAND TO EXECUTE THE RECIPE
            if not recipe:
                recipeText = "is not null"
            else:
                recipeText = f"= '{recipe}'"

            rawGroups = pd.read_sql(
                f"SELECT * FROM raw_frame_sets where recipe_order is not null and complete = 1 and recipe {recipeText} {armText}  order by recipe_order, sof {limitText}",
                con=conn,
            )

        elif reductionTarget.split(".")[-1].lower() == "sof":
            sqlQuery = f"select sof from product_frames where sof = '{reductionTarget}' and complete = 1"

            for _ in range(4):  # Recursively query up to 5 times
                sqlQuery = f"SELECT distinct sof FROM product_frames WHERE file IN (SELECT file FROM sof_map_base WHERE sof in ({sqlQuery})) or sof in ({sqlQuery})"

            sqlQuery = (
                f"SELECT distinct sof, recipe from raw_frame_sets WHERE sof in ({sqlQuery}) order by recipe_order, sof"
            )

            rawGroups = pd.read_sql(sqlQuery, con=conn)

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

        rawGroups["command"] = "soxspipe " + rawGroups["recipe"].str.replace("-obj", "") + " sof/" + rawGroups["sof"]
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

    productPath, qcTable = soxs_recipe.produce_product()
    del soxs_recipe

    log.debug("completed the ``run_recipe`` method")
    return productPath, qcTable


def run_recipe_bulk(log, recipe, sofList, commandList, settings, overwrite, workspaceDirectory, conn, sessionId):
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
    """
    log.debug("starting the ``run_recipe_bulk`` method")

    from fundamentals import fmultiprocess
    from soxspipe.commonutils import data_organiser
    import pandas as pd
    import shutil

    def wrapper(inputDict, log, recipe, settings, overwrite, workspaceDirectory, wrapperTurnOffMP=True):
        import traceback

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
                returnDict["error_message"] = e
        except Exception as e:
            # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
            log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
            log.error(f'\nRecipe Command: {inputDict["command"].replace("-obj ", " ")}\n\n')
            returnDict["status"] = "fail"
            returnDict["error_message"] = e

        return returnDict

    inputDicts = [{"sof": sof, "command": command} for sof, command in zip(sofList, commandList)]

    poolSize = False
    if ("nod" in recipe or "stare" in recipe) and len(inputDicts) < 4:
        turnOffMP = True
        wrapperTurnOffMP = False
    else:
        turnOffMP = False
        wrapperTurnOffMP = True

    if "mflat" in recipe:
        poolSize = 3
        print(
            f"Running {len(inputDicts)} reductions for the {recipe.upper()} recipe in multiprocessing mode with a pool size of {poolSize} to avoid memory issues..."
        )

    log.print(f"Running {len(inputDicts)} reductions for the {recipe.upper()} recipe in multiprocessing mode...")
    results = fmultiprocess(
        log=log,
        function=wrapper,
        inputArray=inputDicts,
        poolSize=poolSize,
        timeout=36000,
        settings=settings,
        overwrite=overwrite,
        recipe=recipe,
        workspaceDirectory=workspaceDirectory,
        wrapperTurnOffMP=wrapperTurnOffMP,
        turnOffMP=turnOffMP,
        mute=True,
        progressBar=True,
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
    passingString = "','".join(passing)
    failingString = "','".join(failing)
    sqlQuery = f"select count(*) from product_frames where (sof in ('{passingString}') and status_{sessionId} = 'fail') or  (sof in ('{failingString}') and (status_{sessionId} != 'fail' or status_{sessionId} is null))"
    c.execute(sqlQuery)
    count = c.fetchone()[0]
    if len(passing) or len(failing):
        sqlQueries = [
            f"update product_frames set status_{sessionId} = 'pass' where sof in ('{passingString}') and (status_{sessionId} != 'pass' or status_{sessionId} is null)",
            f"update product_frames set status_{sessionId} = 'fail' where sof in ('{failingString}') and (status_{sessionId} != 'fail' or status_{sessionId} is null)",
        ]
        for sqlQuery in sqlQueries:
            c.execute(sqlQuery)
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
    return None
