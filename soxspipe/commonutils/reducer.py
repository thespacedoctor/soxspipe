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

        # REQUEST THE WORKSPACE PARAMETERS FROM THE DATA-ORGANISER
        from soxspipe.commonutils import data_organiser

        do = data_organiser(
            log=log,
            rootDir=workspaceDirectory,
        )
        self.sessionId, allSessions = do.session_list(silent=True)

        if self.sessionId is None:
            return None

        self.recipeList = [
            "mbias",
            "mdark",
            "disp_sol",
            "order_centres",
            "mflat",
            "spat_sol",
            "nod-std",
            "stare-std",
            "offset-std",
            "nod",
            "stare",
            "offset",
        ]

        self.sessionPath = workspaceDirectory + "/sessions/" + self.sessionId
        self.sessionDB = workspaceDirectory + "/soxspipe.db"

        return None

    # @profile
    def reduce(self):
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

        i = 0

        if self.reductionTarget != "all":
            self.recipeList = [False]

        for recipe in self.recipeList:
            # rawGroups WILL CONTAIN ONE RECIPE COMMAND PER ENTRY
            rawGroups = self.select_sof_files_to_process(recipe=recipe, reductionTarget=self.reductionTarget)

            for index, row in rawGroups.iterrows():
                recipe = row["recipe"]
                sof = row["sof"]
                startTime = times.get_now_sql_datetime()
                sof = self.sessionPath + "/sof/" + sof

                try:
                    run_recipe(
                        self.log, recipe, sof, settings=self.settings, overwrite=self.overwrite, command=row["command"]
                    )
                    i += 1
                except FileExistsError as e:
                    continue
                except Exception as e:
                    # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
                    self.log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
                    self.log.error(f'\nRecipe Command: {row["command"]}\n\n')

                    if self.quitOnFail:
                        sys.exit(0)

                    do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
                    do.session_refresh()
                    if not self.daemon:
                        print(f"{'='*70}\n")
                    continue

                ## FINISH LOGGING ##
                endTime = times.get_now_sql_datetime()
                runningTime = times.calculate_time_difference(startTime, endTime)
                sys.argv[0] = os.path.basename(sys.argv[0])

                self.log.print(f'\nRecipe Command: {row["command"]} ')
                self.log.print(f"Recipe Run Time: {runningTime}\n\n")
                if not self.daemon:
                    print(f"{'='*70}\n")

        if self.reductionTarget == "all":
            do = data_organiser(log=self.log, rootDir=self.workspaceDirectory)
            incompleteSets = do.get_incomplete_raw_frames_set()
            if len(incompleteSets.index):
                from tabulate import tabulate

                print(
                    "\nSOME CALIBRATION FRAMES ARE NOT PRESENT (OR FAILED TO BE BUILT) FOR THE FOLLOWING DATA SETS AND THEY CANNOT BE REDUCED:"
                )
                print(tabulate(incompleteSets, headers="keys", tablefmt="psql", showindex=False))

        self.log.debug("completed the ``reduce`` method")
        return None

    def select_sof_files_to_process(self, recipe=False, reductionTarget=False):
        """*select all of the SOF files still requiring processing*

        **Key Arguments:**
            - ``recipe`` -- the name of the recipe to filter by (optional)
            - ``reductionTarget`` -- target for reduction: "all", "sof", "ob" (default: False)

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

        conn = sql.connect(self.sessionDB)

        if reductionTarget == "all":
            # GET THE GROUPS OF FILES NEEDING REDUCED, ASSIGN THE CORRECT COMMAND TO EXECUTE THE RECIPE
            if not recipe:
                recipeText = "is not null"
                std = ""
            elif "std" in recipe:
                recipeText = f"= '{recipe.replace("-std","")}'"
                std = " AND `eso dpr type` like '%STD%' "
            else:
                recipeText = f"= '{recipe}'"
                std = " AND `eso dpr type` not like '%STD%' "
            rawGroups = pd.read_sql(
                f"SELECT * FROM raw_frame_sets where recipe_order is not null and complete = 1 and recipe {recipeText} {std} order by recipe_order, sof",
                con=conn,
            )

        elif reductionTarget.split(".")[-1].lower() == "sof":
            sqlQuery = f"'{reductionTarget}'"

            for _ in range(4):  # Recursively query up to 5 times
                sqlQuery = f"SELECT distinct sof FROM product_frames WHERE file IN (SELECT file FROM sof_map_base WHERE sof in ({sqlQuery})) or sof in ({sqlQuery})"

            sqlQuery = (
                f"SELECT distinct sof, recipe from raw_frame_sets WHERE sof in ({sqlQuery}) order by recipe_order, sof"
            )

            rawGroups = pd.read_sql(sqlQuery, con=conn)

        # FILTER DATA FRAME
        # FIRST CREATE THE MASK
        if recipe:
            mask = rawGroups["recipe"] == recipe.replace("-std", "")
            rawGroups = rawGroups.loc[mask]

        rawGroups["command"] = "soxspipe " + rawGroups["recipe"] + " sof/" + rawGroups["sof"]
        if self.pathToSettings:
            rawGroups["command"] += f" -s {self.pathToSettings}"
        conn.close()

        self.log.debug("completed the ``select_sof_files_to_process`` method")
        return rawGroups[["recipe", "sof", "command"]].drop_duplicates()


def run_recipe(log, recipe, sof, settings, overwrite, command=False):
    """*execute a pipeline recipe*

    **Key Arguments:**

    - ``recipe`` -- the name of the recipe tp execute
    - ``sof`` -- path to the sof file containing the files the recipe requires
    - ``command`` -- the command used to run the recipe

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
        )

    if recipe == "mdark":
        from soxspipe.recipes import soxs_mdark

        soxs_recipe = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "disp_sol":
        from soxspipe.recipes import soxs_disp_solution

        soxs_recipe = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "order_centres":
        from soxspipe.recipes import soxs_order_centres

        soxs_recipe = soxs_order_centres(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "mflat":
        from soxspipe.recipes import soxs_mflat

        soxs_recipe = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "spat_sol":
        from soxspipe.recipes import soxs_spatial_solution

        soxs_recipe = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "stare":
        from soxspipe.recipes import soxs_stare

        soxs_recipe = soxs_stare(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    if recipe == "nod":
        from soxspipe.recipes import soxs_nod

        soxs_recipe = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=sof,
            overwrite=overwrite,
            command=command,
        )

    soxs_recipe.produce_product()
    del soxs_recipe

    log.debug("completed the ``run_recipe`` method")
    return None

    # use the tab-trigger below for new method
    # xt-class-method
