#!/usr/bin/env python
# encoding: utf-8
"""
*reduce all the data in a workspace, or target specific obs and files for reduction*

Author
: David Young

Date Created
: January 17, 2024
"""
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class reducer(object):
    """
    *reduce all the data in a workspace, or target specific obs and files for reduction*

    **Key Arguments:**

    - ``log`` -- logger
    - ``workspaceDirectory`` -- path to the root of the workspace
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
            settings=False,
            pathToSettings=False,
            quitOnFail=False,
            overwrite=False,
            daemon=False

    ):
        self.log = log
        log.debug("instantiating a new 'reducer' object")
        self.settings = settings
        self.workspaceDirectory = workspaceDirectory
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

        self.sessionPath = workspaceDirectory + "/sessions/" + self.sessionId
        self.sessionDB = workspaceDirectory + "/soxspipe.db"

        return None

    def reduce(self):
        """
        *reduce the selected data*
        """
        self.log.debug('starting the ``reduce`` method')

        if self.sessionId is None:
            print("Please prepare this workspace using `soxspipe prep` before attempting to reduce the data.")
            return None

        from fundamentals import times
        import traceback

        # rawGroups WILL CONTAIN ONE RECIPE COMMAND PER ENTRY
        # THE ENTRIES ARE SELECTED IN THE ORDER THEY NEED TO RUN
        rawGroups = self.select_sof_files_to_process()

        for index, row in rawGroups.iterrows():
            recipe = row["recipe"]
            sof = row["sof"]
            startTime = times.get_now_sql_datetime()
            try:
                self.run_recipe(recipe, sof)
            except FileExistsError as e:
                continue
            except Exception as e:
                # ONE FAILURE RESET THE SOF FILES SO FUTURE RECIPES DON'T RELY ON FAILED PRODUCTS
                self.log.error(f"\n\nRecipe failed with the following error:\n\n{traceback.format_exc()}")
                self.log.error(f'\nRecipe Command: {row["command"]}\n\n')

                if self.quitOnFail:
                    sys.exit(0)

                from soxspipe.commonutils import data_organiser
                do = data_organiser(
                    log=self.log,
                    rootDir=self.workspaceDirectory
                )
                do.session_refresh()
                if not self.daemon:
                    print(f"{'='*70}\n")
                continue

            ## FINISH LOGGING ##
            endTime = times.get_now_sql_datetime()
            runningTime = times.calculate_time_difference(startTime, endTime)
            sys.argv[0] = os.path.basename(sys.argv[0])

            self.log.print(f'\nRecipe Command: {row["command"]} ')
            self.log.print(f'Recipe Run Time: {runningTime}\n\n')
            if not self.daemon:
                print(f"{'='*70}\n")

        self.log.debug('completed the ``reduce`` method')
        return None

    def select_sof_files_to_process(
            self):
        """*select all of the SOF files still requiring processing*

        **Key Arguments:**
            # -

        **Return:**

        - `rawGroups` -- a dataframe of the containing a list of recipes and sof file paths

        **Usage:**

        ```python
        rawGroups = reducer.select_sof_files_to_process()
        ```
        """
        self.log.debug('starting the ``select_sof_files_to_process`` method')

        import pandas as pd
        import sqlite3 as sql
        conn = sql.connect(
            self.sessionDB)

        # GET THE GROUPS OF FILES NEEDING REDUCED, ASSIGN THE CORRECT COMMAND TO EXECUTE THE RECIPE
        rawGroups = pd.read_sql(
            'SELECT * FROM raw_frame_sets where recipe_order is not null and complete = 1 order by recipe_order', con=conn)
        rawGroups["command"] = "soxspipe " + rawGroups["recipe"] + " sof/" + rawGroups["sof"]
        if self.pathToSettings:
            rawGroups["command"] += f" -s {self.pathToSettings}"
        conn.close()

        self.log.debug('completed the ``select_sof_files_to_process`` method')
        return rawGroups[['recipe', 'sof', 'command']].drop_duplicates()

    def run_recipe(
            self,
            recipe,
            sof):
        """*execute a pipeline recipe*

        **Key Arguments:**

        - ``recipe`` -- the name of the recipe tp execute
        - ``sof`` -- path to the sof file containing the files the recipe requires

        **Usage:**

        ```python
        reducer.run_recipe("mbias", "/path/to/sofs/my_bias_files.sof")
        ```
        """
        self.log.debug('starting the ``run_recipe`` method')

        sof = self.sessionPath + "/sof/" + sof

        if recipe == "mbias":
            from soxspipe.recipes import soxs_mbias
            recipe = soxs_mbias(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            mbiasFrame = recipe.produce_product()

        if recipe == "mdark":
            from soxspipe.recipes import soxs_mdark
            recipe = soxs_mdark(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            mdarkFrame = recipe.produce_product()

        if recipe == "disp_sol":
            from soxspipe.recipes import soxs_disp_solution
            disp_map = soxs_disp_solution(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            ).produce_product()

        if recipe == "order_centres":
            from soxspipe.recipes import soxs_order_centres
            order_table = soxs_order_centres(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            ).produce_product()

        if recipe == "mflat":
            from soxspipe.recipes import soxs_mflat
            recipe = soxs_mflat(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            mflatFrame = recipe.produce_product()

        if recipe == "spat_sol":
            from soxspipe.recipes import soxs_spatial_solution
            disp_map, mapImage2D, res_plots = soxs_spatial_solution(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            ).produce_product()

        if recipe == "stare":
            from soxspipe.recipes import soxs_stare
            recipe = soxs_stare(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            reducedStare = recipe.produce_product()

        if recipe == "nod":
            from soxspipe.recipes import soxs_nod
            recipe = soxs_nod(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            reducedNod = recipe.produce_product()

        self.log.debug('completed the ``run_recipe`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method
