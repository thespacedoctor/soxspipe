#!/usr/bin/env python
# encoding: utf-8
"""
*reduce all the data in a workspace, or target specific obs and files for reduction*

:Author:
    David Young

:Date Created:
    January 17, 2024
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
        - ``overwrite`` -- overwrite existing reductions. Default *False*.

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a reducer object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``reducer`` to documentation
        - create a blog post about what ``reducer`` does
    ```

    ```python
    usage code 
    ```

    """
    # Initialisation
    # 1. @flagged: what are the unique Attributes for each object? Add them
    # to __init__

    def __init__(
            self,
            log,
            workspaceDirectory,
            settings=False,
            pathToSettings=False,
            overwrite=False

    ):
        self.log = log
        log.debug("instantiating a new 'reducer' object")
        self.settings = settings
        self.workspaceDirectory = workspaceDirectory
        self.overwrite = overwrite
        self.pathToSettings = pathToSettings

        from soxspipe.commonutils import data_organiser
        do = data_organiser(
            log=log,
            rootDir=workspaceDirectory,
        )
        self.sessionId, allSessions = do.session_list(silent=True)
        self.sessionPath = workspaceDirectory + "/sessions/" + self.sessionId
        self.sessionDB = self.sessionPath + "/soxspipe.db"

        return None

    # 4. @flagged: what actions does each object have to be able to perform? Add them here
    # Method Attributes
    def reduce(self):
        """
        *reduce the selected data*

        **Return:**
            - ``reducer``

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code 
        ```
        """
        self.log.debug('starting the ``reduce`` method')

        from fundamentals import times

        rawGroups = self.select_sof_files_to_process()

        for index, row in rawGroups.iterrows():
            recipe = row["recipe"]
            sof = row["sof"]
            startTime = times.get_now_sql_datetime()
            try:
                self.run_recipe(recipe, sof)
            except FileExistsError as e:
                continue

            ## FINISH LOGGING ##
            endTime = times.get_now_sql_datetime()
            runningTime = times.calculate_time_difference(startTime, endTime)
            sys.argv[0] = os.path.basename(sys.argv[0])

            self.log.print(f'\nRecipe Command: {row["command"]} ')
            self.log.print(f'Recipe Run Time: {runningTime}\n\n')
            print(f"{'='*70}\n")

        self.log.debug('completed the ``reduce`` method')
        return reducer

    def select_sof_files_to_process(
            self):
        """*select all of the SOF files still requiring processing*

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
        self.log.debug('starting the ``select_sof_files_to_process`` method')

        import pandas as pd
        import sqlite3 as sql

        conn = sql.connect(
            self.sessionDB)

        rawGroups = pd.read_sql(
            'SELECT * FROM raw_frame_sets where recipe_order is not null order by recipe_order', con=conn)
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
        """*run a recipe on an sof *

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

        if recipe == "spat_sol":
            from soxspipe.recipes import soxs_spatial_solution
            disp_map, mapImage2D, res_plots = soxs_spatial_solution(
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

        if recipe == "stare":
            from soxspipe.recipes import soxs_stare
            recipe = soxs_stare(
                log=self.log,
                settings=self.settings,
                inputFrames=sof,
                overwrite=self.overwrite
            )
            reducedStare = recipe.produce_product()

        self.log.debug('completed the ``run_recipe`` method')
        return None

    # use the tab-trigger below for new method
    # xt-class-method

    # 5. @flagged: what actions of the base class(es) need ammending? ammend them here
    # Override Method Attributes
    # method-override-tmpx
