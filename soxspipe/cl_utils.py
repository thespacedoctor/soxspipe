#!/usr/bin/env python
# encoding: utf-8
"""
Documentation for soxspipe can be found here: http://soxspipe.readthedocs.org

Usage:
    soxspipe init
    soxspipe mbias <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>] 
    soxspipe mdark <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    soxspipe mflat <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    soxspipe disp_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    soxspipe order_centres <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]
    soxspipe spat_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile>]

Options:
    init                                   setup the soxspipe settings file for the first time
    mbias                                  the master bias recipe
    mdark                                  the master dark recipe
    mflat                                  the master flat recipe
    disp_sol                               the disp solution recipe
    order_centres                          the order centres recipe
    spat_sol                               the spatial solution recipe

    inputFrames                            path to a directory of frames or a set-of-files file

    -h, --help                             show this help message
    -v, --version                          show version
    -s, --settings <pathToSettingsFile>    the settings file
"""
################# GLOBAL IMPORTS ####################
import sys
import os
os.environ['TERM'] = 'vt100'
import readline
import glob
import pickle
from docopt import docopt
from fundamentals import tools, times
from subprocess import Popen, PIPE, STDOUT


def tab_complete(text, state):
    return (glob.glob(text + '*') + [None])[state]


def main(arguments=None):
    """
    *The main function used when `cl_utils.py` is run as a single script from the cl, or when installed as a cl command*
    """
    # setup the command-line util settings
    su = tools(
        arguments=arguments,
        docString=__doc__,
        logLevel="ERROR",
        options_first=False,
        projectName="soxspipe",
        defaultSettingsFile=True
    )
    arguments, settings, log, dbConn = su.setup()

    # ALIGN ASTROPY LOGGING LEVEL WITH SOXSPIPES
    try:
        from astropy import log as astrolog
        astrolog.setLevel(settings["logging settings"]["root"]["level"])
    except:
        pass

    # tab completion for raw_input
    readline.set_completer_delims(' \t\n;')
    readline.parse_and_bind("tab: complete")
    readline.set_completer(tab_complete)

    # UNPACK REMAINING CL ARGUMENTS USING `EXEC` TO SETUP THE VARIABLE NAMES
    # AUTOMATICALLY
    a = {}
    for arg, val in list(arguments.items()):
        if arg[0] == "-":
            varname = arg.replace("-", "") + "Flag"
        else:
            varname = arg.replace("<", "").replace(">", "")
        a[varname] = val
        if arg == "--dbConn":
            dbConn = val
            a["dbConn"] = val
        log.debug('%s = %s' % (varname, val,))

    ## START LOGGING ##
    startTime = times.get_now_sql_datetime()
    log.info(
        '--- STARTING TO RUN THE cl_utils.py AT %s' %
        (startTime,))

    # set options interactively if user requests
    if "interactiveFlag" in a and a["interactiveFlag"]:

        # load previous settings
        moduleDirectory = os.path.dirname(__file__) + "/resources"
        pathToPickleFile = "%(moduleDirectory)s/previousSettings.p" % locals()
        try:
            with open(pathToPickleFile):
                pass
            previousSettingsExist = True
        except:
            previousSettingsExist = False
        previousSettings = {}
        if previousSettingsExist:
            previousSettings = pickle.load(open(pathToPickleFile, "rb"))

        # x-raw-input
        # x-boolean-raw-input
        # x-raw-input-with-default-value-from-previous-settings

        # save the most recently used requests
        pickleMeObjects = []
        pickleMe = {}
        theseLocals = locals()
        for k in pickleMeObjects:
            pickleMe[k] = theseLocals[k]
        pickle.dump(pickleMe, open(pathToPickleFile, "wb"))

    # PACK UP SOME OF THE CL SWITCHES INTO SETTINGS DICTIONARY
    if a['outputDirectory']:
        settings["intermediate-data-root"] = a['outputDirectory']

    if a["init"]:
        from os.path import expanduser
        home = expanduser("~")
        filepath = home + "/.config/soxspipe/soxspipe.yaml"
        try:
            cmd = """open %(filepath)s""" % locals()
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        except:
            pass
        try:
            cmd = """start %(filepath)s""" % locals()
            p = Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True)
        except:
            pass
        return

    if a["mbias"]:
        from soxspipe.recipes import soxs_mbias
        recipe = soxs_mbias(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        )
        mbiasFrame = recipe.produce_product()
        print("You can find the master bias frame at `%(mbiasFrame)s`" % locals())

    if a["mdark"]:
        from soxspipe.recipes import soxs_mdark
        recipe = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        )
        mdarkFrame = recipe.produce_product()
        print("You can find the master bias frame at `%(mdarkFrame)s`" % locals())

    if a["disp_sol"]:
        from soxspipe.recipes import soxs_disp_solution
        disp_map = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        ).produce_product()
        print(f"\nSingle pinhole first guess dispersion map saved to: {disp_map}")

    if a["order_centres"]:
        from soxspipe.recipes import soxs_order_centres
        order_table = soxs_order_centres(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        ).produce_product()
        print(f"\nThe order centre locations have been saved to an order table: {order_table}")

    if a["spat_sol"]:
        from soxspipe.recipes import soxs_spatial_solution
        disp_map = soxs_spatial_solution(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        ).produce_product()
        print(f"\nFull 2D dispersion map saved to: {disp_map}")

    if a["mflat"]:
        from soxspipe.recipes import soxs_mflat
        recipe = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=a["inputFrames"]
        )
        mflatFrame = recipe.produce_product()
        print(f"\nThe master flat frame has been saved to: {mflatFrame}")

    # CALL FUNCTIONS/OBJECTS

    if "dbConn" in locals() and dbConn:
        dbConn.commit()
        dbConn.close()
    ## FINISH LOGGING ##
    endTime = times.get_now_sql_datetime()
    runningTime = times.calculate_time_difference(startTime, endTime)
    log.info('-- FINISHED ATTEMPT TO RUN THE cl_utils.py AT %s (RUNTIME: %s) --' %
             (endTime, runningTime, ))

    return


if __name__ == '__main__':
    main()
