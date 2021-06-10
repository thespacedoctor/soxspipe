from __future__ import print_function
from builtins import str
import os
import unittest
import shutil
import unittest
import yaml
from soxspipe.utKit import utKit
from fundamentals import tools
from os.path import expanduser
home = expanduser("~")
from astropy.nddata import CCDData
from astropy import units as u
import numpy as np
import math
import numpy.ma as ma

packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings.yaml"
# settingsFile = home + "/.config/soxspipe/soxspipe.yaml"
su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False
)
arguments, settings, log, dbConn = su.setup()

# SETUP AND TEARDOWN FIXTURE FUNCTIONS FOR THE ENTIRE MODULE
moduleDirectory = os.path.dirname(__file__)
utKit = utKit(moduleDirectory)
log, dbConn, pathToInputDir, pathToOutputDir = utKit.setupModule()
utKit.tearDownModule()

try:
    shutil.rmtree(pathToOutputDir)
except:
    pass
# COPY INPUT TO OUTPUT DIR
shutil.copytree(pathToInputDir, pathToOutputDir)

# Recursively create missing directories
if not os.path.exists(pathToOutputDir):
    os.makedirs(pathToOutputDir)

# xt-setup-unit-testing-files-and-folders


class test_soxs_mflat(unittest.TestCase):

    def test_unpack_order_table_function(self):
        orderTablePath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/nir/order_locations_NIR_locations.csv"
        # UNPACK THE ORDER TABLE
        from soxspipe.commonutils.toolkit import unpack_order_table
        orderPolyTable, orderPixelTable = unpack_order_table(
            log=log, orderTablePath=orderTablePath)

        from tabulate import tabulate
        print(tabulate(orderPolyTable, headers='keys', tablefmt='psql'))
        print(tabulate(orderPixelTable.head(100), headers='keys', tablefmt='psql'))

    def test_soxs_mflat_nir_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/nir_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        mflat = this.produce_product()
        print(f"The master flat file has been saved to '{mflat}'")

    def test_soxs_mflat_nir_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/nir_short_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_uvb_dflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/uvb_dflats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    # def test_soxs_mflat_uvb_qflat_function(self):

    #     sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/uvb_qflats.sof"
    #     from soxspipe.recipes import soxs_mflat
    #     this = soxs_mflat(
    #         log=log,
    #         settings=settings,
    #         inputFrames=sofPath
    #     )
    #     this.produce_product()

    def test_soxs_mflat_vis_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/vis_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_vis_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/vis_short_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    # def test_soxs_mflat_function(self):

    #     # utKit.refresh_database() # reset database to database found in
    #     # soxspipe/test/input
    #     from soxspipe.recipes import soxs_mflat
    #     this = soxs_mflat(
    #         log=log,
    #         settings=settings
    #     )
    #     this.get()

    def test_soxs_mflat_function_exception(self):

        from soxspipe.recipes import soxs_mflat
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof//nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_mflat
            this = soxs_mflat(
                log=log,
                settings=settings,
                inputFrames=sofPath
            )
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
