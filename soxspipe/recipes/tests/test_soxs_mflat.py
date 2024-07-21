from __future__ import print_function
import numpy.ma as ma
import math
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData
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

packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings_xsh.yaml"
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

packageDirectory = utKit("").get_project_root()
settingsFile2 = packageDirectory + "/test_settings_soxs_sim.yaml"
# settingsFile = home + "/.config/soxspipe/soxspipe.yaml"
su = tools(
    arguments={"settingsFile": settingsFile2},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False
)
arguments2, settings2, log2, dbConn2 = su.setup()

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

    import pytest

    @pytest.mark.full
    def test_xsh_unpack_order_table_function(self):
        orderTablePath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/nir/20170818T173106_NIR_ORDER_LOCATIONS.fits"
        # UNPACK THE ORDER TABLE
        from soxspipe.commonutils.toolkit import unpack_order_table
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=log, orderTablePath=orderTablePath)

        from tabulate import tabulate
        print(tabulate(orderPolyTable, headers='keys', tablefmt='psql'))
        print(tabulate(orderPixelTable.head(100), headers='keys', tablefmt='psql'))

    def test_xsh_mflat_nir_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/nir_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        mflat = this.produce_product()
        print(f"The master flat file has been saved to '{mflat}'")

    @pytest.mark.full
    def test_soxs_mflat_nir_soxsreal_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/soxs/FLAT/sof/SOXS_NIR_FLATS.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log2,
            settings=settings2,
            inputFrames=sofPath,
            overwrite=True
        )
        mflat = this.produce_product()
        print(f"The master flat file has been saved to '{mflat}'")

    @pytest.mark.full
    def test_xsh_mflat_nir_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/nir_short_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        this.produce_product()

    def test_xsh_mflat_uvb_dflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/uvb_dflats.sof"
        from soxspipe.recipes import soxs_mflat
        try:
            this = soxs_mflat(
                log=log,
                settings=settings,
                inputFrames=sofPath,
                overwrite=True
            )
            this.produce_product()
            assert False
        except Exception as e:
            assert True
            print(str(e))

    @pytest.mark.full
    def test_xsh_mflat_uvb_qflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/uvb_qflats.sof"
        from soxspipe.recipes import soxs_mflat
        try:
            this = soxs_mflat(
                log=log,
                settings=settings,
                inputFrames=sofPath,
                overwrite=True
            )
            this.produce_product()
            assert False
        except Exception as e:
            assert True
            print(str(e))

    def test_xsh_mflat_vis_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/vis_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        this.produce_product()

    # def test_xsh_mflat_vis_short_function(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof/vis_short_flats.sof"
    #     from soxspipe.recipes import soxs_mflat
    #     this = soxs_mflat(
    #         log=log,
    #         settings=settings,
    #         inputFrames=sofPath
    #     )
    #     this.produce_product()

    def test_soxs_mflat_function_exception(self):

        from soxspipe.recipes import soxs_mflat
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mflat/sof//nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_mflat
            this = soxs_mflat(
                log=log,
                settings=settings,
                inputFrames=sofPath,
                overwrite=True
            )
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
