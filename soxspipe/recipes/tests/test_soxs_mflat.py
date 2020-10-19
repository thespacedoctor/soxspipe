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

    def test_soxs_mflat_nir_normalise_function(self):

        # UNPACK THE ORDER TABLE
        from soxspipe.commonutils.toolkit import unpack_order_table
        orders, orderCentres, orderLimits, orderEdgeLow, orderEdgeUp = unpack_order_table(
            log=log, orderTablePath="/Users/Dave/soxspipe-unittests/intermediate/order_locations_NIR.csv")
        print(orders, orderCentres, orderLimits, orderEdgeLow, orderEdgeUp)

        # sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/nir_long_flats.sof"
        # from soxspipe.recipes import soxs_mflat
        # this = soxs_mflat(
        #     log=log,
        #     settings=settings,
        #     inputFrames=sofPath
        # )
        # this.produce_product()

        # # inputFiles = "/path/to/a/directory"
        # # inputFiles = ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
        # inputFiles = '~/xshooter-pipeline-data/unittest_data/xshooter-flats/nir/calibrated/'
        # from soxspipe.commonutils import set_of_files
        # sof = set_of_files(
        #     log=log,
        #     settings=settings,
        #     inputFrames=inputFiles
        # )
        # sofFile, supplementarySof = sof.get()
        # # LIST OF CCDDATA OBJECTS
        # ccds = [c for c in sofFile.ccds(ccd_kwargs={
        #                                 "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        # this.normalise_flats(
        #     inputFlats=ccds,
        #     orderTablePath=supplementarySof["NIR"]["ORDER_CENT"])

        # this.produce_product()

    def test_soxs_mflat_nir_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/nir_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_nir_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/nir_short_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_uvb_dflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/uvb_dflats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_uvb_qflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/uvb_qflats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_vis_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/vis_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_vis_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof/vis_short_flats.sof"
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
            sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-flats/sof//nir_mixed_exptime_darks.sof"
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
