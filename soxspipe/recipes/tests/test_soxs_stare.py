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


class test_soxs_stare(unittest.TestCase):

    import pytest

    # def test_xsh_stare_nir_function(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xsh/SOMEDIRECTORY/sofs/nir_6s_darks.sof"
    #     from soxspipe.recipes import soxs_stare
    #     this = soxs_stare(
    #         log=log,
    #         settings=settings,
    #         inputFrames=sofPath
    #     )
    #     this.produce_product()

    # def test_xsh_stare_uvb_function(self):

    #     sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xsh/SOMEDIRECTORY/sofs/uvb_1x1_dark_3600s.sof"
    #     from soxspipe.recipes import soxs_stare
    #     this = soxs_stare(
    #         log=log,
    #         settings=settings,
    #         inputFrames=sofPath
    #     )
    #     this.produce_product()

    # def test_xsh_stare_vis_function(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xsh/SOMEDIRECTORY/sofs/vis_1x1_dark_3600s.sof"
    #     from soxspipe.recipes import soxs_stare
    #     this = soxs_stare(
    #         log=log,
    #         settings=settings,
    #         inputFrames=sofPath
    #     )
    #     this.produce_product()

    # def test_xsh_stare_function(self):

    #     # utKit.refresh_database() # reset database to database found in
    #     # soxspipe/test/input
    #     from soxspipe.recipes import soxs_stare
    #     this = soxs_stare(
    #         log=log,
    #         settings=settings
    #     )
    #     this.get()

    # @pytest.mark.full
    # def test_soxs_stare_function_soxsreal(self):

    #     sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/sky_subtraction/sof/soxsreal_pseudo_objects.sof"
    #     from soxspipe.recipes import soxs_stare
    #     recipe = soxs_stare(
    #         log=log2,
    #         settings=settings2,
    #         inputFrames=sofPath,
    #         overwrite=True
    #     )
    #     reducedStare = recipe.produce_product()

    def test_soxs_stare_function_exception(self):

        from soxspipe.recipes import soxs_stare
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xsh/SOMEDIRECTORY/sofs/nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_stare
            this = soxs_stare(
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
#ADD A MAIN HERE TO RUN THOSE TESTS
if __name__ == "__main__":
    unittest.main()