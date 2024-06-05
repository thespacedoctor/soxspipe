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
import sys
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


class test_soxs_nod(unittest.TestCase):
    def test_soxs_nod_nir_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-nod/sof/2010.09.14T23.39.32.2932_NIR_NOD_20PT0_XSHOOTER.sof"
        from soxspipe.recipes import soxs_nod
        this = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_nod_uvb_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-nod/sof/2010.09.14T23.39.21.697_UVB_1X1_FAST_NOD_XSHOOTER.sof"
        from soxspipe.recipes import soxs_nod
        this = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_nod_vis_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-nod/sof/2010.09.14T23.39.26.928_VIS_1X1_FAST_NOD_XSHOOTER.sof"
        from soxspipe.recipes import soxs_nod
        this = soxs_nod(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_nod_function_exception(self):
        from soxspipe.recipes import soxs_nod
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xsh/SOMEDIRECTORY/sofs/nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_nod
            this = soxs_nod(
                log=log,
                settings=settings,
                inputFrames=sofPath
            )
        except Exception as e:
            self.fail(f"Exception raised: {str(e)}")

if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(test_soxs_nod("test_soxs_nod_nir_function"))
    runner = unittest.TextTestRunner()
    runner.run(suite)

    #unittest.main()