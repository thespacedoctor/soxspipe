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


class test_soxs_disp_solution(unittest.TestCase):

    def test_soxs_disp_solution_nir_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-disp-solution/sof/20170818_NIR_DISP_SOLUTION.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    def test_soxs_disp_solution_uvb_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-disp-solution/sof/20170818_UVB_DISP_SOLUTION_1x1_fast.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    def test_soxs_disp_solution_vis_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-disp-solution/sof/20170818_VIS_DISP_SOLUTION_1x1_fast.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    # def test_soxs_disp_solution_function(self):

    #     # utKit.refresh_database() # reset database to database found in
    #     # soxspipe/test/input
    #     from soxspipe.recipes import soxs_disp_solution
    #     this = soxs_disp_solution(
    #         log=log,
    #         settings=settings
    #     )
    #     this.get()

    def test_soxs_disp_solution_function_exception(self):

        from soxspipe.recipes import soxs_disp_solution
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mdark/sofs/nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_disp_solution
            this = soxs_disp_solution(
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
