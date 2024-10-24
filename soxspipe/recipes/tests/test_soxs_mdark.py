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


class test_soxs_mdark(unittest.TestCase):

    import pytest

    def test_xsh_mdark_nir_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/nir_6s_darks.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    @pytest.mark.full
    def test_xsh_mdark_nir_function2(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/nir_240s_darks.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    @pytest.mark.full
    def test_xsh_mdark_nir_function3(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/nir_125s_darks.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_xsh_mdark_uvb_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/uvb_1x1_dark_3600s.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    @pytest.mark.full
    def test_xsh_mdark_uvb_function2(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/uvb_2x2_dark_3600s.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    def test_xsh_mdark_vis_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/vis_1x1_dark_3600s.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    @pytest.mark.full
    def test_xsh_mdark_vis_function2(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/vis_2x2_dark_3600s.sof"
        from soxspipe.recipes import soxs_mdark
        this = soxs_mdark(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        )
        productPath = this.produce_product()
        print(f"Here is the final product `{productPath}`")

    # @pytest.mark.full
    # def test_soxs_mdark_e2e_nir_folder_function(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/dark/nir/180"
    #     # sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/mdark/sofs/e2e_nir_180s_no_crh_darks.sof"
    #     from soxspipe.recipes import soxs_mdark
    #     this = soxs_mdark(
    #         log=log2,
    #         settings=settings2,
    #         inputFrames=sofPath,
    #         verbose=True,
    #         overwrite=True
    #     )
    #     productPath = this.produce_product()
    #     print(f"Here is the final product `{productPath}`")
    # def test_soxs_mdark_e2e_nir_function(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/dark/sofs/soxsim_nir_180s_darks.sof"
    #     # sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/mdark/sofs/e2e_nir_180s_no_crh_darks.sof"
    #     from soxspipe.recipes import soxs_mdark
    #     this = soxs_mdark(
    #         log=log2,
    #         settings=settings2,
    #         inputFrames=sofPath,
    #         verbose=True,
    #         overwrite=True
    #     )
    #     productPath = this.produce_product()
    #     print(f"Here is the final product `{productPath}`")
    # @pytest.mark.full
    # def test_soxs_mdark_e2e_nir_function2(self):
    #     sofPath = "~/xshooter-pipeline-data/unittest_data/soxs-sim/dark_CR/sofs/soxsim_nir_300s_wCRH_darks.sof"
    #     from soxspipe.recipes import soxs_mdark
    #     this = soxs_mdark(
    #         log=log2,
    #         settings=settings2,
    #         inputFrames=sofPath,
    #         verbose=True,
    #         overwrite=True
    #     )
    #     productPath = this.produce_product()
    #     print(f"Here is the final product `{productPath}`")
    # def test_soxs_mdark_function(self):
    #     # utKit.refresh_database() # reset database to database found in
    #     # soxspipe/test/input
    #     from soxspipe.recipes import soxs_mdark
    #     this = soxs_mdark(
    #         log=log,
    #         settings=settings
    #     )
    #     this.get()
    @pytest.mark.full
    def test_soxs_mdark_function_exception(self):

        from soxspipe.recipes import soxs_mdark
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_mdark
            this = soxs_mdark(
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
