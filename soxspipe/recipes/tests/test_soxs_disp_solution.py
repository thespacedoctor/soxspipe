from __future__ import print_function
import pytest
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


class test_soxs_disp_solution(unittest.TestCase):

    def test_soxs_real_disp_solution_nir_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/soxs/disp-solution/sof/NIR_DISP_SOLUTION.sof"

        # frame = "~/xshooter-pipeline-data/unittest_data/soxs/disp-solution/NIR/SOXS_NIR_SPH_ARC_OFF_10_002.fits"
        # from astropy.io import fits
        # with fits.open(frame, "readonly") as hdul:
        #     # PRINT SUMMARY OF FITS CONTENT
        #     print(hdul.info())
        #     # READ HEADER INTO MEMORY
        #     hdr = hdul[0].header + hdul[1].header
        #     data = hdul[1].data
        #     hdr["ESO DPR TYPE"] = "WAVE,LAMP"
        # # WRITE OUT A NEW FITS FILE
        # from astropy.io import fits
        # fits.writeto(frame, data, hdr,
        #              output_verify="fix+warn", overwrite=True, checksum=True)

        # frame = "~/xshooter-pipeline-data/unittest_data/soxs/disp-solution/NIR/SOXS_NIR_SPH_ARC_ON_10_004.fits"
        # from astropy.io import fits
        # with fits.open(frame, "readonly") as hdul:
        #     # PRINT SUMMARY OF FITS CONTENT
        #     print(hdul.info())
        #     # READ HEADER INTO MEMORY
        #     hdr = hdul[0].header + hdul[1].header
        #     data = hdul[1].data
        #     hdr["ESO DPR TYPE"] = "WAVE,LAMP"
        # # WRITE OUT A NEW FITS FILE
        # from astropy.io import fits
        # fits.writeto(frame, data, hdr,
        #              output_verify="fix+warn", overwrite=True, checksum=True)

        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log2,
            settings=settings2,
            inputFrames=sofPath,
            overwrite=True
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    def test_xsh_disp_solution_nir_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-disp-solution/sof/20170818_NIR_DISP_SOLUTION.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    def test_xsh_disp_solution_uvb_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-disp-solution/sof/20170818_UVB_DISP_SOLUTION_1x1_fast.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    def test_xsh_disp_solution_vis_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-disp-solution/sof/20170818_VIS_DISP_SOLUTION_1x1_fast.sof"
        from soxspipe.recipes import soxs_disp_solution
        disp_map_path = soxs_disp_solution(
            log=log,
            settings=settings,
            inputFrames=sofPath,
            overwrite=True
        ).produce_product()
        print(f"Here is the final product `{disp_map_path}`")

    @pytest.mark.full
    def test_soxs_disp_solution_function_exception(self):

        from soxspipe.recipes import soxs_disp_solution
        try:
            sofPath = "~/xshooter-pipeline-data/unittest_data/xsh/xshooter-mdark/sofs/nir_mixed_exptime_darks.sof"
            from soxspipe.recipes import soxs_disp_solution
            this = soxs_disp_solution(
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
