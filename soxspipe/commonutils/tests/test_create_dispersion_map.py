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
from astropy.nddata import CCDData
from astropy import units as u
home = expanduser("~")

packageDirectory = utKit("").get_project_root()
settingsFile = packageDirectory + "/test_settings.yaml"
# settingsFile = home + \
#     "/git_repos/_misc_/settings/soxspipe/test_settings.yaml"

su = tools(
    arguments={"settingsFile": settingsFile},
    docString=__doc__,
    logLevel="DEBUG",
    options_first=False,
    projectName=None,
    defaultSettingsFile=False
)
arguments, settings, log, dbConn = su.setup()

# SETUP PATHS TO COMMON DIRECTORIES FOR TEST DATA
moduleDirectory = os.path.dirname(__file__)
pathToInputDir = moduleDirectory + "/input/"
pathToOutputDir = moduleDirectory + "/output/"

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


class test_create_dispersion_map(unittest.TestCase):

    def test_create_dispersion_map_single_nir_function(self):
        frame = "~/xshooter-pipeline-data/unittest_data/create_dispersion_map/single_pinhole_NIR_calibrated.fits"
        from os.path import expanduser
        home = expanduser("~")
        frame = frame.replace("~", home)
        frame = CCDData.read(frame, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                             hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        from soxspipe.commonutils import create_dispersion_map
        mapPath = create_dispersion_map(
            log=log,
            settings=settings,
            pinholeFrame=frame
        ).get()
        print(mapPath)

    def test_create_dispersion_map_multi_nir_function(self):
        frame = "~/xshooter-pipeline-data/unittest_data/create_dispersion_map/20170818T173315_NIR_ARC_MULTIPIN.fits"
        from os.path import expanduser
        home = expanduser("~")
        frame = frame.replace("~", home)
        frame = CCDData.read(frame, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                             hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        from soxspipe.commonutils import create_dispersion_map
        mapPath = create_dispersion_map(
            log=log,
            settings=settings,
            pinholeFrame=frame,
            firstGuessMap="~/xshooter-pipeline-data/unittest_data/create_dispersion_map/20170820T153602_NIR_DISP_MAP.csv"
        ).get()
        print(mapPath)

    def test_create_dispersion_map_function_exception(self):

        from soxspipe.commonutils import create_dispersion_map
        try:
            this = create_dispersion_map(
                log=log,
                settings=settings,
                fakeKey="break the code"
            )
            this.get()
            assert False
        except Exception as e:
            assert True
            print(str(e))

        # x-print-testpage-for-pessto-marshall-web-object

    # x-class-to-test-named-worker-function
