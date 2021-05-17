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
from astropy.io import fits
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
# xt-utkit-refresh-database

class test_filenamer(unittest.TestCase):

    def test_filenamer_function(self):
        from soxspipe.commonutils import filenamer
        utDir = settings["test-data-root"]
        subDirs = [
            "create_dispersion_map",
            "detect_continuum",
            "xshooter-mbias/uvb/1x1/fast_read",
            "xshooter-mbias/uvb/1x1/slow_read",
            "xshooter-mbias/vis",
            "xshooter-mdark/nir/240s",
            "xshooter-detect-order-edges",
            "xshooter-mflat/uvb/qflat",
            "xshooter-mflat/nir/calibrated",
            "xshooter-order-centres/uvb/",
            "xshooter-disp-solution/nir/"
        ]
        filestoname = []
        for dirr in subDirs:
            # GENERATE A LIST OF FILE PATHS
            pathToDirectory = f"{utDir}/{dirr}"
            for d in os.listdir(pathToDirectory):
                filepath = os.path.join(pathToDirectory, d)
                if os.path.isfile(filepath) and os.path.splitext(filepath)[1] == ".fits":
                    filestoname.append(filepath)

        for filepath in filestoname:
            print(f"\nORIG: {filepath}")
            # with fits.open(filepath, "readonly") as hdul:
            frame = CCDData.read(filepath, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                                 hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
            if "SOXSPIPE PRE" in frame.header:
                frame.header["SXSPRE"] = frame.header["SOXSPIPE PRE"]

            filename = filenamer(
                log=log,
                frame=frame,
                settings=settings
            )
            print(filename)

    def test_filenamer_function_exception(self):

        from soxspipe.commonutils import filenamer
        try:
            this = filenamer(
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
