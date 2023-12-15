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
settingsFile = packageDirectory + "/test_settings_xsh.yaml"
# settingsFile = home + \
#     "/git_repos/_misc_/settings/soxspipe/test_settings_xsh.yaml"

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
    import pytest

    @pytest.mark.full
    def test_xsh_create_dispersion_map_single_nir_function(self):
        import pandas as pd
        frame = "~/xshooter-pipeline-data/unittest_data/xsh/create_dispersion_map/single_pinhole_NIR_calibrated.fits"
        from os.path import expanduser
        home = expanduser("~")
        frame = frame.replace("~", home)
        frame = CCDData.read(frame, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                             hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": [],
            "label": []
        })

        from soxspipe.commonutils import create_dispersion_map
        mapPath, mapImagePath, res_plots, qcTable, productsTable = create_dispersion_map(
            log=log,
            settings=settings,
            pinholeFrame=frame,
            qcTable=qc,
            productsTable=products,
            create2DMap=False
        ).get()
        print(mapPath)

    @pytest.mark.full
    def test_xsh_create_dispersion_map_multi_nir_function(self):
        import pandas as pd
        frame = "~/xshooter-pipeline-data/unittest_data/xsh/create_dispersion_map/20170818T173315_NIR_ARC_MULTIPIN.fits"
        from os.path import expanduser
        home = expanduser("~")
        frame = frame.replace("~", home)
        frame = CCDData.read(frame, hdu=0, unit=u.electron, hdu_uncertainty='ERRS',
                             hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')

        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": [],
            "label": []
        })

        from soxspipe.commonutils import create_dispersion_map
        mapPath, mapImagePath, res_plots, qcTable, productsTable = create_dispersion_map(
            log=log,
            settings=settings,
            pinholeFrame=frame,
            firstGuessMap="~/xshooter-pipeline-data/unittest_data/xsh/create_dispersion_map/20170818T172310_NIR_DISP_MAP.fits",
            qcTable=qc,
            productsTable=products,
            create2DMap=False
        ).get()
        print(mapPath)

    @pytest.mark.full
    def test_soxs_create_dispersion_map_function_exception(self):

        from soxspipe.commonutils import create_dispersion_map
        try:
            this, this2, this3 = create_dispersion_map(
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
