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
from astropy.nddata import CCDData
from astropy import units as u
import numpy as np
import math
import numpy.ma as ma

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

    def test_unpack_order_table_function(self):
        orderTablePath = "/Users/Dave/soxspipe-unittests/intermediate/order_locations_NIR.csv"
        # UNPACK THE ORDER TABLE
        from soxspipe.commonutils.toolkit import unpack_order_table
        orderTable = unpack_order_table(
            log=log, orderTablePath=orderTablePath)

        from tabulate import tabulate
        print(tabulate(orderTable, headers='keys', tablefmt='psql'))

        # print(orders)
        # print(orderCentres)
        # print(orderLimits)

    def test_soxs_mflat_nir_normalise_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/nir_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        # this.produce_product()

        # inputFiles = "/path/to/a/directory"
        # inputFiles =
        # ['/path/to/one.fits','/path/to/two.fits','/path/to/three.fits']
        inputFiles = '~/xshooter-pipeline-data/unittest_data/xshooter-mflat/nir/calibrated/'
        from soxspipe.commonutils import set_of_files
        sof = set_of_files(
            log=log,
            settings=settings,
            inputFrames=inputFiles
        )
        sofFile, supplementarySof = sof.get()
        # LIST OF CCDDATA OBJECTS
        ccds = [c for c in sofFile.ccds(ccd_kwargs={
                                        "hdu_uncertainty": 'ERRS', "hdu_mask": 'QUAL', "hdu_flags": 'FLAGS', "key_uncertainty_type": 'UTYPE'})]

        combined_normalised_frame = "/Users/Dave/soxspipe-unittests/intermediate/first_iteration_NIR_master_flat.fits"
        from astropy.nddata import CCDData
        from astropy import units as u
        combined_normalised_frame = CCDData.read(combined_normalised_frame, hdu=0, unit=u.electron,
                                                 hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS', key_uncertainty_type='UTYPE')
        orderTablePath = "/Users/Dave/soxspipe-unittests/intermediate/order_locations_NIR.csv"
        # UNPACK THE ORDER TABLE

        from soxspipe.commonutils.toolkit import unpack_order_table
        orderTable = unpack_order_table(
            log=log, orderTablePath=orderTablePath)

        # GENERATE A MASK CONTAINING ORDER LOCATIONS
        mask = np.ones_like(combined_normalised_frame.data)
        for lower, upper in zip(orderEdgeLow, orderEdgeUp):
            for y, xlow, xhigh in zip(lower[1], lower[0], upper[0]):
                mask[y][math.ceil(xlow):math.floor(xhigh)] = 0

        # SUBTRACT NORMALISED FRAME FROM INDIVIDUAL CALIBRATED FRAMES
        medianNormalisedCCDs = []
        medianNormalisedCCDs[:] = [ccd.divide(np.ma.median(ma.array(ccd.subtract(
            combined_normalised_frame).data, mask=mask))) for ccd in ccds]

        exposure_normalised_flat = this.clip_and_stack(
            frames=medianNormalisedCCDs, recipe="soxs_mflat")

        this._write(exposure_normalised_flat,
                    "/tmp/mflat.fits", overwrite=True)

        # frame = tmpCcds[0]

        # # PLOT ONE OF THE MASKED FRAMES TO CHECK
        # from soxspipe.commonutils.toolkit import quicklook_image
        # for frame in [tmpCcds[0]]:
        #     maskedFrame = ma.array(frame.data, mask=mask)
        #     quicklook_image(log=log, CCDObject=maskedFrame, show=True)

        # print(np.ma.median(maskedFrame))

        import matplotlib.pyplot as plt
        frame = exposure_normalised_flat

        rotatedImg = np.flipud(np.rot90(frame.data, 1))
        std = np.std(frame.data)
        mean = np.mean(frame.data)
        vmax = mean + 3 * std
        vmin = mean - 3 * std
        plt.figure(figsize=(12, 5))

        plt.imshow(rotatedImg, vmin=vmin, vmax=vmax,
                   cmap='gray', alpha=1)
        plt.colorbar()
        plt.xlabel(
            "y-axis", fontsize=10)
        plt.ylabel(
            "x-axis", fontsize=10)

        for l, u in zip(orderEdgeLow, orderEdgeUp):

            plt.fill_between(l[1], l[0], u[0], alpha=0.4)

        plt.gca().invert_yaxis()

        plt.show()

        # EXPOSURE NORMAALISE FRAMES

        # this.produce_product()

    def test_soxs_mflat_nir_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/nir_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_nir_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/nir_short_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_uvb_dflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/uvb_dflats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_uvb_qflat_function(self):

        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/uvb_qflats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_vis_long_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/vis_long_flats.sof"
        from soxspipe.recipes import soxs_mflat
        this = soxs_mflat(
            log=log,
            settings=settings,
            inputFrames=sofPath
        )
        this.produce_product()

    def test_soxs_mflat_vis_short_function(self):
        sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof/vis_short_flats.sof"
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
            sofPath = "~/xshooter-pipeline-data/unittest_data/xshooter-mflat/sof//nir_mixed_exptime_darks.sof"
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
