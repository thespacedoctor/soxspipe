#!/usr/bin/env python
# encoding: utf-8
"""
*fit and subtract background flux from scattered light from frame*

Author
: David Young

Date Created
: June  3, 2021
"""
from soxspipe.commonutils import keyword_lookup
from os.path import expanduser
from soxspipe.commonutils.toolkit import quicklook_image
from soxspipe.commonutils.toolkit import unpack_order_table
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class subtract_background(object):
    """
    *fit and subtract background flux from scattered light from frame*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``frame`` -- the frame to subtract background light from
    - ``recipeName`` -- name of the parent recipe
    - ``sofName`` -- the sof file name given to the parent recipe
    - ``orderTable`` -- the order geometry table
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products
    - ``lamp`` -- needed for UVB flats

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/master/initialisation.html).

    To fit and subtract the background from an image:

    ```python
    from soxspipe.commonutils import subtract_background
    background = subtract_background(
        log=log,
        frame=myCCDDataObject,
        orderTable="/path/to/orderTable",
        settings=settings
    )
    backgroundFrame, backgroundSubtractedFrame = background.subtract()
    ```

    """

    def __init__(
            self,
            log,
            frame,
            orderTable,
            sofName=False,
            recipeName=False,
            settings=False,
            qcTable=False,
            productsTable=False,
            lamp=""
    ):
        self.log = log
        log.debug("instantiating a new 'subtract_background' object")
        self.settings = settings
        self.sofName = sofName
        self.recipeName = recipeName
        self.frame = frame
        self.orderTable = orderTable
        self.qc = qcTable
        self.products = productsTable
        self.lamp = lamp

        from soxspipe.commonutils import detector_lookup

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(
            log=self.log,
            settings=self.settings
        ).get
        kw = self.kw
        self.arm = frame.header[kw("SEQ_ARM")]
        self.dateObs = frame.header[kw("DATE_OBS")]

        self.inst = frame.header[kw("INSTRUME")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(
            log=log,
            settings=settings
        ).get(self.arm)

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext='data', stdWindow=3, surfacePlot=True, title="Initial input frame needing scattered light subtraction")

        return None

    def subtract(self):
        """
        *fit and subtract background light from frame*

        **Return:**

        - ``backgroundSubtractedFrame`` -- a CCDData object of the original input frame with fitted background light subtracted
        """
        self.log.debug('starting the ``subtract`` method')

        import numpy as np
        import pandas as pd

        kw = self.kw
        imageType = self.frame.header[kw("DPR_TYPE")].replace(",", "-")
        imageTech = self.frame.header[kw("DPR_TECH")].replace(",", "-")
        imageCat = self.frame.header[kw("DPR_CATG")].replace(",", "-")

        self.log.print(f"\n# FITTING AND SUBTRACTING SCATTERED LIGHT BACKGROUND FROM {self.arm} {imageCat} {imageTech} {imageType} FRAME")

        binx = 1
        biny = 1
        try:
            binx = self.frame.header[self.kw("WIN_BINX")]
            biny = self.frame.header[self.kw("WIN_BINY")]
        except:
            pass

        # UNPACK THE ORDER TABLE
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=self.orderTable, binx=binx, biny=biny, extend=4)

        originalMask = np.copy(self.frame.mask)

        quicklook_image(
            log=self.log, CCDObject=self.frame.data, show=False, ext=None, surfacePlot=True, title="Initial frame")

        # MASK THE INNER ORDER AREA (AND BAD PIXELS)
        self.mask_order_locations(orderPixelTable)

        quicklook_image(
            log=self.log, CCDObject=self.frame, show=False, ext=None, surfacePlot=True, title="Initial input frame with order pixels masked")

        backgroundFrame = self.create_background_image(
            rowFitOrder=self.settings['background-subtraction']['bspline-deg'], gaussianSigma=self.settings['background-subtraction']['gaussian-blur-sigma'])

        # REPLACE MASK
        self.frame.mask = originalMask

        backgroundSubtractedFrame = self.frame.subtract(backgroundFrame)
        backgroundSubtractedFrame.header = self.frame.header

        # GET FILENAME FOR THE RESIDUAL PLOT
        saveToPath = False
        if self.sofName:
            backgroundQCImage = self.sofName + f"_BKGROUND{self.lamp}.pdf"
            home = expanduser("~")
            self.qcDir = self.settings["workspace-root-dir"].replace("~", home) + f"/qc/{self.recipeName}/"
            self.qcDir = self.qcDir.replace("//", "/")
            # RECURSIVELY CREATE MISSING DIRECTORIES
            if not os.path.exists(self.qcDir):
                os.makedirs(self.qcDir)
            saveToPath = self.qcDir + "/" + backgroundQCImage

        quicklook_image(
            log=self.log, CCDObject=backgroundFrame, show=False, ext='data', stdWindow=9, surfacePlot=False, title="Background Scattered Light", saveToPath=saveToPath)

        if saveToPath and not isinstance(self.products, bool):
            from datetime import datetime
            utcnow = datetime.utcnow()
            utcnow = utcnow.strftime("%Y-%m-%dT%H:%M:%S")
            self.products = pd.concat([self.products, pd.Series({
                "soxspipe_recipe": self.recipeName,
                "product_label": f"BKGROUND{self.lamp}",
                "file_name": backgroundQCImage,
                "file_type": "PDF",
                "obs_date_utc": self.dateObs,
                "reduction_date_utc": utcnow,
                "product_desc": f"Fitted intra-order image background{self.lamp.replace('_',' ')}",
                "file_path": saveToPath,
                "label": "QC"
            }).to_frame().T], ignore_index=True)

        self.log.debug('completed the ``subtract`` method')
        return backgroundFrame, backgroundSubtractedFrame, self.products

    def mask_order_locations(
            self,
            orderPixelTable):
        """*mask the order locations and return the masked frame*

        **Key Arguments:**

        - ``orderPixelTable`` -- the order location in a pandas datafrmae.
        """
        self.log.debug('starting the ``mask_order_locations`` method')

        import numpy as np

        # MASK DATA INSIDE OF ORDERS (EXPAND THE INNER-ORDER AREA IF NEEDED)
        uniqueOrders = np.sort(orderPixelTable['order'].unique())
        if self.axisA == "x":
            axisALen = self.frame.data.shape[1]
            axisBLen = self.frame.data.shape[0]
        else:
            axisALen = self.frame.data.shape[0]
            axisBLen = self.frame.data.shape[1]

        oTop = orderPixelTable['order'].min()
        oBot = orderPixelTable['order'].max()

        for o in uniqueOrders:

            axisBcoord = orderPixelTable.loc[
                (orderPixelTable["order"] == o)][f"{self.axisB}coord"]
            axisAcoord_edgeup = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                f"{self.axisA}coord_edgeup"]
            axisAcoord_edgelow = orderPixelTable.loc[(orderPixelTable["order"] == o)][
                f"{self.axisA}coord_edgelow"]

            if o != oBot:
                next_axisAcoord_edgeup = orderPixelTable.loc[(orderPixelTable["order"] == o + 1)][
                    f"{self.axisA}coord_edgeup"]
                bottomGap = axisAcoord_edgelow.values - next_axisAcoord_edgeup.values
                expandBottom = np.median(bottomGap) / 2 - 3
                if expandBottom < 2:
                    expandBottom = 2
            else:
                expandBottom = expandTop

            if o != oTop:
                previous_axisAcoord_edgelow = orderPixelTable.loc[(orderPixelTable["order"] == o - 1)][
                    f"{self.axisA}coord_edgelow"]
                topGap = previous_axisAcoord_edgelow.values - axisAcoord_edgeup.values
                expandTop = np.median(topGap) / 2 - 3
                if expandTop < 2:
                    expandTop = 2
            else:
                expandTop = expandBottom

            axisAcoord_edgeup += expandTop
            axisAcoord_edgelow -= expandBottom

            axisAcoord_edgelow, axisAcoord_edgeup, axisBcoord = zip(*[(x1, x2, b) for x1, x2, b in zip(axisAcoord_edgelow, axisAcoord_edgeup, axisBcoord) if x1 < axisALen and x2 > 0 and x2 < axisALen and b > 0 and b < axisBLen])
            for b, u, l in zip(axisBcoord, np.ceil(axisAcoord_edgeup).astype(int), np.floor(axisAcoord_edgelow).astype(int)):
                if l < 0:
                    l = 0
                if self.axisA == "x":
                    self.frame.mask[b, l:u] = 1
                else:
                    self.frame.mask[l:u, b] = 1

            if o == oTop:
                # MASK TOP OF FRAME
                mask_bottom = np.array(axisAcoord_edgeup) + 7
                for b, m, in zip(axisBcoord, np.ceil(mask_bottom).astype(int)):
                    if self.axisA == "x":
                        self.frame.mask[b, m:] = 1
                    else:
                        self.frame.mask[m:, b] = 1

            if o == oBot:
                # MASK BOTTOM OF FRAME
                mask_top = np.array(axisAcoord_edgelow) - 7
                for b, m, in zip(axisBcoord, np.ceil(mask_top).astype(int)):
                    if m < 0:
                        m = 0
                    if self.axisA == "x":
                        self.frame.mask[b, :m] = 1
                    else:
                        self.frame.mask[:m, b] = 1

        self.log.debug('completed the ``mask_order_locations`` method')
        return None

    def create_background_image(
            self,
            rowFitOrder,
            gaussianSigma):
        """*model the background image from intra-order flux detected*

        **Key Arguments:**

        - ``rowFitOrder`` -- order of the polynomial fit to flux in each row
        - ``gaussianSigma`` -- the sigma of the gaussian used to blur the final image
        """
        self.log.debug('starting the ``create_background_image`` method')

        from astropy.stats import sigma_clip, mad_std
        import numpy as np
        import pandas as pd
        from astropy.nddata import CCDData
        from scipy.signal import medfilt2d
        from scipy.interpolate import splrep, splev
        import scipy
        import numpy.ma as ma
        import math
        import random
        from soxspipe.commonutils.filenamer import filenamer
        from os.path import expanduser

        maskedImage = np.ma.array(self.frame.data, mask=self.frame.mask)
        # SIGMA-CLIP THE DATA
        clippedDataMask = sigma_clip(
            maskedImage, sigma_lower=10, sigma_upper=50, maxiters=2, cenfunc='median', stdfunc='mad_std')

        # COMBINE MASK WITH THE BAD PIXEL MASK
        mask = (clippedDataMask.mask == 1) | (self.frame.mask == 1)
        maskedImage = np.ma.array(self.frame.data, mask=mask)

        maskedAsNanImage = self.frame.data.copy()
        maskedAsNanImage[mask] = np.nan

        maskedAsNanImage[:20, :] = 0
        maskedAsNanImage[-20:, :] = 0
        maskedAsNanImage[:, -20:] = 0
        maskedAsNanImage[:, :20] = 0

        import skimage
        from skimage.morphology import disk
        maskedAsNanImage = skimage.filters.median(maskedAsNanImage, footprint=disk(4))

        # REMOVE -VE VALUES
        maskedAsNanImage = np.where(maskedAsNanImage < 0, 0, maskedAsNanImage)

        # REGENERATE A MASK OF NANS
        mask = np.isnan(maskedAsNanImage)
        maskedImage = np.ma.array(maskedAsNanImage, mask=mask)

        quicklook_image(
            log=self.log, CCDObject=maskedImage, show=False, ext=True, stdWindow=3, surfacePlot=True, title="Sigma clipped, blurred masked image")

        # PLACEHOLDER ARRAY FOR BACKGROUND IMAGE
        backgroundMap = np.zeros_like(self.frame)

        for idx, row in enumerate(maskedImage):

            # SET X TO A MASKED RANGE ... BLANK DATA BUT WITH MASK FROM IMAGE
            xunmasked = ma.masked_array(np.linspace(
                0, len(row), len(row), dtype=int), mask=row.mask)

            # fit = ma.polyfit(xunmasked, row, deg=rowFitOrder)
            xfit = xunmasked.data

            # yfit = np.polyval(fit, xfit)

            xmasked = xunmasked[~xunmasked.mask]
            xmin = xmasked.min()
            xmax = xmasked.max()
            rowmasked = row[~row.mask]

            window = 7
            hw = math.floor(window / 2)
            # rowmaskedSmoothed = pd.Series(rowmasked).rolling(window=window, center=True).quantile(.1)
            try:
                rowmaskedSmoothed = pd.Series(rowmasked).rolling(window=window, center=True).median()
            except:
                rowmasked = rowmasked.astype(float)
                # rowmasked = rowmasked.byteswap().newbyteorder() ## REMOVE IF ABOVE .astype(float) WORKS
                rowmaskedSmoothed = pd.Series(rowmasked).rolling(window=window, center=True).median()
            rowmaskedSmoothed[:hw] = rowmaskedSmoothed.iloc[hw + 1]
            rowmaskedSmoothed[-hw:] = rowmaskedSmoothed.iloc[-hw - 1]
            rowmasked[:hw] = rowmasked[hw + 1]
            rowmasked[-hw:] = rowmasked[-hw - 1]

            rowmaskedSmoothed = np.where(rowmaskedSmoothed < 0, 0, rowmaskedSmoothed)

            seedKnots = xmasked[1:-1:window * 2]
            tck, fp, ier, msg = splrep(xmasked, rowmaskedSmoothed, t=seedKnots, k=rowFitOrder, full_output=True)
            t, c, k = tck

            if ier == 10:
                self.log.info(f"\t\tpoor fit on columns {idx}.\n")

            yfit = splev(xfit, tck, ext=3)

            # ADD FITTED ROW TO BACKGROUND IMAGE
            backgroundMap[idx, :] = yfit

            if False and random.randint(1, 401) == 42:
                import matplotlib.pyplot as plt
                fig, (ax1) = plt.subplots(1, 1, figsize=(30, 15))
                plt.scatter(xmasked, rowmasked)
                plt.scatter(xmasked, rowmaskedSmoothed)
                plt.plot(xfit, yfit, 'b')
                plt.ylabel("flux")
                plt.xlabel("pixels along row")
                plt.show()

        quicklook_image(
            log=self.log, CCDObject=backgroundMap, show=False, ext=None, surfacePlot=True, title="Scattered light background image")

        backgroundMap = scipy.ndimage.filters.gaussian_filter(backgroundMap, gaussianSigma)

        # SET -VE T0 ZERO
        backgroundMap = np.where(backgroundMap < 0, 0, backgroundMap)

        quicklook_image(
            log=self.log, CCDObject=backgroundMap, show=False, ext=None, surfacePlot=True, title="Scattered light background image with median filtering")

        backgroundMap = CCDData(backgroundMap, unit=self.frame.unit)

        self.log.debug('completed the ``create_background_image`` method')
        return backgroundMap
