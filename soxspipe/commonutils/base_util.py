#!/usr/bin/env python
# encoding: utf-8
"""
*Some common setup methods used by more than one utility*

:Author:
    David Young

:Date Created:
    June 22, 2026
"""
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools

class base_util(object):
    """
    *The base utility class which all other utilities inherit*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary (default: False)
        - ``associatedFrame`` -- the associated frame the utility is working with (default: False)
        - ``dispersionMap`` -- if passed then `read_spectral_format` will be called to give info on the detector format (default: False)
        - ``twoDMapPath`` -- path to the 2D dispersion map. If passed, the map with be opened as a CCDData object (default: False)

    **Usage:**

    To initiate a base_util object, use the following:

    ```python
    from .base_util import base_util
    class new_util(base_util):
        def __init__(
            self,
            log,
            settings,
            associatedFrame=False,
            dispersionMap=False,
            twoDMapPath=False
            # other arguments needed for the new_util class
        ):
        super(new_util, self).__init__(log, settings, associatedFrame=associatedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath)

        ...
    ```
    """
    def __init__(
            self,
            log,
            settings,
            associatedFrame=False,
            dispersionMap=False,
            twoDMapPath=False
    ):
        self.log = log
        log.debug("instansiating a new 'base_util' object")
        self.settings = settings
        self.dispersionMap = dispersionMap
        self.twoDMapPath = twoDMapPath

        import numpy as np
        from astropy.io import fits
        import numpy as np
        import pandas as pd
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe, read_spectral_format
        from soxspipe.commonutils import detector_lookup, keyword_lookup
        from soxspipe.commonutils.toolkit import get_skylines_dataframe
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe

        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        if associatedFrame is not False:
            self.arm = associatedFrame.header[self.kw("SEQ_ARM")]
            self.dateObs = associatedFrame.header[self.kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        self.binx = 1
        self.biny = 1
        try:
            self.binx = int(associatedFrame.header[self.kw("WIN_BINX")])
            self.biny = int(associatedFrame.header[self.kw("WIN_BINY")])
        except:
            pass

        # GET SKYLINES DATAFRAME
        self.skylinesDF = get_skylines_dataframe(
            self.log, self.settings, self.arm, minBrightnessVIS=1, minBrightnessNIR=0.5
        )

        # SET IMAGE ORIENTATION
        self.dispersionAxis = self.detectorParams["dispersion-axis"]
        if self.dispersionAxis == "x":
            self.axisA = "x"
            self.axisB = "y"
        else:
            self.axisA = "y"
            self.axisB = "x"

        if dispersionMap:
            # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
            self.orderNums, self.waveLengthMin, self.waveLengthMax, self.amins, self.amaxs = read_spectral_format(
                log=self.log,
                settings=self.settings,
                arm=self.arm,
                dispersionMap=dispersionMap,
                extended=False,
                binx=self.binx,
                biny=self.biny,
            )

        if self.twoDMapPath:
            self.mapDF, self.interOrderMaskNDArray = twoD_disp_map_image_to_dataframe(
                log=self.log,
                slit_length=self.detectorParams["slit_length"],
                twoDMapPath=twoDMapPath,
                associatedFrame=associatedFrame,
                kw=self.kw,
                dispAxis=self.detectorParams["dispersion-axis"],
            )
            associatedFrame.data[self.interOrderMaskNDArray == 1] = np.nan

        # OPEN AND UNPACK THE 2D IMAGE MAP
        if twoDMapPath:
            self.twoDMap = fits.open(twoDMapPath)

            try:
                dpBinx = self.twoDMap[0].header[self.kw("WIN_BINX")]
                dpBiny = self.twoDMap[0].header[self.kw("WIN_BINY")]
            except:
                dpBinx = 1
                dpBiny = 1

            binxRatio = self.binx / dpBinx
            binyRatio = self.biny / dpBiny

            xdim = int(self.twoDMap[0].data.shape[1] / binxRatio)
            ydim = int(self.twoDMap[0].data.shape[0] / binyRatio)
            xarray = np.tile(np.arange(0, xdim), ydim)
            yarray = np.repeat(np.arange(0, ydim), xdim)

            if binxRatio > 1 or binyRatio > 1:
                from astropy.nddata import block_reduce

                self.twoDMap["WAVELENGTH"].data = block_reduce(
                    self.twoDMap["WAVELENGTH"].data, (binyRatio, binxRatio), func=np.mean
                )
                self.twoDMap["SLIT"].data = block_reduce(self.twoDMap["SLIT"].data, (binyRatio, binxRatio), func=np.mean)
                self.twoDMap["ORDER"].data = block_reduce(self.twoDMap["ORDER"].data, (binyRatio, binxRatio), func=np.mean)

            self.imageMap = pd.DataFrame.from_dict(
                {
                    "x": xarray,
                    "y": yarray,
                    "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(np.float32),
                    "slit_position": self.twoDMap["SLIT"].data.flatten().astype(np.float32),
                    "order": self.twoDMap["ORDER"].data.flatten().astype(np.float32),
                    "flux": associatedFrame.data.flatten().astype(np.float32),
                }
            )
            self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)


        return None

    