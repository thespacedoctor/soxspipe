
#!/usr/bin/env python
# encoding: utf-8
"""
*Some common setup methods used by more than one utility*

:Author:
    David Young

:Date Created:
    June 22, 2026
"""
import os

os.environ['TERM'] = 'vt100'

class base_util:
    """
    *The base utility class which all other utilities inherit*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary (default: False)
        - ``associatedFrame`` -- the associated frame the utility is working with (default: False)
        - ``dispersionMap`` -- if passed then `read_spectral_format` will be called to give info on the detector
          format (default: False)
        - ``twoDMapPath`` -- path to the 2D dispersion map. If passed, the map with be opened as a CCDData
          object (default: False)

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
            super(new_util, self).__init__(
                log,
                settings,
                associatedFrame=associatedFrame,
                dispersionMap=dispersionMap,
                twoDMapPath=twoDMapPath
            )

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

        from astropy.io import fits

        from soxspipe.commonutils import detector_lookup, keyword_lookup
        from soxspipe.commonutils.toolkit import get_skylines_dataframe

        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        if associatedFrame is not False:
            self.arm = associatedFrame.header[self.kw("SEQ_ARM")]
            self.dateObs = associatedFrame.header[self.kw("DATE_OBS")]

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(log=self.log, settings=self.settings).get(self.arm)

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        self.binx, self.biny = self._read_frame_binning(associatedFrame)

        # GET SKYLINES DATAFRAME
        self.skylinesDF = get_skylines_dataframe(
            self.log, self.settings, self.arm, minBrightnessVIS=1, minBrightnessNIR=0.5
        )

        # SET IMAGE ORIENTATION
        self.dispersionAxis = self.detectorParams["dispersion-axis"]
        self.axisA, self.axisB = self._image_orientation_axes(self.dispersionAxis)

        if dispersionMap:
            # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
            (
                self.orderNums,
                self.waveLengthMin,
                self.waveLengthMax,
                self.amins,
                self.amaxs,
            ) = self._read_spectral_format_limits(dispersionMap)

        if self.twoDMapPath:
            self.mapDF, self.interOrderMaskNDArray = self._read_two_d_map_dataframe(twoDMapPath, associatedFrame)
            self._mask_inter_order_pixels(associatedFrame, self.interOrderMaskNDArray)

        # OPEN AND UNPACK THE 2D IMAGE MAP
        if twoDMapPath:
            self.twoDMap = fits.open(twoDMapPath)

            dpBinx, dpBiny = self._read_map_binning()

            binxRatio = self.binx / dpBinx
            binyRatio = self.biny / dpBiny

            xdim = int(self.twoDMap[0].data.shape[1] / binxRatio)
            ydim = int(self.twoDMap[0].data.shape[0] / binyRatio)

            self._rebin_two_d_map(binxRatio, binyRatio)
            self.imageMap = self._build_image_map(associatedFrame, xdim, ydim)

        return

    def _read_frame_binning(
            self,
            associatedFrame):
        """*read the window binning from the associated frame's header, defaulting to 1x1*

        **Key Arguments:**
            - ``associatedFrame`` -- the associated frame the utility is working with, or False

        **Return:**
            - ``binx`` -- the binning along the x-axis
            - ``biny`` -- the binning along the y-axis
        """
        binx = 1
        biny = 1
        try:
            binx = int(associatedFrame.header[self.kw("WIN_BINX")])
            biny = int(associatedFrame.header[self.kw("WIN_BINY")])
        except (KeyError, TypeError, ValueError) as e:
            self.log.debug(f"__init__: `self.binx = int(associatedFrame.header[self.kw('WI...` failed, continuing: {e}")

        return binx, biny

    @staticmethod
    def _image_orientation_axes(
            dispersionAxis):
        """*name the dispersion and cross-dispersion axes for this detector*

        **Key Arguments:**
            - ``dispersionAxis`` -- the detector's dispersion axis, "x" or "y"

        **Return:**
            - ``axisA`` -- the axis running along the dispersion direction
            - ``axisB`` -- the axis running across the dispersion direction
        """
        if dispersionAxis == "x":
            return "x", "y"

        return "y", "x"

    def _read_spectral_format_limits(
            self,
            dispersionMap):
        """*read the spectral format table to determine the limits of the order traces*

        **Key Arguments:**
            - ``dispersionMap`` -- path to the dispersion map solution

        **Return:**
            - ``orderNums`` -- the order numbers
            - ``waveLengthMin`` -- the minimum wavelength of each order
            - ``waveLengthMax`` -- the maximum wavelength of each order
            - ``amins`` -- the minimum pixel position of each order along the dispersion axis
            - ``amaxs`` -- the maximum pixel position of each order along the dispersion axis
        """
        from soxspipe.commonutils.toolkit import read_spectral_format

        return read_spectral_format(
            log=self.log,
            settings=self.settings,
            arm=self.arm,
            dispersionMap=dispersionMap,
            extended=False,
            binx=self.binx,
            biny=self.biny,
        )

    def _read_two_d_map_dataframe(
            self,
            twoDMapPath,
            associatedFrame):
        """*unpack the 2D dispersion map image into a dataframe and an inter-order mask*

        **Key Arguments:**
            - ``twoDMapPath`` -- path to the 2D dispersion map
            - ``associatedFrame`` -- the associated frame the utility is working with

        **Return:**
            - ``mapDF`` -- the 2D dispersion map as a dataframe
            - ``interOrderMaskNDArray`` -- the mask of the pixels lying between the orders
        """
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe

        return twoD_disp_map_image_to_dataframe(
            log=self.log,
            slit_length=self.detectorParams["slit_length"],
            twoDMapPath=twoDMapPath,
            associatedFrame=associatedFrame,
            kw=self.kw,
            dispAxis=self.detectorParams["dispersion-axis"],
        )

    @staticmethod
    def _mask_inter_order_pixels(
            associatedFrame,
            interOrderMaskNDArray):
        """*set the frame's inter-order pixels to NaN, in place*

        **Key Arguments:**
            - ``associatedFrame`` -- the associated frame the utility is working with
            - ``interOrderMaskNDArray`` -- the mask of the pixels lying between the orders
        """
        import numpy as np

        associatedFrame.data[interOrderMaskNDArray == 1] = np.nan

        return

    def _read_map_binning(
            self):
        """*read the window binning from the 2D map's primary header, defaulting to 1x1*

        **Return:**
            - ``dpBinx`` -- the 2D map's binning along the x-axis
            - ``dpBiny`` -- the 2D map's binning along the y-axis
        """
        try:
            dpBinx = self.twoDMap[0].header[self.kw("WIN_BINX")]
            dpBiny = self.twoDMap[0].header[self.kw("WIN_BINY")]
        except KeyError as e:
            self.log.debug(f"__init__: `dpBinx = self.twoDMap[0].header[self.kw('WIN_B...` failed, continuing: {e}")
            dpBinx = 1
            dpBiny = 1

        return dpBinx, dpBiny

    def _rebin_two_d_map(
            self,
            binxRatio,
            binyRatio):
        """*block-reduce the 2D map's planes onto the associated frame's binning*

        **Key Arguments:**
            - ``binxRatio`` -- the frame's x-binning divided by the map's x-binning
            - ``binyRatio`` -- the frame's y-binning divided by the map's y-binning
        """
        import numpy as np

        if binxRatio > 1 or binyRatio > 1:
            from astropy.nddata import block_reduce

            self.twoDMap["WAVELENGTH"].data = block_reduce(
                self.twoDMap["WAVELENGTH"].data, (binyRatio, binxRatio), func=np.mean
            )
            self.twoDMap["SLIT"].data = block_reduce(self.twoDMap["SLIT"].data, (binyRatio, binxRatio), func=np.mean)
            self.twoDMap["ORDER"].data = block_reduce(self.twoDMap["ORDER"].data, (binyRatio, binxRatio), func=np.mean)

        return

    def _build_image_map(
            self,
            associatedFrame,
            xdim,
            ydim):
        """*associate each frame pixel with its wavelength, slit position and order*

        **Key Arguments:**
            - ``associatedFrame`` -- the associated frame the utility is working with
            - ``xdim`` -- the map's x-dimension, rebinned onto the frame's binning
            - ``ydim`` -- the map's y-dimension, rebinned onto the frame's binning

        **Return:**
            - ``imageMap`` -- the image map dataframe, with the inter-order rows removed
        """
        import numpy as np
        import pandas as pd

        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)

        imageMap = pd.DataFrame.from_dict(
            {
                "x": xarray,
                "y": yarray,
                "wavelength": self.twoDMap["WAVELENGTH"].data.flatten().astype(np.float32),
                "slit_position": self.twoDMap["SLIT"].data.flatten().astype(np.float32),
                "order": self.twoDMap["ORDER"].data.flatten().astype(np.float32),
                "flux": associatedFrame.data.flatten().astype(np.float32),
            }
        )
        imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)

        return imageMap
