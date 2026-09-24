#!/usr/bin/env python
# encoding: utf-8
"""
*Using a 2D dispersion image map, transform a SOXS data frame from xy detector pixel space to wavelength-slit position space.*

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
from .base_util import base_util



class image_transformer(base_util):
    """*zoom and rebin 2D frame arrays order-by-order into cross-dispersion slices ready for spectral extraction and image rectification*

    For each order:

    1. Zoom every input array along the spatial (slit) axis by ``zoomFactor`` using
       nearest-neighbour pixel replication (slice into sub-pixels).
    2. Indexes into the zoomed arrays using sub-pixel slit bounds derived from the
       order-pixel table, producing a 2D block of shape
       ``(nSlices, slitHalfLength * 2 * zoomFactor)``.
    3. Rebins that block back to ``(nSlices, slitHalfLength * 2)`` by block-averaging,
       giving finer effective slit sampling.
    4. Optionally sigma-clips the rebinned raw flux to catch additional outlier pixels before extraction.

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- settings dict
    - ``orderPixelTable`` -- dataframe containing per-pixel order-trace metadata (object continuum fit results)
    - ``twoDMapPath`` -- path to the 2D map FITS file (pixel wavelength and slit position values needed for rectification)
    - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
    - ``associatedFrame`` -- an example 2D frame to be rectified. This frame is used to determine detector binning, arm etc.
    - ``slitHalfLength`` -- half-length of the slit in pixels (sets extraction aperture)
    - ``zoomFactor`` -- integer zoom factor applied along the slit axis. Default *11*

    **Return:**

    - ``orderSlices`` -- list of per-order dataframes, each row being one cross-dispersion slice
    - ``wlMinMax`` -- list of ``(wlmin, wlmax)`` tuples in the same order as ``orderSlices``

    **Usage:**

    ```python
    from soxspipe.commonutils.image_transformer import image_transformer

    transformer = image_transformer(
        log=log,
        settings=settings,
        orderPixelTable=orderPixelTable,
        twoDMapPath=twoDMapPath,
        dispersionMap=dispersionMap,
        associatedFrame=skySubtractedFrame,
        slitHalfLength=slitHalfLength,
        zoomFactor=11
    )
    transformer.cache_image(
        imageName="fluxRaw",
        ndarray=skySubtractedFrame,
        associatedMask=badPixelMask
    )
    slices, wlMinMax = transformer.get_order_slices()
    ```
    """
    def __init__(
            self,
            log,
            settings,
            orderPixelTable,
            twoDMapPath,
            dispersionMap,
            associatedFrame,
            slitHalfLength,
            zoomFactor=11
    ):
        super(image_transformer, self).__init__(log, settings, associatedFrame=associatedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath)
        self.orderPixelTable = orderPixelTable
        self.twoDMapPath = twoDMapPath
        self.slitHalfLength = slitHalfLength
        self.zoomFactor = zoomFactor

        import numpy as np
        from soxspipe.commonutils.toolkit import twoD_disp_map_image_to_dataframe

        # ORDERS PRESENT IN THE TRACE TABLE — USED TO SKIP ORDERS WITH NO DETECTED TRACE
        self.uniqueOrders = self.orderPixelTable["order"].unique()
        self.sliceOrders = []

        # AFTER FANCY-INDEXING, THE EXTRACTED BLOCK IS ALWAYS (nSlices, slitHalfLength*2*zoomFactor)
        # SO THE REBIN SHAPE IS ALWAYS (1, zoomFactor) REGARDLESS OF DISPERSION AXIS
        self.rebinShape = (1, self.zoomFactor)

        # SETUP THE ORDER TABLE DATAFRAMES FOR EACH ORDER, INCLUDING SUB-PIXEL INDEXES INTO THE ZOOMED ARRAYS
        self._setup_order_table_dataframes()

        self._cache_image_names = set()

        # CACHE THE WAVELENGTH AND SLIT POSITION MAPS FOR LATER RECTIFICATION
        self.cache_image(
            imageName="wavelength",
            ndarray=self.twoDMap["WAVELENGTH"].data,
        )
        self.cache_image(
            imageName="slit",
            ndarray=self.twoDMap["SLIT"].data,
        )
    
        return None

    def cache_image(
            self,
            imageName,
            ndarray,
            associatedMask=None):
        """
        *Place an image in the transformer's cache. These images can be acted on later.*

        **Key Arguments:**

            - ``imageName`` -- the unique name to give to the image
            - ``ndarray`` -- the 2D frame to rectify and cache
            - ``associatedMask`` -- 2D bad-pixel mask. If passed with the raw flux image, the rebinned raw flux image will be sigma-clipped and new outlying sub-pixels will be added to the rebinned bad-pixel mask (also cached). Default *None*

        **Return:**
        
            - None
        """
        self.log.debug('starting the ``cache_image`` method')

        import numpy as np

        # ZOOM NEW NDARRAY
        zoomedArray = self._zoom_array(ndarray, self.zoomFactor, self.dispersionAxis)

        # ITERATE OVER EACH ORDER, RECTIFYING THE IMAGE INTO WAVELEGTH-SLIT SPACE 
        for order, orderTable, SPIndex in zip(self.sliceOrders, self.orderSlices, self.subPixelIndexes):

            # ZOOM AND REBIN EACH THE NDARRY AND STORE AS A COLUMN IN THE ORDER TABLE
            rebinned = self._rebin_2d(zoomedArray[SPIndex], *self.rebinShape)
            orderTable[imageName] = list(rebinned)
            self._cache_image_names.add(imageName)

            # REBIN THE BAD-PIXEL MASK AND SIGMA-CLIP AGAINST THE RAW FLUX TO FLAG ADDITIONAL BAD PIXELS
            if imageName == "fluxRaw" and associatedMask is not None:
                bpmArray = associatedMask
                if bpmArray is not None:
                    bpmArray = self._zoom_array(bpmArray, self.zoomFactor, self.dispersionAxis)
                    bpmRebinned = self._rebin_2d(bpmArray[SPIndex], *self.rebinShape)
                    bpmRebinned = self._sigma_clip_bpm(rebinned, bpmRebinned, order=order)
                    orderTable["bpMask"] = list(bpmRebinned)
                    self._cache_image_names.add("bpMask")

        self.log.debug('completed the ``cache_image`` method')
        return None


    def _zoom_array(self, arr, zoomFactor, dispersionAxis):
        """Zoom a 2d array along the spatial (slit) axis using nearest-neighbour pixel replication.

        For x-dispersion the slit runs along columns (axis=1); for y-dispersion it runs along rows (axis=0).
        """
        import numpy as np

        if dispersionAxis == "x":
            return np.repeat(arr, zoomFactor, axis=1)
        return np.repeat(arr, zoomFactor, axis=0)


    def _rebin_2d(self, arr, binRows, binCols):
        """Rebin a 2d array to lower resolution by block-averaging groups of pixels.

        ``binRows`` and ``binCols`` must evenly divide the corresponding array dimension.
        """
        new_shape = (arr.shape[0] // binRows, binRows, arr.shape[1] // binCols, binCols)
        return arr.reshape(new_shape).mean(axis=(1, 3))


    def _sigma_clip_bpm(self, rawFluxArray, bpmArray, order=None):
        """Sigma-clip a rebinned raw flux array to flag additional outlier pixels.

        Uses the rebinned raw flux as the data array and clips along the spatial axis
        (axis=1, i.e. across the slit) so that column-wise outliers are caught per
        cross-dispersion slice.
        """
        import numpy as np
        from astropy.stats import sigma_clip

        bpmArray[bpmArray > 0] = 1

        # MASK THE RAW FLUX WITH THE CURRENT BAD-PIXEL MAP THEN SIGMA-CLIP ALONG THE SLIT
        rawFluxMasked = np.ma.array(rawFluxArray, mask=bpmArray)
        rawFluxMasked = sigma_clip(
            rawFluxMasked,
            sigma_lower=2,
            sigma_upper=7,
            maxiters=1,
            cenfunc="mean",
            stdfunc="std",
            axis=1,
        )
        newBpm = rawFluxMasked.mask

        self.log.print(
            f"\t\t{sum(sum(newBpm)) - sum(sum(bpmArray))} additional bad pixels found from initial sigma clipping (order:{order})",
        )

        return newBpm

    def _setup_order_table_dataframes(self):
        """*Setup the individual order dataframes for each order in the order pixel table*
        """
        self.log.debug('starting the ``_setup_order_table_dataframes`` method')
        import numpy as np

        self.orderSlices = []
        self.wlMinMax = []
        self.subPixelIndexes = []

        if self.axisA == "x":
            axisALen = self.twoDMap["WAVELENGTH"].data.shape[1] * self.zoomFactor
            axisBLen = self.twoDMap["WAVELENGTH"].data.shape[0]
        else:
            axisALen = self.twoDMap["WAVELENGTH"].data.shape[0] * self.zoomFactor
            axisBLen = self.twoDMap["WAVELENGTH"].data.shape[1]

        # ITERATE OVER EACH ORDER, CLIPPING TO THE PIXEL BOUNDS (AMIN/AMAX) DEFINED FOR THAT ORDER
        for order, amin, amax, wlmin, wlmax in zip(self.orderNums, self.amins, self.amaxs, self.waveLengthMin, self.waveLengthMax):
            if order not in self.uniqueOrders:
                continue

            # FILTER THE ORDER PIXEL TABLE TO ROWS BELONGING TO THIS ORDER WITHIN THE PIXEL BOUNDS
            orderTable = self.orderPixelTable.loc[
                (self.orderPixelTable["order"] == order)
                & (self.orderPixelTable[f"{self.axisB}coord"] > amin)
                & (self.orderPixelTable[f"{self.axisB}coord"] < amax)
            ].copy()

            # COMPUTE SLIT UPPER AND LOWER SUB-PIXEL BOUNDS FOR EVERY CROSS-DISPERSION SLICE
            axisAstart = (
                np.round((orderTable[f"{self.axisA}coord_centre"] * self.zoomFactor)).astype(int)
                - self.slitHalfLength * self.zoomFactor
            )
            axisAstop = (
                np.round((orderTable[f"{self.axisA}coord_centre"] * self.zoomFactor)).astype(int)
                + self.slitHalfLength * self.zoomFactor
            )

            # MAKE SURE THE DISPERSION (AXIS B) PIXEL COORDINATES ARE INTEGERS
            axisBcoord = orderTable[f"{self.axisB}coord"].round().astype(int)

            validRows = (
                (axisAstart >= 0)
                & (axisAstop <= axisALen)
                & (axisBcoord >= 0)
                & (axisBcoord < axisBLen)
            )

            if not validRows.any():
                continue

            orderTable = orderTable.loc[validRows].copy()
            axisAstart = axisAstart.loc[validRows]
            axisAstop = axisAstop.loc[validRows]
            axisBcoord = axisBcoord.loc[validRows]

            # FOR EACH CROSS-DISPERSION SLICE, BUILD THE ZOOMED SPATIAL PIXEL INDEX LIST
            axisAcoords = list(map(lambda x: list(range(x[0], x[1])), zip(axisAstart, axisAstop)))

            # ASSIGN THE SAME DISPERSION PIXEL COORDINATE TO EVERY PIXEL IN EACH SLICE
            axisBcoords = list(map(lambda x: [x] * self.slitHalfLength * 2 * self.zoomFactor, axisBcoord))

            # SET NUMPY FANCY-INDEX ORDER: ROWS ARE ALWAYS Y, COLUMNS ALWAYS X
            # X-DISPERSION → B=y (rows), A=x (cols); Y-DISPERSION → A=y (rows), B=x (cols)
            if self.dispersionAxis == "x":
                SPIndex = (axisBcoords, axisAcoords)
            else:
                SPIndex = (axisAcoords, axisBcoords)
            self.subPixelIndexes.append(SPIndex)
            self.sliceOrders.append(order)

            # COLLECT THE PROCESSED ORDER SLICE AND ITS WAVELENGTH RANGE FOR DOWNSTREAM EXTRACTION
            self.orderSlices.append(orderTable)
            self.wlMinMax.append((wlmin, wlmax))

        self.log.debug('completed the ``_setup_order_table_dataframes`` method')

        return None

    def get_order_slices(self):
        return self.orderSlices
    
    def get_order_wavelength_ranges(self):
        return self.wlMinMax

    def get_order_rectified(self):
        """*Return the rectified wavelength and slit position images for each order*

        **Return:**

            - ``orderRectifiedImages`` -- list of tuples of the form ``(wlImage, slitImage)`` for each order in ``self.orderSlices``
        """
        import numpy as np 
        self.log.debug('starting the ``get_order_rectified`` method')

        orderRectifiedImages = []
        for orderTable in self.orderSlices:
            order = orderTable["order"].iloc[0]
            rectifiedImageDict = {}
            for imageName in self._cache_image_names:
                rectifiedImageDict[imageName] = np.vstack(orderTable[imageName])

                if False:
                    import matplotlib
                    matplotlib.use("MacOSX")
                    import matplotlib.pyplot as plt
                    fig = plt.figure(
                        num=None,
                        figsize=(135, 1),
                        dpi=None,
                        facecolor=None,
                        edgecolor=None,
                        frameon=True,
                    )
                    fig.suptitle(f"{imageName}, order {order}", fontsize=16)
                    plt.imshow(rectifiedImageDict[imageName].T, interpolation="none", aspect="auto")
                    plt.show()

            orderRectifiedImages.append(rectifiedImageDict)

        self.log.debug('completed the ``get_order_rectified`` method')
        return orderRectifiedImages

