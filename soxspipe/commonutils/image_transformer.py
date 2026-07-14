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
    """*rectify frame arrays, order-by-order, into cross-dispersion slices ready for spectral extraction*

    For each order:

    1. Define the slit-position (s) and wavelength (w) bounds of the order.

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
            mapDF,
            orderPixelTable,
            twoDMapPath,
            dispersionMap,
            associatedFrame,
            slitHalfLength,
            edgeSamples=1
    ):
        super(image_transformer, self).__init__(log, settings, associatedFrame=associatedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath)

        self.log.debug("Starting the image_transformer object")
        self.orderPixelTable = orderPixelTable
        self.twoDMapPath = twoDMapPath
        self.slitHalfLength = slitHalfLength
        self.mapDF = mapDF

        # NUMBER OF BOUNDARY SAMPLE POINTS PER CELL EDGE — FIXED FOR THE LIFE OF THE INSTANCE SO THE
        # RESAMPLING GEOMETRY CAN BE PRECOMPUTED ONCE AND SHARED ACROSS ALL cache_image CALLS
        self.edgeSamples = edgeSamples

        # TODO: MAKE THIS A SETTING IN THE YAML FILE
        self.slitLengthArcsec = 3.0
        self.zoomFactor = 5

        # ORDERS PRESENT IN THE TRACE TABLE — USED TO SKIP ORDERS WITH NO DETECTED TRACE
        # self.orderPixelTable = self.orderPixelTable.loc[self.orderPixelTable["order"] == 16]
        self.uniqueOrders = self.orderPixelTable["order"].unique()

        # DETERMINE THE BOUNDS OF EACH ORDER IN WS-PIXEL SPACE
        self.orderSlitEdges, self.orderWlEdges = self._determine_rectified_image_boundaries()

        # DETECTOR SHAPE — SAME FOR EVERY NDARRAY EVER PASSED TO cache_image, SO ONLY DERIVED ONCE
        self.ny, self.nx = self.twoDMap["WAVELENGTH"].data.shape

        self._cache_image_names = set()

        # PRECOMPUTE THE PIXEL-BOUNDARY/POLYGON-AREA RESAMPLING WEIGHTS ONCE, SHARED BY EVERY cache_image CALL
        self._resamplingWeights = self._precompute_resampling_weights()

        # CACHE THE WAVELENGTH AND SLIT POSITION MAPS FOR LATER RECTIFICATION —
        # BUILT ANALYTICALLY FROM EACH ORDER'S OWN BIN EDGES, NOT RESAMPLED FROM self.twoDMap PIXEL DATA
        self._cache_true_wavelength_slit_images()

        return None

    def cache_image(
            self,
            imageName,
            ndarray,
            associatedMask=None,
            returnCoverage=False,
            debug=False):
        """
        *Place an image in the transformer's cache. These images can be acted on later.*

        **Key Arguments:**

            - ``imageName`` -- the unique name to give to the image
            - ``ndarray`` -- the 2D frame to rectify and cache
            - ``associatedMask`` -- 2D bad-pixel mask. Currently accepted but not applied (no-op), matching the previous behaviour. Default *None*
            - ``returnCoverage`` -- if True, the method will return a list (one per order) of coverage maps of the rectified image. Default *False*
            - ``debug`` -- if True, print and plot the rectified image for each order. Default *False*

        **Return:**

            - ``orderCoverage`` -- list of per-order coverage maps if ``returnCoverage`` is True, else None
        """
        self.log.debug('starting the ``cache_image`` method')

        import numpy as np
        bpmArray = associatedMask

        orderCoverage = [] if returnCoverage else None

        for order, sp_edges, wl_edges, orderTable in zip(self.uniqueOrders, self.orderSlitEdges, self.orderWlEdges, self.orderSlices):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1
            weights = self._resamplingWeights[order]

            # FULLY VECTORIZED WEIGHTED SUM USING THE PRECOMPUTED PIXEL/POLYGON-AREA WEIGHTS
            weighted = ndarray[weights["py"], weights["px"]] * weights["area"]
            flux = np.bincount(weights["flatIdx"], weights=weighted, minlength=n_sp * n_wl).reshape(n_sp, n_wl)

            orderTable[imageName] = list(flux.T)
            # SCALAR BROADCASTS ONCE THE TABLE'S ROW COUNT IS ESTABLISHED — READ BY get_order_rectified()
            orderTable["order"] = order
            self._cache_image_names.add(imageName)

            if returnCoverage:
                orderCoverage.append(weights["coverage"])

            if bpmArray is not None:
                # FULLY VECTORIZED WEIGHTED SUM USING THE PRECOMPUTED PIXEL/POLYGON-AREA WEIGHTS
                weightedBpm = bpmArray[weights["py"], weights["px"]] * weights["area"]
                bpm = np.bincount(weights["flatIdx"], weights=weightedBpm, minlength=n_sp * n_wl).reshape(n_sp, n_wl)
                bpm = bpm > 0
                orderTable[f"bpMask"] = list(bpm.T)
                self._cache_image_names.add("bpMask")

            if debug:
                print(f"Rectifying image '{imageName}' for order {order} with shape {ndarray.shape} into ({n_sp}, {n_wl})")
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
                # ROWS ARE SLIT POSITION (Y-AXIS), COLUMNS ARE WAVELENGTH (X-AXIS)
                plt.imshow(
                    flux,
                    interpolation="none",
                    aspect="auto",
                    origin="lower",
                    extent=[wl_edges[0], wl_edges[-1], sp_edges[0], sp_edges[-1]],
                )
                plt.xlabel("Wavelength (Å)")
                plt.ylabel("Slit Position (arcsec)")
                plt.show()

        self.log.debug('completed the ``cache_image`` method')
        return orderCoverage

    def _cache_true_wavelength_slit_images(self):
        """*Cache the analytic (non-resampled) wavelength and slit-position images for each order*

        Unlike ``cache_image``, this does not resample any detector-pixel-space ndarray. Instead it
        directly computes, for every output (slit, wavelength) cell, the true centre wavelength/slit
        value implied by that order's own ``sp_edges``/``wl_edges`` bin boundaries. The wavelength
        image is constant along each row (varies only with the wavelength bin); the slit image is
        constant along each column (varies only with the slit bin).
        """
        self.log.debug('starting the ``_cache_true_wavelength_slit_images`` method')

        import numpy as np

        cache_image_names = self._cache_image_names

        for order, sp_edges, wl_edges, orderTable in zip(self.uniqueOrders, self.orderSlitEdges, self.orderWlEdges, self.orderSlices):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1

            # BIN-CENTRE VALUES ALONG EACH AXIS, DERIVED PURELY FROM THE EDGE ARRAYS
            wl_centers = (wl_edges[:-1] + wl_edges[1:]) / 2.0
            sp_centers = (sp_edges[:-1] + sp_edges[1:]) / 2.0

            # BROADCAST VIEWS INSTEAD OF MATERIALISING FULL TILED ARRAYS
            wavelengthImage = np.broadcast_to(wl_centers, (n_sp, n_wl))
            slitImage = np.broadcast_to(sp_centers[:, None], (n_sp, n_wl))

            orderTable["wavelength"] = [row for row in wavelengthImage.T]
            orderTable["slit"] = [row for row in slitImage.T]
            orderTable["order"] = order

        cache_image_names.add("wavelength")
        cache_image_names.add("slit")

        self.log.debug('completed the ``_cache_true_wavelength_slit_images`` method')
        return None

    def _precompute_resampling_weights(self):
        """*Precompute, once per instance, the detector-pixel/output-cell polygon-overlap weights used to rectify any cached image*

        For every order and every (slit, wavelength) output cell, determine which detector pixels are overlapped by
        that cell's boundary polygon (converted to detector x,y via ``dispersion_map_to_pixel_arrays``) and by how
        much (clipped polygon area). This geometry is independent of the ndarray content being cached, so it is
        computed once here and reused by every ``cache_image`` call.

        **Return:**

        - ``resamplingWeights`` -- dict keyed by order number, each value a dict of numpy arrays ``i``, ``j``, ``px``, ``py``, ``area``, ``flatIdx`` and a per-cell ``coverage`` grid
        """
        self.log.debug('starting the ``_precompute_resampling_weights`` method')

        import numpy as np
        import pandas as pd
        from .dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays

        ncorners = 4 * self.edgeSamples

        # BUILD ONE FLAT TABLE OF BOUNDARY CORNER POINTS ACROSS ALL ORDERS AND CELLS
        orderChunks, wlChunks, spChunks = [], [], []
        records = []

        for order, sp_edges, wl_edges in zip(self.uniqueOrders, self.orderSlitEdges, self.orderWlEdges):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1
            spBlock, wlBlock = _pixel_boundaries_grid(sp_edges, wl_edges, self.edgeSamples)
            nrows = n_sp * n_wl * ncorners
            spChunks.append(spBlock.reshape(-1))
            wlChunks.append(wlBlock.reshape(-1))
            orderChunks.append(np.full(nrows, order))
            records.append({"order": order, "n_sp": n_sp, "n_wl": n_wl, "nrows": nrows})

        print(11)

        cornersDF = pd.DataFrame({
            "order": np.concatenate(orderChunks),
            "wavelength": np.concatenate(wlChunks),
            "slit_position": np.concatenate(spChunks),
        })

        print(12)

        cornersDF = cornersDF.astype({"order": int, "wavelength": float, "slit_position": float})

        print(13)

        # SINGLE BATCHED CONVERSION OF ALL BOUNDARY CORNER POINTS FROM WAVELENGTH/SLIT/ORDER TO DETECTOR X,Y
        # removeOffDetectorLocation=False: EVERY CORNER (EVEN OFF-DETECTOR) MUST BE KEPT TO PRESERVE POLYGON GEOMETRY AND ROW ORDER
        resultDF = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersionMap,
            orderPixelTable=cornersDF,
            removeOffDetectorLocation=False,
            trimColumns=True,
        )
        fit_x = resultDF["fit_x"].to_numpy()
        fit_y = resultDF["fit_y"].to_numpy()  

        print(14)      

        # REBUILD THE PER-ORDER RESAMPLING WEIGHTS FROM THE FLAT BOUNDARY CORNER TABLE
        resamplingWeights = {}
        offset = 0
        print(1)
        for record in records:
            order = record["order"]
            n_sp = record["n_sp"]
            n_wl = record["n_wl"]
            nrows = record["nrows"]

            print(order)

            xb = fit_x[offset:offset + nrows].reshape(n_sp, n_wl, ncorners)
            yb = fit_y[offset:offset + nrows].reshape(n_sp, n_wl, ncorners)
            offset += nrows

            iList, jList, pxList, pyList, areaList = [], [], [], [], []
            coverage = np.zeros((n_sp, n_wl))
            print(2)
            for i in range(n_sp):
                for j in range(n_wl):
                    print(i, j)
                    xs, ys = xb[i, j], yb[i, j]
                    if not (np.all(np.isfinite(xs)) and np.all(np.isfinite(ys))):
                        continue

                    poly = list(zip(xs, ys))

                    # DETECTOR PIXELS POSSIBLY OVERLAPPED BY THIS POLYGON
                    px_lo = max(int(np.floor(xs.min() + 0.5)), 0)
                    px_hi = min(int(np.floor(xs.max() + 0.5)), self.nx - 1)
                    py_lo = max(int(np.floor(ys.min() + 0.5)), 0)
                    py_hi = min(int(np.floor(ys.max() + 0.5)), self.ny - 1)

                    for py in range(py_lo, py_hi + 1):
                        for px in range(px_lo, px_hi + 1):
                            a = _polygon_area(_clip_to_pixel(poly, px, py))
                            if a > 0.0:
                                iList.append(i)
                                jList.append(j)
                                pxList.append(px)
                                pyList.append(py)
                                areaList.append(a)
                                coverage[i, j] += a

            iArr = np.array(iList, dtype=np.int64)
            jArr = np.array(jList, dtype=np.int64)
            resamplingWeights[order] = {
                "i": iArr,
                "j": jArr,
                "px": np.array(pxList, dtype=np.int64),
                "py": np.array(pyList, dtype=np.int64),
                "area": np.array(areaList, dtype=np.float64),
                "flatIdx": iArr * n_wl + jArr,
                "coverage": coverage,
            }
            print(6)

        self.log.debug('completed the ``_precompute_resampling_weights`` method')
        return resamplingWeights


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

    def _determine_rectified_image_boundaries(self):
        """*Setup the individual order dataframes for each order in the order pixel table*
        """
        self.log.debug('starting the ``_determine_rectified_image_boundaries`` method')
        import numpy as np
        import pandas as pd

        # FOR EACH CONTINUUM-FITTED DATA-POINT, RETURN THE SLIT POSITION 
        aArray = self.orderPixelTable[f"{self.axisA}coord_centre"].round().astype(int)
        bArray = self.orderPixelTable[f"{self.axisB}coord"]

        mapLookup = self.mapDF.set_index([f"{self.axisB}", f"{self.axisA}"])["slit_position"]
        self.orderPixelTable["slit_position"] = mapLookup.reindex(list(zip(aArray, bArray))).to_numpy()
        slitCentreArcsec = np.nanmean(self.orderPixelTable["slit_position"])

        subPixelSize = 0.25 / self.zoomFactor
        slitStart = slitCentreArcsec - self.slitLengthArcsec/2
        slitStop = slitCentreArcsec + self.slitLengthArcsec/2
        slitEdges = np.arange(slitStart, slitStop, subPixelSize)

        orderSlitEdges, orderWlEdges = [], []
        sliceOrders = []
        orderSliceTables = []

        # ITERATE OVER EACH ORDER, CLIPPING TO THE PIXEL BOUNDS (AMIN/AMAX) DEFINED FOR THAT ORDER
        wlMinMax = []

        for order, amin, amax, wlmin, wlmax in zip(self.orderNums, self.amins, self.amaxs, self.waveLengthMin, self.waveLengthMax):
            if order not in self.uniqueOrders:
                continue
            wl_edges = np.arange(wlmin, wlmax, 0.02/self.zoomFactor)
            orderSlitEdges.append(slitEdges)
            orderWlEdges.append(wl_edges)
            sliceOrders.append(order)
            wlMinMax.append((wlmin, wlmax))
            # print(f"Order {order}: Wavelength range: {wlmin:.2f} - {wlmax:.2f} Angstroms, Slit range: {slitStart/10:.1f} - {slitStop/10:.1f} arcsec")
            orderSliceTables.append(pd.DataFrame())

        self.orderSlices = orderSliceTables
        self.uniqueOrders = sliceOrders
        self.wlMinMax = wlMinMax

        self.log.debug('completed the ``_determine_rectified_image_boundaries`` method')
        return orderSlitEdges, orderWlEdges


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
        for orderTable, sp_edges, wl_edges in zip(self.orderSlices, self.orderSlitEdges, self.orderWlEdges):
            order = orderTable["order"].iloc[0]
            rectifiedImageDict = {}
            for imageName in self._cache_image_names:
                rectifiedImageDict[imageName] = np.vstack(orderTable[imageName]).T

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
                    # ROWS ARE SLIT POSITION (Y-AXIS), COLUMNS ARE WAVELENGTH (X-AXIS)
                    plt.imshow(
                        rectifiedImageDict[imageName],
                        interpolation="none",
                        aspect="auto",
                        origin="lower",
                        extent=[wl_edges[0], wl_edges[-1], sp_edges[0], sp_edges[-1]],
                    )
                    plt.xlabel("Wavelength (Å)")
                    plt.ylabel("Slit Position (arcsec)")
                    plt.show()

            orderRectifiedImages.append(rectifiedImageDict)

        self.log.debug('completed the ``get_order_rectified`` method')
        return orderRectifiedImages
    
# ----------------------------------------------------------------------
# Geometry helpers: Sutherland-Hodgman clipping against a unit pixel
# ----------------------------------------------------------------------
 
def _clip_halfplane(poly, axis, value, keep_below):
    """Clip polygon (list of (x, y)) against one axis-aligned half-plane."""
    out = []
    n = len(poly)
    for k in range(n):
        p = poly[k]
        q = poly[(k + 1) % n]
        p_in = (p[axis] <= value) if keep_below else (p[axis] >= value)
        q_in = (q[axis] <= value) if keep_below else (q[axis] >= value)
        if p_in:
            out.append(p)
        if p_in != q_in:                       # edge crosses the boundary
            t = (value - p[axis]) / (q[axis] - p[axis])
            out.append((p[0] + t * (q[0] - p[0]),
                        p[1] + t * (q[1] - p[1])))
    return out
 
    
def _clip_to_pixel(poly, px, py):
    """Clip polygon to the unit square of detector pixel (px, py)."""
    poly = _clip_halfplane(poly, 0, px - 0.5, keep_below=False)
    if not poly:
        return poly
    poly = _clip_halfplane(poly, 0, px + 0.5, keep_below=True)
    if not poly:
        return poly
    poly = _clip_halfplane(poly, 1, py - 0.5, keep_below=False)
    if not poly:
        return poly
    return _clip_halfplane(poly, 1, py + 0.5, keep_below=True)
 
 
def _polygon_area(poly):
    """Unsigned area via the shoelace formula."""
    if len(poly) < 3:
        return 0.0
    a = 0.0
    n = len(poly)
    for k in range(n):
        x0, y0 = poly[k]
        x1, y1 = poly[(k + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * abs(a)
    
def _pixel_boundaries_grid(sp_edges, wl_edges, edge_samples):
    """
    Boundary of every output-pixel cell in (SP, WL) space for a full
    ``sp_edges``/``wl_edges`` grid, traversed counter-clockwise, with
    `edge_samples` points per edge (>=1). Extra points let a curved
    mapping be followed accurately.

    VECTORIZED OVER THE WHOLE (n_sp, n_wl) GRID IN ONE SHOT — REPLACES A
    PER-CELL PYTHON LOOP CALLING THIS ONCE PER CELL, WHICH DOMINATED THE
    COST OF _precompute_resampling_weights AT REALISTIC GRID SIZES
    (n_sp * n_wl IN THE MILLIONS ACROSS ALL ORDERS).

    Returns ``spBlock, wlBlock`` of shape ``(n_sp, n_wl, 4 * edge_samples)``.
    """
    import numpy as np

    sp0, sp1 = sp_edges[:-1, None], sp_edges[1:, None]   # (n_sp, 1)
    wl0, wl1 = wl_edges[None, :-1], wl_edges[None, 1:]     # (1, n_wl)

    t = np.linspace(0.0, 1.0, edge_samples + 1)[:-1]      # exclude endpoint, (edge_samples,)

    n_sp = sp_edges.shape[0] - 1
    n_wl = wl_edges.shape[0] - 1

    shape = (n_sp, n_wl, edge_samples)

    # BROADCAST EACH EDGE OF THE CELL BOUNDARY TO SHAPE (n_sp, n_wl, edge_samples)
    bottom_sp = np.broadcast_to(sp0[:, :, None] + t * (sp1 - sp0)[:, :, None], shape)   # wl = wl0
    bottom_wl = np.broadcast_to(wl0[:, :, None], shape)

    right_sp = np.broadcast_to(sp1[:, :, None], shape)                                  # sp = sp1
    right_wl = np.broadcast_to(wl0[:, :, None] + t * (wl1 - wl0)[:, :, None], shape)

    top_sp = np.broadcast_to(sp1[:, :, None] + t * (sp0 - sp1)[:, :, None], shape)      # wl = wl1
    top_wl = np.broadcast_to(wl1[:, :, None], shape)

    left_sp = np.broadcast_to(sp0[:, :, None], shape)                                   # sp = sp0
    left_wl = np.broadcast_to(wl1[:, :, None] + t * (wl0 - wl1)[:, :, None], shape)

    spBlock = np.concatenate([bottom_sp, right_sp, top_sp, left_sp], axis=2)
    wlBlock = np.concatenate([bottom_wl, right_wl, top_wl, left_wl], axis=2)
    return spBlock, wlBlock
