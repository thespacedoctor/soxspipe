#!/usr/bin/env python
"""
*Using a 2D dispersion image map, transform a SOXS data frame from xy detector pixel space to wavelength-slit position space.*

:Author:
    David Young

:Date Created:
    June 22, 2026
"""

import os

os.environ["TERM"] = "vt100"

import numba

# NUMPY/NUMBA IMPORTED AT MODULE SCOPE (NOT METHOD-LOCAL) SO THE @numba.njit
# DECORATORS BELOW ARE EVALUATED ONCE AT IMPORT TIME
import numpy as np

from .base_util import base_util

# DEGREE OF THE PER-ORDER POLYNOMIAL FIT OF SLIT POSITION VS WAVELENGTH USED TO
# CENTRE THE RECTIFIED SLIT GRID ON THE OBJECT TRACE (CAPPED TO nValid - 1 WHEN
# AN ORDER HAS FEWER TRACE POINTS THAN THIS DEGREE REQUIRES)
SLIT_CENTRE_POLY_DEGREE = 3

# DENSITY OF WAVELENGTH SAMPLES (PER DETECTOR PIXEL ALONG THE ORDER) USED TO TRACE EACH ORDER'S
# PATH ACROSS THE DETECTOR WHEN BUILDING PIXEL-UNIFORM WAVELENGTH BIN EDGES
TRACE_SAMPLES_PER_PIXEL = 4

# SLIT OFFSET (ARCSEC) FROM THE TRACE CENTRE USED TO MEASURE EACH ORDER'S ARCSEC-PER-PIXEL SCALE
SLIT_SCALE_PROBE_ARCSEC = 1.0


class image_transformer(base_util):
    """*rectify frame arrays, order-by-order, into cross-dispersion slices ready for spectral extraction*

    For each order:

    1. Define the slit-position (s) and wavelength (w) bounds of the order, centred on the
       object trace: the slit centre at each wavelength is a per-order polynomial fit of
       slit position vs wavelength through that order's own trace points.

       The wavelength bin edges are spaced evenly in detector pixels along the trace (so they
       are not evenly spaced in wavelength), and the slit bin edges use that order's own
       arcsec-per-pixel scale, both measured from the dispersion map.
    2. Resamples the frame onto an oversampled (slit, wavelength) grid, with the slit axis
       sub-sampled by ``zoomFactorSlit`` and the wavelength axis by ``zoomFactorWavelength``.
    3. Rebins that grid back to detector-pixel resolution by summing each
       ``zoomFactorSlit`` x ``zoomFactorWavelength`` block, so each rectified pixel holds the
       counts of about one detector pixel.
    4. Optionally sigma-clips the rebinned raw flux to catch additional outlier pixels before extraction.

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- settings dict
    - ``orderPixelTable`` -- dataframe containing per-pixel order-trace metadata (object continuum fit results)
    - ``twoDMapPath`` -- path to the 2D map FITS file (pixel wavelength and slit position values needed for rectification)
    - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
    - ``associatedFrame`` -- an example 2D frame to be rectified. This frame is used to determine detector binning, arm etc.
    - ``slitHalfLength`` -- half-length of the slit in detector pixels (sets extraction aperture). May end in .5.
      The rectified grid has ``2 * slitHalfLength`` rows (rounded half up), centred on the trace

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
        slitHalfLength=slitHalfLength
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
        self, log, settings, orderPixelTable, twoDMapPath, dispersionMap, associatedFrame, slitHalfLength, edgeSamples=1
    ):
        super().__init__(
            log, settings, associatedFrame=associatedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath
        )

        self.log.debug("Starting the image_transformer object")
        self.orderPixelTable = orderPixelTable
        self.twoDMapPath = twoDMapPath
        self.slitHalfLength = slitHalfLength

        # NUMBER OF BOUNDARY SAMPLE POINTS PER CELL EDGE — FIXED FOR THE LIFE OF THE INSTANCE SO THE
        # RESAMPLING GEOMETRY CAN BE PRECOMPUTED ONCE AND SHARED ACROSS ALL cache_image CALLS
        self.edgeSamples = edgeSamples
        # SUB-PIXEL SAMPLING FACTORS FOR THE SLIT (Y-AXIS) AND WAVELENGTH (X-AXIS) OF THE RECTIFIED IMAGE
        self.zoomFactorSlit = 5
        self.zoomFactorWavelength = 5

        # ORDERS PRESENT IN THE TRACE TABLE — USED TO SKIP ORDERS WITH NO DETECTED TRACE
        # self.orderPixelTable = self.orderPixelTable.loc[self.orderPixelTable["order"] == 16]
        self.uniqueOrders = self.orderPixelTable["order"].unique()

        # DETECTOR SHAPE — SAME FOR EVERY NDARRAY EVER PASSED TO cache_image, SO ONLY DERIVED ONCE
        self.ny, self.nx = self.twoDMap["WAVELENGTH"].data.shape

        # DETERMINE THE BOUNDS OF EACH ORDER IN WS-PIXEL SPACE
        self.orderSlitEdges, self.orderWlEdges = self._determine_rectified_image_boundaries()
        self._cache_image_names = set()

        # PRECOMPUTE THE PIXEL-BOUNDARY/POLYGON-AREA RESAMPLING WEIGHTS ONCE, SHARED BY EVERY cache_image CALL
        self._resamplingWeights = self._precompute_resampling_weights()

        # CACHE THE WAVELENGTH AND SLIT POSITION MAPS FOR LATER RECTIFICATION —
        # BUILT ANALYTICALLY FROM EACH ORDER'S OWN BIN EDGES, NOT RESAMPLED FROM self.twoDMap PIXEL DATA
        self._cache_true_wavelength_slit_images()

        return

    def cache_image(self, imageName, ndarray, associatedMask=None, returnCoverage=False, debug=False):
        """
        *Place an image in the transformer's cache. These images can be acted on later.*

        **Key Arguments:**

            - ``imageName`` -- the unique name to give to the image
            - ``ndarray`` -- the 2D frame to rectify and cache
            - ``associatedMask`` -- 2D bad-pixel mask. Flagged pixels are left out of the cell sums and each cell is
              renormalised by its good-pixel area.
              A cell whose flagged area exceeds 0.2 is masked in ``bpMask``. Default *None*
            - ``returnCoverage`` -- if True, the method will return a list (one per order) of coverage maps of the rectified image. Default *False*
            - ``debug`` -- if True, print and plot the rectified image for each order. Default *False*

        **Return:**

            - ``orderCoverage`` -- list of per-order coverage maps if ``returnCoverage`` is True, else None
        """
        self.log.debug("starting the ``cache_image`` method")

        bpmArray = associatedMask
        badPixels = None if bpmArray is None else np.asarray(bpmArray).astype(bool)

        orderCoverage = [] if returnCoverage else None

        for order, sp_edges, wl_edges, orderTable in zip(
            self.uniqueOrders, self.orderSlitEdges, self.orderWlEdges, self.orderSlices
        ):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1
            weights = self._resamplingWeights[order]

            # FULLY VECTORIZED WEIGHTED SUM USING THE PRECOMPUTED (ALREADY REBINNED) PIXEL-OVERLAP WEIGHTS,
            # LEAVING OUT FLAGGED PIXELS AND RENORMALISING EACH CELL BY ITS GOOD-PIXEL AREA
            flux = _apply_weights(ndarray, weights, badPixels=badPixels)

            orderTable[imageName] = list(flux.T)
            # SCALAR BROADCASTS ONCE THE TABLE'S ROW COUNT IS ESTABLISHED — READ BY get_order_rectified()
            orderTable["order"] = order
            self._cache_image_names.add(imageName)

            if returnCoverage:
                orderCoverage.append(weights["coverage"])

            if bpmArray is not None:
                bpm = _apply_weights(badPixels, weights) > 0.2
                orderTable["bpMask"] = list(bpm.T)
                self._cache_image_names.add("bpMask")

            if debug:
                print(
                    f"Rectifying image '{imageName}' for order {order} with shape {ndarray.shape} into ({n_sp}, {n_wl})"
                )
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
                plt.ylabel("Slit offset from trace (arcsec)")
                plt.show()

        self.log.debug("completed the ``cache_image`` method")
        return orderCoverage

    def _unzoom(self, arr2d, operation="sum"):
        """Flux-conserving NxM binning (slit rows by zoomFactorSlit, wavelength columns by zoomFactorWavelength)."""
        zs = self.zoomFactorSlit
        zw = self.zoomFactorWavelength
        by = arr2d.shape[0] // zs
        bx = arr2d.shape[1] // zw
        if by == 0 or bx == 0:
            return arr2d
        trimmed = arr2d[: by * zs, : bx * zw]
        if operation == "sum":
            return trimmed.reshape(by, zs, bx, zw).sum(axis=(1, 3))
        if operation == "mean":
            return trimmed.reshape(by, zs, bx, zw).mean(axis=(1, 3))
        raise ValueError("Invalid operation. Use 'sum' or 'mean'.")

    def _cache_true_wavelength_slit_images(self):
        """*Cache the analytic (non-resampled) wavelength and slit-position images for each order*

        Unlike ``cache_image``, this does not resample any detector-pixel-space ndarray. Instead it
        directly computes, for every output (slit, wavelength) cell, the true centre wavelength/slit
        value implied by that order's own ``sp_edges``/``wl_edges`` bin boundaries. The wavelength
        image is constant along each row (varies only with the wavelength bin); the slit image
        varies along both axes, since it is the bin's slit offset plus that order's trace centre
        at the bin's wavelength (``offset + np.polyval(orderSlitCentreCoeffs, wavelength)``).
        """
        self.log.debug("starting the ``_cache_true_wavelength_slit_images`` method")

        import numpy as np

        cache_image_names = self._cache_image_names

        for order, sp_edges, wl_edges, centreCoeffs, orderTable in zip(
            self.uniqueOrders,
            self.orderSlitEdges,
            self.orderWlEdges,
            self.orderSlitCentreCoeffs,
            self.orderSlices,
            strict=True,
        ):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1

            # BIN-CENTRE VALUES ALONG EACH AXIS, DERIVED PURELY FROM THE EDGE ARRAYS
            wl_centers = (wl_edges[:-1] + wl_edges[1:]) / 2.0
            sp_centers = (sp_edges[:-1] + sp_edges[1:]) / 2.0

            wavelengthImage = np.broadcast_to(wl_centers, (n_sp, n_wl))
            # SLIT OFFSET PLUS THE PER-WAVELENGTH TRACE CENTRE — A REAL 2D ARRAY, NOT A BROADCAST VIEW
            slitImage = sp_centers[:, None] + np.polyval(centreCoeffs, wl_centers)[None, :]

            wavelengthImage = self._unzoom(wavelengthImage, operation="mean")
            slitImage = self._unzoom(slitImage, operation="mean")

            orderTable["wavelength"] = [row for row in wavelengthImage.T]
            orderTable["slit"] = [row for row in slitImage.T]
            orderTable["order"] = order

        cache_image_names.add("wavelength")
        cache_image_names.add("slit")

        self.log.debug("completed the ``_cache_true_wavelength_slit_images`` method")
        return

    def _precompute_resampling_weights(self):
        """*Precompute, once per instance, the detector-pixel/output-cell polygon-overlap weights used to rectify any cached image*

        For every order and every (slit, wavelength) output cell, determine which detector pixels are overlapped by
        that cell's boundary polygon (converted to detector x,y via ``dispersion_map_to_pixel_arrays``) and by how
        much (clipped polygon area). This geometry is independent of the ndarray content being cached, so it is
        computed once here and reused by every ``cache_image`` call.

        **Return:**

        - ``resamplingWeights`` -- dict keyed by order number. Each value holds the rebinned-grid weights
          (``px``, ``py``, ``area``, ``flatIdx``, ``shape``) and the zoomed-grid ``coverage`` array
        """
        self.log.debug("starting the ``_precompute_resampling_weights`` method")

        import numpy as np
        import pandas as pd

        from .dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays

        ncorners = 4 * self.edgeSamples

        # BUILD ONE FLAT TABLE OF BOUNDARY CORNER POINTS ACROSS ALL ORDERS AND CELLS
        orderChunks, wlChunks, spChunks = [], [], []
        records = []

        for order, sp_edges, wl_edges, centreCoeffs in zip(
            self.uniqueOrders,
            self.orderSlitEdges,
            self.orderWlEdges,
            self.orderSlitCentreCoeffs,
            strict=True,
        ):
            n_sp = len(sp_edges) - 1
            n_wl = len(wl_edges) - 1
            spBlock, wlBlock = _pixel_boundaries_grid(sp_edges, wl_edges, self.edgeSamples)
            # CONVERT SLIT OFFSETS FROM THE TRACE TO ABSOLUTE SLIT POSITIONS BEFORE THEY ARE
            # HANDED TO THE WAVELENGTH/SLIT -> DETECTOR-PIXEL CONVERTER
            spBlock = spBlock + np.polyval(centreCoeffs, wlBlock)
            nrows = n_sp * n_wl * ncorners
            spChunks.append(spBlock.reshape(-1))
            wlChunks.append(wlBlock.reshape(-1))
            orderChunks.append(np.full(nrows, order))
            records.append({"order": order, "n_sp": n_sp, "n_wl": n_wl, "nrows": nrows})

        cornersDF = pd.DataFrame(
            {
                "order": np.concatenate(orderChunks),
                "wavelength": np.concatenate(wlChunks),
                "slit_position": np.concatenate(spChunks),
            }
        )

        cornersDF = cornersDF.astype({"order": int, "wavelength": float, "slit_position": float})

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

        # REBUILD THE PER-ORDER RESAMPLING WEIGHTS FROM THE FLAT BOUNDARY CORNER TABLE.
        # THE PER-CELL PIXEL-OVERLAP CLIP/AREA WORK ITSELF RUNS IN THE NUMBA-JIT-COMPILED
        # _resample_weights_kernel (MILLIONS OF CELLS ACROSS ALL ORDERS) — EVERYTHING HERE
        # IS EITHER A ONE-OFF VECTORIZED NUMPY PRE-PASS OR O(NUMBER OF ORDERS) BOOKKEEPING
        resamplingWeights = {}
        offset = 0
        for record in records:
            order = record["order"]
            n_sp = record["n_sp"]
            n_wl = record["n_wl"]
            nrows = record["nrows"]

            xb = np.ascontiguousarray(fit_x[offset : offset + nrows].reshape(n_sp, n_wl, ncorners))
            yb = np.ascontiguousarray(fit_y[offset : offset + nrows].reshape(n_sp, n_wl, ncorners))
            offset += nrows

            # VECTORIZED PER-CELL CANDIDATE-PIXEL BOUNDING BOX (REPLACES THE OLD PER-CELL
            # np.isfinite()/.min()/.max() PYTHON LOOP)
            valid = np.all(np.isfinite(xb), axis=2) & np.all(np.isfinite(yb), axis=2)
            xmin = np.where(valid, np.nanmin(xb, axis=2), 0.0)
            xmax = np.where(valid, np.nanmax(xb, axis=2), -1.0)
            ymin = np.where(valid, np.nanmin(yb, axis=2), 0.0)
            ymax = np.where(valid, np.nanmax(yb, axis=2), -1.0)

            px_lo = np.clip(np.floor(xmin + 0.5).astype(np.int32), 0, self.nx - 1)
            px_hi = np.clip(np.floor(xmax + 0.5).astype(np.int32), 0, self.nx - 1)
            py_lo = np.clip(np.floor(ymin + 0.5).astype(np.int32), 0, self.ny - 1)
            py_hi = np.clip(np.floor(ymax + 0.5).astype(np.int32), 0, self.ny - 1)
            # px_hi < px_lo IS THE "SKIP THIS CELL" SENTINEL CONSUMED BY THE KERNEL
            # (MARKS NON-FINITE/OFF-BOUNDARY-DEGENERATE CELLS)
            px_hi = np.where(valid, px_hi, px_lo - 1)

            # SAFE UPPER BOUND ON KERNEL OUTPUT LENGTH: ONLY CANDIDATE PIXELS WITH
            # POSITIVE CLIPPED OVERLAP AREA ARE ACTUALLY KEPT
            candCount = np.where(valid, (px_hi - px_lo + 1) * (py_hi - py_lo + 1), 0).astype(np.int32)
            maxTotal = int(candCount.sum())

            iOut = np.empty(maxTotal, dtype=np.int32)
            jOut = np.empty(maxTotal, dtype=np.int32)
            pxOut = np.empty(maxTotal, dtype=np.int32)
            pyOut = np.empty(maxTotal, dtype=np.int32)
            areaOut = np.empty(maxTotal, dtype=np.float64)
            coverage = np.zeros((n_sp, n_wl))

            nOut = _resample_weights_kernel(
                xb,
                yb,
                px_lo,
                px_hi,
                py_lo,
                py_hi,
                iOut,
                jOut,
                pxOut,
                pyOut,
                areaOut,
                coverage,
            )

            # MERGE THE SUB-CELL OVERLAPS INTO THE REBINNED (DETECTOR-RESOLUTION) GRID, SO EACH CACHED IMAGE
            # IS A SINGLE WEIGHTED SUM OF DETECTOR PIXELS
            resamplingWeights[order] = {
                **_rebin_resampling_weights(
                    iOut[:nOut],
                    jOut[:nOut],
                    pxOut[:nOut],
                    pyOut[:nOut],
                    areaOut[:nOut],
                    nSp=n_sp,
                    nWl=n_wl,
                    zoomSlit=self.zoomFactorSlit,
                    zoomWavelength=self.zoomFactorWavelength,
                    nx=self.nx,
                    ny=self.ny,
                ),
                "coverage": coverage,
            }

        self.log.debug("completed the ``_precompute_resampling_weights`` method")
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

    def _fit_slit_centre_polynomial(self, wavelength, slitPosition, degree, order, fallbackCentreArcsec):
        """*Fit slit position as a polynomial function of wavelength, degrading the degree on a rank-deficient fit*

        **Key Arguments:**

        - ``wavelength`` -- 1D array of finite wavelength values for this order's trace points
        - ``slitPosition`` -- 1D array of finite slit-position values, same length as ``wavelength``
        - ``degree`` -- requested polynomial degree (already capped by the number of distinct wavelengths)
        - ``order`` -- the order number, used only for logging
        - ``fallbackCentreArcsec`` -- constant centre to use if every degree down to 0 is degenerate

        **Return:**

        - ``centreCoeffs`` -- polynomial coefficients (``numpy.polyfit`` convention), highest power first
        - ``isFallback`` -- True if every degree was degenerate and ``centreCoeffs`` is the constant fallback centre
        """
        self.log.debug("starting the ``_fit_slit_centre_polynomial`` method")

        import warnings

        import numpy as np
        from numpy.exceptions import RankWarning

        # STEP THE DEGREE DOWN ON A RANK-DEFICIENT FIT (E.G. REPEATED/NEAR-IDENTICAL WAVELENGTHS)
        # RATHER THAN SILENTLY TRUSTING numpy's LEAST-SQUARES SOLUTION TO AN ILL-POSED SYSTEM
        for candidateDegree in range(degree, -1, -1):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("error", RankWarning)
                    return np.polyfit(wavelength, slitPosition, candidateDegree), False
            except (RankWarning, np.linalg.LinAlgError):
                continue

        # EVERY DEGREE DOWN TO A CONSTANT WAS DEGENERATE — FALL BACK TO THE GLOBAL MEAN
        self.log.warning(
            f"Slit-centre trace fit for order {order} was rank-deficient at every degree down to 0 "
            f"(likely duplicate wavelengths); falling back to the global mean slit position "
            f"({fallbackCentreArcsec:.3f} arcsec) as a constant centre."
        )
        return np.array([fallbackCentreArcsec]), True

    def _determine_rectified_image_boundaries(self):
        """*Setup the individual order dataframes for each order in the order pixel table*"""
        self.log.debug("starting the ``_determine_rectified_image_boundaries`` method")
        import numpy as np
        import pandas as pd

        # FOR EACH CONTINUUM-FITTED DATA-POINT, RETURN THE SLIT POSITION
        aArray = self.orderPixelTable[f"{self.axisA}coord_centre"].round().astype(int)
        bArray = self.orderPixelTable[f"{self.axisB}coord"]

        if self.dispersionAxis == "x":
            mapLookup = self.mapDF.set_index([f"{self.axisA}", f"{self.axisB}"])[["slit_position", "wavelength"]]
            self.orderPixelTable[["slit_position", "wavelength"]] = mapLookup.reindex(
                list(zip(aArray, bArray))
            ).to_numpy()
        else:
            mapLookup = self.mapDF.set_index([f"{self.axisB}", f"{self.axisA}"])[["slit_position", "wavelength"]]
            self.orderPixelTable[["slit_position", "wavelength"]] = mapLookup.reindex(
                list(zip(bArray, aArray))
            ).to_numpy()

        # GLOBAL FALLBACK CENTRE — ONLY USED FOR AN ORDER WITH NO VALID TRACE POINTS OF ITS OWN.
        # ISFINITE (NOT JUST NOTNA) SO A STRAY +/-INF IN THE MAP CAN'T SNEAK THROUGH THE MEAN
        finiteSlitPosition = self.orderPixelTable["slit_position"][np.isfinite(self.orderPixelTable["slit_position"])]
        globalSlitCentreArcsec = finiteSlitPosition.mean() if len(finiteSlitPosition) else np.nan
        if not np.isfinite(globalSlitCentreArcsec):
            raise ValueError(
                "Cannot determine any slit centre: every trace point's slit_position is missing "
                "or non-finite across all orders. Check the 2D map / orderPixelTable inputs."
            )
        # KEPT FOR THE SLIT-DRIFT QC PLOT (THE SINGLE CENTRE USED BEFORE PER-ORDER TRACES)
        self.globalSlitCentreArcsec = globalSlitCentreArcsec

        # from astropy.table import Table
        # t = Table.from_pandas(self.orderPixelTable)
        # t.write("/tmp/table.fits", overwrite=True)

        orderTraces = []
        fallbackFlags = []

        # ITERATE OVER EACH ORDER, CLIPPING TO THE PIXEL BOUNDS (AMIN/AMAX) DEFINED FOR THAT ORDER
        for order, amin, amax, wlmin, wlmax in zip(
            self.orderNums, self.amins, self.amaxs, self.waveLengthMin, self.waveLengthMax
        ):
            if order not in self.uniqueOrders:
                continue
            pixelRange = amax - amin

            # FIT THIS ORDER'S OWN TRACE POINTS: SLIT POSITION AS A POLYNOMIAL FUNCTION OF WAVELENGTH
            orderMask = self.orderPixelTable["order"] == order
            orderTrace = self.orderPixelTable.loc[orderMask, ["slit_position", "wavelength"]]
            # "VALID" MEANS FINITE ON BOTH COLUMNS — dropna() ALONE WOULD LET +/-INF THROUGH
            finiteMask = np.isfinite(orderTrace["slit_position"]) & np.isfinite(orderTrace["wavelength"])
            validTrace = orderTrace.loc[finiteMask]
            nValid = len(validTrace)

            if nValid == 0:
                self.log.warning(
                    f"No valid trace points found for order {order}; falling back to the global "
                    f"mean slit position ({globalSlitCentreArcsec:.3f} arcsec) as a constant centre."
                )
                centreCoeffs = np.array([globalSlitCentreArcsec])
                isFallback = True
            else:
                # DEGREE IS CAPPED BY THE NUMBER OF *DISTINCT* WAVELENGTHS, NOT JUST THE POINT
                # COUNT — REPEATED WAVELENGTHS AT DIFFERENT SLIT POSITIONS MAKE THE FIT
                # RANK-DEFICIENT EVEN WHEN nValid IS LARGE
                nUniqueWl = validTrace["wavelength"].nunique()
                degree = min(SLIT_CENTRE_POLY_DEGREE, nUniqueWl - 1)
                centreCoeffs, isFallback = self._fit_slit_centre_polynomial(
                    validTrace["wavelength"].to_numpy(),
                    validTrace["slit_position"].to_numpy(),
                    degree,
                    order,
                    globalSlitCentreArcsec,
                )

            fallbackFlags.append(isFallback)
            orderTraces.append((order, wlmin, wlmax, pixelRange, centreCoeffs))

        # WAVELENGTH EDGES EVENLY SPACED IN DETECTOR PIXELS ALONG EACH TRACE, PLUS EACH ORDER'S ARCSEC/PIXEL
        orderGeometry = self._measure_order_trace_geometry(orderTraces)

        # SLIT OFFSET EDGES, CENTRED ON ZERO — THE ABSOLUTE SLIT POSITION OF ANY POINT IS THIS OFFSET
        # PLUS THE PER-ORDER, PER-WAVELENGTH TRACE CENTRE (SEE orderSlitCentreCoeffs)
        # THE GRID HOLDS A WHOLE NUMBER OF DETECTOR-PIXEL ROWS (2 x HALF LENGTH, ROUNDED HALF UP), EDGES FROM -H TO +H
        # INCLUSIVE, SO THE ROWS ARE CENTRED ON THE TRACE AND EACH ROW IS EXACTLY zoomFactorSlit CELLS WIDE
        nSlitRows = int(np.floor(2 * self.slitHalfLength + 0.5))
        slitPixelOffsets = np.linspace(-nSlitRows / 2, nSlitRows / 2, nSlitRows * self.zoomFactorSlit + 1)
        orderSlitEdges = [slitPixelOffsets * arcsecPerPixel for _, arcsecPerPixel in orderGeometry]
        orderWlEdges = [wlEdges for wlEdges, _ in orderGeometry]

        self.orderSlices = [pd.DataFrame() for _ in orderTraces]
        self.uniqueOrders = [trace[0] for trace in orderTraces]
        self.wlMinMax = [(trace[1], trace[2]) for trace in orderTraces]
        self.orderSlitCentreCoeffs = [trace[4] for trace in orderTraces]
        self.orderSlitCentreFallback = fallbackFlags

        self.log.debug("completed the ``_determine_rectified_image_boundaries`` method")
        return orderSlitEdges, orderWlEdges

    def _measure_order_trace_geometry(self, orderTraces):
        """*Measure each order's pixel-uniform wavelength edges and arcsec-per-pixel scale from the dispersion map*

        The trace centre (and a point ``SLIT_SCALE_PROBE_ARCSEC`` along the slit from it) is sampled on a dense
        wavelength grid and converted to detector pixels in a single ``dispersion_map_to_pixel_arrays`` call.

        **Key Arguments:**

        - ``orderTraces`` -- list of ``(order, wlmin, wlmax, pixelRange, centreCoeffs)`` tuples

        **Return:**

        - ``orderGeometry`` -- list of ``(wlEdges, arcsecPerPixel)`` tuples, in the same order as ``orderTraces``
        """
        self.log.debug("starting the ``_measure_order_trace_geometry`` method")

        import numpy as np
        import pandas as pd

        from .dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays

        # EACH ORDER CONTRIBUTES ITS DENSE TRACE SAMPLES FOLLOWED BY THE SAME WAVELENGTHS AT THE SLIT PROBE
        chunks = []
        sampleCounts = []
        for order, wlmin, wlmax, pixelRange, centreCoeffs in orderTraces:
            nSamples = max(int(np.ceil(pixelRange * TRACE_SAMPLES_PER_PIXEL)), 2)
            wavelengths = np.linspace(wlmin, wlmax, nSamples)
            centre = np.polyval(centreCoeffs, wavelengths)
            chunks.append(
                pd.DataFrame(
                    {
                        "order": order,
                        "wavelength": np.concatenate([wavelengths, wavelengths]),
                        "slit_position": np.concatenate([centre, centre + SLIT_SCALE_PROBE_ARCSEC]),
                    }
                )
            )
            sampleCounts.append(nSamples)

        samplesDF = pd.concat(chunks, ignore_index=True).astype(
            {"order": int, "wavelength": float, "slit_position": float}
        )
        resultDF = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersionMap,
            orderPixelTable=samplesDF,
            removeOffDetectorLocation=False,
            trimColumns=True,
        )
        # FLOAT64 SO THE CUMULATIVE PATH LENGTH DOES NOT INHERIT THE CONVERTER'S FLOAT32 ROUNDING
        wavelength = resultDF["wavelength"].to_numpy(dtype=float)
        fitX = resultDF["fit_x"].to_numpy(dtype=float)
        fitY = resultDF["fit_y"].to_numpy(dtype=float)
        # SAMPLES OFF THE DETECTOR MUST NOT CONTRIBUTE TO THE TRACE PATH OR THE SLIT SCALE
        offDetector = (fitX < -0.5) | (fitX >= self.nx - 0.5) | (fitY < -0.5) | (fitY >= self.ny - 0.5)
        fitX[offDetector] = np.nan
        fitY[offDetector] = np.nan

        orderGeometry = []
        start = 0
        for nSamples in sampleCounts:
            trace = slice(start, start + nSamples)
            probe = slice(start + nSamples, start + 2 * nSamples)
            start += 2 * nSamples
            wlEdges = _pixel_uniform_wavelength_edges(
                wavelength[trace], fitX[trace], fitY[trace], samplesPerPixel=self.zoomFactorWavelength
            )
            arcsecPerPixel = _arcsec_per_pixel(
                SLIT_SCALE_PROBE_ARCSEC, fitX[trace], fitY[trace], fitX[probe], fitY[probe]
            )
            orderGeometry.append((wlEdges, arcsecPerPixel))

        self.log.debug("completed the ``_measure_order_trace_geometry`` method")
        return orderGeometry

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

        self.log.debug("starting the ``get_order_rectified`` method")

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

                    # Calculate sigma-clipped mean and std
                    from astropy.stats import sigma_clip

                    clipped_data = sigma_clip(rectifiedImageDict[imageName], sigma=3, maxiters=5)
                    clipped_mean = clipped_data.mean()
                    clipped_std = clipped_data.std()

                    plt.imshow(
                        rectifiedImageDict[imageName],
                        interpolation="none",
                        aspect="auto",
                        origin="lower",
                        extent=[wl_edges[0], wl_edges[-1], sp_edges[0], sp_edges[-1]],
                        vmin=clipped_mean - clipped_std,
                        vmax=clipped_mean + 7 * clipped_std,
                    )
                    plt.xlabel("Wavelength (Å)")
                    plt.ylabel("Slit offset from trace (arcsec)")
                    plt.show()

            orderRectifiedImages.append(rectifiedImageDict)

        self.log.debug("completed the ``get_order_rectified`` method")
        return orderRectifiedImages


# ----------------------------------------------------------------------
# Geometry kernels: Sutherland-Hodgman clipping against a unit pixel,
# NUMBA-JIT-COMPILED — SEE _precompute_resampling_weights FOR THE HOT
# LOOP THESE ARE CALLED FROM (MILLIONS OF CELLS PER image_transformer
# INSTANTIATION, SO THIS MUST RUN AT COMPILED SPEED, NOT PURE PYTHON)
# ----------------------------------------------------------------------


@numba.njit(cache=True, inline="always")
def _clip_halfplane_nb(xs_in, ys_in, n_in, axis, value, keep_below, xs_out, ys_out):
    """CLIP AN n_in-VERTEX POLYGON AGAINST ONE AXIS-ALIGNED HALF-PLANE,
    WRITING THE RESULT INTO CALLER-PROVIDED xs_out/ys_out (LENGTH >= n_in + 1).
    RETURNS n_out, THE NUMBER OF OUTPUT VERTICES."""
    n_out = 0
    for k in range(n_in):
        px_ = xs_in[k]
        py_ = ys_in[k]
        kk = k + 1
        if kk == n_in:
            kk = 0
        qx = xs_in[kk]
        qy = ys_in[kk]

        if axis == 0:
            pval, qval = px_, qx
        else:
            pval, qval = py_, qy

        if keep_below:
            p_in = pval <= value
            q_in = qval <= value
        else:
            p_in = pval >= value
            q_in = qval >= value

        if p_in:
            xs_out[n_out] = px_
            ys_out[n_out] = py_
            n_out += 1
        if p_in != q_in:
            t = (value - pval) / (qval - pval)
            xs_out[n_out] = px_ + t * (qx - px_)
            ys_out[n_out] = py_ + t * (qy - py_)
            n_out += 1
    return n_out


@numba.njit(cache=True, inline="always")
def _clip_to_pixel_area_nb(xs, ys, ncorners, px, py, bufA_x, bufA_y, bufB_x, bufB_y):
    """CLIP THE ncorners-VERTEX POLYGON (xs, ys) TO THE UNIT SQUARE OF
    DETECTOR PIXEL (px, py) AND RETURN THE CLIPPED POLYGON'S UNSIGNED AREA.
    USES CALLER-PROVIDED SCRATCH BUFFERS (LENGTH >= ncorners + 4) PING-PONG
    STYLE ACROSS THE 4 HALF-PLANE CLIPS — NO PER-CALL ALLOCATION."""
    n = _clip_halfplane_nb(xs, ys, ncorners, 0, px - 0.5, False, bufA_x, bufA_y)
    if n < 3:
        return 0.0
    n = _clip_halfplane_nb(bufA_x, bufA_y, n, 0, px + 0.5, True, bufB_x, bufB_y)
    if n < 3:
        return 0.0
    n = _clip_halfplane_nb(bufB_x, bufB_y, n, 1, py - 0.5, False, bufA_x, bufA_y)
    if n < 3:
        return 0.0
    n = _clip_halfplane_nb(bufA_x, bufA_y, n, 1, py + 0.5, True, bufB_x, bufB_y)
    if n < 3:
        return 0.0

    # SHOELACE AREA ON THE FINAL CLIPPED POLYGON
    a = 0.0
    for k in range(n):
        kk = k + 1
        if kk == n:
            kk = 0
        a += bufB_x[k] * bufB_y[kk] - bufB_x[kk] * bufB_y[k]
    return 0.5 * abs(a)


@numba.njit(cache=True)
def _resample_weights_kernel(
    xb,
    yb,
    px_lo,
    px_hi,
    py_lo,
    py_hi,
    iOut,
    jOut,
    pxOut,
    pyOut,
    areaOut,
    coverage,
):
    """FOR EVERY (i, j) OUTPUT CELL, CLIP ITS BOUNDARY POLYGON AGAINST EVERY
    CANDIDATE DETECTOR PIXEL IN ITS [px_lo, px_hi] x [py_lo, py_hi] BOUNDING
    BOX AND RECORD THE OVERLAP AREA. CELLS WITH px_hi < px_lo ARE SKIPPED
    (OFF-DETECTOR/NON-FINITE SENTINEL, SET BY THE CALLER).

    RETURNS nOut, THE NUMBER OF (i, j, px, py, area) ENTRIES WRITTEN INTO
    THE PREALLOCATED iOut/jOut/pxOut/pyOut/areaOut[:nOut]. coverage IS
    ACCUMULATED IN PLACE."""
    n_sp, n_wl, ncorners = xb.shape
    MAXV = ncorners + 4
    bufA_x = np.empty(MAXV, dtype=np.float64)
    bufA_y = np.empty(MAXV, dtype=np.float64)
    bufB_x = np.empty(MAXV, dtype=np.float64)
    bufB_y = np.empty(MAXV, dtype=np.float64)

    nOut = 0
    for i in range(n_sp):
        for j in range(n_wl):
            plo = px_lo[i, j]
            phi = px_hi[i, j]
            if phi < plo:
                continue
            qlo = py_lo[i, j]
            qhi = py_hi[i, j]
            xs = xb[i, j]
            ys = yb[i, j]

            for py in range(qlo, qhi + 1):
                for px in range(plo, phi + 1):
                    a = _clip_to_pixel_area_nb(xs, ys, ncorners, px, py, bufA_x, bufA_y, bufB_x, bufB_y)
                    if a > 0.0:
                        iOut[nOut] = i
                        jOut[nOut] = j
                        pxOut[nOut] = px
                        pyOut[nOut] = py
                        areaOut[nOut] = a
                        nOut += 1
                        coverage[i, j] += a
    return nOut


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

    sp0, sp1 = sp_edges[:-1, None], sp_edges[1:, None]  # (n_sp, 1)
    wl0, wl1 = wl_edges[None, :-1], wl_edges[None, 1:]  # (1, n_wl)

    t = np.linspace(0.0, 1.0, edge_samples + 1)[:-1]  # exclude endpoint, (edge_samples,)

    n_sp = sp_edges.shape[0] - 1
    n_wl = wl_edges.shape[0] - 1

    shape = (n_sp, n_wl, edge_samples)

    # BROADCAST EACH EDGE OF THE CELL BOUNDARY TO SHAPE (n_sp, n_wl, edge_samples)
    bottom_sp = np.broadcast_to(sp0[:, :, None] + t * (sp1 - sp0)[:, :, None], shape)  # wl = wl0
    bottom_wl = np.broadcast_to(wl0[:, :, None], shape)

    right_sp = np.broadcast_to(sp1[:, :, None], shape)  # sp = sp1
    right_wl = np.broadcast_to(wl0[:, :, None] + t * (wl1 - wl0)[:, :, None], shape)

    top_sp = np.broadcast_to(sp1[:, :, None] + t * (sp0 - sp1)[:, :, None], shape)  # wl = wl1
    top_wl = np.broadcast_to(wl1[:, :, None], shape)

    left_sp = np.broadcast_to(sp0[:, :, None], shape)  # sp = sp0
    left_wl = np.broadcast_to(wl1[:, :, None] + t * (wl0 - wl1)[:, :, None], shape)

    spBlock = np.concatenate([bottom_sp, right_sp, top_sp, left_sp], axis=2)
    wlBlock = np.concatenate([bottom_wl, right_wl, top_wl, left_wl], axis=2)
    return spBlock, wlBlock


def _pixel_uniform_wavelength_edges(wavelengths, pixelX, pixelY, samplesPerPixel):
    """*Wavelength bin edges spaced evenly in detector pixels along an order's trace*

    **Key Arguments:**

    - ``wavelengths`` -- 1D wavelengths sampled densely along the trace, in increasing order
    - ``pixelX``, ``pixelY`` -- detector position of the trace at each wavelength
    - ``samplesPerPixel`` -- number of bins per detector pixel of path length

    **Return:**

    - ``wlEdges`` -- 1D array of wavelength edges, starting at the first usable wavelength
    """
    import numpy as np

    # KEEP THE LONGEST CONTIGUOUS RUN OF FINITE SAMPLES, SO THE PATH NEVER JUMPS ACROSS A GAP
    finite = np.isfinite(wavelengths) & np.isfinite(pixelX) & np.isfinite(pixelY)
    bounds = np.flatnonzero(np.diff(np.concatenate([[0], finite.astype(np.int8), [0]])))
    runStarts, runEnds = bounds[::2], bounds[1::2]
    if len(runStarts) == 0:
        raise ValueError("Cannot build wavelength edges: no finite trace samples on the detector.")
    longest = np.argmax(runEnds - runStarts)
    run = slice(runStarts[longest], runEnds[longest])
    wavelengths, pixelX, pixelY = wavelengths[run], pixelX[run], pixelY[run]
    if len(wavelengths) < 2:
        raise ValueError("Cannot build wavelength edges: fewer than 2 finite trace samples on the detector.")

    pathLength = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(pixelX), np.diff(pixelY)))])
    # np.interp NEEDS A STRICTLY INCREASING PATH — DROP REPEATED DETECTOR POSITIONS
    keep = np.concatenate([[True], np.diff(pathLength) > 0])
    if keep.sum() < 2:
        raise ValueError("Cannot build wavelength edges: the trace has zero length on the detector.")

    targets = np.arange(0.0, pathLength[keep][-1], 1.0 / samplesPerPixel)
    return np.interp(targets, pathLength[keep], wavelengths[keep])


def _arcsec_per_pixel(slitProbeArcsec, centreX, centreY, offsetX, offsetY):
    """*Median slit scale (arcsec per detector pixel) from trace points and points offset along the slit*

    **Key Arguments:**

    - ``slitProbeArcsec`` -- slit offset (arcsec) between each centre point and its offset point
    - ``centreX``, ``centreY`` -- detector positions of the trace centre
    - ``offsetX``, ``offsetY`` -- detector positions of the same wavelengths offset along the slit

    **Return:**

    - ``arcsecPerPixel`` -- median arcsec per detector pixel along the slit
    """
    import numpy as np

    distance = np.hypot(offsetX - centreX, offsetY - centreY)
    usable = np.isfinite(distance) & (distance > 0)
    if not usable.any():
        raise ValueError("Cannot measure the slit scale (arcsec per pixel): no slit probe landed on the detector.")
    return slitProbeArcsec / np.median(distance[usable])


def _rebin_resampling_weights(i, j, px, py, area, nSp, nWl, zoomSlit, zoomWavelength, nx, ny):
    """*Merge zoomed-grid pixel-overlap weights into the rebinned (detector-resolution) output grid*

    Matches ``image_transformer._unzoom``: sub-cells beyond the last whole zoom block are dropped, and a grid
    smaller than one zoom block is not rebinned. Overlaps of the same detector pixel within one output cell
    are summed into a single weight.

    **Key Arguments:**

    - ``i``, ``j`` -- zoomed-grid slit and wavelength indices of each overlap
    - ``px``, ``py`` -- detector pixel of each overlap
    - ``area`` -- overlap area of each entry (in detector pixels)
    - ``nSp``, ``nWl`` -- zoomed-grid shape
    - ``zoomSlit``, ``zoomWavelength`` -- zoom factors along the slit and wavelength axes
    - ``nx``, ``ny`` -- detector shape

    **Return:**

    - ``rebinned`` -- dict of ``flatIdx``, ``px``, ``py``, ``area`` arrays and the output grid ``shape``
    """
    import numpy as np

    nSpOut, nWlOut = nSp // zoomSlit, nWl // zoomWavelength
    if nSpOut == 0 or nWlOut == 0:
        zoomSlit, zoomWavelength, nSpOut, nWlOut = 1, 1, nSp, nWl

    iOut, jOut = i // zoomSlit, j // zoomWavelength
    keep = (iOut < nSpOut) & (jOut < nWlOut)
    cellIdx = iOut[keep].astype(np.int64) * nWlOut + jOut[keep]
    pixelIdx = py[keep].astype(np.int64) * nx + px[keep]

    uniqueKeys, inverse = np.unique(cellIdx * (nx * ny) + pixelIdx, return_inverse=True)
    summedArea = np.bincount(inverse, weights=area[keep], minlength=len(uniqueKeys))
    uniquePixels = uniqueKeys % (nx * ny)

    return {
        "flatIdx": uniqueKeys // (nx * ny),
        "px": uniquePixels % nx,
        "py": uniquePixels // nx,
        "area": summedArea,
        "shape": (nSpOut, nWlOut),
    }


def _apply_weights(ndarray, weights, badPixels=None):
    """*Weighted sum of detector pixels into the rebinned output grid, leaving out flagged pixels*

    Without ``badPixels`` every pixel is summed with its overlap area. With it, flagged pixels are left out of the
    sum and the cell is scaled by (all overlap area) / (good overlap area), so a cell keeps the flux of one full
    cell and does not carry a flagged pixel's value. A cell with no good pixel keeps its all-pixel sum (it is
    masked by the caller). A cell with no flagged pixel is exactly the plain weighted sum.

    **Key Arguments:**

    - ``ndarray`` -- 2D detector-space image
    - ``weights`` -- rebinned weights dict from ``_rebin_resampling_weights``
    - ``badPixels`` -- optional 2D boolean detector-space array, True where a pixel is flagged. Default *None*

    **Return:**

    - ``rectified`` -- 2D array of shape ``weights["shape"]``
    """
    import numpy as np

    nSpOut, nWlOut = weights["shape"]
    nCells = nSpOut * nWlOut
    flatIdx, area = weights["flatIdx"], weights["area"]
    weighted = ndarray[weights["py"], weights["px"]] * area
    total = np.bincount(flatIdx, weights=weighted, minlength=nCells)
    if badPixels is None:
        return total.reshape(nSpOut, nWlOut)

    isBad = badPixels[weights["py"], weights["px"]]
    goodAreaEach = np.where(isBad, 0.0, area)
    goodArea = np.bincount(flatIdx, weights=goodAreaEach, minlength=nCells)
    allArea = np.bincount(flatIdx, weights=area, minlength=nCells)
    # CELLS WITH NO FLAGGED PIXEL, OR NO GOOD PIXEL, KEEP THE PLAIN SUM
    canRenormalise = (goodArea > 0) & (goodArea < allArea)
    if not canRenormalise.any():
        return total.reshape(nSpOut, nWlOut)

    # WHERE, NOT MULTIPLICATION, SO A NON-FINITE VALUE UNDER THE MASK CANNOT POISON THE CELL
    goodSum = np.bincount(flatIdx, weights=np.where(isBad, 0.0, weighted), minlength=nCells)
    scale = np.divide(allArea, goodArea, out=np.ones_like(allArea), where=canRenormalise)
    return np.where(canRenormalise, goodSum * scale, total).reshape(nSpOut, nWlOut)
