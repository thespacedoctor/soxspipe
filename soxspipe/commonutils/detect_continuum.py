#!/usr/bin/env python
"""
*find and fit the continuum trace across all echelle orders with low-order polynomials.*

Author
: David Young & Marco Landoni

Date Created
: September 10, 2020
"""

################# GLOBAL IMPORTS ####################
from __future__ import annotations

import copy
import math
import os
import sys
from contextlib import suppress
from datetime import datetime
from typing import TYPE_CHECKING

from soxspipe.commonutils import detector_lookup, keyword_lookup
from soxspipe.commonutils.dispersion_map_to_pixel_arrays import (
    dispersion_map_to_pixel_arrays,
)
from soxspipe.commonutils.faint_order_map import make_pixels_at
from soxspipe.commonutils.faint_order_profile import frame_to_nan_array
from soxspipe.commonutils.faint_order_rows import clip_measured_residuals, measured_residuals, transferred_row_mask
from soxspipe.commonutils.faint_order_trace import (
    faint_trace_settings,
    join_slit_mismatches,
    transfer_faint_orders,
)
from soxspipe.commonutils.figure_guard import closes_its_figures
from soxspipe.commonutils.filenamer import filenamer
from soxspipe.commonutils.polynomials import (
    chebyshev_order_xy_polynomials,
    chebyshev_xy_polynomial,
    solve_linear_poly_coefficients,
)
from soxspipe.commonutils.toolkit import (
    QC_PLOT_DECORATION_HEIGHT_INCHES,
    QC_RESIDUAL_ROW_HEIGHT_INCHES,
    QC_STRIP_ROW_HEIGHT_INCHES,
    append_qc,
    cut_image_slice,
    get_calibration_lamp,
    has_sequence_suffixed_traces,
    legend_below_axis,
    nodding_sequence_suffix,
    order_label,
    order_traces_table,
    qc_image_panel_height,
    qc_table_heights,
    read_spectral_format,
)

if TYPE_CHECKING:
    from collections.abc import Sequence

    import numpy as np
    import pandas as pd

os.environ["TERM"] = "vt100"

ORDER_TRACES_EXTNAME = "ORDER_TRACES"

# FALLBACK FOR SETTINGS FILES WRITTEN BEFORE `max-unsampled-trace-fraction` EXISTED; MATCHES THE DEFAULT SETTINGS
DEFAULT_MAX_UNSAMPLED_TRACE_FRACTION = 0.2


def longest_unsampled_stretch(
    sampledWavelengths: Sequence[float] | np.ndarray | pd.Series,
    acceptedWavelengths: Sequence[float] | np.ndarray | pd.Series,
) -> tuple[float, float, float]:
    """*find the longest stretch of an order that no accepted trace sample covers*

    The stretch runs between two neighbouring accepted samples, or from an end of the order to the nearest accepted
    sample. The fitted trace is interpolated or extrapolated across it.

    **Key Arguments:**

    - ``sampledWavelengths`` -- the wavelengths (nm) of every sample position along the order
    - ``acceptedWavelengths`` -- the wavelengths (nm) of the samples kept after the trace fit

    **Return:**

    - ``fraction`` -- the length of the stretch as a fraction of the sampled length of the order
    - ``startNm`` -- the wavelength (nm) at the blue end of the stretch
    - ``endNm`` -- the wavelength (nm) at the red end of the stretch

    **Usage:**

    ```python
    fraction, startNm, endNm = longest_unsampled_stretch(sampled, accepted)
    ```
    """
    import numpy as np

    sampled = np.asarray(sampledWavelengths, dtype=float)
    sampled = sampled[np.isfinite(sampled)]
    accepted = np.asarray(acceptedWavelengths, dtype=float)
    accepted = accepted[np.isfinite(accepted)]

    lowNm, highNm = sampled.min(), sampled.max()
    edges = np.concatenate([[lowNm], np.sort(accepted), [highNm]])
    gaps = np.diff(edges)
    longest = int(np.argmax(gaps))
    orderLengthNm = highNm - lowNm
    fraction = float(gaps[longest] / orderLengthNm) if orderLengthNm > 0 else 0.0

    return fraction, float(edges[longest]), float(edges[longest + 1])


def _order_colours(orders):
    """*assign each order a colour from the matplotlib property cycle, keyed by the order itself*

    **Key Arguments:**

    - ``orders`` -- the orders to colour, in plotting order

    **Return:**

    - ``orderColours`` -- dictionary of order to colour, wrapping around once the cycle is exhausted
    """
    import matplotlib.pyplot as plt

    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    return {order: cycle[i % len(cycle)] for i, order in enumerate(orders)}


class _base_detect:
    # FILE-NAME SUFFIX OF THE NOT-FLATTENED REPEAT OF A TRACE. ONLY DETECT_CONTINUUM SETS IT
    notFlattened = ""

    def fit_order_polynomial(self, pixelList, order, axisBDeg, axisACol, axisBCol, exponentsIncluded=False):
        """*iteratively fit the dispersion map polynomials to the data, clipping residuals with each iteration*

        **Key Arguments:**

        - ``pixelList`` -- data-frame group containing x,y pixel array
        - ``order`` -- the order to fit
        - ``axisBDeg`` -- degree for polynomial to fit
        - ``axisACol`` -- name of columns containing axis to be fitted
        - ``axisBCol`` -- name of columns containing free axis (values known)
        - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**

        - ``coeffs`` -- the coefficients of the polynomial fit
        - ``pixelList`` -- the pixel list but now with fits and residuals included
        """
        self.log.debug("starting the ``fit_order_polynomial`` method")

        import numpy as np
        from astropy.stats import sigma_clip
        from scipy.optimize import curve_fit

        arm = self.arm
        self.axisBDeg = axisBDeg

        clippedCount = 1

        poly = chebyshev_xy_polynomial(
            log=self.log,
            axisBCol=axisBCol,
            axisBDeg=self.axisBDeg,
            exponentsIncluded=exponentsIncluded,
        ).poly

        clippingSigma = self.recipeSettings["poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings["poly-clipping-iteration-limit"]

        iteration = 0
        mask = pixelList["order"] == order
        pixelListFiltered = pixelList.loc[mask]

        coeff = np.ones(self.axisBDeg + 1)
        while clippedCount > 0 and iteration < clippingIterationLimit:
            pixelListFiltered = pixelList.loc[mask]

            startCount = len(pixelListFiltered.index)
            iteration += 1
            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY
            # NOTE X AND Y COLUMN ARE CORRECLY IN xdata AND ydata - WANT TO
            # FIND X (UNKNOWN) WRT Y (KNOWNN)
            try:
                coeff, pcov_x = curve_fit(
                    poly,
                    xdata=pixelListFiltered,
                    ydata=pixelListFiltered[axisACol].values,
                    p0=coeff,
                    maxfev=30000,
                )
            except TypeError:
                # REMOVE THIS ORDER FROM PIXEL LIST
                pixelList = pixelList.loc[~mask]
                coeff = None
                return coeff, pixelList
            except Exception as e:
                raise e

            res, res_mean, res_std, res_median, xfit, poly = self.calculate_residuals(
                orderPixelTable=pixelListFiltered,
                coeff=coeff,
                axisACol=axisACol,
                axisBCol=axisBCol,
            )

            pixelList.loc[mask, "x_fit_res"] = res
            pixelList.loc[mask, "x_fit"] = xfit

            # SIGMA-CLIP THE DATA
            masked_residuals = sigma_clip(
                res,
                sigma_lower=clippingSigma,
                sigma_upper=clippingSigma,
                maxiters=1,
                cenfunc="median",
                stdfunc="mad_std",
            )
            pixelList.loc[mask, "mask"] = masked_residuals.mask

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = pixelList["mask"] == True
            pixelList = pixelList.loc[~removeMask]
            pixelListFiltered = pixelList.loc[mask]
            clippedCount = startCount - len(pixelListFiltered.index)

            sys.stdout.write("\x1b[1A\x1b[2K")
            sys.stdout.flush()
            self.log.print(
                f"\t\tORDER {order:0.0f}: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial"
            )

        self.log.debug("completed the ``fit_order_polynomial`` method")
        return coeff, pixelList

    def fit_global_polynomial(
        self,
        pixelList,
        axisACol="cont_x",
        axisBCol="cont_y",
        orderCol="order",
        exponentsIncluded=False,
        writeQCs=False,
    ):
        """*iteratively fit the global polynomial to the data, fitting axisA as a function of axisB, clipping residuals with each iteration*

        **Key Arguments:**

        - ``pixelList`` -- data-frame group containing x,y pixel array
        - ``exponentsIncluded`` -- the exponents have already been calculated in the dataframe so no need to regenerate. Default *False*

        **Return:**

        - ``coeffs`` -- the coefficients of the polynomial fit
        - ``pixelList`` -- the pixel list but now with fits and residuals included
        - ``allClipped`` -- data that was sigma-clipped
        """
        self.log.debug("starting the ``fit_global_polynomial`` method")

        import numpy as np
        import pandas as pd

        arm = self.arm

        # WORK ON A PRIVATE COPY: CALLERS PASS FILTERED SLICES (E.G. dropna) AND USE THE RETURNED TABLE
        pixelList = pixelList.copy()

        clippedCount = 1

        poly = chebyshev_order_xy_polynomials(
            log=self.log,
            axisBCol=axisBCol,
            orderCol=orderCol,
            orderDeg=self.orderDeg,
            axisBDeg=self.axisBDeg,
            axisB=self.axisB,
            exponentsIncluded=exponentsIncluded,
        ).poly

        clippingSigmaHigh = self.recipeSettings["poly-fitting-residual-clipping-sigma"]
        clippingIterationLimit = self.recipeSettings["poly-clipping-iteration-limit"]

        if axisACol == "stddev":
            clippingSigmaLow = clippingSigmaHigh
            clippingSigmaHigh = clippingSigmaLow * 2
        else:
            clippingSigmaLow = clippingSigmaHigh

        iteration = 0

        allClipped = []

        if "cont_x" in pixelList.columns:
            mask = pixelList["cont_x"].isna() | pixelList["cont_y"].isna()
            pixelList.loc[mask, "pre-clipped"] = True

        # REMOVE FILTERED ROWS FROM DATA FRAME
        if "pre-clipped" in pixelList.columns:
            removeMask = pixelList["pre-clipped"] == True
            allClipped.append(pixelList.loc[removeMask])
            pixelList = pixelList.loc[~removeMask].copy()

        nCoeff = (self.axisBDeg + 1) * (self.orderDeg + 1)

        # PLACEHOLDER SO THE FINAL RESIDUALS CAN BE EVALUATED IF THE CLIPPING LOOP NEVER RUNS
        coeff = np.ones(nCoeff)
        while clippedCount > 0 and iteration < clippingIterationLimit:
            startCount = len(pixelList.index)
            iteration += 1

            # USE LEAST-SQUARED CURVE FIT TO FIT CHEBY POLY

            # TOO FEW POINTS TO CONSTRAIN THE POLYNOMIAL (OR NONE LEFT) - NO SOLUTION
            if len(pixelList.index) < nCoeff:
                coeff = None
                return coeff, pixelList, pixelList

            coeff = solve_linear_poly_coefficients(
                poly=poly,
                pixelList=pixelList,
                yValues=pixelList[axisACol].values,
                nCoeff=nCoeff,
            )

            res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
                orderPixelTable=pixelList,
                coeff=coeff,
                orderCol=orderCol,
                axisACol=axisACol,
                axisBCol=axisBCol,
            )

            pixelList[f"{axisACol}_fit_res"] = res
            pixelList[f"{axisACol}_fit"] = xfit

            # SIGMA-CLIP THE DATA, LEAVING OUT THE SYNTHETIC ROWS OF A TRANSFERRED FAINT TRACE
            pixelList["mask"] = clip_measured_residuals(
                res, transferred_row_mask(pixelList), clippingSigmaLow, clippingSigmaHigh
            )

            # REMOVE FILTERED ROWS FROM DATA FRAME
            removeMask = pixelList["mask"] == True
            allClipped.append(pixelList.loc[removeMask])
            pixelList = pixelList.loc[~removeMask].copy()
            clippedCount = startCount - len(pixelList.index)

            if iteration > 1:
                sys.stdout.flush()
                sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.print(
                f"\t\tGLOBAL FIT: {clippedCount} pixel positions where clipped in iteration {iteration} of fitting the polynomial"
            )

        if len(allClipped):
            allClipped = pd.concat(allClipped, ignore_index=True)

        # A CAP EXIT WITH ROWS STILL BEING CLIPPED LEAVES COEFFICIENTS THAT INCLUDE THE ROWS JUST REJECTED
        if iteration > 0 and clippedCount > 0:
            self.log.info(
                f"fit_global_polynomial: iteration limit of {clippingIterationLimit} reached with {clippedCount} rows "
                "still being clipped; refitting to the surviving rows"
            )
            if len(pixelList.index) < nCoeff:
                return None, pixelList, pixelList
            coeff = solve_linear_poly_coefficients(
                poly=poly,
                pixelList=pixelList,
                yValues=pixelList[axisACol].values,
                nCoeff=nCoeff,
            )

        res, res_mean, res_std, res_median, xfit = self.calculate_residuals(
            orderPixelTable=pixelList,
            coeff=coeff,
            orderCol=orderCol,
            axisACol=axisACol,
            axisBCol=axisBCol,
            writeQCs=writeQCs,
        )
        pixelList[f"{axisACol}_fit_res"] = res
        pixelList[f"{axisACol}_fit"] = xfit

        self.log.debug("completed the ``fit_global_polynomial`` method")
        return coeff, pixelList, allClipped

    def calculate_residuals(self, orderPixelTable, coeff, axisACol, axisBCol, orderCol=False, writeQCs=False):
        """*calculate residuals of the polynomial fits against the observed line postions*

        **Key Arguments:**

        - ``orderPixelTable`` -- data-frame containing pixel list for given order
        - ``coeff`` -- the coefficients of the fitted polynomial
        - ``axisACol`` -- name of x-pixel column
        - ``axisBCol`` -- name of y-pixel column
        - ``orderCol`` -- name of the order column (global fits only)
        - ``writeQCs`` -- write the QCs to dataframe? Default *False*

        **Return:**

        - ``res`` -- x residuals
        - ``mean`` -- the mean of the residuals
        - ``std`` -- the stdev of the residuals
        - ``median`` -- the median of the residuals
        - ``xfit`` -- fitted x values
        """
        self.log.debug("starting the ``calculate_residuals`` method")

        import numpy as np
        import pandas as pd

        arm = self.arm

        poly = chebyshev_order_xy_polynomials(
            log=self.log,
            axisBCol=axisBCol,
            orderCol=orderCol,
            orderDeg=self.orderDeg,
            axisBDeg=self.axisBDeg,
        ).poly

        # CALCULATE RESIDUALS BETWEEN GAUSSIAN PEAK LINE POSITIONS AND POLY
        # FITTED POSITIONS
        thisFit = poly(orderPixelTable, *coeff)
        res = thisFit - orderPixelTable[axisACol].values

        # GET UNIQUE VALUES IN COLUMN
        uniqueorders = len(orderPixelTable["order"].unique())

        # CALCULATE COMBINED RESIDUALS AND STATS, OVER THE MEASURED ROWS ONLY WHEN A FAINT TRACE WAS TRANSFERRED
        isTransferred = transferred_row_mask(orderPixelTable)
        statRes = measured_residuals(res, isTransferred)
        res_mean = np.ma.mean(np.ma.abs(statRes))
        res_std = np.ma.std(statRes)
        res_median = np.ma.median(np.ma.abs(statRes))

        if writeQCs:
            utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

            tag = "continuum"
            if "order-centres" in self.recipeName.lower():
                tag = "order centres"

            self.qc = pd.concat(
                [
                    self.qc,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "qc_name": f"{self.axisA.upper()} RES MIN",
                                "qc_value": round(float(statRes.min()), 3),
                                "qc_comment": f"[px] Minimum residual in {tag} fit along {self.axisA}-axis",
                                "qc_unit": "px",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "to_header": True,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "qc_name": f"{self.axisA.upper()} RES MAX",
                                "qc_value": round(float(statRes.max()), 3),
                                "qc_comment": f"[px] Maximum residual in {tag} fit along {self.axisA}-axis",
                                "qc_unit": "px",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "to_header": True,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )
            self.qc = pd.concat(
                [
                    self.qc,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "qc_name": f"{self.axisA.upper()} RES SD",
                                "qc_value": round(float(res_std), 3),
                                "qc_comment": f"[px] Std-dev of residual {tag} fit along {self.axisA}-axis",
                                "qc_unit": "px",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "to_header": True,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

            self.qc = pd.concat(
                [
                    self.qc,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "qc_name": f"{self.axisA.upper()} RES MEDIAN",
                                "qc_value": round(float(res_mean), 3),
                                "qc_comment": f"[px] Median abolute residual {tag} fit along {self.axisA}-axis",
                                "qc_unit": "px",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "to_header": True,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

            if "order" in self.recipeName.lower():
                c = "Number of order centre traces found"
            else:
                c = "Number of orders containing an object trace"
                if isTransferred is not None and isTransferred.any():
                    # THE COUNT INCLUDES A TRANSFERRED FAINT ORDER ON PURPOSE
                    c += ", a transferred faint order included"

            self.qc = pd.concat(
                [
                    self.qc,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "qc_name": "N ORDERS",
                                "qc_value": uniqueorders,
                                "qc_comment": c,
                                "qc_unit": None,
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "to_header": True,
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

        self.log.debug("completed the ``calculate_residuals`` method")
        return res, res_mean, res_std, res_median, thisFit

    def write_order_table_to_file(self, frame, orderPolyTable, orderMetaTable):
        """*write out the fitted polynomial solution coefficients to file*

        **Key Arguments:**

        - ``frame`` -- the calibration frame used to generate order location data
        - ``orderPolyTable`` -- data-frames containing centre location coefficients (and possibly also order edge coeffs)
        - ``orderMetaTable`` -- extra order meta data to be added in an extra FITS extension

        **Return:**

        - ``order_table_path`` -- path to the order table file
        """

        self.log.debug("starting the ``write_order_table_to_file`` method")

        arm = self.arm
        kw = self.kw

        if not self.sofName:
            filename = filenamer(log=self.log, frame=frame, settings=self.settings)

        else:
            filename = self.sofName + ".fits"
        filename = filename.replace("MFLAT", "FLAT")

        # if self.inst.upper() == "SOXS":
        #     filename = filename.replace("_DLAMP", "")
        #     filename = filename.replace("_QLAMP", "")

        if "mflat" in self.recipeName.lower():
            filename = filename.upper().replace("FLAT", "OLOC").replace(".FITS", ".fits")

        elif "stare" in self.recipeName.lower():
            filename = filename.upper().split(".FITS")[0] + "_OBJTRACE" + self.notFlattened + ".fits"
        elif has_sequence_suffixed_traces(self.recipeName):
            # sequence = "A" if int(frame.header['HIERARCH ESO SEQ CUMOFF Y'] > 0) else "B"
            filename = (
                filename.upper().split(".FITS")[0] + "_OBJTRACE" + self.noddingSequence + self.notFlattened + ".fits"
            )

        if self.lampTag and self.inst.upper() != "SOXS":
            filename = filename.replace(".fits", f"{self.lampTag}.fits")

        order_table_path = f"{self.productDir}/{filename}"

        header = copy.deepcopy(frame.header)
        with suppress(KeyError):
            header.pop(kw("DPR_TECH"))
        with suppress(KeyError):
            header.pop(kw("DPR_CATG"))
        with suppress(KeyError):
            header.pop(kw("DPR_TYPE"))
        with suppress(KeyError):
            header.pop(kw("DET_READ_SPEED"))
        with suppress(KeyError):
            header.pop(kw("CONAD"))
        with suppress(KeyError):
            header.pop(kw("GAIN"))
        with suppress(KeyError):
            header.pop(kw("RON"))

        header["HIERARCH " + kw("PRO_TECH")] = "ECHELLE,SLIT"

        header[kw("SEQ_ARM")] = arm
        header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
        tables = [orderPolyTable, orderMetaTable]
        if (
            "stare" not in self.recipeName.lower()
            and "nod" not in self.recipeName.lower()
            and "offset" not in self.recipeName.lower()
        ):
            header["HIERARCH " + kw("PRO_CATG")] = f"ORDER_TAB_{arm}".upper()
            tables.append(self.build_order_traces_extension(frame, orderPolyTable, orderMetaTable))
        else:
            header["HIERARCH " + kw("PRO_CATG")] = f"OBJECT_TAB_{arm}".upper()

        from soxspipe.commonutils.phase3 import write_fits_table_to_disk

        write_fits_table_to_disk(
            log=self.log,
            settings=self.settings,
            header=header,
            tables=tables,
            filePath=order_table_path,
            qc=self.qc,
        )

        self.log.debug("completed the ``write_order_table_to_file`` method")
        return order_table_path

    def build_order_traces_extension(self, frame, orderPolyTable, orderMetaTable):
        """*sample the order traces at every unbinned pixel, as the ``ORDER_TRACES`` table extension*

        **Key Arguments:**

        - ``frame`` -- the calibration frame used to generate order location data
        - ``orderPolyTable`` -- data-frame containing centre location coefficients (and possibly also order edge coeffs)
        - ``orderMetaTable`` -- order meta data, with the fitted range of each order

        **Return:**

        - ``orderTraces`` -- astropy table named ``ORDER_TRACES``, one row per order and sampled pixel
        """
        self.log.debug("starting the ``build_order_traces_extension`` method")

        from astropy.table import Table

        tableBinx, tableBiny, limitsBinx, limitsBiny = self._order_table_binning(frame)
        orderTraces = Table.from_pandas(
            order_traces_table(
                log=self.log,
                orderPolyTable=orderPolyTable,
                orderMetaTable=orderMetaTable,
                tableBinx=tableBinx,
                tableBiny=tableBiny,
                limitsBinx=limitsBinx,
                limitsBiny=limitsBiny,
            )
        )
        orderTraces.meta["EXTNAME"] = ORDER_TRACES_EXTNAME

        self.log.debug("completed the ``build_order_traces_extension`` method")
        return orderTraces

    def _order_table_binning(self, frame):
        """*find the absolute binning of the pixel scale of the table's coefficients and of its metadata limits*

        An order-centre table is measured on its frame, so ``self.binx`` and ``self.biny`` (the frame's binning)
        apply to both. The master-flat detector gets only the ratio of the flat's binning to the order-centre
        table's, so the order-centre table's binning is the flat's binning (from its header) over that ratio. The
        flat's own binning applies to the metadata limits, which are the flat's frame shape.

        **Key Arguments:**

        - ``frame`` -- the calibration frame used to generate order location data

        **Return:**

        - ``tableBinx``, ``tableBiny`` -- the binning of the pixel scale the coefficients are in
        - ``limitsBinx``, ``limitsBiny`` -- the binning of the frame the metadata limits were measured on
        """
        if "mflat" not in self.recipeName.lower():
            return self.binx, self.biny, self.binx, self.biny

        try:
            limitsBinx = int(frame.header[self.kw("WIN_BINX")])
            limitsBiny = int(frame.header[self.kw("WIN_BINY")])
        except KeyError as e:
            self.log.debug(f"_order_table_binning: no binning in the flat header, assuming none: {e}")
            limitsBinx = limitsBiny = 1

        return round(limitsBinx / self.binx), round(limitsBiny / self.biny), limitsBinx, limitsBiny


class detect_continuum(_base_detect):
    """
    *find and fit the continuum trace across all echelle orders with low-order polynomials.*

    **Key Arguments:**

    - ``log`` -- logger
    - ``traceFrame`` -- calibrated frame containing a source trace (CCDObject)
    - ``dispersion_map`` -- path to dispersion map file containing polynomial fits of the dispersion solution for the frame
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``recipeName`` -- the recipe name as given in the settings dictionary
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products
    - ``sofName`` ---- name of the originating SOF file
    - ``binx`` -- binning in x-axis
    - ``biny`` -- binning in y-axis
    - ``lampTag`` -- add this tag to the end of the product filename. Default *False*
    - ``locationSetIndex`` -- the index of the AB cycle locations, or the offset stack's own index (nodding and
      offset modes only). Default *False*
    - ``orderPixelTable`` -- this is used for tuning the pipeline.  Default *False*
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
    - ``debug`` -- if *True* then extra debugging information is printed. Default *False*
    - ``notFlattened`` -- if *True* the trace is the repeat that a stare, nod or offset standard makes for the
      response curve: its files get a ``_NOTFLAT`` suffix and it adds no product row. Default *False*
    - ``slitOffsetsArcsec`` -- ``{order: offset}`` slit-illumination offsets in arcsec, used to anchor a faint VIS
      order on its bright neighbour. Default *None*

    **Usage:**

    To use the ``detect_continuum`` object, use the following:

    ```python
    from soxspipe.commonutils import detect_continuum
    detector = detect_continuum(
        log=log,
        traceFrame=traceFrame,
        dispersion_map=dispersion_map,
        settings=settings,
        recipeName="soxs-order-centres"
    )
    order_table_path = detector.get()
    ```
    """

    # THE FRAME-TO-MAP BINNING, SET WHEN THE TRACE IS SAMPLED, AND THE SLIT-ILLUMINATION OFFSET OF EACH ORDER (IF KNOWN)
    mapBinning = None
    slitOffsetsArcsec = None
    faintTransferReports = ()

    def __init__(
        self,
        log,
        traceFrame,
        dispersion_map,
        settings=False,
        recipeSettings=False,
        recipeName=False,
        qcTable=False,
        productsTable=False,
        sofName=False,
        binx=1,
        biny=1,
        lampTag=False,
        locationSetIndex=False,
        orderPixelTable=False,
        startNightDate="",
        debug=False,
        notFlattened=False,
        slitOffsetsArcsec=None,
    ):
        self.log = log
        log.debug("instantiating a new 'detect_continuum' object")

        import copy

        self.settings = settings
        self.noddingSequence = nodding_sequence_suffix(
            header=traceFrame.header, locationSetIndex=locationSetIndex, recipeName=recipeName
        )

        self.notFlattened = "_NOTFLAT" if notFlattened else ""
        self.recipeName = recipeName
        self.traceFrame = traceFrame
        self.dispersion_map = dispersion_map
        self.qc = qcTable
        self.products = productsTable
        self.sofName = sofName
        self.binx = binx
        self.biny = biny
        self.recipeSettings = copy.deepcopy(recipeSettings["detect-continuum"])
        self.lampTag = lampTag
        self.orderPixelTable = orderPixelTable
        self.startNightDate = startNightDate
        self.debug = debug
        self.slitOffsetsArcsec = slitOffsetsArcsec

        # KEYWORD LOOKUP OBJECT - LOOKUP KEYWORD FROM DICTIONARY IN RESOURCES
        # FOLDER
        self.kw = keyword_lookup(log=self.log, settings=self.settings).get
        self.arm = traceFrame.header[self.kw("SEQ_ARM")]
        self.dateObs = traceFrame.header[self.kw("DATE_OBS")]
        self.inst = traceFrame.header[self.kw("INSTRUME")]
        self.exptime = traceFrame.header[self.kw("EXPTIME")]

        # if self.exptime < 59 and recipeName != "soxs-stare":
        #     raise Exception("too low")

        # DETECTOR PARAMETERS LOOKUP OBJECT
        self.detectorParams = detector_lookup(log=log, settings=settings).get(self.arm)

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            sliceAxis = "x"
            sliceAntiAxis = "y"
        else:
            sliceAxis = "y"
            sliceAntiAxis = "x"
        self.sliceAxis = sliceAxis
        self.sliceAntiAxis = sliceAntiAxis

        # DEG OF THE POLYNOMIALS TO FIT THE ORDER CENTRE LOCATIONS
        self.axisBDeg = self.recipeSettings["disp-axis-deg"]
        self.orderDeg = self.recipeSettings["order-deg"]

        self.lamp = get_calibration_lamp(log=log, frame=traceFrame, kw=self.kw)

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=recipeName,
            startNightDate=startNightDate,
        )

        # SET IMAGE ORIENTATION
        if self.detectorParams["dispersion-axis"] == "x":
            self.axisA = "x"
            self.axisB = "y"
            self.coeff_dict = {
                "degorder_cent": self.orderDeg,
                "degy_cent": self.axisBDeg,
            }
        else:
            self.axisA = "y"
            self.axisB = "x"
            self.coeff_dict = {
                "degorder_cent": self.orderDeg,
                "degx_cent": self.axisBDeg,
            }

        return

    def get(self):
        """
        *return the order centre table filepath*

        **Return:**

        - ``order_table_path`` -- file path to the order centre table giving polynomial coeffs to each order fit
        """
        self.log.debug("starting the ``get`` method")

        from datetime import datetime

        import numpy as np
        import pandas as pd

        arm = self.arm
        coeff_dict = self.coeff_dict

        if isinstance(self.orderPixelTable, bool):
            orderPixelTable, detectionPercentage = self.sample_trace()
            if detectionPercentage < 10:
                self.log.print(
                    "Could not converge on a good fit to the continuum. Please check the quality of your data."
                )
                return None, self.qc, self.products, None, None, None
            faintSettings = self._faint_trace_settings()
            orderPixelTable = self._transfer_faint_order_traces(orderPixelTable, faintSettings)
        else:
            faintSettings = None
            orderPixelTable = self.orderPixelTable

        # THE SYNTHETIC ROWS OF A TRANSFERRED FAINT TRACE ARE NOT DETECTED SAMPLES
        foundLines = len(orderPixelTable.index)
        isTransferred = transferred_row_mask(orderPixelTable)
        if isTransferred is not None:
            foundLines -= int(isTransferred.sum())

        self.log.print("\n\t## FINDING GLOBAL POLYNOMIAL SOLUTION FOR CONTINUUM TRACES\n")

        # GET UNIQUE VALUES IN COLUMN
        uniqueOrders = orderPixelTable["order"].unique()

        orderLocations = {}
        orderPixelTable[f"cont_{self.axisA}_fit"] = np.nan
        orderPixelTable[f"cont_{self.axisA}_fit_res"] = np.nan

        # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
        fitFound = False
        tryCount = 0
        backupOrderPixelTable = orderPixelTable.copy()
        while not fitFound and tryCount < 5:
            # SETUP EXPONENTS AHEAD OF TIME - SAVES TIME ON POLY FITTING
            for i in range(0, self.axisBDeg + 1):
                orderPixelTable[f"{self.axisB}_pow_{i}"] = orderPixelTable[f"cont_{self.axisB}"].pow(i)
            for i in range(0, self.orderDeg + 1):
                orderPixelTable[f"order_pow_{i}"] = orderPixelTable["order"].pow(i)
            try:
                coeff, orderPixelTable, clippedDataCentre = self.fit_global_polynomial(
                    pixelList=orderPixelTable,
                    axisACol=f"cont_{self.axisA}",
                    axisBCol=f"cont_{self.axisB}",
                    exponentsIncluded=True,
                    writeQCs=True,
                )
                # THE GATES JUDGE THE MEASURED ROWS: SYNTHETIC ROWS LIE ON THE TRACE AND WOULD HIDE A BAD FIT
                mean_res = np.mean(
                    np.abs(
                        measured_residuals(
                            orderPixelTable[f"cont_{self.axisA}_fit_res"].values,
                            transferred_row_mask(orderPixelTable),
                        )
                    )
                )

                if "order" in self.recipeName.lower() and mean_res > 2:
                    orderPixelTable = backupOrderPixelTable
                    raise AttributeError("Failed to continuum trace")
                if mean_res > 10 and not (self.inst == "SOXS" and self.arm == "VIS"):
                    # BAD FIT ... FORCE A FAIL
                    orderPixelTable = backupOrderPixelTable
                    raise AttributeError("Failed to continuum trace")

                n_coeff = 0
                for i in range(0, self.orderDeg + 1):
                    for j in range(0, self.axisBDeg + 1):
                        coeff_dict[f"cent_{i}{j}"] = coeff[n_coeff]
                        n_coeff += 1

                # ITERATIVELY FIT THE POLYNOMIAL SOLUTIONS TO THE DATA
                coeff, orderPixelTable, clippedData = self.fit_global_polynomial(
                    pixelList=orderPixelTable,
                    axisACol="gauss_stddev",
                    axisBCol=f"cont_{self.axisB}",
                    exponentsIncluded=True,
                )

                if len(clippedData) and len(clippedData.index):
                    clippedData = pd.concat([clippedDataCentre, clippedData], ignore_index=True)
                else:
                    clippedData = clippedDataCentre
                fitFound = True

            except Exception as e:
                degList = [self.axisBDeg, self.orderDeg]
                degList[degList.index(max(degList))] -= 1
                self.axisBDeg, self.orderDeg = degList
                coeff_dict["degorder_cent"] = self.orderDeg
                coeff_dict[f"deg{self.axisB}_cent"] = self.axisBDeg
                self.log.print(
                    f"{self.axisB} and Order fitting orders reduced to {self.axisBDeg}, {self.orderDeg} to try and successfully fit the continuum."
                )
                self.recipeSettings["disp-axis-deg"] = self.axisBDeg
                self.recipeSettings["order-deg"] = self.orderDeg
                tryCount += 1
                if tryCount == 5:
                    self.log.print(
                        "Could not converge on a good fit to the continuum. Please check the quality of your data or adjust your fitting parameters."
                    )
                    return None, self.qc, self.products, None, None, None
                if self.settings["tune-pipeline"]:
                    raise e

        if "wavelength" in backupOrderPixelTable.columns:
            self.check_trace_extrapolation(backupOrderPixelTable, orderPixelTable)

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        try:
            nclip = len(clippedData.index)
        except (AttributeError, TypeError) as e:
            self.log.debug(f"get: `nclip = len(clippedData.index)` failed, continuing: {e}")
            nclip = 0
        self.qc = pd.concat(
            [
                self.qc,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "SAMPLES CLIP NUM",
                            "qc_value": nclip,
                            "qc_comment": "Number of continuum sample clipped during solution fitting",
                            "qc_unit": None,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        pclip = nclip / foundLines if foundLines > 0 else 0
        self.qc = pd.concat(
            [
                self.qc,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "SAMPLES CLIP FRAC",
                            "qc_value": round(float(pclip), 3),
                            "qc_comment": "Fraction of detected continuum samples clipped during solution fitting",
                            "qc_unit": None,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        # orderLocations[o] = coeff
        coeff_dict["degorder_std"] = self.orderDeg
        coeff_dict[f"deg{self.axisB}_std"] = self.axisBDeg
        n_coeff = 0
        for i in range(0, self.orderDeg + 1):
            for j in range(0, self.axisBDeg + 1):
                coeff_dict[f"std_{i}{j}"] = coeff[n_coeff]
                n_coeff += 1
        coeffColumns = coeff_dict.keys()

        orderPolyTable = pd.DataFrame([coeff_dict])
        self._warn_on_join_slit_jumps(orderPixelTable, orderPolyTable, faintSettings)

        # HERE IS THE LINE LIST IF NEEDED FOR QC
        orderPixelTable.drop(columns=["mask"], inplace=True)

        # OPTIMISE SUPER: 10s 8 hits
        orderMetaTable = self.build_order_metadata(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable,
        )
        plotPath = self.plot_results(
            orderPixelTable=orderPixelTable,
            orderPolyTable=orderPolyTable,
            clippedData=clippedData,
        )

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        basename = os.path.basename(plotPath)

        if "order" in self.recipeName.lower():
            label = "ORDER_CENTRES_RES"
            product_desc = "Residuals of the order centre polynomial fit"
        else:
            label = "OBJECT_TRACE_RES"
            product_desc = "Residuals of the object trace polynomial fit"

        # THE NOT-FLATTENED REPEAT FINDS THE SAME TRACE, SO IT ADDS NO SECOND ROW FOR THE PDF
        if not isinstance(self.products, bool) and not self.notFlattened:
            self.products = pd.concat(
                [
                    self.products,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "product_label": label + self.noddingSequence,
                                "file_name": basename,
                                "file_type": "PDF",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "product_desc": product_desc,
                                "file_path": plotPath,
                                "label": "QC",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

        # WRITE OUT THE FITS TO THE ORDER CENTRE TABLE
        order_table_path = self.write_order_table_to_file(
            frame=self.traceFrame,
            orderPolyTable=orderPolyTable,
            orderMetaTable=orderMetaTable,
        )

        # mean_res = np.mean(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))
        # std_res = np.std(np.abs(orderPixelTable[f'cont_{self.axisA}_fit_res'].values))

        self.log.debug("completed the ``get`` method")
        return (
            order_table_path,
            self.qc,
            self.products,
            orderPolyTable,
            orderPixelTable,
            orderMetaTable,
        )

    def _faint_trace_settings(self):
        """*the faint-order trace settings, or None where the trace is off*

        The trace is on for the VIS arm of SOXS, in a recipe that is not an order-centre one, when the
        ``faint-order-trace`` setting is on.
        """
        settings = faint_trace_settings.from_block(self.recipeSettings)
        isVisTrace = self.inst.upper() == "SOXS" and self.arm.upper() == "VIS"
        isOn = settings.isEnabled and isVisTrace and "order" not in self.recipeName.lower()
        return settings if isOn else None

    def _make_pixels_at(self):
        """*the function that turns an order, wavelengths and a slit position into pixels of the trace frame*"""
        return make_pixels_at(
            self.log, self.dispersion_map, self.mapBinning or {"x": 1.0, "y": 1.0}, self.axisA, self.axisB
        )

    def _transfer_faint_order_traces(self, orderPixelTable, settings):
        """*replace the sampled trace of each faint order by one anchored to its bright neighbour*

        When the trace is off (``settings`` is None), or no faint order has a bright partner, the table comes back as
        the same object.

        **Key Arguments:**

        - ``orderPixelTable`` -- the table of sampled trace positions
        - ``settings`` -- the ``faint_trace_settings``, or None where the trace is off

        **Return:**

        - ``orderPixelTable`` -- the table to fit
        """
        if settings is None:
            return orderPixelTable
        orderPixelTable, reports = transfer_faint_orders(
            orderPixelTable,
            frame_to_nan_array(self.traceFrame),
            self._make_pixels_at(),
            slitOffsets=self.slitOffsetsArcsec,
            settings=settings,
            axisA=self.axisA,
            axisB=self.axisB,
        )
        self._record_faint_trace_qc(reports)
        return orderPixelTable

    def _fitted_trace(self, orderPolyTable):
        """*the function that evaluates the fitted order-centre polynomial: ``trace(order, axisBPixels)``*"""
        import numpy as np
        import pandas as pd

        poly = chebyshev_order_xy_polynomials(
            log=self.log,
            axisBCol=self.axisB,
            orderCol="order",
            orderDeg=self.orderDeg,
            axisBDeg=self.axisBDeg,
        ).poly
        centCoeff = [float(v) for k, v in orderPolyTable.iloc[0].items() if "cent_" in k]

        def trace(order, axisBPixels):
            axisBPixels = np.asarray(axisBPixels, dtype=float)
            return np.asarray(poly(pd.DataFrame({self.axisB: axisBPixels, "order": float(order)}), *centCoeff))

        return trace

    def _warn_on_join_slit_jumps(self, orderPixelTable, orderPolyTable, settings):
        """*warn when adjacent VIS orders put the fitted trace at different slit positions where they join*

        The check only logs a warning, which the extraction also keeps as a note, and changes nothing else.

        **Key Arguments:**

        - ``orderPixelTable`` -- the fitted trace table
        - ``orderPolyTable`` -- the fitted order-centre polynomial coefficients
        - ``settings`` -- the ``faint_trace_settings``, or None where the trace is off
        """
        if settings is None or "wavelength" not in orderPixelTable.columns:
            return
        mismatches = join_slit_mismatches(
            orderPixelTable,
            self._fitted_trace(orderPolyTable),
            self._make_pixels_at(),
            self.slitOffsetsArcsec,
            settings.maxJoinArcsec,
        )
        for mismatch in mismatches:
            blue = order_label(mismatch.blueOrder, self.arm)
            red = order_label(mismatch.redOrder, self.arm)
            self.log.warning(
                f"The fitted traces of orders {blue} and {red} differ by {mismatch.differenceArcsec:+.1f} arcsec at "
                f"their join, more than {settings.maxJoinArcsec:.1f} arcsec, so they may follow different sources. "
                "Check the object trace."
            )

    def _record_faint_trace_qc(self, reports):
        """*record the faint orders that were tried: a log line and an extraction note each, and a header QC*

        The QC counts only the orders whose rows were replaced. Nothing is written to the header when no order was
        transferred, so a bright frame keeps the headers it always had. A skipped order still gets its log line and
        extraction note, on purpose, as an audit trail of why its trace was left as measured.

        **Key Arguments:**

        - ``reports`` -- the ``transfer_report`` of each faint order that had an anchor
        """
        note = getattr(self.log, "note", None)
        for report in reports:
            line = report.describe(order_label(report.order, self.arm), order_label(report.anchorOrder, self.arm))
            self.log.print(line)
            if callable(note):
                note(line)
        transferred = tuple(report for report in reports if report.is_transferred)
        if not transferred:
            return
        self.faintTransferReports = transferred
        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="N ORDERS FAINT TRACE",
            qcValue=len(transferred),
            qcComment="Number of faint orders whose trace was transferred from a bright neighbour",
            obsDateUtc=self.dateObs,
            reductionDateUtc=datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S"),
            qcUnit=None,
            toHeader=True,
        )

    def create_pixel_arrays(self):
        """*create a pixel array for the approximate centre of each order*

        **Return:**

        - ``orderPixelTable`` -- a data-frame containing lines and associated pixel locations
        - `dmBinx` -- the dispersion map binning in x
        - `dmBiny` -- the dispersion map binning in y
        """
        self.log.debug("starting the ``create_pixel_arrays`` method")

        import numpy as np
        import pandas as pd
        from astropy.io import fits

        # READ THE SPECTRAL FORMAT TABLE TO DETERMINE THE LIMITS OF THE TRACES
        orderNums, waveLengthMin, waveLengthMax, amins, amaxs = read_spectral_format(
            log=self.log,
            settings=self.settings,
            arm=self.arm,
            dispersionMap=self.dispersion_map,
        )

        # READ ORDER SAMPLING RESOLUTION FROM SETTINGS
        sampleCount = self.recipeSettings["order-sample-count"]

        orderPixelRanges = []
        for o, amin, amax in zip(orderNums, amins, amaxs):
            arange = amax - amin
            orderPixelRanges.append(arange)

        smallestRange = min(orderPixelRanges)
        samplePixelSep = int(smallestRange / sampleCount)

        # CREATE THE WAVELENGTH/ORDER ARRAYS TO BE CONVERTED TO PIXELS
        myDict = {
            "order": np.asarray([]),
            "wavelength": np.asarray([]),
            "slit_position": np.asarray([]),
        }

        for o, wmin, wmax, pixelRange in zip(orderNums, waveLengthMin, waveLengthMax, orderPixelRanges):
            orderSampleCount = int(pixelRange / samplePixelSep)
            wrange = wmax - wmin

            wlArray = np.arange(wmin, wmax, (wmax - wmin) / orderSampleCount)
            myDict["wavelength"] = np.append(myDict["wavelength"], wlArray)
            myDict["order"] = np.append(myDict["order"], np.ones(len(wlArray)) * o)
            myDict["slit_position"] = np.append(myDict["slit_position"], np.zeros(len(wlArray)))

        orderPixelTable = pd.DataFrame(myDict)
        orderPixelTable = dispersion_map_to_pixel_arrays(
            log=self.log,
            dispersionMapPath=self.dispersion_map,
            orderPixelTable=orderPixelTable,
        )

        # GRAB HEADER FROM DISPERSION MAP
        with fits.open(self.dispersion_map, memmap=True) as hdul:
            header = hdul[0].header

        try:
            dmBinx = header[self.kw("WIN_BINX")]
            dmBiny = header[self.kw("WIN_BINY")]
        except KeyError as e:
            self.log.debug(f"create_pixel_arrays: `dmBinx = header[self.kw('WIN_BINX')]` failed, continuing: {e}")
            dmBinx = 1
            dmBiny = 1

        self.log.debug("completed the ``create_pixel_arrays`` method")
        return orderPixelTable, dmBinx, dmBiny

    def fit_1d_gaussian_to_slices(self, orderPixelTable, sliceLength, medianStddev=False):
        """Optimized version of Gaussian fitting to slices"""
        import numpy as np
        from astropy.modeling import fitting, models
        from astropy.stats import mad_std
        from scipy.signal import find_peaks

        if not medianStddev:
            medianStddev = 3.0

        # Pre-allocate arrays
        fitx = orderPixelTable["fit_x"].values
        fity = orderPixelTable["fit_y"].values
        n_points = len(fitx)

        contSliceAxis = np.full(n_points, np.nan)
        contSliceAntiAxis = np.full(n_points, np.nan)
        gauss_amplitude = np.full(n_points, np.nan)
        gauss_stddev = np.full(n_points, np.nan)
        gauss_mean = np.full(n_points, np.nan)

        # Create slice array once
        # OPTIMISE: 1s - 6464 hits
        slices = [
            cut_image_slice(
                log=self.log,
                frame=self.traceFrame,
                width=self.sliceWidth,
                length=sliceLength,
                x=x,
                y=y,
                sliceAxis=self.detectorParams["dispersion-axis"],
                median=True,
                debug=self.debug,
            )
            for x, y in zip(fitx, fity)
        ]

        # Initialize fitter once
        fit_g = fitting.LevMarLSQFitter()

        # Process each slice
        for i, s in enumerate(slices):
            slice_data, slice_length_offset, slice_width_centre = s[0], s[1], s[2]

            if slice_data is None:
                continue

            # Replace values < -500 with NaN
            slice_data = np.where(slice_data < -50, np.nan, slice_data)

            # Calculate baseline and std using vectorized operations
            baseline = np.nanpercentile(slice_data, 20)
            std_r = mad_std(slice_data, ignore_nan=True)

            if baseline < 0:
                baseline = 0

            if baseline is None:
                continue

            # Find peaks using vectorized operations
            peaks, _ = find_peaks(slice_data, height=baseline + self.peakSigmaLimit * std_r, width=1)

            if len(peaks) == 0:
                continue

            if len(peaks) > 1:
                # Vectorized closest peak finding
                peak_diffs = np.abs(peaks - sliceLength / 2.0)
                peaks = [peaks[np.argmin(peak_diffs)]]

            # Initialize and fit Gaussian
            # OPTIMISE: 1s
            g_init = models.Gaussian1D(amplitude=slice_data[peaks[0]], mean=peaks[0], stddev=medianStddev)

            try:
                mask = np.isfinite(slice_data)
                x_range = np.arange(len(slice_data))
                # OPTIMISE: 7s - 3185 hits
                g = fit_g(g_init, x_range[mask], slice_data[mask])

                if 0 <= g.mean.value <= sliceLength and g.stddev.value > 0.5:
                    contSliceAxis[i] = g.mean.value + max(0, slice_length_offset)
                    contSliceAntiAxis[i] = slice_width_centre
                    gauss_amplitude[i] = g.amplitude.value
                    gauss_stddev[i] = g.stddev.value
                    gauss_mean[i] = g.mean.value

            except (ValueError, TypeError, RuntimeError, IndexError) as e:
                self.log.debug(f"fit_1d_gaussian_to_slices: `mask = np.isfinite(slice_data)` failed, continuing: {e}")
                continue

        # Assign results back to DataFrame using vectorized operations
        orderPixelTable[f"cont_{self.sliceAxis}"] = contSliceAxis
        orderPixelTable[f"cont_{self.sliceAntiAxis}"] = contSliceAntiAxis
        orderPixelTable["gauss_amplitude"] = gauss_amplitude
        orderPixelTable["gauss_stddev"] = gauss_stddev
        orderPixelTable["gauss_mean"] = gauss_mean

        return orderPixelTable

    def _fitted_order_curves(self, orderPixelTable, orderPolyTable):
        """*evaluate the order-centre polynomial along each order, in detector pixel coordinates*

        **Key Arguments:**

        - ``orderPixelTable`` -- the pixel table with residuals of fits
        - ``orderPolyTable`` -- data-frame of order-location polynomial coeff

        **Return:**

        - ``curves`` -- one ``(order, axisAFit, axisBSamples, axisAStdFit)`` tuple for each order fitted on the detector
        """
        import numpy as np
        import pandas as pd

        if self.axisB == "y":
            axisALength = self.traceFrame.data.shape[1]
            axisBLength = self.traceFrame.data.shape[0]
        elif self.axisB == "x":
            axisALength = self.traceFrame.data.shape[0]
            axisBLength = self.traceFrame.data.shape[1]

        axisBLineList = np.arange(0, axisBLength, 3)

        poly = chebyshev_order_xy_polynomials(
            log=self.log,
            axisBCol=self.axisB,
            orderCol="order",
            orderDeg=self.orderDeg,
            axisBDeg=self.axisBDeg,
        ).poly
        for _index, row in orderPolyTable.iterrows():
            centCoeff = [float(v) for k, v in row.items() if "cent_" in k]
            stdCoeff = [float(v) for k, v in row.items() if "std_" in k]

        df = pd.DataFrame({f"{self.axisB}": axisBLineList})
        curves = []
        for o in orderPixelTable["order"].unique():
            df["order"] = o
            axisAFit = np.asarray(poly(df, *centCoeff))
            axisAStdFit = np.asarray(poly(df, *stdCoeff))

            keep = (axisAFit > 0) & (axisAFit < axisALength - 10)
            if not keep.any():
                self.log.debug(f"_fitted_order_curves: order {o} has no fitted positions on the detector, continuing")
                continue
            curves.append((o, axisAFit[keep], axisBLineList[keep], axisAStdFit[keep]))
        return curves

    def build_order_metadata(self, orderPixelTable, orderPolyTable):
        """*find the extent of each fitted order on the detector*

        The limits come from the order-centre fit alone, so the QC plot settings cannot change them.
        ``{axisA}min`` and ``{axisA}max`` are the smallest and largest fitted axis-A positions in detector pixels.

        **Key Arguments:**

        - ``orderPixelTable`` -- the pixel table with residuals of fits
        - ``orderPolyTable`` -- data-frame of order-location polynomial coeff

        **Return:**

        - ``orderMetaTable`` -- dataframe of order number and the min and max pixel position along both axes
        """
        self.log.debug("starting the ``build_order_metadata`` method")

        import pandas as pd

        curves = self._fitted_order_curves(orderPixelTable, orderPolyTable)
        orderMetaTable = pd.DataFrame(
            {
                "order": [o for o, _, _, _ in curves],
                f"{self.axisB}min": [min(yfit) for _, _, yfit, _ in curves],
                f"{self.axisB}max": [max(yfit) for _, _, yfit, _ in curves],
                f"{self.axisA}min": [min(xfit) for _, xfit, _, _ in curves],
                f"{self.axisA}max": [max(xfit) for _, xfit, _, _ in curves],
            }
        )

        self.log.debug("completed the ``build_order_metadata`` method")
        return orderMetaTable

    def check_trace_extrapolation(self, sampledTable: pd.DataFrame, acceptedTable: pd.DataFrame) -> None:
        """*log when the fitted trace of an order is extrapolated across a long stretch with no accepted samples*

        Records the longest sample-free stretch of each order as the ``SAMPLES GAP MAX O<n>`` QC. The QC is not
        written to the product header and has no acceptable range, so it never fails a recipe.

        **Key Arguments:**

        - ``sampledTable`` -- the table of every sample position, with ``order`` and ``wavelength`` columns
        - ``acceptedTable`` -- the samples kept after the trace fit, with the same columns

        **Usage:**

        ```python
        self.check_trace_extrapolation(backupOrderPixelTable, orderPixelTable)
        ```
        """
        self.log.debug("starting the ``check_trace_extrapolation`` method")

        maxFraction = self.recipeSettings.get("max-unsampled-trace-fraction", DEFAULT_MAX_UNSAMPLED_TRACE_FRACTION)
        reductionDate = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        for order, orderSamples in sampledTable.groupby("order"):
            acceptedWavelengths = acceptedTable.loc[acceptedTable["order"] == order, "wavelength"]
            fraction, startNm, endNm = longest_unsampled_stretch(orderSamples["wavelength"], acceptedWavelengths)
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName=f"SAMPLES GAP MAX O{int(order)}",
                qcValue=round(fraction, 3),
                qcComment="Longest stretch of the order with no accepted trace sample, as a fraction of the order",
                obsDateUtc=self.dateObs,
                reductionDateUtc=reductionDate,
                qcUnit=None,
                toHeader=False,
                qcOrder=int(order),
            )
            if fraction > maxFraction:
                self.log.info(
                    f"The trace of order {int(order)} has no accepted samples between {startNm:.1f} nm and "
                    f"{endNm:.1f} nm ({fraction:.0%} of the order), so the fitted trace is extrapolated there. "
                    "Check the object and order positions."
                )

        self.log.debug("completed the ``check_trace_extrapolation`` method")

    def _scatter_trace_peaks(self, axis, orderPixelTable):
        """*scatter the fitted peak positions, drawing the rows of a transferred faint trace apart in their own marker*

        **Key Arguments:**

        - ``axis`` -- the matplotlib axis to draw on
        - ``orderPixelTable`` -- the fitted trace table, with a ``faint_transfer`` column when an order was transferred
        """
        measured, transferred = orderPixelTable, orderPixelTable.iloc[0:0]
        isTransferred = transferred_row_mask(orderPixelTable)
        if isTransferred is not None:
            measured, transferred = orderPixelTable.loc[~isTransferred], orderPixelTable.loc[isTransferred]
        axis.scatter(
            measured[f"cont_{self.axisB}"],
            measured[f"cont_{self.axisA}"],
            marker="o",
            c="green",
            s=0.3,
            alpha=0.6,
            label="cross-dispersion 1D gaussian peak-position",
        )
        anchors = {report.order: report.anchorOrder for report in self.faintTransferReports}
        for order, rows in transferred.groupby("order"):
            anchor = anchors.get(order)
            source = "a bright neighbouring order" if anchor is None else f"order {int(anchor)}"
            axis.scatter(
                rows[f"cont_{self.axisB}"],
                rows[f"cont_{self.axisA}"],
                marker="^",
                c="tab:blue",
                s=1.5,
                alpha=0.6,
                label=f"trace transferred from {source}",
            )

    @closes_its_figures
    def plot_results(self, orderPixelTable, orderPolyTable, clippedData):
        """*generate a plot of the polynomial fits and residuals*

        **Key Arguments:**

        - ``orderPixelTable`` -- the pixel table with residuals of fits
        - ``orderPolyTable`` -- data-frame of order-location polynomial coeff
        - ``clippedData`` -- the sigma-clipped data

        **Return:**

        - ``filePath`` -- path to the plot pdf
        """
        self.log.debug("starting the ``plot_results`` method")

        import matplotlib.pyplot as plt
        import numpy as np

        from soxspipe.commonutils.toolkit import qc_settings_plot_tables

        arm = self.arm

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotateImage = self.detectorParams["rotate-qc-plot"]
        flipImage = self.detectorParams["flip-qc-plot"]

        # ROTATE THE IMAGE FOR BETTER LAYOUT
        rotatedImg = self.traceFrame.data
        if rotateImage:
            rotatedImg = np.rot90(rotatedImg, rotateImage / 90)
        if flipImage:
            rotatedImg = np.flipud(rotatedImg)
            if not rotateImage:
                # FLIPUD SENDS PIXEL ROW r TO ROW N-1-r. PLOT FROM COPIES SO THE CALLER'S COORDINATES STAY UNFLIPPED
                aLen = rotatedImg.shape[0]
                orderPixelTable = orderPixelTable.assign(
                    **{f"cont_{self.axisA}": aLen - 1 - orderPixelTable[f"cont_{self.axisA}"]}
                )
                clippedData = clippedData.assign(
                    **{
                        f"cont_{self.axisA}": aLen - 1 - clippedData[f"cont_{self.axisA}"],
                        f"fit_{self.axisA}": aLen - 1 - clippedData[f"fit_{self.axisA}"],
                    }
                )

        # REMOVE DUPLICATE ENTRIES IN COLUMN 'qc_name' AND KEEP THE LAST ENTRY
        self.qc = self.qc.drop_duplicates(subset=["qc_name"], keep="last")
        tableSettings = {**self.recipeSettings, **{"exptime": self.exptime}}

        # SIZE EACH IMAGE PANEL TO THE IMAGE ASPECT AND GIVE EACH TABLE ITS OWN FULL-WIDTH ROW SO THEY CANNOT OVERLAP
        imagePanelHeight = qc_image_panel_height(rotatedImg)
        qcTableHeight, settingsTableHeight = qc_table_heights(self.qc, tableSettings)
        heightRatios = [
            imagePanelHeight,
            imagePanelHeight,
            QC_RESIDUAL_ROW_HEIGHT_INCHES,
            QC_STRIP_ROW_HEIGHT_INCHES,
            qcTableHeight,
            settingsTableHeight,
        ]
        fig = plt.figure(
            figsize=(6, sum(heightRatios) + QC_PLOT_DECORATION_HEIGHT_INCHES),
            constrained_layout=True,
        )
        gs = fig.add_gridspec(6, 2, height_ratios=heightRatios)
        toprow = fig.add_subplot(gs[0, :])
        midrow = fig.add_subplot(gs[1, :])
        bottomleft = fig.add_subplot(gs[2, 0])
        bottomright = fig.add_subplot(gs[2, 1])
        fwhmaxis = fig.add_subplot(gs[3, :])
        settingsAx = fig.add_subplot(gs[5, :])
        qcAx = fig.add_subplot(gs[4, :])

        # INTERPOLATION "none" EMBEDS THE IMAGE AT NATIVE RESOLUTION IN THE PDF; RESAMPLING TO THE SAVE DPI
        # OFFSETS THE ORDER BANDS FROM THE VECTOR MARKERS WHEN ZOOMED
        toprow.imshow(rotatedImg, vmin=10, vmax=50, cmap="gray", alpha=0.5, interpolation="none")
        midrow.imshow(rotatedImg, vmin=10, vmax=50, cmap="gray", alpha=0.5, interpolation="none")
        toprow.set_title("1D gaussian peak positions (post-clipping)", fontsize=10)
        self._scatter_trace_peaks(toprow, orderPixelTable)

        # PLOT WHERE A CONTINUUM WAS NOT FOUND - WITH SLICE LENGTH SHOWN
        mask = clippedData["cont_x"].isna() | clippedData["cont_y"].isna()
        for axB, axA in zip(
            clippedData.loc[mask, f"fit_{self.axisB}"],
            clippedData.loc[mask, f"fit_{self.axisA}"],
        ):
            toprow.plot(
                [axB, axB],
                [axA - int(self.sliceLength / 2), axA + int(self.sliceLength / 2)],
                color="black",
                linewidth=0.5,
                alpha=0.3,
                label=("continuum not found" if axB == clippedData.loc[mask, f"fit_{self.axisB}"].iloc[0] else None),
            )
        toprow.scatter(
            clippedData[f"cont_{self.axisB}"],
            clippedData[f"cont_{self.axisA}"],
            marker="x",
            c="red",
            s=5,
            alpha=0.6,
            linewidths=0.5,
            label="peaks clipped during continuum fitting",
        )
        legend_below_axis(toprow, fig)

        # toprow.set_yticklabels([])
        # toprow.set_xticklabels([])

        toprow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        toprow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        toprow.tick_params(axis="both", which="major", labelsize=9)
        toprow.set_xlim([0, rotatedImg.shape[1]])

        if self.detectorParams["dispersion-axis"] == "x":
            toprow.invert_yaxis()
            midrow.invert_yaxis()
        toprow.set_ylim([0, rotatedImg.shape[0]])
        midrow.set_ylim([0, rotatedImg.shape[0]])

        if "order" in self.recipeName.lower():
            midrow.set_title("order-location fit solutions", fontsize=10)
        else:
            midrow.set_title("global polynomal fit of order-centres", fontsize=10)
        uniqueOrders = orderPixelTable["order"].unique()
        orderColours = _order_colours(uniqueOrders)
        labelAdded = None

        for o, xfit, yfit, stdfit in self._fitted_order_curves(orderPixelTable, orderPolyTable):
            lower = xfit - 3 * stdfit
            upper = xfit + 3 * stdfit
            if flipImage and not rotateImage:
                xfit = aLen - 1 - xfit
                lower = aLen - 1 - lower
                upper = aLen - 1 - upper
            if labelAdded == None:
                label1 = "$3\\sigma$ deviation"
                label2 = "polynomial fit"
                labelAdded = True
            else:
                label1 = None
                label2 = None
            midrow.plot(yfit, xfit, linewidth=0.7, label=label2, color=orderColours[o])

            midrow.fill_between(yfit, lower, upper, color=orderColours[o], alpha=0.3, label=label1)
            try:
                midrow.text(
                    yfit[10],
                    xfit[10] + 10,
                    int(o),
                    fontsize=6,
                    c=orderColours[o],
                    verticalalignment="bottom",
                )
            except (IndexError, KeyError) as e:
                self.log.debug(f"plot_results: `midrow.text( yfit[10], xfit[10] + 10, int(...` failed, continuing: {e}")

        # xfit = np.ones(len(xfit)) * \
        #     self.pinholeFrame.data.shape[1] - xfit
        # midrow.scatter(yfit, xfit, marker='x', c='blue', s=4)
        # midrow.set_yticklabels([])
        # midrow.set_xticklabels([])
        midrow.set_ylabel(f"{self.axisA}-axis", fontsize=12)
        midrow.set_xlabel(f"{self.axisB}-axis", fontsize=12)
        midrow.tick_params(axis="both", which="major", labelsize=9)

        legend_below_axis(midrow, fig)

        # PLOT THE FINAL RESULTS:
        for o in uniqueOrders:
            c = orderColours[o]
            mask = orderPixelTable["order"] == o
            bottomleft.scatter(
                orderPixelTable.loc[mask][f"cont_{self.axisA}"].values,
                orderPixelTable.loc[mask][f"cont_{self.axisA}_fit_res"].values,
                alpha=0.2,
                s=1,
                c=c,
            )
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                bottomleft.text(
                    orderPixelTable.loc[mask][f"cont_{self.axisA}"].values[labelIndex],
                    orderPixelTable.loc[mask][f"cont_{self.axisA}_fit_res"].values[labelIndex],
                    int(o),
                    fontsize=8,
                    c=c,
                    verticalalignment="bottom",
                )
            except (IndexError, KeyError) as e:
                self.log.debug(f"plot_results: `bottomleft.text( orderPixelTable.loc[mask]...` failed, continuing: {e}")

        bottomleft.set_xlabel(f"{self.axisA} pixel position", fontsize=10)
        bottomleft.set_ylabel(f"{self.axisA} residual", fontsize=10)
        bottomleft.tick_params(axis="both", which="major", labelsize=9)

        # PLOT THE FINAL RESULTS:
        for o in uniqueOrders:
            c = orderColours[o]
            mask = orderPixelTable["order"] == o
            bottomright.scatter(
                orderPixelTable.loc[mask][f"cont_{self.axisB}"].values,
                orderPixelTable.loc[mask][f"cont_{self.axisA}_fit_res"].values,
                alpha=0.2,
                s=1,
                c=c,
            )
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                bottomright.text(
                    orderPixelTable.loc[mask][f"cont_{self.axisB}"].values[labelIndex],
                    orderPixelTable.loc[mask][f"cont_{self.axisA}_fit_res"].values[labelIndex],
                    int(o),
                    fontsize=8,
                    c=c,
                    verticalalignment="bottom",
                )
            except (IndexError, KeyError) as e:
                self.log.debug(f"plot_results: `bottomright.text( orderPixelTable.loc[mask...` failed, continuing: {e}")

        bottomright.set_xlabel(f"{self.axisB} pixel position", fontsize=10)
        bottomright.tick_params(axis="both", which="major", labelsize=9)
        # bottomright.set_ylabel('x residual')
        bottomright.set_yticklabels([])

        stdToFwhm = 2 * (2 * math.log(2)) ** 0.5
        for o in uniqueOrders:
            c = orderColours[o]
            mask = orderPixelTable["order"] == o
            fwhmaxis.scatter(
                orderPixelTable.loc[mask]["wavelength"].values,
                orderPixelTable.loc[mask]["gauss_stddev_fit"].values * stdToFwhm,
                alpha=0.2,
                s=1,
                c=c,
            )
            if len(orderPixelTable.loc[mask].index) > 10:
                labelIndex = 10
            else:
                labelIndex = 1
            try:
                fwhmaxis.text(
                    orderPixelTable.loc[mask]["wavelength"].values[labelIndex],
                    orderPixelTable.loc[mask]["gauss_stddev_fit"].values[labelIndex] * stdToFwhm,
                    int(o),
                    fontsize=8,
                    c=c,
                    verticalalignment="bottom",
                )
            except (IndexError, KeyError) as e:
                self.log.debug(f"plot_results: `fwhmaxis.text( orderPixelTable.loc[mask]['...` failed, continuing: {e}")
        fwhmaxis.set_xlabel("wavelength (nm)", fontsize=10)
        fwhmaxis.set_ylabel("Cross-dispersion\nFWHM (pixels)", fontsize=10)

        fwhmaxis.set_ylim(
            orderPixelTable["gauss_stddev_fit"].min() * stdToFwhm * 0.5,
            orderPixelTable["gauss_stddev_fit"].max() * stdToFwhm * 1.2,
        )

        qc_settings_plot_tables(
            log=self.log,
            qc=self.qc,
            qcAx=qcAx,
            settings=tableSettings,
            settingsAx=settingsAx,
        )

        # THE TITLE STATISTICS COME FROM THE MEASURED ROWS, SO EXACT SYNTHETIC ROWS DO NOT FLATTER THE FIT
        statRes = measured_residuals(
            orderPixelTable[f"cont_{self.axisA}_fit_res"].values, transferred_row_mask(orderPixelTable)
        )
        mean_res = np.mean(np.abs(statRes))
        std_res = np.std(statRes)
        res_min = np.min(orderPixelTable[f"cont_{self.axisA}_fit_res"].values)
        res_max = np.max(orderPixelTable[f"cont_{self.axisA}_fit_res"].values)
        res_range = res_max - res_min

        subtitle = f"mean res: {mean_res:2.2f} pix, res stdev: {std_res:2.2f}"
        lamp = ""
        if self.lamp:
            lamp = f" {self.lamp} lamp"
        if "order" in self.recipeName.lower():
            fig.suptitle(
                f"traces of order-centre locations - {arm}{lamp} pinhole flat-frame\n{subtitle}",
                fontsize=12,
            )
        else:
            fig.suptitle(f"{arm} object trace locations\n{subtitle}", fontsize=12)

        polyOrders = [self.orderDeg, self.axisBDeg]
        polyOrders[:] = [str(l) for l in polyOrders]
        polyOrders = "".join(polyOrders)

        # plt.show()
        if not self.sofName:
            filename = filenamer(log=self.log, frame=self.traceFrame, settings=self.settings)
            filename = filename.split("FLAT")[0] + "ORDER_CENTRES_residuals.pdf"
        elif "order" in self.recipeName.lower():
            filename = self.sofName + f"_residuals_{polyOrders}.pdf"
        elif has_sequence_suffixed_traces(self.recipeName):
            sequence = self.noddingSequence + self.notFlattened
            filename = self.sofName + "_OBJECT_TRACE_residuals" + sequence + f"_{polyOrders}.pdf"
        else:
            filename = self.sofName + f"_OBJECT_TRACE_residuals{self.notFlattened}_{polyOrders}.pdf"

        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        if not self.settings["tune-pipeline"]:
            plt.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
        if self.debug:
            plt.show()
        plt.close(fig)
        plt.close("all")

        # REDUCE MEMORY USAGE
        cleanUp = [toprow, midrow, bottomleft, bottomright, fwhmaxis, settingsAx, qcAx]
        for item in cleanUp:
            try:
                item.cla()
                del item
            except AttributeError as e:
                self.log.debug(f"plot_results: `item.cla()` failed, continuing: {e}")

        if self.settings["tune-pipeline"]:
            import codecs

            filePath = "residuals.txt"
            exists = os.path.exists(filePath)
            if not exists:
                with codecs.open(filePath, encoding="utf-8", mode="w") as writeFile:
                    writeFile.write("polyOrders,mean_res,std_res,res_min,res_max,res_range \n")
            with codecs.open(filePath, encoding="utf-8", mode="a") as writeFile:
                writeFile.write(
                    f"{polyOrders},{mean_res:2.4f},{std_res:2.4f},{res_min:2.4f},{res_max:2.4f},{res_range:2.4f}\n"
                )

        self.log.debug("completed the ``plot_results`` method")
        return filePath

    @closes_its_figures
    def sample_trace(self):
        """*take many cross-dispersion samples across each order to try and find an object trace*

        **Return:**

        - ``orderPixelTable`` -- the detector locations at which a trace was found
        """
        self.log.debug("starting the ``sample_trace`` method")

        import numpy as np
        import pandas as pd

        # CONVERT WAVELENGTH TO PIXEL POSITIONS AND RETURN ARRAY OF POSITIONS TO
        # SAMPLE THE TRACES
        orderPixelTable, dmBinx, dmBiny = self.create_pixel_arrays()

        binx = 1
        biny = 1
        try:
            binx = self.traceFrame.header[self.kw("WIN_BINX")] / dmBinx
            biny = self.traceFrame.header[self.kw("WIN_BINY")] / dmBiny
        except KeyError as e:
            self.log.debug(f"sample_trace: `binx = self.traceFrame.header[self.kw('WIN_BIN...` failed, continuing: {e}")

        self.mapBinning = {"x": binx, "y": biny}

        # FIT_X AND FIT_Y FROM DISP-SOLUTION
        orderPixelTable["fit_x"] = orderPixelTable["fit_x"] / binx
        orderPixelTable["fit_y"] = orderPixelTable["fit_y"] / biny

        # SLICE LENGTH TO SAMPLE TRACES IN THE CROSS-DISPERSION DIRECTION
        self.sliceLength = self.recipeSettings["slice-length"]
        self.peakSigmaLimit = self.recipeSettings["peak-sigma-limit"]
        self.sliceWidth = self.recipeSettings["slice-width"]

        if (
            self.kw("DPR_TYPE").upper() in self.traceFrame.header
            and "STD" in self.traceFrame.header[self.kw("DPR_TYPE")].upper()
            and self.sliceWidth > 3
        ):
            self.sliceWidth = 3

        # PREP LISTS WITH NAN VALUE IN CONT_X AND CONT_Y BEFORE FITTING
        orderPixelTable[f"cont_{self.axisA}"] = np.nan
        orderPixelTable[f"cont_{self.axisB}"] = np.nan

        # FOR EACH ORDER, FOR EACH PIXEL POSITION SAMPLE, FIT A 1D GAUSSIAN IN
        # CROSS-DISPERSION DIRECTION. RETURN PEAK POSITIONS
        from soxspipe.commonutils.toolkit import quicklook_image

        quicklook_image(
            log=self.log,
            CCDObject=self.traceFrame,
            show=self.debug,
            ext="data",
            stdWindow=3,
            title=False,
            surfacePlot=False,
        )

        if "order" in self.recipeName.lower():
            self.log.print("\n# FINDING & FITTING ORDER-CENTRE CONTINUUM TRACES\n")
        else:
            self.log.print("\n# FINDING & FITTING OBJECT CONTINUUM TRACES\n")

        def find_centre_points(orderPixelTable, medianShift=False, medianStddev=False):
            """*find the central peak of the continuum trace*"""
            import numpy as np
            from astropy.stats import sigma_clip

            # UNIQUE ORDERS
            uniqueOrders = orderPixelTable["order"].unique()
            uniqueOrders = "&".join(map(str, sorted(uniqueOrders)))

            # WORK ON A PRIVATE COPY: THE PROBE CALLS PASS STRIDED SLICES AND DISCARD THE RESULT
            orderPixelTable = orderPixelTable.copy()

            sliceLength = self.sliceLength
            if np.isnan(medianShift):
                sliceLength = sliceLength * 1.2
                shifted = ""
            elif medianShift:
                shifted = f"\nShifted by {medianShift:2.1f} pixels"
                orderPixelTable[f"fit_{self.axisA}"] = orderPixelTable[f"fit_{self.axisA}"] - medianShift
                sliceLength = sliceLength - int(abs(medianShift))
                if sliceLength > medianStddev * 5 + 5:
                    sliceLength = int(medianStddev * 5 + 5)
                if sliceLength < 7:
                    sliceLength = 7

            else:
                shifted = ""

            sliceLengthStr = f"\nSlit Length = {str(sliceLength)} pixels"

            orderPixelTable = self.fit_1d_gaussian_to_slices(
                orderPixelTable=orderPixelTable,
                sliceLength=sliceLength,
                medianStddev=medianStddev,
            )

            if "pre-clipped" not in orderPixelTable.columns:
                orderPixelTable["pre-clipped"] = False

            # FILTER DATA FRAME
            # FIRST CREATE THE MASK

            mask = orderPixelTable["cont_x"] < 0
            orderPixelTable.loc[mask, "cont_x"] = np.nan
            orderPixelTable.loc[mask, "pre-clipped"] = True
            mask = orderPixelTable["cont_y"] < 0
            orderPixelTable.loc[mask, "cont_y"] = np.nan
            orderPixelTable.loc[mask, "pre-clipped"] = True
            mask = orderPixelTable["gauss_stddev"] < 0.3
            orderPixelTable.loc[mask, "pre-clipped"] = True

            increaseSlitLength = False
            for o in orderPixelTable["order"].unique():
                mask = (orderPixelTable["order"] == o) & (orderPixelTable["pre-clipped"] == False)
                if len(orderPixelTable.loc[mask].index) < 5:
                    increaseSlitLength = True
                    continue

            # FIND HOW FAR AWAY FROM THE CENTRE THE CONTINUUM WAS FOUND
            orderPixelTable["centre_shift"] = (
                orderPixelTable[f"fit_{self.axisA}"] - orderPixelTable[f"cont_{self.axisA}"]
            )

            # SIGMA-CLIP THE DATA ON STDDEV
            if "centre" not in self.recipeName:
                mask = orderPixelTable["pre-clipped"] == False
                masked_residuals = sigma_clip(
                    orderPixelTable.loc[mask]["gauss_stddev"],
                    sigma_lower=3,
                    sigma_upper=3,
                    maxiters=5,
                    cenfunc="mean",
                    stdfunc="std",
                )
                orderPixelTable.loc[mask, "pre-clipped"] = masked_residuals.mask
                mask = orderPixelTable["pre-clipped"] == False
                masked_residuals = sigma_clip(
                    orderPixelTable.loc[mask]["centre_shift"],
                    sigma_lower=3,
                    sigma_upper=3,
                    maxiters=5,
                    cenfunc="mean",
                    stdfunc="std",
                )
                orderPixelTable.loc[mask, "pre-clipped"] = masked_residuals.mask
            else:
                orderPixelTable["pre-clipped"] = False

            mask = orderPixelTable["pre-clipped"] == False
            filteredDf = orderPixelTable.loc[mask]
            medianShift = filteredDf["centre_shift"].median()
            medianStddev1 = filteredDf["gauss_stddev"].median()

            from astropy.stats import sigma_clipped_stats

            mean1, medianShift, std1 = sigma_clipped_stats(
                filteredDf["centre_shift"].values,
                sigma=2.5,
                stdfunc="std",
                cenfunc="mean",
                maxiters=5,
            )

            mean2, medianStddev, std2 = sigma_clipped_stats(
                filteredDf["gauss_stddev"].values,
                sigma=2.5,
                stdfunc="std",
                cenfunc="mean",
                maxiters=5,
            )

            if not std1 or ((abs(medianShift) + 2) / std1 < 6 and medianStddev < 1):
                medianShift = False
                medianStddev = False
            elif (abs(medianShift)) < 1.0:
                medianShift = False
            elif medianStddev < 1:
                medianStddev = False

            if self.debug:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(8, 5))
                plt.scatter(
                    orderPixelTable["wavelength"],
                    orderPixelTable["centre_shift"],
                    marker="o",
                    s=30,
                    label=f"Median Shift={medianShift:2.1f}, Median STD={medianStddev:2.1f}",
                )

                # LABEL X and Y AXES
                plt.xlabel("Wavelength (nm)", fontsize=12)
                plt.ylabel(f"{self.axisA} centre shift (pixels)", fontsize=12)
                plt.title(
                    f"{self.arm} continuum trace centre shifts (Orders {uniqueOrders}){sliceLengthStr}{shifted} FINDCENTRE",
                    fontsize=14,
                )

                filteredDf = orderPixelTable.loc[~mask]
                plt.scatter(
                    filteredDf["wavelength"],
                    filteredDf["centre_shift"],
                    marker="x",
                    s=30,
                )
                plt.legend()
                plt.show()

            if self.debug:
                import matplotlib.pyplot as plt

                plt.figure(figsize=(8, 5))
                plt.scatter(
                    orderPixelTable["wavelength"],
                    orderPixelTable["gauss_stddev"],
                    marker="o",
                    s=30,
                    label=f"Median Shift={medianShift:2.1f}, Median STD={medianStddev:2.1f}",
                )

                # LABEL X and Y AXES
                plt.xlabel("Wavelength (nm)", fontsize=12)
                plt.ylabel("Gaussian stddev (pixels)", fontsize=12)
                plt.title(
                    f"{self.arm} continuum trace Gaussian stddevs (Orders {uniqueOrders}){sliceLengthStr}{shifted} FINDCENTRE",
                    fontsize=14,
                )

                mask = orderPixelTable["pre-clipped"] == True
                filteredDf = orderPixelTable.loc[mask]
                plt.scatter(
                    filteredDf["wavelength"],
                    filteredDf["gauss_stddev"],
                    marker="x",
                    s=30,
                )
                plt.legend()
                plt.show()

            # DROP THE MASKED PEAKS?
            if False:
                mask = orderPixelTable["pre-clipped"] == False
                orderPixelTable = orderPixelTable.loc[mask]
                # REMOVE COLUMN FROM DATA FRAME
                # orderPixelTable.drop(columns=['mask'], inplace=True)

            if increaseSlitLength:
                medianShift = np.nan

            return orderPixelTable, medianShift, medianStddev

        def sample_order(orderTable: pd.DataFrame, everyN: int) -> tuple[pd.DataFrame, bool]:
            """*sample one order with its own shift and window*

            Returns the sampled order and whether a denser probe (``everyN`` of 1) is needed first. When it is, the
            returned table is the input, unsampled.
            """
            junk, medianShift, medianStddev = find_centre_points(
                orderPixelTable=orderTable.iloc[::everyN],
                medianShift=False,
                medianStddev=False,
            )

            if np.isnan(medianShift):
                junk, medianShift, medianStddev = find_centre_points(
                    orderPixelTable=orderTable.iloc[::everyN],
                    medianShift=medianShift,
                    medianStddev=False,
                )

            # `find_centre_points` RETURNS False FOR "NO SIGNIFICANT SHIFT" (UNDER 1 PIXEL OR NOISY) AND NaN FOR
            # "TOO FEW SAMPLES"; BOTH MEAN THE SPARSE PROBE FOUND NO USABLE SHIFT. A FLOAT 0.0 IS NEVER RETURNED
            isShiftMissing = medianShift is False or np.isnan(medianShift)
            if isShiftMissing and everyN > 1:
                return orderTable, True

            orderTable, medianShift, medianStddev = find_centre_points(
                orderPixelTable=orderTable,
                medianShift=medianShift,
                medianStddev=medianStddev,
            )
            return orderTable, False

        allLines = len(orderPixelTable.index)
        tmpOrderPixelTable = orderPixelTable.copy()

        everyN = 10
        sanityCheck = False

        while not sanityCheck:
            if self.inst.upper() == "SOXS" and self.arm.upper() == "VIS":
                # EACH VIS ORDER GETS ITS OWN SHIFT AND WINDOW: THE ORDERS SIT AT DIFFERENT OFFSETS FROM THE MAP,
                # SO A SHARED MEDIAN SHIFT LEAVES THE TRACE OF ONE ORDER OUTSIDE THE WINDOW OF ANOTHER (DY-1716)
                orderTables = []
                needsDenserProbe = False
                for order in sorted(tmpOrderPixelTable["order"].unique()):
                    orderTable, needsDenserProbe = sample_order(
                        tmpOrderPixelTable.loc[tmpOrderPixelTable["order"] == order].copy(),
                        everyN,
                    )
                    if needsDenserProbe:
                        break
                    orderTables.append(orderTable)

                if needsDenserProbe:
                    everyN = 1
                    continue

                # RECOMBINE
                orderPixelTable = pd.concat(orderTables, ignore_index=True)
            else:
                tmpOrderPixelTable = orderPixelTable.copy()
                junk, medianShift, medianStddev = find_centre_points(
                    orderPixelTable=tmpOrderPixelTable.iloc[::everyN],
                    medianShift=False,
                    medianStddev=False,
                )

                if np.isnan(medianShift):
                    junk, medianShift, medianStddev = find_centre_points(
                        orderPixelTable=tmpOrderPixelTable.iloc[::everyN],
                        medianShift=medianShift,
                        medianStddev=False,
                    )

                if (not (medianShift) or np.isnan(medianShift)) and everyN > 1:
                    sanityCheck = False
                    everyN = 1
                    # print("Reducing everyN to 1")
                    continue

                if np.isnan(medianShift):
                    self.log.error("FLUX IN THE INPUT FLAT FRAMES IS TOO LOW TO PROCEED. PLEASE CHECK THE RAW FRAMES")
                    raise ValueError("FLUX IN THE INPUT FLAT FRAMES IS TOO LOW TO PROCEED. PLEASE CHECK THE RAW FRAMES")

                orderPixelTable, medianShift, medianStddev = find_centre_points(
                    orderPixelTable=tmpOrderPixelTable,
                    medianShift=medianShift,
                    medianStddev=medianStddev,
                )

            sanityCheck = True

        # MASK orderPixelTable WHERE cont_x or cont_y is NaN
        mask = orderPixelTable["cont_x"].isna() | orderPixelTable["cont_y"].isna()
        foundLines = len(orderPixelTable.loc[~mask].index)
        percent = 100 * foundLines / allLines

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        self.qc = pd.concat(
            [
                self.qc,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "SAMPLES TOT NUM",
                            "qc_value": allLines,
                            "qc_comment": "Total number of samples along orders",
                            "qc_unit": None,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        self.qc = pd.concat(
            [
                self.qc,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "SAMPLES DET NUM",
                            "qc_value": foundLines,
                            "qc_comment": "Number of samples where a continuum is detected",
                            "qc_unit": None,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        self.psamp = percent / 100
        self.qc = pd.concat(
            [
                self.qc,
                pd.DataFrame(
                    [
                        {
                            "soxspipe_recipe": self.recipeName,
                            "qc_name": "SAMPLES DET FRAC",
                            "qc_value": round(float(self.psamp), 3),
                            "qc_comment": "Proportion of samples where a continuum is detected",
                            "qc_unit": None,
                            "obs_date_utc": self.dateObs,
                            "reduction_date_utc": utcnow,
                            "to_header": True,
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )

        self.log.print(f"\tContinuum found in {foundLines} out of {allLines} order slices ({percent:2.0f}%)")

        self.log.debug("completed the ``sample_trace`` method")

        return orderPixelTable, percent

    # use the tab-trigger below for new method
    # xt-class-method
