#!/usr/bin/env python
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

Author
: Marco Landoni & David Young

Date Created
: May 17, 2023
"""

import os
from collections.abc import Mapping
from dataclasses import dataclass
from types import MappingProxyType
from typing import TYPE_CHECKING, NamedTuple

from .base_util import base_util
from .extraction_notes import extraction_note_logger, format_extraction_notes
from .figure_guard import closes_its_figures
from .filenamer import filenamer

if TYPE_CHECKING:
    import pandas as pd

os.environ["TERM"] = "vt100"


class profile_clip_cap_rows(NamedTuple):
    """Slit rows of one order that stopped at the profile clipping iteration limit while still clipping."""

    refitted: int
    kept: int


# WIDTH IN COLUMNS OVER WHICH THE LOCAL FRACTION OF USABLE PIXELS IS MEASURED IN EACH SLIT ROW
PROFILE_SUPPORT_WINDOW = 51
# A MASKED PIXEL IN A ROW WITH LESS LOCAL USABLE FRACTION THAN THIS HAS AN UNCONSTRAINED PROFILE
PROFILE_MIN_LOCAL_SUPPORT = 0.5
# A COLUMN SPILLS IF THE USABLE ROW NEXT TO THE UNSUPPORTED EDGE HOLDS THIS FRACTION OF THE COLUMN PEAK
SLIT_EDGE_PEAK_FRACTION = 0.5
# AN ORDER SPILLS IF MORE THAN THIS FRACTION OF ITS COLUMNS ARE TRUNCATED AT THE SLIT EDGE
SLIT_EDGE_ORDER_FRACTION = 0.2

# SLIT POSITION (ARCSEC) ON WHICH THE BOXCAR FALLBACK IS CENTRED WHEN NO OBJECT TRACE WAS FOUND
SLIT_CENTRE_ARCSEC = 0.0

# SHORT ENOUGH TO FIT ON THE FITS HEADER CARD NEXT TO THE HIERARCH KEYWORD AND THE TEXT VALUE
EXTRACT_METHOD_QC_COMMENT = "Optimal, or boxcar if no trace"
SKY_TUNING_QC_COMMENT = "Wavelength tuning to sky lines"
SKY_TUNING_SKIPPED_QC_VALUE = "skipped: attached arc"

# SKY-LINE MATCH TOLERANCES (nm) FOR THE THREE TUNING ITERATIONS OF AN ORDER MEASURED FROM SCRATCH
SKYLINE_TOLERANCES_NM = (0.5, 0.1, 0.1)
# MATCH TOLERANCE (nm) FOR EACH ITERATION OF THE VIS r ORDER WHEN IT STARTS AT THE SHIFT PREDICTED FROM i (DY-1645)
R_SEEDED_TOLERANCE_NM = 0.05


@dataclass(frozen=True)
class FlatPassHandoff:
    """*what the flattened extraction of a frame hands to the not-flattened extraction of the same frame (DY-1680)*

    The not-flattened (response) extraction runs no skyline peak matching. It applies the per-order shifts that the
    flattened extraction of the same frame or nod pair measured. Both extractions trace the same unflattened frame, so
    the not-flattened extraction also reuses the flattened extraction's object trace instead of detecting it again
    (DY-1650).

    **Key Arguments:**

    - ``skyShiftsNm`` -- ``{order: wavelength shift (nm)}`` that the flattened extraction applied, or *None* when it
      skipped skyline tuning (attached arc, DY-1665)
    - ``objectTracePath`` -- path to the object-trace table the flattened extraction wrote, or *None* when it found
      no object trace
    """

    skyShiftsNm: Mapping[int | str, float] | None
    objectTracePath: str | None

    # COMPARABLE BUT NOT HASHABLE: THE SHIFTS ARE A READ-ONLY MAPPING, WHICH HAS NO HASH
    __hash__ = None

    def __post_init__(self) -> None:
        if self.skyShiftsNm is not None:
            object.__setattr__(self, "skyShiftsNm", MappingProxyType(dict(self.skyShiftsNm)))

    def __reduce__(self):
        # A MAPPINGPROXY CANNOT BE PICKLED, SO PICKLE (AND DEEP-COPY) THROUGH A PLAIN DICT
        shifts = None if self.skyShiftsNm is None else dict(self.skyShiftsNm)
        return (type(self), (shifts, self.objectTracePath))


# TODO: revisit how the wavelength for each slice is calculated ... take from the continuum, or the central 3-5 pixels?


class horne_extraction(base_util):
    """
    *perform optimal source extraction using the Horne method (Horne 1986)*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``recipeSettings`` -- the recipe specific settings
    - ``skySubtractedFrame`` -- path to sky subtracted frame
    - ``unflattenedFrame`` -- path to unflattened frame
    - ``subtractedFrame`` -- path to the frame that has been subtracted in the AB or BA sequence. This is use to extract a sky spectrum. Default *False* (i.e. not used)
    - ``twoDMapPath`` -- path to 2D dispersion map image path
    - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
    - ``qcTable`` -- the data frame to collect measured QC metrics
    - ``productsTable`` -- the data frame to collect output products (if False no products are saved to file)
    - ``dispersionMap`` -- the FITS binary table containing dispersion map polynomial
    - ``sofName`` -- the set-of-files filename
    - ``locationSetIndex`` -- the index of the AB cycle locations, or the offset stack's own index (nodding and
      offset modes only). Default *False*
    - ``startNightDate`` -- YYYY-MM-DD date of the observation night. Default ""
    - ``notFlattened`` -- flag to indicate if the frame is flattened or not. Default *False*
    - ``debug`` -- flag to indicate if debug mode is on (shows plots). Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.
    - ``noTraceFallback`` -- if no object trace is found, fall back to a robust boxcar extraction centred on the slit
      centre (slit position 0 arcsec) instead of giving up. The fallback is flagged by ``traceFallback``, a warning in
      ``extractionNotes``, the ``EXTRACT METHOD`` QC and the plot label. Default *False* (no spectrum is extracted
      when there is no trace)
    - ``scienceHeaders`` -- the headers of the raw science frames of the observation block. The attached-arc test
      measures the arc from the first frame's start and the last frame's end (DY-1666). Default *False* uses the
      header of ``skySubtractedFrame``
    - ``flatPass`` -- the ``FlatPassHandoff`` of the flattened extraction of the same frame, read only when
      ``notFlattened`` is set. The not-flattened extraction applies its sky shifts and runs no skyline peak matching
      (DY-1680), and it reuses its object trace instead of detecting the trace again (DY-1650). Default *None* leaves
      the not-flattened wavelengths untuned and detects the trace again

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (see tutorial here https://fundamentals.readthedocs.io/en/main/initialisation.html).

    To initiate a horne_extraction object, use the following:

    ```python
    from soxspipe.commonutils import horne_extraction
    optimalExtractor = horne_extraction(
        log=log,
        skySubtractedFrame=skySubtractedFrame,
        unflattenedFrame=unflattenedFrame,
        subtractedFrame=subtractedFrame,
        twoDMapPath=twoDMap,
        settings=settings,
        recipeName="soxs-stare",
        qcTable=qc,
        productsTable=products,
        dispersionMap=dispMap,
        sofName=sofName,
        locationSetIndex=locationSetIndex,
        startNightDate=startNightDate,
        debug=debug,
        turnOffMP=turnOffMP

    )
    qc, products = optimalExtractor.extract()
    ```

    """

    # SET TRUE IN __init__ WHEN NO OBJECT TRACE WAS FOUND AND noTraceFallback ASKED FOR A BOXCAR EXTRACTION INSTEAD
    noTraceFallback = False
    traceFallback = False

    # {ORDER: SLIT-ILLUMINATION OFFSET (arcsec)} FOR THE FAINT-ORDER TRACE, OR None WHEN NOT MEASURED
    slitOffsetsArcsec = None

    # BLEND-MERGED VIS CALIBRATION SKYLINES, BUILT ON FIRST USE (SEE ``_merged_skylines``)
    _mergedSkylinesDF = None

    # THE FLATTENED PASS'S HAND-OFF READ BY A NOT-FLATTENED EXTRACTION (__init__), AND THE ONE THIS EXTRACTION BUILT
    flatPass = None
    flatPassHandoff = None
    # {ORDER: SHIFT (nm)} APPLIED BY THE LAST SKYLINE TUNING, OR None WHEN TUNING WAS SKIPPED
    skyShiftsNm = None
    # PATH TO THE OBJECT-TRACE TABLE THIS EXTRACTION USED, OR None WHEN NO TRACE WAS FOUND
    objectTracePath = None
    # THE TRACE'S PIXEL TABLE, OR None WHEN THERE IS NOTHING TO EXTRACT
    orderPixelTable = None

    def __init__(
        self,
        log,
        settings,
        recipeSettings,
        skySubtractedFrame,
        unflattenedFrame,
        twoDMapPath,
        subtractedFrame=False,
        recipeName=False,
        qcTable=False,
        productsTable=False,
        dispersionMap=False,
        sofName=False,
        locationSetIndex=False,
        startNightDate="",
        notFlattened=False,
        debug=False,
        turnOffMP=False,
        noTraceFallback=False,
        scienceHeaders=False,
        flatPass=None,
        slitOffsetsArcsec=None,
    ):
        from astropy import units as u
        from astropy.nddata import CCDData

        from soxspipe.commonutils.toolkit import nodding_sequence_suffix

        # COLLECT THE WARNINGS AND ERRORS OF THE EXTRACTION AND THE HELPERS IT BUILDS, TO PRINT ON THE QC PLOTS
        log = extraction_note_logger(log)

        super().__init__(
            log, settings, associatedFrame=skySubtractedFrame, dispersionMap=dispersionMap, twoDMapPath=twoDMapPath
        )

        log.debug("instantiating a new 'horne_extraction' object")
        self.twoDMapPath = twoDMapPath
        self.products = productsTable
        self.qc = qcTable
        self.slitEdgeOrders = []
        self.recipeName = recipeName
        self.sofName = sofName
        self.noddingSequence = ""
        self.recipeSettings = recipeSettings
        self.startNightDate = startNightDate
        self.debug = debug
        self.unflattenedFrame = unflattenedFrame
        self.turnOffMP = turnOffMP
        self.subtractedFrame = subtractedFrame
        self.noTraceFallback = noTraceFallback
        self.flatPass = flatPass
        self.slitOffsetsArcsec = slitOffsetsArcsec

        if notFlattened:
            self.notFlattened = "_NOTFLAT"
        else:
            self.notFlattened = ""

        if notFlattened:
            self.skySubtractedFrame = unflattenedFrame

        # DETECTING SEQUENCE AUTOMATICALLY
        self.noddingSequence = nodding_sequence_suffix(
            header=skySubtractedFrame.header, locationSetIndex=locationSetIndex, recipeName=recipeName
        )

        # COLLECT SETTINGS FROM SETTINGS FILE
        # HALF LENGTH KEPT AS THE SETTING / 2 (MAY END IN .5) SO THE EXTRACTED SLIT IS AS LONG AS THE SETTING
        self.slitHalfLength = self.recipeSettings["horne-extraction-slit-length"] / 2
        self.clippingSigma = self.recipeSettings["horne-extraction-profile-clipping-sigma"]
        self.clippingIterationLimit = self.recipeSettings["horne-extraction-profile-clipping-iteration-count"]
        self.globalClippingSigma = self.recipeSettings["horne-extraction-profile-global-clipping-sigma"]

        # OPEN THE SKY-SUBTRACTED FRAME
        if isinstance(skySubtractedFrame, CCDData):
            self.skySubtractedFrame = skySubtractedFrame
        else:
            self.skySubtractedFrame = CCDData.read(
                skySubtractedFrame,
                hdu=0,
                unit=u.electron,
                hdu_uncertainty="ERRS",
                hdu_mask="QUAL",
                hdu_flags="FLAGS",
                key_uncertainty_type="UTYPE",
            )

        # THE RAW SCIENCE FRAMES OF THE OBSERVATION BLOCK, WHOSE FIRST START AND LAST END THE ATTACHED-ARC TEST USES
        self.scienceHeaders = list(scienceHeaders) if scienceHeaders else [self.skySubtractedFrame.header]

        # THE INSTRUMENT, READ BEFORE THE TRACE BRANCH SO THE NO-TRACE FALLBACK HAS IT TOO
        self.inst = self.skySubtractedFrame.header.get(self.kw("INSTRUME"), "")

        # USE SUBTRACTED FRAME (NODDING ONLY) TO EXTRACT A SKY SPECTRUM
        if subtractedFrame == False:
            self.subtractedFrame = None
        else:
            # OPEN THE SKY-MODEL FRAME
            if isinstance(subtractedFrame, CCDData):
                self.subtractedFrame = subtractedFrame
            else:
                self.subtractedFrame = CCDData.read(
                    subtractedFrame,
                    hdu=0,
                    unit=u.electron,
                    hdu_uncertainty="ERRS",
                    hdu_mask="QUAL",
                    hdu_flags="FLAGS",
                    key_uncertainty_type="UTYPE",
                )

        # GET A TEMPLATE FILENAME USED TO NAME PRODUCTS
        if self.sofName:
            self.filenameTemplate = self.sofName + ".fits"
        else:
            self.filenameTemplate = filenamer(log=self.log, frame=self.skySubtractedFrame, settings=self.settings)

        from soxspipe.commonutils.toolkit import utility_setup

        self.qcDir, self.productDir = utility_setup(
            log=self.log,
            settings=settings,
            recipeName=recipeName,
            startNightDate=startNightDate,
        )

        # ADJUST SLIT HEIGHT IF BINNING IN USE
        if self.binx > 1 or self.biny > 1:
            if self.detectorParams["dispersion-axis"] == "x":
                self.slitHalfLength /= self.binx
            else:
                self.slitHalfLength /= self.biny

        # REMOVE ZEROS
        mask = (self.imageMap["wavelength"] == 0) & (self.imageMap["slit_position"] == 0)
        self.imageMap = self.imageMap.loc[~mask]

        if not self._find_object_trace(locationSetIndex):
            return

        # CHECK DIFFERENCE IN TIME BETWEEN OBJECT AND MPH FRAME
        mjdDispMap = self.twoDMap["WAVELENGTH"].header[self.kw("MJDOBS")]
        mjdObject = self.skySubtractedFrame.header[self.kw("MJDOBS")]
        self.log.print(f"\nTime difference between object and dispersion map: {mjdObject - mjdDispMap:.4f} days")

        # xpd-update-filter-dataframe-column-values

    def _find_object_trace(self, locationSetIndex) -> bool:
        """*find the object trace in each order, or reuse the trace the flattened extraction found (DY-1650)*

        Sets ``orderPixelTable`` and ``objectTracePath``. A flattened extraction also starts its ``flatPassHandoff``
        with the trace, so the not-flattened extraction of the same frame can reuse it even when this extraction finds
        no trace and stops early. A not-flattened extraction with a ``flatPass`` runs no continuum detection and adds
        no trace files, products or QCs.

        **Key Arguments:**

        - ``locationSetIndex`` -- the index of the AB cycle locations, or the offset stack's own index

        **Return:**

        - ``hasTrace`` -- *False* when no trace was found and there is no slit-centre fallback, so nothing can be
          extracted
        """
        from soxspipe.commonutils import detect_continuum
        from soxspipe.commonutils.toolkit import unpack_order_table

        if self.notFlattened and self.flatPass is not None:
            # BOTH PASSES TRACE THE SAME UNFLATTENED FRAME, SO THE FLATTENED PASS'S TRACE IS READ, NOT FOUND AGAIN
            productPath = self.flatPass.objectTracePath
            self.orderPixelTable = None
        else:
            detector = detect_continuum(
                log=self.log,
                traceFrame=self.unflattenedFrame,
                dispersion_map=self.dispersionMap,
                settings=self.settings,
                recipeSettings=self.recipeSettings,
                sofName=self.sofName,
                recipeName=self.recipeName,
                qcTable=self.qc,
                productsTable=self.products,
                locationSetIndex=locationSetIndex,
                startNightDate=self.startNightDate,
                debug=self.debug,
                notFlattened=bool(self.notFlattened),
                slitOffsetsArcsec=self.slitOffsetsArcsec,
            )
            productPath, self.qc, self.products, _, self.orderPixelTable, _ = detector.get()

        self.objectTracePath = productPath
        if not self.notFlattened:
            self.flatPassHandoff = FlatPassHandoff(skyShiftsNm=None, objectTracePath=productPath)

        if productPath is None:
            # FAILED TO FIND THE TRACE
            if not self.noTraceFallback:
                return False
            self._use_slit_centre_boxcar_fallback()
        else:
            self._unpack_trace(productPath, unpack_order_table)
        return True

    def _use_slit_centre_boxcar_fallback(self):
        """*prepare a robust boxcar extraction at the slit centre, for when no object trace was found*

        Sets ``orderPixelTable`` to a table of every order (the slit centre is then used for all of them), flags the
        fallback and logs a warning that reaches the QC plots.
        """
        import pandas as pd

        self.orderPixelTable = pd.DataFrame({"order": self.orderNums})
        self.traceFallback = True
        self.log.warning(
            "No object trace was found, so optimal extraction is not possible. "
            f"The spectrum is a robust boxcar extraction of a {2 * self.slitHalfLength:g}-pixel aperture "
            "centred on the slit centre (slit position 0 arcsec). Check that the object is on the slit."
        )

    def _record_extract_method_qc(self):
        """*add the ``EXTRACT METHOD`` QC: ``optimal``, or ``boxcar`` when no object trace was found*

        Only the flattened extraction of a recipe that asks for the fallback records it, so the nodding and offset
        recipes and the unflattened re-extraction are unchanged.
        """
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        if not self.noTraceFallback or self.notFlattened or isinstance(self.qc, bool):
            return
        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="EXTRACT METHOD",
            qcValue="boxcar" if self.traceFallback else "optimal",
            qcComment=EXTRACT_METHOD_QC_COMMENT,
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow_string(),
            qcUnit=None,
            toHeader=True,
        )

    def _unpack_trace(self, productPath, unpack_order_table):
        """*read the detected object trace and drop the orders of a SOXS VIS continuum that was badly fitted*

        **Key Arguments:**

        - ``productPath`` -- path to the object-trace order table
        - ``unpack_order_table`` -- the toolkit function, imported by the caller
        """
        # UNPACK THE ORDER TABLE
        orderPolyTable, self.orderPixelTable, orderMetaTable = unpack_order_table(
            log=self.log, orderTablePath=productPath
        )

        # ORDER CHECK IN CASE OF POOR CONTINUUM FITTING
        # GET UNIQUE VALUES IN COLUMN
        if self.arm.upper() == "VIS" and self.inst.upper() == "SOXS":
            keepOrders = []
            uniqueOrders = self.orderPixelTable["order"].unique()
            for o in uniqueOrders:
                mask = self.orderPixelTable["order"] == o
                if self.orderPixelTable.loc[mask][f"{self.axisA}coord_centre"].std() < 100:
                    keepOrders.append(o)
                else:
                    self.log.warning(f"Bad continuum fit to order {o}; this order will not be extracted")
            mask = self.orderPixelTable["order"].isin(keepOrders)
            self.orderPixelTable = self.orderPixelTable.loc[mask]

    def extract(self):
        """*extract the full spectrum order-by-order and return FITS Binary table containing order-merged spectrum*

        **Return:**

        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products
        - ``mergedSpectumDF`` -- path to the FITS binary table containing the merged spectrum
        """
        self.log.debug("starting the ``extract`` method")

        if self.orderPixelTable is None:
            self.log.error("No trace found in the data, optimal extraction cannot be performed.")
            return self.qc, self.products, None, None, None

        self._record_extract_method_qc()

        import copy
        from contextlib import suppress

        # MAKE RELATIVE HOME PATH ABSOLUTE
        from datetime import datetime

        import pandas as pd
        from astropy.table import Table

        from soxspipe.commonutils.phase3 import write_fits_table_to_disk
        from soxspipe.commonutils.toolkit import (
            add_snr_efficiency_qcs,
        )

        kw = self.kw
        arm = self.arm

        extractions = []

        self.log.print("\n# PERFORMING OPTIMAL SOURCE EXTRACTION (Horne Method)")

        self.log.print("\tBuilding wavelength, slit-position, flux, error and bad-pixel arrays")

        # BUILD THE DICTIONARY OF ARRAYS TO ZOOM AND REBIN — KEYS BECOME COLUMN NAMES IN THE ORDER TABLE
        # THIS IS (D-S) IN HORNE 1986, I.E. THE SKY-SUBTRACTED FLUX
        frameArrays = {
            "wavelength": self.twoDMap["WAVELENGTH"].data,
            "fluxRaw": self.skySubtractedFrame.data,
            "slit": self.twoDMap["SLIT"].data,
            "sliceError": self.skySubtractedFrame.uncertainty.array,
        }
        if self.subtractedFrame:
            frameArrays["fluxSky"] = self.subtractedFrame.data

        from soxspipe.commonutils.image_transformer import image_transformer

        # ZOOM AND REBIN ALL ARRAYS ORDER-BY-ORDER, ALSO SIGMA-CLIPPING THE BAD-PIXEL MASK
        transformer = image_transformer(
            log=self.log,
            settings=self.settings,
            orderPixelTable=self.orderPixelTable,
            twoDMapPath=self.twoDMapPath,
            dispersionMap=self.dispersionMap,
            associatedFrame=self.skySubtractedFrame,
            slitHalfLength=self.slitHalfLength,
            fixedSlitCentreArcsec=SLIT_CENTRE_ARCSEC if self.traceFallback else None,
        )
        # THE DRIFT IS MEASURED FROM THE FITTED OBJECT TRACE, WHICH A FALLBACK EXTRACTION DOES NOT HAVE
        if not self.traceFallback:
            self.plot_slit_drift_qc(transformer=transformer)
        transformer.cache_image("fluxRaw", self.skySubtractedFrame.data, associatedMask=self.skySubtractedFrame.mask)
        # VARIANCE IS RECTIFIED WITH THE SAME LINEAR AREA WEIGHTS AS FLUX, NOT THEIR SQUARES. NEIGHBOURING
        # RECTIFIED PIXELS SHARE DETECTOR PIXELS, AND THE EXTRACTION SUMS ACROSS THEM: LINEAR WEIGHTS APPROXIMATE
        # THAT SUM'S VARIANCE BY COUNTING EACH DETECTOR PIXEL'S VARIANCE ABOUT ONCE, WHEREAS SQUARED WEIGHTS
        # DROP THE COVARIANCE AND UNDERESTIMATE THE NOISE (~1.55x INFLATED SNR ON REAL DATA)
        # THE SAME BAD-PIXEL MASK KEEPS FLAGGED PIXELS OUT OF THE VARIANCE SUMS TOO,
        # SO FLUX AND VARIANCE CELLS ARE RENORMALISED ALIKE
        transformer.cache_image(
            "variance", self.skySubtractedFrame.uncertainty.array**2, associatedMask=self.skySubtractedFrame.mask
        )
        if self.subtractedFrame:
            transformer.cache_image("fluxSky", self.subtractedFrame.data)

        orderSlices = transformer.get_order_slices()
        wlMinMax = transformer.get_order_wavelength_ranges()
        orderRectifiedImages = transformer.get_order_rectified()

        # IF NO SKY FRAME WAS PROVIDED, FILL THE SKY COLUMN WITH ZEROS
        if not self.subtractedFrame:
            for orderTable in orderSlices:
                orderTable["fluxSky"] = [0] * len(orderTable)

        self.log.print(f"\tExtracting {len(orderSlices)} orders\n\n")

        from fundamentals import fmultiprocess

        inputArray = [[s, r] for s, r in zip(orderSlices, orderRectifiedImages)]

        # EXTRACTION RUNS IN-PROCESS (turnOffMP=True BELOW). WARNINGS LOGGED INSIDE WORKER PROCESSES WOULD NOT REACH
        # THE PLOT NOTES, BECAUSE EACH WORKER RECORDS INTO ITS OWN COPY OF self.log
        extractions = fmultiprocess(
            log=self.log,
            function=extract_single_order,
            inputArray=inputArray,
            poolSize=False,
            timeout=300,
            funclog=self.log,
            slitHalfLength=self.slitHalfLength,
            clippingSigma=self.clippingSigma,
            clippingIterationLimit=self.clippingIterationLimit,
            globalClippingSigma=self.globalClippingSigma,
            axisA=self.axisA,
            axisB=self.axisB,
            hornePolyOrder=self.recipeSettings["horne-extraction-profile-poly-order"],
            debug=self.debug,
            boxcarOnly=self.traceFallback,
            turnOffMP=True,
            progressBar=True,
        )

        updatedExtractions = []
        for e, wlTuple in zip(extractions, wlMinMax):
            # FILTER DATA FRAME
            # FIRST CREATE THE MASK
            if e is not None:
                mask = (e["wavelengthMean"] > wlTuple[0]) & (e["wavelengthMean"] < wlTuple[1])
                e = e.loc[mask]
                updatedExtractions.append(e)

        extractions = self._summarise_slit_edge_spill(updatedExtractions)

        self.plot_extracted_spectrum_qc(extractions=extractions)

        # MERGE THE ORDER SPECTRA
        extractedOrdersDF = pd.concat(extractions, ignore_index=True)
        extractedOrdersDF = extractedOrdersDF.loc[
            extractedOrdersDF["pixelScaleNm"] < 3
        ]  # FILTER OUT ANY REMAINING BAD PIXELS WITH UNREALISTICALLY LARGE PIXEL SCALE (I.E. WAVELENGTH JUMPS BETWEEN ADJACENT PIXELS)
        extractedOrdersDF = self._correct_wavelengths(extractedOrdersDF, arm)

        mergedSpectumDF, orderJoins = self.merge_extracted_orders(extractedOrdersDF)

        self.qc = add_snr_efficiency_qcs(
            log=self.log,
            spectrumDF=mergedSpectumDF,
            qcTable=self.qc,
            orderJoins=orderJoins,
            recipeName=self.recipeName,
            dateObs=self.dateObs,
            arm=arm,
        )

        if not isinstance(self.products, bool):
            # CONVERT TO FITS BINARY TABLE
            header = copy.deepcopy(self.skySubtractedFrame.header)
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

            # header["HIERARCH " + kw("PRO_TECH")] = header.pop(kw("DPR_TECH"))

            header[kw("SEQ_ARM")] = arm
            header["HIERARCH " + kw("PRO_TYPE")] = "REDUCED"
            header["HIERARCH " + kw("PRO_CATG")] = f"SCI_SLIT_FLUX_{arm}".upper()

            # DISCRETE ORDERS
            filename = self.filenameTemplate.replace(
                ".fits",
                f"_EXTRACTED_ORDERS{self.noddingSequence}{self.notFlattened}.fits",
            )
            filePath = f"{self.productDir}/{filename}"

            write_fits_table_to_disk(
                log=self.log,
                settings=self.settings,
                header=header,
                tables=[extractedOrdersDF],
                filePath=filePath,
                qc=self.qc,
            )

            utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

            if not self.notFlattened:
                self.products = pd.concat(
                    [
                        self.products,
                        pd.DataFrame(
                            [
                                {
                                    "soxspipe_recipe": self.recipeName,
                                    "product_label": f"EXTRACTED_ORDERS_TABLE{self.noddingSequence}{self.notFlattened}",
                                    "file_name": filename,
                                    "file_type": "FITS",
                                    "obs_date_utc": self.dateObs,
                                    "reduction_date_utc": utcnow,
                                    "product_desc": "Table of the extracted source in each order",
                                    "file_path": filePath,
                                    "label": "PROD",
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

            # NOW MERGED SPECTRUM
            filename = self.filenameTemplate.replace(
                ".fits",
                f"_EXTRACTED_MERGED{self.noddingSequence}{self.notFlattened}.fits",
            )
            filePath = f"{self.productDir}/{filename}"

            def _round_scalar_or_quantity(value, decimals):
                if hasattr(value, "value"):
                    value = value.value
                return round(float(value), decimals)

            for col, decimals in [
                ("WAVE", 2),
                ("FLUX_COUNTS", 3),
                ("SNR", 2),
                ("FLUX_DENSITY_COUNTS", 3),
            ]:
                # ONLY THE NOT-FLATTENED SPECTRUM HAS A FLUX DENSITY COLUMN
                if col not in mergedSpectumDF.columns:
                    continue
                mergedSpectumDF[col] = mergedSpectumDF[col].apply(lambda x: _round_scalar_or_quantity(x, decimals))

            mergedTable = Table.from_pandas(mergedSpectumDF)
            write_fits_table_to_disk(
                log=self.log,
                settings=self.settings,
                header=header,
                tables=[mergedTable],
                filePath=filePath,
                qc=self.qc,
            )

            # EXPORTING SPECTRUM IN ASCII FORMAT

            # CHECKING IF WE ARE IN A NODDING SEQUENCE
            if self.noddingSequence or len(self.noddingSequence) > 0 and not self.notFlattened:
                pass
            elif not self.notFlattened:
                # SAVE THE MERGED ASTROPY TABLE TO TXT FILE
                # SAVE THE TABLE stackedSpectrum TO DISK IN ASCII FORMAT
                asciiFilepath = filePath.replace(".fits", ".txt")
                asciiFilename = filename.replace(".fits", ".txt")
                mergedTable["WAVE"] = mergedTable["WAVE"] * 10  # CONVERTING TO ANGSTROMS
                mergedTable["WAVE"].format = "{:.2f}"
                mergedTable.write(asciiFilepath, format="ascii", overwrite=True)

                self.products = pd.concat(
                    [
                        self.products,
                        pd.DataFrame(
                            [
                                {
                                    "soxspipe_recipe": self.recipeName,
                                    "product_label": f"EXTRACTED_MERGED_ASCII{self.notFlattened}",
                                    "file_name": asciiFilename,
                                    "file_type": "TXT",
                                    "obs_date_utc": self.dateObs,
                                    "reduction_date_utc": utcnow,
                                    "product_desc": "Ascii version of extracted source spectrum",
                                    "file_path": asciiFilepath,
                                    "label": "PROD",
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )

                self.products = pd.concat(
                    [
                        self.products,
                        pd.DataFrame(
                            [
                                {
                                    "soxspipe_recipe": self.recipeName,
                                    "product_label": f"EXTRACTED_MERGED_TABLE{self.noddingSequence}{self.notFlattened}",
                                    "file_name": filename,
                                    "file_type": "FITS",
                                    "obs_date_utc": self.dateObs,
                                    "reduction_date_utc": utcnow,
                                    "product_desc": "Table of the extracted, order-merged",
                                    "file_path": filePath,
                                    "label": "PROD",
                                }
                            ]
                        ),
                    ],
                    ignore_index=True,
                )
            extractionFilepath = filePath
        else:
            extractionFilepath = False

        self.log.debug("completed the ``extract`` method")
        return self.qc, self.products, mergedSpectumDF, orderJoins, extractionFilepath

    @property
    def extractionNotes(self) -> list[tuple[str, str]]:  # noqa: N802
        """*the warnings and errors logged during this extraction*

        **Return:**

        - ``notes`` -- a list of ``(level, message)`` pairs, empty if the extractor's log is not a note collector
        """
        if isinstance(self.log, extraction_note_logger):
            return self.log.notes
        return []

    def _summarise_slit_edge_spill(self, extractions):
        """*find the orders where the object spills over the slit edge, report them and drop the per-slice flag*

        Sets ``self.slitEdgeOrders``. Unless this is the unflattened re-extraction, also logs one error naming
        the spilled orders and adds the ``N ORDERS SLIT EDGE`` QC.

        **Key Arguments:**

        - ``extractions`` -- the list of per-order extraction data-frames, each with a ``slitEdgeTruncated`` column

        **Return:**

        - ``extractions`` -- the data-frames without the ``slitEdgeTruncated`` column
        """
        self.slitEdgeOrders = []
        for e in extractions:
            if e.empty:
                continue
            fraction = float(e["slitEdgeTruncated"].mean())
            if fraction > SLIT_EDGE_ORDER_FRACTION:
                self.slitEdgeOrders.append(
                    (
                        int(e["order"].iloc[0]),
                        float(e["wavelengthMean"].min()),
                        float(e["wavelengthMean"].max()),
                        fraction,
                    )
                )

        if not self.notFlattened:
            self._report_slit_edge_spill()

        return [e.drop(columns=["slitEdgeTruncated"]) for e in extractions]

    def _report_slit_edge_spill(self):
        """*log the slit-edge spill error and add the ``N ORDERS SLIT EDGE`` QC*"""
        from soxspipe.commonutils.toolkit import SOXS_VIS_ORDER_BANDS, append_qc, utcnow_string

        if self.slitEdgeOrders:
            bands = SOXS_VIS_ORDER_BANDS if self.arm.upper() == "VIS" else {}
            spilledParts = []
            for order, wlMin, wlMax, fraction in self.slitEdgeOrders:
                band = f"{bands[order]}, " if order in bands else ""
                spilledParts.append(f"{order} ({band}{wlMin:0.1f}-{wlMax:0.1f} nm, {fraction:.0%} of slices)")
            spilled = ", ".join(spilledParts)
            self.log.error(
                f"The object spills over the slit edge in order(s) {spilled}. "
                "Flux beyond the slit edge is lost and the extraction in these orders is unreliable."
            )

        if isinstance(self.qc, bool):
            return
        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName="N ORDERS SLIT EDGE",
            qcValue=len(self.slitEdgeOrders),
            qcComment="Number of orders where the object spills over the slit edge",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow_string(),
            qcUnit=None,
            toHeader=True,
        )

    def weighted_average(self, group):
        import numpy as np

        group["FLUX_COUNTS"] = (group["flux_resampled"] * np.abs(group["snr"])).sum() / np.abs(group["snr"]).sum()
        # group['wf'] = group['flux_resampled'].mean()

        return group["FLUX_COUNTS"]

    def residual_merge(self, group):
        import numpy as np

        # residual = np.abs(group[0]['flux_resampled']) -  np.abs(group[1]['flux_resampled'])
        # if residual <= 0:
        #     group['wf'] = group[0]['flux_resampled']
        # else:
        #     group['wf'] = group[1]['flux_resampled']
        if len(group) > 1:
            if np.abs(group.iloc[0]["flux_resampled"]) >= np.abs(group.iloc[1]["flux_resampled"]):
                group["FLUX_COUNTS"] = group.iloc[1]["flux_resampled"]
            else:
                group["FLUX_COUNTS"] = group.iloc[0]["flux_resampled"]
        else:
            group["FLUX_COUNTS"] = group.iloc[0]["flux_resampled"]

        return group["FLUX_COUNTS"]

    def _arc_attachment(self):
        """*the offsets ``(minutes, rotatorDegrees)`` of an attached SOXS arc from the science block, or None*

        The arc is the dispersion map's primary header and the block is ``self.scienceHeaders``. X-Shooter, an
        unreadable header and an arc that is not attached (``within_attached_arc_limits``) all give *None*.
        """
        from soxspipe.commonutils.ug_wavelength_transfer import attached_arc_offsets, within_attached_arc_limits

        if self.inst.upper() != "SOXS":
            return None
        arcHeader = self._dispersion_map_header()
        if arcHeader is None:
            return None
        offsets = attached_arc_offsets(self.scienceHeaders, arcHeader, self.kw)
        return offsets if within_attached_arc_limits(offsets) else None

    def _correct_wavelengths(self, extractedOrdersDF, arm):
        """*tune the flattened extraction to the sky lines, or give the not-flattened one the flattened pass's shifts*

        The flattened extraction tunes (or skips, ``_tune_or_skip``) and keeps the result as ``flatPassHandoff``. The
        not-flattened extraction runs no skyline peak matching and records no sky QCs, so the reported SKY SHIFT QCs
        are those of the flattened extraction (DY-1680).

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument

        **Return:**

        - ``extractedOrdersDF`` -- the wavelength-corrected data-frame
        """
        if self.notFlattened:
            return self._apply_flat_pass_shifts(extractedOrdersDF, arm)
        # CLEARED FIRST, SO A REPEATED extract() NEVER HANDS ON AN EARLIER CALL'S SHIFTS
        self.skyShiftsNm = None
        tuned = self._tune_or_skip(extractedOrdersDF, arm)
        self.flatPassHandoff = FlatPassHandoff(skyShiftsNm=self.skyShiftsNm, objectTracePath=self.objectTracePath)
        return tuned

    def _apply_flat_pass_shifts(self, extractedOrdersDF, arm):
        """*shift each order by the sky shift the flattened extraction of the same frame applied (DY-1680)*

        Orders the flattened pass did not shift are left as they are. With no hand-off, or a flattened pass that
        skipped tuning, the spectrum is returned unchanged.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument

        **Return:**

        - ``extractedOrdersDF`` -- the shifted data-frame
        """
        if self.flatPass is None:
            self.log.warning(
                "No sky shifts were handed on from the flattened extraction, so the not-flattened extraction is "
                "left untuned"
            )
            return extractedOrdersDF
        if self.flatPass.skyShiftsNm is None:
            self.log.print("\tSkyline tuning was skipped in the flattened extraction, so it is skipped here too")
            return extractedOrdersDF

        return self._shift_orders(
            extractedOrdersDF, self.flatPass.skyShiftsNm.items(), arm, source=" (from the flattened extraction)"
        )

    def _shift_orders(self, extractedOrdersDF, orderShifts, arm, source=""):
        """*add each order's wavelength shift to its ``wavelengthMean`` rows, logging one line per order*

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders (shifted in place and returned)
        - ``orderShifts`` -- ``(order, shift nm)`` pairs; the order ``"all"`` shifts every order
        - ``arm`` -- the arm of the instrument
        - ``source`` -- text naming where the shifts came from, for the log. Default ""

        **Return:**

        - ``extractedOrdersDF`` -- the shifted data-frame
        """
        from soxspipe.commonutils.toolkit import order_label

        for o, s in orderShifts:
            mask = self._order_mask(extractedOrdersDF, o)
            orderText = o if o == "all" else order_label(o, arm)
            self.log.print(f"\t\t{s:0.2f} nm shift{source} applied to order {orderText}")
            extractedOrdersDF.loc[mask, "wavelengthMean"] += s
        return extractedOrdersDF

    def _tune_or_skip(self, extractedOrdersDF, arm):
        """*tune the wavelength solution to the sky lines, unless a SOXS arc is attached to the science frame*

        With an attached arc (``_arc_attachment``) the arc solution is accurate enough that tuning only adds noise
        (DY-1665), so the spectrum is returned unchanged and the skip is logged, noted on the QC plots and recorded
        as the ``SKY TUNING`` QC.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument

        **Return:**

        - ``extractedOrdersDF`` -- the tuned data-frame, or the input if tuning was skipped
        """
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        attachment = self._arc_attachment()
        if attachment is None:
            tuned = self.tune_wavelength_calibration_to_skylines(extractedOrdersDF, arm=arm)
            qcValue = "skylines"
        else:
            message = (
                f"Attached arc (Δt {attachment[0]:0.1f} min, Δrot {attachment[1]:0.1f}°): "
                "skyline wavelength tuning skipped"
            )
            self.log.print(message)
            if isinstance(self.log, extraction_note_logger):
                self.log.note(message)
            tuned = extractedOrdersDF
            qcValue = SKY_TUNING_SKIPPED_QC_VALUE
            self.skyShiftsNm = None

        if not isinstance(self.qc, bool):
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName="SKY TUNING",
                qcValue=qcValue,
                qcComment=SKY_TUNING_QC_COMMENT,
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow_string(),
                qcUnit=None,
                toHeader=True,
            )
        return tuned

    def tune_wavelength_calibration_to_skylines(self, extractedOrdersDF, arm, byOrder=True):
        """Skyline-based wavelength calibration tuning helper.

        Not called for a SOXS arc attached to the science frame (``_tune_or_skip``, DY-1665): its wavelength
        solution is accurate enough that tuning only adds noise.

        In the SOXS VIS arm, tuned order by order, the i order is measured before r, and no order is shifted from
        fewer than two sky lines or from a sky too faint for its lines to be centred (``MIN_SKY_LINE_SNR``, DY-1663;
        see ``_tune_vis_order``). The VIS u and g orders are never measured.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders
        - ``arm`` -- the arm of the instrument
        - ``byOrder`` -- whether to tune calibration order by order, or globally (default: True)

        **Return:**

        - ``extractedOrdersDF`` -- wavelength-corrected dataframe
        """

        self.log.debug("starting the ``tune_wavelength_calibration_to_skylines`` method")

        import numpy as np

        from soxspipe.commonutils.ug_wavelength_transfer import PREDICTED_ORDERS, VIS_I_ORDER

        uniqueOrders = ["all"] if not byOrder else np.sort(extractedOrdersDF["order"].unique())
        isVis = arm.upper() == "VIS"
        # THE DY-1645 VIS BEHAVIOUR IS FOR SOXS ONLY, AND ONLY WHEN TUNING ORDER BY ORDER
        isSoxsVis = byOrder and self._is_soxs_vis(arm)
        dispersions = self._order_dispersions(extractedOrdersDF, uniqueOrders, isSoxsVis) if byOrder else {}
        shiftsNm = dict.fromkeys(uniqueOrders, 0.0)
        # SKY LINES KEPT, SHIFT SCATTER AND (SOXS VIS) SKY-LINE SNRS PER MEASURED ORDER (THE r AND u/g GUARDS)
        orderStats = {}

        # IN VIS i IS MEASURED FIRST, BECAUSE THE r SHIFT IS PREDICTED FROM IT
        sequence = sorted(uniqueOrders, key=lambda o: not (isSoxsVis and o == VIS_I_ORDER))
        for o in sequence:
            orderDF = self._order_slice(extractedOrdersDF, o)
            if orderDF.empty or (isVis and o in PREDICTED_ORDERS):
                continue
            if isSoxsVis:
                shiftsNm[o], stats = self._tune_vis_order(orderDF, o, shiftsNm, orderStats, dispersions)
            else:
                shiftsNm[o], keptShifts = self._measure_order_shift(orderDF, o)
                stats = self._order_stats(keptShifts)
            if stats is not None:
                orderStats[o] = stats

        if byOrder:
            orderShifts = self._finalise_order_shifts(arm, uniqueOrders, shiftsNm, orderStats, dispersions)
        else:
            orderShifts = [shiftsNm[o] for o in uniqueOrders]

        # APPLY PER-ORDER SHIFTS TO WAVELENGTH COLUMN IN FULL DATAFRAME
        self.skyShiftsNm = {
            (o if o == "all" else int(o)): float(s) for o, s in zip(uniqueOrders, orderShifts, strict=True)
        }
        self._shift_orders(extractedOrdersDF, zip(uniqueOrders, orderShifts, strict=True), arm)

        return extractedOrdersDF

    def _is_soxs_vis(self, arm):
        """Return whether the arm is the SOXS VIS arm, the only one with the DY-1645 sky-tuning behaviour."""
        return arm.upper() == "VIS" and self.inst.upper() == "SOXS"

    def _order_dispersions(self, extractedOrdersDF, uniqueOrders, tolerateReadErrors):
        """*the signed dispersion of each order, optionally returning none if the dispersion map cannot be read*

        An unreadable dispersion map only costs the SOXS VIS r prediction, so there (``tolerateReadErrors``) a map
        that is missing, not a FITS file or holds no table is logged and no dispersions are returned. Only that read
        is guarded: any error in the dispersion computation itself goes through, and so does every error on the
        other arms, as it always has.
        """
        from soxspipe.commonutils.dispersion_map_to_pixel_arrays import read_dispersion_map_axes

        if tolerateReadErrors and self.dispersionMap:
            try:
                read_dispersion_map_axes(self.dispersionMap)
            except (OSError, ValueError) as e:
                self.log.info(
                    f"could not read the dispersion map to sign the order dispersions ({type(e).__name__}: {e})"
                )
                return {}
        return self._signed_order_dispersions(extractedOrdersDF, uniqueOrders)

    @staticmethod
    def _order_mask(extractedOrdersDF, order):
        """Return the row mask of one order; the dummy mask for ``"all"`` selects every order."""
        return extractedOrdersDF["order"] < 100 if order == "all" else extractedOrdersDF["order"] == order

    def _order_slice(self, extractedOrdersDF, order):
        """Return the rows of one order (every order for ``"all"``), sorted by wavelength."""
        mask = self._order_mask(extractedOrdersDF, order)
        return extractedOrdersDF.loc[mask].sort_values("wavelengthMean", kind="stable").copy()

    @staticmethod
    def _order_stats(keptShifts):
        """Return the sky lines kept and their shift scatter (nm) for the last iteration, or None if unmeasured."""
        import numpy as np

        if keptShifts is None:
            return None
        keptShifts = np.ravel(keptShifts)
        return {
            "nLines": int(keptShifts.size),
            "scatterNm": float(np.std(keptShifts)) if keptShifts.size else float("nan"),
        }

    def _measure_order_shift(self, orderDF, o, seedNm=0.0, tolerances=SKYLINE_TOLERANCES_NM):
        """*measure one order's wavelength shift from the sky lines*

        **Key Arguments:**

        - ``orderDF`` -- the order's rows, sorted by wavelength
        - ``o`` -- the order
        - ``seedNm`` -- the shift (nm) the first iteration starts from. Default *0*
        - ``tolerances`` -- the catalogue match tolerance (nm) of each iteration

        **Return:**

        - ``orderShift`` -- the total wavelength shift (nm), including the seed
        - ``keptShifts`` -- the per-line shifts kept in the last iteration, or *None* if no iteration found data
        """
        # PRIME THE SHIFT COLUMN WITH THE ORIGINAL VALUES, PLUS THE SEED
        orderDF["wavelength_shifted"] = orderDF["wavelengthMean"] + seedNm
        orderShift = seedNm
        keptShifts = None
        for iteration, tolerance in enumerate(tolerances):
            result = self._measure_shift_iteration(orderDF, o, tolerance, iteration)
            if result is None:
                continue
            medianShift, keptShifts = result
            # APPLY MEDIAN SHIFT TO CURRENT ITERATION
            orderDF["wavelength_shifted"] += medianShift
            orderShift += medianShift

        return orderShift, keptShifts

    def _measure_shift_iteration(self, orderDF, o, tolerance, iteration):
        """*match the sky lines of one tuning iteration and measure the median shift*

        **Key Arguments:**

        - ``orderDF`` -- the order's rows, with the current ``wavelength_shifted`` column
        - ``o`` -- the order
        - ``tolerance`` -- the catalogue match tolerance (nm)
        - ``iteration`` -- the iteration number, for the diagnostic plot

        **Return:**

        - ``(medianShift, keptShifts)`` -- the median shift (nm) and the per-line shifts kept after clipping, or
          *None* if the order holds no valid data
        """
        import numpy as np

        # EXTRACT NUMERIC ARRAYS FROM ORDER DATAFRAME
        wave, sky, objectFlux = self._extract_order_arrays(orderDF)
        valid = wave.notna() & sky.notna()
        if not valid.any():
            return None

        # DETECT PEAKS IN SMOOTHED SKY SPECTRUM
        skyValsOriginal, skyVals, peaks, waveVals, objectVals = self._detect_sky_peaks(wave, sky, objectFlux, valid)
        wmin, wmax = np.nanmin(waveVals), np.nanmax(waveVals)

        # MAP CATALOGUE SKYLINES TO WAVELENGTH COORDINATES FOR THIS ORDER
        localSkylines, calibrationSkylines = self._get_local_skylines_for_order(wmin, wmax, o, "wavelength")

        # MEASURE SUB-SAMPLE LINE CENTRES, SO A FLAT-TOPPED WIDE-SLIT LINE IS NOT CENTRED ON A SHOULDER
        peakCentres = self._measure_peak_centres(waveVals, skyVals, peaks)

        # MATCH LINE CENTRES TO CALIBRATION CATALOGUE SKYLINES
        matchedSkylinePixels, matchedShifts = self._match_peaks_to_skylines(
            peakCentres, calibrationSkylines, tolerance=tolerance
        )

        # SIGMA-CLIP SHIFT DISTRIBUTION AND COMPUTE MEDIAN
        clippedWave, clippedShifts, medianShift, keptShifts = self._compute_clipped_shift(
            matchedSkylinePixels, matchedShifts
        )

        # DIAGNOSTIC PLOT SHOWING SKY, SHIFTS, AND OBJECT SPECTRA
        if self.debug:
            self._plot_skyline_shift_diagnostic(
                waveVals,
                skyValsOriginal,
                skyVals,
                peaks,
                objectVals,
                localSkylines,
                calibrationSkylines,
                matchedSkylinePixels,
                matchedShifts,
                clippedWave,
                clippedShifts,
                medianShift,
                o,
                wmin,
                wmax,
                np.median(orderDF["pixelScaleNm"]),
                "wavelength (nm)",
                iteration,
                "nm",
            )

        return medianShift, keptShifts

    def _tune_vis_order(self, orderDF, o, shiftsNm, orderStats, dispersions):
        """*measure one VIS order's shift, refusing a shift from fewer than two sky lines or from a faint sky*

        The r order may instead be predicted from the already measured i order (``_tune_vis_r_order``).

        **Key Arguments:**

        - ``orderDF`` -- the order's rows, sorted by wavelength
        - ``o`` -- the order
        - ``shiftsNm`` -- ``{order: wavelength shift (nm)}`` of the orders measured so far
        - ``orderStats`` -- ``{order: {"nLines": int, "scatterNm": float, "lineSnr": float, ...}}`` of the orders
          measured so far
        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``

        **Return:**

        - ``shiftNm`` -- the order's wavelength shift (nm)
        - ``stats`` -- the order's sky-line statistics, or *None* if no data
        """
        from soxspipe.commonutils.ug_wavelength_transfer import VIS_R_ORDER

        if o == VIS_R_ORDER:
            return self._tune_vis_r_order(orderDF, shiftsNm, orderStats, dispersions)
        shiftNm, keptShifts = self._measure_order_shift(orderDF, o)
        stats = self._with_line_snr(orderDF, o, self._order_stats(keptShifts))
        reason = self._unusable_shift_reason(o, stats)
        if reason is not None:
            self.log.info(f"{reason}; no shift applied")
            shiftNm = 0.0
        return shiftNm, stats

    def _with_line_snr(self, orderDF, o, stats):
        """*return the order's sky-line statistics with its sky-line SNRs added, or None if it has no data*

        ``lineSnr`` is the combined SNR of the order's calibration lines and ``secondLineSnr`` the SNR of its
        second-brightest one (0 when it has fewer than two).
        """
        from soxspipe.commonutils.sky_line_snr import combine_line_snrs

        if stats is None:
            return None
        snrs = self._order_sky_line_snrs(orderDF, o)
        return {
            **stats,
            "lineSnr": combine_line_snrs(snrs),
            "secondLineSnr": float(snrs[1]) if snrs.size > 1 else 0.0,
        }

    def _order_sky_line_snrs(self, orderDF, o):
        """*the photon-noise SNR of each of the order's calibration sky lines at its arc-solution position, brightest
        first*"""
        import numpy as np
        import pandas as pd

        from soxspipe.commonutils.sky_line_snr import sky_line_snrs

        wave = pd.to_numeric(orderDF["wavelengthMean"], errors="coerce")
        sky = pd.to_numeric(orderDF["skyFlux"], errors="coerce")
        if not wave.notna().any():
            return np.array([])
        _, calibrationSkylines = self._get_local_skylines_for_order(wave.min(), wave.max(), o, "wavelength")
        return sky_line_snrs(wave, sky, calibrationSkylines, self.detectorParams["ron"])

    @staticmethod
    def _unusable_shift_reason(o, stats):
        """Return why an order's measured shift must not be applied, or None if enough sky lines were kept from a
        bright enough sky."""
        from soxspipe.commonutils.ug_wavelength_transfer import MIN_SHIFT_LINES, sky_line_snr_reason

        nLines = 0 if stats is None else stats["nLines"]
        if nLines < MIN_SHIFT_LINES:
            return f"order {o}: {nLines} sky lines kept (min {MIN_SHIFT_LINES})"
        return sky_line_snr_reason(o, stats["lineSnr"])

    def _tune_vis_r_order(self, orderDF, shiftsNm, orderStats, dispersions):
        """*decide the VIS r shift from its own sky lines and the shift predicted from i (DY-1645)*

        When i passes its guards and the science rotator angle is usable, r is predicted from i. r is then measured
        only within ``R_SEEDED_TOLERANCE_NM`` of that prediction. With two or more r lines kept, the measurement and
        the prediction are combined by inverse variance; with fewer, the prediction is used. Otherwise r is
        measured on its own and applied only when two or more lines are kept.

        **Key Arguments:**

        - ``orderDF`` -- the r order's rows, sorted by wavelength
        - ``shiftsNm`` -- ``{order: wavelength shift (nm)}``, holding the i shift
        - ``orderStats`` -- ``{order: {"nLines": int, "scatterNm": float, "lineSnr": float, ...}}``, holding the i
          statistics
        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``

        **Return:**

        - ``shiftNm`` -- the r wavelength shift (nm)
        - ``stats`` -- the r sky-line statistics, flagged ``transferred`` when the i prediction contributed
        """
        from soxspipe.commonutils import ug_wavelength_transfer as ug_transfer

        coeffs = ug_transfer.load_r_transfer_coefficients(self.detectorParams)
        predictedPx, reason = self._predict_r_shift_px(shiftsNm, orderStats, dispersions, coeffs)
        if predictedPx is None:
            return self._measure_r_alone(orderDF, coeffs is not None, reason)

        dispersion = dispersions[ug_transfer.VIS_R_ORDER]
        predictedNm = ug_transfer.detector_px_to_nm_shift(predictedPx, dispersion)
        measuredNm, keptShifts = self._measure_order_shift(
            orderDF,
            ug_transfer.VIS_R_ORDER,
            seedNm=predictedNm,
            tolerances=(R_SEEDED_TOLERANCE_NM,) * len(SKYLINE_TOLERANCES_NM),
        )
        stats = self._order_stats(keptShifts) or {"nLines": 0, "scatterNm": float("nan")}
        stats = self._with_line_snr(orderDF, ug_transfer.VIS_R_ORDER, stats)
        nLines = stats["nLines"]
        snrReason = ug_transfer.sky_line_snr_reason(ug_transfer.VIS_R_ORDER, stats["lineSnr"])
        if nLines < ug_transfer.MIN_SHIFT_LINES:
            source, shiftNm = "predicted", predictedNm
            comment = f"predicted from i and ABSROT; {nLines} sky lines kept"
        elif snrReason is not None:
            source, shiftNm = "predicted", predictedNm
            comment = f"predicted from i and ABSROT; {snrReason}"
        else:
            measuredPx = ug_transfer.nm_shift_to_detector_px(measuredNm, dispersion)
            sigmaPx = stats["scatterNm"] / abs(dispersion) / nLines**0.5
            combinedPx = ug_transfer.combine_r_shift_px(measuredPx, sigmaPx, predictedPx, coeffs["rms_px"])
            source, shiftNm = "combined", ug_transfer.detector_px_to_nm_shift(combinedPx, dispersion)
            comment = f"measured ({nLines} lines) combined with i prediction"
        self._record_shift_source_qc("R SHIFT SOURCE", source, comment)
        return shiftNm, {**stats, "transferred": True}

    def _measure_r_alone(self, orderDF, hasCoefficients, reason):
        """Measure the VIS r order without a prediction, applying it only when two or more sky lines are kept from a
        bright enough sky whose second-brightest line is above the noise."""
        from soxspipe.commonutils.ug_wavelength_transfer import VIS_R_ORDER, r_second_line_snr_reason

        shiftNm, keptShifts = self._measure_order_shift(orderDF, VIS_R_ORDER)
        stats = self._with_line_snr(orderDF, VIS_R_ORDER, self._order_stats(keptShifts))
        # NO STATS GIVES A "0 SKY LINES" REASON, SO THE SECOND-LINE CHECK ONLY EVER READS A MEASURED ORDER
        unusableReason = self._unusable_shift_reason(VIS_R_ORDER, stats) or r_second_line_snr_reason(
            stats["secondLineSnr"]
        )
        if hasCoefficients:
            self.log.info(f"r shift not predicted from i: {reason}")
        if unusableReason is None:
            source, comment = "measured", f"measured from {stats['nLines']} sky lines"
        else:
            self.log.info(f"{unusableReason}; no shift applied")
            source, comment, shiftNm = "none", unusableReason, 0.0
        if hasCoefficients:
            self._record_shift_source_qc("R SHIFT SOURCE", source, comment)
        return shiftNm, stats

    def _predict_r_shift_px(self, shiftsNm, orderStats, dispersions, coeffs):
        """*predict the r detector shift (px) from the i shift, or explain why it cannot be*

        **Key Arguments:**

        - ``shiftsNm`` -- ``{order: wavelength shift (nm)}``, holding the i shift
        - ``orderStats`` -- ``{order: {"nLines": int, "scatterNm": float, "lineSnr": float, ...}}``, holding the i
          statistics
        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``
        - ``coeffs`` -- the ``r-from-i-transfer`` coefficients, or *None* for an arm without them

        **Return:**

        - ``predictedPx`` -- the predicted r detector shift (px), or *None*
        - ``reason`` -- why there is no prediction, or *None*
        """
        from soxspipe.commonutils import ug_wavelength_transfer as ug_transfer

        if coeffs is None:
            return None, "the arm has no r-from-i coefficients"
        missing = [o for o in (ug_transfer.VIS_I_ORDER, ug_transfer.VIS_R_ORDER) if o not in dispersions]
        if missing:
            return None, f"order {missing[0]}: no dispersion"
        iStats = orderStats.get(ug_transfer.VIS_I_ORDER)
        iDispersion = dispersions[ug_transfer.VIS_I_ORDER]
        if iStats is not None:
            iStats = {"nLines": iStats["nLines"], "lineSnr": iStats["lineSnr"]}
        scienceAngle = ug_transfer.usable_rotator_angle(self.skySubtractedFrame.header, self.kw)
        reason = ug_transfer.r_prediction_skip_reason(iStats, scienceAngle)
        if reason is not None:
            return None, reason
        deltaIpx = ug_transfer.nm_shift_to_detector_px(shiftsNm[ug_transfer.VIS_I_ORDER], iDispersion)
        calibrationAngle = ug_transfer.calibration_rotator_angle(self._dispersion_map_header() or {}, self.kw)
        return ug_transfer.predict_r_detector_shift(deltaIpx, scienceAngle, calibrationAngle, coeffs), None

    def _extract_order_arrays(self, orderDF):
        """Extract numeric series for wave, spatial coord, sky, and object flux from an order slice."""
        import numpy as np
        import pandas as pd

        wave = pd.to_numeric(orderDF["wavelength_shifted"], errors="coerce")
        # axisBcoord = pd.to_numeric(orderDF[f"{self.axisB}coord"], errors="coerce")
        sky = pd.to_numeric(orderDF["skyFlux"], errors="coerce")
        if "extractedFluxOptimal" in orderDF.columns:
            objectFlux = pd.to_numeric(orderDF["extractedFluxOptimal"], errors="coerce")
        else:
            objectFlux = pd.Series(np.nan, index=orderDF.index)

        return wave, sky, objectFlux

    def _detect_sky_peaks(self, wave, sky, objectFlux, valid):
        """Smooth sky spectrum with Savitzky-Golay filter and detect peaks above median."""
        import numpy as np
        from astropy import units as u
        from scipy.signal import find_peaks, savgol_filter
        from specutils import Spectrum1D
        from specutils.fitting import fit_generic_continuum

        skyValsOriginal = sky[valid].to_numpy()
        waveVals = wave.loc[valid].to_numpy()
        objectVals = objectFlux.loc[valid].to_numpy()
        skyVals = savgol_filter(sky[valid], window_length=3, polyorder=2)

        spec = Spectrum1D(flux=skyValsOriginal * u.Unit(""), spectral_axis=waveVals * u.nm)
        cont_fit = fit_generic_continuum(spec)
        continuum = cont_fit(waveVals * u.nm)
        skyVals = skyValsOriginal / continuum.value

        peaks, _ = find_peaks(skyVals, height=np.mean(skyVals) + 1 * np.std(skyVals), distance=5)

        waveVals = wave.loc[valid].to_numpy()
        objectVals = objectFlux.loc[valid].to_numpy()
        # axisBVals = axisBcoord.loc[valid].to_numpy()

        return skyValsOriginal, skyVals, peaks, waveVals, objectVals

    def _get_local_skylines_for_order(self, wmin, wmax, order, calibrationCol):
        """Fetch catalogue skylines within wavelength range and project to pixel/wavelength coordinates.

        Returns every catalogue line in range (for plotting) and the calibration lines to match against. In SOXS VIS the
        calibration lines of an order are the blend-merged, isolated centroids (``merge_blended_skylines``, DY-1645);
        everywhere else they are the lines flagged ``ISOLATED``.
        """
        import pandas as pd

        localSkylinesDF = self._project_skylines(self.skylinesDF, wmin, wmax, order)
        localSkylines = pd.to_numeric(localSkylinesDF[calibrationCol], errors="coerce").dropna().to_numpy()

        if order != "all" and self._is_soxs_vis(self.arm):
            calibrationSkylinesDF = self._project_skylines(self._merged_skylines(), wmin, wmax, order)
        elif "ISOLATED" in self.skylinesDF.columns:
            # ISOLATED FLAG MARKS LINES SUITABLE FOR CENTROIDING
            # `== True` IS KEPT: A BARE TRUTH MASK WOULD TREAT A NON-BOOLEAN ISOLATED COLUMN DIFFERENTLY
            calibrationSkylinesDF = localSkylinesDF.loc[localSkylinesDF["ISOLATED"] == True]  # noqa: E712
        else:
            return localSkylines, []
        calibrationSkylines = pd.to_numeric(calibrationSkylinesDF[calibrationCol], errors="coerce").dropna().to_numpy()

        return localSkylines, calibrationSkylines

    def _project_skylines(self, skylinesDF, wmin, wmax, order):
        """Return the catalogue lines in a wavelength range, projected to detector pixels unless ``order`` is "all"."""
        from soxspipe.commonutils import dispersion_map_to_pixel_arrays

        localSkylinesDF = skylinesDF.loc[skylinesDF["WAVELENGTH"].between(wmin, wmax)].copy()
        localSkylinesDF.rename(columns={"WAVELENGTH": "wavelength"}, inplace=True)
        localSkylinesDF["order"] = order
        localSkylinesDF["slit_position"] = 0

        if order != "all":
            localSkylinesDF = dispersion_map_to_pixel_arrays(
                log=self.log, dispersionMapPath=self.dispersionMap, orderPixelTable=localSkylinesDF
            )
        return localSkylinesDF

    def _merged_skylines(self):
        """Return the blend-merged calibration lines, built once from the catalogue without its brightness cut.

        ``self.skylinesDF`` drops the faint lines, but the merge needs them to weight a doublet and to find close
        neighbours, so the full table is loaded here. ``self.skylinesDF`` itself is left as it is.
        """
        from soxspipe.commonutils.toolkit import (
            SKYLINE_MIN_LINE_FLUX,
            get_skylines_dataframe,
            merge_blended_skylines,
        )

        if self._mergedSkylinesDF is None:
            fullDF = get_skylines_dataframe(self.log, self.settings, self.arm, minBrightnessVIS=SKYLINE_MIN_LINE_FLUX)
            self._mergedSkylinesDF = merge_blended_skylines(fullDF)
        return self._mergedSkylinesDF

    def _measure_peak_centres(self, shiftArray, skyVals, peaks):
        """Return one sub-sample centre per sky line holding the detected peaks (see ``skyline_centres``)."""
        from soxspipe.commonutils.skyline_centres import measure_line_centres

        return measure_line_centres(shiftArray, skyVals, peaks)

    def _match_peaks_to_skylines(self, peakPositions, calibrationSkylines, tolerance=15):
        """Match observed sky-line centres to catalogue skyline positions within a tolerance."""
        import numpy as np

        peakPositions = np.asarray(peakPositions, dtype=float)
        matchedSkylinePixels = []
        matchedShifts = []

        for ww in calibrationSkylines:
            if len(peakPositions) == 0:
                continue
            idx = np.where(np.abs(peakPositions - ww) < tolerance)[0]
            for i in idx:
                matchedSkylinePixels.append(ww)
                matchedShifts.append(ww - peakPositions[i])

        return matchedSkylinePixels, matchedShifts

    def _compute_clipped_shift(self, matchedSkylinePixels, matchedShifts):
        """Sigma-clip matched shift distribution and return outlier arrays, median shift and the kept shifts."""
        import numpy as np
        from astropy.stats import sigma_clip

        maskedShifts = sigma_clip(
            matchedShifts, sigma_lower=2.5, sigma_upper=2.5, maxiters=5, cenfunc="mean", stdfunc="std"
        )
        clippedWave = np.asarray(matchedSkylinePixels)[np.asarray(maskedShifts.mask)]
        clippedShifts = np.asarray(maskedShifts.data)[np.asarray(maskedShifts.mask)]
        goodShifts = np.asarray(maskedShifts.data)[~np.asarray(maskedShifts.mask)]

        medianShift = float(np.median(goodShifts)) if len(goodShifts) > 0 else 0.0
        if np.isnan(medianShift):
            medianShift = 0.0

        return clippedWave, clippedShifts, medianShift, goodShifts

    def _record_order_shift_qc(self, order, shiftNm, shiftPx=None):
        """Append the ``SKY SHIFT O<n>`` QC (nm) and, when the signed dispersion is known, ``SKY SHIFT O<n> PX``."""
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        reductionDate = utcnow_string()
        rows = [(f"SKY SHIFT O{int(order)}", round(shiftNm, 3), "nm", "Wavelength shift applied to the order")]
        if shiftPx is not None:
            rows.append(
                (
                    f"SKY SHIFT O{int(order)} PX",
                    round(shiftPx, 3),
                    "px",
                    "Same shift as detector motion along dispersion",
                )
            )
        for name, value, unit, comment in rows:
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName=name,
                qcValue=value,
                qcComment=comment,
                obsDateUtc=self.dateObs,
                reductionDateUtc=reductionDate,
                qcUnit=unit,
                toHeader=True,
                qcOrder=int(order),
            )

    def _record_line_snr_qc(self, order, lineSnr):
        """Append the ``SKY LINE SNR O<n>`` QC: the order's combined sky-line signal-to-noise (DY-1663)."""
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName=f"SKY LINE SNR O{int(order)}",
            qcValue=round(lineSnr, 1),
            qcComment="Combined sky-line signal-to-noise",
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow_string(),
            qcUnit=None,
            toHeader=True,
            qcOrder=int(order),
        )

    def _finalise_order_shifts(self, arm, uniqueOrders, shiftsNm, orderStats, dispersions):
        """*add the modelled VIS u/g shifts to the measured per-order shifts and record the sky-shift QCs*

        **Key Arguments:**

        - ``arm`` -- the arm of the instrument
        - ``uniqueOrders`` -- the orders, in the order the shifts are returned
        - ``shiftsNm`` -- ``{order: measured wavelength shift (nm)}``
        - ``orderStats`` -- per measured order, the sky lines kept and their shift scatter (nm) in the last iteration
        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``, from ``_signed_order_dispersions``

        **Return:**

        - ``orderShifts`` -- a list of per-order wavelength shifts (nm), with u/g modelled when the guards pass
        """
        from soxspipe.commonutils.ug_wavelength_transfer import nm_shift_to_detector_px

        if arm.upper() == "VIS":
            shiftsNm = self._apply_ug_transfer(shiftsNm, orderStats, dispersions)
        for o in uniqueOrders:
            shiftPx = nm_shift_to_detector_px(shiftsNm[o], dispersions[o]) if o in dispersions else None
            self._record_order_shift_qc(o, shiftsNm[o], shiftPx)
            lineSnr = orderStats.get(o, {}).get("lineSnr")
            if lineSnr is not None:
                self._record_line_snr_qc(o, lineSnr)
        return [shiftsNm[o] for o in uniqueOrders]

    def _apply_ug_transfer(self, shiftsNm, orderStats, dispersions):
        """*predict the VIS u and g shifts from the r/i sky shifts and the rotator angles (DY-1608)*

        u and g have no usable sky lines. When the r and i shifts pass the guards and the science rotator angle is
        usable, their shifts come from the ``ug-wavelength-transfer`` model; otherwise they are left unshifted and
        the reason is logged. Arms without the model's coefficients are returned unchanged and get no QC.

        **Key Arguments:**

        - ``shiftsNm`` -- ``{order: measured wavelength shift (nm)}``
        - ``orderStats`` -- per measured order, ``{"nLines": int, "scatterNm": float, "lineSnr": float}``; an r order
          whose shift came from i also carries ``"transferred": True`` and passes its guards
        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``

        **Return:**

        - ``shiftsNm`` -- a new dictionary with the u/g shifts replaced when the model applies
        """
        from soxspipe.commonutils import ug_wavelength_transfer as ug_transfer

        coeffs = ug_transfer.load_ug_transfer_coefficients(self.detectorParams)
        if coeffs is None:
            return shiftsNm

        statsPx = {
            o: {
                "nLines": orderStats[o]["nLines"],
                "lineSnr": orderStats[o].get("lineSnr"),
                "shiftPx": ug_transfer.nm_shift_to_detector_px(shiftsNm[o], dispersions[o]),
                "transferred": orderStats[o].get("transferred", False),
            }
            for o in ug_transfer.MEASURED_ORDERS
            if o in orderStats and o in dispersions
        }
        scienceAngle = ug_transfer.usable_rotator_angle(self.skySubtractedFrame.header, self.kw)
        reason = ug_transfer.ug_skip_reason(statsPx, scienceAngle)
        if reason is None:
            missing = [o for o in ug_transfer.PREDICTED_ORDERS if o in shiftsNm and o not in dispersions]
            if missing:
                reason = f"order {missing[0]}: no dispersion"
        if reason is not None:
            self.log.info(f"u/g orders left unshifted: {reason}")
            self._record_ug_source_qc("none", reason)
            return shiftsNm

        calibrationAngle = ug_transfer.calibration_rotator_angle(self._dispersion_map_header() or {}, self.kw)
        predictedPx = ug_transfer.predict_ug_detector_shifts(
            statsPx[ug_transfer.VIS_R_ORDER]["shiftPx"],
            statsPx[ug_transfer.VIS_I_ORDER]["shiftPx"],
            scienceAngle,
            calibrationAngle,
            coeffs,
        )
        updated = dict(shiftsNm)
        for o in ug_transfer.PREDICTED_ORDERS:
            if o in updated:
                updated[o] = ug_transfer.detector_px_to_nm_shift(predictedPx[o], dispersions[o])
        self._record_ug_source_qc("model", "predicted from r/i shifts and ABSROT")
        return updated

    def _record_ug_source_qc(self, source, comment):
        """Append the ``UG SHIFT SOURCE`` QC: ``model`` or ``none``, with the reason as the comment."""
        self._record_shift_source_qc("UG SHIFT SOURCE", source, comment)

    def _record_shift_source_qc(self, qcName, source, comment):
        """Append a text QC naming where a VIS shift came from, with the reason as the comment."""
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        self.qc = append_qc(
            self.qc,
            recipeName=self.recipeName,
            qcName=qcName,
            qcValue=source,
            qcComment=comment,
            obsDateUtc=self.dateObs,
            reductionDateUtc=utcnow_string(),
            qcUnit=None,
            toHeader=True,
        )

    def _signed_order_dispersions(self, extractedOrdersDF, orders):
        """*the signed dispersion d(wavelength)/d(pixel) of each order along the dispersion axis (nm/px)*

        ``pixelScaleNm`` is an absolute value, so the sign comes from the dispersion map: two wavelengths per order
        (the 25th and 75th percentiles of the order's wavelengths) are mapped to detector pixels on the slit centre.
        An order whose sample points fall off the detector, or map to one pixel, has no entry.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- the extracted orders, with ``order`` and ``wavelengthMean`` columns
        - ``orders`` -- the orders to sample

        **Return:**

        - ``dispersions`` -- ``{order: signed dispersion (nm/px)}``
        """
        import numpy as np
        import pandas as pd

        from soxspipe.commonutils import dispersion_map_to_pixel_arrays

        rows = []
        for o in orders:
            wavelengths = extractedOrdersDF.loc[extractedOrdersDF["order"] == o, "wavelengthMean"]
            if wavelengths.empty:
                continue
            for w in wavelengths.quantile([0.25, 0.75]):
                rows.append({"order": o, "wavelength": float(w), "slit_position": 0.0})
        if not rows:
            return {}

        pixels = dispersion_map_to_pixel_arrays(
            log=self.log, dispersionMapPath=self.dispersionMap, orderPixelTable=pd.DataFrame(rows)
        )
        dispersions = {}
        for o, group in pixels.groupby("order"):
            if len(group) != 2:
                continue
            group = group.sort_values("wavelength")
            deltaPixel = float(group[f"fit_{self.axisB}"].iloc[1]) - float(group[f"fit_{self.axisB}"].iloc[0])
            deltaWavelength = float(group["wavelength"].iloc[1]) - float(group["wavelength"].iloc[0])
            if deltaPixel != 0 and deltaWavelength != 0 and np.isfinite(deltaWavelength / deltaPixel):
                dispersions[o] = deltaWavelength / deltaPixel
        return dispersions

    def _dispersion_map_header(self):
        """*the dispersion map's primary header (a copy of the multi-pinhole frame header), or None if unreadable*"""
        from astropy.io import fits

        if not self.dispersionMap:
            self.log.debug("there is no dispersion map to read the multi-pinhole header from")
            return None
        try:
            return fits.getheader(self.dispersionMap, 0)
        except (OSError, KeyError, ValueError, TypeError) as e:
            self.log.debug(f"could not read the multi-pinhole header from the dispersion map ({type(e).__name__}: {e})")
            return None

    @closes_its_figures
    def _plot_skyline_shift_diagnostic(
        self,
        shiftArray,
        skyValsOriginal,
        skyVals,
        peaks,
        objectVals,
        localSkylines,
        calibrationSkylines,
        matchedSkylinePixels,
        matchedShifts,
        clippedWave,
        clippedShifts,
        medianShift,
        order,
        wmin,
        wmax,
        pixelScaleMedian,
        shiftLabel,
        iteration,
        shiftUnits,
    ):
        """Three-panel diagnostic plot: sky spectrum with skyline markers, shift scatter, object spectrum."""
        import matplotlib
        import matplotlib.pyplot as plt

        from soxspipe.commonutils.toolkit import order_label

        matplotlib.use("MacOSX")
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 9), dpi=180, gridspec_kw={"height_ratios": [2, 1, 2]})
        smoothColor = "tab:orange" if iteration == 0 else "tab:green"

        # SKY SPECTRUM WITH DETECTED PEAKS AND CATALOGUE SKYLINE POSITIONS
        # ax1.plot(shiftArray, skyValsOriginal, color="black", alpha=0.7, linewidth=0.7, label="skyFlux - original")
        ax1.plot(shiftArray, skyVals, color=smoothColor, linewidth=0.7, label="skyFlux - smoothed")
        for ww in localSkylines:
            ax1.axvline(ww, color="grey", alpha=0.25, linewidth=0.5)
        if len(localSkylines):
            ax1.plot([], [], color="grey", alpha=0.7, linewidth=0.8, label="skyline")
        for ww in calibrationSkylines:
            ax1.axvline(ww, color="blue", alpha=0.5, linewidth=0.6)
        if len(calibrationSkylines):
            ax1.plot([], [], color="blue", alpha=0.7, linewidth=0.8, label="calibration skyline")
        ax1.plot(shiftArray[peaks], skyVals[peaks], "x", color="red", label="skyFlux peaks")
        ax1.set_ylabel("sky flux ($e^{-}$)")
        ax1.legend(loc="best", fontsize=8)

        # MATCHED PEAK SHIFTS WITH SIGMA-CLIPPED OUTLIERS AND MEDIAN INDICATOR
        ax2.axhline(0, color="k", linestyle="--", linewidth=0.7)
        if matchedShifts:
            ax2.scatter(matchedSkylinePixels, matchedShifts, c="tab:green", s=30)
            ax2.scatter(clippedWave, clippedShifts, c="tab:red", s=30, marker="x")
            ax2.axhline(
                medianShift,
                color="tab:orange",
                linestyle=":",
                linewidth=0.8,
                label=f"median={medianShift:.3f} {shiftUnits} (catalogue - observed)",
            )
            ax2.legend(loc="best", fontsize=8)
        ax2.set_ylabel(f"shift ({shiftUnits})")
        ax2.set_xlim(ax1.get_xlim())

        # OBJECT SPECTRUM FOR REFERENCE
        ax3.plot(shiftArray, objectVals, color="tab:purple", linewidth=0.2, label="objectFlux")
        ax3.set_xlabel(shiftLabel)
        ax3.set_ylabel("object flux ($e^{-}$)")
        ax3.legend(loc="best", fontsize=8)
        ax3.set_xlim(ax1.get_xlim())

        if order == "all":
            fig.suptitle(
                f"\nGlobal skyline matching and shifting. WL {wmin:0.2f}-{wmax:0.2f}nm. Mean pixel {pixelScaleMedian:0.2f} nm."
            )
        else:
            fig.suptitle(
                f"\nOrder {order_label(order, self.arm)} skyline matching and shifting. "
                f"WL {wmin:0.2f}-{wmax:0.2f}nm. Mean pixel {pixelScaleMedian:0.2f} nm."
            )
        plt.show()
        plt.close(fig)

    @staticmethod
    def _orders_red_to_blue(extractedOrdersDF: "pd.DataFrame") -> list[int]:
        """*list the extracted orders from the reddest to the bluest, by their median wavelength*

        Orders whose median wavelength is NaN are left out.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- the extraction table, with ``order`` and ``wavelengthMean`` columns

        **Return:**

        - ``orders`` -- the order numbers as integers, reddest first. Empty if the table has no usable order

        **Raises:**

        - ``ValueError`` -- if an order id is not a whole number
        """
        medianWavelengths = extractedOrdersDF.groupby("order")["wavelengthMean"].median().dropna()
        for order in medianWavelengths.index:
            if not float(order).is_integer():
                raise ValueError(f"Expected a whole-number order id in the extraction table, found {order!r}.")
        return [int(o) for o in medianWavelengths.sort_values(ascending=False, kind="stable").index]

    @staticmethod
    def _resample_grid(minWavelength, maxWavelength, ratio):
        """*build the merged-spectrum wavelength grid: multiples of the output step lying inside [min, max]*

        **Key Arguments:**

        - ``minWavelength`` -- the shortest extracted wavelength (nm)
        - ``maxWavelength`` -- the longest extracted wavelength (nm)
        - ``ratio`` -- the reciprocal of the output step size

        **Return:**

        - ``grid`` -- the wavelength grid (nm), first bin >= ``minWavelength`` and last bin <= ``maxWavelength``
        """
        import numpy as np

        firstIndex = int(np.ceil(minWavelength * ratio))
        lastIndex = int(np.floor(maxWavelength * ratio))
        # GUARD AGAINST FLOATING-POINT ROUNDING PUSHING AN END BIN OUTSIDE THE DATA
        if firstIndex / ratio < minWavelength:
            firstIndex += 1
        if lastIndex / ratio > maxWavelength:
            lastIndex -= 1
        if lastIndex < firstIndex:
            raise ValueError("No output wavelength bin lies inside the extracted wavelength range.")
        return np.arange(firstIndex, lastIndex + 1) / ratio

    def merge_extracted_orders(self, extractedOrdersDF):
        """*merge the extracted order spectra in one continuous spectrum*

        **Key Arguments:**

        - ``extractedOrdersDF`` -- a data-frame containing the extracted orders

        **Return:**

        - None
        """

        self.log.debug("starting the ``merge_extracted_orders`` method")

        # THINGS TO TRY
        # - experiment with FluxConserving, Linear and Spline resampling
        # - need to dynamically define the ends of each order (possibly from order centre traces) and don't extract beyong these points
        # - set a S/N threshold, below which the data point is ignored
        # - run some kind of median smoothing to remove obvious spikes (with higher resolution than spectrograph)

        self.log.print("\n# MERGING ORDERS INTO SINGLE SPECTRUM")

        # ASTROPY HAS RESET LOGGING LEVEL -- FIX
        import logging

        import astropy.units as u
        import numpy as np
        import pandas as pd
        from astropy.nddata import VarianceUncertainty
        from specutils import Spectrum1D

        from soxspipe.commonutils.toolkit import plot_merged_spectrum_qc

        logging.getLogger().setLevel(logging.INFO + 5)

        kw = self.kw
        arm = self.arm

        # PARAMETERS FROM INPUT FILE
        # THIS IS THE STEP SIZE IN NM (0.06 nm IS SIMILAR TO XSHOOTER EXTRACTION)
        if self.arm.upper() == "NIR":
            stepWavelengthOrderMerge = 0.06
        elif self.arm.upper() == "UVB" or self.arm.upper() == "VIS":
            stepWavelengthOrderMerge = 0.02

        ratio = 1 / stepWavelengthOrderMerge
        order_list = []

        # JOIN NEIGHBOURS IN WAVELENGTH, WHATEVER THE ARM NUMBERS ITS ORDERS. EACH PAIR IS (RED ORDER, BLUE ORDER)
        orderSequence = self._orders_red_to_blue(extractedOrdersDF)
        orderPairs = list(zip(orderSequence[:-1], orderSequence[1:], strict=True))

        orderJoins = {}
        orderGaps = {}
        for red, blue in orderPairs:
            redMin = extractedOrdersDF.loc[extractedOrdersDF["order"] == red]["wavelengthMean"].min()
            blueMax = extractedOrdersDF.loc[extractedOrdersDF["order"] == blue]["wavelengthMean"].max()
            orderJoins[f"{red}-{blue}"] = (blueMax + redMin) / 2.0
            orderGaps[f"{red}-{blue}"] = blueMax - redMin

        # A FIX FOR THE UV XSHOOTER DATA (FLATS TAKEN WITH DIFFERENT LAMPS)
        # GET UNIQUE VALUES IN COLUMN
        if self.arm.upper() == "UVB":
            orderJoins["20-21"] = 363.5

        stepRatio = 20

        for red, blue in orderPairs:
            thisKey = f"{red}-{blue}"
            gap = orderGaps[thisKey]
            if gap > stepWavelengthOrderMerge * stepRatio * 2.1:
                mask = (extractedOrdersDF["order"] == red) & (
                    extractedOrdersDF["wavelengthMean"] <= orderJoins[thisKey] - stepWavelengthOrderMerge * stepRatio
                )
                extractedOrdersDF = extractedOrdersDF.loc[~mask]
                mask = (extractedOrdersDF["order"] == blue) & (
                    extractedOrdersDF["wavelengthMean"] > orderJoins[thisKey] + stepWavelengthOrderMerge * stepRatio
                )
                extractedOrdersDF = extractedOrdersDF.loc[~mask]

        if self.arm.upper() == "UVB":
            # CLIP DICHROICH REGION
            mask = extractedOrdersDF["wavelengthMean"] > 556
            extractedOrdersDF = extractedOrdersDF.loc[~mask]

        # SORT BY COLUMN NAME
        extractedOrdersDF.sort_values(["wavelengthMean"], ascending=[True], inplace=True, kind="stable")

        # DEFINE THE WAVELENGTH ARRAY
        # ENSURE THE COLUMN IS NOT EMPTY AND CONTAINS VALID NUMERIC VALUES
        if extractedOrdersDF["wavelengthMean"].notna().any():
            min_wavelength = np.min(extractedOrdersDF["wavelengthMean"])
            max_wavelength = np.max(extractedOrdersDF["wavelengthMean"])

            # ENSURE THE RANGE IS VALID
            if min_wavelength < max_wavelength and stepWavelengthOrderMerge > 0:
                # DEFINE THE WAVELENGTH ARRAY, ALIGNED TO STEP MULTIPLES AND INSIDE THE EXTRACTED RANGE
                wave_resample_grid = self._resample_grid(min_wavelength, max_wavelength, ratio)
            else:
                raise ValueError("Invalid range or step size for wavelength resampling.")
        else:
            raise ValueError("The 'wavelengthMean' column is empty or contains invalid values.")

        # ADD UNITS TO THE VARIOUS COLUMNS
        extractedOrdersDF["extractedFluxOptimal"] = extractedOrdersDF["extractedFluxOptimal"].values * u.electron
        extractedOrdersDF["extractedFluxBoxcar"] = extractedOrdersDF["extractedFluxBoxcar"].values * u.electron
        extractedOrdersDF["extractedFluxBoxcarRobust"] = (
            extractedOrdersDF["extractedFluxBoxcarRobust"].values * u.electron
        )
        extractedOrdersDF["pixelScaleNm"] = extractedOrdersDF["pixelScaleNm"].values * u.nm
        extractedOrdersDF["wavelengthMean"] = extractedOrdersDF["wavelengthMean"].values * u.nm

        flux_orig = extractedOrdersDF["extractedFluxOptimal"].values
        spectrum_orig = Spectrum1D(
            flux=flux_orig,
            spectral_axis=extractedOrdersDF["wavelengthMean"].values,
            uncertainty=VarianceUncertainty(extractedOrdersDF["varianceSpectrum"].values),
            bin_specification="center",
        )

        sky_orig = extractedOrdersDF["skyFlux"].values * u.electron
        spectrumSky_orig = Spectrum1D(
            flux=sky_orig,
            spectral_axis=extractedOrdersDF["wavelengthMean"].values,
            bin_specification="center",
        )

        # Fast linear resampling — np.interp is O(N log M) vs FluxConservingResampler's
        # O(N·M), and is accurate for output grids with similar resolution to the input.
        wave_in = spectrum_orig.spectral_axis.to_value(u.nm)
        flux_in = spectrum_orig.flux.to_value(u.electron)
        var_in = spectrum_orig.uncertainty.array
        sky_in = spectrumSky_orig.flux.to_value(u.electron)

        # Deduplicate wavelengths that may overlap at order join boundaries
        _, unique_idx = np.unique(wave_in, return_index=True)
        flux_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            flux_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)
        variance_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            var_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)
        sky_resampled = np.interp(
            wave_resample_grid,
            wave_in[unique_idx],
            sky_in[unique_idx],
            left=0.0,
            right=0.0,
        ).astype(np.float32)

        # flux_resampled = median_smooth(flux_resampled, width=3)
        merged_orders = pd.DataFrame()
        merged_orders["WAVE"] = wave_resample_grid * u.nm
        merged_orders["FLUX_COUNTS"] = flux_resampled * u.electron
        merged_orders["VARIANCE"] = variance_resampled * u.electron
        merged_orders["SKY_COUNTS"] = sky_resampled * u.electron
        merged_orders["SNR"] = merged_orders["FLUX_COUNTS"].values / np.sqrt(merged_orders["VARIANCE"].values)

        # merged_orders = calculate_rolling_snr(dataframe=merged_orders, flux_column="FLUX_COUNTS", window_size=300)

        # THE MASTER FLAT KEEPS THE PIXEL SCALE, SO A FLAT-FIELDED SPECTRUM HAS ALREADY BEEN DIVIDED BY IT. DIVIDING
        # AGAIN WOULD COUNT THE DISPERSION TWICE, SO ONLY THE NOT-FLATTENED SPECTRUM (USED FOR THE EFFICIENCY) GETS THE
        # FLUX DENSITY COLUMN, AND ONLY IT PAYS FOR COMPUTING ONE
        if self.notFlattened:
            merged_orders["FLUX_DENSITY_COUNTS"] = (
                self._add_and_resample_flux_density(
                    extractedOrdersDF, wave_in[unique_idx], unique_idx, wave_resample_grid
                )
                * u.electron
                / u.nm
            )

        if False:
            self.products, filePath = plot_merged_spectrum_qc(
                merged_orders=merged_orders,
                products=self.products,
                log=self.log,
                qcDir=self.qcDir,
                filenameTemplate=self.filenameTemplate,
                noddingSequence=self.noddingSequence,
                dateObs=self.dateObs,
                arm=self.arm,
                recipeName=self.recipeName,
                orderJoins=orderJoins,
                debug=self.debug,
            )

        return merged_orders, orderJoins

    def _add_and_resample_flux_density(self, extractedOrdersDF, uniqueWave, uniqueIdx, waveGrid):
        """*add the flux density columns to the order rows and return the optimal one resampled onto the merge grid*

        The order rows gain the ``extractedFluxDensityOptimal``, ``extractedFluxDensityBoxcar`` and
        ``extractedFluxDensityBoxcarRobust`` columns: each extracted flux divided by the pixel size in wavelength.
        The columns are added in place, as they were before DY-1719, so the not-flattened ``_EXTRACTED_ORDERS``
        table keeps them unchanged when the merge works on the caller's frame.

        **Key Arguments:**

        - ``extractedOrdersDF`` -- the extracted orders, with flux and pixel-scale quantity columns
        - ``uniqueWave`` -- the distinct order wavelengths (nm) as plain values
        - ``uniqueIdx`` -- the row index of each of those wavelengths
        - ``waveGrid`` -- the merge grid (nm)

        **Return:**

        - ``fluxDensity`` -- the resampled flux density, as float32 values
        """
        import astropy.units as u
        import numpy as np

        pixelScale = extractedOrdersDF["pixelScaleNm"].values
        for method in ("Optimal", "Boxcar", "BoxcarRobust"):
            extractedOrdersDF[f"extractedFluxDensity{method}"] = (
                extractedOrdersDF[f"extractedFlux{method}"].values / pixelScale
            )
        fluxDensityIn = extractedOrdersDF["extractedFluxDensityOptimal"].values.to_value(u.electron / u.nm)
        return np.interp(waveGrid, uniqueWave, fluxDensityIn[uniqueIdx], left=0.0, right=0.0).astype(np.float32)

    def _pinhole_offset_labels(self):
        """*the legend labels for the science frame's obs-time and rotator-angle offsets from the multi-pinhole frame*

        The multi-pinhole header is read from the primary header of the dispersion map, which is a copy of the
        multi-pinhole frame header. If the map cannot be read there are no labels and the plot is still made.

        **Return:**

        - ``labels`` -- list of legend labels, or *None* if the dispersion map header cannot be read
        """
        from soxspipe.commonutils.slit_drift_qc import pinhole_offset_labels

        pinholeHeader = self._dispersion_map_header()
        if pinholeHeader is None:
            return None
        return pinhole_offset_labels(self.skySubtractedFrame.header, pinholeHeader, self.kw)

    def plot_slit_drift_qc(self, transformer):
        """*plot the object trace's slit position against wavelength, before rectification*

        Nothing is written for the not-flattened re-extraction.

        **Key Arguments:**

        - ``transformer`` -- the image transformer holding the per-order slit-centre fits

        **Usage:**

        ```python
        optimalExtractor.plot_slit_drift_qc(transformer=transformer)
        ```

        """
        self.log.debug("starting the ``plot_slit_drift_qc`` method")

        if self.notFlattened:
            return

        from soxspipe.commonutils.slit_drift_qc import (
            plot_slit_drift_qc,
            slit_drift_arcsec,
            slit_drift_series,
            slit_zero_crossings,
        )
        from soxspipe.commonutils.toolkit import append_product, utcnow_string

        series = slit_drift_series(
            transformer.orderPixelTable,
            transformer.uniqueOrders,
            transformer.orderSlitCentreCoeffs,
            transformer.wlMinMax,
            transformer.orderSlitCentreFallback,
        )
        driftArcsec = slit_drift_arcsec(series)
        crossings = slit_zero_crossings(series)
        filename = self.filenameTemplate.replace(".fits", f"_SLIT_DRIFT_QC_PLOT{self.noddingSequence}.pdf")
        filePath = f"{self.qcDir}/{filename}"
        plot_slit_drift_qc(
            series,
            transformer.globalSlitCentreArcsec,
            f"Slit position of the object trace before rectification ({self.arm.upper()})",
            filePath,
            driftArcsec=driftArcsec,
            zeroWavelength=crossings[0] if crossings else None,
            extraLabels=self._pinhole_offset_labels(),
        )
        self._record_slit_drift_qcs(driftArcsec, crossings)

        if not isinstance(self.products, bool):
            self.products = append_product(
                self.products,
                recipeName=self.recipeName,
                productLabel=f"SLIT_DRIFT_QC_PLOT{self.noddingSequence}",
                fileName=filename,
                filePath=filePath,
                productDesc="Slit position of the object trace vs wavelength before rectification",
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow_string(),
                fileType="PDF",
                label="QC",
            )

        self.log.debug("completed the ``plot_slit_drift_qc`` method")

    def _record_slit_drift_qcs(self, driftArcsec, crossings):
        """*add the ``SLIT DRIFT`` and ``SLIT ZERO WLEN`` QCs*

        **Key Arguments:**

        - ``driftArcsec`` -- the slit-position drift in arcsec, or *None* if every order fell back to the mean
        - ``crossings`` -- the wavelengths (nm) where the trace crosses slit position 0, blue to red. Only the bluest
          is recorded
        """
        from soxspipe.commonutils.toolkit import append_qc, utcnow_string

        if driftArcsec is None:
            self.log.warning("Every order fell back to the mean slit position, so no SLIT DRIFT QC can be measured.")
        if len(crossings) > 1:
            others = ", ".join(f"{wavelength:0.1f} nm" for wavelength in crossings[1:])
            self.log.info(
                f"The object trace also crosses slit position 0 at {others}; "
                f"only the bluest crossing ({crossings[0]:0.1f} nm) is recorded."
            )
        if isinstance(self.qc, bool):
            return

        utcnow = utcnow_string()
        if driftArcsec is not None:
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName="SLIT DRIFT",
                qcValue=round(driftArcsec, 3),
                qcComment="Slit position drift of the object trace, reddest minus bluest wavelength",
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow,
                qcUnit="arcsec",
                toHeader=True,
            )
        if crossings:
            self.qc = append_qc(
                self.qc,
                recipeName=self.recipeName,
                qcName="SLIT ZERO WLEN",
                qcValue=round(crossings[0], 2),
                qcComment="Wavelength where the object trace crosses slit position 0",
                obsDateUtc=self.dateObs,
                reductionDateUtc=utcnow,
                qcUnit="nm",
                toHeader=True,
            )

    @closes_its_figures
    def plot_extracted_spectrum_qc(self, extractions):
        """*plot extracted spectrum QC plot*

        **Key Arguments:**

        - ``extractions`` -- dataframes hosting order extractions

        **Usage:**

        ```python
        optimalExtractor.plot_extracted_spectrum_qc(extractions)
        ```

        """
        self.log.debug("starting the ``plot_extracted_spectrum_qc`` method")

        from datetime import datetime

        # DO NOT PLOT IF PRODUCT TABLE HAS NOT BEEN PASSED
        # if isinstance(self.products, bool) and not self.products:
        #     return
        import matplotlib.pyplot as plt
        import pandas as pd
        from astropy.stats import sigma_clip, sigma_clipped_stats

        from soxspipe.commonutils.toolkit import add_notes_below_figure, order_label

        fig = plt.figure(figsize=(14, 12), constrained_layout=True, dpi=180)
        # THE FALLBACK SPECTRUM IS THE BOXCAR, SO ITS DUPLICATE PANEL GETS ALMOST NO HEIGHT AND IS REMOVED BELOW
        gs = fig.add_gridspec(4, 1, height_ratios=[1, 0.01, 1, 1] if self.traceFallback else None)
        toprow = fig.add_subplot(gs[0:1, :])
        secondrow = fig.add_subplot(gs[1:2, :], sharex=toprow)
        thirdrow = fig.add_subplot(gs[2:3, :], sharex=toprow)
        fourthrow = fig.add_subplot(gs[3:4, :], sharex=toprow)
        addedLegend = True

        allExtractions = pd.concat(extractions, ignore_index=True)

        # SIGMA-CLIP THE DATA TO FIND PLOT AXIS LIMITS
        arrayMask = sigma_clip(
            allExtractions["extractedFluxBoxcarRobust"],
            sigma_lower=30000,
            sigma_upper=5.0,
            maxiters=3,
            cenfunc="mean",
            stdfunc="std",
        )
        mean, median, std = sigma_clipped_stats(
            allExtractions["extractedFluxBoxcarRobust"],
            sigma=5.0,
            stdfunc="std",
            cenfunc="mean",
            maxiters=3,
        )
        maxFlux = arrayMask.max() + 0.5 * std

        for df in extractions:
            if not len(df["order"].values):
                continue
            orderLabel = order_label(df["order"].values[0], self.arm)

            extracted_wave_spectrum = df["wavelengthMean"]
            extracted_spectrum = df["extractedFluxOptimal"]
            extracted_spectrum_nonopt = df["extractedFluxBoxcarRobust"]
            extracted_variance_spectrum = df["varianceSpectrum"]
            extracted_sky_spectrum = df["skyFlux"]
            extracted_snr = df["snr"]

            try:
                if "PAE" in self.settings and self.settings["PAE"] and True:
                    line = toprow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum[10:-10],
                        zorder=2,
                        linewidth=0.1,
                    )
                    secondrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum_nonopt[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    thirdrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_snr[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    fourthrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_sky_spectrum[10:-10],
                        color=line[0].get_color(),
                        zorder=2,
                        linewidth=0.1,
                    )
                    toprow.text(
                        extracted_wave_spectrum[10:-10].mean(),
                        extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std(),
                        orderLabel,
                        fontsize=10,
                        c=line[0].get_color(),
                        verticalalignment="bottom",
                    )

                else:
                    if addedLegend:
                        label_opt = (
                            "Robust Boxcar Extraction (no object trace)" if self.traceFallback else "Optimal Extraction"
                        )
                        label_nonopt = "Robust Boxcar Extraction"
                        label_snr = "SNR"
                        label_sky = "Sky Flux"
                    else:
                        label_opt = None
                        label_nonopt = None
                        label_snr = None
                        label_sky = None
                    line = toprow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum[10:-10],
                        alpha=0.8,
                        zorder=2,
                        label=label_opt,
                        linewidth=0.5,
                    )
                    secondrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_spectrum_nonopt[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=1,
                        label=label_nonopt,
                        linewidth=0.5,
                    )
                    thirdrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_snr[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=2,
                        label=label_snr,
                        linewidth=0.5,
                    )
                    fourthrow.plot(
                        extracted_wave_spectrum[10:-10],
                        extracted_sky_spectrum[10:-10],
                        color=line[0].get_color(),
                        alpha=0.8,
                        zorder=2,
                        label=label_sky,
                        linewidth=0.5,
                    )
                    if extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std() < maxFlux:
                        toprow.text(
                            extracted_wave_spectrum[10:-10].mean(),
                            extracted_spectrum[10:-10].mean() + 1.4 * extracted_spectrum[10:-10].std(),
                            orderLabel,
                            fontsize=10,
                            c=line[0].get_color(),
                            verticalalignment="bottom",
                        )
                    addedLegend = False
            except (ValueError, IndexError, KeyError, TypeError):
                self.log.warning(f"Order skipped: {orderLabel}")

        if (
            self.skylinesDF is not None
            and not self.skylinesDF.empty
            and {"WAVELENGTH", "FLUX"}.issubset(self.skylinesDF.columns)
        ):
            xlim = fourthrow.get_xlim()
            for wavelength in self.skylinesDF["WAVELENGTH"].to_numpy():
                fourthrow.axvline(
                    wavelength,
                    ymin=0,
                    ymax=1,
                    color="grey",
                    alpha=0.2,
                    linewidth=0.5,
                    zorder=0,
                )
            fourthrow.set_xlim(xlim)

        toprow.set_title(f"Extracted Object Spectrum QC ({self.arm.upper()})", fontsize=11)
        toprow.set_ylabel("flux ($e^{-}$)", fontsize=10)
        secondrow.set_ylabel("boxcar flux ($e^{-}$)", fontsize=10)
        thirdrow.set_ylabel("SNR", fontsize=10)
        fourthrow.set_ylabel("sky flux ($e^{-}$)", fontsize=10)

        # toprow.set_yscale("log")
        # secondrow.set_yscale("log")
        fourthrow.set_xlabel("wavelength (nm)", fontsize=10)

        toprow.legend(fontsize=8, loc="best")
        secondrow.legend(fontsize=8, loc="best")
        thirdrow.legend(fontsize=8, loc="best")
        fourthrow.legend(fontsize=8, loc="best")

        toprow.set_ylim(-200, maxFlux)
        secondrow.set_ylim(-200, maxFlux)
        if self.traceFallback:
            secondrow.remove()
        # PRINT EVERY WARNING AND ERROR OF THE EXTRACTION UNDER THE PLOT
        add_notes_below_figure(fig, format_extraction_notes(self.extractionNotes))
        filename = self.filenameTemplate.replace(
            ".fits",
            f"_EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}{self.notFlattened}.pdf",
        )
        filePath = f"{self.qcDir}/{filename}"
        # plt.tight_layout()
        # import matplotlib

        # matplotlib.use("MacOSX")
        # plt.show()
        plt.savefig(filePath, dpi=120, format="pdf", bbox_inches="tight")
        plt.close(fig)
        plt.close("all")

        utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")

        if not self.notFlattened and not isinstance(self.products, bool):
            self.products = pd.concat(
                [
                    self.products,
                    pd.DataFrame(
                        [
                            {
                                "soxspipe_recipe": self.recipeName,
                                "product_label": f"EXTRACTED_ORDERS_QC_PLOT{self.noddingSequence}{self.notFlattened}",
                                "file_name": filename,
                                "file_type": "PDF",
                                "obs_date_utc": self.dateObs,
                                "reduction_date_utc": utcnow,
                                "product_desc": "QC plot of extracted source",
                                "file_path": filePath,
                                "label": "QC",
                            }
                        ]
                    ),
                ],
                ignore_index=True,
            )

        self.log.debug("completed the ``plot_extracted_spectrum_qc`` method")
        return


def extract_single_order(
    inputData,
    funclog,
    slitHalfLength,
    clippingSigma,
    clippingIterationLimit,
    globalClippingSigma,
    axisA,
    axisB,
    hornePolyOrder=3,
    debug=False,
    boxcarOnly=False,
):
    """
    *extract the object spectrum for a single order*

    **Key Arguments:**

    - ``inputData`` -- the order's cross-dispersion slices dataframe and rectified images
    - ``funclog`` -- logger
    - ``slitHalfLength`` -- half-length of the extraction aperture in detector pixels
    - ``clippingSigma`` -- sigma used to clip outliers when fitting the object profile
    - ``clippingIterationLimit`` -- iteration limit for the profile clipping
    - ``globalClippingSigma`` -- sigma used to clip the profile globally
    - ``axisA`` -- the dispersion axis
    - ``axisB`` -- the cross-dispersion axis
    - ``hornePolyOrder`` -- order of the polynomial fitted to the object profile along the dispersion axis. Default *3*
    - ``debug`` -- show the rectified images. Default *False*
    - ``boxcarOnly`` -- skip the object-profile fit and report the robust boxcar sum as the extracted spectrum. Used
      when no object trace was found. Default *False*

    **Return:**

    - ``crossDispersionSlicesDF`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    """
    # log.debug('starting the ``extract_single_order`` method')

    import numpy as np

    crossDispersionSlicesDF, orderRectifiedImages = inputData[0], inputData[1]

    # 1) SELECTING THE ORDER FROM THE ORDER PIXEL TABLE - THIS IS THE CONTINUUM OF THE OBJECT
    order = crossDispersionSlicesDF["order"].values[0]

    # SLICE SINGLE IMAGE INTO CROSS-DISPERSION SLICES
    crossDispersionSlicesDF, orderRectifiedImages = generate_masks(
        crossDispersionSlicesDF=crossDispersionSlicesDF, orderRectifiedImages=orderRectifiedImages
    )

    # RETURN IF NO SLICES WERE CREATED
    if not len(crossDispersionSlicesDF.index):
        return None

    if boxcarOnly:
        # NO OBJECT TRACE: NOTHING TO FIT A PROFILE TO, AND NO SLIT-EDGE SPILL TO MEASURE
        crossDispersionSlicesDF["slitEdgeTruncated"] = False
        extractions = compute_boxcar_extractions(
            crossDispersionSlicesDF=crossDispersionSlicesDF, orderRectifiedImages=orderRectifiedImages
        )
    else:
        extractions, capRows = _extract_optimal_order(
            crossDispersionSlicesDF=crossDispersionSlicesDF,
            orderRectifiedImages=orderRectifiedImages,
            slitHalfLength=slitHalfLength,
            clippingSigma=clippingSigma,
            clippingIterationLimit=clippingIterationLimit,
            hornePolyOrder=hornePolyOrder,
            axisB=axisB,
            order=order,
            debug=debug,
        )
        _log_profile_clip_cap(log=funclog, order=order, clippingIterationLimit=clippingIterationLimit, capRows=capRows)

    if "skyFlux" not in extractions.columns:
        extractions["skyFlux"] = np.nan

    return extractions[
        [
            "order",
            "wavelengthMean",
            "pixelScaleNm",
            "varianceSpectrum",
            "snr",
            "extractedFluxOptimal",
            "extractedFluxBoxcar",
            "extractedFluxBoxcarRobust",
            "skyFlux",
            "slitEdgeTruncated",
        ]
    ]


def _extract_optimal_order(
    crossDispersionSlicesDF,
    orderRectifiedImages,
    slitHalfLength,
    clippingSigma,
    clippingIterationLimit,
    hornePolyOrder,
    axisB,
    order,
    debug,
):
    """*fit the object profile of one order and extract it optimally*

    **Return:**

    - ``extractions`` -- the order's slices dataframe with the optimal and boxcar extraction columns added
    - ``capRows`` -- a ``profile_clip_cap_rows`` count of the slit rows that stopped at the profile clipping iteration
      limit while still clipping. The caller logs it; this function logs nothing
    """
    import matplotlib.pyplot as plt

    # 2) DETERMINING LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS - FITTING OF THE FRACTIONAL FLUX
    # ITERATE FIRST PIXEL IN EACH SLICE AND THEN MOVE TO NEXT
    crossDispersionSlicesDF, orderRectifiedImages, capRows = fit_object_profile(
        crossDispersionSlicesDF=crossDispersionSlicesDF,
        orderRectifiedImages=orderRectifiedImages,
        slitHalfLength=slitHalfLength,
        clippingSigma=clippingSigma,
        clippingIterationLimit=clippingIterationLimit,
        hornePolyOrder=hornePolyOrder,
        axisB=axisB,
        order=order,
        debug=debug,
        plt=plt,
    )

    # FLAG THE WAVELENGTH SLICES WHERE THE OBJECT IS CUT OFF BY THE SLIT EDGE
    # (MUST BE SET BEFORE THE SLICES ARE SORTED AND FILTERED)
    crossDispersionSlicesDF["slitEdgeTruncated"] = _slit_edge_truncated_columns(
        orderRectifiedImages["fluxRaw"], orderRectifiedImages["mask"], orderRectifiedImages["unsupported"]
    )

    # PLOT THE RECTIFIED IMAGES
    if debug:
        plot_rectified_images(orderRectifiedImages=orderRectifiedImages, order=order)

    extractions = compute_extractions(
        crossDispersionSlicesDF=crossDispersionSlicesDF, orderRectifiedImages=orderRectifiedImages, order=order
    )
    return extractions, capRows


def _log_profile_clip_cap(log, order: int, clippingIterationLimit: int, capRows: profile_clip_cap_rows) -> None:
    """*log, once for the order, how many slit rows stopped at the profile clipping iteration limit*

    **Key Arguments:**

    - ``log`` -- logger
    - ``order`` -- the order number
    - ``clippingIterationLimit`` -- iteration limit for the profile clipping
    - ``capRows`` -- counts of slit rows that hit the limit while still clipping, as ``refitted`` and ``kept``

    Logs nothing when no row hit the limit.
    """
    totalRows = capRows.refitted + capRows.kept
    if not totalRows:
        return
    log.info(
        f"fit_object_profile: order {order}: iteration limit of {clippingIterationLimit} reached with points still "
        f"being clipped in {totalRows} slit rows; {capRows.refitted} refitted to the surviving points, "
        f"{capRows.kept} kept their last fit because too few points survived"
    )


def compute_extractions(crossDispersionSlicesDF, orderRectifiedImages, order):
    """Compute optimal and boxcar extracted spectra from cross-dispersion slices."""
    import numpy as np

    # WORK ON A COPY SO THE CALLER'S SLICES FRAME IS NEVER MODIFIED
    crossDispersionSlicesDF = crossDispersionSlicesDF.copy()

    # CALCULATE HORNE 86 NUMERATOR (EQU 8)
    orderRectifiedImages["horneNumerator"] = np.ma.masked_array(
        orderRectifiedImages["fluxRaw"] * orderRectifiedImages["objectProfile"] / orderRectifiedImages["variance"],
        mask=orderRectifiedImages["mask"],
    )
    horneNumeratorSum = orderRectifiedImages["horneNumerator"].T.sum(axis=1)

    # CALCULATE HORNE 86 DENOMINATOR (EQU 8)
    orderRectifiedImages["horneDenominator"] = np.ma.masked_array(
        np.power(orderRectifiedImages["objectProfile"], 2) / orderRectifiedImages["variance"],
        mask=orderRectifiedImages["mask"],
    )
    horneDenominatorSum = orderRectifiedImages["horneDenominator"].T.sum(axis=1)

    orderRectifiedImages["optimalExtraction"] = np.ma.masked_array(
        orderRectifiedImages["horneNumerator"] / orderRectifiedImages["horneDenominator"],
        mask=orderRectifiedImages["mask"],
    )

    # plot_rectified_images(orderRectifiedImages=orderRectifiedImages, order=order)

    wavelengthMasked = np.ma.masked_array(orderRectifiedImages["wavelength"], mask=orderRectifiedImages["mask"])
    crossDispersionSlicesDF["wavelengthMean"] = np.ma.mean(wavelengthMasked.T, axis=1)

    # CALCULATE THE FINAL EXTRACTED SPECTRA
    crossDispersionSlicesDF["varianceSpectrum"] = 1 / horneDenominatorSum
    crossDispersionSlicesDF["extractedFluxOptimal"] = horneNumeratorSum / horneDenominatorSum
    crossDispersionSlicesDF["extractedFluxBoxcar"] = orderRectifiedImages["fluxRaw"].T.sum(axis=1)
    if "fluxSky" in orderRectifiedImages.keys():
        crossDispersionSlicesDF["skyFlux"] = orderRectifiedImages["fluxSky"].T.mean(axis=1)
    crossDispersionSlicesDF["extractedFluxBoxcarRobust"] = (
        np.ma.masked_array(orderRectifiedImages["fluxRaw"].T, mask=orderRectifiedImages["mask"].T)
        .sum(axis=1)
        .astype(float)
    )
    crossDispersionSlicesDF["snr"] = crossDispersionSlicesDF["extractedFluxOptimal"] / np.power(
        crossDispersionSlicesDF["varianceSpectrum"], 0.5
    )

    return _sort_and_clean_extractions(crossDispersionSlicesDF)


def compute_boxcar_extractions(crossDispersionSlicesDF, orderRectifiedImages):
    """*compute the robust boxcar spectrum of an order, reported in the optimal-extraction columns*

    The boxcar sums the unmasked rectified flux across the slit. It needs no object profile, so it is used when no
    object trace was found.

    **Key Arguments:**

    - ``crossDispersionSlicesDF`` -- dataframe of the order's cross-dispersion slices (with pixel scale)
    - ``orderRectifiedImages`` -- the rectified images of the order, including the final ``mask`` and ``variance``

    **Return:**

    - ``crossDispersionSlicesDF`` -- a filtered copy of the slices dataframe with the extracted columns added, sorted by
      wavelength. The input dataframe is not modified.
      ``extractedFluxOptimal`` is a copy of ``extractedFluxBoxcarRobust``, and ``varianceSpectrum`` is the sum of the
      variance over the same unmasked pixels
    """
    import numpy as np

    # WORK ON A COPY SO THE CALLER'S SLICES FRAME IS NEVER MODIFIED
    crossDispersionSlicesDF = crossDispersionSlicesDF.copy()
    mask = orderRectifiedImages["mask"]

    wavelengthMasked = np.ma.masked_array(orderRectifiedImages["wavelength"], mask=mask)
    crossDispersionSlicesDF["wavelengthMean"] = np.ma.mean(wavelengthMasked.T, axis=1).filled(np.nan)

    crossDispersionSlicesDF["extractedFluxBoxcar"] = orderRectifiedImages["fluxRaw"].T.sum(axis=1)
    if "fluxSky" in orderRectifiedImages:
        crossDispersionSlicesDF["skyFlux"] = orderRectifiedImages["fluxSky"].T.mean(axis=1)
    crossDispersionSlicesDF["extractedFluxBoxcarRobust"] = (
        np.ma.masked_array(orderRectifiedImages["fluxRaw"].T, mask=mask.T).sum(axis=1).astype(float).filled(np.nan)
    )
    crossDispersionSlicesDF["varianceSpectrum"] = (
        np.ma.masked_array(orderRectifiedImages["variance"].T, mask=mask.T).sum(axis=1).astype(float).filled(np.nan)
    )
    crossDispersionSlicesDF["extractedFluxOptimal"] = crossDispersionSlicesDF["extractedFluxBoxcarRobust"].copy()
    crossDispersionSlicesDF["snr"] = crossDispersionSlicesDF["extractedFluxOptimal"] / np.sqrt(
        crossDispersionSlicesDF["varianceSpectrum"]
    )

    return _sort_and_clean_extractions(crossDispersionSlicesDF)


def _sort_and_clean_extractions(crossDispersionSlicesDF):
    """*sort the extracted slices by wavelength and drop the unusable ones*

    **Key Arguments:**

    - ``crossDispersionSlicesDF`` -- the extracted slices. This dataframe is not modified

    **Return:**

    - ``crossDispersionSlicesDF`` -- a filtered copy, without zero-wavelength or incomplete slices
    """
    # SORT BY COLUMN NAME
    crossDispersionSlicesDF = crossDispersionSlicesDF.sort_values(["wavelengthMean"], ascending=[True], kind="stable")

    # REMOVE 0 WAVELENGTH (A LOCAL COPY, SO THE IN-PLACE dropna BELOW NEVER TARGETS A SLICE)
    crossDispersionSlicesDF = crossDispersionSlicesDF.loc[crossDispersionSlicesDF["wavelengthMean"] > 0].copy()

    crossDispersionSlicesDF.dropna(
        how="any",
        subset=[
            "pixelScaleNm",
            "wavelengthMean",
            "extractedFluxOptimal",
            "snr",
            "varianceSpectrum",
            "extractedFluxBoxcarRobust",
        ],
        inplace=True,
    )

    return crossDispersionSlicesDF


@closes_its_figures
def plot_rectified_images(orderRectifiedImages, order):
    """Plot available rectified order image layers for debug inspection."""

    import matplotlib.pyplot as plt
    from astropy.stats import sigma_clipped_stats

    if not isinstance(orderRectifiedImages, dict) or len(orderRectifiedImages) == 0:
        return

    import matplotlib

    matplotlib.use("MacOSX")

    for key, value in orderRectifiedImages.items():
        mean, median, std = sigma_clipped_stats(value, sigma=5.0, stdfunc="std", cenfunc="mean", maxiters=3)
        fig = plt.figure(
            num=None,
            figsize=(135, 1),
            dpi=None,
            facecolor=None,
            edgecolor=None,
            frameon=True,
        )
        fig.suptitle(f"{key} for Order {order}", fontsize=16)
        if "mask" in key.lower():
            plt.imshow(value, interpolation="none", aspect="auto", vmin=0, vmax=1, cmap="viridis")
        else:
            plt.imshow(
                value, interpolation="none", aspect="auto", vmin=mean - 2 * std, vmax=mean + 2 * std, cmap="viridis"
            )
        plt.show()

    return


def generate_masks(crossDispersionSlicesDF, orderRectifiedImages):
    """This function is used to create masks for the cross-dispersion slices and to calculate the pixel scale in wavelength space.

    **Key Arguments:**

    - ``crossDispersionSlicesDF`` -- the order dataframe
    - ``orderRectifiedImages`` -- the rectified images for the order

    **Return:**

    - ``crossDispersionSlicesDF`` -- dataframe containing metadata for each cross-dispersion slice (single data-points in extracted spectrum)
    - ``orderRectifiedImages`` -- the rectified images for the order with updated masks
    """

    import numpy as np
    from astropy.stats import sigma_clip

    # UNPACK RETIFIED ORDER IMAGES
    bpMask = orderRectifiedImages["bpMask"]
    fluxRaw = orderRectifiedImages["fluxRaw"]
    wavelength = orderRectifiedImages["wavelength"]

    # CALCULATE THE PIXEL SCALE BEFORE ANY CLIPPING OCCURS
    crossDispersionSlicesDF["pixelScaleNm"] = np.ma.mean(wavelength, axis=0)
    this = (
        crossDispersionSlicesDF["pixelScaleNm"].values[2:] - crossDispersionSlicesDF["pixelScaleNm"].values[:-2]
    ) / 2
    this = np.insert(this, 0, np.nan)
    this = np.append(this, np.nan)
    crossDispersionSlicesDF["pixelScaleNm"] = np.abs(this)

    ## REMOVE BAD PIXELS AND COSMIC RAYS FROM THE FLUX ARRAY
    fluxRawMask = _sigma_clip_and_mask(fluxRaw, bpMask)

    # FAIL SAFE FOR BAD WAVELENGTH VALUES - SOME ODD WAVELENGTHS FROM DISPERSION SOLUTION CAN CAUSE PROBLEMS WITH PROFILE FITTING
    wlMask = wavelength.copy()
    wlMask[wlMask > 0] = 3
    wlMask[wlMask < 0.1] = 1
    wlMask[wlMask > 2] = 0
    wavelengthMasked = np.ma.array(wavelength, mask=wlMask)
    wavelengthMasked = sigma_clip(
        wavelengthMasked,
        sigma_lower=1,
        sigma_upper=1,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )

    orderRectifiedImages["mask"] = fluxRawMask | (wlMask.astype(bool))

    # IF THERE IS MORE THAN 1 PIXEL MASKED IN THE CROSS-DISPERSION DIRECTION, THEN FLAG THE ENTIRE COLUMN AS BAD
    fullColumnMaskFlags = np.sum(orderRectifiedImages["mask"], axis=1)
    # fullColumnMaskFlags = fullColumnMaskFlags > 1
    # orderRectifiedImages["mask"][fullColumnMaskFlags] = True
    crossDispersionSlicesDF["mask"] = [x for x in orderRectifiedImages["mask"].T]

    return crossDispersionSlicesDF, orderRectifiedImages


def fit_object_profile(
    crossDispersionSlicesDF,
    orderRectifiedImages,
    slitHalfLength,
    clippingSigma,
    clippingIterationLimit,
    hornePolyOrder,
    axisB,
    order,
    debug,
    plt,
):
    """*fit the object profile of one order along the dispersion axis, row by row*

    **Key Arguments:**

    - ``crossDispersionSlicesDF`` -- the order's cross-dispersion slices dataframe. Gains an ``objectProfile`` column
    - ``orderRectifiedImages`` -- the order's rectified images. Gains the ``objectProfile`` and ``unsupported`` layers
    - ``slitHalfLength`` -- half-length of the extraction aperture in detector pixels
    - ``clippingSigma`` -- sigma used to clip outliers when fitting the object profile
    - ``clippingIterationLimit`` -- iteration limit for the profile clipping
    - ``hornePolyOrder`` -- order of the polynomial fitted to each slit row along the dispersion axis
    - ``axisB`` -- the cross-dispersion axis
    - ``order`` -- the order number
    - ``debug`` -- plot each row's fit. Default *False*
    - ``plt`` -- matplotlib pyplot module

    **Return:**

    - ``crossDispersionSlicesDF`` -- the slices dataframe with the object profile added
    - ``orderRectifiedImages`` -- the rectified images with the object profile added
    - ``capRows`` -- a ``profile_clip_cap_rows`` count of the slit rows that stopped at ``clippingIterationLimit`` while
      still clipping: ``refitted`` to the surviving points, or ``kept`` at their last fit because too few points
      survived. The function itself logs nothing, since a message per row would fire thousands of times
    """
    import numpy as np
    from astropy.stats import sigma_clip

    crossSlitProfiles = []
    # SLIT ROWS THAT STOPPED AT THE CLIPPING ITERATION LIMIT WHILE STILL CLIPPING
    refittedRows = 0
    keptRows = 0

    shape = orderRectifiedImages["fluxRaw"].shape
    ww = shape[1]
    ss = shape[0]

    fluxRawMasked = np.ma.masked_array(orderRectifiedImages["fluxRaw"], mask=orderRectifiedImages["mask"])
    # RETURN THE SUM OF THE ARRAY ELEMENTS OVER THE GIVEN AXIS. MASKED ELEMENTS ARE SET TO 0 INTERNALLY.
    fluxRawMaskedSum = fluxRawMasked.sum(axis=0)

    ## THIS IS THE NORMALISED FLUX USED FOR FITTING THE DISPERSION PROFILES - THIS IS THE FRACTIONAL FLUX IN HORNE 1986 PAPER
    fluxRawNormalisedMasked = fluxRawMasked / fluxRawMaskedSum[np.newaxis, :]

    dispersionAxisPixelsOrignal = range(0, ww)

    # DETERMINE LOW-ORDER POLYNOMIALS FOR FITTING THE PROFILE ALONG THE WAVELENGTH AXIS
    for slitPixelIndex in range(0, ss):
        iteration = 0
        clipped_count = 1

        fractions = fluxRawNormalisedMasked[slitPixelIndex, :]
        dispersionAxisPixels = np.ma.masked_array(range(0, ww), mask=fractions.mask)

        mask = orderRectifiedImages["mask"][slitPixelIndex, :]

        # fractions MAY STILL CONTAIN BAD-PIXEL/CRHs SO DROP PIXELS MASKED IN STEP 1 ABOVE
        a = [fractions, dispersionAxisPixels]
        fractions, dispersionAxisPixels = [np.ma.compressed(np.ma.masked_array(i, mask)) for i in a]

        coeff = []

        # STOP AT THE ITERATION LIMIT, OR AS SOON AS A PASS CLIPS NOTHING (A FURTHER PASS WOULD REPEAT THE SAME FIT)
        while (iteration < clippingIterationLimit) and (clipped_count > 0):
            # FIT A POLYNOMIAL TO THE FRACTIONAL FLUXES
            if len(dispersionAxisPixels):
                coeff = np.polyfit(dispersionAxisPixels, fractions, deg=hornePolyOrder)
            else:
                coeff = []
                break
            residuals = fractions - np.polyval(coeff, dispersionAxisPixels)
            passCount = len(fractions)

            # REMOVE REMAINING OUTLIERS
            masked_residuals = sigma_clip(
                residuals,
                sigma_lower=clippingSigma,
                sigma_upper=clippingSigma,
                maxiters=1,
                cenfunc="mean",
                stdfunc="std",
            )
            # REDUCE ARRAYS TO NON-MASKED VALUES
            a = [fractions, dispersionAxisPixels]
            fractions, dispersionAxisPixels = [
                np.ma.compressed(np.ma.masked_array(i, masked_residuals.mask)) for i in a
            ]
            clipped_count = passCount - len(fractions)
            iteration = iteration + 1

        # A CAP EXIT WITH POINTS STILL BEING CLIPPED LEAVES COEFFICIENTS THAT INCLUDE THE POINTS JUST REJECTED.
        # REFIT ONLY WHEN THE SURVIVORS CONSTRAIN EVERY COEFFICIENT; OTHERWISE KEEP THE LAST FIT
        if len(coeff) and clipped_count > 0:
            if len(dispersionAxisPixels) > hornePolyOrder:
                coeff = np.polyfit(dispersionAxisPixels, fractions, deg=hornePolyOrder)
                refittedRows += 1
            else:
                keptRows += 1

        # GENERATE THE FINAL FITTING PROFILE FOR THIS SLIT POSITION
        if len(coeff):
            profile = np.polyval(coeff, dispersionAxisPixelsOrignal)
            profile[profile < 0] = 0
        else:
            profile = np.zeros_like(dispersionAxisPixelsOrignal)
        crossSlitProfiles.append(profile)

        if debug:
            import matplotlib
            import matplotlib.pyplot as plt

            matplotlib.use("MacOSX")
            plt.scatter(dispersionAxisPixels, fractions, alpha=0.2)
            plt.plot(dispersionAxisPixelsOrignal, profile, color="red")
            plt.title(f"Fitted Profile for Order {order}")
            plt.ylim([-1, 1])
            plt.show()

    crossSlitProfiles = np.array(crossSlitProfiles)
    transposedProfiles = crossSlitProfiles.tolist()
    crossDispersionProfile = np.array([np.array(t) for t in transposedProfiles])

    # OFF-ORDER PIXELS (NON-FINITE FLUX) CANNOT HOLD OBJECT FLUX, BUT THEIR EXTRAPOLATED
    # ROW POLYNOMIALS WOULD TAKE PROFILE WEIGHT FROM THE ON-ORDER ROWS.
    # ISOLATED MASKED PIXELS (CRHs, BAD PIXELS) IN WELL-SAMPLED ROWS KEEP THEIR WEIGHT
    # SO HORNE CORRECTS FOR THEIR MISSING FLUX. MASKED PIXELS IN ROWS WITH SPARSE LOCAL
    # SUPPORT (E.G. OBJECT ON THE SLIT EDGE) GET NO WEIGHT: THEIR ROW POLYNOMIAL IS AN
    # UNCONSTRAINED EXTRAPOLATION THAT WOULD INFLATE THE FLUX RESCALING.
    # THIS ASSUMES ON-ORDER BAD PIXELS ARRIVE FINITE AND FLAGGED IN THE MASK, NOT AS NaN
    unsupported = _unsupported_pixel_mask(orderRectifiedImages["fluxRaw"], orderRectifiedImages["mask"])
    crossDispersionProfile[unsupported] = 0
    orderRectifiedImages["unsupported"] = unsupported

    crossDispersionProfileSums = np.array([x.sum() for x in crossDispersionProfile.T])
    # A COLUMN WITH NO ON-ORDER WEIGHT NORMALISES TO NaN AND IS DROPPED BY compute_extractions
    with np.errstate(divide="ignore", invalid="ignore"):
        orderRectifiedImages["objectProfile"] = crossDispersionProfile / crossDispersionProfileSums[np.newaxis :]
    crossDispersionSlicesDF["objectProfile"] = [x for x in orderRectifiedImages["objectProfile"].T]

    return crossDispersionSlicesDF, orderRectifiedImages, profile_clip_cap_rows(refittedRows, keptRows)


def _unsupported_pixel_mask(fluxRaw, mask):
    """*flag pixels that cannot hold object flux or whose row profile is an unconstrained extrapolation*

    **Key Arguments:**

    - ``fluxRaw`` -- the rectified flux array (rows are slit positions, columns are wavelength). Off-order pixels
      are non-finite
    - ``mask`` -- the boolean bad-pixel mask, True where a pixel is masked

    **Return:**

    - ``unsupported`` -- boolean array, True where the pixel is non-finite or is masked in a slit row with sparse
      local usable pixels
    """
    import numpy as np
    from scipy.ndimage import uniform_filter1d

    isFinite = np.isfinite(fluxRaw)
    usable = ~mask & isFinite
    # FRACTION OF USABLE PIXELS IN A WINDOW ALONG THE DISPERSION AXIS, MEASURED IN EACH SLIT ROW
    localSupport = uniform_filter1d(usable.astype(float), size=PROFILE_SUPPORT_WINDOW, axis=1, mode="nearest")
    return ~isFinite | (mask & (localSupport < PROFILE_MIN_LOCAL_SUPPORT))


def _slit_edge_truncated_columns(fluxRaw, mask, unsupported):
    """*flag wavelength columns where the object is cut off by the slit edge*

    A column is truncated when the unsupported rows touch the top or bottom of the window and the nearest supported
    row still holds a large share of the column's peak flux, i.e. the object continues beyond the usable slit.

    **Key Arguments:**

    - ``fluxRaw`` -- the rectified flux array (rows are slit positions, columns are wavelength)
    - ``mask`` -- the boolean bad-pixel mask, True where a pixel is masked
    - ``unsupported`` -- boolean array from ``_unsupported_pixel_mask``

    **Return:**

    - ``truncated`` -- boolean array with one value per wavelength column
    """
    import numpy as np

    rowCount = fluxRaw.shape[0]
    columnIndex = np.arange(fluxRaw.shape[1])
    usable = ~mask & np.isfinite(fluxRaw)
    usableFlux = np.where(usable, fluxRaw, np.nan)
    peak = np.max(np.where(usable, fluxRaw, -np.inf), axis=0)

    # NEAREST ROW WITH USABLE FLUX TO THE TOP AND BOTTOM EDGES OF THE WINDOW (AN ISOLATED MASKED
    # PIXEL IN A WELL-SAMPLED ROW IS SUPPORTED BUT HAS NO FLUX TO COMPARE, SO LOOK PAST IT)
    supported = usable & ~unsupported
    hasSupport = supported.any(axis=0)
    firstSupported = np.argmax(supported, axis=0)
    lastSupported = rowCount - 1 - np.argmax(supported[::-1], axis=0)

    threshold = SLIT_EDGE_PEAK_FRACTION * peak
    with np.errstate(invalid="ignore"):
        topSpill = unsupported[0] & (usableFlux[firstSupported, columnIndex] >= threshold)
        bottomSpill = unsupported[-1] & (usableFlux[lastSupported, columnIndex] >= threshold)
    return hasSupport & (peak > 0) & (topSpill | bottomSpill)


def _sigma_clip_and_mask(fluxRaw, bpMask):
    import numpy as np
    from astropy.stats import sigma_clip

    # BUILDS A 2D DISTANCE-FROM-CENTER WEIGHTING ARRAY (SAME SHAPE AS fluxRaw) .. ONLY USED FOR SIGMA-CLIPPING THE FLUX ARRAY, NOT FOR WEIGHTING THE PROFILE FITTING
    # EACH ROW GETS IDENTICAL CENTER-DISTANCE WEIGHTS, LATER USED TO BIAS SIGMA-CLIPPING SO OUTLIERS FARTHER FROM THE CENTER ARE PENALIZED MORE.
    sliceIndexArray = np.tile(
        np.abs(np.arange(fluxRaw.shape[1]) - fluxRaw.shape[1] // 2) * 50 + 1,
        (fluxRaw.shape[0], 1),
    )
    bpMask[bpMask > 1] = 1
    fluxRawMasked = np.ma.array(fluxRaw, mask=np.isnan(fluxRaw) | bpMask)
    fluxRawMasked = sigma_clip(
        fluxRawMasked * sliceIndexArray,
        sigma_lower=3000,
        sigma_upper=7,
        maxiters=3,
        cenfunc="mean",
        stdfunc="std",
        axis=1,
    )

    ## COMPUTE NUMBER OF NEWLY MASKED PIXELS IN EACH ROW AFTER SIGMA-CLIPPING
    newlyMaskedPixels = np.zeros(fluxRaw.shape[0])
    for i in range(fluxRaw.shape[0]):
        newlyMaskedPixels[i] = np.sum(fluxRawMasked.mask[i]) - np.sum(bpMask[i])
        # print(f"Row {i}: {np.sum(fluxRawMasked.mask[i])} pixels masked, {np.sum(bpMask[i])} pixels were already masked, {newlyMaskedPixels[i]} newly masked")
    # print(f"Total newly masked pixels across all slices: {int(np.sum(newlyMaskedPixels))}")

    return fluxRawMasked.mask
