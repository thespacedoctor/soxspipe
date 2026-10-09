#!/usr/bin/env python
"""
*The recipe for creating master-bias frames *

Author
: David Young & Marco Landoni

Date Created
: January 22, 2020
"""

################# GLOBAL IMPORTS ####################

import os
import sys

from soxspipe.commonutils.toolkit import generic_quality_checks, utcnow_string

from .base_recipe import base_recipe

os.environ["TERM"] = "vt100"

BIAS_HISTOGRAM_HALF_WIDTH_SIGMA = 5
BIAS_HISTOGRAM_BINS = 100
# HALF-WIDTH (e-) OF THE HISTOGRAM RANGE WHEN EVERY PLOTTED PIXEL HAS THE SAME VALUE
BIAS_HISTOGRAM_FLAT_HALF_WIDTH = 0.5
RAW_SERIES_COLOUR = "C0"
MASTER_SERIES_COLOUR = "C1"
# OPACITY OF THE FILLED HISTOGRAMS, LOW ENOUGH FOR THE OVERLAP OF THE TWO SERIES TO SHOW
HISTOGRAM_FILL_ALPHA = 0.4
STATS_TABLE_FONT_SIZE = 7
STATS_TABLE_EDGE_COLOUR = "0.8"
STATS_TABLE_EDGE_WIDTH = 0.4
# LEFT PADDING OF THE LABEL COLUMN, AS A FRACTION OF THE CELL WIDTH
STATS_TABLE_LABEL_PAD = 0.04
# TABLE POSITION IN AXES COORDINATES (LEFT, BOTTOM, WIDTH, HEIGHT): BELOW THE AXES AND THE X-AXIS LABEL
STATS_TABLE_BBOX = (0.0, -0.45, 1.0, 0.28)
STATS_TABLE_COLUMNS = ("", "mean (e-)", "median (e-)", "RON (e-)")
STATS_TABLE_NOT_APPLICABLE = "\u2014"


def _finite_unmasked_pixels(data, *masks):
    """*return the finite pixel values that no mask flags*

    **Key Arguments:**

    - ``data`` -- the pixel array
    - ``masks`` -- any number of boolean masks the same shape as `data` (`True` = exclude). `None` masks are skipped

    **Return:**

    - ``pixels`` -- 1D array of the finite, unmasked pixel values

    **Usage:**

    ```python
    pixels = _finite_unmasked_pixels(frame.data, frame.mask, badPixelMask)
    ```
    """
    import numpy as np

    data = np.asarray(data)
    keep = np.isfinite(data)
    for mask in masks:
        if mask is not None:
            keep &= ~np.asarray(mask, dtype=bool)

    return data[keep]


def _histogram_bin_edges(rawPixels, masterPixels, centre, frameRon):
    """*return the bias histogram bin edges: the centre plus or minus five frame RONs*

    If the frame RON is zero or not finite, the edges span the range of all plotted pixels instead.

    **Key Arguments:**

    - ``rawPixels`` -- the plotted raw-frame pixel values
    - ``masterPixels`` -- the plotted master-bias pixel values
    - ``centre`` -- the bias level to centre the histogram on
    - ``frameRon`` -- the RON of the plotted raw frame

    **Return:**

    - ``binEdges`` -- `BIAS_HISTOGRAM_BINS + 1` strictly increasing bin edges
    """
    import numpy as np

    if np.isfinite(frameRon) and frameRon > 0:
        low = centre - BIAS_HISTOGRAM_HALF_WIDTH_SIGMA * frameRon
        high = centre + BIAS_HISTOGRAM_HALF_WIDTH_SIGMA * frameRon
    else:
        allPixels = np.concatenate([rawPixels, masterPixels])
        low, high = float(np.min(allPixels)), float(np.max(allPixels))
        if low == high:
            low, high = low - BIAS_HISTOGRAM_FLAT_HALF_WIDTH, high + BIAS_HISTOGRAM_FLAT_HALF_WIDTH

    return np.linspace(low, high, BIAS_HISTOGRAM_BINS + 1)


def _unclipped_raw_pixels(noiseFrame, badPixelMask, meanLevel):
    """*return the flux of the raw pixels that were neither sigma-clipped nor flagged bad*

    **Key Arguments:**

    - ``noiseFrame`` -- the mean-subtracted raw frame, whose mask is the sigma-clip mask
    - ``badPixelMask`` -- the bad-pixel mask the raw frame had before it was sigma-clipped. `None` if it had none
    - ``meanLevel`` -- the mean bias level that was subtracted from the frame

    **Return:**

    - ``pixels`` -- the flux (e-) of the surviving finite pixels with the mean bias level added back

    **Usage:**

    ```python
    pixels = _unclipped_raw_pixels(noiseFrame, badPixelMask, meanLevel)
    ```
    """
    return _finite_unmasked_pixels(noiseFrame.data, noiseFrame.mask, badPixelMask) + meanLevel


def bias_distribution_summary(
    rawPixels,
    frameRon,
    masterPixels,
    rawRon,
    masterRon,
):
    """*summarise the raw and master bias pixel distributions for the QC plot*

    **Key Arguments:**

    - ``rawPixels`` -- the flux (e-) of the unclipped pixels of the plotted raw frame
    - ``frameRon`` -- the sigma-clipped standard deviation (e-) of the plotted raw frame
    - ``masterPixels`` -- the flux (e-) of the unmasked master-bias pixels
    - ``rawRon`` -- the mean of the per-frame RON values (the RAW RON QC)
    - ``masterRon`` -- the standard deviation (e-) of the stacked, mean-subtracted frame (the MASTER RON QC)

    **Return:**

    - ``summary`` -- dictionary with the keys `rawMean`, `rawMedian`, `masterMean`, `masterMedian`, `binEdges`
      and `statsTable`. `statsTable` is a dictionary with `columns` (the four column headers) and `rows` (one list
      of cell strings per row: the plotted raw frame, the master bias and the RAW RON QC). Numbers are formatted
      to two decimal places. The mean and median cells of the RAW RON QC row hold a dash because that value is
      a mean of per-frame RONs and not a pixel statistic. `masterRon` is carried by the master bias row

    Raises `ValueError` if either sample has no finite, unmasked pixels. A bias frame with every pixel masked is
    unusable, so the recipe stops rather than writing a plot of NaN values.

    **Usage:**

    ```python
    summary = bias_distribution_summary(rawPixels, frameRon, masterPixels, rawRon, masterRon)
    ```
    """
    import numpy as np

    for name, pixels in (("raw frame", rawPixels), ("master bias", masterPixels)):
        if len(pixels) == 0:
            raise ValueError(f"The {name} has no finite, unmasked pixels to plot in the bias distribution QC plot")

    rawMean = float(np.mean(rawPixels))
    rawMedian = float(np.median(rawPixels))
    masterMean = float(np.mean(masterPixels))
    masterMedian = float(np.median(masterPixels))
    binEdges = _histogram_bin_edges(rawPixels, masterPixels, rawMean, frameRon)

    statsTable = {
        "columns": list(STATS_TABLE_COLUMNS),
        "rows": [
            ["raw frame", f"{rawMean:.2f}", f"{rawMedian:.2f}", f"{frameRon:.2f}"],
            ["master bias", f"{masterMean:.2f}", f"{masterMedian:.2f}", f"{masterRon:.2f}"],
            ["RAW RON (QC)", STATS_TABLE_NOT_APPLICABLE, STATS_TABLE_NOT_APPLICABLE, f"{rawRon:.2f}"],
        ],
    }

    return {
        "rawMean": rawMean,
        "rawMedian": rawMedian,
        "masterMean": masterMean,
        "masterMedian": masterMedian,
        "binEdges": binEdges,
        "statsTable": statsTable,
    }


def _style_stats_table(table):
    """*apply the house style to the statistics table: thin grey edges, a bold header and a left-aligned label column*

    **Key Arguments:**

    - ``table`` -- the `matplotlib.table.Table` to style in place
    """
    table.auto_set_font_size(False)
    for (rowIndex, columnIndex), cell in table.get_celld().items():
        cell.set_fontsize(STATS_TABLE_FONT_SIZE)
        cell.set_edgecolor(STATS_TABLE_EDGE_COLOUR)
        cell.set_linewidth(STATS_TABLE_EDGE_WIDTH)
        if rowIndex == 0:
            cell.set_text_props(fontweight="bold")
        if columnIndex == 0:
            cell.set_text_props(ha="left")
            cell.PAD = STATS_TABLE_LABEL_PAD


def _draw_bias_distribution(ax, rawPixels, masterPixels, summary, arm):
    """*draw the filled raw and master bias histograms, their mean and median lines, and the statistics table onto `ax`*

    Each histogram is a translucent fill in the series colour with no outline. The statistics table sits below the
    x-axis label, in axes coordinates, so it stays inside a figure saved with `bbox_inches="tight"`.

    **Key Arguments:**

    - ``ax`` -- the matplotlib axes to draw on
    - ``rawPixels`` -- the plotted raw-frame pixel values (e-)
    - ``masterPixels`` -- the plotted master-bias pixel values (e-)
    - ``summary`` -- the dictionary returned by `bias_distribution_summary`; its `statsTable` fills the table
    - ``arm`` -- the arm name, used in the title

    **Usage:**

    ```python
    fig, ax = plt.subplots()
    _draw_bias_distribution(ax, rawPixels, masterPixels, summary, "VIS")
    ```
    """
    for label, pixels, mean, median, colour in (
        ("raw frame (earliest by MJD-OBS)", rawPixels, summary["rawMean"], summary["rawMedian"], RAW_SERIES_COLOUR),
        ("master bias", masterPixels, summary["masterMean"], summary["masterMedian"], MASTER_SERIES_COLOUR),
    ):
        ax.hist(
            pixels,
            bins=summary["binEdges"],
            histtype="stepfilled",
            color=colour,
            alpha=HISTOGRAM_FILL_ALPHA,
            linewidth=0,
            label=label,
        )
        ax.axvline(mean, color=colour, linestyle="-", linewidth=0.8)
        ax.axvline(median, color=colour, linestyle=":", linewidth=0.8)

    ax.set_yscale("log")
    ax.set_xlabel("flux (e-)")
    ax.set_ylabel("pixel count")
    ax.set_title(f"{arm} bias pixel-flux distribution (solid: mean, dotted: median)", fontsize=9)
    ax.legend(fontsize=8, loc="upper right")

    statsTable = summary["statsTable"]
    table = ax.table(
        cellText=statsTable["rows"],
        colLabels=statsTable["columns"],
        cellLoc="right",
        colLoc="right",
        bbox=STATS_TABLE_BBOX,
    )
    _style_stats_table(table)


class soxs_mbias(base_recipe):
    """
    *The* `soxs_mbias` *recipe is used to generate a master-bias frame from a set of input raw bias frames. The recipe is used only for the UV-VIS arm as NIR frames have bias (and dark current) removed by subtracting an off-frame of equal exposure length.*

    **Key Arguments**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary
    - ``inputFrames`` -- input fits frames. Can be a directory, a set-of-files (SOF) file or a list of fits frame paths.
    - ``verbose`` -- verbose. True or False. Default *False*
    - ``overwrite`` -- overwrite the product file if it already exists. Default *False*
    - ``command`` -- the command called to run the recipe
    - ``debug`` -- debug mode. True or False. Default *False*
    - ``turnOffMP`` -- turn off multiprocessing. True or False. Default *False*. If True, multiprocessing will be turned off and the recipe will run in serial. This is useful for debugging.

    **Usage**

    ```python
    from soxspipe.recipes import soxs_mbias
    mbiasFrame = soxs_mbias(
        log=log,
        settings=settings,
        inputFrames=fileList
    ).produce_product()
    ```
    """

    # INITIALISATION

    def __init__(
        self,
        log,
        settings=False,
        inputFrames=[],
        verbose=False,
        overwrite=False,
        command=False,
        debug=False,
        turnOffMP=False,
    ):
        # INHERIT INITIALISATION FROM  base_recipe
        this = super().__init__(
            log=log,
            settings=settings,
            inputFrames=inputFrames,
            overwrite=overwrite,
            recipeName="soxs-mbias",
            command=command,
            debug=debug,
            verbose=verbose,
            turnOffMP=turnOffMP,
        )
        log.debug("instantiating a new 'soxs_mbias' object")
        self.settings = settings
        self.inputFrames = inputFrames
        self.verbose = verbose

        # INITIAL ACTIONS
        # CONVERT INPUT FILES TO A CCDPROC IMAGE COLLECTION (inputFrames >
        # imagefilecollection)
        from soxspipe.commonutils.set_of_files import set_of_files

        sof = set_of_files(
            log=self.log,
            settings=self.settings,
            inputFrames=self.inputFrames,
            ext=self.settings["data-extension"],
        )
        self.inputFrames, self.supplementaryInput = sof.get()

        # VERIFY THE FRAMES ARE THE ONES EXPECTED BY SOXS_MBIAS - NO MORE, NO LESS.
        # PRINT SUMMARY OF FILES.
        self.log.print("# VERIFYING INPUT FRAMES")
        self.verify_input_frames()
        sys.stdout.flush()
        sys.stdout.write("\x1b[1A\x1b[2K")
        self.log.print("# VERIFYING INPUT FRAMES - ALL GOOD")

        # SORT IMAGE COLLECTION
        self.inputFrames.sort(["MJD-OBS"])

        # PREPARE THE FRAMES - CONVERT TO ELECTRONS, ADD UNCERTAINTY AND MASK
        # EXTENSIONS
        self.inputFrames = self.prepare_frames(save=self.settings["save-intermediate-products"])

        return

    def verify_input_frames(self):
        """*verify the input frame match those required by the soxs_mbias recipe*

        If the fits files conform to the required input for the recipe, everything will pass silently; otherwise, an exception will be raised.
        """
        self.log.debug("starting the ``verify_input_frames`` method")

        kw = self.kw

        # BASIC VERIFICATION COMMON TO ALL RECIPES
        imageTypes, imageTech, imageCat = self._verify_input_frames_basics()

        error = False

        # MIXED INPUT IMAGE TYPES ARE BAD
        if len(imageTypes) > 1:
            error = "Input frames are a mix of %(imageTypes)s" % locals()
        # NON-BIAS INPUT IMAGE TYPES ARE BAD
        elif imageTypes[0] != "BIAS":
            error = "Input frames not BIAS frames"

        if error:
            sys.stdout.flush()
            sys.stdout.write("\x1b[1A\x1b[2K")
            self.log.error("# VERIFYING INPUT FRAMES - **ERROR**\n")
            self.log.print(self.inputFrames.summary)
            self.log.print("")
            raise TypeError(error)

        self.imageType = imageTypes[0]

        self.log.debug("completed the ``verify_input_frames`` method")
        return

    def produce_product(self):
        """*generate a master bias frame*

        **Return:**

        - ``productPath`` -- the path to the master bias frame
        """
        self.log.debug("starting the ``produce_product`` method")

        from soxspipe.commonutils import toolkit

        arm = self.arm
        kw = self.kw
        dp = self.detectorParams

        (
            combined_bias_mean,
            masterMedianBiasLevel,
            rawRon,
            masterRon,
            rawFrameSample,
        ) = self._combine_bias_frames()

        # OPTIMISE: 24%
        self.qc_periodic_pattern_noise(frames=self.inputFrames)

        self.qc_ron(
            frameType="MBIAS",
            frameName="master bias",
            masterFrame=combined_bias_mean,
            rawRon=rawRon,
            masterRon=masterRon,
        )

        self.qc_bias_structure(combined_bias_mean)

        # ADD QUALITY CHECKS
        self.qc = generic_quality_checks(
            log=self.log,
            frame=combined_bias_mean,
            settings=self.settings,
            recipeName=self.recipeName,
            qcTable=self.qc,
        )

        medianFlux = self.qc_median_flux_level(
            frame=combined_bias_mean,
            frameType="MBIAS",
            frameName="master bias",
            medianFlux=masterMedianBiasLevel,
        )

        self.update_fits_keywords(frame=combined_bias_mean)

        # WRITE TO DISK
        toolkit.frame_to_32(combined_bias_mean)
        productPath = self._write(
            frame=combined_bias_mean,
            filedir=self.workspaceRootPath,
            filename=False,
            overwrite=True,
        )
        filename = os.path.basename(productPath)

        utcnow = utcnow_string()

        self.dateObs = combined_bias_mean.header[self.kw("DATE_OBS")]

        # PLOT BEFORE THE MBIAS ROW SO THE MASTER BIAS STAYS THE LAST PRODUCT
        self.plot_bias_distribution_qc(
            rawFrameSample=rawFrameSample,
            masterFrame=combined_bias_mean,
            rawRon=rawRon,
            masterRon=masterRon,
            productPath=productPath,
        )

        self.add_product(
            productLabel="MBIAS",
            fileName=filename,
            filePath=productPath,
            productDesc=f"{self.arm} Master bias frame",
            reductionDateUtc=utcnow,
            fileType="FITS",
            label="PROD",
        )

        qcTable = self.report_output()
        self.clean_up()

        self.log.debug("completed the ``produce_product`` method")
        return productPath, qcTable

    def _combine_bias_frames(self):
        """*stack the raw bias frames into the master-bias frame*

        **Return:**

        - ``combined_bias_mean`` -- the stacked master-bias frame
        - ``masterMedianBiasLevel`` -- the median of the per-frame mean bias levels
        - ``rawRon`` -- the mean of the per-frame read-out noise values
        - ``masterRon`` -- the standard deviation of the stacked noise frame
        - ``rawFrameSample`` -- the earliest raw frame (by MJD-OBS) for the QC plot, a dictionary with the keys
          `pixels` (flux of its unclipped, unflagged pixels) and `frameRon` (its RON)

        **Usage:**

        ```python
        combined_bias_mean, masterMedianBiasLevel, rawRon, masterRon, rawFrameSample = self._combine_bias_frames()
        ```
        """
        import numpy as np

        # LIST OF CCDDATA OBJECTS
        # OPTIMISE: 9%
        ccds = list(
            self.inputFrames.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        )

        # COPY THE BAD-PIXEL MASK OF THE FIRST FRAME: subtract_mean_flux_level REPLACES IT WITH THE CLIP MASK
        badPixelMask = None if ccds[0].mask is None else np.array(ccds[0].mask, dtype=bool)

        # OPTIMISE: 33%
        # `strict=False` IS THE CURRENT BEHAVIOUR MADE EXPLICIT, NOT A CHANGE
        meanBiasLevels, rons, noiseFrames = zip(*[self.subtract_mean_flux_level(c) for c in ccds], strict=False)
        masterMeanBiasLevel = np.mean(meanBiasLevels)
        masterMedianBiasLevel = np.median(meanBiasLevels)
        rawRon = np.mean(rons)

        # SAMPLE BEFORE STACKING IN CASE THE STACK TOUCHES THE INPUT FRAMES
        rawFrameSample = {
            "pixels": _unclipped_raw_pixels(noiseFrames[0], badPixelMask, meanBiasLevels[0]),
            "frameRon": float(rons[0]),
        }

        # OPTIMISE: 19%
        combined_noise = self.clip_and_stack(
            frames=list(noiseFrames),
            recipe="soxs_mbias",
            ignore_input_masks=False,
            post_stack_clipping=True,
        )

        # ESTIMATE THE RON FROM THE PIXELS THAT SURVIVED CLIPPING, AS soxs_mdark DOES
        masterRon = np.std(np.ma.array(combined_noise.data, mask=combined_noise.mask))

        # USE COMBINED NOISE MASK AS MBIAS MASK
        combined_noise.data = (
            np.ma.array(combined_noise.data, mask=combined_noise.mask, fill_value=0).filled() + masterMeanBiasLevel
        )
        combined_noise.uncertainty = np.ma.array(
            combined_noise.uncertainty.array,
            mask=combined_noise.mask,
            fill_value=rawRon,
        ).filled()
        combined_bias_mean = combined_noise
        combined_bias_mean.mask = combined_noise.mask

        return combined_bias_mean, masterMedianBiasLevel, rawRon, masterRon, rawFrameSample

    def plot_bias_distribution_qc(self, rawFrameSample, masterFrame, rawRon, masterRon, productPath):
        """*plot the raw and master bias pixel-flux distributions and register the plot as a QC product*

        **Key Arguments:**

        - ``rawFrameSample`` -- the earliest raw frame's `pixels` and `frameRon`
          (see `_combine_bias_frames`)
        - ``masterFrame`` -- the master bias frame
        - ``rawRon`` -- the mean of the per-frame RON values (the RAW RON QC)
        - ``masterRon`` -- the MASTER RON QC value
        - ``productPath`` -- the path of the master bias product, used to name the plot

        **Return:**

        - ``filePath`` -- the path to the PDF plot

        **Usage:**

        ```python
        filePath = self.plot_bias_distribution_qc(rawFrameSample, masterFrame, rawRon, masterRon, productPath)
        ```
        """
        self.log.debug("starting the ``plot_bias_distribution_qc`` method")

        import matplotlib.pyplot as plt

        masterPixels = _finite_unmasked_pixels(masterFrame.data, masterFrame.mask)

        summary = bias_distribution_summary(
            rawPixels=rawFrameSample["pixels"],
            frameRon=rawFrameSample["frameRon"],
            masterPixels=masterPixels,
            rawRon=rawRon,
            masterRon=masterRon,
        )

        filename = os.path.basename(productPath).replace(".fits", "_BIAS_DISTRIBUTION_QC_PLOT.pdf")
        filePath = os.path.join(self.qcDir, filename)

        fig, ax = plt.subplots()
        try:
            _draw_bias_distribution(ax, rawFrameSample["pixels"], masterPixels, summary, self.arm)
            fig.savefig(filePath, format="pdf", bbox_inches="tight")
        finally:
            plt.close(fig)

        self.add_product(
            productLabel="BIAS_DISTRIBUTION_QC_PLOT",
            fileName=filename,
            filePath=filePath,
            productDesc=f"{self.arm} raw vs master bias pixel-flux distribution",
            reductionDateUtc=utcnow_string(),
            fileType="PDF",
            label="QC",
        )

        self.log.debug("completed the ``plot_bias_distribution_qc`` method")
        return filePath

    def qc_bias_structure(self, combined_bias_mean):
        """*calculate the structure of the bias*

        **Key Arguments:**

        - ``combined_bias_mean`` -- the mbias frame

        **Return:**

        - ``structx`` -- slope of BIAS in X direction
        - ``structx`` -- slope of BIAS in Y direction

        **Usage:**

        ```python
        structx, structy = self.qc_bias_structure(combined_bias_mean)
        ```
        """
        self.log.debug("starting the ``qc_bias_structure`` method")

        import numpy as np

        plot = False

        collaps_ax1 = np.nansum(combined_bias_mean, axis=0)
        collaps_ax2 = np.nansum(combined_bias_mean, axis=1)

        x_axis = np.linspace(0, len(collaps_ax1), len(collaps_ax1), dtype=int)
        y_axis = np.linspace(0, len(collaps_ax2), len(collaps_ax2), dtype=int)

        # FIT WITH A LINE AND COLLECT THE SLOPE
        coeff_ax1 = np.polyfit(x_axis, collaps_ax1, deg=1)
        coeff_ax2 = np.polyfit(y_axis, collaps_ax2, deg=1)

        if plot == True:
            import matplotlib.pyplot as plt

            plt.plot(x_axis, collaps_ax1)
            plt.plot(x_axis, np.polyval(coeff_ax1, x_axis))
            plt.xlabel("x-axis")
            plt.ylabel("Summed Pixel Values")
            plt.show()

            plt.plot(y_axis, collaps_ax2)
            plt.plot(y_axis, np.polyval(coeff_ax2, y_axis))
            plt.xlabel("y-axis")
            plt.ylabel("Summed Pixel Values")
            plt.show()

        # ONE TIMESTAMP DELIBERATELY SHARED ACROSS BOTH QC ROWS BELOW
        utcnow = utcnow_string()

        self.add_qc(
            qcName="STRUCTX",
            qcValue=coeff_ax1[0],
            qcComment="Slope of BIAS in X direction",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.add_qc(
            qcName="STRUCTY",
            qcValue=coeff_ax2[0],
            qcComment="Slope of BIAS in Y direction",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.log.debug("completed the ``qc_bias_structure`` method")
        return coeff_ax1[0], coeff_ax2[0]

    def _periodic_noise_ratio(self, frame):
        """*return the periodic-pattern-noise ratio of a single raw bias frame*

        **Key Arguments:**

        - ``frame`` -- a single raw bias frame (CCDData)

        **Return:**

        - ``ratio`` -- the standard deviation of the sigma-clipped 2D FFT divided by its median absolute deviation

        **Usage:**

        ```python
        ratio = self._periodic_noise_ratio(frame)
        ```
        """
        import numpy as np
        from astropy.stats import sigma_clip
        from ccdproc import block_reduce
        from scipy.stats import median_abs_deviation

        from soxspipe.commonutils.toolkit import quicklook_image

        # FORCE CONVERSION OF CCDData OBJECT TO NUMPY ARRAY
        maskedDataArray = np.ma.array(frame.data, mask=frame.mask)
        # BIN THE FRAME TO INCREASE SPEED
        maskedDataArray = block_reduce(maskedDataArray, 5, np.mean)
        dark_image_grey_fourier = np.fft.fftshift(np.fft.fft2(maskedDataArray.filled(np.median(frame.data))))

        # SIGMA-CLIP THE DATA
        masked_dark_image_grey_fourier = sigma_clip(
            dark_image_grey_fourier,
            sigma_lower=100,
            sigma_upper=100,
            maxiters=1,
            cenfunc="mean",
        )
        goodData = np.ma.compressed(masked_dark_image_grey_fourier)

        frame_mad = median_abs_deviation(goodData, axis=None)
        frame_std = np.std(goodData)

        quicklook_image(
            log=self.log,
            CCDObject=abs(masked_dark_image_grey_fourier),
            show=False,
            ext=None,
            stdWindow=0.1,
        )
        quicklook_image(
            log=self.log,
            CCDObject=abs(dark_image_grey_fourier),
            show=False,
            ext=None,
            stdWindow=0.1,
        )

        return frame_std / frame_mad

    def qc_periodic_pattern_noise(self, frames):
        """*calculate the periodic pattern noise based on the raw input bias frames*

        A 2D FFT is applied to each of the raw bias frames and the standard deviation and median absolute deviation calcualted for each result. The maximum std/mad is then added as the ppnmax QC in the master bias frame header.

        **Key Arguments:**

        - ``frames`` -- the raw bias frames (imageFileCollection)

        **Return:**

        - ``ppnmax``

        **Usage:**

        ```python
        self.qc_periodic_pattern_noise(frames=self.inputFrames)
        ```
        """
        self.log.debug("starting the ``qc_periodic_pattern_noise`` method")

        # LIST OF CCDDATA OBJECTS
        ccds = [
            c
            for c in frames.ccds(
                ccd_kwargs={
                    "hdu_uncertainty": "ERRS",
                    "hdu_mask": "QUAL",
                    "hdu_flags": "FLAGS",
                    "key_uncertainty_type": "UTYPE",
                }
            )
        ]

        ratios = [self._periodic_noise_ratio(frame) for frame in ccds]

        utcnow = utcnow_string()

        ppnmax = max(ratios)

        self.add_qc(
            qcName="FPN FRACMAX",
            qcValue=ppnmax,
            qcComment="Max periodic pattern noise ratio in raw bias frames",
            reductionDateUtc=utcnow,
            qcUnit=None,
            toHeader=True,
        )

        self.log.debug("completed the ``qc_periodic_pattern_noise`` method")
        return ppnmax

    # USE THE TAB-TRIGGER BELOW FOR NEW METHOD
    # xt-class-method
