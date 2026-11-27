"""
*Combine the order-merged extractions of the nod positions or cycles of a recipe into one spectrum*

Author
: David Young

Date Created
: November 23, 2026
"""

from __future__ import annotations

import math
from typing import Any

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.on_target_frames import with_frame_accounting, without_pointing
from soxspipe.commonutils.spectrum_names import SPECTRUM_ASCII_LABEL, SPECTRUM_SUFFIX, SPECTRUM_TABLE_LABEL
from soxspipe.commonutils.toolkit import utcnow_string

# UP TO THIS MANY SPECTRA, THE MEDIAN IS THE MEAN, SO ITS VARIANCE IS SUM(VAR) / N**2
MEDIAN_IS_MEAN_MAX_SPECTRA = 2
# THE LARGE-N VARIANCE OF THE MEDIAN OF GAUSSIAN SAMPLES IS THIS FACTOR TIMES THE VARIANCE OF THEIR MEAN
MEDIAN_VARIANCE_FACTOR = math.pi / 2


def median_stack_variance(spectra):
    """*the variance of the per-wavelength median of a set of stacked spectra*

    A row contributes at its wavelength when it has a flux, as for the flux median, so N is counted per
    wavelength. For N <= 2 the median is the mean and the variance is sum(var) / N**2. For larger N it is
    (pi / 2) * mean(var) / N, the large-N variance of the median of Gaussian samples. The variance is NaN where
    no row contributes, or where a contributing row has no variance.

    **Key Arguments:**

    - ``spectra`` -- the concatenated spectra (dataframe) with rounded ``WAVE``, ``FLUX_COUNTS`` and ``VARIANCE``

    **Return:**

    - ``variance`` -- the variance of the median at each ``WAVE`` (series indexed by ``WAVE``, sorted ascending)

    **Usage:**

    ```python
    variance = median_stack_variance(mergedDataframe)
    ```
    """
    import numpy as np
    import pandas as pd

    # WORK ON PLAIN FLOAT ARRAYS IN FRESH SERIES. COPYING THE FRAME CONSOLIDATES ITS ASTROPY QUANTITY COLUMNS,
    # WHICH FAILS ON THEIR MIXED UNITS, AND CONCATENATED SPECTRA REPEAT INDEX LABELS, WHICH BREAKS ALIGNMENT
    flux = np.asarray(spectra["FLUX_COUNTS"].to_numpy(), dtype=np.float64)
    contributing = pd.Series(~np.isnan(flux))
    variance = pd.Series(np.asarray(spectra["VARIANCE"].to_numpy(), dtype=np.float64)).where(contributing)
    byWave = pd.Series(np.asarray(spectra["WAVE"].to_numpy(), dtype=np.float64))
    contributorCount = contributing.groupby(byWave).sum()
    hasMissingVariance = (contributing & variance.isna()).groupby(byWave).any()
    varianceSum = variance.groupby(byWave).sum(min_count=1).mask(hasMissingVariance)
    meanVariance = varianceSum / contributorCount**2
    medianVariance = MEDIAN_VARIANCE_FACTOR * (varianceSum / contributorCount) / contributorCount
    return meanVariance.where(contributorCount <= MEDIAN_IS_MEAN_MAX_SPECTRA, medianVariance)


class extraction_stacking:
    """
    Mixin that gives a recipe `stack_extractions`, which combines per-position or per-cycle extractions.

    The nod recipe combines its cycles with it. The stare recipe, reducing a nodding SOF in stare mode, combines the
    extractions of the A and B stacks with it. The recipe supplies ``log``, ``recipeName``, ``productStem``,
    ``settings``, ``arm``, ``qc``, ``productDir``, ``masterHeaderFrame``, ``add_product`` and ``update_fits_keywords``.
    It may also set ``onTargetFrames``, the on-target frames of the product, which give the merged spectrum its
    frame-accounting cards (DY-2356).
    """

    # THE ON-TARGET FRAMES THAT REACH THE PRODUCT BEING WRITTEN, OR None WHEN THE RECIPE HAS NOT BUILT THE LIST. AN
    # IMMUTABLE DEFAULT, REBOUND (NEVER MUTATED) AS FRAMES LEAVE THE LIST
    onTargetFrames = None

    def _merged_spectrum_header(self) -> Any:
        """*the header of the merged spectrum: the master header with the frame-accounting cards of the on-target list*

        ``stack_extractions`` stamps ``masterHeaderFrame`` with every raw frame of the recipe, so the calibrated
        spectrum lists them all. The header of one input frame lists the raw frames of its own cycle only, or none. The
        cards go on a copy, because flux calibration and the response divide the counts by the exposure time of one
        frame, which is the ``EXPTIME`` of the master header. With no on-target list a copy of the master header comes
        back with only its pointing ``RA`` and ``DEC`` removed.

        **Return:**

        - ``header`` -- a copy of the header of ``masterHeaderFrame``, with the frame-accounting and target position
          cards when there is an on-target list
        """
        return with_frame_accounting(self.masterHeaderFrame.header, self.onTargetFrames, self.log)

    def stack_extractions(self, dataFrameList, notFlattened=False, orderJoins=None):
        """*merge the extractions of the nodding cycles, or of the nod positions of a stare-mode nodding SOF*

        The spectra are grouped by wavelength and combined by their median, which is the mean for two spectra.

        **Key Arguments:**

        - ``dataFrameList`` -- a list of order-merged spectrum dataframes
        - ``notFlattened`` -- if True, the extraction was performed on non-flattened data and the products
          are written with a ``_NOTFLAT`` suffix. Default *False*
        - ``orderJoins`` -- the order joins the extraction reported. Default *None*

        **Return:**

        - ``stackedSpectrum`` -- the combined spectrum, as an astropy table
        - ``filePath`` -- the path the combined spectrum was written to

        **Usage:**

        ```python
        stackedSpectrum, filePath = soxs_nod.stack_extractions(
            [mergedSpectrumDF_A, mergedSpectrumDF_B], orderJoins=orderJoins)
        ```
        """
        self.log.debug("starting the ``stack_extractions`` method")

        import numpy as np
        import pandas as pd
        from astropy import units as u
        from astropy.table import Table

        from soxspipe.commonutils.fits_headers import set_header_card
        from soxspipe.commonutils.phase3 import PER_ARM_MODE, with_flux_calibration_cards, write_fits_table_to_disk
        from soxspipe.commonutils.specutils_calls import spectrum1d_quietly
        from soxspipe.commonutils.toolkit import (
            add_snr_efficiency_qcs,
        )

        postfix = "_NOTFLAT" if notFlattened else ""

        # AN EXTRACTION THAT FOUND NO OBJECT TRACE RETURNS NO SPECTRUM, AND WITH NONE TO STACK THERE IS NO PRODUCT
        if all(spectrum is None for spectrum in dataFrameList):
            where = (
                "in the offset frame of any cycle"
                if "offset" in self.recipeName
                else "in either nod position of any nodding cycle"
            )
            error = (
                f"No object trace was found {where}, so there is no extracted spectrum to stack. "
                "Please check that the object is in the slit and bright enough to be traced."
            )
            self.log.error(error)
            raise Exception(error)

        # MERGE THE PANDAS DATAFRAMES MERDGED_ORDERS_A AND MERGEDSPECTRUMDF_B INTO A SINGLE DATAFRAME, THEN
        # GROUP BY WAVE AND TAKE THE MEDIAN OF EACH COLUMN

        merged_dataframe = pd.concat(dataFrameList)
        # BEFORE GROUPING, WE NEED TO TRUNCATE THE WAVELENGTH TO THE 4 DIGITS
        merged_dataframe["WAVE"] = merged_dataframe["WAVE"].apply(lambda x: round(float(x.value), 4))
        groupedDataframe = merged_dataframe.groupby(by="WAVE", as_index=False).median()
        # THE FLUX IS THE MEDIAN, BUT THE MEDIAN OF THE INPUT VARIANCES IS THE ERROR OF ONE SPECTRUM, NOT OF THE STACK
        # KEEP THE EXTRACTION'S FLOAT32 VARIANCE COLUMN TYPE; INTEGER OR FLOAT64 INPUT GIVES FLOAT64
        varianceType = np.promote_types(np.asarray(merged_dataframe["VARIANCE"].to_numpy()).dtype, np.float32)
        groupedDataframe["VARIANCE"] = (
            groupedDataframe["WAVE"].map(median_stack_variance(merged_dataframe)).astype(varianceType)
        )

        self.filenameTemplate = self.productStem + ".fits"

        # PREPARING THE HEADER
        kw = keyword_lookup(log=self.log, settings=self.settings).get

        # SELECTING HEADER A_MINUS_B (IS THIS THE SAME?)
        self.update_fits_keywords(frame=self.masterHeaderFrame)
        # THE NOT-FLATTENED SPECTRUM FEEDS THE RESPONSE ONLY, SO IT CARRIES NO FRAME ACCOUNTING
        # ITS HEADER IS THE MASTER HEADER WITHOUT THE POINTING RA AND DEC (DY-2357)
        header = without_pointing(self.masterHeaderFrame.header) if notFlattened else self._merged_spectrum_header()

        set_header_card(header, kw("PRO_TYPE"), "REDUCED")
        set_header_card(header, kw("PRO_CATG"), f"SCI_SLIT_FLUX_{self.arm}".upper())
        # THE STACK IS IN COUNTS. FLUX CALIBRATION REWRITES THE FILE WITH THE ABSOLUTE CARDS (DY-2357)
        header = with_flux_calibration_cards(header, isCalibrated=False)

        flux_orig = groupedDataframe["FLUX_COUNTS"].values * u.electron
        # THE RESULT IS UNUSED, BUT BUILDING IT VALIDATES THE GROUPED FLUX AND WAVELENGTH ARRAYS AGAINST EACH
        # OTHER AND RAISES WHEN THEY DISAGREE. DELETING IT WOULD REMOVE THAT CHECK.
        spectrum_orig = spectrum1d_quietly(  # noqa: F841
            flux=flux_orig,
            spectral_axis=groupedDataframe["WAVE"].values * u.nm,
            bin_specification="center",
        )

        groupedDataframe["SNR"] = groupedDataframe["FLUX_COUNTS"].values / np.sqrt(groupedDataframe["VARIANCE"].values)

        # PREPARING THE HDU. WAVE IS ALREADY ROUNDED TO 4 DECIMALS BY THE GROUPING ABOVE. ROUNDING IT TO 2 WOULD MOVE A
        # LOG-GRID POINT BY UP TO HALF OF A BIN, SO ONLY THE LINEAR X-SHOOTER GRIDS KEEP THE 2-DECIMAL WAVE (DY-2348)
        roundingDecimals = [("FLUX_COUNTS", 3), ("SNR", 2), ("FLUX_DENSITY_COUNTS", 3)]
        if str(self.settings.get("instrument", "soxs")).lower() != "soxs":
            roundingDecimals.insert(0, ("WAVE", 2))
        for col, decimals in roundingDecimals:
            # ONLY A NOT-FLATTENED STACK HAS A FLUX DENSITY COLUMN (DY-1622)
            if col not in groupedDataframe.columns:
                continue
            # `decimals` IS BOUND IN A DEFAULT ARGUMENT SO THE LAMBDA CANNOT SEE A LATER LOOP VALUE. `.apply`
            # CALLS IT WITHIN THE SAME ITERATION, SO THIS CHANGES NOTHING; IT ONLY MAKES THAT EXPLICIT.
            groupedDataframe[col] = groupedDataframe[col].apply(lambda x, decimals=decimals: round(float(x), decimals))
        stackedSpectrum = Table.from_pandas(groupedDataframe, index=False)

        self.utcnow = utcnow_string()
        self.dateObs = header[kw("DATE_OBS")]

        self.qc = add_snr_efficiency_qcs(
            log=self.log,
            spectrumDF=groupedDataframe,
            qcTable=self.qc,
            orderJoins=orderJoins,
            recipeName=self.recipeName,
            dateObs=self.dateObs,
            arm=self.arm,
        )

        # WRITE PRODUCT TO DISK
        filename = self.filenameTemplate.replace(".fits", SPECTRUM_SUFFIX + postfix + ".fits")
        filePath = f"{self.productDir}/{filename}"

        write_fits_table_to_disk(
            log=self.log,
            settings=self.settings,
            header=header,
            tables=[stackedSpectrum],
            filePath=filePath,
            qc=self.qc,
            phase3=PER_ARM_MODE,
            extname="SPECTRUM",
        )

        self.add_product(
            productLabel=SPECTRUM_TABLE_LABEL,
            fileName=filename,
            filePath=filePath,
            productDesc="Table of the extracted, order-merged spectrum. All nodding cycles combined.",
            reductionDateUtc=self.utcnow,
            fileType="FITS",
            label="PROD",
        )

        # THE FLATTENED ASCII SPECTRUM IS WRITTEN AFTER FLUX CALIBRATION (SEE _write_ascii_spectrum). ONLY THE
        # NOT-FLATTENED STACK, WHICH IS NEVER FLUX CALIBRATED, GETS ITS ASCII FILE HERE
        if notFlattened:
            asciiFilename = filename.replace(".fits", ".txt")
            asciiFilePath = f"{self.productDir}/{asciiFilename}"
            stackedSpectrum2 = stackedSpectrum.copy()
            stackedSpectrum2["WAVE"] = stackedSpectrum2["WAVE"] * 10
            stackedSpectrum2["WAVE"].format = "{:.2f}"  # CONVERTING TO ANGSTROMS
            stackedSpectrum2.write(asciiFilePath, format="ascii", overwrite=True)
            self.add_product(
                productLabel=SPECTRUM_ASCII_LABEL,
                fileName=asciiFilename,
                filePath=asciiFilePath,
                productDesc="Ascii version of extracted source spectrum",
                reductionDateUtc=self.utcnow,
                fileType="TXT",
                label="PROD",
            )

        self.log.debug("completed the ``stack_extractions`` method")
        return stackedSpectrum, filePath
