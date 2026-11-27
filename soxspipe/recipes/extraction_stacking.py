"""
*Combine the order-merged extractions of the nod positions or cycles of a recipe into one spectrum*

Author
: David Young

Date Created
: November 23, 2026
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING, Any

from soxspipe.commonutils import keyword_lookup
from soxspipe.commonutils.log_wavelength_grid import grid_index, grid_step_kms, grid_wavelength
from soxspipe.commonutils.on_target_frames import without_pointing
from soxspipe.commonutils.response_basis import dispersion_binning
from soxspipe.commonutils.spectrum_names import SPECTRUM_ASCII_LABEL, SPECTRUM_SUFFIX, SPECTRUM_TABLE_LABEL
from soxspipe.commonutils.toolkit import utcnow_string
from soxspipe.commonutils.wavelength_keywords import with_merged_spectrum_cards

if TYPE_CHECKING:
    import numpy as np
    import pandas as pd

# UP TO THIS MANY SPECTRA, THE MEDIAN IS THE MEAN, SO ITS VARIANCE IS SUM(VAR) / N**2
MEDIAN_IS_MEAN_MAX_SPECTRA = 2
# THE LARGE-N VARIANCE OF THE MEDIAN OF GAUSSIAN SAMPLES IS THIS FACTOR TIMES THE VARIANCE OF THEIR MEAN
MEDIAN_VARIANCE_FACTOR = math.pi / 2
# X-SHOOTER, WHICH KEEPS LINEAR GRIDS, GROUPS THE BINS OF A STACK BY THEIR WAVELENGTH TO THIS MANY DECIMALS (NM)
LINEAR_GRID_GROUPING_DECIMALS = 4
# THE INPUT WAVELENGTHS OF A SOXS STACK ARE GRID POINTS ALREADY. ONE THAT IS FARTHER FROM ITS GRID POINT THAN THIS
# (RELATIVE) WAS MADE ON ANOTHER GRID, SUCH AS ONE WITH A DIFFERENT BINNING
GRID_POINT_TOLERANCE = 1.0e-6


def median_stack_variance(spectra: pd.DataFrame, groupKeys: np.ndarray) -> pd.Series:
    """*the variance of the per-bin median of a set of stacked spectra*

    A row contributes to its bin when it has a flux, as for the flux median, so N is counted per bin. For N <= 2 the
    median is the mean and the variance is sum(var) / N**2. For larger N it is (pi / 2) * mean(var) / N, the large-N
    variance of the median of Gaussian samples. The variance is NaN where no row contributes, or where a contributing
    row has no variance.

    **Key Arguments:**

    - ``spectra`` -- the concatenated spectra (dataframe) with ``FLUX_COUNTS`` and ``VARIANCE``
    - ``groupKeys`` -- the key of the bin of each row of ``spectra``, in the order of its rows: the integer bin index on
      the log-wavelength grid for SOXS, the rounded wavelength for X-Shooter

    **Return:**

    - ``variance`` -- the variance of the median in each bin (series indexed by the group key, sorted ascending)

    **Usage:**

    ```python
    variance = median_stack_variance(mergedDataframe, groupKeys)
    ```
    """
    import numpy as np
    import pandas as pd

    # WORK ON PLAIN FLOAT ARRAYS IN FRESH SERIES. COPYING THE FRAME CONSOLIDATES ITS ASTROPY QUANTITY COLUMNS,
    # WHICH FAILS ON THEIR MIXED UNITS, AND CONCATENATED SPECTRA REPEAT INDEX LABELS, WHICH BREAKS ALIGNMENT
    flux = np.asarray(spectra["FLUX_COUNTS"].to_numpy(), dtype=np.float64)
    contributing = pd.Series(~np.isnan(flux))
    variance = pd.Series(np.asarray(spectra["VARIANCE"].to_numpy(), dtype=np.float64)).where(contributing)
    byBin = pd.Series(np.asarray(groupKeys))
    contributorCount = contributing.groupby(byBin).sum()
    hasMissingVariance = (contributing & variance.isna()).groupby(byBin).any()
    varianceSum = variance.groupby(byBin).sum(min_count=1).mask(hasMissingVariance)
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

    # THE PATH TO THE DISP_TAB THE RECIPE REDUCED ITS FRAMES WITH, OR None WHEN IT HAS NOT FOUND ONE. IT GIVES THE
    # MERGED SPECTRUM ITS LAMRMS, LAMNLIN AND SPEC_ERR CARDS (DY-2358)
    dispMap = None

    def _merged_spectrum_header(self, header: Any = None) -> Any:
        """*the header of the merged spectrum: a frame header with the frame-accounting, target and wavelength cards*

        ``stack_extractions`` stamps ``masterHeaderFrame`` with every raw frame of the recipe, so the calibrated
        spectrum lists them all. The header of one input frame lists the raw frames of its own cycle only, or none. The
        cards go on a copy, because flux calibration and the response divide the counts by the exposure time of one
        frame, which is the ``EXPTIME`` of the master header. With no on-target list the copy has only its pointing
        ``RA`` and ``DEC`` removed. The wavelength-calibration and resolution cards are added either way (DY-2358).

        **Key Arguments:**

        - ``header`` -- the header to stamp. Default *None*, the header of ``masterHeaderFrame``

        **Return:**

        - ``header`` -- a stamped copy of the header, with the frame-accounting and target position cards when there is
          an on-target list, and ``SPEC_RES``, ``LAMRMS``, ``LAMNLIN`` and ``SPEC_ERR`` where they can be computed
        """
        kw = keyword_lookup(log=self.log, settings=self.settings).get
        return with_merged_spectrum_cards(
            self.masterHeaderFrame.header if header is None else header,
            self.onTargetFrames,
            self.dispMap,
            self.arm,
            kw,
            self.settings,
            self.log,
        )

    def _is_soxs(self) -> bool:
        """*whether the recipe reduces SOXS data, whatever the case of the ``instrument`` setting*

        **Return:**

        - ``isSoxs`` -- *True* for SOXS, *False* for X-Shooter, which keeps its linear grids
        """
        return str(self.settings.get("instrument", "soxs")).lower() == "soxs"

    def _stack_grid_step_kms(self) -> float | None:
        """*the log-wavelength step of the grid the spectra to stack sit on (km/s), or *None* for X-Shooter*

        It is the step of the arm per unbinned pixel times the binning of the frames, as the merge used (DY-2348).

        **Return:**

        - ``stepKms`` -- the step per bin, or *None* when the instrument has no log-wavelength grid
        """
        if not self._is_soxs():
            return None
        kw = keyword_lookup(log=self.log, settings=self.settings).get
        binning = dispersion_binning(self.masterHeaderFrame.header, self.arm, kw)
        return grid_step_kms(self.detectorParams, binning)

    def _median_by_bin(self, mergedDataframe: pd.DataFrame) -> pd.DataFrame:
        """*the median of each column of the concatenated spectra in each bin, with the variance of that median*

        A SOXS bin is the integer index on the log-wavelength grid of the arm, and its ``WAVE`` is the grid wavelength
        of the index. An X-Shooter bin is the wavelength rounded to 4 decimals. The frame is not changed.

        **Key Arguments:**

        - ``mergedDataframe`` -- the concatenated spectra (dataframe) with ``WAVE`` as plain floats in nm, and a
          ``VARIANCE``

        **Return:**

        - ``medians`` -- one row per bin, sorted by wavelength, with ``WAVE`` and the median of every other column,
          except ``VARIANCE``, which is the variance of the median (see ``median_stack_variance``)

        **Raises:**

        - ``ValueError`` -- if a SOXS wavelength is not on the log-wavelength grid of the stack, which happens when the
          spectra were merged with another binning than the one the stack reads
        """
        import numpy as np

        waves = mergedDataframe["WAVE"].to_numpy(dtype=float)
        stepKms = self._stack_grid_step_kms()
        if stepKms is None:
            groupKeys = np.round(waves, LINEAR_GRID_GROUPING_DECIMALS)
        else:
            groupKeys = grid_index(waves, stepKms)
            if not np.allclose(grid_wavelength(groupKeys, stepKms), waves, rtol=GRID_POINT_TOLERANCE, atol=0.0):
                raise ValueError(
                    f"The wavelengths to stack are not on the {stepKms:.3f} km/s log-wavelength grid the binning of "
                    "the master header gives, so the spectra were merged with another binning. Please check the "
                    "binning cards of the frames."
                )
        medians = mergedDataframe.groupby(by=groupKeys).median()
        binKeys = medians.index.to_numpy()
        medians = medians.reset_index(drop=True)
        medians["WAVE"] = binKeys if stepKms is None else grid_wavelength(binKeys, stepKms)
        # THE FLUX IS THE MEDIAN, BUT THE MEDIAN OF THE INPUT VARIANCES IS THE ERROR OF ONE SPECTRUM, NOT OF THE STACK
        # KEEP THE EXTRACTION'S FLOAT32 VARIANCE COLUMN TYPE; INTEGER OR FLOAT64 INPUT GIVES FLOAT64
        varianceType = np.promote_types(np.asarray(mergedDataframe["VARIANCE"].to_numpy()).dtype, np.float32)
        medians["VARIANCE"] = (
            median_stack_variance(mergedDataframe, groupKeys).reindex(binKeys).to_numpy().astype(varianceType)
        )
        return medians

    def stack_extractions(self, dataFrameList, notFlattened=False, orderJoins=None):
        """*merge the extractions of the nodding cycles, or of the nod positions of a stare-mode nodding SOF*

        The spectra are grouped by bin and combined by their median, which is the mean for two spectra. A SOXS bin is
        its integer index on the log-wavelength grid of the arm, so no bin is split or merged by rounding, and the
        ``WAVE`` of the stack is the grid wavelength of the index (DY-2349). X-Shooter groups by the wavelength to 4
        decimals.

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
        # GROUP BY BIN AND TAKE THE MEDIAN OF EACH COLUMN

        merged_dataframe = pd.concat(dataFrameList)
        # PLAIN FLOATS IN NM. COPYING THE FRAME WOULD CONSOLIDATE ITS ASTROPY QUANTITY COLUMNS, WHICH FAILS ON THEIR
        # MIXED UNITS, SO THE CONCATENATED FRAME ITSELF TAKES THE FLOAT WAVELENGTHS
        merged_dataframe["WAVE"] = [float(x.value) for x in merged_dataframe["WAVE"]]
        groupedDataframe = self._median_by_bin(merged_dataframe)

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

        # PREPARING THE HDU. WAVE IS A GRID POINT (SOXS) OR ROUNDED TO 4 DECIMALS (X-SHOOTER) BY THE GROUPING ABOVE.
        # ROUNDING IT TO 2 WOULD MOVE A LOG-GRID POINT BY UP TO HALF OF A BIN, SO ONLY THE LINEAR X-SHOOTER GRIDS KEEP
        # THE 2-DECIMAL WAVE (DY-2348)
        roundingDecimals = [("FLUX_COUNTS", 3), ("SNR", 2), ("FLUX_DENSITY_COUNTS", 3)]
        if not self._is_soxs():
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
