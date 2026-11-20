# DY-2327: multi-arm keyword rules in the SDP standard

This note records what the ESO Science Data Products standard (version 8, released 2022-03-15) and ESO practice say about the header keywords of one 1D spectrum that joins the arms of a multi-arm spectrograph. It is for the people who will decide the design of the joined SOXS UV-VIS + NIR spectrum (Linear DY-2328). It reports rules and options only. It makes no decisions.

## Contents

- [Sources](#sources)
- [Summary](#summary)
- [What the standard says about multi-arm products](#what-the-standard-says-about-multi-arm-products)
- [Precedents](#precedents)
- [Keyword by keyword](#keyword-by-keyword)
- [Overlap and bin-width changes](#overlap-and-bin-width-changes)
- [Associations for per-arm ancillary files](#associations-for-per-arm-ancillary-files)
- [What the standard leaves open](#what-the-standard-leaves-open)

## Sources

- [SDP] ESO Science Data Products Standard, version 8, as Markdown in the P32.11 workspace at `output/p3sdpstd_md/p3sdpstd.md`. Section numbers below refer to it.
- [CHK] ESO Phase 3 checks summary, dated 2025-03-27: <https://www.eso.org/sci/observing/phase3/phase3_checks_summary.txt> (read 2026-11-20).
- [FAQ] ESO Phase 3 FAQ: <https://www.eso.org/sci/observing/phase3/faq.html> (read 2026-11-20). It has no entry on multi-arm or joined spectra.
- [XSH-RD] ESO Phase 3 Data Release Description, collection XSHOOTER, release 1, document version 2.9, 2026-07-21: <https://www.eso.org/rm/api/v1/public/releaseDescriptions/70>.
- [XQ-RD] ESO Phase 3 Data Release Description, collection XQ-100, release 1, 2016-06-22: <https://www.eso.org/sci/observing/phase3/data_releases/xq100_dr1.pdf>.
- [ARCH] ESO archive TAP service (`ivoa.ObsCore`) and header service (`https://archive.eso.org/hdr?DpId=<ADP id>`), queried 2026-11-20. Headers read: XQ-100 joint spectrum `ADP.2016-07-15T08:22:40.175` (`J1042+1957_rescale.fits`), its UVB and NIR arm spectra `ADP.2016-07-15T08:22:40.177` and `ADP.2016-07-15T08:22:40.423`, and an ESO X-Shooter NIR IDP `ADP.2025-12-19T10:57:26.336`.

## Summary

- The standard has exactly one rule written for multi-arm spectra. Section 5.1, `EXPTIME`: "When collating spectra from a multi-arm spectrograph, e.g. X-SHOOTER, EXPTIME shall be set to the median exposure time of the individual arms." With two arms, the median is the mean of the two arm values.
- No other keyword has a multi-arm rule. For most keywords the general definition still gives one answer when applied to the joined spectrum (for example `MJD-OBS`, `MJD-END`, `NCOMBINE`, `WAVELMIN`, `SNR`, `SPEC_VAL`). For some it does not (`TEXPTIME`, `SPEC_RES`, `LAMRMS`/`LAMNLIN`, `SPEC_ERR`/`SPEC_SYE`, `SPEC_BIN`, the barycentric velocity, `APERTURE`, `PROVi`).
- ESO itself does not make joined X-Shooter spectra. [XSH-RD]: "XSHOOTER products always come separately for each arm: UVB, VIS, or NIR." The only joined multi-arm spectra in the ESO Phase 3 archive are the 100 XQ-100 "joint" spectra (an external collection, made in 2016 under an earlier standard version). Their headers are the only working precedent.
- One automated check matters for the choice of `TEXPTIME`: [CHK] requires `1:TELAPSE >= TEXPTIME (tol=0.1s)` and `TEXPTIME <= MJD-END - MJD-OBS (tol. 0.01s)`. The SOXS arms expose at the same time, so a `TEXPTIME` that adds the two arms would be about twice the elapsed time and would fail.
- One automated check matters for `SPEC_BIN`: [CHK] requires `median(WAVE[i+1]-WAVE[i])/2 <= SPEC_BIN <= median(WAVE[i+1]-WAVE[i])*3/2`. On a joined grid with two bin widths, a mean bin width can fall outside this range.
- Section 8.1 allows a change of bin width ("Different binnings are allowed.") but requires the wavelength array to be "strictly increasing". The overlap can therefore not be two concatenated arm arrays. It must be one array.

## What the standard says about multi-arm products

A search of [SDP] for "multi-arm", "X-SHOOTER", "XSHOOTER", "arm", "collat", "join", "overlap" and "merge" finds these passages. Nothing else in the standard mentions arms.

| Section | Exact wording | Relevance |
|---|---|---|
| 5.1, `EXPTIME` | "When collating spectra from a multi-arm spectrograph, e.g. X-SHOOTER, EXPTIME shall be set to the median exposure time of the individual arms." | The only multi-arm rule. |
| 8.1 | "For spectral data in the spectrum binary table format, the wavelength axis shall have values between TDMIN and TDMAX. Different binnings are allowed. The sequence of elements of the wavelength vector for 1D spectra must be in ascending order, i.e. strictly increasing." | Bin-width changes are allowed. The overlap must give one monotonic array. |
| 8.1 | "we require that each FITS file contain only one science spectrum. [...] 2D spectral frames may be submitted in addition as associated files (see section 4.2)." | One spectrum per file. Per-arm 2D files go in as associations. |
| 8.2, Table 5 | `EXPOSURE`: "The exposure array for e.g. combined spectra of different wavelength coverage." Req-10: "If TTYPEi = 'EXPOSURE', TUCDi shall be set to 'time.duration;obs.exposure' and TUTYPEi to an empty string ('')." | An optional per-bin exposure column, meant for combined spectra whose parts have different exposure times. |
| 10.2, footnote 24 (IFS cubes) | "The requirement of one unique science data HDU applies to the main science file but not to the associated ancillary files, which may contain multiple data extensions, referring for example to individual instrument arms, dispersion orders, or exposures." | Written for cubes, but it shows that ESO accepts one ancillary FITS file with one extension per arm. |
| 5.13 | "science data products may carry over information from the original raw observational data provided the keywords apply to the data product as a whole." | Per-arm raw cards (detector, instrument) do not apply to the joined product as a whole. |
| 5.1, `INSTRUME` | "If any SDP originates from more than one instrument, the value shall be INSTRUME = 'MULTI' [...]" | SOXS is one instrument, so `MULTI` does not apply. |
| 5.1, `FILTER` | "To identify a band-merged data product [...] FILTER = 'MULTI'" | Imaging only. Table 8 marks `FILTER` N/A for `SCIENCE.SPECTRUM`, and [CHK] applies `FILTER=MULTI` "only for OBSTECH=*IMAGE*". |

## Precedents

### ESO X-Shooter IDPs (per arm, not joined)

[XSH-RD] states: "There is one product file for each set of input raw files. XSHOOTER products always come separately for each arm: UVB, VIS, or NIR." and "All XSHOOTER spectra come separately by arm, like the raw data. Typically there are three spectra (from the UVB, VIS and NIR arm). The three products have essentially no spectral gap." The only cross-arm item ESO made was a preview: "for the VIS arm, the second graphical file contains an overview plot of all three arms, discontinued for data taken from December 2025 onwards."

Other X-Shooter IDP facts that bear on the joined design, from [XSH-RD] and the NIR IDP header in [ARCH]:

- Wavelengths are topocentric: "no corrections for barycentric or heliocentric motion have been applied. The corresponding values have been calculated by the pipeline and are stored in the header (HIERARCH.ESO.QC.VRAD.BARYCOR and …HELICOR". The header has `SPECSYS = 'TOPOCENT'`.
- "The header key SPEC_RES contains the nominal resolving power, as derived from the arc lamp calibrations of the same slit width."
- `SNR`: "There is a column “SNR” in the product table [...] Its mean value is written as header key SNR."
- The arm is named by `DISPELEM = 'NIR'` and `HIERARCH ESO SEQ ARM = 'NIR'`. `DISPELEM` is not an SDP keyword.
- `PROVi` lists the raw `XSHOO.*.fits` files of that arm. `HIERARCH ESO PRO REC1 RAWi NAME/CATG` lists the same raw files. The 2D file is `ASSON1` with no `ASSOC1` (it is FITS), and the PNG preview is `ASSON2`/`ASSOC2 = 'ANCILLARY.PREVIEW'`/`ASSOM2`.

### XQ-100 joint spectra (the only joined precedent)

[XQ-RD] describes the method: "We also produced joint spectra of the three arms together, by cutting and pasting the arms at fixed wavelengths (560 nm for the UVB-VIS interface; 1012.5 nm for the VIS-NIR interface). [...] Median flux values in the superposition regions were used to equalize the VIS and NIR arms to the flux level of the UVB arm. The equalization was only meant to avoid occasional jumps across the arms [...]". Each arm had its own velocity bin ("20 km s–1 for the UVB, 11 km s–1 for the VIS, and 19 km s–1 for the NIR"), so the joint spectrum changes bin width at each cut. The joint file has only `WAVE`, `FLUX` and `ERR_FLUX` columns.

Header values of the joint spectrum `J1042+1957_rescale.fits` against its arm spectra [ARCH]:

| Keyword | UVB arm | VIS arm | NIR arm | Joint | What the joint value is |
|---|---|---|---|---|---|
| `EXPTIME` | 3520 | 3320 | 3600 | 3480, comment "(average)" | Mean of the arms. The median would be 3520. |
| `TEXPTIME` | 3520 | (not read) | 3600 | 3600, comment "(maximum)" | Maximum of the arms. |
| `MJD-OBS` | 56710.198500591 | 56710.19856055 | 56710.19859834 | 56710.198500591 | Earliest arm start. |
| `MJD-END` | 56710.24111651 | 56710.24059799 | 56710.24144450 | 56710.24144450 | Latest arm end. |
| `TELAPSE` | 3682.0 | | 3701.9 | 3710.4 | `MJD-END - MJD-OBS` of the joint. |
| `NCOMBINE` | 4 | | 4 | 12 | Sum over arms: every raw file of every arm. |
| `PROVi` | 4 raw `XSHOO.*` files | | 4 raw files | 3 entries: the three arm product `ORIGFILE`s | Points to the arm products, not to raw files. |
| `OBID1` | 661984 | | 661984 | 661984 | Same OB. |
| `SPEC_RES` | 5100 | 8800 | 5300 | 5700, comment "(average)" | Not a plain mean (6400) or median (5300). The method is not documented. |
| `SPEC_BIN` | 0.0282 | | 0.1025 | 0.0554, comment "(averaged)" | Equals `SPEC_BW / NELEM`, the mean bin. |
| `SNR` | 15.2 | 32.0 | 21.2 | 26.0 | Computed on the joint array. |
| `APERTURE` | 1.0 arcsec | 0.9 arcsec | 0.9 arcsec | 0.9333 arcsec, comment "(average)" | Mean slit width. |
| `WAVELMIN`/`WAVELMAX` | 309.98/560.02 | | 982.00/2480.30 | 309.98/2480.30 | Full joint range. |
| `INSTRUME` | `XSHOOTER` | | `XSHOOTER` | `XSHOOTER` | Not `MULTI`. |
| `DISPELEM` | `UVB` | | `NIR` | `UVB,VIS,NIR` | Comma list of arms. |
| `SPECSYS` | `HELIOCEN` | | `HELIOCEN` | `HELIOCEN` | One frame. |
| `FLUXERR` | 10.0 | | 10.0 | 10.0 | Same as arms. |
| `ASSONi` | none | | `ASSON1` = NIR nodded spectrum | none | Joint carries no associations. |

Cautions on this precedent: XQ-100 is an external collection (`P3ORIG = 'EDP'`), made in 2016, before version 8 of the standard. Its `EXPTIME` is the mean of three arms, not the median that version 8 asks for. The ESO TAP query found no other `SCIENCE.SPECTRUM` in the archive that spans both below 600 nm and above 1200 nm, so ESO has no in-house joined-spectrum practice to copy.

## Keyword by keyword

"Rule" quotes the standard. "Multi-arm" says whether the standard gives a multi-arm rule. "Options" lists choices where the rule does not give one answer. The per-arm rulings referred to are those in the map DY-2291.

### `EXPTIME` and `TEXPTIME`

- Rule, `EXPTIME` (5.1): "Total integration time per data element (pixel, spectral bin, spaxel)" and the multi-arm sentence quoted above: "the median exposure time of the individual arms". Each arm value is the per-arm sum of on-target exposure times (DY-2295). With two arms, the median is `(EXPTIME_UVVIS + EXPTIME_NIR) / 2`.
- Rule, `TEXPTIME` (5.1): "Arithmetic sum of the integration time of all exposures included in this product (in seconds)." and "TEXPTIME shall be less than or equal to (MJD-END–MJD-OBS)."
- Multi-arm: only for `EXPTIME`. The standard is silent for `TEXPTIME`.
- Checks [CHK]: `EXPTIME < (TEXPTIME+0.01)`, `TEXPTIME <= MJD-END - MJD-OBS (tol. 0.01s)`, `1:TELAPSE >= TEXPTIME (tol=0.1s)`.
- Conflict: read literally, "sum of the integration time of all exposures" adds the exposures of both arms. The arms expose at the same time, so this sum is about twice the elapsed time and fails the `TELAPSE` check.
- Options for `TEXPTIME`: (a) the maximum of the arm sums (what XQ-100 did; passes both checks); (b) the median of the arm sums, equal to `EXPTIME` (the 5.1 note "in the case where there is stacking and no mosaic mode, the value of TEXPTIME is equal to EXPTIME"); (c) the sum over both arms (fails the checks; not viable unless ESO says otherwise).
- Optional extra: an `EXPOSURE` column (8.2, Table 5, Req-10) that holds the per-bin exposure time, so each arm's true value is in the data.

### `MJD-OBS` and `MJD-END`

- Rule (5.1): `MJD-OBS` is "the start of the earliest exposure if the data product results from the combination of multiple exposures" and `MJD-END` is "the end of the latest exposure". Both "shall include on-target exposures only".
- Multi-arm: no rule, but the general rule gives one answer: the earliest on-target start over both arms and the latest on-target end over both arms. XQ-100 did this.
- Derived: `TELAPSE = MJD-END - MJD-OBS` and `TMID = (MJD-OBS + MJD-END)/2` (5.6.1), checked by [CHK].

### `PROVi`, `NCOMBINE`, `OBIDi`, `PROG_ID`/`PROGIDi`

- `PROVi` rule (5.2.1): "Processing provenance in terms of the list of science files originating this data product. The list of files must not contain duplicates. [...] The PROVi keywords refer to the precursor products and not to the raw data, unless the raw data are the immediate precursors." Req-3: "intermediate products to be referenced must be submitted to ESO not later than the products resulting thereof."
- `PROVi` multi-arm: no rule. Options: (a) the raw `ARCFILE` names of the on-target frames of both arms, in one list (valid if the per-arm spectra are not submitted, so the raw frames are the immediate precursors in the archive); (b) the per-arm Phase 3 spectra by `ORIGFILE` (what XQ-100 did), which needs the per-arm spectra to be submitted to the same collection first or at the same time. The current ruling is that the per-arm merged spectra stay internal, which points to (a).
- `NCOMBINE` rule (5.1): "Number of raw science data files that were combined to generate this data product. Calibration and acquisition data files do not contribute to this count. [...] shall include on-target exposures only." Each arm writes its own raw file, so the general rule gives the sum of on-target raw files over both arms. XQ-100 did this (4 + 4 + 4 = 12).
- `OBIDi` rule (5.1): "the complete set of original observations should be listed" and "a given Observation Block ID must be listed in all the products that are based on this observation". Both arms run in the same OB, so this is the distinct list of OB IDs over both arms, normally `OBID1` only.
- `PROG_ID` rule (5.1): "If any SDP originates from more than one observing run or observing programme, the value shall be PROG_ID = 'MULTI'". Distinct over both arms; normally one value.

### `SPEC_RES`

- Rule (5.11): "Average spectral resolving power [...] In case the spectral resolving power varies within the spectral range, SPEC_RES is defined at the central wavelength of the data (WAVELMIN+WAVELMAX)/2. In case of data combination from different IFUs, slices, exposures (e.g. MUSE) the average spectral resolving power should be assigned to SPEC_RES."
- Multi-arm: no rule. The two sentences point to different answers. For SOXS the central wavelength `(WAVELMIN+WAVELMAX)/2` lies in the NIR arm, so the first sentence gives the NIR value. The second sentence gives an average.
- Options: (a) the R of the arm that contains `(WAVELMIN+WAVELMAX)/2` (the NIR nominal value from the DY-2306 table); (b) the mean of the two arm values; (c) a mean weighted by wavelength span or by number of bins; (d) the lower of the two values (conservative). XQ-100 wrote a value labelled "(average)" whose method is not documented.
- Check [CHK]: `SPEC_RES>0` only.

### `SPEC_BIN`

- Rule (5.6): "The average spectral coordinate bin size in [nm]. It characterises the spectral sampling in the case of the spectrum binary table format for 1D spectra."
- Multi-arm: no rule. 8.1 allows different binnings.
- Check [CHK]: `median(WAVE[i+1]-WAVE[i])/2 <= SPEC_BIN <= median(WAVE[i+1]-WAVE[i])*3/2`. The check uses the median bin width over the whole joined array.
- Consequence: with a fine UV-VIS step and a coarser NIR step, the median bin width is the step of whichever arm has more bins. A mean bin width (`SPEC_BW / NELEM`, as XQ-100 wrote) can exceed 1.5 times that median. Example with assumed numbers only: 25 000 bins of 0.02 nm and 20 000 bins of 0.06 nm give a median of 0.02 nm, an allowed range of 0.01 to 0.03 nm, and a mean of 0.038 nm, which fails. The real bin counts depend on the join point chosen in DY-2328.
- Options: (a) the median bin width (passes the check by construction); (b) the mean bin width (may fail); (c) the bin width at `SPEC_VAL`.

### `LAMRMS`, `LAMNLIN`, `SPEC_ERR`, `SPEC_SYE`

- `LAMRMS` rule (5.6): "The root-mean-square (RMS) of the residuals of the wavelength solution along the wavelength axis [nm]. [...] LAMRMS = sqrt( sum_i(R_i2)/N ) where R_i= residual of wavelength solution for i-th arc line, N = numbers of arc lines = LAMNLIN."
- `LAMNLIN` rule (5.6): "The number of arc lines used in the fit of the wavelength solution."
- `SPEC_ERR` rule (5.6): "Statistical error in spectral coordinate [nm] [...] In most cases, including UVES and XSHOOTER, the formula SPEC_ERR = LAMRMS / sqrt(LAMNLIN) is a good approximation [...] Note that the value of LAMRMS / sqrt(LAMNLIN) is a lower limit to SPECERR."
- `SPEC_SYE` rule (5.6): "aims at capturing the systematic error in spectral coordinate and is expressed in nanometers [nm]. [...] If systematic uncertainties are known to represent an insignificant contribution [...] then SPECSYE can be set to zero."
- Table 8 for `SCIENCE.SPECTRUM`: `LAMRMS` and `LAMNLIN` are O (optional); `SPEC_ERR` and `SPEC_SYE` are R (recommended).
- Multi-arm: no rule for any of the four. The definitions assume one wavelength solution.
- Options for `LAMRMS`/`LAMNLIN`: (a) pool the lines of both solutions: `LAMNLIN = N_UVVIS + N_NIR` and `LAMRMS = sqrt((N_UVVIS·RMS_UVVIS² + N_NIR·RMS_NIR²) / LAMNLIN)`, which is the 5.6 formula applied to all lines; (b) the worse arm's values; (c) leave both out, which Table 8 allows because they are optional.
- Options for `SPEC_ERR`: (a) the pooled `LAMRMS / sqrt(LAMNLIN)`; (b) the larger of the two arm values (conservative). Note that pooling lowers `SPEC_ERR` because `LAMNLIN` grows, although no single wavelength is better known.
- Options for `SPEC_SYE`: the maximum over both arms extends the DY-2307 ruling (maximum over orders) without change.

### `WAVELMIN`, `WAVELMAX`, `SPEC_VAL`, `SPEC_BW`, `TDMIN1`, `TDMAX1`

- Rule (5.6): `WAVELMIN`/`WAVELMAX` "should describe the minimum and maximum wavelengths in nanometers outside of which the described spectrum does not carry any valuable scientific information" and "It must be true that: SPECTRUM.WAVE[1] <= WAVELMIN < WAVELMAX <= SPECTRUM.WAVE[NELEM]". "WAVELMIN / WAVELMAX need to be specified in either air or vacuum wavelength consistent with the definition of the spectral axis."
- Rule (5.6.1): `SPEC_VAL` "is set to (WAVELMAX+WAVELMIN)/2.0"; `SPEC_BW` "is set to WAVELMAX-WAVELMIN".
- Rule (8.1, Table 8 footnote 63): `TDMIN1`/`TDMAX1` are "Mandatory for the TTYPE1 array only (start/stop spectral coordinates)".
- Multi-arm: no rule, and none is needed. The definitions apply to the joined array: `WAVELMIN` from the blue end of UV-VIS, `WAVELMAX` from the red end of NIR. XQ-100 did this.
- Checks [CHK]: `WAVELMIN greater than data[0][0] (tol=bin_max/2)`, `WAVELMAX smaller than data[0][-1] (tol=bin_max/2)`, `TDMIN1 = data[0][0] (tol = bin_max/2)`, `TDMAX1 = data[0][-1] (tol = bin_max/2)`, and `SPEC_VAL`/`SPEC_BW` to `tol=1e-3`.
- Both arms use air wavelengths (DY-2292), so one `TUCD1 = 'em.wl;obs.atmos'` covers the joined array.

### `SPECSYS` and the barycentric velocity

- Rule (5.6): `SPECSYS` is the "Frame of reference for spectral coordinates", one value from Table 2. Table 8: mandatory, primary header.
- Multi-arm: no rule. The standard has no keyword for the velocity applied, and it says nothing on which epoch the correction uses. `HIERARCH ESO QC VRAD BARYCOR`/`HELICOR` are ESO pipeline cards, not SDP keywords.
- What the rule does require: one frame for the whole array. If both arms are corrected, `SPECSYS = 'BARYCENT'` holds for the joined spectrum.
- Options: (a) one velocity at the mean on-target midpoint over both arms, applied to the joined `WAVE`; (b) a per-arm velocity at each arm's own mean midpoint, applied before the join (the DY-2303 ruling run per arm), with the two values recorded in arm-labelled QC cards; (c) as (b) but record only one value. My estimate, not from a source: the two arms of one SOXS observation share a start and differ in end time by at most one NIR exposure, so their mean midpoints differ by minutes and their velocities by well under 0.1 km/s, far below the ±30.4 km/s total from DY-2302.

### `SNR`, `FLUXERR`, `FLUXCAL`, `TOT_FLUX`

- `SNR` rule (5.11): "Median(*) signal to noise ratio per pixel [...] (*)The median is to be computed on all "good" data elements, i.e. excluding any QUAL!=0 data." Multi-arm: no rule, but the general rule gives one answer: the median over all good bins of the joined array. The result leans towards the arm with more bins. XQ-100 computed it on the joint array.
- `FLUXCAL` rule (5.1): "For spectroscopic data, it represents the type of flux calibration, whether 'ABSOLUTE' or 'UNCALIBRATED'". One value per file. The standard is silent on a product in which one arm is calibrated and the other is not. Options: `ABSOLUTE` only when both arms had a response applied (the DY-2294 rule applied to both arms); otherwise `UNCALIBRATED`, or no joined product.
- `FLUXERR` rule (5.1 and Table 8 footnote 43): "the fractional uncertainty on the flux scale expressed as a percentage" with "-1" for `UNCALIBRATED` and "-2 [...] when the flux error cannot be determined". Silent on multi-arm. The NIR scaling to VIS adds an inter-arm term; with `FLUXERR = -2` (DY-2294) nothing changes now. A measured value later would need to cover both arms, for example the larger arm value.
- `TOT_FLUX` rule (5.6): "set to T (true) to indicate that the flux data represent the total source flux". Silent on multi-arm. Scaling NIR to VIS does not correct slit losses, so the DY-2294 ruling `TOT_FLUX = F` is unchanged.

### `INSTRUME`, `PROCSOFT`, `PRO REC1` cards and other identification keywords

- `INSTRUME` (5.1): one instrument, so `INSTRUME = 'SOXS'`. `MULTI`/`INSTRi` apply only to "more than one instrument". XQ-100 kept `INSTRUME = 'XSHOOTER'`.
- `TELESCOP` (5.1): one telescope; `MULTI` does not apply.
- `PROCSOFT` (5.1): "Indicates the reduction software system including its version number used to produce this data product." One value. Silent on multi-arm; both arms are reduced by the same `soxspipe` version, so DY-2300 holds.
- `OBSTECH` (5.1): "Normally, OBSTECH should be adopted from the keyword HIERARCH ESO DPR TECH of the original data." The arms share the template, so one value; the standard is silent if the arms disagree.
- `DISPELEM`: not an SDP keyword. ESO X-Shooter IDPs use it to name the arm (`'NIR'`); XQ-100 wrote `'UVB,VIS,NIR'` for the joint. Optional.
- `HIERARCH ESO PRO REC1 ...` cards: the standard does not define them. 5.13 allows carried-over keywords "provided the keywords apply to the data product as a whole" and asks providers to ensure "user-defined keywords do not conflict with the SDP keywords". A joined product made from two per-arm recipe runs plus a join step has no single `PRO REC1`. Options: (a) `PRO REC1` describes the join step and lists both arm products (or all raw files) as inputs; (b) keep the per-arm recipe chains as `PRO REC1` and `PRO REC2`, the DICD convention for chained recipes; (c) drop per-arm `PRO REC` cards from the joined file and keep them in the per-arm merged spectra. Per-arm `HIERARCH ESO DET`, `INS` and `SEQ ARM` cards fail the "whole product" test unless they are arm-labelled or dropped. This touches the DY-2300 ruling "keep all `HIERARCH ESO` cards".

### `APERTURE`, `RA_ERR`, `DEC_ERR`

- `APERTURE` (5.6.1): "It shall be set to the width of the slit or the diameter of the fiber." One value. Silent on arms with different slit widths. XQ-100 wrote the mean, labelled "(average)". Options: mean, the narrower slit, or the slit of the arm that holds `SPEC_VAL`. If both SOXS arms always use the same width, there is no choice to make.
- `RA_ERR`/`DEC_ERR` (5.6): "could be set to half the width of the slit". The same choice follows.

## Overlap and bin-width changes

- Bin width: 8.1 says "Different binnings are allowed." A joined spectrum may change step at the join. The only constraint is the `SPEC_BIN` median check above.
- Monotonic array: 8.1 says the wavelength vector "must be in ascending order, i.e. strictly increasing", and [CHK] tests "the data in the first column is strictly increasing". The two arms cannot both keep their bins in the overlap. The overlap must become one sequence.
- Equal-length arrays: 8.2 Req-1, "All data arrays in the first row of the binary table shall have the same number of points", and [CHK] "TFORMi must all be the same". Any per-arm column (for example an arm flag or `EXPOSURE`) must span the whole joined array.
- No NaN or Inf: [CHK] "Inf data values not allowed"; the DY-2293 ruling already forbids NaN.
- The standard says nothing on how to treat the overlap: no rule on cut versus average, on which arm wins, on flux scaling between arms, or on how to flag the join. Options from practice: (a) cut at one fixed wavelength and keep each arm's own bins on either side (XQ-100); (b) cut at a per-spectrum point, for example where the two arms' S/N cross; (c) resample both arms in the overlap onto one grid and combine them with inverse-variance weights. A `QUAL` bit could mark bins in the overlap or bins carried by the scaled arm; the standard allows provider-defined bits ("other positive integers greater than 1 may be used to flag data which is bad or dubious for specific reasons").
- The flux scaling of NIR to VIS is not addressed anywhere in the standard. The release description is where ESO expects such processing to be described (4.2: "Phase 3 generally mandates that the data documentation is to be provided in the release description PDF file").

## Associations for per-arm ancillary files

- Rule (4.2): "`ASSONi`: The list of files (FITS / non-FITS) associated to this data product. i is a sequential number starting from 1." "`ASSOCi`: Specifies the product category of the associated non-FITS file given by ASSONi." "`ASSOMi`: Specifies the md5sum of the associated non-FITS file given by ASSONi. Applies to non-FITS files only."
- Req-1 (4.2): "Each ancillary file must be associated to at least one science file and can be associated to more than one science file."
- Req-2 (4.2): "An associated FITS file cannot contain other associations in turn, i.e. nested associations are not supported."
- Table 8 footnotes 31 and 32: `ASSOCi`/`ASSOMi` are "Mandatory for the cases where ancillary files are provided in a non-FITS format"; `ASSONi` is "Mandatory for the cases where ancillary files are provided in association with the scientific data."
- Check [CHK]: "ASSOC* must be in ANCILLARY.[PREVIEW, HARPSTAR, FEROSTAR, README, TAR, MOSSPECTRA.PREVIEW, SRCTBL, 2DECHELLE.TAR, FORS2TAR, ANCILLARY.CUBE.LMV, CRIRESPTAR]". `ANCILLARY.2DSPECTRUM` is a recognised category but not in this `ASSOC` list, as DY-2301 already found. Because `ASSOCi` is for non-FITS files, a FITS 2D file needs only `ASSONi` and carries its own `PRODCATG`; the ESO X-Shooter IDP does this (`ASSON1` = the 2D FITS file, no `ASSOC1`).
- Multi-arm: no rule. Nothing limits the number of associations, so one joined product may list both per-arm 2D files.
- Options: (a) `ASSON1` = the joined preview PNG (with `ASSOC1`/`ASSOM1`), `ASSON2` = UV-VIS 2D, `ASSON3` = NIR 2D, each 2D file with its own `PRODCATG = 'ANCILLARY.2DSPECTRUM'`; (b) one 2D FITS file with one extension per arm (section 10.2 footnote 24 shows ESO accepts per-arm extensions in an ancillary file); (c) per-arm previews in place of, or in addition to, one joined preview. Each 2D file must itself hold no `ASSONi` (Req-2).

## What the standard leaves open

These items have no rule in [SDP], no answer in [CHK] or [FAQ], and no ESO in-house precedent. They go to DY-2328.

1. `TEXPTIME` for arms that expose at the same time: maximum of arms, median of arms, or sum (the sum fails the `TELAPSE` check).
2. `SPEC_RES` for two arms with different R: central-wavelength arm, mean, weighted mean, or minimum.
3. `LAMRMS`/`LAMNLIN`: pool both solutions, take the worse arm, or omit (both optional).
4. `SPEC_ERR`: pooled formula or worse arm. `SPEC_SYE`: maximum over arms is the natural extension of DY-2307.
5. `SPEC_BIN` for a grid with two steps: median, mean or value at `SPEC_VAL`. Only the median is sure to pass the check.
6. Barycentric velocity: one value for the joined spectrum or one per arm, and how to record it.
7. `PROVi`: raw files of both arms, or per-arm products (which must then be submitted too).
8. `PRO REC` cards and per-arm `HIERARCH ESO` cards in the joined file, given 5.13's "apply to the data product as a whole".
9. `APERTURE`, `RA_ERR`, `DEC_ERR` if the arm slit widths differ.
10. `FLUXCAL` when only one arm has a response.
11. Overlap treatment: where and how to join, how to scale, and whether to flag the join or the scaled arm in `QUAL`.
12. Whether to add an `EXPOSURE` column.
13. How the joined product links the per-arm 2D files and previews, and whether ESO accepts `ANCILLARY.2DSPECTRUM` (already on the first Test Mode batch from DY-2301).

Fully determined by the general rules, with no choice left: `EXPTIME` (median, so mean, of the two arms), `MJD-OBS` (earliest), `MJD-END` (latest), `TELAPSE`, `TMID`, `NCOMBINE` (sum of on-target raw files over both arms), `OBIDi` and `PROG_ID` (distinct lists), `WAVELMIN`/`WAVELMAX`/`SPEC_VAL`/`SPEC_BW`/`TDMIN1`/`TDMAX1` (from the joined array), `SNR` (median over good bins of the joined array), `SPECSYS` (one frame for the whole array), `INSTRUME = 'SOXS'` (no `MULTI`), `TOT_FLUX = F`, and a single strictly increasing `WAVE` array with a step change allowed.
