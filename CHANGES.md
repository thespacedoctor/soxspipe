# Release Notes

* **FEATURE**: new `--no-stack` option on `reduce all`, `reduce sof`, `stare`, `stare_std`, `nod`, `nod_std` and `offset` extracts each frame or cycle on its own for one run; standards still write the stacked response first, then also extract each frame (`_F<N>`) or cycle (`_CYCLES`, inspection only) (DY-1723).
* **ENHANCEMENT**: **RESP product format change**: the RESP product no longer holds the legacy table of polynomial coefficients as its first extension; its only table is `RESPONSE`, which is now HDU 1, so any external code that reads HDU 1 by number gets the curve instead of the coefficients and must read `RESPONSE` by name. The `legacy_poly_order` setting is removed and an old workspace `poly_order` is ignored (DY-1719).
* **ENHANCEMENT**: stare, nod and offset now check the `ESO PRO RESP BASIS` card as soon as they find a RESP product in the SOF, so an old-basis RESP fails before any extraction with the same message as flux calibration, which still checks it again (DY-1719).
* **ENHANCEMENT**: Horne extraction no longer computes and resamples the flux density of flat-fielded spectra, which never used it; products are unchanged, except that a flat-fielded `_EXTRACTED_ORDERS` table whose orders do not overlap no longer carries the three `extractedFluxDensity*` columns (DY-1719).
* **ENHANCEMENT**: the offline test suite runs on several CPU cores with `pytest-xdist` (new in the `tests` and `dev` extras): use `-n 4` locally; the required CI `tests` job runs with `-n auto` (2 workers on the hosted runner) and its time limit rises from 20 to 30 minutes; results and coverage match a serial run (DY-1673).
* **ENHANCEMENT**: stacked nod pairs shift the B-A difference onto the A-B difference along the slit using the measured object positions and average them, then extract once, so a science nod makes 1 Horne extraction and a standard 2 instead of 2 and 4; the sky frame is built the same way; falls back to separate A-B and B-A extractions, with a warning, when the positions are not measured or the nod throw is shorter than the Horne aperture; new `NOD COMBINED` and `NOD COMBINE FALLBACK NUM` QCs (DY-1674).
* **ENHANCEMENT**: the SOXS attached-arc test now measures the gap from the end of the earlier exposure to the start of the later, over every raw frame of the stare, nod or offset block, so an arc taken right after a long sequence is attached; NIR exposures end after `ESO DET SEQ1 EXPTIME`, not one DIT (DY-1666).
* **ENHANCEMENT**: SOXS stare, nod and offset skip skyline wavelength tuning when the arc is attached (taken within 7 minutes and 5 degrees of rotator angle of the science frame), because the arc solution is more accurate; new `SKY TUNING` QC records `skylines` or `skipped: attached arc` (DY-1665).
* **ENHANCEMENT**: the three package `__init__.py` files are exempt from the import sorter (`I001`), and a test fails if their import order breaks (DY-258).
* **ENHANCEMENT**: `_dataframe_to_sqlite` (data organiser and recipes) now writes each dataframe in one `BEGIN IMMEDIATE` transaction, so a failed attempt leaves no partial rows, a `replace` never exposes an empty table, and inserts are about 100 times faster on the autocommit workspace connections; both copies share `sqlite_dataframe_writer` (DY-1515).
* **ENHANCEMENT**: fit the response curve as a log-space spline with telluric and stellar-line masks, so the VIS u band is no longer 30-50% wrong; flux calibration reads the RESP `RESPONSE` table and new `RESP` QCs flag poorly constrained curves (DY-1614, DY-1388).
* **ENHANCEMENT**: nod and offset align frames along the slit through the 2D map and stack them before extraction by default; new `stack_before_extraction` setting and `SLIT ALIGN FALLBACK NUM` QC (DY-1612).
* **ENHANCEMENT**: stare SOFs hold all frames of one OB per arm; new `stack_before_extraction` setting (default `True`) stacks them or extracts each on its own (DY-1611).
* **ENHANCEMENT**: skyline wavelength tuning now predicts the SOXS VIS u and g shifts, which have no usable sky lines, from the r and i sky shifts and the change in `ESO ADA ABSROT` since calibration (coefficients in the `ug-wavelength-transfer` block of the VIS detector parameters; leave-one-night-out rms 0.96 and 1.14 px), so VIS stare, nod and offset wavelength-calibrated u/g output changes; `SKY SHIFT O<N>` is now in nm, and `SKY SHIFT O<N> PX` and `UG SHIFT SOURCE` are new; `flag_poor_data` no longer crashes on text QC values (DY-1608).
* **ENHANCEMENT**: `soxs_stare` now falls back to a robust boxcar extraction at the slit centre when no object trace is found, so a spectrum is still produced; new `EXTRACT METHOD` QC and plot warning flag it, and a standard star with no trace fails with a clear message instead of making a response curve (DY-1623).
* **ENHANCEMENT**: Horne extraction now logs one info line per order with the number of slit rows that stopped at the profile clipping iteration limit while still clipping, split into rows refitted to the surviving points and rows that kept their last fit; nothing is logged when every row converges and extraction output is unchanged (DY-1606).
* **ENHANCEMENT**: the not-flattened extraction that stare, nod and offset standards run for the response curve reuses the object trace of the flattened extraction instead of detecting it again, so it no longer writes `_NOTFLAT` object-trace files; output is unchanged (DY-1650).
* **REFACTOR**: the rule "run serially when multiprocessing is off or debug is on" is now one `runSerially` property on `base_recipe` and `create_dispersion_map` instead of six copies, and the unused copy in Horne extraction is gone; behaviour is unchanged (DY-1607).
* **DOCS**: the UV-VIS instrument page lists the VIS order numbers, bands and wavelength ranges, and the QC pages say what `<N>` means (DY-1563).
* **FIXED**: `reduce all` no longer reduces a SOF again once its session status is `pass`, so a stare SOF reduced frame by frame (`stack_before_extraction: False` or `--no-stack`), which never writes the predicted `_EXTRACTED_MERGED.fits`, is no longer overwritten on every run; failed or unfinished SOFs, including one with any product row not `pass`, are still retried and `-x` still reduces everything again; a SOF that passed but lost its product files is skipped until `-x` (DY-1705).
* **FIXED**: flux calibration and the response curve of a stare, nod or offset now stop with an error naming each frame and its `EXPTIME` when the contributing frames (all stacked stare frames, nod A and B frames, offset ON frames) do not share one exposure time, because they divided the counts of the combined spectrum by the first frame's `EXPTIME` without checking; equal exposures give unchanged output (DY-1714).
* **FIXED**: the atmospheric extinction correction now holds the last table value (about 0.01 mag/airmass) above 900 nm instead of returning NaN, so X-Shooter VIS response and flux-calibrated spectra are no longer blanked above 900 nm; SOXS VIS, which ends near 853 nm, is unchanged (DY-1688).
* **FIXED**: the response curve and flux calibration now divide the flat-fielded counts by the detector binning along the dispersion axis (`ESO DET BINY` on VIS), so a binned VIS science spectrum is no longer calibrated 2x too bright against an unbinned standard; the RESP product carries `ESO PRO RESP BASIS = PER_UNBINNED_PIXEL` and flux calibration refuses a RESP product without it (rebuild the response), and flat-fielded merged spectra no longer carry `FLUX_DENSITY_COUNTS`, which counted the dispersion twice (DY-1622).
* **FIXED**: the VIS spatial solution no longer misnumbers an order's multi-pinhole spots by one when an end pinhole is missing: the first-iteration shift is now no shift or one pinhole spacing, whichever centres the pattern between the slit edges; VIS accepts 8 pinholes and an 8-pinhole solution writes its 2D map image (DY-1618).
* **FIXED**: the UVB and VIS extinction correction of a response curve now uses the same mean airmass as flux calibration, over every raw frame of the standard (all stacked stare frames, nod A and B frames, offset ON frames only), not the start and end midpoint of the first frame; a UVB or VIS standard with no readable airmass now stops with a clear error, and NIR response is unchanged (DY-1713).
* **FIXED**: the UVB and VIS extinction correction of a flux-calibrated stare, nod or offset spectrum now uses the exposure-time-weighted mean airmass of every frame in the product (each frame the midpoint of its `ESO TEL AIRM START` and `END`; nod counts A and B frames, offset only ON frames), not the end airmass of the first frame; the value is written to the new `ESO PRO AIRM MEAN` header card, and NIR flux is unchanged (DY-1695).
* **FIXED**: the not-flattened extraction that stare, nod and offset standards run for the response curve no longer tunes its own skyline shifts, which were unstable on NIR; it applies the per-order shifts the flattened extraction of the same frame or nod pair measured, and records no sky QCs, so the reported `SKY SHIFT` QCs now come from the flattened extraction; flattened science output is unchanged (DY-1680).
* **FIXED**: shifting a frame along the slit now keeps the object flux and lands it on its target where the 2D map's pixel area changes along the slit, so the stacked nod shift-and-add no longer adds up to 10% of flux at the NIR K-band red end (2000-2030 nm, response -6.5%) (DY-1679).
* **FIXED**: when nod or offset combines several extracted spectra (each cycle extracted on its own, or a stacked nod that falls back to separate A-B and B-A extractions), the reported `VARIANCE` is now the variance of the per-wavelength median flux, not the median of the input variances: sum(var)/N² for N ≤ 2 spectra and (π/2)·mean(var)/N above, where N counts only the spectra with a value at that wavelength; `SNR` and the SNR QCs rise by up to √N, and flux is unchanged (DY-1675).
* **FIXED**: the UVB and VIS extinction correction now interpolates the La Silla extinction table linearly instead of stepping to the next table value, so the response curve and every flux-calibrated UVB/VIS spectrum no longer jump by up to about 5% at each table wavelength (every 10 nm below 400 nm, every 20 nm above); outside the table the old behavior is kept, and NIR is unchanged (DY-1684).
* **FIXED**: slit alignment no longer raises `TypeError` when no object trace is found in a frame; every order now falls back to the median shift (DY-1674).
* **FIXED**: each session now groups every raw frame when it is created, switched to or prepared, so a session created after the first `prep` gets its SOFs, a late frame reaches every session, and a set with unchanged frames keeps its SOF name and status; `raw_frames.processed` no longer decides what a session groups, and a rejoin resets the status of the active session only (DY-274).
* **FIXED**: SOXS VIS skyline tuning no longer shifts an order to noise when its sky lines are faint: an r or i shift is applied only if the combined photon-noise SNR of its catalogue lines is at least 5.4, r measured alone also needs its second-brightest line at SNR 2.5 or more, and the same check guards the r prediction and the u/g model; the 1.5 px scatter guard is removed because it could never fire; new `SKY LINE SNR O<N>` QC; VIS output changes for one of 25 checked standards, NIR and X-Shooter are unchanged (DY-1663).
* **FIXED**: `soxspipe reduce all -m` workers no longer write the shared log file themselves; they send their records to the parent process, which writes and rotates the log through its own handlers, so simultaneous rollovers can no longer lose records or raise logging errors; worker warnings and errors still reach the terminal and the file handler levels are unchanged (DY-1534).
* **FIXED**: the not-flattened repeat that stare, nod and offset standards run for the response curve no longer writes a second object-trace residual PDF and `OBJECT_TRACE_RES` product row under the flattened pass's names; its trace files take a `_NOTFLAT` suffix and it adds no row, and no science product changes (DY-1234).
* **FIXED**: the X-Shooter UVB and VIS `ThAr.fits` line atlas now names its source column `source` (it was `Source`), so `bootstrap_dispersion_solution: true` no longer raises `KeyError` when the static line list is built; only the column name changed (DY-1605).
* **FIXED**: the NIR response curve no longer oscillates at the telluric band wings. NIR telluric masks are wider and now include the 890-980 nm water band, bins in the raised wing next to each mask are dropped (new `telluric_wing_max_nm`, `telluric_wing_ref_nm` and `telluric_wing_threshold` settings), and NIR uses a penalised spline that crosses masked gaps with little bending (new `smoothing_penalty` setting, NIR 1 with 30 nm knots; NIR `max_gap_nm_max` 270). `RESP BIN CLIP FRAC` no longer counts wing bins and the response plot marks them. NIR response output changes; VIS and UVB are unchanged (DY-1647).
* **FIXED**: SOXS VIS skyline wavelength tuning no longer shifts the r order to a noise peak. i is measured first, r is predicted from the i shift and the change in `ESO ADA ABSROT` (coefficients in the new `r-from-i-transfer` block of the VIS detector parameters; leave-one-night-out rms 0.36 px), r lines are searched only within 0.05 nm of that prediction, and a measured r is combined with the prediction. Unresolved sky-line doublets in the atlas are merged into flux-weighted lines. No VIS order is shifted from fewer than 2 lines. The new `R SHIFT SOURCE` QC says how r was set. VIS wavelength-calibrated output changes; NIR and X-Shooter are unchanged (DY-1645).
* **FIXED**: a recipe forced to fail with a message no longer loses that message when a QC row also fails, so the reason reaches the database `error_message` (DY-1623).
* **FIXED**: nod and nod_std cycles keep the spectrum and order joins of whichever nod position found an object trace, so the SNR and efficiency QCs are no longer skipped when only one position found one (DY-1624).
* **FIXED**: nod, nod_std and offset recipes now fail with a clear "no object trace was found" error when no cycle extracted a spectrum, instead of a pandas `ValueError`; per-cycle stacking now skips cycles with no trace and keeps the order joins of the last cycle that had one (DY-1617).
* **FIXED**: nod and nod-std recipes no longer crash when one of the A or B extractions finds no object trace; the SNR and efficiency QCs are skipped with a warning and the product is still written (DY-1616).
* **FIXED**: nod and offset slit alignment no longer crashes reduction workers with a segmentation fault; pixels where the 2D map gives an infinite source position (zero determinant) were sent to `scipy.ndimage.map_coordinates`, which read memory before the image. These pixels are still flagged bad, and products are unchanged (DY-1615).
* **FIXED**: `soxspipe reduce all -m` no longer hangs for up to 10 hours when a worker process dies; each SOF runs as its own task, a crashed SOF is recorded as failed, the other SOFs that were running are re-run once on their own, and the results are written to the database. A batch that reaches the timeout now stops its workers and records the unfinished SOFs as failed (DY-1615).
* **FIXED**: SOXS VIS orders now use one numbering everywhere (r=1, u=2, g=3, i=4) in `_EXTRACTED_ORDERS` products, QC cards, plots and messages; extraction products numbered them 1=i, 2=r, 3=g, 4=u before (DY-1563).
* **FIXED**: per-order QCs now write their own header cards (`SNR MEDIAN O<N>`, `EFF MEDIAN O<N>`, `X RES MEDIAN O<N>` and the rest of the dispersion family) on every arm; before, every order wrote the same card, so the unsuffixed card held the last order's value, not the all-order value (DY-1563).
* **FIXED**: the SOXS VIS single- and multi-pinhole line lists now mark the order-1 Th lines at 657.77, 658.46 and 659.35 nm as deleted; each is blended with a stronger Th neighbour and sat 0.8 to 2.7 px off the dispersion and spatial solutions in every pinhole. VIS dispersion and spatial solutions, and the products built on them, change slightly (DY-1578).
* **FIXED**: `compute_extractions` and `compute_boxcar_extractions` no longer sort the caller's slices dataframe in place or add extraction columns to it; they work on a copy, with identical output (DY-1481).
* **REFACTOR**: the disabled dispersion-map coefficient cache is removed (`get_cached_coeffs`, which is no longer exported from `soxspipe.commonutils`, its file-name helper and the dead cache write in `write_map_to_file`, which referenced an undefined `hduList`); the x/y fit uses the default coefficients directly and output is unchanged (DY-27).
* **FIXED**: Horne object-profile fitting now refits each slit row on the points it keeps (when they are more than the polynomial degree) when clipping stops at the iteration cap while still rejecting points, so the profile no longer includes the outliers the last pass removed; rows that converge are unchanged, and on synthetic frames with cosmic rays the profile moves by up to 0.07% of its peak (DY-1554).
* **FIXED**: the parameter-tuning pools of `soxs_disp_solution`, `soxs_order_centres` and `soxs_spatial_solution` now run serially when multiprocessing is turned off (`turnOffMP`), as in a bulk-reduction worker, instead of trying to start a child pool from a daemonic worker (DY-1535).
* **FIXED**: `prepare()` now sets the active session ID together with the session path when it reads `.sessionid`, so status reads and writes use the session named in the file (DY-1223).
* **FIXED**: the dispersion and spatial solution now stops with a clear "no usable arc lines survived pre-fit clipping" error, with the detected and clipped counts, when the pre-fit QC clip drops every line, instead of crashing in the QC histogram or reporting a failed fit; both xy residual histograms show a "no lines" note when empty (DY-1551).
* **FIXED**: the published order-centre (`ORDER_TAB`) and object-trace (`OBJECT_TAB`) tables now hold `{axisA}min` and `{axisA}max` as the smallest and largest fitted cross-dispersion pixel on the detector; for UVB, VIS and X-Shooter they held the mirrored values (axis length minus the fit), so these columns change, while the SOXS NIR values, the order list and the `{axisB}` columns do not (DY-1219).
* **FIXED**: the SOXS NIR `line-atlas` now points at a packaged `ArHgNeXe.fits` (a copy of the X-Shooter atlas, the only NIR atlas with the `wave`, `ion`, `source` and `amplitude` columns the static line list reads), and the duplicate VIS `line-atlas` key is removed (`thar_murphy.fits`, the value already in effect); tests check that the detector-parameters yaml files repeat no key, name only packaged static calibrations, and name line atlases with the columns the code reads (DY-1579).
* **FIXED**: continuum tracing outside SOXS VIS now retries with every row sampled when the strided probes find no shift (NaN), instead of stopping with "flux too low" first; the error is raised only if the every-row pass also fails (DY-1459).
* **FIXED**: NIR recipes whose SOF names a master dark now use the read noise measured on that master (`QC RAW RON`) instead of the fixed 7.3 e⁻ from the detector settings, for the prepared-frame errors and the nod/offset difference errors; the stare sky fit also adds the master-dark term, RON·√(1+1/N). A master dark from one raw frame, or an older one without `ESO PRO DATANCOM`, falls back to the settings value with a warning, and a value more than 15% from it is used but logged. NIR errors, sky weights and extraction weights change (DY-1581).
* **FIXED**: `detrend` no longer warns that a SOXS dark has a different exposure time from the frame, because SOXS darks are never scaled. X-Shooter keeps one warning, now only when the times differ by more than the 1 ms header tolerance and the dark is scaled (DY-1553).
* **FIXED**: UVB/VIS recipes whose SOF names a master bias now use the read noise measured on that master (`QC RAW RON`) instead of the nominal header value, for the prepared-frame errors, the bias-pedestal correction and the nod/offset difference errors; the stare sky fit also adds the master-bias term, RON·√(1+1/N). Older masters without `ESO PRO DATANCOM` fall back to the header value with a warning, and a value more than 15% from nominal is used but logged. VIS errors, sky weights and extraction weights change; NIR is unchanged (DY-1458).
* **FIXED**: the `RAW RON` QC of `soxs_mbias` and `soxs_mdark` is now measured from the difference of the first two raw frames divided by √2, so fixed-pattern noise no longer inflates it; a master built from one raw frame writes no `RAW RON`, and master biases and darks now record the number of combined frames in `ESO PRO DATANCOM` (DY-1458).
* **FIXED**: skyline wavelength tuning now measures each matched sky line at the midpoint of its half-maximum crossings, not at its brightest sample, so the flat-topped, dipped sky lines of wide-slit frames are no longer measured on a shoulder. On a VIS 5″ standard the r and i shifts go from +0.10 and -0.10 nm to -0.02 and -0.05 nm; wavelength-calibrated output changes (DY-1599).
* **FIXED**: `soxspipe reduce` no longer crashes with `KeyError: 'workspace-root-dir'`; the packaged advanced settings again default the workspace root to the current directory, which DY-1225 had removed (DY-1601).
* **FIXED**: `soxs_nod`, `soxs_nod_std` and `soxs_offset` now read the UVB/VIS master bias their SOF names and remove the bias-pedestal photon noise from the A and B errors before differencing; VIS nod and offset ERRS, SNR and extraction weights change, and a SOF with no master bias logs a warning (DY-1390).
* **FIXED**: `detrend` now raises an error when a UVB or VIS frame gets both a master bias and a master dark, because UVB/VIS master darks still carry the bias pedestal and the pedestal was subtracted twice. This includes `soxs_mflat` given a master bias with lamp-off or raw dark frames. NIR is unchanged (DY-1580).
* **FIXED**: stare sky-subtraction rolling-window clipping now sorts each order by wavelength itself, so for distinct wavelengths the sky model and clip flags no longer depend on the row order of the input. Output is unchanged for wavelength-sorted input, as production supplies today (DY-1528).
* **ENHANCEMENT**: the response-function product (`RESP_TAB_<ARM>`) now also holds the fitted response curve in a second table extension, `RESPONSE`, with `WAVE` (nm) and `RESPONSE` columns sampled on the extracted standard-star wavelengths; the coefficient table stays the first extension, and flux calibration now reads it by HDU index (DY-1378).
* **FIXED**: the scattered-light background mask now finds each order's neighbours by their position on the detector, not by order number. SOXS orders now grow into half of each inter-order gap (they grew by 2 px only), and X-Shooter no longer masks the whole frame beyond its first and last orders. Background-subtracted output changes in `soxs_mflat`, `soxs_stare` and `soxs_stare_std` (DY-1518).
* **FIXED**: the `INNER ORDER PIX MEAN`/`SUM` QCs and the 38 dispersion-solution residual, FWHM and R QCs now store numbers rounded to 3 and 5 decimals, not formatted strings, in the QC table, FITS headers and database; all three QC header writers now share `phase3.add_qc_to_header`, which leaves a NaN or infinite QC out of the header instead of crashing the write (DY-1259).
* **FIXED**: Horne profile fitting now stops a slit row's clipping loop once a pass clips nothing, instead of always running to the cap, and `horne-extraction-profile-clipping-iteration-count` now gives that many passes, not one fewer. Rows that still clip at the cap get one more fitting pass, so extracted spectra can change by a very small amount (DY-1257).
* **FIXED**: SOXS darks are now always subtracted unscaled, and X-Shooter darks are scaled whenever their exposure time differs from the frame's, where both used to compare whole seconds only (0.4 s and 0.9 s matched). SOXS NIR output changes where the matched dark's whole-second exposure time differs (DY-1263).
* **FIXED**: the `STRUCTX` and `STRUCTY` master-bias QCs are now fitted against a 0 to N-1 pixel axis; the old axis skipped N-1 and ended at N, so the slopes change, most on small frames (DY-1258).
* **DOCS**: the `soxs_stare`, `soxs_stare_std`, `soxs_nod`, `soxs_nod_std` and `soxs_offset` QC tables now list the `N ORDERS SLIT EDGE` QC; only a flux-standard stare reduction that builds a response curve fails on a spill (DY-1507).
* **DOCS**: the recipe and user-manual recipe pages now give each figure and table its own label (QC figures use `<recipe>_qc_fig`), and the `soxs_nod_std` and `soxs_stare_std` QC-metric links now point to their own tables (DY-1508).
* **FIXED**: multiprocessing reductions no longer hang on macOS when a forked worker crashes in Apple Accelerate; the bulk recipe pool now uses the spawn start method (DY-1532).
* **ENHANCEMENT**: `soxspipe reduce all` has a new `-c, --cores <cores>` option that sets the number of worker processes in each recipe's pool when you use `-m`, up to the number of CPUs; without it, pool sizes are unchanged (DY-1377).
* **ENHANCEMENT**: warnings and errors logged during optimal extraction are now printed below the per-order extracted-spectrum QC plot and the order-merged (stare) and stacked (nod, offset) spectrum QC plots, ERROR in red and WARNING in orange, at most 20 lines; nod and offset notes start with their nodding position, e.g. `[A1]` (DY-1559).
* **ENHANCEMENT**: the slit drift QC plot (`SLIT_DRIFT_QC_PLOT`) legend now also gives the science minus multi-pinhole observation time and `ABSROT` difference, assuming 0 deg when the multi-pinhole frame has no `ABSROT`; `read_rotator_angle` moves to `slit_illumination_qc` (DY-1531).
* **FIXED**: the stare sky model no longer carries the wing light of bright objects: each order is refitted without a model of the wings (settings `object_wing_correction`, `object_wing_gate_fraction`). Sky-subtracted output and bright-object fluxes change (DY-1279).
* **FEATURE**: `soxs_stare` and `soxs_stare_std` measure the along-slit offset of the sky image from the master flat for each order, mask the slit ends not lit in both from the sky model (or flag them when sky subtraction is off), and add the `SLIT ILLUM OFFSET O<N>`, `SLIT ILLUM APPLIED` and `SLIT ILLUM N ORDERS` QCs, the `SLIT_ILLUMINATION_QC_PLOT` product and an `ABSROT` spread warning; set by the new `slit-illumination` settings block (on for SOXS, off for X-Shooter), and sky-subtracted output changes (DY-1280).
* **FIXED**: the dispersion and spatial solution x and y polynomials are now solved by exact column-scaled linear least squares, not by scipy `curve_fit`, which stopped short of the minimum; `XY RES MEDIAN` drops about 50% on VIS and 2-7% on NIR DSOL, and no longer swings on tiny input changes; dispersion maps move by a median of 0.1-0.4 px, and more at order ends beyond the last arc line. Downstream output changes (DY-1525).
* **FIXED**: stare sky subtraction now leaves out the pixels where the master flat is not finite, not positive, or below 2% of the 95th percentile of its order, weights the sky fit by the flat-fielded noise, and tests for noisy regions in noise units, so the out-of-band order ends no longer distort the sky model or remove the knots of low-flat regions. Sky-subtracted output changes (DY-1517).
* **ENHANCEMENT**: sky subtraction no longer clips the slit edges of an order that the slit-illumination mask covers. Before, such an order lost another 7% (VIS) or 2% (NIR) of its remaining slit at each end. Other orders are unchanged. Sky-subtracted output changes (DY-1280).
* **FIXED**: SOXS VIS recipes now read the read noise of each readout mode from `ESO DET OUT1 RON` (3.3, 3.5 or 5.4 e-), not the 3.8 e- default, and only raw frames are checked for a read-noise mix, so a calibration product from another readout mode no longer stops a recipe. NIR keeps the 7.3 e- default. VIS ERRS and SNR change (DY-1457).
* **FIXED**: the response function now matches the points kept when its clipping loop stops at `max_iteration` while still rejecting points, and a warning is logged; before, it was fitted before the last rejection and could be off by several percent. Flux-calibrated output changes (DY-1245).
* **FIXED**: the response function fit now drops non-finite raw response samples (from zero extracted flux) before smoothing and fitting, and logs how many; before, one such sample in a NIR standard could give all-NaN coefficients and NaN flux-calibrated output with no error. A response with no finite samples now fails the recipe (DY-1255).
* **FIXED**: sigma-clipped polynomial fits (order centres, object traces, order edges, dispersion solution) that stop at the clipping-iteration cap now refit to the rows they keep, so the coefficients and RES QCs match the final set. Reduction output changes (DY-1246).
* **FIXED**: the rectified slit grid has the set number of rows, centred on the object trace, where it was one row short and half a pixel low; an odd `horne-extraction-slit-length` (default 15) now extracts the full length. Extracted spectra change (DY-1247).
* **FIXED**: merged 1D spectra no longer start or end with a made-up bin of zero flux, zero variance and NaN SNR outside the extracted wavelength range (DY-1248).
* **FIXED**: the MASTER RON QC no longer counts clipped and masked pixels, so its value drops by a few percent (DY-1249).
* **FIXED**: order-centre and object-trace positions on the cross-dispersion axis no longer sit 0.5 px too high. Reduction output changes (DY-1250).
* **FIXED**: the 2D dispersion map is now block-reduced by the frame binning divided by the map's own binning, in the map dataframe and in the stare slit-illumination measurement, so a binned map used with a frame of the same binning no longer crashes with `All arrays must be of the same length`; an unbinned map is unchanged (DY-1260).
* **FIXED**: rectification scales dispersion-solution pixel positions to the frame binning when the two differ, for example a 2x2 frame with a 1x1 solution; output for matching binning is unchanged (DY-1251).
* **FIXED**: rectification leaves flagged pixels out of each cell's flux and variance and rescales the cell to its good-pixel area, so unmasked cells no longer carry part of a hot pixel's value. Extracted spectra change (DY-1252).
* **FIXED**: sky-fit weights in noisy low-flux regions are bounded by the pixel noise, so the sky model is no longer pulled toward zero there. Sky-subtracted output changes (DY-1253).
* **FIXED**: sky subtraction no longer crashes on an order with fewer unclipped pixels than the residual window, no longer turns off pandas SettingWithCopy warnings for the whole process, and logs the measured residual floor instead of a constant 5; a knot rule that never fired and the unused `_rectify_order` and `calculate_residuals` methods are removed, and the SKYMODEL and SKYSUB product descriptions say their ERRS is the object frame's uncertainty (DY-1286).
* **FIXED**: extraction (stare, nod, offset) now records the running recipe in the `soxspipe_recipe` column of its products-table rows, instead of always `soxs-stare` (DY-907).
* **FIXED**: `soxs_offset` gives each cycle and the stacked pair their own object-trace table, residual PDF, extraction products and `_ONOFF` frame (suffixes `_1`, `_2`, `_STACK`), because real offset frames have no `CUMOFF Y` and every cycle used to overwrite one file, so none overwrites another (DY-906).
* **FIXED**: every recipe now runs from a list of frame paths or a directory as well as a `.sof` file, instead of failing with `AttributeError`; the night date is the earliest `MJD-OBS` minus 15 hours, the QC `sof_name` is `<recipe>_<ARM>_<night>.sof`, and input with no frames or no `MJD-OBS` raises a clear `ValueError` (DY-90).
* **FIXED**: the sky image-comparison QC plot rotates and flips its three panels, and labels their axes, with the arm's `rotate-qc-plot` and `flip-qc-plot` detector parameters, as the other QC plots do (DY-699).
* **FIXED**: reading a preserved database (QC history and session statuses restored on a refresh) now goes through a scratch copy, so no file in `backups/` is created or changed; `open_backup_read_only` is now a context manager (DY-273).
* **FIXED**: `data_organiser.session_switch` now updates `sessionId` as well as `sessionPath`, so a later `prepare()` on the same instance reads and writes the status column of the session switched to (DY-910).
* **FIXED**: `detrend` no longer counts the electrons of the bias pedestal (about 650 e- in the VIS arm) as signal in the photon-noise term of a frame's ERRS; before, ERRS was 3-6 times too large in sky-dominated pixels of VIS stare, spatial-solution, order-centre, dispersion-solution and master-flat reductions, so extracted-spectrum variance falls and SNR rises by roughly 8-20% (DY-1283).
* **FIXED**: sky subtraction sets QUAL on every pixel it zero-fills and no longer writes a 0 uncertainty, and the SKYSUB bad-pixel QC ignores these flags (DY-1284).
* **FIXED**: the sky fit now weights each pixel by 1/σ from a noise model of the previous fit (σ² = model + RON²) instead of a rolling noise estimate and, in VIS, the pixel's own flux; this removes a bias that put the VIS sky model too high at skyline peaks. The DY-1253 noisy-region reweight is removed and the sky model is no longer clipped at zero. Sky-subtracted output changes (DY-1282).
* **FIXED**: `soxs_stare` now gives the un-flattened sky-subtracted frame an uncertainty of the flat times the flat-corrected uncertainty, where it kept the flat-corrected one; the `_NOTFLAT` extraction of a flux standard and its response curve change, and other stare products do not (DY-1391).
* **FIXED**: continuum tracing outside the SOXS VIS arm (SOXS NIR and X-Shooter order centres, master flats and object traces) no longer crashes with `IndexError` when its first sampling pass finds too few points and the wider retry finds the trace; the retry's trace width was stored in place of the pixel table (DY-1294).
* **FIXED**: the sky b-spline fit now clips outliers on both sides after each fit, at `bspline_fitting_residual_clipping_sigma` times the noise model, and refits with the final knots until no new pixel is clipped (at most 5 refits); cosmic-ray residuals and cold pixels no longer pull the sky model, and clipped pixels are flagged `flagged_bspline_clipped`. Sky-subtracted output changes (DY-1281).
* **FIXED**: the dispersion solution no longer crashes with a `TypeError` when fewer arc lines are left than the polynomial has coefficients, at the start or after clipping, or when the degrees reach zero; it lowers the polynomial degrees and retries, and raises `ArithmeticError` only if no reduced-degree fit succeeds (DY-1256).
* **FIXED**: `clip_and_stack` now leaves a pixel out of the frames that flag it, where only a pixel flagged in every frame was left out, and takes the uncertainty from the frames that contributed at the pixel; `ignore_input_masks=True` now differs from `False`, and the caller's frame masks are no longer cleared. Master bias and master dark change at pixels clipped in some frames, and stacked nod and offset frames change at cosmic-ray pixels (DY-1254).
* **ENHANCEMENT**: astropy is upgraded from 6.1.2 to 7.2.2 and `uncompresspy` 0.4.1 is added, so astropy can open LZW-compressed `.fits.Z` frames. Reduction output changes slightly, because astropy 7 no longer computes float32 sigma-clipping statistics with bottleneck, which lost precision (DY-694).
* **FEATURE**: `soxs_mbias` writes `<product>_BIAS_DISTRIBUTION_QC_PLOT.pdf`, which compares the pixel flux histograms of the earliest raw bias frame and the master bias, with the RON values (DY-267).
* **FEATURE**: extraction (stare, nod, offset) writes `<template>_SLIT_DRIFT_QC_PLOT<nodding sequence>.pdf`, which plots the object slit position against wavelength before rectification, with per-order fits and residuals; it is not written for the `_NOTFLAT` re-extraction (DY-270).
* **FEATURE**: extraction (stare, nod, offset) adds the `SLIT DRIFT` (arcsec; fitted slit position 5% in from the red end of the reddest order minus that 5% in from the blue end of the bluest) and `SLIT ZERO WLEN` (nm; bluest wavelength where a fitted order curve crosses slit position 0, never across the gap between orders) QCs to the FITS header, and the slit-drift QC plot shows the drift in its legend and marks the zero crossing; neither QC is written for the `_NOTFLAT` re-extraction (DY-1506).
* **ENHANCEMENT**: optimal extraction logs an ERROR naming each order where the object spills over the slit edge, and adds the `N ORDERS SLIT EDGE` QC; when a flux standard spills, `soxs-stare` fails the recipe so its response curve is not used to flux-calibrate other frames (DY-802).
* **ENHANCEMENT**: the `soxs_mbias` bias-distribution QC plot now draws filled, translucent histogram bars with no outline, and shows its mean, median and RON values in a table under the plot instead of free text (DY-909).
* **ENHANCEMENT**: `.fits.Z` frames are now read in place by `prep` and the recipes, so the `uncompress` command is no longer needed; when `X.fits` and `X.fits.Z` both exist the compressed frame wins and `X.fits` is deleted (DY-310).
* **ENHANCEMENT**: continuum-trace fitting (object traces, order centres, order edges) now works on its own copies of the pixel tables, so it no longer writes into pandas slices of a caller's table and runs with `mode.chained_assignment = "raise"`; reduction output is unchanged (DY-1294).
* **ENHANCEMENT**: `prep` no longer turns off pandas SettingWithCopy warnings for the whole process, so every `soxspipe reduce` run shows them again; order-table unpacking, optimal extraction and the dispersion-map fit now work on their own copies of their tables, so they no longer write into pandas slices and run with `mode.chained_assignment = "raise"`; reduction output is unchanged (DY-1462).
* **ENHANCEMENT**: `soxspipe prep` reads only the primary header of `.fits.Z` frames instead of decompressing each whole file, so indexing new frames is about 20x faster (DY-1509)
* **ENHANCEMENT**: `build_sof_files` refreshes a `calibration_matches` table that the `cal_<type>` views read, from one `CALIBRATION_RULES` rule per calibration type in `calibration_matching.py`, instead of cross-joining every incomplete product with every upstream product; the first SOF-building pass on a 141,391-frame workspace falls from about 1024 s to about 50 s with the same SOFs, ties at equal distance are broken by upstream `mjd-obs` then `filepath`, and a workspace from an older soxspipe makes `prep` and `reduce` stop with a message to run `soxspipe prep --refresh` (DY-1513).
* **ENHANCEMENT**: raw `.fits.Z` frames are now decompressed with the system `gzip`, about 20x faster than `uncompresspy` (DY-1510)
* **REFACTOR**: `soxs_mflat` keeps each lamp's products in one record keyed by lamp tag instead of parallel lists, and the UV stitch reads the D and QTH lamps by name; a stitch without one of them now raises a `ValueError` naming the missing lamp (DY-131).
* **REFACTOR**: the disabled slit-illumination and tilt corrections in sky subtraction (`cross_dispersion_flux_normaliser`, `adjust_tilt`) and the unused `slit_illumination_order` setting are removed; output is unchanged (DY-1280).
* **FIXED**: the global order-centre, object-trace and order-edge polynomial fits now use a column-scaled linear least-squares solve instead of `curve_fit`, so they no longer depend on floating-point noise; before, an unrelated code change could move the SOXS NIR order edges and so change the master flat and every later NIR product. Reduction output changes slightly (DY-1222).
* **FIXED**: the lower order-edge fit in `soxs_mflat` now skips rows with no lower-edge position, as the upper-edge fit already did, so an order whose edges could not be measured no longer fails the master flat (DY-1222).
* **FIXED**: the order-centre and object-trace residual plots, the dispersion-map residual plot and the sky-subtraction QC plots now embed their images at full resolution, so the markers line up with the frame when you zoom in. In the residual and dispersion-map plots the QC and settings tables are on separate rows, and on arms that are flipped but not rotated, such as SOXS NIR, the order-centre markers and fit lines are no longer one pixel off their row (DY-1218).
* **FIXED**: in the `soxs_mflat` order-edge QC plot (`_ORD_LOC.pdf`), the flat image is now embedded at full resolution, so the edge markers line up with the order bands when you zoom in. The QC and settings tables are on separate rows and no longer overlap. On arms that are flipped but not rotated, such as SOXS NIR, the edge markers and fit lines are no longer one pixel off their row (DY-1217).
* **FIXED**: `soxs_mflat` records exactly one `MFLAT` product row, for the final master flat; the per-lamp master flats are still written but no longer get `MFLAT`, `MFLAT_DLAMP`, `MFLAT_QLAMP` or `MFLAT_DOME` rows (DY-117).
* **FIXED**: optimal (Horne) extraction no longer inflates the flux up to about 2.6x when the object core sits on masked slit-edge rows; masked pixels in slit rows with sparse local support now get no profile weight (DY-802).
* **FIXED**: a SOF listing `X.fits.Z` no longer reads `X.fits` instead or fails, and prepared frames from `.Z` inputs no longer crash on a `.Z` output name (DY-310).
* **FIXED**: optimal (Horne) extraction no longer inflates the flux by up to 1e6x when the object trace sits on an order edge; the object profile now gives no weight to pixels off the order (DY-752).
* **FIXED**: the master flat now flags partly illuminated order-edge pixels, which flat-fielding amplified about 6x; set `order-edge-min-flat-fraction` (default 0.5, 0 turns it off) and re-make existing master flats to gain the flags; these pixels are reported as `N ORDER EDGE` and are not counted in the `COLDPIX` QC (DY-753).
* **FIXED**: `soxs_stare` now shows the raw-stack surface plot only in debug mode (`-d`), where before it opened a window on every run and waited until the window was closed.
* **FIXED**: a poor b-spline sky fit (FITPACK `ier` of 10 or more, including 50) on the first iteration now raises a named `ValueError` instead of `UnboundLocalError`, and after a revert the returned knots match the reverted spline (DY-602).
* **FIXED**: noisy-region knot pruning in `subtract_sky` removes only the knots that bound noisy pixels, where it always dropped the reddest knot as well (DY-601).
* **FIXED**: `subtract_sky` anchors the order ends with sigma-clipped, clamped straight-line fits over each end window instead of the whole-order median, and always uses the fitted spline, so the flux-error ratio on a perfectly flat sky is now 0 (DY-593).
* **FIXED**: pixels released when object clipping is reset (more than 85% clipped, or the VIS retry) now rejoin the sky fit (DY-595).
* **FIXED**: aggressive object masking in `subtract_sky` records the full slit range of each object, including objects that start at the slit edge (DY-596).
* **FIXED**: the sky QC plot now applies `rotate-qc-plot` and `flip-qc-plot` like the other QC plots, so SOXS NIR image panels are no longer upside down relative to them, and the image axis labels name the right detector axes (DY-698).
* **FIXED**: the clipped-pixel panel of the sky QC plot no longer hides parts of the order on rotated arms such as SOXS VIS; the order-outline layer is now rotated like the other layers.
* **FIXED**: sky QC plots place pixels correctly for y-dispersion arms such as SOXS NIR and no longer crash on non-square frames (DY-597).
* **FIXED**: the Read the Docs build installs the package with `pip` instead of `python setup.py install`, as the repository no longer has a `setup.py`.
* **FIXED**: a raw frame added after `prep` grouped its set now rejoins that set on the next plain `prep`: the set keeps its SOF name, its status is reset in every session, and its stale product is deleted (DY-263).
* **FIXED**: `quicklook_image` now shows the inter-order mask from a dispersion map instead of silently dropping it, and no longer modifies the caller's frame (DY-81).
* **FIXED**: `quicklook_image` no longer crashes with `RecursionError` when it is given a `CCDData` with `ext=False` and a dispersion map overlay, as in `soxspipe stare -xd`.
* **FIXED**: `quicklook_image` called without `settings` now skips the dispersion map and skyline overlays with a warning instead of raising (DY-290).
* **FIXED**: the real-data workflow's `NPY_DISABLE_CPU_FEATURES` now names numpy's `AVX512_*` dispatch groups, so AVX512 runners no longer run the AVX512_SKX/CLX kernels (DY-40)
* **FIXED**: before each b-spline sky fit, `subtract_sky` drops knots that would leave a knot interval with no sample (duplicate knots, knots outside the sampled wavelength range, and knots with no sample since the previous knot), which used to make FITPACK reject the fit (`ier=30`) and revert (or raise on the first iteration); dropped and duplicate knots no longer count as new knots in the stop check (DY-697).
* **FIXED**: each recipe now gets its own scratch directory under `<workspace>/tmp` from `tempfile.mkdtemp` instead of a random number, so two recipes running at the same time can no longer share one and delete each other's files (DY-49).
* **FIXED**: `prep` and directory SOFs now read their directory listings in name order, so which frames a capped `prep` pass takes, and which map wins when an arm has two, no longer depend on the filesystem (DY-48).
* **FIXED**: the order-centre QC plot from `detect_continuum` keeps each order's own colour in every panel and no longer drops the last order from the lower panels when an earlier order's fitted trace is off the detector (DY-604).
* **FIXED**: a recipe whose QC rows still fail to write to the database after seven attempts now raises the original error, so the SOF is marked `fail`, where before the rows were dropped silently (DY-603).
* **FIXED**: in a named session, an out-of-range QC value now marks the SOF `fail` in that session's status, where before the session showed `pass` and the failure was written to the `base` session instead (DY-272).
* **FIXED**: VIS sky-subtraction object clipping now restarts its convergence check after a retry, so an equal first pass no longer ends the loop; its retries stop lowering the sigma limit at 1.5 with a warning, and a non-finite `percentile_clipping_sigma` raises a clear error; a bright VIS order skips sky subtraction on its own (data passed through, sky model 0 and flagged in QUAL) instead of turning sky subtraction off for the whole frame; and the arm test is an exact match. Sky-subtracted output changes where a retry or a bright order occurs (DY-1285).
* **FIXED**: a `prep` interrupted after marking a product SOF `complete = -1` no longer leaves that SOF to be marked complete by the next `prep` without its calibrations being checked (DY-1513).
* **FIXED**: the 2D dispersion map reader and the response function now close the FITS files they open, so a reduction no longer leaves them to the garbage collector with a `ResourceWarning`; products are unchanged (DY-1709).
* **TEST**: the real-data gate also checks the NIR nod-standard response at 2020 nm, an extrapolated point inside the DY-1647 telluric masks, held to 5% of the centre recorded from six identical CI runs (DY-1712)
* **TEST**: the real-data nod-standard, stare and NIR offset baselines are re-recorded from six identical CI runs after DY-1674, DY-1675, DY-1679, DY-1680 and DY-1684; both nod-standard `SNR MEDIAN` values rise by about √2 (VIS 207.34 to 292.84, NIR 50.06 to 70.14), and the response, efficiency and flux-calibrated medians move by 0.6% or less (DY-1690, DY-1691)
* **TEST**: the real-data NIR response, stare and offset flux-calibration baselines are re-recorded from six identical CI runs after the NIR response fit change; the NIR response median rises 2.8% and both NIR flux-calibrated medians rise 0.6% (DY-1648)
* **TEST**: the real-data stare, NIR offset and calibration baselines are re-recorded from six identical CI runs after DY-1525, DY-1254, DY-1279 and DY-1518; VIS `XY RES MEDIAN` falls about 50% and the VIS stare sky model median falls from 61 to 9 e-. The stare and offset checks now report every value before they assert (DY-1521, DY-1524)
* **TEST**: the real-data NIR offset and stare baselines are re-recorded from CI after the astropy 7.2.2 upgrade (DY-694)
* **TEST**: the real-data VIS and NIR stare baselines are re-recorded from six identical CI runs after the sky-model rework; the VIS sky model 90th percentile falls from 875 to 550 (DY-1514)
* **TEST**: characterization tests raise subtract_sky.py coverage from 67.7% to 95.4% (DY-40)
* **TEST**: the subtract_sky end-anchoring characterization test pins the spline derivative instead of a flux-error-ratio sum that differed between AVX2 and AVX512 CPUs (DY-40)
* **TEST**: the real-data workflow now also reduces a VIS and an NIR stare SOF, and `tests/real_data/test_stare_baseline.py` checks the stare product set, merged and flux-calibrated spectra, sky model levels and QC values within bands (DY-696)
* **TEST**: the real-data workflow now runs only on pull requests into `main` and on manual dispatch; it no longer runs on pull requests into `develop` or on a weekly schedule (DY-840)
* **TEST**: `tools/check_changelog.py` rejects a bullet repeated within one release block, in pre-commit and CI (DY-851)
* **TEST**: CI now fails when any file is not formatted with the pinned `ruff format`; the whole repository was reformatted once to start clean (DY-850)
* **TEST**: the real-data gate now fails when any recipe in its reduction chain records a QC `fail` in the `quality_control` table, or when the table is empty; an expected failure must be listed with a reason in `tests/real_data/test_qc_verdict.py` (DY-908)
* **TEST**: the real-data gate now holds all 23 calibration SOFs and both standard-star `nod` SOFs to baseline QC values within bands, checks each standard star's response curve on a fixed wavelength grid, and fails when it reduces a SOF with no baseline or stops reducing one (DY-257)

## v0.18.0 - September 30, 2026

* **FEATURE**: added the shared QC, product, timestamp and plot helpers — `utcnow_string`, `append_qc`, `append_product` and `save_qc_plot` in `soxspipe/commonutils/toolkit.py`, plus thin `add_qc`/`add_product` delegators on `base_recipe` that add nothing but a fallback to `self.recipeName` and `self.dateObs`.
* **FEATURE**: added a `[tool.ruff]` house-rule configuration to `pyproject.toml` (line length 120, camelCase-friendly naming ignores, `max-statements = 50`) and applied the 705 safe ruff fixes ahead of it. Import sorting is deliberately withheld from the three package `__init__.py` files, whose import order works around a circular import.
* **ENHANCEMENT**: widened the `eff-median` QC acceptable ranges for `soxs-nod-std` in `soxs_default_settings.yaml` (upper limits 0.30 to 0.60 and 0.2 to 0.5).
* **ENHANCEMENT**: `image_transformer` now has separate `zoomFactorSlit` and `zoomFactorWavelength` sub-sampling factors in place of a single `zoomFactor`, so the slit and wavelength axes of the rectified image can be oversampled independently.
* **ENHANCEMENT**: the missing-calibrations report now says why each calibration is missing (`failed QC`, `failed run`, `not yet reduced`, `not observed` or `no match`), and `get_blocking_calibration_sets` names the raw sof and failure message responsible, printed as a second table by `soxspipe prep` and `soxspipe reduce all`; for `mbias`/`mflat` the reason and blocking sof are scoped to the science set's own binning and readout speed, matching the real `cal_mbias`/`cal_mflat`.
* **ENHANCEMENT**: Raise the whole-package branch-coverage threshold to 70% and add broad deterministic unit and integration coverage for pipeline utilities and recipes.
* **ENHANCEMENT**: refine skyline-based wavelength-shift calibration
* **ENHANCEMENT**: tune skyline detection, iterative matching, clipping, and VIS-order shift fallback behavior.
* **ENHANCEMENT**: update calibration data and default dispersion, spatial, and sky-subtraction thresholds.
* **ENHANCEMENT**: using numba JIT to speed up image rectification
* **REFACTOR**: Widen `structx` and `structy` acceptance limits for SOXS master-bias quality checks.
* **REFACTOR**: Raised the supported Python version floor from 3.11 to 3.12.
* **REFACTOR**: huge refactor of entire code-base. Unit-test coverage is now approaching 90%.
* **REFACTOR**: Limit multiprocessing pool size for stare, nod, and offset recipes to avoid memory issues.
* **REFACTOR**: Update predicted paths for standard response products (now looks for RESP function instead of merge spectra when deciding if the data has been reduced yet).
* **REFACTOR**: correct rectified-image diagnostic orientation and add local CodeGraph metadata exclusions.
* **REFACTOR**: optimising converting dispersion map to pixel arrays
* **REFACTOR**: vectorising calculation of rectification weights
* **DOCS**: adding to doc FAQs
* **FIXED**: a recipe marked as failed because QC values are outside the acceptable limits is now logged as a warning instead of an error.
* **FIXED**: `session_switch` now points the shared `sof_map` database view at the session being switched to.
* **FIXED**: a workspace database locked or busy in another process now makes `data_organiser._get_or_create_db_connection` refuse with `DatabasePreservationError` after the first failed attempt, instead of retrying 50 times (up to about four hours with the 300 s busy timeout)
* **FIXED**: `soxs_mflat` NIR verification now rejects mixed input image types.
* **FIXED**: `base_recipe.qc_ron` now records a MASTER RON QC row when it measures the master frame itself, instead of returning the value without recording it.
* **FIXED**: `soxspipe prep --refresh` and the automatic rebuild of a database that will not open no longer lose quality-control history: the old `soxspipe.db` is first saved to a new `backups/` directory in the workspace root (a verified single-file snapshot, or a byte-for-byte copy if it will not open), and its `quality_control` rows are restored into the rebuilt database before the QC acceptable-range checks run; if the database cannot be saved it is left in place, nothing is rebuilt, and `soxspipe prep` exits with status 1. The automatic rebuild of an unopenable database now restore each session's per-SOF pass/fail status from the backup when the SOF still holds the same frames; SOFs whose frames differ are queued again.
* **FIXED**: `image_transformer` wavelength bins are now evenly spaced in detector pixels along each order's trace instead of evenly spaced in nm, removing an 11–21% blue-to-red flux ramp within every order
* **FIXED**: `image_transformer` slit bins now use each order's arcsec-per-pixel scale measured from the dispersion map instead of a hardcoded value, so the Horne aperture set by `horne-extraction-slit-length` is now that many detector pixels.
* **FIXED**: `image_transformer._determine_rectified_image_boundaries` now fits per-order polynomials of slit position vs wavelength from trace points, so the rectified window follows the object trace along each order (the object trace seems to drift through slit-positions as wavelength increases).
* **FIXED**: `base_recipe._write(..., maskToZero=True)` wrote `1` into masked pixels instead of `0`, contradicting its argument name, docstring and log message, all of which already said `0`.
* **FIXED**: `soxspipe/commonutils/detect_continuum.py`'s `get()` reports a failed continuum fit by returning `None` as its product path, both when the sampled trace detects under 10% of the order and after five failed fitting attempts, and `soxs_order_centres._fit_order_centres` did not check for this, passing the `None` straight to `os.path.basename` so a user whose data could not be fitted saw a bare `TypeError: expected str, bytes or os.PathLike object, not NoneType` traceback instead of a recipe failure.
* **FIXED**: `soxs_order_centres.verify_input_frames` had two rejection messages that never fired as written. The NIR mixed-image-type message was assigned to `erorr`, a misspelling, so it was never raised.
* **FIXED**: `soxs_mflat.mask_low_sens_pixels` sampled the median order flux from the wrong axis when the dispersion axis is y. The `axisA == "y"` branch unmasked the vertical band `interOrderMask[l:u, b]` but then read the flux with `frame.data[b, l:u]`, the horizontal band the `x` branch reads, so the returned `medianOrderFluxDF` came from a different region from the one just unmasked, and the read raised `IndexError` whenever a sampled `b` fell outside the frame's row range, which a frame with more columns than rows can reach. The read is now `frame.data[l:u, b]`, matching the mask.
* **FIXED**: `fit_bspline_curve_to_sky` no longer recomputes the knot-quantile expression that just failed, which had both lost the original traceback and silently retried a transient failure. A failed tuning iteration in `soxs_disp_solution` now propagates rather than leaving the tuning grid quietly recording nothing.
* **FIXED**: the reduction is no longer order-dependent on the CPU it runs on. `create_dispersion_map` sorted its line table on `wavelength` alone, which is not a total key because wavelength repeats across orders and slit positions. Pandas defaults to an unstable quicksort and NumPy dispatches float sorts to SIMD kernels whose permutation of tied rows varies by CPU microarchitecture, so the same input produced different row orders on different runners, changing summation order in the polynomial fits. The sort now uses the full `(wavelength, order, slit_index)` key with `kind="stable"`, and all 20 `sort_values` call sites in the reduction path were made stable.
* **FIXED**: the sigma-clipping decision in `create_dispersion_map.fit_polynomials` no longer flips on last-bit floating-point noise. A residual a fraction of an ULP either side of the threshold changed which lines survived, which changed the next `curve_fit`, which moved the merged spectrum's red end by 12 bins. Residuals are now quantised to 9 decimals as they are computed, well below any physically meaningful difference and well above the noise floor.
* **FIXED**: QC and product rows are now built with `pd.DataFrame([row])` instead of `pd.Series(row).to_frame().T`, so appending a row no longer upcasts the whole `qc_value` column to `object` and numeric QC values keep their dtype through to the FITS headers.
* **FIXED**: the "['skyFlux'] not in index" reduction fails.
* **FIXED**: Detrend master-dark flat frames once, preserve slit-position data in dispersion-map images, and skip optional efficiency products when no estimate is available.
* **FIXED**: Guard merged spectrum QC sky plotting against non-positive sky counts
* **FIXED**: Handle empty and prepared-extension file inputs safely, omit NaN Phase 3 QC headers, and invoke FITS decompression without a shell.
* **FIXED**: Improve flux standard matching via 3 FITS header keywords
* **FIXED**: Restore editable installs by using the supported default setuptools-scm parser for `v`-prefixed tags.
* **FIXED**: Skip efficiency-product creation when the optional response efficiency estimate is unavailable.
* **FIXED**: fixing database to collect the correct response curves
* **FIXED**: pinning the normalisation factor using by ccdproc for flat correction
* **FIXED**: prevent dispersion-map transformations from running outside debug mode and handle invalid sky-plot statistics.
* **FIXED**: `soxs_mflat` order median fluxes now ignore masked pixels, and a fully masked order raises a `ValueError` naming the order (#527).
* **FIXED**: `toolkit.spectroscopic_image_quality_checks` now defaults a missing `WIN_BINX` or `WIN_BINY` header to 1 binning for non-NIR arms instead of failing with an unbound variable, and always uses 1 for NIR arms (DY-78, #523).
* **TEST**: Add a verified, opt-in real-data NIR-offset acceptance workflow with immutable archive inventory checks and approved scalar baselines.
* **TEST**: Add an isolated synthetic test foundation, branch-coverage ratchet, and repository-owned required CI workflow.
* **TEST**: Add synthetic file-contract coverage for SOF inputs, FITS frame preparation, and Phase 3 products.
* **TEST**: Add synthetic master-flat, order-edge, and response workflow coverage, and raise the branch-coverage ratchet to 46%.
* **TEST**: Characterize dispersion conversion, order tables, stacking, subtraction, extraction, response, flux calibration, continuum fitting, and order-edge detection with deterministic synthetic data.
* **TEST**: Characterize recipe construction, validation, orchestration, product metadata, and failure behavior with isolated synthetic inputs.
* **TEST**: Limit the legacy integration workflow to pull requests targeting `main` and a weekly scheduled run.

## v0.17.4 - July 9, 2026

* **FIXED**: standard extraction fail when no sky have been subtracted.
* **FIXED**: Ignore flux standards missing from the static library when creating SOF files
* **REFACTOR**: Broaden master bias structure QC acceptable ranges.
* **REFACTOR**: Increase decimal precision for many QC metrics.

## v0.17.3 - July 9, 2026

* **FIXED**: issue with stare efficiencies being reported lower than expected (issue with unflattening sky-subtracted data).
* **FIXED**: database to collect the correct response curves for science objects.
* **FIXED**: matching against STD akas.
* **FIXED**: other bugs squashed

## v0.17.2 - July 1, 2026

* **FEATURE**: Added `image_transformer` module to help with image rectification.
* **ENHANCEMENT**: Add order-level rectified image handling for flux, variance, sky, masks, wavelength, and object profile data.
* **ENHANCEMENT**: Split extraction into dedicated mask generation, profile fitting, extraction computation, and debug plotting helpers.
* **ENHANCEMENT**: finding of standard name in FITS headers more robust
* **ENHANCEMENT**: Refreshed static skyline list resource.
* **ENHANCEMENT**: added a 'write_fits_table_to_disk' function to consolidate code writting FITS binary tables to file.
* **ENHANCEMENT**: added a phase3 module (yet to be completed)
* **ENHANCEMENT**: There are now separate recipes for spectroscopic standard star reductions (`soxs_nod_std` and `soxs_stare_std`). Closes #422.
* **ENHANCEMENT**: Improve sky residual modelling and debug plots
* **ENHANCEMENT**: added debugging plots to reveal individual order sky lines
* **ENHANCEMENT**: Logging parallel recipe runtimes
* **ENHANCEMENT**: Expand extracted-order QC plots to show optimal flux, boxcar flux, SNR, and sky flux panels.
* **ENHANCEMENT**: Sky is now plotted in the nodding and offset data extraction QC plots
* **ENHANCEMENT**: Reintroduced writing of the modelled scatter background as a QC FITS image product (alongside PDF) 
* **ENHANCEMENT:** adding quickstart data link to the docs (thanks Markus)
* **ENHANCEMENT:** moving to github actions (away from jenkins)
* **REFACTOR**: Move shared detector, skyline, map, binning setup etc into `base_util` to remove duplicate code.
* **REFACTOR**: clarifying stacked vs frame clipping in settings and docs
* **REFACTOR**: Normalised sky residuals by error.
* **REFACTOR**: Simplify residual floor and skyline flagging logic in determine_residual_floor.
* **REFACTOR**: adjusted laCosmic settings so as not to clip bright skylines
* **REFACTOR**: Disable post-stack clipping when combining normalised mflat frames
* **REFACTOR**: Reuse the shared skyline loader in quicklook and merged-spectrum QC plotting.
* **REFACTOR**: Replaced the `skyModelFrame` parameter in `horne_extraction` with `subtractedFrame` for extracting a sky spectrum in nodding mode -- the subtracted (B or A) frame is now used directly instead of a separate sky model frame
* **REFACTOR**: Tune default clipping iterations and sigma thresholds for master bias and dark recipes.
* **REFACTOR:** updated SNR calculation on extracted spectra. The variance using is now from the propagated errors instead of the SNR being derived from the extracted spectrum itself.
* **REFACTOR:** allowing for NODDING to continue if no trace found in single AB sequence.
* **REFACTOR:** forcing sky-subtraction to be switched off for relatively bright sources in VIS
* **REFACTOR:** information tables are now printed to the terminal in markdown format
* **REFACTOR:** moving to `pyproject.toml` and away from older `setup.py` installation
* **REFACTOR:** renaming `master` branch to `main`
* **REFACTOR:** updating NODDING setting to catch continuum traces closer to the edges (spatially) of the orders 
* **DOCS**: updating docs to include offset mode.
* **DOCS**: updating docs to reflect the decoupling of standard star and science objects in the stare and nodding recipes.
* **FIXED**: soxspipe version checking does not preemptively uncompress files. Thanks @kalabartainaf
* **FIXED:** A fix to the uncompression of last few .Z frames in a batch. Thanks @kalabartainaf
* **FIXED**: Update UTC timestamp handling to avoid deprecated `utcnow` usage.
* **FIXED**: fixing output directory cl-switch bug.
* **FIXED**: Respect post-stack clipping setting during frame combination 
* **FIXED**: Correct top and bottom frame masking in background subtraction
* **FIXED**: Add separate per-frame clipping settings for pre-stack sigma clipping and apply stacked-frame clipping after combination.
* **FIXED**: pass `order` argument to `initial_sigma_clipping` call so per-order clipping is correctly parameterised
* **FIXED**: Replace the zero readout-noise error with a clearer corrupted-frame diagnostic.
* **FIXED:** fixing reduction of VIS binned object data - response curve binning verification was too strict

## v0.17.1 - May 1, 2026

* **FIXED:** case sensitivity issue with the name of the `skylines.fits` file on ubuntu.

## v0.17.0 - May 1, 2026

* **FEATURE:** Offset mode has been integrated into the pipeline.
* **ENHANCEMENT:** Print a table of exported files when running the `raw` command.
* **ENHANCEMENT:** Add new parameters to fit and clip noise from the sky-model spectrum.
* **ENHANCEMENT:** Sky spectrum added to the stare extraction QC plots.
* **ENHANCEMENT:** Add SNR values to the bottom panel of the extracted spectra QC plots.
* **ENHANCEMENT:** Sky flux and variance added to the extracted spectra FITS tables (stare so far).
* **REFACTOR:** Don’t attempt sky-subtraction on *very* bright sources in VIS.
* **REFACTOR:** If continuum fitting fails on a single NOD cycle, still continue and attempt to fit other cycles.
* **REFACTOR:** Add FITS verification.
* **REFACTOR:** Improve sky-subtraction.


## v0.16.0 - April 20, 2026

* **FEATURE:** Added a multiprocessing option when reducing all frames.
* **FEATURE:** Added guardrails to recipes using a new `qc-acceptable-ranges` settings block for each recipe. Lower and upper acceptable ranges can now be set for each QC metric, and a forced failure is triggered if a QC falls outside its defined range.
* **FEATURE:** You can now export all the raw frames needed to reduce a science SOF file (see the `soxspipe raw sof` command and the FAQs in docs).
* **ENHANCEMENT:** Added `ESO TPL NEXP` to the data-organiser database to allow incomplete OB datasets to be filtered out.
* **ENHANCEMENT:** Added an information line at the end of a multiprocessing batch reduction run with details of passes and failures.
* **ENHANCEMENT:** Added checks for corrupt frames so the pipeline warns about them before attempting to reduce the data without the offending files.
* **ENHANCEMENT:** Added database queries to enforce QC guardrails.
* **ENHANCEMENT:** Added guardrail values to the QC database table.
* **ENHANCEMENT:** Added gain to the data organiser.
* **ENHANCEMENT:** Added a `GOODLINES FRAC` QC for the proportion of good, unclipped lines in a pinhole frame.
* **ENHANCEMENT:** Added global and order-specific median efficiency QCs (`EFF MEDIAN`).
* **ENHANCEMENT:** Added SNR QC metrics for all extracted spectra, including global and order-specific values (`SNR MEDIAN`).
* **ENHANCEMENT:** Added `PINHOLE COUNT MIN` to the spatial solution QCs.
* **ENHANCEMENT:** Made `X DIFF SD` and `Y DIFF SD` use the absolute values of the pixel shifts, since VIS pseudo-bands move in different directions.
* **ENHANCEMENT:** Added instrument temperature and AFC numbers to the data-organiser database.
* **ENHANCEMENT:** Added many more FITS keywords to the `raw_frames` table in `soxspipe.db` and the data organiser.
* **ENHANCEMENT:** Added median pinhole FWHM and resolution to FITS header QCs.
* **ENHANCEMENT:** Added more pinhole checks to the spatial solution.
* **ENHANCEMENT:** Added order numbers and error messages to the QC and products tables in the database.
* **ENHANCEMENT:** Added Python files to the allow-list of file extensions that are not moved to the misc directory when `prep` is run.
* **ENHANCEMENT:** Added separate settings for nodding standards versus objects.
* **ENHANCEMENT:** Allowed the response function code to work in multiprocessing mode.
* **ENHANCEMENT:** Creating useful views on the `soxspipe.db` (sorting products and raw frames by recipe).
* **ENHANCEMENT:** Refreshed the pipeline session after each reduction to update the status of each product in the database.
* **ENHANCEMENT:** Updated some QC names and added several more.
* **ENHANCEMENT:** Updated the `prep` command so that the `-r` flag no longer removes previous ERROR logs.
* **REFACTOR:** Renamed QC metrics as suggested in RIX 119.
* **REFACTOR:** Added try/except handling when creating QC and product directories to avoid race condition issues.
* **REFACTOR:** Changed the grid size from 2x2 to 3x3 when determining line shifts in NIR dispersion solutions.
* **REFACTOR:** Cleaned up standards data reduction.
* **REFACTOR:** Optimised the spatial solution for speed.
* **REFACTOR:** Preserved `detector_x` and `detector_y` as the original line-list values in the multipinhole dispersion solution. The X DIFF SD and Y DIFF SD QCs now show the true spread in pixel shifts.
* **REFACTOR:** Reduced the size of dispersion map image files.
* **REFACTOR:** Reduced the size of MBIAS and MDARK products by using `float32` instead of `float64`.
* **REFACTOR:** Reduced the typical `mflat` recipe memory footprint from approximately 7 GB to approximately 2 GB.
* **REFACTOR:** Removed cruft from the data-organiser module.
* **REFACTOR:** Removed previous guardrails that forced recipe failures when QC values strayed beyond set limits. This is now handled with the `qc-acceptable-ranges` setting.
* **REFACTOR:** Removed VIS darks from the data organiser so they are not automatically written to SOF files.
* **REFACTOR:** Sped up Horne extraction.
* **REFACTOR:** Sped up the SNR calculation for extracted spectra.
* **REFACTOR:** Updated extracted spectra in ASCII files to use wavelength in Angström instead of nm, as requested by the consortium. FITS table extractions remain in nm.
* **REFACTOR:** Updated the detect continuum algorithm to fail if fewer than 30% of lines are detected, increased from 10%.
* **DOCS:** Added `horne-extraction-profile-poly-order` and `use_lacosmic` to the `soxs-nod` and `soxs-stare` parameter documentation.
* **DOCS:** Added a note in the documentation about the `soxs_sof_map.yaml` file.
* **DOCS:** Added a QC plot for the response curve and efficiency.
* **DOCS:** Added a section to the documentation about reducing data in multiprocessing mode.
* **DOCS:** Added acceptable ranges to all recipe QCs.
* **DOCS:** Added details of how the efficiency is calculated.
* **DOCS:** Added a note about QC metric guardrails in the `qc-acceptable-ranges` section of the settings.
* **DOCS:** Added details of the `max_iteration` and `poly_order` response function parameters to the `soxs-nod` and `soxs-stare` documentation.
* **DOCS:** Added documentation for targeting and reducing individual science SOF files.
* **DOCS:** Added `GOODLINES FRAC`, `FWHM PIN SD` & `PINHOLE COUNT MIN` to the QC descriptions for dispersion solutions.
* **DOCS:** Added `EFF MEDIAN` & `SNR MEDIAN` to the QC descriptions for standard star reductions.
* **DOCS:** Updated some QC names and added several more.
* **DOCS:** Added stare and nodding QC plots.
* **DOCS:** Noted that the pipeline can now reduce entire datasets in parallel using multiprocessing.
* **DOCS:** Noted that the pipeline can now specifically target and reduce individual SOF files.
* **DOCS:** Replaced Xshooter example QC plots and images with SOXS versions.
* **DOCS:** Updated the quickstart guide to use SOXS data instead of Xshooter.
* **FIXED:** Fixed a bug in the SOXS stare recipe settings.
* **FIXED:** Fixed a substantial bug in `mflat` where the original input flats were being modified during normalisation instead of a copy. This bug was introduced during an attempt to reduce memory usage.
* **FIXED:** Fixed logging for the watch command.
* **FIXED:** Fixed spectroscopic QC measurements where the mask was incorrect.
* **FIXED:** Fixed the inter-order mask used to isolate orders in the `mflat` recipe.
* **FIXED:** fixing slit width matching for std response curve files and science frames.
* **FIXED:** science frames are now being reduced again with the reduce all command.

## v0.15.2 - January 26, 2026

* conda release issue fixed

## v0.15.0 - January 26, 2026

* **FEATURE:** Add a `calculate_rolling_snr` function to robustly compute the SNR ratio across the entire spectrum wavelength range. Algorithm is the [presented here](https://esahubble.org/static/archives/stecfnewsletters/pdf/hst_stecf_0042.pdf).
* **FEATURE:** Added a `--refresh` flag to the prep command. Using this flag will completely refresh the `soxspipe.db` database, delete all ERROR logs, and rebuild all SOF files from scratch.
* **FEATURE:** Show progress when uncompressing `.Z` files.
* **ENHANCEMENT:** A single QC file is written for each recipe. This is a plain text file containing all computed QC values for the recipe.
* **ENHANCEMENT:** Added true NIR gain and RON to detector settings file.
* **ENHANCEMENT:** Checks are in place to determine if darks, pinhole, order trace, and flat frames are 'good'. If not, an error is forced and the pipeline removes the faulty calibration files from downstream SOF files.
* **ENHANCEMENT:** Dark scaling warning is thrown once per recipe as opposed to once per raw frame to be detrended (resulting in multiple duplicate warnings).
* **ENHANCEMENT:** Integrating flux-calibration into the dataorganiser.
* **ENHANCEMENT:** Standard star response function now recorded in the products table of the DO database.
* **REFACTOR:** Cleaned up the text printed after `soxspipe prep`.
* **REFACTOR:** Complete refactoring of the data-organiser. Now much faster, and data-organiser dictionaries are located in one YAML file.
* **REFACTOR:** Corrected validation of response curve files by not checking binning with NIR data (previously this was causing needless fails).
* **REFACTOR:** If the ICS is in simulation mode, data is ignored (e.g. NISE is forced to stay at the 5" slit).
* **REFACTOR:** Massive refactoring of `create_dispersion_map` to make it more readable and manageable.
* **REFACTOR:** Move to adjusting the NIR predicted line-locations in detector quadrants rather than order by order.
* **REFACTOR:** Move to detect continuum on the unflattened frame in stare mode (to mirror what has been done in nodding).
* **REFACTOR:** Moving management of matplotlib backend to a single location in the `base_recipe` class (easier to manage).
* **REFACTOR:** Moving module level imports into classes and methods to reduce memory footprint.
* **REFACTOR:** Raw frames are now sorted into UTC dated folders instead of Chilean time (to match La Silla setup).
* **REFACTOR:** Removed extra stare reductions at the end of `soxspipe reduce all`.
* **FIXED:** Added `jinja2` to `setup.py`.
* **FIXED:** Fixed an issue with the macOS matplotlib backend tripping up Ubuntu reductions (again).

## v0.14.1 - October 20, 2025

* **FEATURE:** Dome flats are now recognised by the data-organiser and are combined into master flats.
* **FEATURE:** Efficiency estimates now also done in stare-mode (alongside nodding).
* **FEATURE:** Efficiency estimates now written out as a FITS binary table.
* **FEATURE:** SNR calculated on merged spectra.
* **ENHANCEMENT:** Spectrum plots now have 3 panels for linear flux, log flux and SNR vs wavelength.
* **ENHANCEMENT:** Order join locations now plotted in extraction QC plots.
* **ENHANCEMENT:** If there is no or little flux in the raw flat frames, an ERROR is raised, and the resulting master flat is removed from the reduction cascade.
* **ENHANCEMENT:** The data organiser will selectively find the nearest dome master flat. If none exists, it will look for a QTH master flat.
* **ENHANCEMENT:** After each reduction run, the pipeline will now report a table of frames that could not be reduced due to a lack of calibration frames. Closes #389 
* **ENHANCEMENT:** Pushing VIS extractions to the very edge of the detector to increase order overlap regions and overall wavelength coverage.
* **REFACTOR:** optimisations to stare and nodding receipes (now ~2x faster)
* **REFACTOR:** A new algorithm has been added to 'dither' the line-lists for each order (both arms) to attempt to successfully find lines if a significant shift is seen between the detector coordinates in the static line-list and those observed on the detector.
* **REFACTOR:** Adjusted all line-lists to reflect the formats now seen in La Silla.
* **REFACTOR:** Nodding parameters changed for robustness. Fixes #390
* **REFACTOR:** Data organiser needs to rebuild the cache of SOF it needs to process after each recipe's specific set is completed. Closes #388 
* **REFACTOR:** Data-organiser and recipes can now determine if CONAD or GAIN keyword values should be used as e-/ADU. Basically `gain = max(CONAD,GAIN)`. Fixes #385
* **REFACTOR:** u-band object tracing is now very robust. Fixes #360
* **REFACTOR:** clipping the NaN flux value before modelling the sky (was making ~1 in 300 images fail)
* **REFACTOR:** if a standard star observed in NODDING mode is not in the static calibration std-star library, the pipeline will issue a warning rather than failing.
* **REFACTOR:** la-cosmic parameters changed to differ for stare vs nod. This improves object tracing in nodding mode. Fixes #390
* **FIXED:** discontinuity issue when stitching stare orders together. Fixes #392
* **FIXED:** master flat not getting used during nodding reductions. 
* **FIXED:** numbering of orders in the unmerged VIS spectrum plots was wrong. Orders are now labelled correctly.
* **FIXED:** tighter clipping when stacking the flats helps reduce noise in the VIS u-band.
* **FIXED:** standard star names can be resolved from OBJECT *or* TARGET keywords
* **FIXED:** The data organiser is now differentiating between stare, nod, and offset modes again. Fixes #383 
* **FIXED:** If one of the cross-dispersion profiles used in the Horne extraction is completely masked, the recipe can skip this profile and continue.
* **DOCS:** Added a warning for the user:

    > A typical user will not need to use the watch command. It is included in the pipeline to allow for automatic/autonomous reductions as data is flowing into a directory. Most users will only need to use the reduce command.


## v0.13.6 - August 30, 2025

* **FIXED:** Final extraction was switched off accidentally. Fixed now.

## v0.13.5 - August 29, 2025

* **REFACTOR:** do not scale SOXS darks if exptimes are not equivalent (dark does not scale linearly). Instead, do not subtract dark and warn the user. fixes #376
* **REFACTOR:** huge speed gain on data organiser prep command
* **REFACTOR:** ignore SOXS Deut flats by default (too many lines seen in the orders)
* **REFACTOR:** improved order edge detection (especially VIS u-band). fixes #378
* **REFACTOR:** Increasing the default SOXS slice height to detect continuum in VIS
* **REFACTOR:** much improved dispersion SOXS VIS solution
* **REFACTOR:** object tracing improved
* **ENHANCEMENT:** added a few more keyword to the data organiser database
* **ENHANCEMENT:** If a previously failed recipe is rerun and succeeds, update SOF files to re-add the product. closes #380
* **ENHANCEMENT:** many updates to the data-organiser to start handling ACQ camera data
* **ENHANCEMENT:** object clipping during VIS sky-subtraction now adjusts some variables dynamically to get an improved result
* **ENHANCEMENT:** The pipeline can now detect and fail gracefully whenever fewer than 9 pinholes are found in any given order. This will stop the pipeline from propagating bad dispersion solutions if there is a problem with the data. fixes #379


## v0.13.4 - April 29, 2025

* **REFACTOR:** small performance gain in prep command
* **ENHANCEMENT:** pipeline no longer attempts to rerun failed recipes (unless forced)
* **ENHANCEMENT:** recipes are now run in chronological order

## v0.13.3 - April 28, 2025

* **FIXED:** soxspipe prep command with the `--vlt` option was having issues when the workspace being prepared was on a separate mount from the data

## v0.13.2 - April 24, 2025

* **FEATURE:** response curves and efficiency plots generated for flux standards in nodding mode
* **REFACTOR:** optimisation of spatial solution recipe
* **REFACTOR:** improved XSH UVB mflat when binning is used
* **REFACTOR:** optimisation of master flat code
* **REFACTOR:** optimisation of background light fitting code
* **REFACTOR:** some optimisations of continuum detection code
* **REFACTOR:** some optimisations of the horne extraction
* **REFACTOR:** optimisation of stare mode (>2 times faster)
* **REFACTOR:** if the object trace produces residuals with a mean > 10 pixels, then a fitting fail is forced.
* **REFACTOR:** Created a `utility_setup` tool to create the QC and product directories needed for recipes and utils (to stop duplicating code)
* **ENHANCEMENT:** Added a lookup table for standard star aliases sometimes used in the FITS header naming.
* **ENHANCEMENT:** Added fill_value="extrapolate" to `interp1d` as some standard stars contained extracted flux outside of the database-stored absolute flux.
* **FIXED:** issue with data-organiser tripping on duplicate files
* **FIXED:** Small bug fixes

## v0.13.1 - April 10, 2025

* **ENHANCEMENT:** can reduce all SOXS UVVIS binnings
* **ENHANCEMENT:** strategically added a few more lines to the SOXS UVVIS line lists.
* **ENHANCEMENT:** added a parameter in the settings file to turn plotting of the sky-model QC plot on or off (off by default as this is a time expensive plot).
* **ENHANCEMENT:** included a parameter in the settings file to turn on-frame sky-subtraction (stare-mode) on or off be default
* **ENHANCEMENT:** added a '--vlt' flag to the prep command. If used, the pipeline will opt to use the `/data/raw` and `/data/reduced` folders found in a typical VLT environment workstation.
* **REFACTOR:** changed parameters sent to daostarfinder to better detect pinhole lines
* **REFACTOR:** data-organiser can now read the exptime from the ACQ camera images (found in the `HIERARCH ESO DET3 EXPO TIME` keyword)	
* ... and much more

## v0.13.0 - April 2, 2025

* **ENHANCEMENT:** added settings and refactored code for SOXS nodding .. first SOXS on sky data can now be reduced.
* **REFACTOR:** rebuilt VIS line lists
* **REFACTOR:** changed some internal file names
* **REFACTOR:** retuned some SOXS VIS recipe parameters
* **FIXED:** gain is now getting read from the FITS headers for non-NIR frames
* **FIXED:** fixed the gain keyword in the SOXS keyword lookup table
* **FIXED:** don't shift slices for 2nd iteration of object trace detection for SOXS VIS (u&g shift in opposite direction from r&i)
* **FIXED:** a SQLite query containing "fail". It is now 'fail' in single quotes.
* **FIXED:** an issue with mixed slit width in SOXS stare mode. We had a PAE setting to allow order centre traces to be reduced in stare mode, but this setting was tripping up true stare-mode data.
* ... and much more

## v0.12.3 - February 25, 2025

* **FIXED:** fixing string literal bug resulting from new sqlite3 release

## v0.12.2 - February 21, 2025

* **ENHANCEMENT:** Pipeline name and version added to product headers.
* **ENHANCEMENT:** recipe ID added to product headers.
* **ENHANCEMENT:** added comments to all settings in the SOXS default settings file
* **ENHANCEMENT:** there are now 2 rounds of object trace detection. The first round helps to locate the slit-position and a typical standard-deviation of the object profile, and the second round uses this information to preform a better measurement and fit of the trace. This is improving the fitting for faint objects.
* **REFACTOR:** adjusted default 'slice' lengths used to detect a continuum in nodding frames.
* **REFACTOR:** many adjustments to improve robustness and extractions
* **FIXED:** recipe command was getting written to the bottom of the log file of the previous recipe instead of the top of the current recipe log

## v0.12.1 - February 10, 2025

* **ENHANCEMENT:** resolution plots added to the `soxs-spatial-solution` QC PDF (if a through-slit arc frame file is provided).

## v0.12.0 - January 31, 2025

* **ENHANCEMENT:** Major updates to the documentation following PAE.

## v0.11.10 - December 19, 2024

* **ENHANCEMENT:** Raw frames & their categories are reported in the FITS headers of the recipe products.
* **ENHANCEMENT:** Intermediate calibration products are reported in the recipe product FITS headers alongside their MD5 hash (DATAMD5)
* **ENHANCEMENT:** Recipe settings and values are reported in the FITS headers of the recipe products.
* **ENHANCEMENT:** The recipe command is now printed before AND after execution to the terminal and the recipe log file.
* **REFACTOR:** The length of SOF and product filenames have been compressed so they can be reported within the FITS header character limit.
* **REFACTOR:** reduced products are now stored in the ESO-compliant nested-folder scheme `/reduced/YYYY-MM-DD/`
* **REFACTOR:** QCs are also stored in a similar nested-folder scheme `/qc/YYYY-MM-DD/`
* **FIXED:** All FITS file products pass fitsverify checks

## v0.11.9 - November 25, 2024

* **REFACTOR:** updating the data-organiser to fix a performance issue.
* **FIXED:** Allow the object trace polynomial orders to dynamically change if a fit to the continuum is not found (ala the order-centre tracing)

## v0.11.8 - November 21, 2024

* **FEATURE:** SOXS UVVIS line-list (first draft) now ships with the code.
* **FEATURE:** pipeline can now 'watch' a folder and automatically reduce raw data added to it. This 'watch' feature can also be run as a system daemon.
* **ENHANCEMENT:** Python 3.11 is the mimimum python version now
* **ENHANCEMENT:** pinning the main packages used by soxspipe
* **ENHANCEMENT:** now recording 'ESO ADA ABSROT END' in the soxspipe.db database
* **ENHANCEMENT:** adding a check to see if the continuum fit is good in each order ... remove bad orders (SOXS VIS only so far)
* **ENHANCEMENT:** improving SOXS extractions (using order centre traces with flat lamps)
* **REFACTOR:** Data-organiser can now work with new SOXS DPR keywords.
* **REFACTOR:** changing order centre clipping to mean and std (not median and mad)
* **REFACTOR:** updated NIR line list after SOXS format change
* **FIXED:** a few bugs

## v0.11.6 - October 15, 2024

* **FIXED:** a regression in the detect continuum code

## v0.11.4 - October 9, 2024

* **ENHANCEMENT:** LaCosmic is run on images prior to source extraction to help remove cosmic ray hits.
* **ENHANCEMENT:** Dichroic region is now clipped from the Xshooter UVB order merge spectrum
* **REFACTOR:** Improved source continuum tracing making nodding reductions more robust.
* **FIXED:** a bug where binning was not taken into consideration when reading the detector format before performing source extraction.

## v0.11.2 - October 7, 2024

* **ENHANCEMENT:** New raw data is to be 'streamed' into the workspace's root folder. When 'soxspipe prep .' is run, the new data is found and filed automatically into the correct `/raw/YYYY-MM-DD/` folder. Adding new raw frames directly into the `raw/YYYY-MM-DD/` nested folder structure is also possible.
* **ENHANCEMENT:** "STD,TELLURIC" frames now recognised in stare mode.
* **REFACTOR:** Made big speed gains in the `order_to_image` method. This speeds up the `soxs_spatial_solution` dramatically.
* **REFACTOR:** exptimes recorded in the SQLite database are now rounded to 2 decimal places. There are occasions where exposure times in the FITS headers are given to 4 dp, and a set of (e.g. flat) frames have exposure times that differ by 0.0001 secs. The pipeline divided these frames into two sets when they should have been grouped together. This was causing failure on some mflat recipes. 
* **FIXED:** Fixed numpy 2.0 compatibility issues so soxspipe can now run with numpy v2.0 and greater
* **FIXED:** Instrument = "SHOOT" now recognised as Xshooter data
* **FIXED:** bug in the detect continuum code that would fail to fit a gaussian in the cross-dispersion slices if NaN values were present.
* **DOCS:** complete update of the SOXS documentation. See https://soxspipe.readthedocs.io
* **DOCS:** 4 installation methods are now reported in the documentation (anaconda is not required to install the pipeline).[main version](https://soxspipe.readthedocs.io/en/main/)user_manual/installation.html

## v0.11.1 - August 15, 2024

* **REFACTOR:** interpolated wavelength resolutions set to match that of the Xshooter pipeline in order-merged spectra (NIR were the same, but now UVB and VIS arms also match Xshooter).
* **FIXED:** bug in the Horne extraction causing artificial 'undulations' in the extracted spectra.
* **FIXED:** bad pixel treatment in Horne extraction (was severely affecting NIR extraction in particular)
* **FIXED:** filename case-sensitivity issue when working on a case-sensitive file system
* **FIXED:** issue where the prep command was looking for a settings file in "~/.config/soxspipe"

## v0.11.0 - July 21, 2024

* **FEATURE:** ascii exports
* **FEATURE:** nodding mode
* **FEATURE:** allowing order-trace frames to be reduced in stare mode for PAE
* **FEATURE:** pipeline is now reducing soxs UVVIS data robustly.
* **FEATURE:** adding code to help tune the pipeline. Add the setting `tune-pipeline: True` to the setttings file and run a recipe command.
* **ENHANCEMENT:** pipeline parameters (in default settings file) are now optimally tuned from UVVIS up to flats and NIR up to spatial solution.
* **ENHANCEMENT:** This release includes many robustness updates

## v0.10.2 - April 23, 2024

* **ENHANCEMENT:** the calibration lamp name is now added to the sof filenames, and hence the product file names
* **ENHANCEMENT:** file summary now shows which calibration lamps are used
* **ENHANCEMENT:** adding bad-pixel maps for SOXS detectors (currently blank)
* **ENHANCEMENT:** pipeline can select default settings for either soxs or xsh (previously there was only one default settings file)
* **FIXED**: SOXS VIS darks are now getting split by EXPTIME in data-organiser

## v0.10.1 - April 11, 2024

* **FEATURE:** the data-organiser has been 'plumbed' to work with SOXS data (will now work with Xshooter or SOXS data).
* **ENHANCEMENT:** clipping of entire MPH set based on their combined RMS scatter from their predicted locations. MPH sets with large scatter are consider poor and removed before polynomial fitting.
* **ENHANCEMENT:** option added to relevant recipes settings to allow toggling of fitting and subtracting of intra-order scattered background light (`subtract_background`)
* **REFACTOR**: added arm and lamp to QC plot titles.
* **REFACTOR** recipe settings can now be set independently for each arm.
* **REFACTOR:** fitting of the scatter background light is now much more robust.
* **REFACTOR:** The scattered light background images are now saved as QC PDFs instead of FITS frames.
* **FIXED**: fixed issue where logs were getting duplicated.
* **FIXED**: the scaling and stitching together of the UVB D2 and QTH lamp flats.

## v0.10.0 - February 20, 2024

* a `bootstrap_dispersion_solution` has been added to the advanced settings. It this is set to True, the pipeline will attempt to bootstrap the initial dispersion solution (using the static line list) with lines from a line-atlas. The line-atlas contains more lines and lines with more precise wavelength measurements.
* **FEATURE**: a new 'reducer' module and terminal command replace the old `_reduce_all/sh` script. This allows the data-organiser to dynamically self-correct if a recipe fails.
* **ENHANCEMENT:** robustness fixes and updates.
* `pinhole_fwhm_px_min` and `pinhole_fwhm_px_max` settings added to `soxs-spatial-solution`. Detected pinholes with a FWHM below/above these values get clipped.
* **FIXED**: The bad-pixel mask from the noise map of the mbias frame is now injected mbias product. The Xshooter UVB electron trap is now clearly visible in the master bias quality extension.
* `mph_line_set_min` setting added to `soxs-spatial-solution`. Full multi-pinholes sets (same arc line) with fewer than mph_line_set_min lines detected get clipped.

## v0.9.9 - January 24, 2024

* **FIXED**: bug fix logger

## v0.9.8 - January 19, 2024

* **FIXED**: bug fix in collecting settings files from the default location


## v0.9.7 - December 7, 2023

* **ENHANCEMENT:** the instrument name is now included in the SOF & product filename.
* **ENHANCEMENT:** setting bad pixels to zero in sky-subtracted product frames.
* **FIXED**: blocking filters now taken into account when building the master-flats and determining the order-edges.
* **FIXED**: master flats taken with blocking filters are no longer matched with multi-pinhole frames fro the spatial solution recipe.
* **FIXED**: master flats with identical slit-widths now matched against science frames within the data-organiser when building SOF files.
* **FIXED**: the order of the columns in the extracted & merged spectrum tables is now WAVE, FLUX (was FLUX, WAVE).
* **FIXED**: specutils dependency added to conda-forge requirements.

## v0.9.4 - December 5, 2023

* **REFACTOR:** orders are now clipped so that only the pixels deemed to be within the flux receiving regions of the order are extracted (according to the static calibration spectral format table).
* **REFACTOR:** `soxspipe prep` will now warn user if no FITS files are found in the workspace directory and quit before moving any files (safer).
* **REFACTOR:** `soxspipe session` commands will look for a sessions directory before creating any new files and folders (cleaner).
* **REFACTOR:** `read_spectral_format` function can now return limits to the usable region of each spectral order if a dispersion map is given.
* **FIXED**: fixes to make detect_continuum more robust.

## v0.9.2 - November 29, 2023

* **ENHANCEMENT:** intra-order background (scattered light) fits are now being written to FITS image files in the QC directories and reported at the end of a recipe run.
* **ENHANCEMENT:** added a `create_dispersion_solution_grid_lines_for_plot` function to allow adding dispersion solution grid to QC plots.  This is extremely useful for quickly diagnosing problems with the fits.
* **REFACTOR:** All product FITS files now pass fitverify without error or warnings. All issues were due to using '-' instead of underscores in FITS binary table column names.
* **REFACTOR:** bad-pixel values set to 0 in data extensions of products
* **REFACTOR:** nans have been replaced by zero in FITS image product
* **FIXED**: a mismatch between daofind results and the original input pixel table was causing dispersion solution to break (a recent bug introduced during code optimisations)
* **FIXED**: the internal soxspipe logger was being interfered with by astropy so that logs were sometimes getting redirected to the wrong place

## v0.9.0 - October 11, 2023

* **FEATURE:** added a `predict_product_path` function to determine the product path from a recipe's input sof file
* **FEATURE:** Merging of individual order extracted spectra from object frame into a single spectrum for each arm
* **FEATURE:** Object spectra are now extracted from the sky-subtracted frames using the Horne 86 method
* **FEATURE:** Real SOXS data is now included in the unit-test suite (starting to replace simulated data unit-tests). `soxs-disp-solu` recipe so far.
* **FEATURE:** SOXS NIR Xe line-lists added to static-calibration suite (single and multi pinhole).
* **FEATURE:** when running a recipe, `soxspipe` writes informative logs to stdoutAND to a log file adjacent to the recipe's product file(s). Error logs are also written if a recipe fails (see docs).
* **ENHANCEMENT:** recipe timing added to the end of the logs
* **ENHANCEMENT:** fitted lines from the dispersion solution are written out to file as a QC product
* **ENHANCEMENT:** flux (and other daostarfinder metrics) are now recorded in the detected line-list QC file. This will help measure degradation of arc-lamps over time.
* **ENHANCEMENT:** FWHM and pixel-scale added to fitted lines from the dispersion solution
* **ENHANCEMENT:** legends added to many of the QC plots
* **ENHANCEMENT:** OB ids now getting add to the data-organiser database tables.
* **ENHANCEMENT:** object trace FITS binary table added to stare-mode products (alongside complimentary QC plot)
* **ENHANCEMENT:** products and QC outputs are differentiated in the table reported upon recipe completion (see label column).
* **ENHANCEMENT:** verifying the master flat used to calibrate object/std spectra has the same slit-witdh as used to take the science frames
* **REFACTOR:** `init` command has been subsumed into the prep command. The `prep` command will generate a settings file to live within the prepared workspace.
* **REFACTOR:** `misc/` directory created by data-organiser even if empty
* **REFACTOR:** close matplotlib plot after writing plots to file
* **REFACTOR:** command-line startup speeds improved
* **REFACTOR:** continuum fitting code made more robust against edge cases (orders of the fit are automatically reduced if fit does not converge)
* **REFACTOR:** soxspipe now has a 'full' and a 'lite' test-suite. Using the lite suite will speed up deploying of new releases.
* **DOCS:** updated docs with a more robust SOXSPIPE upgrade path (users having issue with `conda update ...`)
* **FIXED**: sky-subtraction code and data-organiser fixed to work with binned data


## v0.8.0 - May 18, 2023

* **FEATURE:** we now have a data-organiser to sort data, prepare the required SOF files and generate reduction scripts.
* **ENHANCEMENT:** '.db', '.yaml', '.sh' and '.log' extensions skipped when moving items to the misc folder
* **ENHANCEMENT:** move information printed to STDOUT when preparing a workspace to inform the user of how the data is organised
* **ENHANCEMENT:** code can automatically adjust polynomial fitting parameters to find a dispersion solution if those provided in the settings file fail.
* **ENHANCEMENT:** uncompression of fits.Z files (if any) occurs before data-organising
* **REFACTOR:** speed & robustness improvements to dispersion solution to 2D image map conversion.
* **REFACTOR:** much fast check for product existence so recipes are quickly skipped if they have already run.
* **REFACTOR:** removed the `intermediate-data-root` setting renamed to a more accurate `workspace-root-dir`
* **REFACTOR:** removed the `reduced-data-root` setting.
* **REFACTOR:** updating all depreciated pandas commands so pipeline is now compatible with 1.X and 2.X versions of pandas 
* **FIXED** pandas 1.X and pandas 2.X were doing different things when renaming columns in data-frames. Both 1.X and 2.X now work in the pipeline.

## v0.7.2 - March 3, 2023

* **REFACTOR:** Big improvements on sky-subtraction  
* **REFACTOR:** UV order-edge detection more robust  
* **REFACTOR:** changed quickstart guide compress to gzipped tar  
* **REFACTOR:** updated default settings to be more robust  

## v0.7.1 - November 4, 2022

* **FEATURE:** UV D-Lamp and QTH-Lamp master flats now being stitched together  
* **FEATURE:** errors in error maps now being treated correctly and propagating to combined images  
* **FEATURE:** Pipeline can now 'remember' where it left off in the reduction cascade. If it has run a recipe before it will exit with a message to the user informing them how to force the recipe to rerun.  
* **FEATURE:** added a `twoD_disp_map_image_to_dataframe` function to toolkit  
* **ENHANCEMENT:** PRO CATG now written to product FITS header  
* **ENHANCEMENT:** Handling of binned images when generating flats and order-locations  
* **ENHANCEMENT:** Where possible, product files are given the same name as the SOF file used to generate them (replacing `.sof` extension with `.fits`)  
* **ENHANCEMENT:** SOF files can now contain a file 'tag' to allow users to read the SOF file contents and know exactly which files are being passed to the recipe (e.g. `MASTER_BIAS_UVB`, `LAMP,DORDERDEF_UVB` ... )  
* **ENHANCEMENT:** dispersion solution now working with simulated NIR SOXS data  
* **ENHANCEMENT:** quicklook now renders dispersion solution grid  
* **ENHANCEMENT:** \~40% speed gain in combining images.  
* **REFACTOR:** 2D Map generation now ~6-8 times faster (seeding solutions with nearest neighbour with cubic spline method)  
* **REFACTOR:** SOF filenames reworked to contain the UTC observation date instead of MJD (more in-line with ESO ecosystems)  
* **REFACTOR:** updated workflow for master bias combination  
* **REFACTOR:** updated workflow for master dark combination  
* **REFACTOR:** QC PDF plots now added to their own directory separate from the products    
* **REFACTOR:** products now sub-divided into recipe directories (e.g. `./products/soxs-mbias/`)    
* **DOCS:** mflat docs brought up-to-date    
* **FIXED:** mflat recipe now exits if flat frames are not of a consistent exptime.    


## v0.6.2 - April 13, 2022

* **ENHANCEMENT:** quickstart guide added for calibration recipes  
* **FEATURE:** QCs added for dispersion solution and order centre recipes  
* **REFACTOR:** clean up of stdout information  

## v0.6.1 - April 11, 2022

* **FEATURE:** shipping static calibration files with the code (one less thing for end-users to install and set-up)

## v0.6.0 - April 10, 2022

This is only a summary of some of the updates included in this release:

* **ENHANCEMENT:** All CSV files moved to FITS binary tables - metadata very useful for developing data organiser
* **FEATURE:** 2D image map now created by create_dispersion_solution
`subtract_calibrations` util renamed to `detrend` and added ability to flat correct
* **FEATURE:** 2D image map of wavelength values, slit-position values and order values written alongside polynomial solutions of full dispersion solution
* **FEATURE:** soxspipe now on conda
* **FEATURE:** QCs now being written to FITS header
* **FEATURE:** adding QC and product collection in mbias recipe
* **ENHANCEMENT** RON and bias structure QCs now reported by mbias
* **ENHANCEMENT** nan ignored when scaling quicklook images
* **ENHANCEMENT** RON and bias structure QCs now reported by mdark
* **ENHANCEMENT:** QCs have an option to *NOT* (`to_header`) write to FITS header (default is to write)
* **REFACTOR:** better treatment of masked pixels when stacking images (e.g. in mbias and mdark)
* **REFACTOR:** removed raw frame reports and neater QC table
* **REFACTOR:** fits header keywords neatly sorted before writing to file
* **FIX:** Correct management of mask when determining RON on bias and darks

## v0.5.1 - September 29, 2021

* **FEATURE:** recipes now have a `qc` and `products` attribute. These are pandas data frames used to collect QCs and generated products throughout the life-time of the recipe. They are printed to STDOUT at the end of the recipe (can be used in the future to send post request to health monitor API with JSON content in request body).
* **ENHANCEMENT** added code-base to conda-forge
* **ENHANCEMENT** added bottleneck to the install requirement (makes image combination more efficient)
* **ENHANCEMENT** masked pixel now coloured red in quicklook plots (easier to differentiate from good pixels)
* **ENHANCEMENT** low-sensitivity pixels in lamp-flats now identified and added to bad-pixel mask
* **ENHANCEMENT** add a verbosity flag to the command-line and a verbose parameter to each recipe
* **REFACTOR** inter-order pixel value in flats now set to unity (instead of running background fitting and subtraction)
* **REFACTOR:** recipes now have their recipe name as a `recipeName` attribute

## v0.5.0 - June 10, 2021

* **FEATURE** Added a new `filenamer` module that implements a strict intermediate and reduced file-naming scheme
* **FEATURE:** `soxs_mflat` recipe now included
* **FEATURE:** `soxs_spatial_solution` recipe is now included
* **FEATURE:** `subtract_background` utility added
* **FEATURE:** added a `detect_order_edges` object
* **FEATURE:** Added a `dispersion_map_to_pixel_arrays` function to convert from order-based and wavelength arrays to pixel arrays (first guess dispersion map only so far)
* **FEATURE:** added a quicklook function in toolkit to quickly visualise a frame
* **FEATURE:** added a toolkit module for small functions used throughout soxspipe 
* **FEATURE:** added function in toolkit to unpack an order table into lists of coordinates, one list per order.
* **FEATURE:** added image slice tool to toolkit
* **ENHANCEMENT** Added a `-o <outputDirectory>` switch to the command-line to optionally override the 'intermediate-data-root' setting in the settings file.
* **ENHANCEMENT:** added a fraction of a second tolerance when matching exptimes between darks and science/calibration frames 
* **ENHANCEMENT:** y limits now added to the order table to show limits of order locations on detector
* **REFACTOR:** Change the "SOXSPIPE PRE" date stamp keyword to "SXSPRE" to future-proof for phase III (8 character keyword limit)
* **REFACTOR:** Pandas tables are now used through-out code to pass line-lists between methods
* **REFACTOR:** refactoring of polynomial fitting has made creation of dispersion maps ~50 times faster
* **REFACTOR:** removed OBID from file names and added readout mode. This information is more helpful at the glance.
* **FIXED:** correct binning reported in product file names
* **FIXED:** lines in a sof file beginning with a `#` are considered as comments and therefore ignored by the pipeline.

## v0.4.1 - September 15, 2020

* **FEATURE:** add command-line util for soxs order_centres recipe
* **FEATURE** added the `detect_continuum` utility to fit order centre locations in single pinhole flat frames.
* **ENHANCEMENT:** added a supplementary file list for non-fits input files in set-of-file util
* **ENHANCEMENT:** adding more information residual plots & visualisation of fitting for disp solution
* **ENHANCEMENT:** check that files in the sof files exist before proceeding.
* **ENHANCEMENT:** added spectral format table lookup to detector settings file
* **REFACTOR:** moved chebyshev order/wavelength polynomials into its own class - decoupled from create_dispersion_map class

## v0.4.0 - September 3, 2020

* **FEATURE:** added create_dispersion_map class to be used in `soxs_disp_solution` and `soxs_spatial_solution`
* **FEATURE:** added a `detrend` method to subtract calibration frames (bias and dark) from an input frame
* **FEATURE:** added the dispersion solution recipe and unit tests
* **FEATURE:** added the disp_solution command-line tool
* **DOCS:** major docs overhaul
* **ENHANCEMENT:** added predicted lines lists to detector parameter file
* **ENHANCEMENT:** DPR CATG and DPR TECH added to metadata of sof imagefilecollection objects
* **ENHANCEMENT:** wcs copied from a single frame into the combined frames during clip and stack
* **REFACTOR:** bad-pixel map paths abstracted to detector settings files
* **REFACTOR:** renaming of unit-testing test data directories
* **REFACTOR:** only filenames reported by sof summaries when files are found in the same directory (easier to read on terminal) 
* **FIXED:** fixed detector science pixels for UVB

## v0.3.1 - August 25, 2020

* **FEATURE:** recipe & command-line tool for master dark creation (`mdark`)
* **ENHANCEMENT:** default binning add to detector settings file
* **ENHANCEMENT:** added mixed exposure time unit test for dark-frames
* **ENHANCEMENT:** added default values for gain and ron in the detector settings files. Default values can be overwritten if correct GAIN and RON are found in fits-headers (overwritten for UVB and VIS but not NIR for XShooter)
* **ENHANCEMENT:** can now interrupt "~" as home directory in sof file path
* **FIXED:** binning factor used when trimming frames
* **FIXED:** the trimming dimensions of NIR frames - bad-pixel map now aligns correctly with data frame
* **FIXED:** science pixels for all 3 xshooter detectors in parameters file

## v0.3.0 - August 18, 2020

* **FEATURE:** added a `write` method to the `_base_recipe` to write frames to disk (renames extensions to ESO preferred naming scheme)
* **FEATURE:** detector lookup class added alongside yaml files to host detector specific parameters (rotation, science-pixels etc). Code has been updated to remove hard-wired detector values.
* **FEATURE:** added a cleanup method to remove intermediate file once recipe completes
* **ENHANCEMENT:** parameters for clip and stack method added to the settings files
* **ENHANCEMENT:** added strict typing of data and variables with astropy units to avoid silent mistakes in frame arithmetic 
* **ENHANCEMENT:** added mixing of readout speeds to input frame verification checks
* **ENHANCEMENT:** added readnoise and gain to the list of keyword values to check during frame verification
* **ENHANCEMENT:** inject a 'SOXSPIPE PRE' keyword with timestamp value into prepared frames
* **ENHANCEMENT:** check frames for 'SOXSPIPE PRE' keyword before preparing - raises exception if found
* **ENHANCEMENT:** ron and gain are added to the recipe's detector lookup dictionary during frame verification (so they don't need read again later)
* **REFACTOR:** moved stacking code to it own `clip_and_stack` method hosted in the `_base_recipe`
* **REFACTOR:** moved basic input frame verifications to the `_base_recipe` - so not to repeat code
* **REFACTOR:** removed python 2.7 support - not feasible with CCDProc
* **DOCS:** added workflow diagrams to the documentation for many of the methods implemented (`prepare_frames()`, ``clip_and_stack``)

## v0.2.0 - February 27, 2020

* **FEATURE** added keyword lookups - abstracting exact keyword names from code
