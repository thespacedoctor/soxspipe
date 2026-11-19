# DY-2297: sources for SPEC_RES, LAMRMS and SPEC_ERR

These notes answer the research question in DY-2297 (map DY-2291): what `soxspipe` already measures that can fill `SPEC_RES`, `LAMRMS`, `LAMNLIN`, `SPEC_ERR` and `SPEC_SYE` in a Phase 3 `SCIENCE.SPECTRUM` header, and how well those values describe the science slit. They are for whoever designs the Phase 3 keyword code in `soxspipe/commonutils/phase3.py`.

## Answer

The only wavelength calibration a science SOF carries is the multi-pinhole `DISP_TAB` (and its `DISP_IMAGE`) from `soxs-spat-solution`. Its header (since v0.16.0) and its `_FITTED_LINES.fits` QC table hold everything needed for `LAMRMS`, `LAMNLIN` and `SPEC_ERR`, but no existing QC is `LAMRMS` as defined: `X RES MEDIAN` and `X RES SD` are the median and standard deviation of the absolute residual in pixels, not an RMS in nm, and the dispersion axis is x for VIS but y for NIR. `LAMRMS` must be computed from the per-line table as the RMS of the along-dispersion residual times the per-line `pixelScaleNm`, over unclipped rows; this gives about 0.0035 nm (VIS) and 0.0066 nm (NIR) on current products. `R PIN MEDIAN` is not usable for `SPEC_RES`: the pinhole resolving power (about 9400 VIS and 12400 NIR) is close to that of the 0.5 arcsec slit and is 1.5 to 2 times higher than the 1.0 arcsec slit value. Arc lines through the science slits and night-sky lines in real science frames give R of about 10600 to 10900 (NIR) and 8900 (VIS) at 0.5 arcsec, 6300 to 6500 (NIR) and 4900 to 5400 (VIS) at 1.0 arcsec, about 4200 to 4300 (NIR) and 3300 to 3600 (VIS) at 1.5 arcsec, and roughly 1000 to 1400 at 5.0 arcsec. SOXS takes slit arcs at all four slit widths in its routine calibration plan, but `soxspipe` never uses them in a science cascade. ESO's own X-Shooter products write a nominal per-slit-width `SPEC_RES` derived from slit arc calibrations, and compute `SPEC_ERR` as `LAMRMS / LAMNLIN`, not `LAMRMS / sqrt(LAMNLIN)`. `SPEC_ERR` from the arc fit is tiny (about 1e-4 nm) compared with the skyline shifts the pipeline applies to science frames (0.06 to 0.5 nm), so the honest wavelength error budget sits in `SPEC_SYE`, which nothing in `soxspipe` measures yet.

**Confidence: high** for the keyword sources, the pixel-to-nm conversion, the axis trap and the X-Shooter practice (read from code, product headers and ESO archive headers). **Medium-high** for the slit resolving powers at 0.5, 1.0 and 1.5 arcsec: two independent methods agree within about 5% (NIR) and 10% (VIS), but they rest on one night of arcs and a handful of science frames, and the fit window moves wide-slit values by up to 20%. **Low** for 5.0 arcsec (few usable lines). **Medium** for what `SPEC_SYE` should hold: that is a design choice for Dave.

## Contents

- [What the standard asks](#what-the-standard-asks)
- [Which products of a science SOF carry the values](#which-products-of-a-science-sof-carry-the-values)
- [Converting a pixel residual to nm](#converting-a-pixel-residual-to-nm)
- [LAMRMS, LAMNLIN and SPEC_ERR on current products](#lamrms-lamnlin-and-spec_err-on-current-products)
- [Pinhole resolving power against the science slits](#pinhole-resolving-power-against-the-science-slits)
- [What X-Shooter Phase 3 products do](#what-x-shooter-phase-3-products-do)
- [SPEC_SYE](#spec_sye)
- [Consequences for other tickets on the map](#consequences-for-other-tickets-on-the-map)
- [Method and reproducibility](#method-and-reproducibility)
- [Sources](#sources)

## What the standard asks

Table 8 of the ESO Science Data Products Standard [S1] sets, for `SCIENCE.SPECTRUM`, all in the primary header:

| Keyword | Status | Definition (section 5.6 or 5.11 of [S1]) |
|---|---|---|
| `SPEC_RES` | Mandatory | Average resolving power lambda / delta lambda; "could be estimated by measuring the FWHM of the ARC-lamp emission lines"; defined at the central wavelength (`WAVELMIN` + `WAVELMAX`) / 2 if it varies. |
| `SPEC_ERR` | Recommended | Statistical error in the spectral coordinate [nm]. "In most cases, including UVES and XSHOOTER, the formula SPEC_ERR = LAMRMS / sqrt(LAMNLIN) is a good approximation", and that value "is a lower limit to SPEC_ERR". |
| `SPEC_SYE` | Recommended | Systematic error in the spectral coordinate [nm], including "any residual offset of the wavelength calibration and any error related to broadening of the lines". May be zero if insignificant or already corrected. |
| `LAMRMS` | Optional | RMS of the residuals of the wavelength solution [nm], sqrt(sum R_i^2 / N), with N = `LAMNLIN`. |
| `LAMNLIN` | Optional | Number of arc lines used in the fit of the wavelength solution. |

`soxspipe` writes none of these today. `soxspipe/resources/xsh_keywords.yaml` (lines 251 to 255) has aliases for them; `soxs_keywords.yaml` and `phase3.py` do not.

## Which products of a science SOF carry the values

A science SOF lists `DISP_TAB_<ARM>` and `DISP_IMAGE_<ARM>` from `soxs-spat-solution` (the multi-pinhole 2D solution), for example `sof/20260128T072009_NIR_3_STARE-OBJ_SLIT1_0_600_0S_SOXS_PAQS135626+042348.sof` in `/Volumes/SOXS/reductions/2026_jan`. The single-pinhole `soxs-disp-solution` product is not in the science SOF; it appears only as `ESO PRO REC1 CAL3` in the spatial-solution header, and its line set is a subset (one pinhole, `slit_index` 4).

`create_dispersion_map.calculate_residuals` (`soxspipe/commonutils/create_dispersion_map.py`, lines 2341 to 2700) writes these header QCs to the `DISP_TAB` primary header, overall and per order:

| QC | What it is | Unit |
|---|---|---|
| `DETLINES TOT` | lines in the input line list (all pinholes) | count |
| `DETLINES NUM` | lines detected (all pinholes), before clipping | count |
| `DETLINES CLIP NUM` | detected lines clipped during the fit | count |
| `X RES MIN/MEDIAN/MAX/SD`, `Y RES ...` | min, median, max and std of the absolute fit residual along x and y | px |
| `XY RES ...` | the same for the combined residual | px |
| `FWHM PIN MEDIAN/SD` | median and std of the pinhole-spot FWHM (photutils `IRAFStarFinder`, a round-profile FWHM, not FWHM along dispersion) | px |
| `R PIN MEDIAN/SD` | median and std of lambda / (`pixelScaleNm` x `fwhm_pin_px`) | none |

The QCs `FWHM PIN` and `R PIN` came in with commit `45837574` (2026-02-23, first tag v0.16.0). Products made with v0.15.2 or earlier (for example all of `/Volumes/SOXS/reductions/2026_jan`) only carry `TLINE`, `NLINE`, `PLINE` and `CLINE`, so Phase 3 export from such products cannot read these values from the header.

The per-line table `<stem>_FITTED_LINES.fits` (written to `qc/`, not listed in the science SOF) has one row per pinhole detection with `wavelength`, `order`, `slit_index`, `residuals_x`, `residuals_y`, `sigma_clipped`, `fwhm_pin_px`, `R_pin` and `pixelScaleNm` (column list `FITTED_LINE_COLUMNS`, lines 178 to 210). It also has `R_slit` and `fwhm_slit_px`, but these are zero unless a slit arc was in the spatial-solution SOF (see below).

Check on 2026-09-22 VIS products (`/Volumes/SOXS/reductions/2026_sept_oct_commissioning_dev`): `DETLINES NUM` = 1206 and `DETLINES CLIP NUM` = 128 in the header; the `FITTED_LINES` table has 1206 rows of which 1078 = 1206 - 128 have `sigma_clipped` false. So the unclipped count is `DETLINES NUM - DETLINES CLIP NUM`.

## Converting a pixel residual to nm

- The dispersion axis is set per arm by `dispersion-axis` in `soxspipe/resources/soxs_detector_parameters.yaml`: `x` for VIS (line 36) and `y` for NIR (line 108). So the along-dispersion residual is `X RES` for VIS but `Y RES` for NIR. Reading `X RES` for both arms would give the cross-dispersion residual for NIR.
- `pixelScaleNm` is computed per line (lines 2436 to 2458) by evaluating the fitted x and y polynomials at lambda + 2 nm and lambda - 2 nm and dividing 4 nm by the distance between the two points. It is the local nm per pixel along the order trace. On current products it is 0.019 to 0.052 nm/pix (VIS, median 0.034 nm/pix) and 0.034 to 0.093 nm/pix (NIR, median 0.048 nm/pix).
- A wavelength residual per line is then delta lambda_i = r_i x `pixelScaleNm`_i, where r_i is the residual along the dispersion axis. Because the order trace is not exactly parallel to the detector axis (the UV-VIS spatial axis is rotated relative to the spectral axis by design [S5]), a cleaner value projects the (x, y) residual onto the local trace direction, or evaluates the inverse map at the observed position. The difference is small (the trace tilt is a few degrees) but the projection is what the code should do.
- The header median and SD of the absolute residual cannot be turned into an RMS exactly. With the per-line table available, compute `LAMRMS` directly.

## LAMRMS, LAMNLIN and SPEC_ERR on current products

Computed from the `FITTED_LINES` tables, unclipped rows only, along-dispersion residual times `pixelScaleNm`:

| Product | Arm | Unclipped rows | Distinct lines (order, wavelength) | RMS residual | `LAMRMS` | `LAMRMS / sqrt(rows)` | `LAMRMS / sqrt(distinct)` |
|---|---|---|---|---|---|---|---|
| `20260122T101652_VIS_1X1_1_SSOL` (v0.18.0) | VIS | 1051 | 117 | 0.100 pix | 0.0036 nm | 1.1e-4 nm | 3.3e-4 nm |
| `20260923T093959_VIS_1X1_1_SSOL` (dev) | VIS | 1078 | 120 | 0.095 pix | 0.0033 nm | 1.0e-4 nm | 3.1e-4 nm |
| `20260122T102149_NIR_3_SSOL` (v0.18.0) | NIR | 2361 | 267 | 0.131 pix | 0.0066 nm | 1.4e-4 nm | 4.0e-4 nm |
| `20260128T103034_NIR_3_SSOL` (v0.15.2) | NIR | 2415 | 272 | 0.144 pix | 0.0069 nm | 1.4e-4 nm | 4.2e-4 nm |

The first two rows are in `P32.15_soxspipe_oct_26_release/sessions/base/qc/2026-01-21/soxs-spat-solution/` and `2026_sept_oct_commissioning_dev/sessions/base/qc/2026-09-22/soxs-spat-solution/`. The single-pinhole solution of 2026-09-22 gives the same `LAMRMS` (0.0031 nm from 125 lines), so the multi-pinhole fit does not degrade the along-dispersion residual.

`LAMRMS` / lambda is about 4e-6 in both arms, or 1.2 km/s.

What `LAMNLIN` should count is a choice. The standard says "number of arc lines used in the fit". Each arc line appears once per pinhole (9 pinholes), so rows outnumber distinct lines by about 9. The nine detections of one line are not independent wavelength measurements, so dividing by sqrt(rows) understates the error by a factor of 3. ESO's X-Shooter products count rows (see below), so counting rows matches the house precedent; counting distinct lines matches the meaning of the formula.

## Pinhole resolving power against the science slits

### Why they differ

The pinhole is a 0.5 arcsec aperture: the NIR slit stage holds "0.5'', 1.0'', 1.5'', 5'' [slits], single pinhole 0.5'' and multi pinhole" [S6]. A science slit wider than the pinhole projects a wider line on the detector, so `R PIN` measures the resolving power of a 0.5 arcsec aperture, not of the slit the science frame used. The published design values are all for the 1 arcsec slit [S2, S3, S4, S5].

### Measurements

Two independent measurements, both made without running a reduction (method below):

1. Arc lines (ThAr for VIS, ArNeHgXe for NIR) in raw slit arcs from 2026-01-22 around 10:12 to 10:21 UT, one per slit width, mapped through the `DISP_IMAGE` of the same morning (v0.18.0), within 3 arcsec of the slit centre.
2. Night-sky emission lines in `SKYSUB + SKYMODEL` of real stare science frames (v0.18.0 and develop), 2 to 5 arcsec from the slit centre so that the target does not contribute.

Median R over isolated lines (interquartile range in brackets):

| Aperture | NIR arc | NIR sky | VIS arc | VIS sky | Pinhole `R PIN MEDIAN` |
|---|---|---|---|---|---|
| Pinhole (0.5 arcsec) | | | | | NIR 12401, VIS 9431 (2026-01-21) |
| 0.5 arcsec slit | 10563 (8988 to 11593) | 10928 (9384 to 11646), 1 frame | 8965 (8010 to 10094) | 8900 to 11600, 2 frames, 3 to 5 lines each | |
| 1.0 arcsec slit | 6312 (5888 to 6581) | 6537 and 6516, 2 frames | 4906 (4595 to 5380) | 5373 and 5480, 2 frames | |
| 1.5 arcsec slit | 4292 (4120 to 4392) | 4211 and 4204, 2 frames | 3271 (3160 to 3546) | 3595 and 3513, 2 frames | |
| 5.0 arcsec slit | about 1350, 9 lines, unreliable | not measured | about 1000, 6 lines, unreliable | not measured | |
| Design, 1.0 arcsec [S3, S4, S5] | 5000 to 5300 | | 4600 to 5000 per quasi-order (4700 average) | | |

Frames used for the sky lines: NIR 0.5 arcsec `20260128T024829` (LTT3218, 150 s), NIR 1.0 arcsec `20260122T063210` and `20260121T045115`, NIR 1.5 arcsec `20260118T024705` and `20260109T021758`, VIS 1.0 arcsec `20260122T063212` and `20251010T054324`, VIS 1.5 arcsec `20250912T040511` and `20251012T015521` (all in `P32.15_soxspipe_oct_26_release`), VIS 0.5 arcsec `20260923T083849` and `20260925T061929` (in `2026_sept_oct_commissioning_dev`). VIS sky values come mostly from orders 1 (r) and 4 (i); the u and g orders have too few sky lines.

What this shows:

- R scales close to 1 / slit width for 1.0 and 1.5 arcsec, as a slit-limited spectrograph should. R at 0.5 arcsec is a little below the pinhole value, consistent with a 0.5 arcsec aperture plus the optics.
- At 1.0 arcsec, the pinhole R is about 1.9 times (NIR) and 1.7 to 1.9 times (VIS) the slit R. Writing `R PIN MEDIAN` as `SPEC_RES` would overstate the resolving power of most science spectra by that factor.
- Arcs and sky agree within about 5% in NIR. In VIS the arc values are 8 to 10% lower than the sky values; the sky values come only from the r and i orders, and the dense ThAr spectrum has more blends. Either way, arc frames through the science slit are a sound source, and the systematic difference between methods is about 10%.
- The measured 1.0 arcsec values are above the design values (NIR 6300 to 6500 against 5000 to 5300). This agrees with the first-light statement that "SOXS meets or exceeds its design specifications in spectral resolution" [S7]. It is not proof the measurement has no bias: the fit uses a single Gaussian, and some OH and ThAr lines are blends, which would push R down, not up.
- R varies along and between orders (NIR 1.0 arcsec arc: 5900 to 6800 by order). `SPEC_RES` at the central wavelength needs either the value of the order covering (`WAVELMIN` + `WAVELMAX`) / 2 or one average per arm and slit.
- Sky lines fill the slit, so they give the slit-limited value. A point source observed in seeing better than the slit has a narrower line spread function. Both ESO release descriptions call the arc value a lower limit for that reason [S8, S9].

### Slit arcs exist but are not used

The definitive raw workspace `/Volumes/SOXS/raw_archive/soxs_commissioning_data/soxspipe.db` has slit arcs (`ECHELLE,SLIT`, `LAMP,WAVE`) at all four slit widths in both arms from 2025-05 to 2026-10: about 720 per slit width in VIS and about 177 per slit width in NIR, as many as the pinhole arcs. But:

- `soxspipe/resources/soxs_sof_map.yaml` has no entry for slit arcs, so the data organiser never puts them in an SOF (the `P32.15` workspace database has no slit arcs at all).
- `soxs_spatial_solution._read_slit_arc_frame` accepts one slit arc if a hand-written SOF carries it, and `_calculate_resolution_on_slit` then fills `R_slit` per line, for the QC plot only. No header QC is written, the slit width is that of the one arc, and the function prints `COME BACK HERE AND TRASH SOME LINES` for every line. It is a prototype, not a product.

## What X-Shooter Phase 3 products do

The XSHOOTER release description v2.9 (2026-07-21) [S8] says: "The header key SPEC_RES contains the nominal resolving power, as derived from the arc lamp calibrations of the same slit width. The actual spectral resolving power of the science data has not been measured but could actually be higher, depending on seeing conditions and exposure time." The UVES release description [S9] says the same: a nominal value based on the slit width, derived from QC monitoring of arc lamp exposures, and a lower limit.

Headers of three current X-Shooter IDPs from one 2024-07-09 observation (pipeline `xshoo/3.6.5a`, read from the ESO archive header service [S10]):

| Product | Arm, slit | `SPEC_RES` | `LAMRMS` | `LAMNLIN` | `SPEC_ERR` | `SPEC_SYE` |
|---|---|---|---|---|---|---|
| `ADP.2024-07-16T13:56:15.104` | UVB, 1.0 arcsec | 5453 | 0.002963 nm | 3363 | 8.81e-7 nm | 5.17e-4 nm |
| `ADP.2024-07-16T13:56:15.108` | VIS, 0.9 arcsec | 8935 | 0.003071 nm | 4157 | 7.39e-7 nm | 3.10e-4 nm |
| `ADP.2024-07-16T13:56:15.113` | NIR, 0.9 arcsec | 5573 | 0.010336 nm | 1859 | 5.56e-6 nm | 9.28e-5 nm |

- `SPEC_RES` is one fixed value per arm and slit width (same values on every product of that setup in the archive query).
- In all three, `SPEC_ERR` = `LAMRMS / LAMNLIN` exactly (0.002963 / 3363 = 8.81e-7), not `LAMRMS / sqrt(LAMNLIN)` as the standard recommends. That looks like an error in ESO's own processing; `soxspipe` should follow the standard, not this precedent.
- `LAMNLIN` is in the thousands, so it counts multi-pinhole detections, not distinct lines.
- The source of `SPEC_SYE` is not documented in [S8].

## SPEC_SYE

The science reduction already measures a systematic offset of the arc solution: skyline tuning writes `SKY SHIFT O<n>` (nm) and `SKY SHIFT O<n> PX` to `EXTRACTED_MERGED` (`horne_extraction.py`, around line 1991). On `20260122T063210_NIR` the applied shifts are 0.19 to 0.49 nm (4 to 7 pix); on `20260122T063212_VIS` they are 0.06 to 0.12 nm (2 to 3 pix), and the u and g shifts are predicted from r and i (`UG SHIFT SOURCE model`), not measured. These shifts are 100 to 1000 times `SPEC_ERR`. Once applied, the remaining systematic error is the uncertainty of the shift, which `soxspipe` does not report. Choices for `SPEC_SYE`:

- the scatter of skyline centroids about the tuned solution, per order, combined over the spectrum (needs new code);
- for orders whose shift is modelled rather than measured (VIS u and g), a value from the validation of that model (the r/i to u/g transfer is not yet validated, see the 2026-10-08 skyline tuning study);
- zero only if Dave rules that the tuning residual is insignificant, which the numbers above do not support yet.

## Consequences for other tickets on the map

- The Phase 3 keyword ticket must read `LAMRMS` from `FITTED_LINES`, not from a header QC, or a new header QC (for example `LAMRMS`, in nm, over unclipped rows along the dispersion direction) must be added to `DISP_TAB`. Since Phase 3 export reads the science SOF, the second is simpler. The QC must use the per-arm dispersion axis.
- `SPEC_RES` needs a new source. Two options: (a) a fixed table per arm and slit width, measured once from slit arcs and maintained like X-Shooter's nominal values; (b) put slit arcs in the cascade (new `soxs_sof_map.yaml` rule, measure per-order R on the slit arc matching the science slit, write a QC). Option (a) matches ESO practice and needs no cascade change; (b) tracks instrument changes. This is a decision for Dave.
- `R PIN MEDIAN` stays useful as an instrument health QC, not as a science keyword.
- `SPEC_SYE` should come from skyline tuning; that couples the Phase 3 work to the skyline tuning code and its u/g open question.
- Exports from products older than v0.16.0 cannot get these values from headers.

## Method and reproducibility

All scripts are in the session scratchpad and are not committed; the logic is short enough to restate.

- `LAMRMS`: read `_FITTED_LINES.fits`; keep `sigma_clipped == False`; r = `residuals_x` (VIS) or `residuals_y` (NIR); `LAMRMS = sqrt(mean((r x pixelScaleNm)^2))`.
- Line widths: for each order in the `DISP_IMAGE` (`WAVELENGTH`, `SLIT`, `ORDER` extensions), select pixels in the slit range, sort by wavelength, build a median profile at half-pixel steps, find peaks with prominence 20 (sky) or 30 (arc) times the robust noise, keep peaks with no neighbour within 8 pix, fit a Gaussian plus constant to the individual pixels within 5 pix, and keep fits with 0.8 pix < FWHM < 8 pix. R = centre / FWHM. Raw arcs are trimmed to the `science-pixels` region of `soxs_detector_parameters.yaml`, pixels above 40000 ADU are dropped, and no bias, dark or lamp-off frame is subtracted (the fit has a constant term). Because each pixel already integrates over its width, these FWHM values are comparable to the pinhole FWHM from `IRAFStarFinder`, which also works on pixel data.
- The fit window and isolation distance scale with an expected-width factor `WIN` (1 for 0.5 arcsec, 1.5 for 1.0 arcsec, 2 for 1.5 arcsec, 5 for 5.0 arcsec on arcs; 1 for all sky frames). A first arc pass with `WIN` = 1 and no saturation check gave 6250 (NIR 1.0), 3954 (NIR 1.5), 4762 (VIS 1.0) and 2688 (VIS 1.5), so wide-slit values depend on the window by up to 20%. Lines with any saturated pixel (VIS above 60000 ADU, NIR above 40000 ADU) in the window are rejected. At 5.0 arcsec only 6 to 9 isolated lines survive and the result is rough.
- No reductions were run.

## Sources

- [S1] ESO, Science Data Products Standard, ESO-044286 issue 8, 2022-03-15: sections 5.6 and 5.11, Table 8. Local Markdown copy: `P32.11_eso_phase_3/output/p3sdpstd_md/p3sdpstd.md`.
- [S2] Schipani, P. et al. 2018, "SOXS: a wide band spectrograph to follow up transients", Proc. SPIE 10702, 107020F, doi:10.1117/12.2307349 (2018SPIE10702E..0FS): Resolution-Slit product about 4500; Table: UV-VIS 3500 to 7000 (about 4500 average), NIR 5000, both for a 1 arcsec slit; slit widths 0.5, 1, 1.5, 5 arcsec; detector scales 0.28 and 0.25 arcsec/pix.
- [S3] Zanmar Sanchez, R. et al. 2018, "Optical design of the SOXS spectrograph for ESO NTT", Proc. SPIE 10702, 107022Y (2018SPIE10702E..27Z), arXiv:1809.01521: Table 2, average R about 4700 (UV-VIS) and about 5300 (NIR) for a 1 arcsec slit.
- [S4] Vitali, F. et al. 2018, "The NIR spectrograph for the new SOXS instrument at the NTT", Proc. SPIE 10702 (2018SPIE10702E..28V), arXiv:1809.01511: R = 5000 for a 1 arcsec slit.
- [S5] Rubin, A. et al. 2018, "MITS: the multi-imaging transient spectrograph for SOXS", Proc. SPIE 10702 (2018SPIE10702E..2ZR), arXiv:1809.01522: Table 4, average R per quasi-order u 4995, g 4634, r 4633, i 4601 (1 arcsec slit); spatial axis slightly rotated relative to the spectral axis.
- [S6] Genoni, M. et al. 2024, "SOXS NIR: optomechanical integration and alignment, optical performance verification before full instrument integration", Proc. SPIE 13096 (2024SPIE13096E..2TG), arXiv:2407.17244: section 5, single pinhole 0.5 arcsec; NIR R measured on 1 arcsec slit lab arcs compliant with R >= 3500.
- [S7] Schipani, P. et al. 2026, "SOXS at the NTT: completion and first on-sky performance of the new transient spectroscopic facility", Proc. SPIE 14149, doi:10.1117/12.3101083 (2026SPIE14149E..25S), abstract.
- [S8] ESO, Phase 3 Data Release Description, XSHOOTER release 1, version 2.9, 2026-07-21 (Hanuschik, Hummel, Wolff, Marsset), section Data Quality, paragraph SPEC_RES. `https://www.eso.org/rm/api/v1/public/releaseDescriptions/70`.
- [S9] ESO, Phase 3 Data Release Description, UVES, 2018-01-15, paragraph Resolving power. `https://www.eso.org/rm/api/v1/public/releaseDescriptions/71`.
- [S10] ESO archive header service, `https://archive.eso.org/hdr?DpId=<dp_id>`, for the three `ADP.2024-07-16T13:56:15.*` products; found with the ObsCore TAP service `https://archive.eso.org/tap_obs`.
- [C1] `soxspipe` at `origin/develop` `ad73aeb0` (2026-11-19): `soxspipe/commonutils/create_dispersion_map.py`, `soxspipe/recipes/soxs_spatial_solution.py`, `soxspipe/resources/soxs_detector_parameters.yaml`, `soxspipe/resources/soxs_sof_map.yaml`, `soxspipe/resources/xsh_keywords.yaml`, `soxspipe/commonutils/horne_extraction.py`.
- [D1] Real products under `/Volumes/SOXS/reductions/` (`2026_jan`, `P32.15_soxspipe_oct_26_release`, `2026_sept_oct_commissioning_dev`) and raw frames under `/Volumes/SOXS/raw_archive/soxs_commissioning_data/raw/2026-01-22/`.
