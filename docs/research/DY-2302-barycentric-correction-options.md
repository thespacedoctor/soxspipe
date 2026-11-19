# DY-2302: How to apply a barycentric correction in soxspipe

These are research notes for the Linear ticket DY-2302 on the map DY-2291 ("soxspipe science spectra that meet the ESO Phase 3 standard"). They collect the facts and options for a barycentric wavelength correction in `soxspipe`. They make no decision. The grilling ticket DY-2303 ("Where the barycentric correction sits in the reduction") decides. The notes were checked against `develop` at `ad73aeb0` on 2026-11-19.

## Summary

The correction is small but measurable for SOXS. Over a year it is at most ±30.4 km/s at La Silla. That is about 0.45 of a resolution element (about 67 km/s at R = 4500), 1.6 native UV-VIS pixels (about 18.5 km/s each), or 2.3 native NIR pixels (about 13.3 km/s each). How you apply it matters much less than whether and where you apply it. `astropy` gives the velocity to the mm/s level from inputs that every SOXS raw frame already has: `ESO TEL GEOLAT`/`GEOLON`/`GEOELEV`, `MJD-OBS` with the exposure time, and the target coordinates. The error budget is dominated by the wavelength calibration, not by the correction. A whole observing sequence drifts by at most about 0.13 km/s per hour, so one correction at the mean epoch of a stack is good enough for every SOXS sequence. The flux calibration is not affected in practice: shifting the wavelengths by 30 km/s changes the response factor by less than 0.06 % and the extinction factor by less than 0.04 %. Two facts may change the "settled" note on the map. First, the Phase 3 standard does not require a barycentric scale: `SPECSYS` is mandatory, but `TOPOCENT` is one of its listed values. Second, ESO's own X-Shooter Phase 3 spectra are topocentric, with the correction stored in `ESO QC VRAD BARYCOR`. `soxspipe` has no barycentric code today. The SOXS keyword map has no `SPECSYS` alias and no `TEL GEO*` aliases.

## Options

Each option ends in the Phase 3 product carrying `SPECSYS` and the velocity that was (or was not) applied.

| Option | What it does | For | Against |
|---|---|---|---|
| A. Do not apply; record only | Keep the wavelengths topocentric, write `SPECSYS = 'TOPOCENT'`, and write the barycentric (and optionally heliocentric) velocity to header cards, as X-Shooter does. | No change to the wavelength chain, the stacking or the flux calibration. Matches the ESO X-Shooter Phase 3 precedent and is valid under the standard. Telluric and sky features stay where users expect them. | Every user must apply the correction. It contradicts the map note "`SPECSYS = 'BARYCENT'`" and needs Dave to reopen that ruling. |
| B. Scale the per-pixel wavelengths before the order merge | In `horne_extraction.merge_extracted_orders`, after the skyline tuning, multiply `wavelengthMean` by `1 + v/c` and then resample onto the usual fixed grid. This is what ESPRESSO does (correct the wavelength solution, then rebin and merge). | One interpolation only. The merged grid stays on round multiples of 0.02 nm or 0.06 nm, so the nod and offset per-cycle stack (which groups rows by exact `WAVE`) still works, and each cycle can get its own velocity. | The correction reaches every downstream product: `_EXTRACTED_MERGED`, the `_NOTFLAT` efficiency pass, the QC plots and, unless you exclude them, the flux standards that build the response. The extraction needs the velocity, so the recipes must compute it and pass it in. A 2D stack still needs one velocity for the whole stack. |
| C. Scale the wavelength column of the final product only | After the flux calibration, at the Phase 3 writer, multiply `WAVE` by `1 + v/c` with no resampling, and write `SPECSYS = 'BARYCENT'`. | The flux calibration stays on the observed scale, the same scale the response was measured on. One place in the code, and the internal products stay topocentric. No second interpolation. | The bins are no longer round numbers (step `0.02 nm × (1 + v/c)`). The current rounding of `WAVE` to 2 decimals must go, or it adds up to ±0.005 nm (±3 km/s at 500 nm). A stack gets one velocity at its mean epoch. |
| D. Option C, then resample back onto the round grid | As C, then interpolate back onto the original grid. | The bins stay at round values. | A second interpolation smooths the spectrum and correlates the noise between bins. It brings nothing that C does not already give. |
| E. Correct each frame, then stack the 1D spectra | Extract every frame or cycle on its own, correct each at its own midpoint, resample to a common grid and combine, as PypeIt does. | Exact per-frame velocities. | It means dropping the default 2D stack before extraction (`stack_before_extraction: True`), which helps faint targets. The gain is at most about 0.3 km/s for a 3 h sequence, which is less than 2 % of a pixel. |

Points that apply to every option that applies the correction (B to E):

- Use the barycentric correction (`kind="barycentric"`), not the heliocentric one. The two differ by at most a few m/s for SOXS targets (2.95 m/s for V818 Sco, see below), so the choice follows the convention. `BARYCENT` is the frame that ESPRESSO uses.
- Apply the factor as `λ_bary = λ_obs × (1 + v_corr/c)`. This follows the `astropy` convention that the correction is multiplicative (`rv + vcorr + rv × vcorr/c`). The relativistic form PypeIt uses, `sqrt((1 + β)/(1 − β))`, differs by about 1.5 m/s at 30 km/s (the β²/2 term).
- The air wavelength scale (DY-2292) is not a problem. Scaling air wavelengths directly instead of converting to vacuum, shifting and converting back adds less than 0.9 m/s.
- Telluric absorption and sky residuals move by the correction in the product. This is expected (telluric correction is out of scope on the map).

## Evidence

### What the Phase 3 standard asks for

- Section 5.6 lists `SPECSYS` as "(S)", "Frame of reference for spectral coordinates", with values from Table 2 (`output/p3sdpstd_md/p3sdpstd.md`, lines 762 and 782 to 797). Table 2 includes `TOPOCENT` (topocentric), `BARYCENT` (solar system barycentric) and `HELIOCEN` (heliocentric).
- In Table 8, `SPECSYS` is mandatory for `SCIENCE.SPECTRUM` (`p3sdpstd.md`, line 1992: row 41, "M" in the spectrum columns).
- The standard does not say that a spectrum must be barycentric. It requires that the frame is declared.
- `TMID` is "Exposure midpoint (MJD). It shall be set to (MJD-OBS + MJD-END)/2.0" and `TELAPSE` is `MJD-END − MJD-OBS` (`p3sdpstd.md`, section 5.6.1). For a stack with gaps, `TMID` is the midpoint of the whole span, not the mean of the frame midpoints.

### What ESO does for X-Shooter

- ESO's release description for the X-Shooter Phase 3 collection says: "The wavelength scale is topocentric; no corrections for barycentric or heliocentric motion have been applied. The corresponding values have been calculated by the pipeline and are stored in the header (HIERARCH.ESO.QC.VRAD.BARYCOR and …HELICOR, in km/s)." ([XSHOOTER_ECHELLE release description v1.4](https://www.eso.org/sci/observing/phase3/data_releases/xshooter_echelle_v1.4.pdf), page 7, "Wavelength scale").
- Archive products confirm this. Three `SCIENCE.SPECTRUM` ADP files of AT2017gfo (`/Volumes/SOXS/xshooter/objects/xsh_nodding_at2017gfo_kilanova/ADP.2017-09-22T07:51:47.{341,346,350}.fits`, `xsh_scired_slit_nod`, `xshoo/2.7.0b`) carry `SPECSYS = 'TOPOCENT'`, `ESO QC VRAD BARYCOR = -24.004` and `ESO QC VRAD HELICOR = -24.007`. The ancillary 2D spectra carry the same QC cards.
- For that product, `astropy` gives −24.066 km/s at `MJD-OBS` and −24.091 km/s at the midpoint of the six-frame span. The X-Shooter value differs by 0.06 to 0.09 km/s. That is less than 1 % of a SOXS pixel. I did not read the X-Shooter CPL source to find the algorithm or the epoch it uses. The X-Shooter pipeline manual (issue 12.17, [ftp.eso.org](https://ftp.eso.org/pub/dfs/pipelines/instruments/xshooter/xshoo-pipeline-manual-12.17.pdf)) does not mention the barycentric correction. It only says that the reference spectrum of the flux standard is moved to the radial velocity of the observed star.
- The existing `soxspipe` X-Shooter keyword map already has `SPECSYS: "SPECSYS"` and `SPECSYS_VALUE: "TOPOCENT"` (`soxspipe/resources/xsh_keywords.yaml`, lines 216 to 217), plus `GEOLAT` and `GEOLON` (lines 21 to 22). No code reads `SPECSYS` today. The SOXS map `soxspipe/resources/soxs_keywords.yaml` has none of these aliases.

### How ESPRESSO does it

From the ESPRESSO pipeline user manual, issue 3.6.0, 2026-07-27 ([ftp.eso.org](https://ftp.eso.org/pub/dfs/pipelines/instruments/espresso/espdr-pipeline-manual-3.6.0.pdf)):

- Section 11.1.9 "Barycentric Correction" (pages 65 to 66): a dedicated code computes the observer's velocity projected on the line of sight "with the accuracy of at least 1 cm/s", the barycentric arrival time, and a yearly upper bound (`QC BERV`, `QC BJD`, `QC BERVMAX`). The inputs are the ICRS target position, its proper motion, the UTC date of observation, and the observer's longitude, latitude and altitude.
- The science recipe steps (page 72 to 73) are, in order: "Compute barycentric correction using FITS header information", "Correct wavelength solution from BERV, shifting it to the barycenter of the SS", "Rebin and merge orders", then "Flux calibrate S1D spectrum ... using absolute efficiency curve". So ESPRESSO shifts the wavelength solution before resampling, and its flux calibration runs on the barycentric scale. This is option B.
- ESPRESSO reduces one exposure at a time, so it has no stack problem. The manual also measures a "flux-weighted mid-exposure time on exposure-meter data" for QC. It does not say which epoch feeds the BERV. SOXS has no exposure meter.

### How PypeIt does it

From the PypeIt source at commit `c1798aa8` (branch `release`):

- The parameter is `calibrations.wavelengths.refframe`, with options `observed`, `heliocentric` and `barycentric`. The default is `heliocentric` (`pypeit/par/pypeitpar.py`, lines 3529 to 3534 and 3680), with the note "TODO: What should the default be? None or 'heliocentric'?".
- `refframe_correct` in `pypeit/pypeit_steps.py` (lines 1171 to 1228) computes one velocity per frame and multiplies both the 1D `BOX_WAVE`/`OPT_WAVE` arrays (`SpecObj.apply_helio`, `pypeit/specobj.py`, lines 647 to 667) and the 2D wavelength image by the factor. It runs after extraction and before coadding, so each exposure is corrected before the 1D spectra are combined. This is option E.
- The epoch is `fitstbl.construct_obstime(frame)`, which is `Time(self['mjd'][row], format='mjd')` (`pypeit/metadata.py`, lines 522 to 533). That is the frame's `mjd` metadata, normally the exposure start, not the midpoint.
- The velocity comes from PypeIt's own `geomotion_velocity` (`pypeit/core/wave.py`): the Earth's barycentric velocity plus the observatory's GCRS velocity, projected on the target direction. A code comment asks why it does not use `SkyCoord.radial_velocity_correction`. The factor is relativistic: `vel_corr = np.sqrt((1. + vel/299792.458) / (1. - vel/299792.458))`.
- `collate_1d` has a separate `refframe` option to correct spec1d files before coadding, and it refuses to correct a spectrum twice (`pypeit/collate.py`, lines 616 to 653).

### The astropy tool

`SkyCoord.radial_velocity_correction(kind='barycentric', obstime=None, location=None)` in `astropy` 7.2.2, the version `soxspipe` pins (`pyproject.toml`: `astropy==7.2.2`). The source is `astropy/coordinates/sky_coordinate.py`, lines 1520 to 1723. Facts from its docstring and source:

- It returns the correction to add to a measured radial velocity. The correction is multiplicative: `rv = rv + vcorr + rv * vcorr / c`, in the optical convention `zb = vcorr / c` (Wright and Eastman 2014). The docstring says to use barycentric rather than heliocentric "if m/s precision is required".
- `obstime` and `location` come from the arguments or the `SkyCoord` frame attributes. It raises `ValueError` if a value is given in both places, and `TypeError` if it is in neither.
- The barycentric branch is equation 28 of Wright and Eastman (2014): Lorentz factor, gravitational redshift at the site, no Shapiro delay. If the `SkyCoord` carries proper motion and an `obstime`, the method moves the target to the observation epoch and adds the star's own motion terms.
- It uses the built-in solar system ephemeris by default. JPL ephemerides are optional through `solar_system_ephemeris`.

Inputs and where SOXS raw frames have them (V818 Sco, `/Volumes/SOXS/raw_archive/soxs_phase3_dataset/soxs_offset_V818Sco_raw/SOXS.2025-05-13T04:08:15.092.fits.Z` and the NIR frame at `04:08:50.132`):

| Input | Header cards | V818 Sco value | Note |
|---|---|---|---|
| Site | `ESO TEL GEOLAT`, `ESO TEL GEOLON`, `ESO TEL GEOELEV` | −29.2584 deg, −70.7345 deg, 2377.0 m | No SOXS keyword alias yet. `EarthLocation.of_site("La Silla Observatory")` gives −29.2567 deg, −70.73 deg, 2347 m. The velocity changes by 0.03 m/s. |
| Start time | `MJD-OBS` (`MJDOBS` alias) | 60808.1723969 | UTC. `MJD-END` is not in the raw frames (`# MJDEND` is commented out in `soxs_keywords.yaml`, line 2). |
| Exposure length | `EXPTIME` (UV-VIS); `ESO DET SEQ1 EXPTIME` (NIR, alias `SEQ_EXPTIME`) | 60.0 s; 63.272 s | NIR `EXPTIME` is one DIT (DY-1666 comment in `soxs_keywords.yaml`, line 10). |
| Target | `ESO TEL TARG ALPHA`/`DELTA` (packed), `EQUINOX`, `ESO TEL TARG PMA`/`PMD`, `ESO TEL TARG EPOCH` | 161955.069, −153825.017, J2000, −0.007185, −0.012332 | DY-2296 chose `TEL TARG ALPHA`/`DELTA` over the `RA`/`DEC` pointing cards. For the velocity it does not matter: a 104 arcsec error changes a 30 km/s correction by at most 15 m/s. |

Results for this frame (exposure midpoint, ICRS position from `RA`/`DEC`): barycentric 7.3811 km/s, heliocentric 7.3782 km/s, difference 2.95 m/s.

### How large the correction is, compared with SOXS

These numbers come from `astropy` 7.2.2 for the La Silla site above (script run in the session scratchpad). The pixel scales come from real `soxspipe` products.

- Annual range: over a grid of 36 × 13 sky positions (RA every 10 deg, Dec −90 deg to +30 deg) and 74 dates across a year, the correction ranges from −30.40 km/s to +30.41 km/s.
- Diurnal term: the site's rotation speed is 0.406 km/s. The fastest change is about 0.107 km/s/h from rotation plus about 0.02 km/s/h from the orbit. For V818 Sco the largest change over 3 h was 0.337 km/s, with a peak rate of 0.123 km/s/h.
- Real sequence: the V818 Sco offset set in `soxs_phase3_dataset` (2 science frames per arm, 0.03 h apart) has a correction of 7.377 km/s to 7.381 km/s. It drifts by 3 m/s, and the correction at the mean epoch equals the mean of the per-frame corrections to better than 0.01 m/s.
- Timing: a 1 min error in the epoch changes the correction by about 2 m/s. Using the exposure start instead of the midpoint of a 3600 s exposure moves it by at most about 0.065 km/s.
- Resolving power: SOXS "is a medium resolution spectrograph (R ∼ 4500 for a 1 arcsecond slit)" ([Rubin et al. 2020, "Progress on the UV-VIS arm of SOXS"](https://arxiv.org/abs/2012.12681), section 1). The optical design paper gives "R>4500" ([Zanmar Sanchez et al. 2018](https://arxiv.org/abs/1809.01521), abstract). `soxspipe` does not measure the resolving power today, so these are design values. A resolution element at R = 4500 is c/R = 66.6 km/s. Narrower slits give a higher R and a smaller element.
- Native pixels, as the median `pixelScaleNm` per order in real `_EXTRACTED_ORDERS` products:
  - UV-VIS (`20250907T050927_VIS_1X1_1_STARE_SLIT1.5_900.0S_SOXS_SAXJ1808_EXTRACTED_ORDERS.fits`): 0.023 nm to 0.047 nm, 17.3 km/s to 19.0 km/s.
  - NIR (`20260128T072009_NIR_3_STARE-OBJ_SLIT1_0_600_0S_SOXS_PAQS135626+042348_EXTRACTED_ORDERS.fits`): 0.037 nm to 0.087 nm, 13.3 km/s to 13.6 km/s.
- Merged-spectrum bins: 0.02 nm for UVB/VIS and 0.06 nm for NIR (`horne_extraction.py`, lines 2403 to 2407). That is 17.1 km/s to 7.1 km/s per bin across 350 nm to 850 nm, and 22.5 km/s to 8.9 km/s per bin across 800 nm to 2030 nm.

| Quantity | UV-VIS | NIR |
|---|---|---|
| Maximum correction | 30.4 km/s | 30.4 km/s |
| As wavelength | 0.035 nm at 350 nm, 0.085 nm at 850 nm | 0.081 nm at 800 nm, 0.20 nm at 2000 nm |
| Resolution elements (R = 4500) | 0.46 | 0.46 |
| Native pixels | about 1.6 | about 2.3 |
| Merged bins | 1.8 to 4.3 | 1.4 to 3.4 |
| Drift over a 3 h sequence | ≤ 0.34 km/s (≤ 2 % of a pixel) | ≤ 0.34 km/s (≤ 3 % of a pixel) |

So a skipped correction is a systematic shift of up to half a resolution element. It is large for radial-velocity work, and it is the error that `SPEC_SYE` would have to describe for a topocentric product. The precision of the correction itself (mm/s to m/s) is far below anything SOXS can measure.

### Time scales and offline use

- The site velocity needs UT1, so `astropy` reads the IERS Earth rotation tables. With `iers.conf.auto_download = False`, a correction for a date after the IERS-B table that `soxspipe`'s environment has (last MJD 61182, 2026-05-22) raises `ValueError: interpolating from IERS_Auto using predictive values that are more than 30.0 days old`. The test date 2025-05-13 worked offline.
- Setting `astropy.utils.iers.conf.auto_max_age = None` lets the call run and gives `ErfaWarning` "dubious year" and an `AstropyWarning` about polar motion. The `iers_degraded_accuracy` setting is the documented way to choose between error, warn and ignore. The error this causes is tiny: a UT1 error of 0.9 s moves the diurnal term by about 0.03 m/s.
- So for ESO's in-house processing, a build must either make sure that current IERS data are present or accept degraded IERS accuracy on purpose. A network call must not decide whether a reduction fails.

### What soxspipe already has that could carry the correction

- The wavelength chain: per-pixel wavelengths (`wavelengthMean`) are tuned to the sky lines in `_correct_wavelengths` (`horne_extraction.py`, line 929, called at line 637). This happens on the observed scale, as it must. `merge_extracted_orders` (line 2364) then builds a fixed grid with `_resample_grid` (line 2338: multiples of the step inside the data range) and resamples with `np.interp` (line 2504). Option B would act between those two steps.
- Rounding: the written merged table rounds `WAVE` to 2 decimals (`horne_extraction.py`, line 729), and so does the nod stack (`soxs_nod.py`, line 2043). On the round 0.02 nm and 0.06 nm grids this rounding loses nothing. With option C, it would.
- The stare stack: `soxs_stare._reduce_object_frames` stacks the 2D frames with `clip_and_stack` (`soxs_stare.py`, line 673) before one extraction. The stacked frame takes the header of the first frame (`base_recipe.py`, line 1980), so its `MJD-OBS` is the start of the first frame. The stack is a sigma-clipped mean, so its effective epoch is the mean of the frame midpoints. The recipes already require equal exposure times for a response curve (`sequence_airmass.require_equal_exposure_times`).
- The nod and offset stacks: by default the A and B frames are aligned and stacked before extraction (`stack_before_extraction: True` in `soxs_default_settings.yaml`; `soxs_nod._should_align_and_stack`, line 449). With `--no-stack`, each cycle is extracted and `stack_extractions` (`soxs_nod.py`, line 1945) concatenates the merged spectra and takes the median per exact `WAVE` value (lines 2002 to 2003). That works only because every cycle lands on the same round grid. A per-cycle correction must either act before the merge resampling (option B) or resample to a common grid before this group-by. `soxs_offset` inherits the same code (`soxs_offset.py`, line 31).
- The per-frame headers that a mean epoch needs are kept: `soxs_nod.produce_product` copies `self.scienceHeaders` from all object frames. The airmass code already uses the same pattern: `sequence_airmass.mean_airmass` takes the mean over every frame's headers.
- Flux calibration: `flux_calibration.calibrate` evaluates the extinction factor and the response at the spectrum's own `WAVE` values (`flux_calibration.py`, lines 222 and 229; the response interpolation is `_response_factors_from_table`, line 28). The response is a least-squares cubic spline with 30 nm to 40 nm knots and masked telluric regions (`response_fit.py`; `knot_spacing_nm` in `soxs_default_settings.yaml`, lines 1186 and 1221). I measured the effect of evaluating it 30 km/s off:
  - VIS response (`/Volumes/SOXS/reductions/P32.11_response_fit_ab/v6/.../20251026T070308_VIS_1X1_1_NOD_STD_FLUX_SLIT5_0_300_0S_SOXS_RESP.fits`): at most 0.061 %, median 0.028 %.
  - NIR response (`.../dy1647_ps_40_1/.../20260115T014023_NIR_3_NOD_STD_FLUX_SLIT5_0_60_0S_SOXS_RESP.fits`): at most 0.035 %, median 0.017 %.
  - La Silla extinction curve (`extinction_lasilla.dat`): at most 0.02 % at airmass 1 and 0.04 % at airmass 2, both at 360 nm.

  These are 100 times smaller than a typical `FLUXERR`. So the order of the flux calibration and the barycentric correction is a question of clarity, not of accuracy. Calibrating first on the observed scale (options A, C and D) keeps the response and the science spectrum on the same scale by construction. The response itself should be built from standards on the observed scale: it is an instrument property.
- Header writing: `phase3.basic_header_scrubbing` and `write_fits_table_to_disk` (`phase3.py`, lines 13 and 137) are where the Phase 3 cards are set, and where `SPECSYS` and the velocity card would go.
- Telluric: no science telluric correction exists in the code. The response fit masks telluric regions only. If one is ever added, it must run on the observed scale, before any barycentric shift.

## Open points for DY-2303

- Whether the map ruling "`SPECSYS = 'BARYCENT'`" stands, now that the standard allows `TOPOCENT` and the ESO X-Shooter collection uses it (option A).
- Where the correction goes if it is applied: before the order merge (B), or on the final product (C or D).
- Which epoch to use for a stack: the mean of the frame midpoints (recommended by the numbers above), `TMID`, or the first frame's `MJD-OBS`. The difference is at most about 0.2 km/s.
- Which header cards record the velocity. Options: reuse the ESO `QC VRAD BARYCOR`/`HELICOR` names that X-Shooter writes, or add a soxspipe card. Whether to record the epoch used.
- Whether the internal `_EXTRACTED_MERGED` and efficiency products stay topocentric. They do under options A, C and D.
- How to handle IERS data offline: ship current IERS data or accept degraded accuracy deliberately.
- Whether to use the header site cards or a fixed La Silla location. They differ by 0.03 m/s, but the header cards need no network and follow any update of the site values.
