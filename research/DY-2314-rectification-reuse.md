# DY-2314: what the existing rectification code can produce for a 2D spectrum

This note answers the research ticket DY-2314 for the design ticket DY-2315 ("Design of the rectified 2D spectrum"). It is for the person who designs the `ANCILLARY.2DSPECTRUM` product. All `file:line` references are to `develop` at `6bdef45d`. The timing figures come from one small run on real SOXS commissioning products (see [Timing and memory](#timing-and-memory)).

## Contents

- [Short answer](#short-answer)
- [What `get_order_rectified()` returns in Horne extraction](#what-get_order_rectified-returns-in-horne-extraction)
- [The grid](#the-grid)
- [Sky subtraction](#sky-subtraction)
- [How long the layers live](#how-long-the-layers-live)
- [Stare, nod and offset: which frame is rectified](#stare-nod-and-offset-which-frame-is-rectified)
- [The `soxs_straighten` recipe](#the-soxs_straighten-recipe)
- [Putting the orders on one grid](#putting-the-orders-on-one-grid)
- [Flux calibration of a 2D array](#flux-calibration-of-a-2d-array)
- [Timing and memory](#timing-and-memory)
- [Verdict on reuse](#verdict-on-reuse)
- [Open questions for DY-2315](#open-questions-for-dy-2315)

## Short answer

The Horne layers can be kept after extraction at almost no memory cost, but as they are they make a poor 2D spectrum:

- They cover only the extraction aperture (10 to 15 detector pixels, about 3 to 4 arcsec), not the 12 arcsec slit.
- Each order has its own irregular wavelength grid and its own slit scale, and the slit axis follows the object trace.
- Their wavelength layer does not include the skyline wavelength correction that the 1D spectrum gets.
- In the default nod mode there are two sets of layers (A−B and B−A), each on its own grid.

The rectification machinery (`image_transformer`) is the right tool. A second transformer on a full-slit, fixed-centre grid gives a proper 2D spectrum, but its set-up memory grows with the slit length. The 15-row Horne grid already peaks at 6.3 GB on VIS and 12.8 GB on NIR. A full-slit grid peaked at 13.7 GB on VIS and would need about 40 GB on NIR, more than the Mac has. The memory use of `_precompute_resampling_weights()` must come down before a full-slit grid can run in every science reduction.

## What `get_order_rectified()` returns in Horne extraction

`horne_extraction.extract()` builds one `image_transformer` (`soxspipe/commonutils/horne_extraction.py:561-570`) and caches these images on it:

| Layer | Source | Where |
|---|---|---|
| `fluxRaw` | the sky-subtracted frame, with its bad-pixel mask | `horne_extraction.py:574` |
| `variance` | the square of the frame's uncertainty, with linear (not squared) area weights | `horne_extraction.py:581-583` |
| `bpMask` | made by `cache_image()` when a mask is given: True where more than 0.2 of the cell area is flagged | `image_transformer.py:238-241` |
| `fluxSky` | the "subtracted frame": the sky model (stare) or the other nod or OFF frame (nod, offset). No mask | `horne_extraction.py:584-585` |
| `wavelength`, `slit` | computed from each order's bin edges, not resampled from the 2D map | `image_transformer.py:291-336` |

`get_order_rectified()` (`image_transformer.py:778-836`) returns a list with one dict per order. Each dict maps a layer name to a 2D `numpy` array of shape (slit rows, wavelength columns). It stacks the per-column lists that `cache_image()` stored in each order's `DataFrame` (`image_transformer.py:794`).

The extraction then adds more layers to the same dicts, in place:

- `mask`: `bpMask` or non-finite flux or a bad wavelength (`generate_masks()`, `horne_extraction.py:3555`).
- `objectProfile`: the normalised Horne profile (`fit_object_profile()`, `horne_extraction.py:3721`).
- `unsupported`: pixels with no support for the profile fit (`horne_extraction.py:3714-3716`).

Because the extraction runs serially (`turnOffMP=True`, `horne_extraction.py:621`; the serial path of `fundamentals.fmultiprocess` calls the function on the same objects), these additions are visible in the caller's `orderRectifiedImages` list after the loop. A boxcar fallback (no object trace) adds `mask` only.

The `variance` layer sums each detector pixel's variance with its linear overlap area. This is a deliberate approximation that is correct for sums across neighbouring cells (the comment at `horne_extraction.py:575-580`). For a single rectified cell it overestimates the variance, because neighbouring cells share detector pixels. A 2D product would carry that per-cell variance as it is, or document it.

## The grid

The grid is set per order, and it is not fixed:

- **Wavelength axis.** The bin edges are evenly spaced in detector pixels along the trace, so they are not evenly spaced in wavelength (`image_transformer.py:41-49`, `_pixel_uniform_wavelength_edges()` at `image_transformer.py:1015-1049`). Each rectified column holds about one detector pixel. In the VIS run the step was 0.019 to 0.026 nm in order 2 and 0.037 to 0.053 nm in order 4. Wavelengths are in nm (the 2D map's units). The debug plot labels at `image_transformer.py:269` and `:829` say "Å", which is wrong.
- **Slit axis.** There are `2 × slitHalfLength` rows (rounded half up), one detector pixel each (`image_transformer.py:680-684`). The slit edges are in arcsec, scaled by that order's own arcsec-per-pixel, measured from the dispersion map (`_arcsec_per_pixel()`, `image_transformer.py:1052-1071`). So the arcsec step differs between orders.
- **Slit centre.** Each order is centred on a polynomial fit (degree 3) of the object trace's slit position against wavelength (`_order_slit_centre()`, `image_transformer.py:599-643`). The `slit` layer is the offset plus that trace centre (`image_transformer.py:324`), so a row is a fixed offset from the object, not a fixed slit position. With `fixedSlitCentreArcsec` set, every order is centred on one slit position instead (`image_transformer.py:570-571`, `:615-616`). Horne extraction uses that only for its no-trace fallback (`SLIT_CENTRE_ARCSEC = 0.0`, `horne_extraction.py:54`, `:569`).
- **Aperture.** `slitHalfLength` is half of `horne-extraction-slit-length` (`horne_extraction.py:306`), divided by the binning along the slit (`horne_extraction.py:364-369`). The default settings are 15 pixels for `soxs-stare` (`soxspipe/soxs_default_settings.yaml:546`, `:666`), 11 for `soxs-nod` (`:903`) and 10 for `soxs-offset` (`:1006`). The slit is 12 arcsec long (`soxspipe/resources/soxs_detector_parameters.yaml`, `slit_length`), so the layers hold about a quarter to a third of it.
- **Oversampling.** The resampling uses a 5 × 5 sub-grid per output cell (`zoomFactorSlit`, `zoomFactorWavelength`, `image_transformer.py:131-132`), then sums each block back to detector resolution (`_rebin_resampling_weights()`, `image_transformer.py:1074-1115`).

## Sky subtraction

Yes, `fluxRaw` is sky-subtracted, because it is rectified from `skySubtractedFrame` (`horne_extraction.py:574`):

- **Stare:** the frame after `_subtract_sky()` (`soxspipe/recipes/soxs_stare.py:778-780`, `:802`). The sky model is the `fluxSky` layer.
- **Nod and offset:** the pair difference (A−B, B−A, or ON−OFF). There is no sky model: `fluxSky` is the other frame (`soxspipe/recipes/soxs_nod.py:1862`, `:1871`), and `carrySkyCounts=False` keeps it out of the 1D product (`soxs_nod.py:1929`).

The science 1D spectrum that is flux-calibrated comes from the flat-fielded frame (`soxs_stare.py:820-829`), and so do these layers.

## How long the layers live

The layers live only inside one call of `horne_extraction.extract()`. `orderSlices`, `orderRectifiedImages` and `transformer` are local variables (`horne_extraction.py:587-600`). Nothing stores them on `self` or writes them to disk, so they are freed when `extract()` returns. `extract_single_order()` returns only the 1D columns (`horne_extraction.py:3094-3107`).

Keeping them is cheap. For a 15-pixel aperture the five cached layers take 7.3 MB on VIS and 14.8 MB on NIR (measured). The added `mask`, `objectProfile` and `unsupported` layers take about the same again. The cost is in building the transformer, not in holding its output (see [Timing and memory](#timing-and-memory)).

## Stare, nod and offset: which frame is rectified

| Recipe and mode | Frames rectified by Horne extraction | Where |
|---|---|---|
| `soxs_stare` | one sky-subtracted, flat-fielded stack | `soxs_stare.py:800-820` |
| `soxs_nod`, stacked (default `stack_before_extraction: True`), no shift-and-add plan | two: A−B (object at A) and B−A (object at B), each on its own grid | `soxs_nod.py:1862-1877` |
| `soxs_nod`, stacked, with a shift-and-add plan | one: B−A shifted onto A−B and averaged | `soxs_nod.py:1598-1630`, `:1862` |
| `soxs_nod`, `--no-stack` or setting off | two per nodding cycle | `soxs_nod.py:1059-1167` |
| `soxs_offset`, stacked | one: stacked ON − stacked OFF | `soxspipe/recipes/soxs_offset.py:627-655`; B−A is skipped because the recipe name has no "nod" (`soxs_nod.py:1866`) |
| `soxs_offset`, per cycle | one ON−OFF per cycle | `soxs_offset.py:497-606` |

A standard star always takes the stacked path (`soxs_nod.py:450-472`). A plan is used only when both nod positions are measured and the shift clears the Horne aperture (`soxs_nod.py:648-655`).

The 1D products are combined after extraction. `stack_extractions()` takes the median per wavelength of the A and B spectra, and of the cycles (`soxs_nod.py:2024`). A 2D product would need the same median of 2D layers, which needs one shared grid first. A−B and B−A grids differ slightly, because each grid follows its own trace across the detector.

With an 11-pixel nod aperture centred on the positive trace, the negative beams fall outside the layers. A full-slit 2D nod product shows the usual positive and negative traces.

The extraction also runs a second time, on the not-flattened frame, when a response is built (`soxs_stare.py:1604-1624`, `soxs_nod.py:1243-1263`, `soxs_offset.py:668-683`). Only standard stars do this. The science 2D spectrum would come from the flat-fielded pass only.

## The `soxs_straighten` recipe

It produces nothing. `produce_product()` returns `productPath = None` after writing the QC report (`soxspipe/recipes/soxs_straighten.py:194-224`), and its input check is mostly commented out (`soxs_straighten.py:155-176`). No recipe, the reducer, the CLI or the SOF maps call it. Only its import (`soxspipe/recipes/__init__.py:12`), its docs page (`docs/source/recipes/soxs_straighten.md`) and a settings block in `soxspipe/test_settings_xsh.yaml:222` refer to it. The settings `straighten_grid_res_wavelength` and `straighten_grid_res_split` in its docs page are not read by any code.

The only other caller of `get_order_rectified()` is a debug block in `create_dispersion_map.map_to_image()` (`soxspipe/commonutils/create_dispersion_map.py:3514-3528`), which rectifies a multi-pinhole frame with a 20-pixel half-length and discards the result.

## Putting the orders on one grid

The 1D merge is `merge_extracted_orders()` (`horne_extraction.py:2351-2554`):

1. It orders the orders red to blue and puts each join at the middle of the overlap (`horne_extraction.py:2401-2408`). UVB order join 20-21 is fixed at 363.5 nm (`:2413-2414`).
2. Where the overlap is wider than 0.84 nm (VIS) or 2.52 nm (NIR), it drops each order's points more than 20 output steps past the join (`horne_extraction.py:2418-2430`). Points in the remaining overlap from both orders stay.
3. It drops points whose variance is far above the local median (`_variance_keep_mask()`, `horne_extraction.py:2440`).
4. It sorts all points by wavelength, removes duplicate wavelengths and linearly interpolates flux and variance onto a fixed grid: 0.02 nm for UVB and VIS, 0.06 nm for NIR, aligned to multiples of the step (`horne_extraction.py:2390-2394`, `:2325-2349`, `:2491-2504`). There is no weighting between orders in the overlap: the interpolation uses whichever points fall next to each output wavelength.

This logic can be applied row by row to 2D layers, because each slit row is a 1D spectrum. It needs two things first:

- **One slit grid for all orders.** The Horne layers do not have one: the arcsec step differs per order and the rows follow each order's trace polynomial. A transformer run with `fixedSlitCentreArcsec` and one common arcsec-per-pixel would give one, or each order must be resampled along the slit.
- **A 2D form of the variance clip.** `_variance_keep_mask()` works on a single 1D sequence. Per row it is the same function, but a clipped pixel would then leave a hole in that row only.

The merge grid does not match the native grids. On VIS the 0.02 nm step is close to the native step of order 2 (0.019 to 0.026 nm) and finer than orders 1, 3 and 4 (up to 0.053 nm). On NIR the 0.06 nm step is coarser than the blue orders (0.033 nm) and finer than the red ones (0.094 nm). So a merged 2D product oversamples some orders and undersamples others, as the 1D product already does.

## Flux calibration of a 2D array

Yes, the factors apply column by column in the same way. The 1D formula is `counts / exptime × extinction × response × 1e-17` (`_calculate_flux_calibration()`, `soxspipe/commonutils/flux_calibration.py:22-29`), where:

- `extinction` and `response` are functions of wavelength only (`flux_calibration.py:274-284`, `load_response_factors()` at `:109-136`). The extinction factor applies to UVB and VIS only.
- `counts` are flat-fielded electrons per unbinned pixel: the merged `FLUX_COUNTS` divided by the binning along the dispersion axis (`counts_per_unbinned_pixel()`, `soxspipe/commonutils/response_basis.py:67-81`).
- The merged `FLUX_COUNTS` are per native pixel, not per nm: they are interpolated pixel sums (`horne_extraction.py:2491-2509`). Only the not-flattened spectrum gets a per-nm column (`horne_extraction.py:2529-2539`).

A rectified cell also holds the counts of about one detector pixel (`image_transformer.py:51-53`). So for column `j` at wavelength `λj` the same scalar `exptime^-1 × extinction(λj) × response(λj) × 1e-17 / binning` multiplies every slit row, and the error is `sqrt(variance)` times the same factor (`flux_calibration.py:296-301`). The result is in erg/cm²/s/Å per slit pixel. The sum over the aperture rows matches the boxcar 1D spectrum, not the optimal one.

Two things carry over from the 1D code:

- The response is NaN outside its wavelength range (`flux_calibration.py:81-106`). The same columns become NaN in 2D.
- For nod, the exposure time is that of one A frame and the airmass is the mean of all A and B frames (`soxs_nod.py:1328-1329`). The 2D A and B layers must be combined by median, as the 1D is, for that exposure time to stay correct.

The calibration should be applied after the wavelengths get the skyline correction. `_correct_wavelengths()` shifts only the 1D `wavelengthMean` column after extraction (`horne_extraction.py:643`, `:916-941`). The 2D wavelength layer does not get the shift, so the 2D product must apply the same per-order shift, or the NIR global-fit correction (`apply_sky_corrections`, `horne_extraction.py:974-992`), before it is resampled.

## Timing and memory

I timed `image_transformer` alone on real products, with no reduction. The script builds the transformer with a fixed slit centre (so no object trace is needed), caches `fluxRaw` and `variance`, and calls `get_order_rectified()`. It ran with the `soxspipe` conda environment (Python 3.12) against the worktree code.

Inputs, from `/Volumes/SOXS/reductions/2026_jan_commissioning/reduced/`:

- VIS: frame `2026-01-10/soxs-stare/20260111T014559_VIS_1X1_1_STARE_STD_FLUX_SLIT5_0_300_0S_SOXS_SKYSUB.fits` (4096 × 808 pixels) and map `2026-01-10/soxs-spat-solution/20260111T095127_VIS_1X1_1_SSOL_MULTPIN_30_0S_SOXS`.
- NIR: frame `2026-01-10/soxs-stare/20260111T014557_NIR_3_STARE_STD_FLUX_SLIT5_0_150_0S_SOXS_SKYSUB.fits` and map `2026-01-11/soxs-spat-solution/20260112T100055_NIR_3_SSOL_MULTPIN_15_0S_SOXS`.

| Arm | Slit rows | Build transformer | Cache 2 images | Peak RSS | Layer memory (5 layers) |
|---|---|---|---|---|---|
| VIS | 15 | 8.8 s | 0.1 s | 6.3 GB | 7.3 MB |
| VIS | 44 (about 12 arcsec) | 30.5 s | 0.2 s | 13.7 GB | 21.4 MB |
| NIR | 15 | 18.5 s | 0.1 s | 12.8 GB | 14.8 MB |

VIS order shapes at 15 rows: (15, 3589), (15, 3670), (15, 3607) and (15, 3873) for orders 1 to 4. NIR has 15 orders (10 to 24), each about 2040 columns (order 24 has 1242), with a native step from 0.033 nm (order 23) to 0.094 nm (order 10).

Findings:

- Caching another image costs almost nothing once the transformer exists. The cost is in `_precompute_resampling_weights()` (`image_transformer.py:339-484`), which converts every sub-cell corner of every order to detector pixels in one `DataFrame` (`image_transformer.py:384-401`). That table has `rows × 5 × columns × 5 × 4` entries per order.
- Time and peak memory grow about linearly with the slit length (VIS: 15 to 44 rows took 3.5 times the time and 2.2 times the peak). NIR at 15 rows already peaks at 12.8 GB, because it has about twice as many columns in total as VIS.
- The 15-row peaks (6.3 GB VIS, 12.8 GB NIR) are already paid inside every Horne extraction. The NIR figure matches the known 13 GB peak of one NIR stare reduction (soxspipe General `CLAUDE.md`, "Machine limits"), so the Horne transformer is very likely that peak.
- I did not run a full-slit NIR grid. Scaling the NIR figure by the VIS ratio gives about 40 GB, which the Mac does not have.
- The VIS 44-row run was not wrapped in `realdata_slot.sh`; the NIR run was. Nothing else was running at the time of the VIS run.

## Verdict on reuse

**Reuse the machinery, not the Horne layers as they are.**

Reusing the Horne layers as they are gives an aperture cut-out (3 to 4 arcsec) that follows the object, on per-order irregular grids, with uncorrected wavelengths, and with two sets of layers per nod. Keeping them costs a few MB and no run time. But turning them into a merged 2D product needs a slit resampling per order, the skyline shift, a resampling onto the merge grid, and a median of A−B and B−A. It also leaves out the nod negatives and any extended emission beyond the aperture. It is acceptable only if DY-2315 decides that the 2D spectrum is a "what the extraction saw" product.

A proper full-slit 2D spectrum needs a second `image_transformer` per extracted frame, with `fixedSlitCentreArcsec` and a slit length that covers the slit. This gives straight rows at fixed slit positions and one arcsec scale per order. It also uses the existing `cache_image()` path for flux, variance and mask, and needs no new resampling code. Its costs and blockers:

- **Blocker: memory.** `_precompute_resampling_weights()` must build its corner table per order (or in chunks) instead of all orders at once, or use a smaller sub-grid for this grid. Without that change a full-slit grid does not fit on NIR at all (about 40 GB). The same change would also lower the existing Horne peak.
- **Cost: run time.** About 30 s extra per VIS frame at full slit with the current code, and roughly 60 s per NIR frame by the same scaling (not measured). That repeats for each nod difference and each cycle in `--no-stack` mode.
- **Needed in both cases:** the per-order skyline shift on the wavelength layer, row-by-row order merging onto the 1D grid, a median of nod differences, and column-wise flux calibration.

## Open questions for DY-2315

1. Full slit or aperture? Full slit shows the nod negatives and extended emission but needs the memory fix. The aperture is free but follows the object.
2. Which grid? The 1D merge grid (0.02 nm or 0.06 nm, linear, merged orders) matches the 1D product. Per-order native grids avoid a resampling but need one image extension per order.
3. Nod default mode: median of A−B and B−A on the shared grid (as the 1D), or the A−B difference only, which shows both beams?
4. `--no-stack` mode: is a 2D product made per cycle and then combined, or only in stacked mode?
5. Variance: keep the linear-weight variance of each cell (it overestimates single-cell noise), or compute a squared-weight variance for the 2D product?
6. Is a 2D spectrum written for a run without a response (counts only), as the 1D preview plots uncalibrated runs in their own units?
7. Where should the 2D spectrum be built: inside `horne_extraction.extract()` (where the trace, the sky shifts and the frame are at hand) or in a new module called by each recipe after its stack?
