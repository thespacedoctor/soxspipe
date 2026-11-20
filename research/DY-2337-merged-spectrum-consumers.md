# DY-2337: code that reads `_EXTRACTED_MERGED` and `_FLUXCAL` products

This note lists every place in `soxspipe` (develop at `8b039737`, 2026-11-20) that writes or reads the order-merged spectrum (`<stem>_EXTRACTED_MERGED.fits`) or the flux-calibrated spectrum (`<stem>_FLUXCAL.fits`), their column names, their suffixes, or their `PRO CATG` and SOF tags. It says what must change when one merged spectrum (`<SOF stem>_SPECTRUM.fits`) replaces both, with the DY-2293 column set and the DY-2326 fixed log-λ grid. It is for whoever splits the Phase 3 build into implementation issues (DY-2338). All paths are relative to the repository root.

## Contents

- [Short answer](#short-answer)
- [What the two files hold today](#what-the-two-files-hold-today)
- [Response, efficiency and SNR: the DY-2326 claim](#response-efficiency-and-snr-the-dy-2326-claim)
- [Consumers that must change](#consumers-that-must-change)
- [Consumers that need no change](#consumers-that-need-no-change)
- [The data organiser and `product_frames`](#the-data-organiser-and-product_frames)
- [Tests and fixtures](#tests-and-fixtures)
- [Docs](#docs)
- [Scripts](#scripts)
- [Unverified](#unverified)
- [Sources](#sources)

## Short answer

- The DY-2326 claim holds for the logic of the response, efficiency and SNR code: none of it depends on equal bins, and every one of them takes its file paths from the caller. It does not hold in three details. (1) The response reads `FLUX_COUNTS` from the standard star's file on disk, so the merged spectrum of an uncalibrated standard-star run must keep a column named `FLUX_COUNTS`. The DY-2293 rule "uncalibrated runs put the counts in `FLUX` and omit `*_COUNTS`" must apply only to the `_ph3` file. (2) The efficiency reads `FLUX_DENSITY_COUNTS` from the `_EXTRACTED_MERGED_NOTFLAT.fits` file, which is not in the DY-2293 column set; that file must stay, or the column must move. (3) Flux-conserving summation changes what a `FLUX_COUNTS` value means (counts per output bin, not per native pixel), so `RESPONSE_BASIS` in `soxspipe/commonutils/response_basis.py:22` must change value so that old RESP products are refused. The density columns (`FLUX_DENSITY_COUNTS`, and arguably the flat-fielded counts) need normalized weights, not plain summed weights.
- The writers are three places: `horne_extraction.extract` (`_EXTRACTED_MERGED`), `flux_calibration.calibrate` (`_FLUXCAL`) and `soxs_nod.stack_extractions` (stacked `_EXTRACTED_MERGED`, which `soxs_offset` inherits). The readers that use the suffix or the columns are the response function, the QC-plot re-read of `_FLUXCAL` in the three recipes, the ASCII writer, `predict_product_path` (skip-if-exists and the recipe log name), and the data organiser's predicted product name.
- `stack_extractions` must group by bin index instead of `WAVE` rounded to 4 decimals, and the 2-decimal and 3-decimal rounding of `WAVE`, `FLUX_COUNTS`, `SNR` and `FLUX_DENSITY_COUNTS` on write must go in both `horne_extraction.py:735` and `soxs_nod.py:2063`.
- `product_frames` registers one `OBJECT_TAB_<ARM>` row per science SOF, with `file = <SOF stem>_EXTRACTED_MERGED.fits`. `OBJECT_TAB` is not the `PRO CATG` the file carries (`SCI_SLIT_FLUX_<ARM>`), and `OBJECT_TAB_<ARM>` is the real `PRO CATG` of the object-trace table (`_OBJTRACE.fits`). No code reads the `OBJECT_TAB` rows today. The `filepath` of every `*_obj` row has a wrong directory, built at `soxspipe/commonutils/data_organiser.py:4088`.

## What the two files hold today

The table below lists the files and their columns as the code writes them.

| File | Writer | `PRO CATG` | Columns | Rounded on write |
|---|---|---|---|---|
| `<stem>_EXTRACTED_MERGED.fits` (stare, single extraction) | `soxspipe/commonutils/horne_extraction.py:720-768` | `SCI_SLIT_FLUX_<ARM>` (`horne_extraction.py:677`) | `WAVE`, `FLUX_COUNTS`, `VARIANCE`, `SKY_COUNTS` (stare only, `horne_extraction.py:2510-2521`), `SNR` | `WAVE` 2, `FLUX_COUNTS` 3, `SNR` 2 decimals (`horne_extraction.py:735`) |
| `<stem>_EXTRACTED_MERGED_NOTFLAT.fits` (standard stars) | same writer, `notFlattened=True` | `SCI_SLIT_FLUX_<ARM>` | as above plus `FLUX_DENSITY_COUNTS` (e-/nm, `horne_extraction.py:2528-2537`) | `FLUX_DENSITY_COUNTS` 3 decimals |
| `<stem>_EXTRACTED_MERGED.fits` (nod and offset stack) | `soxspipe/recipes/soxs_nod.py:1966-2133` (`stack_extractions`) | `SCI_SLIT_FLUX_<ARM>` (`soxs_nod.py:2042`) | `WAVE`, `FLUX_COUNTS`, `VARIANCE`, `SNR` (+ `FLUX_DENSITY_COUNTS` for `_NOTFLAT`) | `WAVE` 2, `FLUX_COUNTS` 3, `SNR` 2 (`soxs_nod.py:2063-2074`) |
| `<stem>_EXTRACTED_MERGED<seq>.fits` (one nod position, one cycle) | `horne_extraction.py:720` | `SCI_SLIT_FLUX_<ARM>` | as stare, no `SKY_COUNTS` | as stare; not registered as a product (`horne_extraction.py:754-755`) |
| `<stem>_FLUXCAL.fits` | `soxspipe/commonutils/flux_calibration.py:340-371` | `SCI_SLIT_FLUX_<ARM>` (`flux_calibration.py:331`) | `WAVE`, `FLUX_CALIBRATED`, `ERR_CALIBRATED`, `FLUX_COUNTS`, `VARIANCE`, `SKY_COUNTS` (when present), `SNR` (`flux_calibration.py:58-77`) | none |

The merge grid is built at `horne_extraction.py:2390-2394` (0.06 nm NIR, 0.02 nm UV-VIS) and `horne_extraction.py:2325` (`_resample_grid`), and the interpolation is `np.interp` at `horne_extraction.py:2491-2521` and `horne_extraction.py:2622` (flux density).

## Response, efficiency and SNR: the DY-2326 claim

DY-2326 decision 10 says: "No code change to the response, efficiency or SNR code." The evidence below confirms the logic and refutes three details.

### What confirms it

- The response function takes the paths of the two standard-star files as arguments (`stdExtractionPath`, `stdNotFlatExtractionPath`, `soxspipe/commonutils/response_function.py:134-151`). The callers pass whatever the extractor or the stack wrote: `soxs_stare.py:1637-1641`, `soxs_nod.py:1269-1273`, `soxs_offset.py:690-694`. A file rename touches the callers, not `response_function.py`.
- The response reads only `WAVE` (`response_function.py:263`), `FLUX_COUNTS` (through `counts_per_unbinned_pixel`, `response_function.py:265`, `soxspipe/commonutils/response_basis.py:67-80`) and `FLUX_DENSITY_COUNTS` (`response_function.py:267`). It sorts by `WAVE` (`response_function.py:183`) and interpolates the catalogue flux at each `WAVE` (`response_function.py:285-290`). Nothing assumes equal bins.
- The spline fit works in nm: knot spacing, telluric regions and line half-widths are all in nm (`soxspipe/commonutils/response_fit.py:52-56`, `response_fit.py:256-275`).
- The efficiency's Savitzky-Golay filter is 21 rows wide (`response_function.py:364`). On a log-λ grid 21 rows is a constant velocity window, so the "constant-velocity window" of DY-2326 needs no code.
- The SNR and efficiency QCs (`add_snr_efficiency_qcs`, `soxspipe/commonutils/toolkit.py:3468-3630`) take medians per order, with orders assigned by the order-join wavelengths. They work with any grid.
- Flux calibration interpolates the response in log space at each `WAVE` (`flux_calibration.py:80-104`). It works with any grid.

### What refutes it, in part

1. **The standard star's file must keep `FLUX_COUNTS`.** The response reads the counts from the file on disk, not from memory. DY-2293 says the internal product holds "the same columns" as the `_ph3` file, and that an uncalibrated `_ph3` file puts the counts in `FLUX` and omits the `*_COUNTS` columns. A standard-star run is always uncalibrated. If the internal file follows that rule, `counts_per_unbinned_pixel` raises `KeyError: 'FLUX_COUNTS'`. Decide that the internal merged spectrum always carries `FLUX_COUNTS`, calibrated or not, and that the uncalibrated layout rule applies only to `_ph3`.
2. **The efficiency needs `FLUX_DENSITY_COUNTS` from the not-flattened file.** That column is not in the DY-2293 column set, and the `_NOTFLAT` file is not named in the "replaces both files" rule. Decide whether `_EXTRACTED_MERGED_NOTFLAT.fits` stays as an internal standard-star file (simplest), or becomes `_SPECTRUM_NOTFLAT.fits` with `FLUX_DENSITY_COUNTS` kept.
3. **The response basis changes meaning.** Today `np.interp` samples the flat-fielded counts of the nearest native pixels, so a `FLUX_COUNTS` value is "flat-fielded electrons per native pixel" (the basis that `response_basis.py:1-25` documents). DY-2326 decision 3 makes each output value "the overlap-weighted sum of input pixels". An output bin narrower than a native pixel then holds only part of that pixel's counts, so the value becomes "electrons per output bin", which is the per-pixel value times (bin width / pixel width). The grid step is the smallest per-order median native step, so that ratio is at or below 1 and changes from order to order. Two consequences follow.
   - A RESP product made before the change, applied to a science spectrum made after it, is wrong by that ratio. The response code needs no logic change, but `RESPONSE_BASIS = "PER_UNBINNED_PIXEL"` (`response_basis.py:22`) must take a new value, so that `require_response_basis` (`response_basis.py:83-101`) refuses the old products, as it did for DY-1580. The tests that pin the value must follow.
   - The flat field has already removed the pixel scale from the flat-fielded counts (`response_basis.py:4-9`), so they are a density, not a count. Summing them with overlap fractions puts the pixel-scale ratio back in, with steps at the order joins where the native pixel width jumps. A science spectrum and its response on the same grid cancel the ratio, but the response spline (knots in nm) may not follow a step at a join. Weighted means (weights normalized to sum to 1) avoid that; the variance then takes the normalized weights squared. This is a design choice for the DY-2326 build issue, not a finding about today's code.
4. **`FLUX_DENSITY_COUNTS` is a density in e-/nm** (`horne_extraction.py:2616-2622`: flux divided by `pixelScaleNm`). DY-2326 says it "uses the same weights" as the flux. If those weights are summed, not normalized, the efficiency (`response_function.py:342-350`) is wrong by the bin-to-pixel ratio. Resample it with normalized weights, or compute it after the merge as the not-flattened counts divided by the output bin width in nm.

The SNR QC thresholds `snr-median: [50,1000]` and `[7,200]` (`soxspipe/soxs_default_settings.yaml:860` and `:905`, both under `soxs-nod-std`) need no change in code. The per-bin SNR changes a little, because bins differ from today's grid; the real-data A/B in the build issue should check them.

## Consumers that must change

The table below lists each place, what it reads, and the change. "Writer" rows are where the suffix or the column names are made; the readers depend on them.

| Place | What it reads or writes | Change needed |
|---|---|---|
| `soxspipe/commonutils/horne_extraction.py:720-768` | writes `_EXTRACTED_MERGED{seq}{notflat}.fits`, product label `EXTRACTED_MERGED_TABLE` | Write `_SPECTRUM` (or hand the merged frame to the recipe, which writes one file after flux calibration). Keep a file on disk before `response_function` runs on a standard star. Decide the per-nod-position and `_NOTFLAT` names. |
| `soxspipe/commonutils/horne_extraction.py:735` | rounds `WAVE` to 2, counts to 3, `SNR` to 2 decimals on write | Remove (DY-2293: no rounding; DY-2326 decision 7: `WAVE` at full float64). The response reads this rounded file today. |
| `soxspipe/commonutils/horne_extraction.py:2390-2537` (`merge_extracted_orders`), `:2325` (`_resample_grid`), `:2622` | linear grid and `np.interp` | Log-λ grid, flux-conserving weights, hard cut at joins, `ERR`/`QUAL` per bin (DY-2326, DY-2293). This is the producer, listed so the column names stay `WAVE`, `FLUX_COUNTS`, `VARIANCE`, `SKY_COUNTS`, `SNR`, `FLUX_DENSITY_COUNTS` for the in-memory readers below. X-Shooter keeps the old path. |
| `soxspipe/commonutils/flux_calibration.py:17-18` | requires `WAVE`, `FLUX_COUNTS`, `VARIANCE`, `SNR` in memory | No change if the in-memory merge frame keeps those names. If `VARIANCE` becomes `ERR_COUNTS`, change the list and the error formula at `flux_calibration.py:296-301`. |
| `soxspipe/commonutils/flux_calibration.py:58-77`, `:331-371` | writes `_FLUXCAL.fits` with `FLUX_CALIBRATED`, `ERR_CALIBRATED`; label `EXTRACTED_FLUXCAL_SPECTRUM` | Return the calibrated columns instead of writing a file, or write `_SPECTRUM` with `FLUX`, `ERR`, `QUAL`, `BGFLUX`, `FLUX_COUNTS`, `ERR_COUNTS`, `SNR`, `SKY_COUNTS`. One file per pass, so the merge writer and this writer must be combined. |
| `soxspipe/recipes/soxs_nod.py:1966-2133` (`stack_extractions`), `:40-80` (`median_stack_variance`) | concatenates in-memory spectra; groups by `round(WAVE, 4)` (`soxs_nod.py:2023-2024`); rounds on write (`:2063-2074`); writes `_EXTRACTED_MERGED{_NOTFLAT}.fits`, label `EXTRACTED_MERGED_TABLE` (`:2091-2110`); `_NOTFLAT` ASCII (`:2117-2131`) | Group by integer bin index (DY-2326 decision 7), in both functions; drop the rounding; carry `QUAL` (bit 8 for short stacks); write `_SPECTRUM`. `soxs_offset` calls it at `soxs_offset.py:598`, `:655`, `:676`. |
| `soxspipe/recipes/soxs_stare.py:1688-1711`, `soxs_nod.py:1373-1397`, `soxs_offset.py:775-799` | re-read `_FLUXCAL.fits` from disk and copy `FLUX_CALIBRATED` into `FLUX_COUNTS` for the calibrated QC plot | Read the `FLUX` column of the merged spectrum (or the in-memory frame). |
| `soxspipe/recipes/soxs_stare.py:859-883`, `soxs_nod.py:587-603` (also used by `soxs_offset.py:231`, `:261`) | builds the ASCII path from the merged path; `soxs_nod.py:600` hard-codes `_EXTRACTED_MERGED.fits` | Pass the one `_SPECTRUM` path. |
| `soxspipe/commonutils/spectrum_ascii.py:19-40`, `:63-103` | chooses `_FLUXCAL.txt` or `_EXTRACTED_MERGED.txt`; label `EXTRACTED_MERGED_ASCII`; writes `WAVE` with 2 decimals in Angstrom (`:59`) | One source file, `<stem>_SPECTRUM.txt`. Decide the label. Two decimals in Angstrom is 0.01 Å, which is finer than any log-λ bin; it can stay. |
| `soxspipe/commonutils/toolkit.py:2102-2105` (`predict_product_path`) | predicts `<stem>_EXTRACTED_MERGED.fits` for `_STARE` and `_NOD` SOFs | Predict `_SPECTRUM.fits`. The prediction drives skip-if-exists (`base_recipe.py:296-323`, `cl_utils.py:231`, `stare_mode.py:107`) and names the recipe log `<stem>_EXTRACTED_MERGED.log` and `_ERROR.log` (`toolkit.py:2184-2185`). `_OFFSET` SOFs match no branch and get `<stem>.fits`; fix that at the same time. |
| `soxspipe/resources/soxs_sof_map.yaml:131-134`, `:161-164`, `:188-191`, `:222-225`, `:244-247`, `:273-276`, `:298-301`, `:327-330`, `:351-354` | predicted product `OBJECT_TAB` with `to: _EXTRACTED_MERGED.fits` for all nine stare, nod and offset filters | Change the suffix to `_SPECTRUM.fits`. Decide the predicted `eso pro catg` (see below). |
| `soxspipe/commonutils/data_organiser.py:4088-4095` | builds `filepath` as `./reduced/<night>/soxs-<recipe with _ replaced by ->/<file>` | Strip `_obj` as `utility_setup` does (`toolkit.py:2662`), so `stare_obj` gives `soxs-stare`, not `soxs-stare-obj`. Needed before DY-2328 reads these rows. |
| `soxspipe/commonutils/toolkit.py:2946-2968` | QC plot names `_EXTRACTED_MERGED_QC_PLOT` and `_EXTRACTED_MERGED_FLUXCALIBRATED_QC_PLOT` | Optional rename; nothing reads them. |
| `soxspipe/commonutils/response_basis.py:20-23` | `RESPONSE_BASIS = "PER_UNBINNED_PIXEL"` | New value when the merge sums overlap weights (see above). |

## Consumers that need no change

- `soxspipe/commonutils/response_function.py` logic: reads paths it is given; columns `WAVE`, `FLUX_COUNTS`, `FLUX_DENSITY_COUNTS`; header cards behind the keyword aliases `SEQ_ARM`, `DATE_OBS`, `EXPTIME` and `WIN_BINY` (`response_function.py:171-176`, `response_basis.py:47`). Holds only if the new file keeps those columns and the raw-frame header cards.
- `soxspipe/commonutils/toolkit.py:3468-3630` (`add_snr_efficiency_qcs`): uses `WAVE`, `SNR` or `EFFICIENCY` and the order joins. Called in memory from `horne_extraction.py:647`, `soxs_nod.py:2080`, `response_function.py:414`.
- `soxspipe/commonutils/toolkit.py:2835-3260` (`plot_merged_spectrum_qc` and its panels): reads `WAVE`, `FLUX_COUNTS`, `SNR`, `SKY_COUNTS` from the frame it is given (`toolkit.py:3023-3059`, `:3102-3110`, `:3222-3231`). Only its callers' re-read of `_FLUXCAL` changes.
- `soxspipe/commonutils/reducer.py`: no suffix, column or tag. It reads `product_frames.sof` and `complete` (`reducer.py:357`, `:436-439`) and the run status. It is affected only through `predict_product_path`.
- `soxspipe/commonutils/calibration_matching.py:142-147`: the only science-chain rule (`std_flux`) matches `RESP_TAB%` products of the `*_std` recipes. No rule reads `OBJECT_TAB`. `tests/unit/test_calibration_matching.py:430` pins that a `SCI_SLIT_FLUX_VIS` product is not matched.
- `soxspipe/commonutils/phase3.py`: writes headers and tables, reads no spectrum.
- `soxspipe/resources/xsh_keywords.yaml:241`: `FLUXCAL` is the Phase 3 header card, not the product suffix.
- `soxspipe/commonutils/detect_continuum.py:637`: writes `OBJECT_TAB_<ARM>` on the object-trace table (`_OBJTRACE.fits`, `detect_continuum.py:596`). It is not the merged spectrum. Leave it.

## The data organiser and `product_frames`

### What it registers for science and standard products

`predict_product_frames` (`soxspipe/commonutils/data_organiser.py:4035-4097`) writes one row per entry in the `products:` list of each `soxs_sof_map.yaml` filter. It copies the raw-group columns, sets `eso pro type`, `eso pro tech` and `eso pro catg + "_" + <ARM>` from the yaml, sets `file` from the SOF name with the yaml `replace` rule, and builds `filepath`.

| Recipe | Rows per SOF | `eso pro catg` | `file` | `filepath` directory |
|---|---|---|---|---|
| `stare_obj` | 1 | `OBJECT_TAB_<ARM>` | `<SOF stem>_EXTRACTED_MERGED.fits` | `soxs-stare-obj` (wrong; real `soxs-stare`) |
| `nod_obj` | 1 | `OBJECT_TAB_<ARM>` | `<SOF stem>_EXTRACTED_MERGED.fits` | `soxs-nod-obj` (wrong; real `soxs-nod`) |
| `offset_obj` | 1 | `OBJECT_TAB_<ARM>` | `<SOF stem>_EXTRACTED_MERGED.fits` | `soxs-offset-obj` (wrong; real `soxs-offset`) |
| `stare_std` (flux) | 2 | `OBJECT_TAB_<ARM>`, `RESP_TAB_<ARM>` | `_EXTRACTED_MERGED.fits`, `_RESP.fits` | `soxs-stare-std` (right) |
| `nod_std`, `offset_std` (flux) | 2 | as above | as above | `soxs-nod-std`, `soxs-offset-std` (right) |
| telluric `*_std` | 1 | `OBJECT_TAB_<ARM>` | `_EXTRACTED_MERGED.fits` | right |

I checked this against the definitive workspace database (`/Volumes/SOXS/raw_archive/soxs_commissioning_data/soxspipe.db`, opened read-only). For example, `stare_obj | OBJECT_TAB_NIR | ./reduced/2026-01-24/soxs-stare-obj/20260125T034741_NIR_3_STARE_OBJ_SLIT0_5_5_0S_SOXS_CD-325599_EXTRACTED_MERGED.fits`. Under `/Volumes/SOXS/reductions/`, 129 product folders are named `soxs-stare` and 36 `soxs-stare-std`, and none is named `soxs-stare-obj`.

The `product_frames` columns are: `uuid`, `instrume`, `eso seq arm`, `file`, `sof`, `mjd-obs`, `mjd-date`, `night start date`, `night start mjd`, `eso pro catg`, `eso pro type`, `eso pro tech`, `exptime`, `binning`, `rospeed`, `slit`, `slitmask`, `lamp`, `recipe`, `recipe_order`, `template`, `filter`, `object`, `eso tpl name`, `eso tpl start`, `eso tpl nexp`, `eso obs name`, `eso obs start`, `eso obs id`, `set_first_file`, `complete`, `status`, `status_<session>`, `filepath`, `gain`, `error_message`, temperature and AFC columns, `absrot`. DY-2328 can match arms on `eso obs id`, `eso tpl start` and `object` from these.

### Where the `filepath` is built, and why it is wrong

`soxspipe/commonutils/data_organiser.py:4088-4095` concatenates `"./reduced/" + night + "/soxs-" + recipe.replace("_", "-") + "/" + file`. The recipe value is the yaml `recipe:` key (`stare_obj`, `nod_obj`, `offset_obj`). The recipes write to `utility_setup`'s directory, which drops `-obj` (`soxspipe/commonutils/toolkit.py:2662`, pinned by `tests/unit/test_toolkit.py:536-548`). The organiser does not. The fix is one `.replace("_obj", "")` before the hyphen swap, or one shared helper that both call.

The row is also wrong in three other ways that DY-2328 must handle:

- `eso pro catg` is `OBJECT_TAB_<ARM>`, but the file carries `SCI_SLIT_FLUX_<ARM>`, and `OBJECT_TAB_<ARM>` is the real catg of the `_OBJTRACE.fits` table. A SOF built from `product_frames` would tag the merged spectrum with the trace-table tag. Choose one catg for the merged spectrum (for example keep `SCI_SLIT_FLUX_<ARM>`, or a new one) and use it in both the yaml and the writer.
- The stare-mode reduction of a nod SOF (`--stare`, `reduce-as-stare`) writes `<SOF stem>_STAREMODE_...` under `soxs-stare` (`soxspipe/commonutils/stare_mode.py:11`, `docs/source/soxs_file_naming_scheme.md:61`), but its row still names the nod file under `soxs-nod-obj` (`tests/integration/test_data_organiser_nod_stare.py:222`).
- A SOF reduced one frame at a time writes `<sof>_F<N>_...` and never the predicted file (`tests/integration/test_reduce_all_skips_per_frame_stare.py:1-5`).

No code reads the science `OBJECT_TAB` rows or their `filepath` today. The `filepath` join at `data_organiser.py:3425-3434` (failed-product cascade) compares `product_frames.filepath` with `sof_map.filepath`; science products never enter a `sof_map`, so the wrong directory has no effect yet.

## Tests and fixtures

These must change with the rename or the new columns. Counts are matching lines.

- Suffix in product names or predicted paths: `tests/unit/test_toolkit.py:243`, `:264`, `:294`, `:360`, `:392-453` (`predict_product_path`); `tests/unit/test_base_recipe_output_directory.py:20`, `:97`; `tests/unit/test_cli_dispatch.py:1463`, `:1659`; `tests/integration/test_data_organiser_nod_stare.py:222`; `tests/integration/test_reduce_all_skips_per_frame_stare.py:3` (docstring).
- `_FLUXCAL` file and `FLUX_CALIBRATED` / `ERR_CALIBRATED` columns: `tests/unit/test_flux_calibration_columns.py` (15 lines, 111-298); `tests/integration/test_flux_calibration_contracts.py:105-352` (label `EXTRACTED_FLUXCAL_SPECTRUM` at `:111`); `tests/integration/test_response_flux_basis.py:135`.
- ASCII naming and label: `tests/unit/test_spectrum_ascii.py` (22 lines, 21-171).
- Product labels and merged-file names: `tests/unit/test_horne_extraction_helpers.py:897-1033`; `tests/integration/test_soxs_nod_cycle_characterization.py:638-733`; `tests/integration/test_observing_recipe_orchestration.py:184-1789`; `tests/integration/test_soxs_stare_reduction_characterization.py:76`, `:827`, `:873`, `:1253`, `:1285`; `tests/integration/test_soxs_nod_reduction_characterization.py:741-825`; `tests/integration/test_soxs_offset_reduction_characterization.py:436-954`; `tests/integration/test_stare_no_sky_model_chain.py:73-110`; `tests/unit/test_toolkit_characterization.py:590-591`; `tests/unit/test_merged_spectrum_qc_characterization.py:407`.
- `OBJECT_TAB` as the predicted catg: `tests/integration/test_calibration_matches_table.py:46-48`; `tests/integration/test_data_organiser_nod_stare.py:222`.
- Grid, interpolation and stacking: `tests/unit/test_merge_variance_mask.py`, `tests/unit/test_nod_stack_variance.py` (groups by rounded `WAVE`, `:57-109`), `tests/unit/test_numerical_defect_pins.py:221`, `tests/integration/test_soxs_nod_shift_and_add.py`, `tests/integration/test_soxs_nod_host_contamination.py`, `tests/integration/test_soxs_nod_no_trace.py`, `tests/integration/test_response_basis_early_check.py`.
- Fixture: `tests/factories/dataframes.py:252-273` (`merged_spectrum`) builds the in-memory merge frame with `WAVE`, `FLUX_COUNTS`, `VARIANCE`, `SKY_COUNTS`, `SNR`. Add the new columns there.
- Real data (opt-in): `tests/real_data/test_stare_baseline.py:183-272` and `tests/real_data/test_nir_offset_baseline.py:52-96` assert the file names, the `_EXTRACTED_MERGED.log` name and the `_FLUXCAL` column list, and hold median baselines that the new grid will move.
- Response basis: `tests/integration/test_response_flux_basis.py:167` pins `PER_UNBINNED_PIXEL`; it changes with the value.

## Docs

- `docs/source/recipes/output/soxs_stare.inc:8-24`, `soxs_stare_std.inc:8-26`, `soxs_nod.inc:8-19`, `soxs_nod_std.inc:8-21`, `soxs_offset.inc:8-19`: product tables list `EXTRACTED MERGED TABLE`, `EXTRACTED FLUXCAL SPECTRUM` with the old column list, the ASCII names and the QC plot names.
- `docs/source/soxs_file_naming_scheme.md:61-63`: `<SOF>_EXTRACTED_MERGED.fits` and the stare-mode and per-frame variants.
- `docs/source/recipes/soxs_nod.md:97`: `<sof name>_STAREMODE_EXTRACTED_MERGED.fits`.
- `docs/source/data_organiser_and_reducer/data_organiser.md:27`, `:56`: predicted products and `product_frames`; add the `filepath` fix if it lands.
- `docs/source/apidocs/soxspipe/soxspipe.commonutils.toolkit.md`: generated; it regenerates from the docstrings of `predict_product_path` and `plot_merged_spectrum_qc`.
- `CHANGES.md`: three historic entries mention the old names; leave them.

## Scripts

- The repository has no script that reads these files (`tools/` holds only CI and lint helpers; no match).
- `../P32.14_soxspipe_general/output/scripts/compare_products.py` (workspace, not repo) compares every `*.fits` under `reduced/` by relative path and by column names (`compare_products.py:79-87`, `:109-114`). It does not name any suffix, so it needs no change. An A/B across the rename reports the old and new files as `ONLY IN A` and `ONLY IN B`, and the same file with new columns as `DIFFERS`; the build's A/B must compare the old `_FLUXCAL.fits` with the new `_SPECTRUM.fits` by hand or with a rename map. `drive_recipe.py` and `stage_sof.py` print the returned product path and name no suffix.
- `output/scripts/DY-2297/lamrms.py` and `linewidth.py` in this workspace read arc and line tables, not spectra.

## Unverified

- The size of the step at order joins that summed overlap weights would put in the flat-fielded counts. I reasoned it from the code; no reduction was run.
- Whether `ESO TPL START` is filled in `product_frames` for every science row (the column exists in the definitive database; I did not count nulls).
- That nothing outside the repository (for example ESO tooling or Dave's notebooks) reads `_FLUXCAL.fits` by name.

## Sources

- Repository `soxspipe`, branch `develop`, commit `8b039737` (2026-11-20). Every `file:line` above is from that commit.
- Linear DY-2293 resolution comment (2026-11-19), DY-2326 resolution comment (2026-11-20), and ADR 0003 in the ESO Phase 3 workspace (`docs/adr/0003-fixed-log-wavelength-grid-per-arm.md`).
- `/Volumes/SOXS/raw_archive/soxs_commissioning_data/soxspipe.db`, `product_frames` table, read-only query on 2026-11-20; folder listing of `/Volumes/SOXS/reductions/*/sessions/*/reduced/*/`.
