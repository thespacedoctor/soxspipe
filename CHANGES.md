# Release Notes


* **REFACTOR**: `soxspipe/recipes/base_recipe.py`, part 3 of 4 of its tier-2 conversion under the house-rule refactor (DY-87, following DY-85 and DY-86), splits the largest functions left between 50 and 200 lines, leaving docstrings, comments, lint and the template revision for DY-88. `clip_and_stack` (198 lines) gains `_ccds_for_stacking`, `_clip_individual_frames` and `_mean_combine_frames`; `detrend` (155) gains `_subtract_dark_frame` and `_subtract_scattered_light`; `qc_ron` (131) gains `_measure_raw_frame_ron` and the `_image_stats` static method that replaces its nested closure; `update_fits_keywords` (114) gains `_stamp_raw_frame_records`, `_stamp_calibration_frame_records` and `_stamp_recipe_parameter_records`; `report_output` (107) gains `_qc_report_columns` and `_format_qc_for_display`. The largest function in the module is now `clip_and_stack` at 138 lines. `prepare_frames`, `clean_up`, `_write`, `flag_poor_data`, `qc_median_flux_level`, `add_qc` and `add_product` are left unsplit deliberately. New characterization tests in `tests/unit/test_base_recipe_output_characterization.py`, 8 tests, land before the split and cover the module's largest uncovered block, `qc_ron`'s raw read-out-noise measurement, plus `update_fits_keywords` and `_write`. Module coverage goes from 88.01% to 92.4% line and branch, the module's whole-file ruff count goes from 64 to 61, and the offline suite reaches 983 tests. Two defects found while characterizing are pinned as they stand and filed, not fixed: DY-91, `_write` writes 1, not 0, into masked pixels when `maskToZero=True`, and DY-92, `qc_ron` records no MASTER RON QC row when it measures the master frame itself. A worktree A/B run against `develop` shows the stacked frame, its mask and uncertainty, the detrended frame and both read-out-noise values bit-identical.
* **REFACTOR**: `soxspipe/recipes/base_recipe.py`, part 2 of 4 of its tier-2 conversion under the house-rule refactor (DY-86, following DY-85), splits the three functions over 200 lines, leaving the remaining splits for DY-87 and docstrings and lint for DY-88. `__init__` (205 lines) gains `_resolve_product_path`, `_open_session_database`, `_advanced_settings` and `_empty_qc_and_product_tables`; `_prepare_single_frame` (225) gains `_add_uncertainty_map`, `_bad_pixel_mask`, `_clean_cosmic_rays` and `_write_prepared_frame`; `_verify_input_frames_basics` (232) gains `_verify_single_instrument_and_arm`, `_verify_single_binning`, `_verify_single_read_speed`, `_verify_single_gain`, `_verify_single_slit_width`, `_verify_single_readnoise`, `_collect_frame_classifications` and a shared `_report_verification_error` that replaces the seven identical copies of the error banner each guard printed before raising. The module now has no function over 200 lines and ruff reports `PLR0915` against none of the three. New characterization tests in `tests/unit/test_base_recipe_init_characterization.py`, `tests/unit/test_base_recipe_prepare_frame_characterization.py` and `tests/unit/test_base_recipe_verify_frames_characterization.py`, 36 tests, land before the split and reach the branches the offline suite never did: the constructor's standard-star renaming, both refusals to overwrite an existing product and the session-database status bookkeeping; frame preparation's corrupted-frame and binary-table early returns, the dummy bad-pixel map written to a missing path, cosmic-ray cleaning, the response-curve flag and the SOXS temperature keywords, with the gain conversion and both uncertainty maps pinned numerically; and every "input frames are a mix of X" guard plus the NIR paths that read gain and readnoise from the detector parameter file and the exemption that lets a SOXS NIR nodding recipe carry two slit widths. Module coverage goes from 73.97% to 88.01% line and branch and the suite from 938 to 975 tests, and all 37 new tests pass unchanged against the branch point. `_prepare_single_frame` is the module's one hot-path function, called once per input frame; an A/B benchmark of old and new on a 256x256 frame, 40 reps alternated across two rounds, gives medians of 9.45 ms and 9.01 ms for the base revision against 9.72 ms and 8.97 ms for the split, a ratio of about 1.00 and well inside the 1.2x threshold, with the prepared data, uncertainty and mask sums bit-identical. No extracted method emits an entry or exit log the base revision did not emit, no side effect is reordered (the arm values are still read before the recipe settings are resolved, which cost one argument on `_verify_single_instrument_and_arm`), and the public surface is unchanged. Three function-local imports move with the code they serve, `soxspipe.commonutils.toolkit` into `_resolve_product_path`, `sqlite3 as sql` into `_open_session_database` and `yaml` into `_advanced_settings`; all three are already in `sys.modules` before this module finishes importing, so no load time and no `ImportError` surfacing moves, and none is monkeypatched by module-level name. `pandas` is not already loaded, so its import stays at constructor entry and the module is passed into `_empty_qc_and_product_tables`, which code review asked for. The existing QC column-order test now reads the literal from `_empty_qc_and_product_tables`, where it moved to, and a new test asserts the constructor still assigns both tables in that order, so the wiring is pinned as well as the literal. The module's whole-file ruff count goes from 78 to 75, three `PLR0915` findings cleared by the splits themselves; ten pre-existing findings that the code motion re-dated onto changed lines are suppressed on their own lines with reasons rather than rewritten, since rewriting them would be a behaviour edit and the module's lint sweep is DY-88's, and the changed-line ratchet is clean with 64 pre-existing findings in the same files left alone. The count of functions over 50 lines rises from 15 to 20, the expected effect of each extracted method carrying a house-format docstring. Two defects are pinned as they stand and filed, not fixed: DY-89, the mixed-arm guard raises `KeyError` because it interpolates `%(imageTypes)s` before that name exists, and DY-90, the constructor cannot be given anything but a `.sof` path because `startNightDate` is set on one branch only.
* **REFACTOR**: `soxspipe/recipes/base_recipe.py` is the second and last tier-2 module converted under the house-rule refactor, and this is the first of four pull requests against it (DY-85, with DY-86, DY-87 and DY-88 to follow): characterization tests plus adoption of the shared QC helpers, leaving the three functions over 200 lines, the remaining splits, docstrings and lint for later pull requests. Nine of the module's ten duplication hits are converted: the `SXSPRE` pre-tag timestamp in `_prepare_single_frame` now calls `utcnow_string(microseconds=True)`, and `flag_poor_data`, `qc_ron` and `qc_median_flux_level` build their five QC rows through the `add_qc` delegator, each function still minting one timestamp and sharing it across its rows. Contrary to the ticket's instruction, all five rows adopt the helper despite listing `qc_comment` before `qc_unit`, the reverse of the helper's key order: `__init__` creates `self.qc` with all nine columns up front and `pd.concat` aligns an appended row to the columns the table already has, so the row's key order never reaches the table, and the template's rule is missing that precondition. New characterization tests in `tests/unit/test_base_recipe_helper_adoption.py`, 12 tests, pin the `__init__` column literal, the resulting column order after every append, each row's values, the two rendered timestamp formats, and, under a clock stub that returns a new value per call, that `flag_poor_data` and `qc_ron` each read the clock exactly once; a mutation minting a second timestamp for the CPATH TEMP row fails them. The integration suite gains an assertion pinning the `SXSPRE` card's microsecond format. Coverage goes from 73.93% to 73.97% line and branch and the suite from 936 to 938 tests, the module from 2252 to 2197 lines, and duplication hits from 10 to 1 — the survivor being the `utcnow` assigned and never read in `detrend`, which is dead but is an `F841` rather than a conversion and goes to DY-88. The four now-unused function-local `from datetime import datetime` imports and three `import pandas as pd` imports are removed, no import is hoisted to module level, the module's whole-file ruff count stays at 78, no hot-path function was created, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/toolkit.py`, part 4 of 4 of its tier-2 conversion under the house-rule refactor (DY-82, following DY-77, DY-79 and DY-80), takes the remaining function splits: six functions split, chosen as the five ruff reported `PLR0915` against plus `read_spectral_format`. `unpack_order_table` gains `_order_table_axes`, `_order_table_pixel_grid`, `_order_polynomial_coefficients`, `_evaluate_order_polynomial` and `_rescale_order_tables_for_binning`, with the four near-identical inline Chebyshev evaluations collapsing onto one helper and the standard-deviation polynomial keeping its deliberate reuse of the centre-trace degrees; `read_spectral_format` gains `_dispersion_axis_pixel_limits`; `twoD_disp_map_image_to_dataframe` gains `_associated_frame_binning`, `_open_disp_map_hdus`, `_disp_map_pixel_columns` and `_add_pixel_scale_column`, and loses the `PLR0915` suppression DY-80 deferred to this ticket; `quicklook_image` gains `_quicklook_header_lookups`, `_apply_inter_order_mask` and `_draw_detector_image`; `spectroscopic_image_quality_checks` gains `_inner_order_mask`; and `plot_merged_spectrum_qc` gains `_draw_merged_flux_panels` and `_draw_merged_snr_panel`. New characterization tests in `tests/integration/test_toolkit_split_characterization.py`, 6 tests, pin the previously uncovered branches inside the functions being split — the x-axis order table with its edge polynomials, axis-A binning, the block-reduced binned dispersion-map path, the home-relative map path, and the dispersion-map pixel limits with their clamping and binning — taking module coverage from 93.99% to 97.6% line and branch and the suite from 920 to 926 tests. `plot_merged_spectrum_qc` had gone back over 200 lines at the branch point (217), because DY-80's docstrings lengthened it; it is now 142, and no function is over 200 lines. The module's whole-file ruff count goes from 4 to 0, the docstring checker stays at 0 findings, and the count of functions over 50 lines rises from 18 to 21, the expected effect of each extracted helper carrying a house-format docstring. No extracted helper emits an entry or exit log the base revision did not emit, and the public surface is unchanged. Function-local imports move into the new helpers — `math`, `numpy` and `pandas` into `_order_table_pixel_grid`; `numpy` and `pandas` into `_dispersion_axis_pixel_limits`; `astropy.io.fits` and `astropy.nddata` into `_open_disp_map_hdus`; `numpy` into `_disp_map_pixel_columns`, `_add_pixel_scale_column` and `_inner_order_mask`; `astropy.stats` into `_draw_merged_flux_panels` and `_draw_merged_snr_panel`; and `soxspipe.commonutils`'s `detector_lookup` and `keyword_lookup` into `_quicklook_header_lookups` — of which `math`, `numpy` and `soxspipe.commonutils` are rebinds, already in `sys.modules` by the time `toolkit` imports, while `pandas` and the astropy submodules are deferred rather than dropped; `_quicklook_header_lookups` deliberately keeps the lookup import function-local, because the module-level names are bound at import and would not see a lookup patched on the package, as four existing tests proved when the import was first hoisted. The code review blocked on the relocated imports, and the repo owner accepted them on the record, under template v4's rule that a relocated function-local import is stated in the pull request rather than forbidden.
* **REFACTOR**: `soxspipe/commonutils/toolkit.py`, part 3 of 3 of its tier-2 conversion under the house-rule refactor (DY-80, following part 1's DY-77 and part 2's DY-79), closes out the docstring drift and the named ruff findings, with no duplication to adopt and no functions left over 200 lines to split. New characterization tests in `tests/unit/test_toolkit_lint_characterization.py`, 35 tests, pin the behaviour the lint fixes touch: UP031's disagreement on `np.ma.masked` formatting, E712's `== False` vs `.eq(False)` equivalence, the SIM108 ternary branches, the unpack_order_table clamps, the spectroscopic mask loop, the read_spectral_format dispersionMap branch, MaxFilter.filter's `None` return, and the `KeyError` each retained F841 lookup raises. Docstrings for all 15 drifted functions are reconciled to the house format, argument bullets move to double backticks, missing arguments and Return sections are added, `plot_merged_spectrum_qc`, `get_skylines_dataframe` and `calculate_rolling_snr` are converted from Google style, `extinction_correction_factor` gains a docstring, and 13 prose comments are uppercased, taking the docstring checker to 0 findings. Eight SIM108 ternaries, `strict=False` on seven zip calls, an `l` rename in seven places, `.eq(False)`, an explicit `return None` in `MaxFilter.filter` and an inlined extinction return clear their findings, two pieces of pure dead code are deleted, and seven new noqa comments, each with a reason, cover the retained F841 lookups, the masked-constant UP031 pair, a SIM401 the SIM108 rewrite exposed, N802 and PLR0915 on the still-oversized `twoD_disp_map_image_to_dataframe`, and the deliberately disabled debug plot in `cut_image_slice`; the module's whole-file ruff count goes from 54 to 4, all PLR0915 left for DY-82, and toolkit coverage goes from 90.9% to 93.99% line and branch. Two defects are pinned and filed, not fixed: DY-83, `predict_product_path`'s duplicated elif makes the `_EXTRACTED_MERGED` branch unreachable, and DY-84, `cut_image_slice` returns the built-in `slice` type when `median=False`. No hot-path function was created, no function-local import moved, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/toolkit.py`, part 2 of 3 of its tier-2 conversion under the house-rule refactor (DY-79, following part 1's DY-77), splits the two functions over 200 lines. `quicklook_image` is split into `_quicklook_frame_array()`, `_quicklook_instrument()`, `_rotate_for_display()`, `_draw_surface_plot()`, `_draw_dispersion_grid_lines()` and a shared `_label_detector_axes()` that replaces two identical axis-label blocks, cutting it from 270 to 160 lines; `plot_merged_spectrum_qc` is split into `_set_wavelength_xlim()`, which replaces four identical try/except blocks, `_snr_orders_and_values()`, `_annotate_snr_values()`, `_draw_sky_panel()`, `_mark_order_joins()` and `_mark_skylines()`, cutting it from 284 to 181 lines; every helper is a contiguous block moved unchanged, none logs on entry or exit, and panel creation and drawing order are untouched. New characterization tests in `tests/unit/test_quicklook_image_characterization.py` and `tests/unit/test_merged_spectrum_qc_characterization.py`, 41 tests pinning real `Figure` objects and their artists, pass unchanged against `develop`'s module, and each of five trial mutations, reverted afterwards, failed at least one of them; toolkit coverage goes from 86.98% to 90.87% line and branch, with every line and branch of both functions and all 12 new helpers covered, and the suite goes from 844 to 885 tests. Lint findings re-dated onto moved lines are cleared (`F841`, `E741`, `B905`, `B007`), and one new `# noqa: E712` with a reason is added where a bare truth mask would treat a non-boolean isolated column differently; the module's whole-file ruff count goes from 59 to 54. One defect is pinned and filed, not fixed: DY-81, `quicklook_image` never applies the inter-order mask in the dispersion-map branch. Docstrings and the module's remaining lint are left for DY-80. No hot-path function was created, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/toolkit.py` is the first tier-2 module converted under the house-rule refactor, and this is the first of several pull requests against it (DY-77): characterization tests plus adoption of the module's own shared helpers, leaving the two functions over 200 lines, the remaining oversized functions, docstrings and lint for later pull requests. `generic_quality_checks`, `spectroscopic_image_quality_checks` and `add_snr_efficiency_qcs`, previously naive `datetime.utcnow()` sites, now call `utcnow_string()` alongside the already timezone-aware `plot_merged_spectrum_qc`, each function still minting one timestamp and sharing it across all its rows; the two now-unused function-local `from datetime import datetime` imports are removed. `quicklook_image` and `plot_merged_spectrum_qc` now save through `save_qc_plot()`, whose defaults (`dpi=120`, `format="pdf"`, `bbox_inches="tight"`) match both old calls exactly, leaving the third `savefig` the survey counts as the one inside `save_qc_plot` itself. The QC row builders are not converted to `append_qc`, because the inline rows put `qc_comment` before `qc_unit` and the helper reverses that order, which would change the column order. New characterization tests in `tests/unit/test_toolkit_characterization.py`, 15 tests, pin each converted function's rows by column and column order, the exact rendered timestamp under a frozen clock, one shared timestamp per call under a ticking clock, and the `savefig` path and arguments plus the `clf`/`close` ordering, and all pass unchanged against `develop`'s module; `spectroscopic_image_quality_checks`, previously essentially uncovered, is now pinned for both dispersion axes and the NIR binning fallback. Toolkit coverage goes from 80.88% to 86.98% line and branch, and the suite goes from 829 to 844 tests. One defect is pinned and filed, not fixed: DY-78, an `UnboundLocalError` for a non-NIR frame with no `WIN_BINX`. The module's whole-file ruff count stays at 59, no hot-path function was created, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/set_of_files.py`, the ninth and last tier-1 module converted under the house-rule refactor, splits its 254-line `get()`, tier 1's only function over 200 lines, by input shape into `_collection_from_frames()`, `_frames_from_sof_file()`, `_common_location()` and `_supplementary_files_in_directory()`, cutting `get()` to 74 lines and collapsing the three repeated copies of the ext>0 primary-header fallback into one; the helpers emit no entry or exit logs, and the function-local `import codecs` moves into `_frames_from_sof_file()`. Characterization tests in `tests/unit/test_set_of_files.py` and `tests/unit/test_set_of_files_headers.py`, 43 new tests, take coverage from 84.4% to 100% line and branch, pinning input order across all three shapes, the ext>0 fallback join, session path rewriting, CRLF SOF lines, and the constructor defaults, and all pass unchanged against `develop`'s module. The module's whole-file ruff count goes from 15 to 0 and the docstring checker from 5 to 0; `B006` and `S101` stay under `noqa` with reasons. The only behaviour change is that the SOF read now closes its handle when `read()` raises. Four defects are pinned and filed, not fixed: DY-73, DY-74, DY-75 and DY-76. No hot-path function was created and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/polynomials.py`, the eighth tier-1 module converted under the house-rule refactor, is unchanged code, confirmed by an AST comparison with docstrings stripped matching `develop`. New characterization tests in `tests/unit/test_polynomials.py` pin all three `poly` methods at full degree with non-integer inputs at `rel=1e-12, abs=0`, the precomputed power-column names, how `exponentsIncluded` of `None`, `0` and `np.False_` branch, the wrong-coefficient-count `ValueError`, and the log levels, taking coverage from 85.1% to 100%. The two `self.exponentsIncluded == False` tests stay under `# noqa: E712` with a reason, since `not` would send `None` down the other branch. Docstrings gain the `coeff` varargs documented by bare name across all three `poly` methods, argument order matching the signature on `chebyshev_xy_polynomial`, corrected usage examples, wrapped long lines and a removed dead `line_profiler` import, taking the module's whole-file ruff count from 7 to 0. The 66-line `chebyshev_order_wavelength_polynomials.poly` is left unsplit as the optimiser's hot path, so no hot-path function was created, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/filenamer.py`, the seventh tier-1 module converted under the house-rule refactor, splits its 191-line `filenamer()` into `_binning_fragment()`, `_readout_fragment()`, `_frame_type()`, `_lamp_fragment()` and `_mask_slit()`, each a contiguous block moved unchanged. Header reads keep their order, so the early return for processed products, both `TypeError` paths and the `LookupError` for an unknown arm behave as before, and the helpers emit no entry or exit logs. Characterization tests in `tests/unit/test_filenamer.py` take coverage from 86.0% to 100% and pin the WCS binning (which truncates towards zero), the readout fragment for SOXS and X-Shooter, further frame-type and mask combinations, the diagnostic prints before each `TypeError`, and the two lookup arguments. They also pin that a missing `ESO OBS ID` raises `KeyError` and an unknown arm raises `LookupError`, so the unused `obid` and `dp` locals stay under `# noqa: F841` with a reason; the dead `obsmode` local, dead commented-out binning code and a docstring typo are removed, and `SIM108` and `SIM102` are cleared. The module's whole-file ruff count goes from 11 to 0, the public surface is unchanged, and no hot-path function was created.
* **REFACTOR**: `soxspipe/commonutils/uncompress.py`, the sixth tier-1 module converted under the house-rule refactor, splits its 81-line `uncompress()` into `_build_batches()`, which lists the `.Z` archives and groups them 25 to a batch, and `_run_batch()`, which runs one batch and reports progress. Control flow is unchanged, including the return after an `OSError` that lets later batches run and both `sys.exit(0)` calls on the missing-binary path; the helpers emit no entry or exit logs, and the function-local `Popen` import now runs once per batch. Existing tests already covered the missing-binary and non-zero-exit paths, so the new characterization tests in `tests/unit/test_uncompress.py` pin the gaps: batch boundaries at 25, the progress line and cursor-up overwrite, and carry-on after an `OSError` or a non-zero exit, taking coverage from 91.0% to 100%. `S603` is suppressed with a reason, since the command is a fixed argument list and never a shell string, and the module's whole-file ruff count goes from 1 to 0. No hot-path function was created and the public surface is unchanged. The `sys.exit(0)` defect is filed as DY-72 and not fixed here.
* **REFACTOR**: `soxspipe/commonutils/getpackagepath.py`, the fifth tier-1 module converted under the house-rule refactor, drops an unused local and an assign-then-return, clearing ruff F841 and RET504, and gains a `**Return:**` section on `getpackagepath()`, clearing the one docstring checker finding. The returned string is unchanged. Survey found 100% coverage, and `tests/unit/test_package_path.py` already asserts on the returned path rather than completion, so no characterization tests were warranted and steps 1, 2 and 3 of the template were skipped. The module's whole-file ruff count goes from 2 to 0, no hot-path function was created, and the public surface is unchanged.
* **DOCS**: `soxspipe/commonutils/combiner.py`, the fourth tier-1 module converted under the house-rule refactor, gets a house-format docstring with a `**Return:**` section on `average_combine()`, which had an empty one, and its three inline comments are now uppercase. Survey found zero ruff findings, no duplication and no oversized functions, so commits 2, 3 and 5 were skipped. No test called `average_combine()` directly, so new characterization tests in `tests/unit/test_combiner.py` pin the NaN-ignoring per-pixel mean, the mask on pixels blank in every frame, `NCOMBINE`, the unit and the combiner dtype at `rtol=1e-12, atol=0`. The private `_combination_setup` call and the function-local imports are untouched, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/dispersion_map_to_pixel_arrays.py`, the third tier-1 module converted under the house-rule refactor and the first with a real hot path, splits its two over-50-line functions into five small, contiguous helpers: `dispersion_map_to_pixel_arrays()` loses its polynomial-loading block to `_load_axis_polynomials()` and its filter/trim block to `_filter_and_trim_pixel_table()`; `get_cached_coeffs()` loses its cache-filename construction to `_dispersion_cache_file_path()`, its FITS cache read to `_load_cached_coefficients()`, and its ones-array fallback to `_default_coefficients()`. An A/B benchmark (old module from a worktree at the branch point vs new, 200 reps) shows both functions at or under 1x median, well inside the 1.2x threshold. `reset == False` becomes `not reset` and an unused loop index is renamed, both on lines the split had already touched; `E741` and a `B007` on untouched lines are left as pre-existing debt. Docstrings gain `**Return:**` sections and house-style double-backtick argument markup, clearing all 8 checker findings. Two review rounds caught that the split had silently changed import timing for `math`/`numpy`/`astropy.table.Table` and dropped a cache-metadata validation; both are restored exactly, with a new regression test pinning the validation and confirmed to fail without it. 100% coverage, public surface unchanged.
* **DOCS**: `soxspipe/commonutils/detector_lookup.py`, the second tier-1 module converted under the house-rule refactor, gets a `**Return:**` section on `get()`, naming the arm-scoped sub-dictionary it returns rather than the same-named but differently-scoped full-instrument attribute. Survey found zero ruff findings, no duplication and no oversized functions, so commits 2, 3 and 5 were skipped; existing tests in `tests/unit/test_lookups.py` already assert on returned detector parameters rather than completion, so no characterization tests were warranted. Coverage stays at 100%, and the public surface is unchanged.
* **REFACTOR**: `soxspipe/commonutils/keyword_lookup.py` is the first tier-1 module converted under the house-rule refactor. Characterization tests in `tests/unit/test_lookups.py` pin `get()`'s index formatting — a falsy index adds no suffix, a float truncates toward zero, the sign sits outside the zero padding, 100 gives `PROV100`, and Decimal and numpy scalars work, while str, bytes, bytearray, memoryview, list and complex raise `TypeError`, NaN raises `ValueError`, and inf raises `OverflowError` — and the instrument selection, covering `settings["instrument"]`, an explicit `instrument` argument winning over settings, the fallback to soxs, and the name being lower-cased. The index if/else collapses to one ternary, clearing ruff SIM108; UP031 is suppressed with `# noqa: UP031` and a comment explaining that no f-string or hand-written conversion reproduces `%d` exactly, after a private-helper attempt was dropped because reviewers found inputs where it diverged. Docstring typos are fixed, the class docstring's argument order now matches `__init__`, and two over-long lines are wrapped. The module's whole-file ruff count goes from 4 to 0, coverage stays at 100%, and the public surface is unchanged.
* **ENHANCEMENT**: exception sweep phase 2 removes the swallow from the handlers the classification marked as swallowed bugs, so those failures now propagate. `utility_setup` creates its QC and product directories with `os.makedirs(..., exist_ok=True)`, absorbing only the benign concurrent-creation race and surfacing a permission or disk-full failure where it happens instead of as a confusing write error later. `fit_bspline_curve_to_sky` no longer recomputes the knot-quantile expression that just failed, which had both lost the original traceback and silently retried a transient failure. A failed tuning iteration in `soxs_disp_solution` now propagates rather than leaving the tuning grid quietly recording nothing for that permutation. Each change carries a unit test that pins the propagation. The matching handlers in `soxs_order_centres` and `soxs_spatial_solution` keep their phase-1 log-and-swallow form, because their `parameterTuning` workers raise `NameError` on every call and cannot run at all; that defect is filed separately and its fix removes those two swallows.
* **ENHANCEMENT**: the real-data reduction gate now runs automatically on pull requests into `develop`, so a module pull request no longer depends on a human remembering to dispatch it. The trigger is limited to changes under `soxspipe/**`, so documentation and tooling pull requests skip it, and a job-level draft guard plus the `opened`, `ready_for_review`, `synchronize` and `reopened` event list keeps draft pull requests off the runner while still starting the first run for a pull request opened ready. The path filter also names `tests/real_data/**`, `pyproject.toml` and the workflow file itself, so a change to the baselines, the archive manifest or the pinned dependencies still runs the gate it belongs to. A pull request from a fork is skipped, because the archive URL and its hashes come from the manifest in the diff and a fork could otherwise point a 120-minute runner at any archive it liked; such a pull request is reached with `workflow_dispatch` after reading the diff. A `concurrency` group cancels a superseded run rather than leaving it to burn up to 120 minutes, keyed by pull request for a pull request and by run id otherwise, so a manual dispatch and the weekly schedule cannot cancel each other. `workflow_dispatch` and the weekly schedule are unchanged, and `develop` gains no branch protection, so the check is advisory.
* **FIXED**: the hard-rule gate tests no longer depend on `PATH`. The required-suite step runs pytest through `env` with an explicit variable list that omits `PATH`, so the environment's `bin` directory is invisible and the four tests in `tests/integration/test_lint_ratchet_gate.py` could not find `ruff`. `tools/lint_ratchet.py` now runs ruff as `sys.executable -m ruff` whenever ruff is installed beside the running interpreter, falling back to `PATH` only when it is not, which also pins the run to the pinned `ruff==0.16.5` rather than whichever copy `PATH` reaches first.
* **REFACTOR**: `soxspipe/recipes/soxs_mbias.py` is the pilot module for the house-rule refactor and is now converted end to end. Its seven duplication hits become `utcnow_string()`, `self.add_qc()` and `self.add_product()` calls, with `qc_unit=None`, `file_type` and `label` passed explicitly so no column the inline rows carried is dropped, and the two `qc_bias_structure` rows still sharing one timestamp. The frame stacking moves into `_combine_bias_frames` and the per-frame FFT ratio into `_periodic_noise_ratio`, cutting `produce_product` from 119 lines to 82 and `qc_periodic_pattern_noise` from 98 to 51 without touching the public surface. A worktree A/B against the branch point shows the QC values and rows bit-identical, and the new per-frame call costs 1.002x the old inline loop.
* **ENHANCEMENT**: `E722`, `S110` and `S112` are now hard repo-wide gates rather than ratcheted rules. The phase-1 exception sweep left no bare except and no silent swallow anywhere in `soxspipe`, `tests` or `tools`, so there is no pre-existing debt for a ratchet to tolerate, and a new one is a silent swallow rather than a style lapse. `tools/lint_ratchet.py` checks those three rules from the repository root — so a Python file outside `soxspipe`, `tests` and `tools`, such as `docs/source/conf.py`, is gated too — and fails on a finding wherever it sits, reporting them separately from the changed-line findings. The scan runs ruff `--isolated --ignore-noqa`, so neither a `# noqa` comment nor a `per-file-ignores` entry added later can silence a rule that is supposed to be absolute; every other rule stays ratcheted to the diff, `SIM105` included. The CI step and the pre-commit hook both pick this up without changing, because both already call the tool.
* **REFACTOR**: every bare `except:` and every silently swallowing handler in the `soxspipe` package is now narrowed to a specific exception type and reports what it caught. All 117 handlers were converted from the approved classification table, so the package holds no bare except and no silent `pass` body outside the tests. Control flow is untouched — a handler that swallowed still swallows, and propagation is a later phase — which an AST comparison of every changed file verifies. Handlers whose `try` body is the progress logger or `stdout` itself stay silent, as `contextlib.suppress`, because reporting through the resource that just failed can raise again. Eight `os.makedirs` and `os.symlink` handlers catch `OSError` rather than the table's `FileExistsError`, so a permission or disk-full failure keeps being swallowed as it was before.
* **ENHANCEMENT**: the real-data workflow now tees the reduction's output to a file and uploads it as a `reduction-log` run artifact, kept for 90 days and uploaded even when the run fails. The per-recipe run times the pipeline already prints are otherwise lost when the run log ages out, so this captures a coarse performance trend at no new infrastructure cost.
* **FEATURE**: added `tools/survey_module.py`, which reports the numbers a per-module refactor ticket opens with: line count, the module's coverage at the branch point, functions over 50 and over 200 lines, duplication hits, bare excepts, remaining whole-file ruff findings, and docstring drift. Output is Markdown bullets that paste into a ticket or a pull-request description unedited, and the duplication definition — naive `utcnow` sites, single-row QC and products `pd.concat` row builders, and `savefig` calls — is the set the shared helpers in `toolkit.py` replace.
* **FEATURE**: added `tools/check_docstrings.py`, a static checker that reports where a docstring no longer agrees with the signature it documents. It parses with `ast` and never imports the module under test, reads `__init__` argument documentation from the enclosing class docstring, and reports undocumented and phantom arguments, missing `**Key Arguments:**` sections, malformed argument bullets, and undocumented or phantom returns.
* **FEATURE**: added the shared QC, product, timestamp and plot helpers — `utcnow_string`, `append_qc`, `append_product` and `save_qc_plot` in `soxspipe/commonutils/toolkit.py`, plus thin `add_qc`/`add_product` delegators on `base_recipe` that add nothing but a fallback to `self.recipeName` and `self.dateObs`. No existing call site is converted yet. An `OMITTED` sentinel keeps "column not passed" distinct from "column passed as None", the reduction timestamp is always supplied by the caller so a batch of rows can share one, and `save_qc_plot(bboxInches=None)` omits the `bbox_inches` keyword entirely, so every call site the helpers will later replace can be reproduced exactly.
* **FEATURE**: added the lint ratchet, `tools/lint_ratchet.py`, and the CI step that runs it. A pull request now fails when a ruff finding lands on a line it changed, while the roughly 1,380 pre-existing findings stay ignored, mirroring the `diff-cover --fail-under=80` coverage gate. A `.pre-commit-config.yaml` hook runs the same check locally, and `ruff` is pinned so a new ruff release cannot fail a pull request that changed nothing.
* **FEATURE**: added a `[tool.ruff]` house-rule configuration to `pyproject.toml` (line length 120, camelCase-friendly naming ignores, `max-statements = 50`) and applied the 705 safe ruff fixes ahead of it. Import sorting is deliberately withheld from the three package `__init__.py` files, whose import order works around a circular import.
* **FIXED**: the real-data reduction is now bit-identical whichever runner GitHub allocates. The nondeterminism left after the stable-sort and residual-quantisation fixes was hardware, not code: six runs of one commit showed the five on AMD EPYC processors producing byte-identical products while the one on an AVX512-capable Intel Xeon diverged in 1,110 of 1,411 product digests, first at the master flat and then through everything derived from it, moving the merged spectrum by five rows and the flux-calibrated median by 0.33 per cent. `NPY_DISABLE_CPU_FEATURES` pins numpy's SIMD dispatch to the AVX2 path, which left only the response fit differing between Zen3 and Zen4 runners, because `np.polyfit` reaches LAPACK through `np.linalg.lstsq` and OpenBLAS selects its kernels independently of numpy. `OPENBLAS_CORETYPE`, `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` pin that too. Eight runs across four processor models, two of them Intel Xeons, then produced one identical set of product digests, matching the numbers the AVX2 runners had been producing all along, so no baseline was re-recorded. The numpy setting is a denylist, so an instruction set it does not name will diverge again.
* **FIXED**: the reduction is no longer order-dependent on the CPU it runs on. `create_dispersion_map` sorted its line table on `wavelength` alone, which is not a total key because wavelength repeats across orders and slit positions. Pandas defaults to an unstable quicksort and NumPy dispatches float sorts to SIMD kernels whose permutation of tied rows varies by CPU microarchitecture, so the same input produced different row orders on different runners, changing summation order in the polynomial fits. The sort now uses the full `(wavelength, order, slit_index)` key with `kind="stable"`, and all 20 `sort_values` call sites in the reduction path were made stable.
* **FIXED**: the sigma-clipping decision in `create_dispersion_map.fit_polynomials` no longer flips on last-bit floating-point noise. A residual a fraction of an ULP either side of the threshold changed which lines survived, which changed the next `curve_fit`, which moved the merged spectrum's red end by 12 bins. Residuals are now quantised to 9 decimals as they are computed, well below any physically meaningful difference and well above the noise floor.
* **FIXED**: the real-data snapshot unit tests now use a dedicated `cache_root` fixture and look the cached snapshot up by its manifest digest, instead of sharing `tmp_path/cache` with the autouse `isolated_runtime` fixture and taking the first `iterdir()` entry. The old arrangement passed locally and failed on CI, where directory order could put the matplotlib or numba cache first.
* **FIXED**: QC and product rows are now built with `pd.DataFrame([row])` instead of `pd.Series(row).to_frame().T`, so appending a row no longer upcasts the whole `qc_value` column to `object` and numeric QC values keep their dtype through to the FITS headers.
* **REFACTOR**: removed the deprecated inline test tree (`soxspipe/recipes/tests/`, `soxspipe/commonutils/tests/`, `soxspipe/tests/`), 33 files and 3,892 lines. The top-level `tests/` tree is canonical. `soxspipe/utKit.py` is retained, since `tests/unit/` still imports it.
* **TEST**: the real-data NIR offset baseline is asserted within tolerance bands rather than as exact scalars, and the band centres were re-recorded from agreeing post-fix runs. The previous approved values were one face of a coin flip and matched no code state. Bit-identical reproduction is verified on Zen 3 and Zen 4 only, so the bands also absorb cross-architecture variation that has not yet been measured.
* **TEST**: Add a verified, opt-in real-data NIR-offset acceptance workflow with immutable archive inventory checks and approved scalar baselines.
* **ENHANCEMENT**: Raise the whole-package branch-coverage threshold to 70% and add broad deterministic unit and integration coverage for pipeline utilities and recipes.
* **ENHANCEMENT**: refine skyline-based wavelength-shift calibration
* **ENHANCEMENT**: tune skyline detection, iterative matching, clipping, and VIS-order shift fallback behavior.
* **ENHANCEMENT**: update calibration data and default dispersion, spatial, and sky-subtraction thresholds.
* **ENHANCEMENT**: using numba JIT to speed up image rectification
* **REFACTOR**: Limit multiprocessing pool size for stare, nod, and offset recipes to avoid memory issues.
* **REFACTOR**: Update predicted paths for standard response products (now looks for RESP function instead of merge spectra when deciding if the data has been reduced yet).
* **REFACTOR**: correct rectified-image diagnostic orientation and add local CodeGraph metadata exclusions.
* **REFACTOR**: optimising converting dispersion map to pixel arrays
* **REFACTOR**: vectorising calculation of rectification weights
* **DOCS**: adding to doc FAQs
* **FIXED**: Detrend master-dark flat frames once, preserve slit-position data in dispersion-map images, and skip optional efficiency products when no estimate is available.
* **FIXED**: Guard merged spectrum QC sky plotting against non-positive sky counts
* **FIXED**: Handle empty and prepared-extension file inputs safely, omit NaN Phase 3 QC headers, and invoke FITS decompression without a shell.
* **FIXED**: Improve flux standard matching via 3 FITS header keywords
* **FIXED**: Restore editable installs by using the supported default setuptools-scm parser for `v`-prefixed tags.
* **FIXED**: Skip efficiency-product creation when the optional response efficiency estimate is unavailable.
* **FIXED**: fixing database to collect the correct response curves
* **FIXED**: pinning the normalisation factor using by ccdproc for flat correction
* **FIXED**: prevent dispersion-map transformations from running outside debug mode and handle invalid sky-plot statistics.
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
