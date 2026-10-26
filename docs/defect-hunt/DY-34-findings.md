# DY-34: static-analysis defect hunt (2 of 5)

This is an investigation only. No production code was changed and nothing was filed in Linear. The repository owner approves this table before any issue is filed.

## Summary

**Command** (ruff 0.15.20, run on `develop` at `4fb7e85`; no `--fix`, no import sorting):

```bash
ruff check soxspipe --select F821,F507,F811,B --output-format json > /tmp/ruff-dy34.json
```

**Environment:** Python 3.12.3 in a fresh venv, `pip install -e ".[tests]"`.

**Baseline:** `python -m pytest tests/unit tests/integration -q` gave **1689 passed** (2325 warnings, 155 s), so the baseline is green. The real-data tests were not run.

### Counts per rule before triage (122 findings)

| Rule | Count | Meaning |
|---|---|---|
| B905 | 82 | `zip()` without `strict=` |
| F821 | 12 | undefined name |
| B006 | 10 | mutable argument default |
| B007 | 9 | unused loop variable |
| F507 | 4 | `%`-format placeholder/argument mismatch |
| B904 | 2 | `raise` inside `except` without `from` |
| B018 | 1 | useless expression |
| F811 | 1 | redefinition of an unused name |
| B023 | 1 | closure binds a loop variable |

The counts have moved from the ticket's (16 `F821`, 21 `F507`). Earlier sweeps fixed most of them.

### Counts per outcome

| Outcome | Findings | Proposed issues |
|---|---|---|
| Reachable defect, **wrong science** | **0** | 0 |
| Reachable defect, crash or silent failure (tier 1) | 5 | 5 (I-1 … I-5) |
| Reachable defect, cosmetic (tier 2) | 3 | 1 (I-6) |
| Reachable, not a defect | 101 | 0 (one housekeeping caution) |
| Dormant (X-Shooter) | 0 | 0 |
| Already filed | 2 | DY-27 |
| Dead code | 11 | housekeeping list |
| **Total** | **122** | **6** |

**No wrong-science defect was found.** None of the reachable findings changes a numerical output without raising an error. Two findings came close, and the reasoning is recorded so it can be checked:

- `subtract_background.py:279` (I-3): in multi-order tables the last order silently reuses the *previous* order's `expandTop`. Both self-consistent detector geometries yield no numerical change. If orders run toward higher axis-A, every inter-order gap is negative, so every expansion clamps to 2 px, including the stale one. If orders run the other way, the "mask bottom of frame" block masks the whole region below the first order, which hides the stale width. This is tier 2 reasoning, not a test. See I-3.
- `base_recipe.py:2610` (I-5): QC rows can be lost silently. This changes pipeline *control* (a failed QC may not block downstream recipes), not a computed value.

### Already-filed items from the brief

| Item | Status in the current ruff output |
|---|---|
| Undefined `hduList` in the dispersion-map cache path (DY-27) | Present as rows 9 and 10, marked **already filed**. |
| `parameterTuning` / undefined `self` in `soxs_order_centres`, `soxs_spatial_solution` (DY-61, fixed) | **No longer reported.** Ruff emits no finding in either file under these rules, so there is no row. |
| `utKit.refresh_database()` (DY-32, dead code) | **No longer reported** under F821/F507/F811/B, so there is no row. |

### Method

Each finding was traced to an entry point by reading call sites (`grep` for callers, then the guard conditions on the path), not by re-reading the flagged line. Many `B905` rows zip together columns of one dataframe, or lists appended in the same loop. Those are recorded as "equal length by construction", with the construction named. Lines behind `if False:`, `if True: … else:`, `… and False`, or in functions with no caller are **dead**. None of the findings sits in an X-Shooter-only branch. The one `xsh` branch near a finding (`base_recipe.py:1023`) is the *SOXS* `else` arm, so there are no "dormant (X-Shooter)" rows.

**Entry-point key** used in the table:

| Key | Path |
|---|---|
| E1 | `soxspipe reduce` → `soxs_disp_solution` / `soxs_spatial_solution` → `create_dispersion_map.get()` |
| E2 | `soxs_spatial_solution` with multi-pinhole first-guess map (`minpin == 9`) → `create_dispersion_map.map_to_image()` |
| E3 | `soxs_order_centres` → `detect_continuum.get()` |
| E4 | `soxs_mflat` → `detect_order_edges.get()` |
| E5 | `soxs_stare` (`subtractSky`) → `subtract_sky.subtract()` |
| E6 | `soxs_stare` / `soxs_nod` / `soxs_offset` → `horne_extraction.extract()` (and its `image_transformer`) |
| E7 | every recipe → `base_recipe` product write / QC report |
| E8 | `soxspipe prep` → `data_organiser.prepare()` |
| E9 | `soxspipe reduce all` → `reducer` |
| E10 | `soxs_mflat` / `base_recipe.detrend` background step → `subtract_background.subtract()` |
| E11 | public Python API only (documented usage, not reachable from the CLI) |
| E12 | only with `--debug` |
| E13 | only with a non-default setting (named in the row) |
| E-std | standard-star reductions → `response_function.get()` |

## Findings table

Reachable defects come first, then already-filed rows, then reachable non-defects. Dead code is grouped in [Housekeeping](#housekeeping-dead-code-low-severity) at the end.

| # | Rule | `path:line` | Description | Status | Entry point traced from | Tier | Wrong science | Proposed issue title |
|---|---|---|---|---|---|---|---|---|
| 1 | F821 | `soxspipe/commonutils/data_organiser.py:207` | `directory` used before assignment when `rootDir` starts with `~`; the expansion would also be discarded (assigned to `directory`, not `rootDir`) | reachable | E11 `data_organiser(rootDir="~/…")` (documented Python API; the CLI `chdir`s and passes `.`) | 1 | no | **I-1** Workspace path beginning with `~` crashes the data organiser with UnboundLocalError |
| 2 | B905 | `soxspipe/commonutils/detect_continuum.py:1314` | `zip(uniqueOrders, colors)`: `colors` skips any order whose fit falls wholly off-detector (`continue` at :1253), so later orders get the wrong colour and the last order is dropped from this panel | reachable | E3 `detect_continuum.get()` → `plot_results()` (QC PDF only) | 2 | no | **I-6** Order-centre QC plot mislabels order colours and drops orders when one order's fit is off-detector |
| 3 | B905 | `soxspipe/commonutils/detect_continuum.py:1345` | `zip(uniqueOrders, colors)`: `colors` skips any order whose fit falls wholly off-detector (`continue` at :1253), so later orders get the wrong colour and the last order is dropped from this panel | reachable | E3 `detect_continuum.get()` → `plot_results()` (QC PDF only) | 2 | no | **I-6** Order-centre QC plot mislabels order colours and drops orders when one order's fit is off-detector |
| 4 | B905 | `soxspipe/commonutils/detect_continuum.py:1376` | `zip(uniqueOrders, colors)`: `colors` skips any order whose fit falls wholly off-detector (`continue` at :1253), so later orders get the wrong colour and the last order is dropped from this panel | reachable | E3 `detect_continuum.get()` → `plot_results()` (QC PDF only) | 2 | no | **I-6** Order-centre QC plot mislabels order colours and drops orders when one order's fit is off-detector |
| 5 | F821 | `soxspipe/commonutils/horne_extraction.py:179` | `filenamer` never imported, so the no-SOF filename fallback raises NameError | reachable | E11 → `soxs_stare(inputFrames=[…])` (`sofName=False`) → E6 `horne_extraction.__init__` | 1 | no | **I-2** Optimal extraction crashes with NameError when a recipe is run from a frame list instead of a SOF file |
| 6 | F821 | `soxspipe/commonutils/subtract_background.py:279` | `expandTop` read before assignment when the first order iterated is the reddest (`o == oBot`), i.e. a single-order order table | reachable | E10 `subtract_background.subtract()` → `mask_order_locations()` | 1 | no | **I-3** Background subtraction crashes with UnboundLocalError on a single-order order table |
| 7 | F821 | `soxspipe/commonutils/subtract_sky.py:1489` | `tck_previous` read before assignment when FITPACK reports a poor fit (`ier` 10/30) on the first b-spline pass | reachable | E5 `subtract_sky.subtract()` → `fit_bspline_curve_to_sky()` | 1 | no | **I-4** Sky subtraction crashes with UnboundLocalError instead of reverting when the first b-spline fit is poor |
| 8 | B904 | `soxspipe/recipes/base_recipe.py:2610` | Flagged for missing `from`, but the `raise` itself is unreachable (`keepTrying > 6` is never true inside `while keepTrying < 7`), so a seventh failed insert is silently swallowed | reachable | E7 `base_recipe` QC report → `_dataframe_to_sqlite(…, "quality_control")` | 1 | no (silent loss of QC rows, not a numerical change) | **I-5** Failed QC database inserts are silently dropped after seven retries, so QC failures may not block downstream recipes |
| 9 | F821 | `soxspipe/commonutils/create_dispersion_map.py:1818` | Undefined `hduList` in the dispersion-map cache write (inside `if False:`) | already filed (DY-27) | E1 `write_map_to_file()`; block is currently disabled | — | no | already filed: DY-27 |
| 10 | F821 | `soxspipe/commonutils/create_dispersion_map.py:1819` | Same as above (`hduList.writeto`) | already filed (DY-27) | as above | — | no | already filed: DY-27 |
| 11 | B905 | `soxspipe/commonutils/create_dispersion_map.py:522` | `zip()` without `strict=`. Equal length by construction: four parallel QC lists built together | reachable | E1 `_write_qc_metrics()` | 2 | no | — (no defect) |
| 12 | B007 | `soxspipe/commonutils/create_dispersion_map.py:1356` | Unused `index` in `iterrows()`; cosmetic | reachable | E12 `_plot_predicted_line_windows()` (`--debug` only) | 2 | no | — (no defect) |
| 13 | B905 | `soxspipe/commonutils/create_dispersion_map.py:1606` | `zip()` without `strict=`. Equal length by construction: four arrays derived from the same `xArray`/`yArray` | reachable | E1 `detect_pinhole_arc_lines()` | 2 | no | — (no defect) |
| 14 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2085` | `zip()` without `strict=`. Equal length by construction: two literal 4-element lists | reachable | E1 `calculate_residuals()` | 2 | no | — (no defect) |
| 15 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2170` | `zip()` without `strict=`. Equal length by construction: parallel QC lists built together | reachable | E1 `calculate_residuals()` | 2 | no | — (no defect) |
| 16 | F507 | `soxspipe/commonutils/create_dispersion_map.py:2449` | `"""…""" % locals()` with no placeholders. With a mapping on the right, `%` does not raise and the message is printed verbatim (checked: `'''curvefit x''' % {'a': 1}` → `curvefit x`) | reachable | E1 `fit_polynomials()` | 2 | no | — (no defect) |
| 17 | F507 | `soxspipe/commonutils/create_dispersion_map.py:2464` | `"""…""" % locals()` with no placeholders. With a mapping on the right, `%` does not raise and the message is printed verbatim (checked: `'''curvefit x''' % {'a': 1}` → `curvefit x`) | reachable | E1 `fit_polynomials()` | 2 | no | — (no defect) |
| 18 | F507 | `soxspipe/commonutils/create_dispersion_map.py:2478` | `"""…""" % locals()` with no placeholders. With a mapping on the right, `%` does not raise and the message is printed verbatim (checked: `'''curvefit x''' % {'a': 1}` → `curvefit x`) | reachable | E1 `fit_polynomials()` | 2 | no | — (no defect) |
| 19 | F507 | `soxspipe/commonutils/create_dispersion_map.py:2501` | `"""…""" % locals()` with no placeholders. With a mapping on the right, `%` does not raise and the message is printed verbatim (checked: `'''curvefit x''' % {'a': 1}` → `curvefit x`) | reachable | E1 `fit_polynomials()` | 2 | no | — (no defect) |
| 20 | B018 | `soxspipe/commonutils/create_dispersion_map.py:2739` | Bare `xlen, ylen` expression; no effect | reachable | E2 `create_placeholder_images()` | 2 | no | — (no defect) |
| 21 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2772` | `zip()` without `strict=`. Equal length by construction: columns of one order's pixel table | reachable | E2 `create_placeholder_images()` | 2 | no | — (no defect) |
| 22 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2775` | `zip()` without `strict=`. Equal length by construction: columns of one order's pixel table | reachable | E2 `create_placeholder_images()` | 2 | no | — (no defect) |
| 23 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2780` | `zip()` without `strict=`. Equal length by construction: tuple-unzipped from one comprehension | reachable | E2 `create_placeholder_images()` | 2 | no | — (no defect) |
| 24 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2792` | `zip()` without `strict=`. Equal length by construction: tuple-unzipped from one comprehension | reachable | E2 `create_placeholder_images()` | 2 | no | — (no defect) |
| 25 | B905 | `soxspipe/commonutils/create_dispersion_map.py:2881` | `zip()` without `strict=`. Equal length by construction: columns of the spectral-format table | reachable | E2 `map_to_image()` | 2 | no | — (no defect) |
| 26 | B905 | `soxspipe/commonutils/create_dispersion_map.py:3164` | `zip()` without `strict=`. Equal length by construction: columns of one dataframe | reachable | E2 `convert_and_fit()` | 2 | no | — (no defect) |
| 27 | B007 | `soxspipe/commonutils/create_dispersion_map.py:3789` | Unused group `name`; cosmetic | reachable | E1 `_create_dispersion_map_qc_plot()` | 2 | no | — (no defect) |
| 28 | B905 | `soxspipe/commonutils/create_dispersion_map.py:4441` | `zip()` without `strict=`. Equal length by construction: columns of the spectral-format table | reachable | E13 `create_new_static_line_list()` (`bootstrap_dispersion_solution: True`, default False) | 2 | no | — (no defect) |
| 29 | B905 | `soxspipe/commonutils/create_dispersion_map.py:4464` | `zip()` without `strict=`. Equal length by construction: `[4]`/`[0.0]` or `range(len(x))`/`x` | reachable | E13 as above | 2 | no | — (no defect) |
| 30 | F811 | `soxspipe/commonutils/data_organiser.py:705` | Second function-local `import pandas as pd` (harmless; local re-imports are load-bearing per DY-34 brief, do not auto-remove) | reachable | CLI `soxspipe list ob` → `list_obs()` | 2 | no | — (no defect) |
| 31 | B905 | `soxspipe/commonutils/data_organiser.py:860` | `zip()` without `strict=`. Equal length by construction: two columns of `rawFrames` | reachable | E8 `_sync_raw_frames()` | 2 | no | — (no defect) |
| 32 | B905 | `soxspipe/commonutils/data_organiser.py:1255` | `zip()` without `strict=`. Equal length by construction: two literal 7-element lists | reachable | E8 `_populate_raw_frames_extra_columns()` | 2 | no | — (no defect) |
| 33 | B007 | `soxspipe/commonutils/detect_continuum.py:933` | Unused `o`; cosmetic | reachable | E3 `create_pixel_arrays()` | 2 | no | — (no defect) |
| 34 | B905 | `soxspipe/commonutils/detect_continuum.py:933` | `zip()` without `strict=`. Equal length by construction: columns of the spectral-format table | reachable | E3 `create_pixel_arrays()` | 2 | no | — (no defect) |
| 35 | B905 | `soxspipe/commonutils/detect_continuum.py:947` | `zip()` without `strict=`. Equal length by construction: format-table columns plus a list appended once per row in the loop above | reachable | E3 `create_pixel_arrays()` | 2 | no | — (no defect) |
| 36 | B905 | `soxspipe/commonutils/detect_continuum.py:1014` | `zip()` without `strict=`. Equal length by construction: two outputs of one fit | reachable | E3 `fit_1d_gaussian_to_slices()` | 2 | no | — (no defect) |
| 37 | B905 | `soxspipe/commonutils/detect_continuum.py:1159` | `zip()` without `strict=`. Equal length by construction: two columns under one mask | reachable | E3 `plot_results()` | 2 | no | — (no defect) |
| 38 | B007 | `soxspipe/commonutils/detect_continuum.py:1222` | Unused `index`; loop keeps only the last row's coefficients, but `orderPolyTable` has exactly one row (`pd.DataFrame([coeff_dict])`, :837) | reachable | E3 `plot_results()` | 2 | no | — (no defect) |
| 39 | B905 | `soxspipe/commonutils/detect_continuum.py:1244` | `zip()` without `strict=`. Equal length by construction: tuple-unzipped from one comprehension | reachable | E3 `plot_results()` | 2 | no | — (no defect) |
| 40 | B905 | `soxspipe/commonutils/detect_continuum.py:1247` | `zip()` without `strict=`. Equal length by construction: polynomial evaluated on the same `df` rows as `axisBlinelist` | reachable | E3 `plot_results()` | 2 | no | — (no defect) |
| 41 | B007 | `soxspipe/commonutils/detect_order_edges.py:714` | Unused `index`; same single-row poly-table pattern as above | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 42 | B905 | `soxspipe/commonutils/detect_order_edges.py:741` | `zip()` without `strict=`. Equal length by construction: tuple-unzip | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 43 | B905 | `soxspipe/commonutils/detect_order_edges.py:742` | `zip()` without `strict=`. Equal length by construction: polynomial evaluated on the same rows as `axisBlinelist` | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 44 | B905 | `soxspipe/commonutils/detect_order_edges.py:744` | `zip()` without `strict=`. Equal length by construction: tuple-unzip | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 45 | B905 | `soxspipe/commonutils/detect_order_edges.py:745` | `zip()` without `strict=`. Equal length by construction: as :742 | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 46 | B905 | `soxspipe/commonutils/detect_order_edges.py:754` | `zip()` without `strict=`. Equal length by construction: tuple-unzip | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 47 | B905 | `soxspipe/commonutils/detect_order_edges.py:757` | `zip()` without `strict=`. Equal length by construction: same slice of two equal-length arrays | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 48 | B905 | `soxspipe/commonutils/detect_order_edges.py:766` | `zip()` without `strict=`. Equal length by construction: tuple-unzip | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 49 | B905 | `soxspipe/commonutils/detect_order_edges.py:769` | `zip()` without `strict=`. Equal length by construction: same slice of two equal-length arrays | reachable | E4 `plot_results()` | 2 | no | — (no defect) |
| 50 | B007 | `soxspipe/commonutils/dispersion_map_to_pixel_arrays.py:48` | Unused `index`; cosmetic | reachable | E1/E6 `_read_dispersion_map_axes()` | 2 | no | — (no defect) |
| 51 | B905 | `soxspipe/commonutils/horne_extraction.py:345` | `zip()` without `strict=`. Equal length by construction: both lists come from one `image_transformer` (all built from the same `orderTraces`) | reachable | E6 `extract()` | 2 | no | — (no defect) |
| 52 | B905 | `soxspipe/commonutils/horne_extraction.py:367` | `zip()` without `strict=`. Equal length by construction: `fmultiprocess` returns one result per input, in order (`map`) | reachable | E6 `extract()` | 2 | no | — (no defect) |
| 53 | B023 | `soxspipe/commonutils/horne_extraction.py:481` | Lambda captures loop variable `decimals`, but `.apply` consumes it in the same iteration, so no late-binding effect | reachable | E6 `extract()` | 2 | no | — (no defect) |
| 54 | B905 | `soxspipe/commonutils/horne_extraction.py:709` | `zip()` without `strict=`. Equal length by construction: one shift appended per order; the only `continue` (empty order) cannot fire because `uniqueOrders` comes from the same frame | reachable | E6 `tune_wavelength_calibration_to_skylines()` | 2 | no | — (no defect) |
| 55 | B905 | `soxspipe/commonutils/horne_extraction.py:714` | `zip()` without `strict=`. Equal length by construction: as :709 | reachable | E6 `tune_wavelength_calibration_to_skylines()` | 2 | no | — (no defect) |
| 56 | B905 | `soxspipe/commonutils/image_transformer.py:155` | `zip()` without `strict=`. Equal length by construction: all three lists built from the same `orderTraces` (:578-581) | reachable | E6 `cache_image()` | 2 | no | — (no defect) |
| 57 | B905 | `soxspipe/commonutils/image_transformer.py:507` | `zip()` without `strict=`. Equal length by construction: two arrays of one grid | reachable | E6 `_determine_rectified_image_boundaries()` | 2 | no | — (no defect) |
| 58 | B905 | `soxspipe/commonutils/image_transformer.py:512` | `zip()` without `strict=`. Equal length by construction: two arrays of one grid | reachable | E6 `_determine_rectified_image_boundaries()` | 2 | no | — (no defect) |
| 59 | B905 | `soxspipe/commonutils/image_transformer.py:532` | `zip()` without `strict=`. Equal length by construction: columns of the spectral-format table | reachable | E6 `_determine_rectified_image_boundaries()` | 2 | no | — (no defect) |
| 60 | B905 | `soxspipe/commonutils/image_transformer.py:679` | `zip()` without `strict=`. Equal length by construction: built from the same `orderTraces` | reachable | E6 `get_order_rectified()` | 2 | no | — (no defect) |
| 61 | B905 | `soxspipe/commonutils/phase3.py:90` | `zip()` without `strict=`. Equal length by construction: three lists from one header walk (ruff reports this line twice: inner and outer `zip`) | reachable | E7 `_write()` → `sort_keywords()` | 2 | no | — (no defect) |
| 62 | B905 | `soxspipe/commonutils/phase3.py:90` | `zip()` without `strict=`. Equal length by construction: three lists from one header walk (ruff reports this line twice: inner and outer `zip`) | reachable | E7 `_write()` → `sort_keywords()` | 2 | no | — (no defect) |
| 63 | B905 | `soxspipe/commonutils/phase3.py:93` | `zip()` without `strict=`. Equal length by construction: tuple-unzipped from one sort | reachable | E7 `_write()` → `sort_keywords()` | 2 | no | — (no defect) |
| 64 | B905 | `soxspipe/commonutils/phase3.py:130` | `zip()` without `strict=`. Equal length by construction: four columns of the QC frame | reachable | E1/E3/E6 `write_fits_table_to_disk()` | 2 | no | — (no defect) |
| 65 | B007 | `soxspipe/commonutils/reducer.py:222` | Unused `index`; cosmetic | reachable | E9 `reducer.reduce()` | 2 | no | — (no defect) |
| 66 | B905 | `soxspipe/commonutils/reducer.py:648` | `zip()` without `strict=`. Equal length by construction: lists built together | reachable | E9 `run_recipe_bulk()` | 2 | no | — (no defect) |
| 67 | B905 | `soxspipe/commonutils/reducer.py:747` | `zip()` without `strict=`. Equal length by construction: appended together (:715-716) | reachable | E9 `run_recipe_bulk()` | 2 | no | — (no defect) |
| 68 | B905 | `soxspipe/commonutils/response_function.py:78` | `zip()` without `strict=`. Equal length by construction: `polyval` over the same wavelength array | reachable | E-std `_fit_response_polynomial()` | 2 | no | — (no defect) |
| 69 | B905 | `soxspipe/commonutils/response_function.py:709` | `zip()` without `strict=`. Equal length by construction: four columns of the QC frame | reachable | E-std `write_response_function_to_file()` | 2 | no | — (no defect) |
| 70 | B905 | `soxspipe/commonutils/subtract_background.py:296` | `zip()` without `strict=`. Equal length by construction: tuple-unzip | reachable | E10 `mask_order_locations()` | 2 | no | — (no defect) |
| 71 | B905 | `soxspipe/commonutils/subtract_background.py:299` | `zip()` without `strict=`. Equal length by construction: columns of one order's pixel table | reachable | E10 `mask_order_locations()` | 2 | no | — (no defect) |
| 72 | B905 | `soxspipe/commonutils/subtract_background.py:306` | `zip()` without `strict=`. Equal length by construction: tuple-unzipped from one comprehension | reachable | E10 `mask_order_locations()` | 2 | no | — (no defect) |
| 73 | B905 | `soxspipe/commonutils/subtract_background.py:324` | `zip()` without `strict=`. Equal length by construction: derived from the same unzipped tuple | reachable | E10 `mask_order_locations()` | 2 | no | — (no defect) |
| 74 | B905 | `soxspipe/commonutils/subtract_background.py:336` | `zip()` without `strict=`. Equal length by construction: derived from the same unzipped tuple | reachable | E10 `mask_order_locations()` | 2 | no | — (no defect) |
| 75 | B007 | `soxspipe/commonutils/subtract_sky.py:260` | Unused `o`; cosmetic | reachable | E5 `subtract()` | 2 | no | — (no defect) |
| 76 | B905 | `soxspipe/commonutils/subtract_sky.py:260` | `zip()` without `strict=`. Equal length by construction: one dataframe per `uniqueOrders` entry (:229-239) | reachable | E5 `subtract()` | 2 | no | — (no defect) |
| 77 | B905 | `soxspipe/commonutils/subtract_sky.py:276` | `zip()` without `strict=`. Equal length by construction: appended once per order in the loop at :260 | reachable | E5 `subtract()` | 2 | no | — (no defect) |
| 78 | B905 | `soxspipe/commonutils/subtract_sky.py:503` | `zip()` without `strict=`. Equal length by construction: two columns of one dataframe | reachable | E12/E13 `plot_sky_sampling()` (`--debug` or `sky_model_qc_plot: True`) | 2 | no | — (no defect) |
| 79 | B905 | `soxspipe/commonutils/subtract_sky.py:583` | `zip()` without `strict=`. Equal length by construction: four literal lists of equal length | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 80 | B905 | `soxspipe/commonutils/subtract_sky.py:669` | `zip()` without `strict=`. Equal length by construction: four literal lists of equal length | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 81 | B905 | `soxspipe/commonutils/subtract_sky.py:715` | `zip()` without `strict=`. Equal length by construction: four literal lists of equal length | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 82 | B905 | `soxspipe/commonutils/subtract_sky.py:718` | `zip()` without `strict=`. Equal length by construction: two columns under one mask | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 83 | B905 | `soxspipe/commonutils/subtract_sky.py:829` | `zip()` without `strict=`. Equal length by construction: three columns of one dataframe | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 84 | B905 | `soxspipe/commonutils/subtract_sky.py:867` | `zip()` without `strict=`. Equal length by construction: three columns of one dataframe | reachable | E12/E13 `plot_sky_sampling()` | 2 | no | — (no defect) |
| 85 | B904 | `soxspipe/commonutils/subtract_sky.py:1480` | Re-raise without `from`; Python still chains the context implicitly. Cosmetic | reachable | E5 `fit_bspline_curve_to_sky()` | 2 | no | — (no defect) |
| 86 | B905 | `soxspipe/commonutils/subtract_sky.py:1701` | `zip()` without `strict=`. Equal length by construction: three columns of one dataframe | reachable | E5 `add_data_to_placeholder_images()` | 2 | no | — (no defect) |
| 87 | B905 | `soxspipe/commonutils/subtract_sky.py:1711` | `zip()` without `strict=`. Equal length by construction: three columns of one dataframe | reachable | E5 `add_data_to_placeholder_images()` | 2 | no | — (no defect) |
| 88 | B905 | `soxspipe/commonutils/subtract_sky.py:1720` | `zip()` without `strict=`. Equal length by construction: three columns of one dataframe | reachable | E5 `add_data_to_placeholder_images()` | 2 | no | — (no defect) |
| 89 | B905 | `soxspipe/commonutils/subtract_sky.py:1768` | `zip()` without `strict=`. Equal length by construction: two columns of one dataframe | reachable | E13 `plot_image_comparison()` (`sky_model_qc_plot: True`) | 2 | no | — (no defect) |
| 90 | B905 | `soxspipe/commonutils/subtract_sky.py:2113` | `zip(bins[e:-e], result[e:-e])`: lengths differ by one *by design* (100 edges vs 99 `value_counts` intervals); pairs each interval with its left edge, correctly. Adding `strict=True` here would break object masking | reachable | E5 `clip_object_slit_positions()` (`aggressive_object_masking`) | 2 | no | — (no defect; do **not** add `strict=True`) |
| 91 | B905 | `soxspipe/commonutils/subtract_sky.py:2554` | `zip()` without `strict=`. Equal length by construction: matplotlib's own handles/labels pair | reachable | E12 `refresh_and_plot()` (`--debug` only) | 2 | no | — (no defect) |
| 92 | B905 | `soxspipe/recipes/base_recipe.py:1023` | `zip()` without `strict=`. Equal length by construction: both lists cut by the same `valid` mask | reachable | E7 `_verify_single_gain()` (SOXS branch) | 2 | no | — (no defect) |
| 93 | B905 | `soxspipe/recipes/base_recipe.py:1351` | `zip()` without `strict=`. Equal length by construction: four columns of the QC frame | reachable | E7 `_write()` | 2 | no | — (no defect) |
| 94 | B905 | `soxspipe/recipes/base_recipe.py:2384` | `zip()` without `strict=`. Equal length by construction: three columns of one table | reachable | E7 `_stamp_raw_frame_records()` | 2 | no | — (no defect) |
| 95 | B905 | `soxspipe/recipes/base_recipe.py:2425` | `zip()` without `strict=`. Equal length by construction: four columns of one table | reachable | E7 `_stamp_calibration_frame_records()` | 2 | no | — (no defect) |
| 96 | B006 | `soxspipe/recipes/soxs_disp_solution.py:56` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_disp_solution.__init__` | 2 | no | — (no defect) |
| 97 | B006 | `soxspipe/recipes/soxs_mbias.py:57` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_mbias.__init__` | 2 | no | — (no defect) |
| 98 | B006 | `soxspipe/recipes/soxs_mdark.py:58` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_mdark.__init__` | 2 | no | — (no defect) |
| 99 | B905 | `soxspipe/recipes/soxs_mdark.py:240` | `zip()` without `strict=`. Equal length by construction: `zip(*…)` unzip of 3-tuples | reachable | `soxs_mdark` → `_combine_dark_frames()` | 2 | no | — (no defect) |
| 100 | B006 | `soxspipe/recipes/soxs_mflat.py:71` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_mflat.__init__` | 2 | no | — (no defect) |
| 101 | B905 | `soxspipe/recipes/soxs_mflat.py:364` | `zip()` without `strict=`. Equal length by construction: four literal 4-element lists | reachable | `soxs_mflat` → `produce_product()` | 2 | no | — (no defect) |
| 102 | B905 | `soxspipe/recipes/soxs_mflat.py:1115` | `zip()` without `strict=`. Equal length by construction: two columns of the order table | reachable | `soxs_mflat` → `_order_centre_mask()` | 2 | no | — (no defect) |
| 103 | B905 | `soxspipe/recipes/soxs_mflat.py:1121` | `zip()` without `strict=`. Equal length by construction: two columns of the order table (axis-swapped names are correct for `axisA == "y"`) | reachable | `soxs_mflat` → `_order_centre_mask()` | 2 | no | — (no defect) |
| 104 | B905 | `soxspipe/recipes/soxs_mflat.py:1473` | `zip()` without `strict=`. Equal length by construction: four columns of the order table | reachable | `soxs_mflat` → `mask_low_sens_pixels()` | 2 | no | — (no defect) |
| 105 | B006 | `soxspipe/recipes/soxs_nod.py:64` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_nod.__init__` | 2 | no | — (no defect) |
| 106 | B006 | `soxspipe/recipes/soxs_offset.py:59` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_offset.__init__` | 2 | no | — (no defect) |
| 107 | B905 | `soxspipe/recipes/soxs_offset.py:367` | `zip()` without `strict=`. Equal length by construction: appended together (:269-270) | reachable | `soxs_offset` → `_split_on_off_frames()` | 2 | no | — (no defect) |
| 108 | B006 | `soxspipe/recipes/soxs_order_centres.py:59` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_order_centres.__init__` | 2 | no | — (no defect) |
| 109 | B006 | `soxspipe/recipes/soxs_spatial_solution.py:51` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_spatial_solution.__init__` | 2 | no | — (no defect) |
| 110 | B006 | `soxspipe/recipes/soxs_stare.py:55` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result before any in-place `.sort`) | reachable | CLI `soxspipe reduce` → `soxs_stare.__init__` | 2 | no | — (no defect) |
| 111 | B006 | `soxspipe/recipes/soxs_straighten.py:59` | `inputFrames=[]` default; never mutated (rebound to `sof.get()` result) | reachable | E11 `soxs_straighten` (Python API only; not wired into the CLI or reducer) | 2 | no | — (no defect) |

## Proposed issues

Every reproducing test below was run from the repository root in the Python 3.12 venv. The test files live outside the repository, in the session scratchpad (`$SCRATCH` in the commands), and are **not** part of the suite. None is marked `xfail` or `skip`. `tests.conftest.RecordingLogger` and `tests.factories` are imported only to build inputs.

### I-1: Workspace path beginning with `~` crashes the data organiser with UnboundLocalError

**Symptom**
`data_organiser(log=log, rootDir="~/my-workspace")` raises `UnboundLocalError: cannot access local variable 'directory'`. The class docstring documents this constructor as the way to build an organiser. The home-expansion branch assigns `directory = directory.replace(...)` instead of updating `rootDir`, so a fix to the name alone would still leave `~` unexpanded. The CLI is not affected because it `chdir`s into the workspace and passes `"."`.

**Evidence** (tier 1)

```python
from tests.conftest import RecordingLogger

from soxspipe.commonutils.data_organiser import data_organiser


def test_tilde_root_dir_is_expanded():
    # THE DOCUMENTED PYTHON API ACCEPTS rootDir; A "~"-PREFIXED PATH SHOULD BE EXPANDED, NOT CRASH
    do = data_organiser(log=RecordingLogger(), rootDir="~/dy34-workspace", dbConnect=False)
    assert not do.rootDir.startswith("~")
```

```bash
python -m pytest $SCRATCH/dy34/test_repro_data_organiser_tilde.py -p no:cacheprovider -q
```

```text
E           UnboundLocalError: cannot access local variable 'directory' where it is not associated with a value
soxspipe/commonutils/data_organiser.py:207: UnboundLocalError
```

**Location**
`soxspipe/commonutils/data_organiser.py:207`

**Not fixed here**
This finding is filed, not fixed, by DY-34.

### I-2: Optimal extraction crashes with NameError when a recipe is run from a frame list instead of a SOF file

**Symptom**
The recipe docstrings say `inputFrames` can be "a directory, a set-of-files (SOF) file or a list of fits frame paths". When a list is passed, `base_recipe` sets `sofName = False`. `soxs_stare` and `subtract_sky` both handle that with a `filenamer` fallback (DY-111). `horne_extraction.__init__` has the same fallback but never imports `filenamer`, so extraction dies with `NameError`. (Related, not a ruff finding: `base_recipe.py:1983` would later fail on `self.sofName + ".sof"` on the same path.)

**Evidence** (tier 1)

```python
import numpy as np
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.horne_extraction import horne_extraction
from tests.conftest import RecordingLogger
from tests.factories import instrument_header


def test_horne_extraction_names_products_without_a_sof_name(monkeypatch, tmp_path):
    # SKIP THE DISPERSION-MAP/2D-MAP LOADING IN base_util; ONLY THE FILENAME-TEMPLATE BRANCH IS UNDER TEST
    def _fake_base_init(self, log, settings, **kwargs):
        self.log = log
        self.settings = settings

    monkeypatch.setattr(horne_extraction.__mro__[1], "__init__", _fake_base_init)
    frame = CCDData(
        np.ones((8, 8)),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros((8, 8), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((8, 8)), unit=u.electron),
    )
    recipeSettings = {
        "horne-extraction-slit-length": 20,
        "horne-extraction-profile-clipping-sigma": 3,
        "horne-extraction-profile-clipping-iteration-count": 5,
        "horne-extraction-profile-global-clipping-sigma": 25,
    }

    # soxs_stare SETS sofName=False WHEN IT IS HANDED A LIST OF FRAMES RATHER THAN A .sof FILE
    extractor = horne_extraction(
        log=RecordingLogger(),
        settings={},
        recipeSettings=recipeSettings,
        skySubtractedFrame=frame,
        unflattenedFrame=frame,
        twoDMapPath=None,
        sofName=False,
    )

    assert extractor.filenameTemplate.endswith(".fits")
```

```bash
python -m pytest $SCRATCH/dy34/test_repro_horne_no_sof_name.py -p no:cacheprovider -q
```

```text
E           NameError: name 'filenamer' is not defined
soxspipe/commonutils/horne_extraction.py:179: NameError
```

**Location**
`soxspipe/commonutils/horne_extraction.py:179`

**Not fixed here**
This finding is filed, not fixed, by DY-34.

### I-3: Background subtraction crashes with UnboundLocalError on a single-order order table

**Symptom**
`mask_order_locations()` iterates orders in ascending order. For the reddest order (`o == oBot`) it sets `expandBottom = expandTop` *before* computing this order's `expandTop`. With a single-order table, `expandTop` has never been assigned and the call raises `UnboundLocalError`. The existing unit tests that use a one-order table monkeypatch `mask_order_locations` away (`tests/unit/test_background_subtraction.py:135`, `:201`), so the suite does not catch this. In multi-order tables the same line reads the *previous* order's `expandTop`. By the geometry argument in the summary, that stale read changes no pixel mask (tier 2, not tested). Whoever fixes this should confirm the real detector orientation.

**Evidence** (tier 1)

```python
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.nddata import CCDData, StdDevUncertainty

from soxspipe.commonutils.subtract_background import subtract_background
from tests.conftest import RecordingLogger
from tests.factories import instrument_header


def test_mask_order_locations_handles_a_single_order_table():
    worker = object.__new__(subtract_background)
    worker.log = RecordingLogger()
    worker.frame = CCDData(
        np.full((32, 32), 10.0, dtype=np.float32),
        unit=u.electron,
        meta=instrument_header(),
        mask=np.zeros((32, 32), dtype=bool),
        uncertainty=StdDevUncertainty(np.ones((32, 32)), unit=u.electron),
    )
    worker.axisA, worker.axisB = "x", "y"
    orderPixels = pd.DataFrame({"order": [10], "ycoord": [5], "xcoord_edgeup": [24.0], "xcoord_edgelow": [20.0]})

    # ONE ORDER IS BOTH oTop AND oBot, SO THE "o == oBot" BRANCH READS expandTop BEFORE IT IS EVER ASSIGNED
    worker.mask_order_locations(orderPixels)

    assert worker.frame.mask[5, 20:24].all()
```

```bash
python -m pytest $SCRATCH/dy34/test_repro_background_single_order.py -p no:cacheprovider -q
```

```text
E               UnboundLocalError: cannot access local variable 'expandTop' where it is not associated with a value
soxspipe/commonutils/subtract_background.py:279: UnboundLocalError
```

**Location**
`soxspipe/commonutils/subtract_background.py:279`

**Not fixed here**
This finding is filed, not fixed, by DY-34.

### I-4: Sky subtraction crashes with UnboundLocalError instead of reverting when the first b-spline fit is poor

**Symptom**
When `scipy.interpolate.splrep` returns `ier` 10 or 30, `fit_bspline_curve_to_sky()` logs "Reverting to last iteration" and sets `tck = tck_previous`. On the first pass `tck_previous` does not exist yet, so the stare reduction dies with `UnboundLocalError` instead of reporting the poor fit. The test below uses the **real** `splrep` (no mocking): coincident quantile-placed starter knots are enough to trigger it. scipy 1.18.1 returns `ier=30` for coincident knots rather than raising. Side note, not a ruff finding: the same probe returned `ier=50` for other bad knot vectors, which this branch does not handle at all, so the fitted `tck` is used silently. That deserves a look when this is fixed.

**Evidence** (tier 1)

```python
import numpy as np
import pandas as pd

from soxspipe.commonutils.subtract_sky import subtract_sky
from tests.conftest import RecordingLogger


def test_poor_first_bspline_fit_is_reported_not_crashed():
    worker = object.__new__(subtract_sky)
    worker.log = RecordingLogger()
    worker.arm = "NIR"
    worker.binx = worker.biny = 1
    worker.debug = False
    worker.bspline_order = 3
    worker.recipeSettings = {
        "sky-subtraction": {
            "min_points_per_knot": 5,
            "bspline_fitting_residual_clipping_sigma": 3,
            "bspline_iteration_limit": 7,
            "residual_floor_percentile": 50,
            "starting_points_per_knot": 25,
        }
    }
    rng = np.random.default_rng(1)
    # 900 OF 1200 PIXELS SHARE ONE WAVELENGTH, SO THE QUANTILE-PLACED STARTER KNOTS COINCIDE AND
    # SCIPY'S splrep (UNMOCKED) RETURNS A "POOR FIT" ier CODE ON THE VERY FIRST PASS
    wavelength = np.concatenate([np.full(900, 1000.0), np.linspace(1000.1, 1010.0, 300)])
    imageMapOrder = pd.DataFrame(
        {
            "order": 12,
            "wavelength": wavelength,
            "flux": rng.normal(100.0, 5.0, wavelength.size),
            "residual_windowed_std": np.full(wavelength.size, 5.0),
            "flagged_all_clipped": False,
            "flagged_noisy_region": False,
        }
    )

    # THE METHOD HANDLES ier IN (10, 30) BY "REVERTING TO LAST ITERATION", BUT ON THE FIRST PASS
    # THERE IS NO LAST ITERATION
    worker.fit_bspline_curve_to_sky(imageMapOrder)
```

```bash
python -m pytest $SCRATCH/dy34/test_repro_sky_bspline_first_pass.py -p no:cacheprovider -q
```

```text
E               UnboundLocalError: cannot access local variable 'tck_previous' where it is not associated with a value
soxspipe/commonutils/subtract_sky.py:1489: UnboundLocalError
```

**Location**
`soxspipe/commonutils/subtract_sky.py:1489`

**Not fixed here**
This finding is filed, not fixed, by DY-34.

### I-5: Failed QC database inserts are silently dropped after seven retries, so QC failures may not block downstream recipes

**Symptom**
Ruff flags `raise Exception(e)` at `base_recipe.py:2610` for missing `from`. Tracing it shows the line can never run. The loop is `while keepTrying < 7` and the guard is `if keepTrying > keepTryingMax - 1` (`> 6`), which is never true inside the loop. After seven failed `to_sql` attempts, for example on a database locked by parallel `reduce all -m` workers, the method returns normally and the QC rows are lost. The docstring promises "Exception if the insert fails after seven attempts". `data_organiser._dataframe_to_sqlite` has the same loop and does raise (`if keepTrying > 5: raise`). The `quality_control` table drives `product_frames.status_base = "fail"` (`data_organiser.py:692-693`), so a lost row can let a QC-failed product feed downstream recipes. The B904 complaint itself (missing `from`) is cosmetic. The defect is the unreachable raise.

**Evidence** (tier 1)

```python
import sqlite3
import time

import pandas as pd
import pytest

from soxspipe.recipes.base_recipe import base_recipe
from tests.conftest import RecordingLogger


def test_insert_that_fails_every_attempt_raises(monkeypatch):
    recipe = object.__new__(base_recipe)
    recipe.log = RecordingLogger()
    # A CLOSED CONNECTION MAKES EVERY to_sql ATTEMPT FAIL
    recipe.conn = sqlite3.connect(":memory:")
    recipe.conn.close()
    sleeps = []
    monkeypatch.setattr(time, "sleep", sleeps.append)

    # THE DOCSTRING PROMISES "Exception if the insert fails after seven attempts"
    with pytest.raises(Exception):
        recipe._dataframe_to_sqlite(pd.DataFrame({"qc_name": ["X"]}), "quality_control")

    assert len(sleeps) == 7
```

```bash
python -m pytest $SCRATCH/dy34/test_repro_sqlite_insert_swallowed.py -p no:cacheprovider -q
```

```text
E       Failed: DID NOT RAISE Exception
```

**Location**
`soxspipe/recipes/base_recipe.py:2610`

**Not fixed here**
This finding is filed, not fixed, by DY-34.

### I-6: Order-centre QC plot mislabels order colours and drops orders when one order's fit is off-detector

**Symptom**
In `detect_continuum.plot_results()`, the middle panel loops over `uniqueOrders` and `continue`s (`:1253`) past any order whose fitted trace lies wholly off the detector. It appends a colour only for orders it plots, and it records them in `foundOrders`. The three lower panels then `zip(uniqueOrders, colors)`. After a skipped order, every later order is drawn in its neighbour's colour, and the last order is silently dropped from those panels. This affects the QC PDF only. No product values change.

**Evidence** (tier 2)
No cheap test. `plot_results()` needs a full order-pixel table, poly table, clipped-data frame, detector parameters and a dispersion map, and the defect only shows up in matplotlib artist colours in a written PDF. The argument: `colors.append(...)` (`:1271`) sits after the `continue` (`:1253`) in the loop over `uniqueOrders` (`:1238`), so `len(colors) < len(uniqueOrders)` whenever any order is skipped. `zip` then truncates and pairs index *i* of `uniqueOrders` with the colour of the *i*-th *plotted* order. `foundOrders` already holds the matching order list. The same pattern in `detect_order_edges.py:835` is dead (inside `if False:`).

**Location**
`soxspipe/commonutils/detect_continuum.py:1314` (also `:1345`, `:1376`)

**Not fixed here**
This finding is filed, not fixed, by DY-34.

## Housekeeping (dead code, low severity)

These lines cannot run in any reduction today. They are proposed for a housekeeping list, not as bugs. None is an X-Shooter path.

| # | Rule | `path:line` | Description | Status | Why it cannot run | Tier | Wrong science | Proposed |
|---|---|---|---|---|---|---|---|---|
| 112 | B905 | `soxspipe/commonutils/create_dispersion_map.py:4337` | `zip()` without `strict=` | dead | `update_static_line_list_detector_positions()` has no caller (only a commented-out call at :995) | 2 | no | housekeeping |
| 113 | B905 | `soxspipe/commonutils/create_dispersion_map.py:4778` | `zip()` without `strict=` | dead | `straighten_mph_sets()` only called under `if self.firstGuessMap and False:` (:929) | 2 | no | housekeeping |
| 114 | B007 | `soxspipe/commonutils/data_organiser.py:500` | Unused `dirs` | dead | inside `if False:` (error-log/SOF deletion) | 2 | no | housekeeping |
| 115 | B905 | `soxspipe/commonutils/detect_order_edges.py:835` | `zip(uniqueOrders, colors)` (same misalignment pattern as I-6) | dead | inside `if False:` at :833 | 2 | no | housekeeping |
| 116 | F821 | `soxspipe/commonutils/horne_extraction.py:662` | Undefined `axisBVals` in the pixel-shift branch | dead | `calibrationCol` is pinned by `if True:` at :602, so the `else` branch is never taken | 2 | no | housekeeping |
| 117 | F821 | `soxspipe/commonutils/subtract_sky.py:1400` | Undefined `tck` passed to `adjust_tilt` | dead | guarded by `… and False` at :1398 | 2 | no | housekeeping |
| 118 | B905 | `soxspipe/commonutils/subtract_sky.py:1872` | `zip()` without `strict=` | dead | `subtract_sky._rectify_order()` has no caller | 2 | no | housekeeping |
| 119 | F821 | `soxspipe/commonutils/subtract_sky.py:2584` | Undefined `plot` (`plot.show()`) | dead | `refresh_and_plot()` is only ever called with `pause=0.2`, so the `else` is never reached; the enclosing plot is `--debug` only | 2 | no | housekeeping |
| 120 | F821 | `soxspipe/commonutils/subtract_sky.py:2940` | Undefined `sky_flux_rolling_median` | dead | inside `if False:` | 2 | no | housekeeping |
| 121 | F821 | `soxspipe/commonutils/subtract_sky.py:2951` | Undefined `sky_d1` | dead | inside `if False:` | 2 | no | housekeeping |
| 122 | F821 | `soxspipe/commonutils/subtract_sky.py:2963` | Undefined `sky_d3` | dead | inside `if False:` | 2 | no | housekeeping |

### Housekeeping notes for the reachable non-defects

- **F507 ×4** (`create_dispersion_map.py:2449/2464/2478/2501`): drop the redundant `% locals()`. Behaviour is unchanged either way.
- **F811** (`data_organiser.py:705`): the duplicate local import is harmless. Do not let an auto-fixer remove function-local imports (see the DY-34 brief).
- **B905**: do **not** blanket-apply `strict=True`. `subtract_sky.py:2113` zips 100 bin edges against 99 interval counts on purpose, and would raise.
