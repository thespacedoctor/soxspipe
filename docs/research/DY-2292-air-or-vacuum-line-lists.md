# DY-2292: air or vacuum wavelengths in the SOXS arc line lists

These notes answer Linear issue DY-2292 for the ESO Phase 3 map (DY-2291). They are for the people and agents who will set `TUCD1`, `WAVELMIN` and `WAVELMAX` in the `soxspipe` Phase 3 products.

## Answer

The SOXS arc line lists hold standard dry-air wavelengths, in nm, in both arms, so every `soxspipe` SOXS spectrum has an air wavelength scale. The Phase 3 products must therefore use `TUCD1 = 'em.wl;obs.atmos'` (standard section 8.2, Req-8), and `WAVELMIN`/`WAVELMAX` must be air values. The evidence is a line-by-line match against the NIST Atomic Spectra Database (ASD): every line in the four VIS pinhole tables (Th, Ar) agrees with the NIST air value to within 0.1 pm, and 90 to 93% of the lines in the two NIR pinhole tables (Ar, Ne, Xe, Hg) agree with NIST air to within 6 pm. No more of them sit near a vacuum value than chance allows, and the vacuum values are 100 to 550 pm away. The sky-line list that the pipeline uses to tune the science wavelength scale is also in air, and no code converts between air and vacuum, so the product scale stays in air. The project docs already say "The final wavelength solution is provided in air". The UV-VIS and NIR arms agree. The only vacuum list is `ThAr.fits`, which only the X-Shooter settings use.

**Confidence: high.** The match is unambiguous for every line list the SOXS settings use. One small detail is still open and does not change the answer: the NIR lists that come from the PypeIt-style NIST files sit 1 to 4 pm above NIST air. That is less than 3% of the air-to-vacuum difference (see [A small offset in the NIR lists](#a-small-offset-in-the-nir-lists)).

## Contents

- [What the standard asks for](#what-the-standard-asks-for)
- [Which line lists the SOXS settings use](#which-line-lists-the-soxs-settings-use)
- [Method](#method)
- [Result for every list](#result-for-every-list)
- [Strong-line cross-check](#strong-line-cross-check)
- [The sky-line list is also in air](#the-sky-line-list-is-also-in-air)
- [Repository history and docs](#repository-history-and-docs)
- [A small offset in the NIR lists](#a-small-offset-in-the-nir-lists)
- [Lists that differ from the others](#lists-that-differ-from-the-others)
- [What this means for other tickets](#what-this-means-for-other-tickets)
- [Sources](#sources)

## What the standard asks for

Section 8.2, Req-8 of the ESO Phase 3 standard says: "To discriminate wavelength measured in dry air (at standard temperature and pressure) from the wavelength measured in vacuum, TUCD1 (see section 5.18) must be set to 'em.wl;obs.atmos' and 'em.wl', respectively" (`output/p3sdpstd_md/p3sdpstd.md`, line 1323). The examples in the standard add the comment `/ em.wl: vacuum, em.wl;obs.atmos: air` (lines 1279 and 1304).

NIST ASD air wavelengths are standard dry air. The ASD help page says the refractive index comes from "the five-parameter formula given by E.R. Peck and K. Reeder" (J. Opt. Soc. Am. 62, 958, 1972), and that air wavelengths are vacuum wavelengths divided by that index [S1]. That is the same convention Req-8 calls "dry air (at standard temperature and pressure)".

## Which line lists the SOXS settings use

`soxspipe/resources/soxs_detector_parameters.yaml` names the line lists for each arm.

| Arm | Setting | File |
|---|---|---|
| VIS | `predicted pinhole lines: single: 1x1` | `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` |
| VIS | `predicted pinhole lines: single: 2x1` | `SOXS_SPH_TAB_VIS_ThAr_2x1.fits` |
| VIS | `predicted pinhole lines: multi: 1x1` | `SOXS_MPH_TAB_VIS_ThAr_1x1.fits` |
| VIS | `predicted pinhole lines: multi: 2x1` | `SOXS_MPH_TAB_VIS_ThAr_2x1.fits` |
| VIS | `line-atlas` | `thar_murphy.fits` |
| NIR | `predicted pinhole lines: single: 1x1` | `ArHgNeXe_clean_within_2.0pixel.fits` |
| NIR | `predicted pinhole lines: multi: 1x1` | `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` |
| NIR | `line-atlas` | `ArHgNeXe.fits` |

The pinhole tables set the wavelength scale. `create_dispersion_map.get_predicted_line_list` loads them and fits the dispersion solution to their `wavelength` column. The `line-atlas` file is read only by `create_new_static_line_list`, which runs only on the `bootstrap_dispersion_solution` path (`soxspipe/commonutils/create_dispersion_map.py`, lines 1565 to 1568 and 4887 to 5023). The atlases were checked too, because the NIR pinhole tables were built from `ArHgNeXe.fits`.

## Method

1. I downloaded the NIST ASD line tables for Ar I, Ar II, Ne I, Ne II, Hg I, Xe I, Xe II, Kr I, Th I and Th II from 280 nm to 2600 nm. I took each table twice: once with the default wavelength type (vacuum below 200 nm, air from 200 nm to 2000 nm, vacuum above 2000 nm [S1]) and once in vacuum for all wavelengths [S2]. I kept both the observed and the Ritz wavelengths.
2. Above 2000 nm the default output is vacuum. For that range I made the air values myself by dividing the NIST vacuum values by the Peck and Reeder index [S1, S3].
3. For each unique wavelength in a SOXS list, I found the nearest NIST air line and the nearest NIST vacuum line. I counted a match when the distance was 6 pm (0.006 nm) or less.
4. To measure the chance rate, I shifted each list by +0.1 nm and counted again. This matters for the dense Th lists, where a random wavelength often falls near some Th line.

I read the FITS files with the `soxspipe` conda env Python and astropy, from `origin/develop` at commit `ad73aeb0`.

## Result for every list

The air-to-vacuum difference is 114 pm at 405 nm, 278 pm at 1014 nm and 552 pm at 2026 nm (from the NIST pairs). A 6 pm window separates the two conventions clearly.

| File | Used for | Unique lines | Within 6 pm of NIST air | Within 6 pm of NIST vacuum | Chance level |
|---|---|---|---|---|---|
| `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` | VIS single pinhole 1x1 | 132 | 100% | 38% | 34% |
| `SOXS_SPH_TAB_VIS_ThAr_2x1.fits` | VIS single pinhole 2x1 | 132 | 100% | 38% | 34% |
| `SOXS_MPH_TAB_VIS_ThAr_1x1.fits` | VIS multi pinhole 1x1 | 132 | 100% | 38% | 34% |
| `SOXS_MPH_TAB_VIS_ThAr_2x1.fits` | VIS multi pinhole 2x1 | 132 | 100% | 38% | 34% |
| `thar_murphy.fits` | VIS line atlas (bootstrap only) | 3070 | 94% | 36% | 37% |
| `ArHgNeXe_clean_within_2.0pixel.fits` | NIR single pinhole | 551 | 93% | 7% | 8% |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | NIR multi pinhole | 507 | 90% | 7% | 8% |
| `ArHgNeXe.fits` | NIR line atlas (bootstrap only) | 5283 | 97% | 26% | 26% |
| `SOXS_SPH_TAB_NIR_Xe.fits` | not used (commented out) | 101 | 95% | 10% | 10% |
| `SOXS_MPH_TAB_NIR_Xe.fits` | not used (commented out) | 101 | 94% | 10% | 10% |
| `SOXS_SPH_TAB_VIS_ThAr_1x1_pre_20241001.fits` | not used | 96 | 100% | 38% | 46% |
| `SOXS_MPH_TAB_VIS_ThAr_1x1_old.fits` | not used | 115 | 37% | 40% | 33% |
| `Ar.fits`, `Ne.fits`, `Hg.fits`, `Xe.fits` | not used | 250 to 850 | 82 to 98% | 7 to 12% | 7 to 12% |
| `ThAr.fits` | X-Shooter line atlas only | 17099 | 38% | 98% | 38% |

The vacuum match rate of each SOXS list equals its chance rate, so the vacuum hits are random. With a tighter 1 pm window, all VIS pinhole lines still match NIST air (119 of 132 match air only, and none match vacuum only), and the median difference from NIST air is 0.00 to 0.01 pm.

Most of the NIR lines that do not match NIST air within 6 pm came from the Gemini GNIRS, Keck MOSFIRE and X-Shooter lists (the `filename` column in `ArHgNeXe.fits`). These lists give wavelengths to 0.01 nm or 0.001 nm, and some of their lines are blends or have no single NIST counterpart. They do not match NIST vacuum either, so they do not point toward vacuum.

## Strong-line cross-check

The table compares strong lines from the VIS and NIR pinhole tables that the settings use with their NIST values. "NIST type" says if the NIST value is the observed or the Ritz wavelength.

| SOXS table | Ion | SOXS wavelength (nm) | NIST air (nm) | NIST vacuum (nm) | SOXS minus air (pm) | SOXS minus vacuum (pm) | NIST type |
|---|---|---|---|---|---|---|---|
| `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` | Ar II | 405.2921 | 405.2921 | 405.4065 | +0.0 | -114.4 | observed |
| `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` | Th I | 545.2219 | 545.2219 | 545.3734 | +0.0 | -151.5 | Ritz |
| `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` | Th I | 764.7379 | 764.7379 | 764.9484 | +0.0 | -210.5 | Ritz |
| `SOXS_SPH_TAB_VIS_ThAr_1x1.fits` | Th I | 841.6727 | 841.6727 | 841.9040 | -0.0 | -231.3 | observed |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Hg I | 1013.9787 | 1013.9788 | 1014.2567 | -0.1 | -278.0 | Ritz |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Ne I | 1117.7538 | 1117.7525 | 1118.0585 | +1.4 | -304.7 | observed |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Ar I | 1530.1900 | 1530.1881 | 1530.6062 | +1.9 | -416.2 | Ritz |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Hg I | 1693.3600 | 1693.3595 | 1693.8220 | +0.5 | -462.0 | Ritz |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Ne I | 1957.7118 | 1957.7114 | 1958.2459 | +0.4 | -534.1 | Ritz |
| `SOXS_MPH_TAB_NIR_ArHgNeXe.fits` | Xe I | 2026.2250 | 2026.2243 | 2026.7774 | +0.7 | -552.4 | observed |

The last row is above 2000 nm, where the NIST default output is vacuum. The SOXS value is the air value, so the NIR list builders converted these lines too, and did not copy the NIST default output.

## The sky-line list is also in air

The science recipes tune the wavelength scale to sky lines when `tune-wavelength-to-skylines: True` (`soxspipe/soxs_default_settings.yaml`), and the NIR sky refit matches observed OH lines to the same catalogue. Both arms use `skylines.fits`. If that catalogue were in vacuum, the tuning would move the science scale toward vacuum. It is in air:

- VIS: [O I] 557.7347 nm in the list, NIST air 557.7339 nm, vacuum 557.8887 nm. Na I D2 588.9959 nm in the list, air 588.9950 nm, vacuum 589.1583 nm. [O I] 630.0309 nm in the list, air 630.0304 nm.
- NIR: I compared the 523 lines from 980 nm to 2290 nm with the OH line table of Oliva et al. (2015), which gives vacuum wavelengths [S4], converted to air with the Peck and Reeder index. Sample lines: 1202.0960 nm in the list, Oliva in air 1202.0950 nm, Oliva in vacuum 1202.4240 nm. 1623.1406 nm in the list, air 1623.1422 nm, vacuum 1623.5856 nm. 2245.3777 nm in the list, air 2245.3822 nm, vacuum 2245.9950 nm. The median difference from air is -1.1 pm.

I searched `soxspipe/` for any air-to-vacuum or vacuum-to-air conversion (`air_to_vac`, `airtovac`, `vactoair`, `refractive`, and similar names) and found none. `soxspipe/commonutils/phase3.py` does not write `TUCD1` yet.

## Repository history and docs

- `docs/source/utils/create_dispersion_map.md`, line 3: "The final wavelength solution is provided in air." This sentence came in with commit `4af272c0` (2025-01-22, "Many doc fixes").
- `ArHgNeXe.fits` (commit `27e86b55`, 2024-02-20, "add new ArHgNeXe line list") is a merged list. Its `filename` column names the input files, and `XSH_NIR_line_list_air.dat` states its convention in its name. The other inputs are named like the PypeIt arc line files (`ArI_lines.dat`, `*_DeVeny_lines.dat`, `*_lines_NIST_dry.dat`). Their `source` values are `NIST`, PypeIt solution names, or the GNIRS and MOSFIRE web pages. No column states air or vacuum.
- `ThAr.fits` (commit `790fa21c`, 2025-03-27, "fixing wavelength from A to nm") has `TCOMM2 = 'Expression: $2*0.1'`, which is a unit change only. Its values are vacuum.
- `thar_murphy.fits` (commit `a5addfe9`, 2025-08-22) was made from `thar_MM201006.txt` with `TCOMM1 = 'Expression: $2/10.'`, which is again only a unit change from angstrom to nm. Its values match NIST air.
- Commit `505bdef4` (DY-2084) moved the NIR multi-pinhole Xe I and Hg I lines to NIST Ritz values. The new values are also air, as the table above shows.

## A small offset in the NIR lists

The lines that come from the `*_lines_NIST_dry.dat`, `XeI_lines_NIST.dat` and `*_DeVeny_lines.dat` inputs are not exact NIST air copies. They sit above NIST air by an amount that falls with wavelength:

| Wavelength range (nm) | Median SOXS minus NIST air (pm) | As a fraction of the air-to-vacuum difference |
|---|---|---|
| 300 to 500 | +3.3 to +3.8 | 2.6 to 3.3% |
| 500 to 800 | +2.2 to +2.7 | 1.1 to 1.6% |
| 800 to 1100 | +1.6 to +1.9 | 0.6 to 0.8% |
| 1100 to 1500 | +1.2 | 0.35% |
| 1500 to 2000 | +0.8 to +0.9 | 0.2% |

The relative offset goes roughly as the inverse square of the wavelength (about 9 ppm at 400 nm and 0.5 ppm at 1750 nm). Converting NIST vacuum values with the Peck and Reeder, Edlén (1953), Ciddor (1996) or Morton (IAU) formulas does not remove it, because these formulas agree with each other to better than 0.1 pm here. The cause is not known. It could be a different conversion or rounding step in the upstream PypeIt-style files.

The offset does not change the answer. In the NIR wavelength range that the SOXS pinhole tables cover (about 800 nm to 2070 nm), it is 0.5 to 1.9 ppm, which is 0.15 to 0.6 km/s. The air-to-vacuum difference is about 274 ppm. The Xe I and Hg I lines that DY-2084 moved to NIST Ritz values have little or none of this offset (for example -0.1 pm at 1013.98 nm and +0.5 pm at 1693.36 nm).

## Lists that differ from the others

- `ThAr.fits` is in vacuum: 98% of its 17099 lines match NIST vacuum. Only `soxspipe/resources/xsh_detector_parameters.yaml` names it, as the X-Shooter UVB and VIS `line-atlas`, and it is read only on the bootstrap path. It does not affect SOXS products. If the X-Shooter path is used again, this list is not in the same convention as the X-Shooter pinhole tables. I did not check those tables.
- `SOXS_SPH_TAB_NIR_Xe.fits` and `SOXS_MPH_TAB_NIR_Xe.fits` are in air below 2000 nm, but their three lines above 2000 nm (2019.271, 2026.777 and 2028.797 nm) are the NIST vacuum values. This is the NIST default output above 2000 nm. The settings no longer use these two tables; their entries in `soxs_detector_parameters.yaml` are commented out.
- `SOXS_MPH_TAB_VIS_ThAr_1x1_old.fits` matches neither convention better than chance, so its wavelengths are probably shifted or wrongly identified. The settings do not use it.
- The UV-VIS and NIR arms do not differ: both use air.

## What this means for other tickets

- The Phase 3 spectrum writer must set `TUCD1 = 'em.wl;obs.atmos'`, with the comment `em.wl: vacuum, em.wl;obs.atmos: air`, on the wavelength column of every SOXS VIS and NIR spectrum. `phase3.py` does not write `TUCD1` yet.
- `WAVELMIN`, `WAVELMAX` and any other wavelength keywords (for example `SPEC_VAL` and `SPEC_BW`) can come straight from the product wavelength array. No conversion is needed, because the array is already in air.
- If the map decides to deliver vacuum wavelengths, that is a new conversion step, which would be the first one in the code. Nothing in the current pipeline makes vacuum wavelengths.
- Anyone who adds a new arc or sky list must convert lines above 2000 nm from the NIST default (vacuum) to air.

## Sources

- [S1] NIST Atomic Spectra Database, lines help page, "Choice of wavelength type" and the air-to-vacuum conversion: <https://physics.nist.gov/PhysRefData/ASD/Html/lineshelp.html> (read 2026-11-19).
- [S2] NIST Atomic Spectra Database lines form, queried 2026-11-19 for Ar I, Ar II, Ne I, Ne II, Hg I, Xe I, Xe II, Kr I, Th I and Th II from 280 nm to 2600 nm, CSV output, default wavelength type (`show_av=2`) and vacuum for all wavelengths (`show_av=3`): <https://physics.nist.gov/cgi-bin/ASD/lines1.pl>.
- [S3] E. R. Peck and K. Reeder, "Dispersion of air", J. Opt. Soc. Am. 62, 958 (1972).
- [S4] E. Oliva et al., "Lines and continuum sky emission in the near infrared", A&A 581, A47 (2015), VizieR catalogue J/A+A/581/A47, table 1 ("Vacuum wavelength of the first line of the doublet"): <https://vizier.cfa.harvard.edu/viz-bin/VizieR?-source=J/A+A/581/A47>.
- [S5] ESO Phase 3 standard, section 8.2, Req-8, as Markdown at `output/p3sdpstd_md/p3sdpstd.md` (line 1323) in the P32.11 workspace.
- [S6] `soxspipe` repository, `origin/develop` at `ad73aeb0`: `soxspipe/resources/soxs_detector_parameters.yaml`, `soxspipe/resources/xsh_detector_parameters.yaml`, `soxspipe/commonutils/create_dispersion_map.py`, `docs/source/utils/create_dispersion_map.md`, and the FITS files in `soxspipe/resources/static_calibrations/soxs/`.
