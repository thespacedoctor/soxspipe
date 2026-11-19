# DY-2296: where the target RA and DEC come from

These notes answer Linear ticket DY-2296 for the ESO Phase 3 map (DY-2291). They are for whoever sets the `RA` and `DEC` keywords of soxspipe Phase 3 spectra.

## Answer

The target position is in `HIERARCH ESO TEL TARG ALPHA` and `HIERARCH ESO TEL TARG DELTA`. These cards hold the coordinates typed into the observation block (OB), in packed sexagesimal form (`HHMMSS.TTT` and `±DDMMSS.TTT`), mean place, equinox and epoch J2000 (`TEL TARG COORDTYPE = 'M'`, `TEL TARG EQUINOX = 2000`, `TEL TARG EPOCH = 2000`, `TEL TARG EPOCHSYSTEM = 'J'` on every frame). Decoded to decimal degrees, they match SIMBAD J2000 positions to 0.01 arcsec or better for 12 of the 14 targets checked; the two misses (0.6 arcsec and 0.9 arcsec) are OB entry errors, not card errors. The card values do not change within an OB: nod, jitter and sky offsets move the telescope, not `TEL TARG`, so one value serves a whole stacked product. The plain `RA`/`DEC` cards are the telescope pointing (their own comment says `RA (J2000) pointing`). On 2211 SOXS science frames the pointing is a median 4.5 arcsec from `TEL TARG` (5 to 95 percent range 1.4 arcsec to 14 arcsec on stare frames), and up to 104 arcsec, because it also holds the acquisition and pointing-model residual, the nod or sky offset, and the proper motion applied by the telescope. So the pointing is not a usable target position, and footnote 10 of the standard forbids it anyway. ESO's own X-Shooter internal data products (pipeline `xshoo/3.6.0`) do not follow footnote 10: they copy the pointing `RA`/`DEC` of the first raw file unchanged, so there is no ESO precedent for `TEL TARG`. Two open points are for Dave: the epoch for high proper motion stars (the J2000 catalogue position can be more than 1 arcmin from where the star was when observed), and what to do when `TEL TARG` is missing.

**Confidence: high** that `TEL TARG ALPHA`/`DELTA` is the right source and is accurate to the precision of the OB entry. **Medium** on the edge cases (proper motion epoch, moving targets, blind or manual offsets during acquisition), which the headers do not fully describe.

## Contents

- [What the standard asks for](#what-the-standard-asks-for)
- [Survey method](#survey-method)
- [Candidate cards and their format](#candidate-cards-and-their-format)
- [Card values against pointing on sampled frames](#card-values-against-pointing-on-sampled-frames)
- [How accurate the target cards are](#how-accurate-the-target-cards-are)
- [Why the pointing differs from the target](#why-the-pointing-differs-from-the-target)
- [Offset cards by observing template](#offset-cards-by-observing-template)
- [What ESO does for X-Shooter internal data products](#what-eso-does-for-x-shooter-internal-data-products)
- [Edge cases and open points](#edge-cases-and-open-points)
- [What soxspipe does today](#what-soxspipe-does-today)
- [Sources](#sources)

## What the standard asks for

Section 5.1 of the ESO Science Data Products Standard defines `RA` and `DEC` as "Equatorial coordinates in decimal degrees (J2000)", and for spectra as the "spectroscopic target position" with footnote 10. Footnote 10 reads: "Estimates like the center of the slit or the telescope pointing positions shall not be used." `EQUINOX` is mandatory and must be 2000 for FK5, and `RADESYS` is ICRS or FK5 [S1]. The standard does not name an epoch for the position, only the equinox.

## Survey method

- Raw data: the definitive SOXS raw workspace `/Volumes/SOXS/raw_archive/soxs_commissioning_data`. The `raw_frames` table of its `soxspipe.db` was opened read-only to list every frame with `eso dpr catg = 'SCIENCE'` and `eso seq arm` VIS or NIR [S2].
- That gave 2213 frames from 115 nights (2025-05-06 to 2026-10-12) and 262 `OBJECT` values: 1468 stare, 678 nod (`AutoNodOnSlit`), 50 generic offset and 16 fixed sky offset frames, plus 1 synchronous stare frame.
- Only the primary header of each frame was read, by streaming the first 576000 bytes of the decompressed `.fits.Z` file into `astropy.io.fits.Header.fromstring`. No pixel data was read and no reduction was run.
- Two frames (`SOXS.2025-12-02T06:17:01.195` NIR and `SOXS.2025-12-02T06:17:03.063` VIS, target HV 11471) have no `RA`, `DEC` or `TEL TARG` cards at all. The statistics below use the other 2211 frames.
- Angular separations were computed with `astropy.coordinates.SkyCoord.separation`. Offsets are given as (ΔRA·cos δ, ΔDec) in arcsec, pointing minus target.
- SIMBAD positions came from the CDS Sesame name resolver (`nph-sesame/-oI/SN`), which returns ICRS J2000 coordinates and Gaia DR3 proper motions [S6].

## Candidate cards and their format

Every one of the 2211 frames carries the full `TEL TARG` block. The format column comes from the VLT TCS keyword dictionary `ESO-VLT-DIC.TCS` [S3].

| Card | Dictionary format | Values seen on SOXS frames | Meaning |
|---|---|---|---|
| `HIERARCH ESO TEL TARG ALPHA` | `%.3f`, unit `HHMMSS.TTT`, range 0 to 240000 | for example `24240.704` (M 77), `161955.069` (V818 Sco) | Target right ascension from the OB |
| `HIERARCH ESO TEL TARG DELTA` | `%.3f`, unit `DDMMSS.TTT`, range -900000 to 900000 | for example `-47.83` (M 77, that is -00:00:47.83), `155313.228` (GD 71) | Target declination from the OB |
| `HIERARCH ESO TEL TARG COORDTYPE` | `M` (mean) or `A` (apparent) | `M` on all 2211 frames | Mean place |
| `HIERARCH ESO TEL TARG EQUINOX` | years | `2000.0` on all frames | Equinox |
| `HIERARCH ESO TEL TARG EPOCH`, `EPOCHSYSTEM` | 1950 or 2000; `J` or `B` | `2000.0`, `J` on all frames | Epoch of the position |
| `HIERARCH ESO TEL TARG PMA`, `PMD` | `%.6f`, arcsec/yr, range -10 to 10 | nonzero on 817 frames (stars); zero for transients and galaxies | Proper motion from the OB |
| `HIERARCH ESO TEL TARG PARALLAX` | arcsec | `0.0` on all frames | Parallax |
| `HIERARCH ESO TEL TARG RADVEL` | km/s | `0.0` on the frames inspected | Radial velocity |
| `RA`, `DEC` | decimal degrees | comment `[deg] hh:mm:ss.s RA (J2000) pointing` on all frames; 6 and 5 decimals | Telescope pointing |
| `EQUINOX`, `RADECSYS` | | `2000.0`, `FK5` on all frames | Frame of the pointing cards |
| `HIERARCH ESO ADA GUID RA`, `DEC` | decimal degrees | present on 2155 frames | Guide star, not the target |
| `HIERARCH ESO INS ADC1 RA`, `DEC` (and `ADC2`) | packed sexagesimal | for example `24402.235` for M 77 | Telescope position given to the ADC, of date (not J2000); not the target |
| `HIERARCH ESO OBS TARG NAME` | string | present on 2189 frames | Target name only, no coordinates |

Decoding rule: for `ALPHA = a`, RA in degrees is `15 × (hh + mm/60 + ss/3600)` with `hh = floor(a/10000)`, `mm = floor((a − 10000·hh)/100)`, `ss = a − 10000·hh − 100·mm`. For `DELTA = d`, take the sign of `d` and decode `|d|` the same way without the factor 15. The sign must come from the number itself, because a declination between -1 deg and 0 deg is a small negative number such as `-47.83`. The same packed format is used for X-Shooter (`101421.290`, `-424830.900` in [S5]).

The precision of the cards is 0.001 s in RA (0.015 arcsec at the equator) and 0.01 arcsec or better in declination. This is far below the accuracy needed for `RA`/`DEC`.

No SOXS science frame has a card for blind offsets, acquisition offsets, or non-sidereal tracking rates. The acquisition frames of the two worst outliers (below) have no such card either [S2].

## Card values against pointing on sampled frames

Target = decoded `TEL TARG ALPHA`/`DELTA`. Pointing = `RA`/`DEC`. Offset is pointing minus target, in arcsec, as (ΔRA·cos δ, ΔDec). PM is `TEL TARG PMA`, `PMD` in arcsec/yr. Frame names drop the `.fits.Z` suffix; all are under `raw/<night start date>/` [S2].

| Frame | Arm | Template | DPR TYPE | OBJECT | TARG ALPHA, DELTA (card) | Target (deg) | Pointing RA, DEC (deg) | PM | Offset cards | Offset (arcsec) | Separation (arcsec) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| `SOXS.2025-06-02T23:55:25.034` | VIS | stare | OBJECT | SN2025kzr | 102936.43, -240633.61 | 157.401792, -24.109336 | 157.402294, -24.109610 | 0, 0 | none | +1.7, -1.0 | 1.9 |
| `SOXS.2025-08-07T03:35:33.462` | NIR | stare | OBJECT | LTT7379 | 183625.95, -441836.907 | 279.108125, -44.310252 | 279.104521, -44.311730 | -0.170, -0.160 | none | -9.3, -5.3 | 10.7 |
| `SOXS.2025-06-03T03:14:51.403` | VIS | stare | OBJECT | LAWD 37 | 114542.916, -645029.461 | 176.428817, -64.841517 | 176.472271, -64.843900 | +2.662, -0.345 | none | +66.5, -8.6 | 67.1 |
| `SOXS.2026-01-18T02:47:05.965` | NIR | stare | OBJECT | NEO500080 | 61901.1, 45208.0 | 94.754583, +4.868889 | 94.752433, +4.867060 | 0, 0 | none | -7.7, -6.6 | 10.1 |
| `SOXS.2026-09-26T07:05:03.778` | VIS | stare | OBJECT | SN2025pht | 44128.87, -25155.66 | 70.370292, -2.865461 | 70.386721, -2.889130 | 0, 0 | none | +59.1, -85.2 | 103.7 |
| `SOXS.2025-05-13T04:19:57.654` | VIS | nod | OBJECT | V V818 Sco | 161955.069, -153825.017 | 244.979454, -15.640282 | 244.978396, -15.640110 | -0.007, -0.012 | CUMOFF RA, DEC +3.2, +0.0 | -3.7, +0.6 | 3.7 |
| `SOXS.2025-05-13T04:23:08.106` | NIR | nod | OBJECT | V V818 Sco | 161955.069, -153825.017 | 244.979454, -15.640282 | 244.976316, -15.640110 | -0.007, -0.012 | CUMOFF RA, DEC -4.0, +0.0 | -10.9, +0.6 | 10.9 |
| `SOXS.2025-05-13T04:08:15.092` | VIS | offset (fixed sky) | OBJECT | V V818 Sco | 161955.069, -153825.017 | 244.979454, -15.640282 | 244.977593, -15.640110 | -0.007, -0.012 | FIXOFF RA, DEC +0.0, -0.4 | -6.5, +0.6 | 6.5 |
| `SOXS.2025-05-13T04:10:26.213` | NIR | offset (fixed sky) | OBJECT | V V818 Sco | 161955.069, -153825.017 | 244.979454, -15.640282 | 244.983376, -15.634460 | -0.007, -0.012 | FIXOFF RA, DEC +20.0, +20.4 | +13.6, +21.0 | 25.0 |
| `SOXS.2025-11-30T05:27:27.385` | NIR | offset (fixed sky) | OBJECT | M 77 | 24240.704, -47.83 | 40.669600, -0.013286 | 40.671201, -0.015400 | 0, 0 | FIXOFF RA, DEC +0.0, -0.4 | +5.8, -7.6 | 9.5 |
| `SOXS.2025-11-30T05:32:02.820` | NIR | offset (fixed sky) | SKY | M 77 | 24240.704, -47.83 | 40.669600, -0.013286 | 40.670636, -0.009730 | 0, 0 | FIXOFF RA, DEC -2.0, +20.4 | +3.7, +12.8 | 13.3 |
| `SOXS.2025-10-25T07:20:41.198` | VIS | offset (generic) | OBJECT | GD 71 | 55227.619, 155313.228 | 88.115079, +15.887008 | 88.115883, +15.882490 | +0.077, -0.173 | none | +2.8, -16.3 | 16.5 |
| `SOXS.2025-10-25T07:26:14.359` | VIS | offset (generic) | SKY | GD 71 | 55227.619, 155313.228 | 88.115079, +15.887008 | 88.116459, +15.883040 | +0.077, -0.173 | none | +4.8, -14.3 | 15.1 |
| `SOXS.2026-01-24T05:33:26.555` | NIR | offset (generic) | OBJECT | EG21 | 31031.019, -683603.38 | 47.629246, -68.600939 | 47.626162, -68.602270 | +0.040, -0.103 | none | -4.1, -4.8 | 6.3 |
| `SOXS.2026-01-24T06:50:03.625` | NIR | offset (generic) | OBJECT | CD-32 5613 | 84132.429, -325632.915 | 130.385121, -32.942476 | 130.377130, -32.934550 | -1.061, +1.346 | none | -24.1, +28.5 | 37.4 |

Separation of pointing from target over all 2211 frames, by template and `DPR TYPE` [S2]:

| Template | DPR TYPE | Frames | Median (arcsec) | 90th percentile (arcsec) | Max (arcsec) |
|---|---|---|---|---|---|
| Stare (`SOXS_slt_obs_Stare`) | OBJECT | 1466 | 4.5 | 12.5 | 103.7 |
| Stare, synchronous | OBJECT | 1 | 1.7 | 1.7 | 1.7 |
| Nod (`SOXS_slt_obs_AutoNodOnSlit`) | OBJECT | 678 | 8.0 | 29.0 | 49.9 |
| Fixed sky offset (`SOXS_slt_obs_FixedSkyOffset`) | OBJECT | 14 | 10.4 | 33.7 | 33.8 |
| Fixed sky offset | SKY | 2 | 13.3 | 13.3 | 13.3 |
| Generic offset (`SOXS_slt_obs_GenericOffset`) | OBJECT | 39 | 37.2 | 42.4 | 44.4 |
| Generic offset | SKY | 11 | 15.1 | 15.4 | 36.7 |

## How accurate the target cards are

The decoded `TEL TARG` position was compared with the SIMBAD J2000 position for 14 targets [S6]. Each target had one `TEL TARG` value across all its frames, except M 77, which had two values 0.08 arcsec apart (two OBs).

| OBJECT | Frames | TEL TARG minus SIMBAD (arcsec) | Median pointing minus TEL TARG (arcsec) | Note |
|---|---|---|---|---|
| GD 71 | 69 | 0.01 | 12.6 | |
| Feige 110 | 90 | 0.01 | 4.3 | |
| EG21 | 34 | 0.00 | 7.0 | |
| LTT7379 | 7 | 0.00 | 10.7 | |
| M 77 | 18 | 0.00 and 0.08 | 9.6 | two OBs |
| CD-32 5613 | 108 | 0.00 | 38.4 | PM 1.7 arcsec/yr |
| 108 Vir | 28 | 0.01 | 4.7 | |
| V V818 Sco | 27 | 0.01 | 6.5 | |
| LTT7987 | 21 | 0.01 | 16.1 | |
| SN 2025pht | 18 | 0.00 | 7.2 | transient |
| PSR J1023+0038 | 6 | 0.00 | 7.8 | |
| G 24-9 | 8 | 0.01 | 19.5 | |
| CD-32d9927 | 6 | 0.60 | 2.6 | OB has the PM of LAWD 37 (+2.6616, -0.34493); SIMBAD gives -0.003, +0.007 |
| HR5501 | 6 | 0.86 | 5.2 | OB has the PM of LTT7379 (-0.16979, -0.15981) |

Findings:

- `TEL TARG` reproduces the catalogue position that the observer typed in. For the stars it is the Gaia-based SIMBAD position at epoch J2000 (it is labelled FK5; the FK5 to ICRS difference is below about 0.05 arcsec, which does not matter here).
- `TEL TARG PMA` equals the SIMBAD `pmRA` (which is μα·cos δ) divided by 1000 for every star checked, for example GD 71 `0.076728` against 76.728 mas/yr.
- The accuracy of `TEL TARG` is the accuracy of the OB entry. It cannot be checked from the header, and two of 14 targets had copied-in errors. An OB with a wrong position is invisible to the pipeline.

## Why the pointing differs from the target

1. **Acquisition and pointing-model residual.** After the acquisition template puts the target on the slit, the pointing still differs from the target by a few arcsec. On the V818 Sco nod pair, pointing minus `CUMOFF` gives (-6.9, +0.6) arcsec at both nod positions, and the fixed sky offset frames of the same star 10 minutes earlier give (-6.5, +0.6) arcsec. On stare frames with no proper motion the offset is (-0.6, -2.6) arcsec median, with a 5 to 95 percent range of about ±6 arcsec in RA and -8 arcsec to +4 arcsec in Dec, and the monthly median moves by several arcsec [S2].
2. **Proper motion.** The telescope applies `TEL TARG PMA`/`PMD` to the J2000 position, so `RA`/`DEC` follow the star. For LAWD 37 (PM 2.66 arcsec/yr) the pointing is 67 arcsec from `TEL TARG` in 2025; after propagating the proper motion to `MJD-OBS` the residual is 1.2 arcsec. On the 408 stare frames of stars with a nonzero proper motion, propagating it cuts the 95th percentile of the separation from 38 arcsec to 9 arcsec [S2].
3. **Nod, jitter and sky offsets.** These move the pointing by `SEQ CUMOFF` or `SEQ FIXOFF`, or by unrecorded amounts in the generic offset template (see next section). `TEL TARG` stays fixed.
4. **Unrecorded moves during acquisition.** For SN2025mqc (2025-06-09, 71 arcsec; the 2025-06-11 frames show the same 69 arcsec offset) and SN2025pht (2026-09-26, 104 arcsec), the earlier acquisition frames of the 2025-06-09 and 2026-09-26 visits point at `TEL TARG`, and the last acquisition frame jumps away from it. No header card records why. In these cases neither `TEL TARG` nor the pointing is a proven position of what was on the slit [S2].

## Offset cards by observing template

| Template | Cards present | Frames | Behavior seen |
|---|---|---|---|
| `SOXS_slt_obs_Stare` | none | 1468 | Pointing fixed during the template |
| `SOXS_slt_obs_AutoNodOnSlit` | `SEQ CUMOFF RA/DEC/X/Y`, `SEQ RELOFF RA/DEC`, `SEQ NOD THROW`, `SEQ JITTER WIDTH` | 678 | Pointing = target + acquisition residual + `CUMOFF RA/DEC` (arcsec on sky) |
| `SOXS_slt_obs_FixedSkyOffset` | `SEQ FIXOFF RA/DEC`, `SEQ JITTER WIDTH` | 16 | `FIXOFF` is the step from the previous exposure, not the total offset: CD-44 12736 steps (0, -0.8), (+30, +30.8), (0, +0.1), (-30, -30.1) arcsec |
| `SOXS_slt_obs_GenericOffset` | none | 50 | Pointing moves between exposures but no card records the offset; `DPR TYPE` is the only object or sky flag |

For Phase 3 this does not matter: every on-target and sky frame in a product has the same `TEL TARG`, so the product `RA`/`DEC` does not depend on the offsets.

## What ESO does for X-Shooter internal data products

The XSHOOTER_ECHELLE release description (v1.4) and the PHOENIX processing page do not say how `RA`/`DEC` are set [S7][S8]. So eight X-Shooter internal data products and their `PROV1` raw files were read from the ESO archive header service [S4][S5]:

| ADP | PRO CATG | DPR TECH | ADP RA, DEC (deg) | PROV1 raw pointing RA, DEC (deg) | ADP minus raw pointing (arcsec) | ADP minus raw TEL TARG (arcsec) |
|---|---|---|---|---|---|---|
| `ADP.2024-10-24T13:23:40.647` | `SCI_SLIT_FLUX_IDP_UVB` | STARE | 266.211617, -24.59633 | 266.211617, -24.59633 | 0.000 | 2.54 |
| `ADP.2024-10-24T13:23:40.716` | `SCI_SLIT_FLUX_IDP_UVB` | STARE | 1.554877, -39.24183 | 1.554877, -39.24183 | 0.000 | 5.89 |
| `ADP.2024-10-24T13:23:40.686` | `SCI_SLIT_FLUX_IDP_UVB` | STARE | 15.185919, -56.77299 | 15.185919, -56.77299 | 0.000 | 11.30 |
| `ADP.2024-10-24T13:23:40.776` | `SCI_SLIT_FLUX_IDP_UVB` | NODDING | 3.606011, -30.38052 | 3.606011, -30.38052 | 0.000 | 2.14 |
| `ADP.2024-10-24T13:23:40.660` | `SCI_SLIT_FLUX_IDP_UVB` | STARE | 71.721433, 16.99981 | 71.721433, 16.99981 | 0.000 | 1.49 |
| `ADP.2024-10-24T14:34:50.728` | `SCI_SLIT_FLUX_IDP_UVB` | NODDING | 40.353935, -36.89031 | 40.353935, -36.89031 | 0.000 | 3.88 |
| `ADP.2024-10-30T16:03:52.164` | `SCI_SLIT_FLUX_IDP_UVB` | STARE | 122.475745, -47.36162 | 122.475745, -47.36162 | 0.000 | 0.88 |
| `ADP.2023-03-03T09:05:21.303` | `TELL_SLIT_FLUX_IDP_UVB` | NODDING | 153.587563, -42.80826 | 153.587563, -42.80826 | 0.000 | 3.24 |

In every case the product `RA`/`DEC` equal the pointing of the first raw file to the last digit, including the nod offset of the first nod frame (for example `SEQ CUMOFF RA, DEC = -0.644, +2.416` arcsec on the last row). Only the comment changes, from `RA (J2000) pointing` to `[deg] Spectroscopic target position (J2000)`. The products carry `PROCSOFT = 'xshoo/3.6.0'` and still have `HIERARCH ESO TEL TARG ALPHA`/`DELTA` in the header [S5].

So ESO's in-house X-Shooter products do not follow footnote 10 of the standard they validate against. This means:

- There is no ESO precedent for using `TEL TARG`. The standard [S1] is the only authority.
- The archive accepts pointing-based `RA`/`DEC` for X-Shooter, so the Phase 3 validator probably does not check them against the target. This is inferred, not checked.

## Edge cases and open points

- **Epoch for high proper motion stars (for Dave).** `TEL TARG` is the J2000 position at epoch J2000. For LAWD 37 that is 67 arcsec from the star in 2025, and 38 arcsec for CD-32 5613. The standard says "J2000" but does not name an epoch [S1]. Option 1: write `TEL TARG` unchanged (catalogue position, the same value across epochs, as the standard asks for `OBJECT` in footnote 8). Option 2: propagate `PMA`/`PMD` to `MJD-OBS` (the sky position of the target at the time of the observation; the pointing cards show the telescope does this). Option 2 depends on correct proper motions in the OB, which failed for 2 of 14 targets checked.
- **Missing cards.** Two frames (HV 11471, 2025-12-02) have no `TEL TARG`, `RA` or `DEC`. The Phase 3 writer needs a rule for this, for example to fail the Phase 3 export with a clear message.
- **Moving targets.** NEO500080 has a different `TEL TARG` on each night and no tracking-rate cards. `TEL TARG` is then the ephemeris position the observer entered for that night, which is the best value in the header.
- **Unrecorded acquisition moves.** In the SN2025mqc and SN2025pht cases, the telescope moved 70 arcsec to 104 arcsec away from `TEL TARG` during acquisition. If the observer moved to a different source, `TEL TARG` is wrong for that product, and nothing in the header says so.
- **Extended targets.** For M 77 `TEL TARG` is the nucleus position from the OB, not the slit position on the galaxy. The standard asks for the target position, so this is acceptable.
- **Frame label.** The `TEL TARG` values are Gaia-based ICRS positions in practice, but every header says `RADECSYS = 'FK5'`, `EQUINOX = 2000`. Keeping `RADESYS = 'FK5'` and `EQUINOX = 2000.0` satisfies section 5.1 [S1].

## What soxspipe does today

- `soxspipe/resources/soxs_keywords.yaml` maps the `RA` and `DEC` aliases to the pointing cards `RA` and `DEC`. There is no alias for `TEL TARG ALPHA`/`DELTA`.
- `basic_header_scrubbing` in `soxspipe/commonutils/phase3.py` renames `RADECSYS` to `RADESYS` and `HIERARCH ESO TEL TARG EQUINOX` to `EQUINOX`, but does not touch `RA`/`DEC`. So a product today keeps the pointing of the raw frame its header came from, which footnote 10 forbids.
- Side finding for the offset work: `soxspipe/recipes/soxs_offset.py` splits on and off frames by `kw("OFFSET_RA")` and `kw("OFFSET_DEC")`, which map to `HIERARCH ESO SEQ FIXOFF RA/DEC`. Those cards exist only in the fixed sky offset template (16 frames). The 50 generic offset frames have no offset cards, and `FIXOFF` is a step from the previous exposure, not a position. This was not tested against the recipe.

## Sources

- [S1] ESO Science Data Products Standard, section 5.1 table (`RA`, `DEC`, `EQUINOX`, `RADESYS`) and footnotes 8 to 10. Local Markdown copy: `P32.11_eso_phase_3/output/p3sdpstd_md/p3sdpstd.md`, lines 518 to 538. Original: <https://www.eso.org/sci/observing/phase3/p3sdpstd.pdf>.
- [S2] SOXS raw science frame headers, `/Volumes/SOXS/raw_archive/soxs_commissioning_data/raw/<night start date>/SOXS.*.fits.Z`, listed from the `raw_frames` table of `soxspipe.db` in the same folder (opened read-only). Survey run on 2026-11-19: 2213 frames, primary headers only.
- [S3] VLT TCS keyword dictionary `ESO-VLT-DIC.TCS` (version 1.38, 2004), entries `TEL TARG ALPHA`, `DELTA`, `EPOCH`, `EPOCHSYSTEM`, `EQUINOX`, `PMA`, `PMD`, `PARALLAX`, `COORDTYPE`. Copy at <https://research.ast.cam.ac.uk/vdfs/internal/vista/WP2/modules/dicTCS/config/ESO-VLT-DIC.TCS>. The SOXS frames give `ESO-VLT-DIC.TCS` as their TEL dictionary (`HIERARCH ESO TEL DID`).
- [S4] ESO archive ObsCore TAP service, <https://archive.eso.org/tap_obs>, table `ivoa.ObsCore`, query `obs_collection = 'XSHOOTER' AND dataproduct_type = 'spectrum'`.
- [S5] ESO archive header service, `https://archive.eso.org/hdr?DpId=<id>`, for the eight ADPs in the table above and their `PROV1` raw files (for example `XSHOO.2023-02-25T04:20:05.725`).
- [S6] CDS Sesame name resolver (SIMBAD), mirror <https://vizier.idia.ac.za/viz-bin/nph-sesame/-oI/SN?>, queried 2026-11-19 for each target name in the accuracy table.
- [S7] ESO Phase 3 release description XSHOOTER_ECHELLE v1.4, <https://www.eso.org/sci/observing/phase3/data_releases/xshooter_echelle_v1.4.pdf>. It gives no rule for `RA`/`DEC`.
- [S8] ESO PHOENIX X-Shooter processing page, <https://eso.org/qc/PHOENIX/XSHOOTER/processing.html>. It gives no rule for `RA`/`DEC`.
