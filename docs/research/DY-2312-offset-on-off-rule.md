# DY-2312: how reliable the soxs-offset ON/OFF frame rule is

These notes answer Linear ticket DY-2312 for the ESO Phase 3 map (DY-2291). They are for whoever changes how `soxs-offset` splits its science frames into on-target (ON) and sky (OFF) frames, and for whoever computes the Phase 3 `EXPTIME`, `NCOMBINE`, `MJD-OBS`, `MJD-END` and `PROVi` keywords from the ON list.

## Answer

The current rule, ON when `HIERARCH ESO SEQ FIXOFF RA + DEC < 0`, gives the right answer for all 16 frames it can see in the archive, but half of those answers are luck. `SEQ FIXOFF RA/DEC` is the telescope step from the previous exposure (including a random jitter), not the position. So a frame taken at the same position as the one before it has a step of a few tenths of an arcsec, and its sign is the sign of the jitter. Eight of the 16 frames are such jitter-only steps. The rule also fails for any sky position whose RA and DEC offsets add up to a negative number, which the template allows (each offset may be -100 to +100 arcsec). The rule cannot work at all for the second spectroscopic offset template, `SOXS_slt_obs_GenericOffset` (50 frames), because those frames carry no offset cards; the recipe would raise a `KeyError`. Today this does not happen only because the data organiser never routes those frames to the recipe (`DPR TECH = ECHELLE,SLIT,GENOFFSET` is not in `soxs_sof_map.yaml`). A third, separate problem: the data organiser sets `ignore = 1` on every frame with `DPR TYPE = SKY`, so for the newer offset observations the OFF frames never reach the recipe, and the M 77 offset SOF holds one ON frame and no OFF frame.

Recommended rule:

1. For `DPR TECH = ECHELLE,SLIT,OFFSET` (`SOXS_slt_obs_FixedSkyOffset`, and `SOXS_slt_cal_SpecphotStdOffset`, which has the same parameters), classify each frame by its position, not by its step. Take the frames of one template execution (one `TPL START`), and measure the angular distance from each frame's pointing (`RA`, `DEC`) to the pointing of the template's first exposure (lowest `TPL EXPNO`). A frame within a small tolerance of that first position (the jitter box width `SEQ JITTER WIDTH` plus 1 arcsec) is ON; a frame farther away is OFF. The template always starts on the object (object-sky-sky-object cycles, per the SOXS user manual). As a cross-check, the running sum of `SEQ FIXOFF RA/DEC` in exposure order per arm gives the same positions to within 0.7 arcsec on every frame found; stop with an error if the two disagree. This rule classifies all 16 FixedSkyOffset frames correctly, as confirmed by the trace signal in the raw frames.
2. For `DPR TECH = ECHELLE,SLIT,GENOFFSET` (`SOXS_slt_obs_GenericOffset`), `DPR TYPE` is the only per-frame flag. It comes from the observer's list `SEQ.OBS.TYPE` (O or S for each offset). Use `OBJECT` as ON and `SKY` as OFF, but only accept an `OBJECT` frame as ON when its pointing is within 1 arcsec (about one slit width) of the template's first `OBJECT` frame. Refuse, with a clear message, an `OBJECT` frame that the telescope moved away from that position, and refuse a set with no `SKY` frame. This rule is correct for the 10 frames of the three well-formed GenericOffset observations (SN2025kzr and two M 77 visits). It cannot be made correct from the headers for the GD 71 and CD-32 5613 test observations of October 2025, where the labels and the light on the slit disagree.
3. Do not use `DPR TYPE` to split FixedSkyOffset frames. The May 2025 frames are all `OBJECT`, including the sky position; only the November 2025 frames mark the sky position `SKY`.
4. Stop the data organiser from ignoring `DPR TYPE = SKY` frames whose `DPR TECH` is an offset technique, and add the `GENOFFSET` technique to the SOF map if GenericOffset science is to be reduced.

Frames the recommended rule cannot classify: 24 frames, all from two October 2025 GenericOffset test OBs (GD 71, OB 1000603202, 12 frames; CD-32 5613, OB 1000603272, 12 frames). Their headers say one thing and the slit holds another (details below). In addition, the two January 2026 GenericOffset OBs (EG21 and CD-32 5613, 16 frames) can be classified (all ON) but cannot be reduced as offset data, because the telescope never moved off the target and no sky frame exists.

**Confidence: high** for the FixedSkyOffset findings (every frame checked against both the pointing and the trace signal). **Medium** for GenericOffset: the rule rests on the observer typing correct O/S labels, and only three of the six GenericOffset observations used the template as intended.

## Contents

- [The current rule](#the-current-rule)
- [Survey method](#survey-method)
- [What the templates do, per the user manual](#what-the-templates-do-per-the-user-manual)
- [Offset frames in the archive](#offset-frames-in-the-archive)
- [FixedSkyOffset frame by frame](#fixedskyoffset-frame-by-frame)
- [GenericOffset frame by frame](#genericoffset-frame-by-frame)
- [Why the current rule is fragile](#why-the-current-rule-is-fragile)
- [How the data organiser filters offset frames](#how-the-data-organiser-filters-offset-frames)
- [Recommended rule and what it cannot classify](#recommended-rule-and-what-it-cannot-classify)
- [Sources](#sources)

## The current rule

`soxs_offset._split_on_off_frames` reads `kw("OFFSET_RA")` and `kw("OFFSET_DEC")`, which `soxs_keywords.yaml` maps to `HIERARCH ESO SEQ FIXOFF RA` and `HIERARCH ESO SEQ FIXOFF DEC`. A frame is ON when their sum is negative, and OFF otherwise. `_check_on_off_balance` then refuses a set with unequal ON and OFF counts, or with no ON frame [S1].

## Survey method

- The raw-frame table of the definitive SOXS workspace, `/Volumes/SOXS/raw_archive/soxs_commissioning_data/soxspipe.db`, was opened read-only and searched for any frame whose `eso dpr tech` or template name mentions an offset, or whose `eso dpr type` mentions sky [S2].
- That gave 66 spectroscopic offset frames: 16 with `DPR TECH = ECHELLE,SLIT,OFFSET` and 50 with `ECHELLE,SLIT,GENOFFSET`. The other offset-like frames are imaging (`IMAGE,OFFSET,POINT`, 86 frames) and nod frames (`ECHELLE,SLIT,NODDING`), which `soxs-offset` never sees.
- The primary header of all 66 frames was read in full (by streaming the decompressed `.fits.Z` file into `astropy.io.fits.Header.fromstring`). Every card containing `SEQ`, `TPL`, `OFF`, `TEL`, `DPR`, `OBS` or `SKY` was listed.
- Truth was taken from the light on the slit. Each frame was decompressed one at a time. For VIS, the median of rows 1800 to 2200 gives a profile across the four orders; for NIR, the median of columns 950 to 1100 does the same. The profile minus a 41-pixel running median keeps narrow point-source peaks and removes slit-filling light (sky, galaxy halo). The highest remaining peak, divided by the exposure time, is the trace signal in the tables below (counts per second, ADU for VIS, detector units for NIR). Two frames were also inspected by eye as images. No reduction was run.
- Pointing offsets are (ΔRA·cos δ, ΔDec) in arcsec of the `RA`/`DEC` pointing cards from the decoded `TEL TARG ALPHA`/`DELTA` target position (see DY-2296).

## What the templates do, per the user manual

The SOXS User Manual (SOXS-MAN-0003, issue 2.1, 2026-07-07) describes two spectroscopic offset templates [S3]:

- `SOXS_slt_obs_FixedSkyOffset` (section 4.2.6 and the template table on page 42): "alternating between an object (O) and sky position (S), with the possibility of adding a small jittering around the object and the sky position. One cycle is a pair of OS or SO observations. Cycles are repeated in OSSO sequences." Parameters include `SEQ.FIXOFF.RA`, `SEQ.FIXOFF.DEC` (the sky offset, -100 to 100 arcsec), `SEQ.JITTER.WIDTH` (0 to 2 arcsec) and `SEQ.NCYCLE`. The fixed values list `DPR.TYPE = OBJECT` and `DPR.TECH = ECHELLE,SLIT,OFFSET`. The calibration template `SOXS_slt_cal_SpecphotStdOffset` (page 47) has the same offset parameters with `DPR.TYPE = STD,FLUX`.
- `SOXS_slt_obs_GenericOffset` (section 4.2.7 and page 43): "define any pattern by providing a list of (cumulative) telescope offsets". Parameters include `SEQ.OBS.TYPE` ("List of TYPE of images at each offsets", O or S, default `O S O S`), `SEQ.OFFSET.COORDS` (SKY or SLIT), and `SEQ.RELOFF1`/`SEQ.RELOFF2` (the offset lists). `DPR.TECH = ECHELLE,SLIT,GENOFFSET`.
- Section 6.1: "Templates use cumulative offsets: the position at a given time is derived from the sum of all offsets specified so far in the template." This matches what the headers show: `SEQ FIXOFF` in a frame header is the step to that frame's position, not the position.

The SOXS Template Manual (SOXS-MAN-0001), which the user manual cites, was not found online.

## Offset frames in the archive

| Template (`TPL ID`) | `DPR TECH` | OB | OBS NAME | OBJECT | Night | Template runs | Frames per run (VIS + NIR) | `DPR TYPE` values |
|---|---|---|---|---|---|---|---|---|
| `SOXS_slt_obs_FixedSkyOffset` | `ECHELLE,SLIT,OFFSET` | 1000542321 | SOXS_ScoX-1_offset_test | V V818 Sco | 2025-05-12 | 1 | 2 + 2 | all OBJECT |
| `SOXS_slt_obs_FixedSkyOffset` | `ECHELLE,SLIT,OFFSET` | 1000542880 | SOXS_LTT7379-STD_FixedSkyOffset | CD-44 12736 | 2025-05-13 | 1 | 4 + 4 | all OBJECT |
| `SOXS_slt_obs_FixedSkyOffset` | `ECHELLE,SLIT,OFFSET` | 1000605643 | SOXS_M77_test_OB_slt_FixedSkyOffset | M 77 | 2025-11-29 | 1 | 2 + 2 | OBJECT, SKY |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000550097 | SOXS_SN2025kzr_GenericOffset | SN2025kzr | 2025-06-05 | 1 | 2 + 0 | OBJECT, SKY |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000603202 | SOXS_GD71_General_Offset | GD 71 | 2025-10-24 | 3 | 4 + 0 | OBJECT, SKY |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000603272 | SOXS_LTT3218_Generic_Offset | CD-32 5613 | 2025-10-25 | 3 | 4 + 0 | all OBJECT |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000605678 | SOXS_M77_test_OB_slt_FixedSkyOffset | M 77 | 2025-11-30 | 1 | 2 + 2 | OBJECT, SKY |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000606271 | SOXS_M77_test_OB_slt_FixedSkyOffset | M 77 | 2025-12-09 | 1 | 2 + 2 | OBJECT, SKY |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000609277 | EG21-STD_stare_15 | EG21 | 2026-01-23 | 1 | 4 + 4 | all OBJECT |
| `SOXS_slt_obs_GenericOffset` | `ECHELLE,SLIT,GENOFFSET` | 1000609325 | LTT3218-STD_stare_15 | CD-32 5613 | 2026-01-23 | 1 | 4 + 4 | all OBJECT |

`TPL NAME` is "Alternating Object and Sky position" for FixedSkyOffset and "Define patterns via list of offsets" for GenericOffset. `TPL VERSION` is `$Revision: 401845 $` (FixedSkyOffset) and `$Revision: 407031 $` (GenericOffset) on every frame, and `TPL PRESEQ` is `SOXS_slt_obs_AllNodOffset.seq` on all 66 frames. No `SpecphotStdOffset` frame exists in the archive.

Offset-related cards present [S2]:

| Card | FixedSkyOffset (16 frames) | GenericOffset (50 frames) |
|---|---|---|
| `SEQ FIXOFF RA`, `SEQ FIXOFF DEC` | yes | no |
| `SEQ JITTER WIDTH` | yes (1.0 or 2.0) | no |
| `SEQ CUMOFF *`, `SEQ RELOFF *`, `SEQ OBS TYPE`, `SEQ NCYCLE`, `TEL OFFS *` | no | no |
| `TEL TARG ALPHA`, `DELTA`, `PMA`, `PMD` | yes, constant within an OB | yes, constant within an OB |
| `RA`, `DEC` (pointing) | yes | yes |
| `TPL EXPNO`, `TPL NEXP`, `TPL START` | yes | yes |
| `DPR TYPE` | OBJECT, or SKY (November 2025 only) | OBJECT or SKY, from `SEQ.OBS.TYPE` |

`TPL EXPNO` counts the exposures of both arms together: in a FixedSkyOffset run with one exposure per arm per position, the VIS and NIR frames at one position have consecutive numbers (for example NIR 1, VIS 2 at the object), and the two arms at one position share the same pointing and the same `FIXOFF` step.

## FixedSkyOffset frame by frame

Step is `SEQ FIXOFF RA, DEC`. Position is the running sum of the steps for that arm. Pointing is the `RA`/`DEC` offset from `TEL TARG`. Rule is what the current sign test decides. Trace is the trace signal (counts per second). Truth comes from the trace.

| Frame | Arm | EXPNO/NEXP | DPR TYPE | Step (arcsec) | Position (arcsec) | Pointing (arcsec) | Rule | Trace | Truth |
|---|---|---|---|---|---|---|---|---|---|
| `SOXS.2025-05-13T04:08:15.092` | VIS | 2/4 | OBJECT | 0.0, -0.4 | 0.0, -0.4 | -6.5, +0.6 | ON | 2.55 | ON |
| `SOXS.2025-05-13T04:08:50.132` | NIR | 1/4 | OBJECT | 0.0, -0.4 | 0.0, -0.4 | -6.5, +0.6 | ON | 1.60 | ON |
| `SOXS.2025-05-13T04:09:51.115` | VIS | 4/4 | OBJECT | +20.0, +20.4 | +20.0, +20.0 | +13.6, +21.0 | OFF | 0.08 | OFF |
| `SOXS.2025-05-13T04:10:26.213` | NIR | 3/4 | OBJECT | +20.0, +20.4 | +20.0, +20.0 | +13.6, +21.0 | OFF | 0.20 | OFF |
| `SOXS.2025-05-14T04:02:38.213` | VIS | 2/8 | OBJECT | 0.0, -0.8 | 0.0, -0.8 | -9.6, -3.9 | ON | 7.81 | ON |
| `SOXS.2025-05-14T04:03:11.831` | NIR | 1/8 | OBJECT | 0.0, -0.8 | 0.0, -0.8 | -9.6, -3.9 | ON | 3.91 | ON |
| `SOXS.2025-05-14T04:08:13.629` | VIS | 4/8 | OBJECT | +30.0, +30.8 | +30.0, +30.0 | +20.5, +26.8 | OFF | 0.03 | OFF |
| `SOXS.2025-05-14T04:08:47.248` | NIR | 3/8 | OBJECT | +30.0, +30.8 | +30.0, +30.0 | +20.5, +26.8 | OFF | 0.03 | OFF |
| `SOXS.2025-05-14T04:13:47.944` | VIS | 6/8 | OBJECT | 0.0, +0.1 | +30.0, +30.1 | +20.5, +26.9 | OFF | 0.03 | OFF |
| `SOXS.2025-05-14T04:14:21.554` | NIR | 5/8 | OBJECT | 0.0, +0.1 | +30.0, +30.1 | +20.5, +26.9 | OFF | 0.04 | OFF |
| `SOXS.2025-05-14T04:19:24.078` | VIS | 8/8 | OBJECT | -30.0, -30.1 | 0.0, 0.0 | -9.6, -3.2 | ON | 34.02 | ON |
| `SOXS.2025-05-14T04:19:57.675` | NIR | 7/8 | OBJECT | -30.0, -30.1 | 0.0, 0.0 | -9.6, -3.2 | ON | 11.18 | ON |
| `SOXS.2025-11-30T05:27:27.385` | NIR | 1/4 | OBJECT | 0.0, -0.4 | 0.0, -0.4 | +5.8, -7.6 | ON | 6.53 | ON |
| `SOXS.2025-11-30T05:27:29.300` | VIS | 2/4 | OBJECT | 0.0, -0.4 | 0.0, -0.4 | +5.8, -7.6 | ON | 4.28 | ON |
| `SOXS.2025-11-30T05:32:02.820` | NIR | 3/4 | SKY | -2.0, +20.4 | -2.0, +20.0 | +3.7, +12.8 | OFF | 0.03 | OFF |
| `SOXS.2025-11-30T05:32:04.714` | VIS | 4/4 | SKY | -2.0, +20.4 | -2.0, +20.0 | +3.7, +12.8 | OFF | 0.01 | OFF |

Findings:

- The change in pointing between consecutive positions equals the `FIXOFF` step to within 0.7 arcsec (for example CD-44 12736: pointing moves by (+30.1, +30.7), then (0.0, +0.1), then (-30.1, -30.1)). So `FIXOFF` is the step, and the running sum is the position.
- CD-44 12736 is the documented object-sky-sky-object (OSSO) cycle. Its ON trace is 100 to 1000 times its OFF trace in both arms.
- Every one of the 16 frames is classified correctly by the current rule. But 8 of them (the six first-position frames, steps (0, -0.4) and (0, -0.8), and the two repeated-sky frames, step (0, +0.1)) have a step that is only jitter, and the rule's answer depends on the sign of that jitter.
- The first step is (0, -0.4 × `JITTER WIDTH`) in all three runs. That looks systematic, but three runs are too few to rely on, and the repeated-sky step (+0.1) does not follow the same pattern.
- `DPR TYPE` marks the sky position `SKY` in the November 2025 run but not in the May 2025 runs, although the template version is the same. The manual lists `DPR.TYPE = OBJECT` as fixed.

## GenericOffset frame by frame

There are no offset cards, so the current rule raises `KeyError` on every frame. Pointing is the offset from `TEL TARG`; Moved is the distance from the pointing of the template's first exposure.

| OB, template start | Frame | Arm | EXPNO | DPR TYPE | Pointing (arcsec) | Moved (arcsec) | Trace | Truth |
|---|---|---|---|---|---|---|---|---|
| SN2025kzr, 2025-06-06T00:17:35 | `SOXS.2025-06-06T00:17:55.956` | VIS | 1/2 | OBJECT | -4.5, -3.6 | 0.0 | 0.102 | ON (faint) |
| | `SOXS.2025-06-06T00:48:31.368` | VIS | 2/2 | SKY | +25.5, +26.4 | 42.4 | 0.002 | OFF |
| GD 71, 2025-10-25T07:20:19 | `SOXS.2025-10-25T07:20:41.198` | VIS | 1/4 | OBJECT | +2.8, -16.3 | 0.0 | 5.46 | ON |
| | `SOXS.2025-10-25T07:26:14.359` | VIS | 2/4 | SKY | +4.8, -14.3 | 2.8 | 4.96 | ON (star still on slit) |
| | `SOXS.2025-10-25T07:31:50.662` | VIS | 3/4 | OBJECT | +8.8, -10.3 | 8.5 | 0.007 | OFF (no star) |
| | `SOXS.2025-10-25T07:37:26.056` | VIS | 4/4 | SKY | +14.8, -4.2 | 17.0 | 0.007 | OFF |
| GD 71, 2025-10-25T07:43:00 | `SOXS.2025-10-25T07:43:23.158` | VIS | 1/4 | OBJECT | +2.8, -16.3 | 0.0 | 4.55 | ON |
| | `SOXS.2025-10-25T07:48:56.690` | VIS | 2/4 | SKY | +4.8, -14.3 | 2.8 | 0.20 | OFF, with a weak star wing (4 % of ON) |
| | `SOXS.2025-10-25T07:54:36.013` | VIS | 3/4 | OBJECT | +8.8, -10.3 | 8.5 | 0.003 | OFF (no star) |
| | `SOXS.2025-10-25T08:00:12.200` | VIS | 4/4 | SKY | +14.8, -4.2 | 17.0 | 0.003 | OFF |
| GD 71, 2025-10-25T08:05:46 | `SOXS.2025-10-25T08:06:03.230` | VIS | 1/4 | OBJECT | +2.8, -16.3 | 0.0 | 3.51 | ON |
| | `SOXS.2025-10-25T08:11:36.288` | VIS | 2/4 | SKY | +4.8, -14.3 | 2.8 | 0.16 | OFF, with a weak star wing (5 %) |
| | `SOXS.2025-10-25T08:17:14.525` | VIS | 3/4 | OBJECT | +8.8, -10.3 | 8.5 | 0.003 | OFF (no star) |
| | `SOXS.2025-10-25T08:22:50.705` | VIS | 4/4 | SKY | +14.8, -4.2 | 17.0 | 0.005 | OFF |
| CD-32 5613, 2025-10-26T07:43:35 | `SOXS.2025-10-26T07:43:52.968` | VIS | 1/4 | OBJECT | -31.1, +28.1 | 0.0 | 20.8 | ON |
| | `SOXS.2025-10-26T07:49:26.187` | VIS | 2/4 | OBJECT | -29.1, +30.1 | 2.8 | 21.0 | ON (star still on slit) |
| | `SOXS.2025-10-26T07:55:02.270` | VIS | 3/4 | OBJECT | -25.1, +34.1 | 8.5 | 0.03 | OFF (no star) |
| | `SOXS.2025-10-26T08:00:38.567` | VIS | 4/4 | OBJECT | -19.2, +40.1 | 17.0 | 0.08 | OFF (no star) |
| CD-32 5613, 2025-10-26T08:06:12 | `SOXS.2025-10-26T08:06:33.191` | VIS | 1/4 | OBJECT | -31.1, +28.1 | 0.0 | 16.2 | ON |
| | `SOXS.2025-10-26T08:12:05.759` | VIS | 2/4 | OBJECT | -29.1, +30.1 | 2.8 | 2.05 | partial (13 % of ON) |
| | `SOXS.2025-10-26T08:17:40.762` | VIS | 3/4 | OBJECT | -25.1, +34.1 | 8.5 | 0.01 | OFF (no star) |
| | `SOXS.2025-10-26T08:23:16.104` | VIS | 4/4 | OBJECT | -19.2, +40.1 | 17.0 | 0.07 | OFF (no star) |
| CD-32 5613, 2025-10-26T08:28:48 | `SOXS.2025-10-26T08:29:03.235` | VIS | 1/4 | OBJECT | -31.1, +28.1 | 0.0 | 14.6 | ON |
| | `SOXS.2025-10-26T08:34:36.771` | VIS | 2/4 | OBJECT | -29.1, +30.1 | 2.8 | 1.41 | partial (10 % of ON) |
| | `SOXS.2025-10-26T08:40:16.082` | VIS | 3/4 | OBJECT | -25.1, +34.1 | 8.5 | 0.01 | OFF (no star) |
| | `SOXS.2025-10-26T08:45:52.227` | VIS | 4/4 | OBJECT | -19.2, +40.1 | 17.0 | 0.04 | OFF (no star) |
| M 77, 2025-12-01T02:43:49 | `SOXS.2025-12-01T02:44:10.355` | NIR | 1/4 | OBJECT | +5.4, -7.9 | 0.0 | 9.23 | ON |
| | `SOXS.2025-12-01T02:44:12.268` | VIS | 2/4 | OBJECT | +5.4, -7.9 | 0.0 | 5.51 | ON |
| | `SOXS.2025-12-01T02:46:47.588` | NIR | 3/4 | SKY | +3.4, +12.1 | 20.1 | 0.08 | OFF |
| | `SOXS.2025-12-01T02:46:49.459` | VIS | 4/4 | SKY | +3.4, +12.1 | 20.1 | 0.01 | OFF |
| M 77, 2025-12-10T02:18:56 | `SOXS.2025-12-10T02:19:12.454` | NIR | 1/4 | OBJECT | +3.2, -7.7 | 0.0 | 22.5 | ON |
| | `SOXS.2025-12-10T02:19:14.369` | VIS | 2/4 | OBJECT | +3.2, -7.7 | 0.0 | 4.65 | ON |
| | `SOXS.2025-12-10T02:21:49.284` | NIR | 3/4 | SKY | +1.2, +12.3 | 20.1 | 0.19 | OFF |
| | `SOXS.2025-12-10T02:21:51.149` | VIS | 4/4 | SKY | +1.2, +12.3 | 20.1 | 0.02 | OFF |
| EG21, 2026-01-24T05:32:58 | 8 frames, `SOXS.2026-01-24T05:33:26.555` to `05:50:30.271` | VIS, NIR | 1-8/8 | OBJECT | -4.1 to -3.6, -4.8 | 0.5 or less | VIS 11.2 to 12.6, NIR 0.29 to 0.37 | all ON |
| CD-32 5613, 2026-01-24T06:49:38 | 8 frames, `SOXS.2026-01-24T06:50:03.625` to `07:07:12.278` | VIS, NIR | 1-8/8 | OBJECT | -24.1, +28.3 to +29.3 | 1.0 or less | VIS 6.6 to 24.1, NIR 0.23 to 2.76 | all ON |

Findings:

- SN2025kzr and the two M 77 GenericOffset runs use the template as intended (object, then sky at 20 to 42 arcsec). `DPR TYPE` is right for all 10 frames. The M 77 sky frames still hold diffuse galaxy light across the slit (higher slit-filling level than the SN sky frame), but no nucleus trace.
- The GD 71 and October CD-32 5613 OBs are test patterns with cumulative steps of (+2, +2), (+4, +4), (+6, +6) arcsec. The slit position angle was 0 deg (north-south), so the first 2 arcsec RA step should take the star off the 1 arcsec slit. In the first run of each OB, the second exposure still has the full star on the slit, although its pointing card says the telescope moved. In the later runs the same exposure has 4 to 13 % of the star. No header card explains the difference. The GD 71 `SKY` label at position 2 is wrong in run 1, and the `OBJECT` label at position 3 is wrong in all three runs. The October CD-32 5613 OB has no `SKY` label at all, and half its `OBJECT` frames hold no star.
- The January 2026 EG21 and CD-32 5613 OBs (named `..._stare_15`) used GenericOffset without moving the telescope: all 16 frames are on the star and labelled `OBJECT`. They are stare data in an offset template.

## Why the current rule is fragile

1. **It reads a step as a position.** `FIXOFF` is the move from the previous exposure. At a repeated position (the second sky exposure of an OSSO cycle, or the first exposure, which only jitters around the acquisition position) the step is a jitter of a few tenths of an arcsec with either sign. In this archive the signs happened to fall the right way for all 8 such frames.
2. **It assumes the sky offset is positive.** The template accepts `SEQ.FIXOFF.RA` and `SEQ.FIXOFF.DEC` from -100 to +100 arcsec. With a sky offset whose components sum to a negative number, such as (-20, +10), every step out to the sky is called ON and every step back to the object is called OFF. M 77 (-2.0, +20.4) is still positive.
3. **It needs cards that only one template writes.** GenericOffset frames have no `FIXOFF`, so the recipe would raise `KeyError` before it could report anything useful.
4. **It counts on equal ON and OFF numbers.** That holds for OSSO, but not for a GenericOffset list such as `O O S` or for a set with no sky.

## How the data organiser filters offset frames

- `soxs_sof_map.yaml` routes only `DPR TECH = ECHELLE,SLIT,OFFSET` to `offset_obj` (`science_offset`, `DPR TYPE` OBJECT or OBJECT,ASYNC) and `offset_std` (`flux_standard_offset`, `telluric_standard_offset`). There is no entry for `ECHELLE,SLIT,GENOFFSET`, so the 50 GenericOffset frames are never put in a SOF [S1].
- `data_organiser._flag_files_to_ignore` sets `ignore = 1` on every frame whose `DPR TYPE` is not in the list of reducible types. `SKY` is not in that list, so all 13 `SKY` offset frames (2 FixedSkyOffset, 11 GenericOffset) carry `ignore = 1` in the archive database [S1][S2].
- Result in the archive: the V818 Sco and CD-44 12736 offset SOFs are complete (their sky frames are labelled `OBJECT`). The M 77 FixedSkyOffset SOFs (`20251130T052727_NIR_..._M77.sof` and `20251130T052729_VIS_..._M77.sof`) hold one ON frame and no OFF frame each, so `_check_on_off_balance` refuses them [S2].

## Recommended rule and what it cannot classify

For `ECHELLE,SLIT,OFFSET` (FixedSkyOffset, SpecphotStdOffset):

1. Group the frames by `TPL START` (one template run). Both arms share one pointing sequence.
2. Take as the object position the pointing (`RA`, `DEC`) of the frame with the lowest `TPL EXPNO` in the run. The manual's OSSO order guarantees that frame is on the object.
3. For each frame, compute the angular distance from that position. ON if the distance is at most `SEQ JITTER WIDTH` + 1 arcsec; OFF otherwise. (With the observed sky offsets of 20 to 43 arcsec and jitter widths of 1 to 2 arcsec, the margin is large.)
4. Cross-check with the running sum of `SEQ FIXOFF RA/DEC` per arm in `TPL EXPNO` order, and stop with a clear error if a frame's two positions differ by more than about 1 arcsec, or if any ON or OFF group is empty.
5. Do not use the sign of `FIXOFF`, and do not use `DPR TYPE` (it is `OBJECT` at the sky position in the May 2025 runs).

This works for all 16 frames. It needs the first exposure of the run in the SOF, which the organiser's template grouping already provides. The pointing cards are used only for differences within one run, so the acquisition residual and proper motion that make them a poor target position (DY-2296) cancel out.

For `ECHELLE,SLIT,GENOFFSET` (GenericOffset), if it is to be reduced:

1. Route it in `soxs_sof_map.yaml`, and stop ignoring its `SKY` frames.
2. OFF = `DPR TYPE` contains `SKY`. ON = `DPR TYPE` is `OBJECT` and the pointing is within 1 arcsec of the template run's first `OBJECT` frame.
3. Refuse, by name, any `OBJECT` frame farther than that (it is a mapping position, not an ON frame), and refuse a run with no `SKY` frame, pointing the user to the stare recipe.

This is right for SN2025kzr and both M 77 runs. It refuses the January 2026 EG21 and CD-32 5613 runs (no sky), which is correct for an offset reduction.

Frames no header-based rule can classify correctly (24, all test OBs):

- GD 71, OB 1000603202 (12 VIS frames). Position-3 `OBJECT` frames have no star (refused by the rule above, which is the right outcome). The run-1 position-2 `SKY` frame `SOXS.2025-10-25T07:26:14.359` holds the full star and would be used as sky; the run-2 and run-3 position-2 `SKY` frames hold a 4 to 5 % star wing. Nothing in the header separates these cases.
- CD-32 5613, OB 1000603272 (12 VIS frames). No `SKY` labels; position-2 frames hold 10 to 100 % of the star with identical headers. The rule refuses the whole OB (no sky), which is right for an offset reduction, but it cannot say which frames are truly on the target.

For Phase 3, the frame accounting (`EXPTIME`, `NCOMBINE`, `MJD-OBS`, `MJD-END`, `PROVi`; DY-2295) should come from the ON list produced by the recommended rule. Only the OB-level target position `TEL TARG` is the same for ON and OFF frames, so `RA`/`DEC` (DY-2296) do not depend on this choice.

## Sources

- [S1] soxspipe source at `develop` 6bdef45d: `soxspipe/recipes/soxs_offset.py` (`_split_on_off_frames`, `_check_on_off_balance`), `soxspipe/resources/soxs_keywords.yaml` (`OFFSET_RA`, `OFFSET_DEC`), `soxspipe/resources/soxs_sof_map.yaml` (`flux_standard_offset`, `telluric_standard_offset`, `science_offset`), `soxspipe/commonutils/data_organiser.py` (`_flag_files_to_ignore`).
- [S2] SOXS raw frames and database, `/Volumes/SOXS/raw_archive/soxs_commissioning_data`: `raw_frames` table of `soxspipe.db` (opened read-only), the six `sof/*OFFSET_OBJ*.sof` files, and the primary headers and pixel data of the 66 frames listed above (`raw/<night>/SOXS.*.fits.Z`). Survey run on 2026-11-19.
- [S3] SOXS User Manual, SOXS-MAN-0003, issue 2.1, 2026-07-07, sections 4.2.6, 4.2.7 and 6.1, template tables on pages 42, 43 and 47. <https://soxs.inaf.it/download/soxs-man-0003_user-manual_2.1.pdf>.
- [S4] DY-2296 research notes, branch `research/DY-2296-target-ra-dec`, `docs/research/DY-2296-target-ra-dec.md` (decoding of `TEL TARG`, pointing cards, first observation that `FIXOFF` is a step).
