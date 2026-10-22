:::{table} Output files for the `soxs_stare` recipe and their respective ESO PRO keywords. Some output file may not be generated depending on the user's pipeline settings.
:name: soxs_stare_output
:class: longtable
:widths: 15, 35, 10, 20, 10, 10

| Label                     | Content                                      | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------------- | -------------------------------------------- | --------- |---------- |---------- |---------- |
| `EXTRACTED ORDERS TABLE`    | Table of the extracted source in each order  | FITS      |     `SCI_SLIT_FLUX_<ARM>`     |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `EXTRACTED MERGED TABLE`    | Table of the extracted, order-merged         | FITS      |      `SCI_SLIT_FLUX_<ARM>`      |     `REDUCED`      |      `ECHELLE,SLIT`     |
| `EXTRACTED FLUXCAL SPECTRUM` | Table of the flux calibrated extracted spectrum | FITS | `SCI_SLIT_FLUX_<ARM>` | `REDUCED` | `ECHELLE,SLIT` |
| `SKY SUBTRACTED OBJECT`     | The sky-subtracted object                    | FITS      |       -       |    -       |     -      |
| `SKY MODEL`                 | The sky background model                     | FITS      |         -       |    -       |     -      |
| `SKY SUB RESIDUALS`         | The sky subtraction residuals                | FITS      |        -       |    -       |     -      |
| `EXTRACTED MERGED ASCII`    | Ascii version of extracted source spectrum   | TXT       |       -       |    -       |     -      |
| `BKGROUND`                  | Fitted intra-order image background          | PDF     |     -       |    -       |     -      |
| `SKY MODEL QC PLOTS`        | QC plots for the sky-background modelling    | PDF       |     -       |    -       |     -      |
| `SKY SUBTRACTION QUICKLOOK` | Sky-subtraction quicklook                    | PDF       |     -       |    -       |     -      |
| `LACOSMIC_MASK` | QC product: the pixels L.A.Cosmic flagged as cosmic rays, without the default bad-pixel mask. One `CRMASK` image extension (1 = flagged, 0 = not) per raw frame in observation time order. Each header gives the raw file name (`ORIGFILE`), `MJD-OBS`, `DATE-OBS` and the flagged-pixel count (`NCRPIX`). Written to the QC directory as `_LACOSMIC_MASK.fits`, only if `use_lacosmic` is on | FITS      |     -       |    -       |     -      |
| `OBJECT TRACE RES`          | Residuals of the object trace polynomial fit | PDF       |     -       |    -       |     -      |
| `EXTRACTED ORDERS QC PLOT`  | QC plot of extracted source                  | PDF       |     -       |    -       |     -      |
| `SLIT DRIFT QC PLOT`      | QC plot of the object slit position against wavelength before rectification. The legend gives the slit drift, the science minus multi-pinhole observation time, and the science minus multi-pinhole `ABSROT` (the multi-pinhole angle is assumed to be 0 deg when its frame has no `ABSROT`). A vertical line marks where the trace crosses slit position 0 | PDF       |     -       |    -       |     -      |
| `SLIT_ILLUMINATION_QC_PLOT` | QC plot of the skyline slit profile of each order against the master flat's lit slit range, with the fitted slit offset. Written as `_SLIT_ILLUM_QC_PLOT.pdf`, only if `slit-illumination: measure-offset` is on | PDF       |     -       |    -       |     -      |
| `EXTRACTED MERGED QC PLOT`  | QC plot of extracted order-merged source     | PDF       |     -       |    -       |     -      |
| `EXTRACTED MERGED FLUXCALIBRATED QC PLOT` | QC plot of extracted order-merged, flux-calibrated source | PDF | - | - | - |

:::

