:::{table} Output files for the `soxs_offset` recipe and their respective ESO PRO keywords. Some output file may not be generated depending on the user's pipeline settings.
:name: soxs_offset_output
:widths: 15, 35, 10, 20, 10, 10


| Label                     | Content                                      | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------------- | -------------------------------------------- | --------- |---------- |---------- |---------- |
| `EXTRACTED ORDERS TABLE`    | Table of the extracted source in each order  | FITS      |     `SCI_SLIT_FLUX_NIR`      |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `SPECTRUM TABLE` | Table of the extracted, order-merged spectrum, `<SOF stem>_SPECTRUM.fits`: `WAVE`, `FLUX_COUNTS`, `ERR_COUNTS`, `SNR`. When a response product exists, the flux calibration replaces it with the table of `WAVE`, `FLUX`, `ERR`, `FLUX_COUNTS`, `ERR_COUNTS`, `SNR`. `ERR_COUNTS` is the error of `FLUX_COUNTS` in electrons, and `ERR` is `ERR_COUNTS` scaled by the same factors as `FLUX`. A standard star also writes `<SOF stem>_SPECTRUM_NOTFLAT.fits`, with a `FLUX_DENSITY_COUNTS` column, for the response curve. | FITS | `SCI_SLIT_FLUX_<ARM>` | `REDUCED` | `ECHELLE,SLIT` |
| `SPECTRUM ASCII`    | ASCII spectrum of the `<SOF stem>_SPECTRUM.fits` table, `WAVE` in Å: `<SOF stem>_SPECTRUM.txt`  | TXT       |     -       |    -       |     -      |
| `LACOSMIC_MASK` | QC product: the pixels L.A.Cosmic flagged as cosmic rays, without the default bad-pixel mask. One `CRMASK` image extension (1 = flagged, 0 = not) per raw frame in observation time order. Each header gives the raw file name (`ORIGFILE`), `MJD-OBS`, `DATE-OBS` and the flagged-pixel count (`NCRPIX`). Written to the QC directory as `_LACOSMIC_MASK.fits`, only if `use_lacosmic` is on | FITS      |     -       |    -       |     -      |
| `OBJECT TRACE RES`          | Residuals of the object trace polynomial fit | PDF       |     -       |    -       |     -      |
| `BKGROUND`                  | Fitted intra-order image background          | PDF       |     -       |    -       |     -      |
| `EXTRACTED ORDERS QC PLOT`  | QC plot of extracted source                  | PDF       |     -       |    -       |     -      |
| `SLIT DRIFT QC PLOT`      | QC plot of the object slit position against wavelength before rectification. The legend gives the slit drift, the science minus multi-pinhole observation time, and the science minus multi-pinhole `ABSROT` (the multi-pinhole angle is assumed to be 0 deg when its frame has no `ABSROT`). A vertical line marks where the trace crosses slit position 0 | PDF       |     -       |    -       |     -      |
| `SLIT ALIGNMENT QC PLOT` | QC plot of the slit shift applied to each frame (or to each order of a frame) before stacking. Written only when more than one frame is aligned | PDF | - | - | - |
| `EXTRACTED MERGED QC PLOT`  | QC plot of extracted order-merged source     | PDF       |     -       |    -       |     -      |
| `EXTRACTED MERGED FLUXCALIBRATED QC PLOT` | QC plot of extracted order-merged, flux-calibrated source | PDF | - | - | - |

:::

