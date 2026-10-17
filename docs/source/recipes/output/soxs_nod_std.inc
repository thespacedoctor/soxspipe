:::{table} Output files for the `soxs_nod_std` recipe and their respective ESO PRO keywords. Some output file may not be generated depending on the user's pipeline settings.
:name: soxs_nod_std_output
:widths: 15, 35, 10, 20, 10, 10


| Label                     | Content                                      | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------------- | -------------------------------------------- | --------- |---------- |---------- |---------- |
| `EXTRACTED ORDERS TABLE`    | Table of the extracted source in each order  | FITS      |     `SCI_SLIT_FLUX_NIR`      |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `EXTRACTED MERGED TABLE`    | Table of the extracted, order-merged         | FITS      |     `SCI_SLIT_FLUX_NIR`      |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `EXTRACTED FLUXCAL SPECTRUM` | Table of the flux calibrated extracted spectrum | FITS | `SCI_SLIT_FLUX_<ARM>` | `REDUCED` | `ECHELLE,SLIT` |
| `RESPONSE_FUNC` | Response function coefficients | FITS | `RESP_TAB_<ARM>` | `REDUCED` | `ECHELLE,SLIT` |
| `EFFICIENCY` | Efficiency estimate | FITS | `EFFICIENCY_TAB_<ARM>` | `REDUCED` | `ECHELLE,SLIT` |
| `EXTRACTED MERGED ASCII`    | Ascii version of extracted source spectrum   | TXT       |     -       |    -       |     -      |
| `OBJECT TRACE RES`          | Residuals of the object trace polynomial fit | PDF       |     -       |    -       |     -      |
| `BKGROUND`                  | Fitted intra-order image background          | PDF       |     -       |    -       |     -      |
| `EXTRACTED ORDERS QC PLOT`  | QC plot of extracted source                  | PDF       |     -       |    -       |     -      |
| `SLIT DRIFT QC PLOT`      | QC plot of the object slit position against wavelength before rectification. The legend gives the slit drift, the science minus multi-pinhole observation time, and the science minus multi-pinhole `ABSROT` (the multi-pinhole angle is assumed to be 0 deg when its frame has no `ABSROT`). A vertical line marks where the trace crosses slit position 0 | PDF       |     -       |    -       |     -      |
| `SLIT ALIGNMENT QC PLOT` | QC plot of the slit shift applied to each frame (or to each order of a frame) before stacking. Written only when more than one frame is aligned | PDF | - | - | - |
| `EXTRACTED MERGED QC PLOT`  | QC plot of extracted order-merged source     | PDF       |     -       |    -       |     -      |
| `EXTRACTED MERGED FLUXCALIBRATED QC PLOT` | QC plot of extracted order-merged, flux-calibrated source | PDF | - | - | - |
| `RESPONSE QC PLOT` | Response curve QC plot. | PDF | - | - | - |

:::

