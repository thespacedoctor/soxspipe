:::{table} Output files for the `soxs_nod` recipe and their respective ESO PRO keywords. Some output file may not be generated depending on the user's pipeline settings.
:name: soxs_nod_output
:widths: 15, 35, 10, 20, 10, 10

| Label                     | Content                                      | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------------- | -------------------------------------------- | --------- |---------- |---------- |---------- |
| `EXTRACTED ORDERS TABLE`    | Table of the extracted source in each order  | FITS      |     `SCI_SLIT_FLUX_NIR`      |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `EXTRACTED MERGED TABLE`    | Table of the extracted, order-merged         | FITS      |     `SCI_SLIT_FLUX_NIR`      |      `REDUCED`     |     `ECHELLE,SLIT`      |
| `EXTRACTED MERGED ASCII`    | Ascii version of extracted source spectrum   | TXT       |     -       |    -       |     -      |
| `OBJECT TRACE RES`          | Residuals of the object trace polynomial fit | PDF       |     -       |    -       |     -      |
| `BKGROUND`                  | Fitted intra-order image background          | PDF       |     -       |    -       |     -      |
| `EXTRACTED ORDERS QC PLOT`  | QC plot of extracted source                  | PDF       |     -       |    -       |     -      |
| `EXTRACTED MERGED QC PLOT`  | QC plot of extracted order-merged source     | PDF       |     -       |    -       |     -      |

:::


