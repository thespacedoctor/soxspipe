:::{table} Output files for the `soxs_stare` recipe. Some output file may not be generated depending on the user's pipeline settings.
:name: soxs_stare_output

| Label                     | Content                                      | Data Type |
| ------------------------- | -------------------------------------------- | --------- |
| `SKY_SUBTRACTED_OBJECT`     | The sky-subtracted object                    | FITS      |
| `SKY_MODEL`                 | The sky background model                     | FITS      |
| `SKY_SUB_RESIDUALS`         | The sky subtraction residuals                | FITS      |
| `EXTRACTED_ORDERS_TABLE`    | Table of the extracted source in each order  | FITS      |
| `EXTRACTED_MERGED_ASCII`    | Ascii version of extracted source spectrum   | TXT       |
| `EXTRACTED_MERGED_TABLE`    | Table of the extracted, order-merged         | FITS      |
| `BKGROUND`                  | Fitted intra-order image background          | PDF       |
| `SKY_MODEL_QC_PLOTS`        | QC plots for the sky-background modelling    | PDF       |
| `SKY SUBTRACTION QUICKLOOK` | Sky-subtraction quicklook                    | PDF       |
| `OBJECT_TRACE_RES`          | Residuals of the object trace polynomial fit | PDF       |
| `EXTRACTED_ORDERS_QC_PLOT`  | QC plot of extracted source                  | PDF       |
| `EXTRACTED_MERGED_QC_PLOT`  | QC plot of extracted order-merged source     | PDF       |

:::

