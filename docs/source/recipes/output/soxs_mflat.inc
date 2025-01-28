:::{table} Output files for the `soxs_mflat` recipe and their respective ESO PRO keywords (NIR and VIS).
:name: soxs_mflat_output
:widths: 15, 35, 10, 20, 10, 10

| Label         | Content                                                       | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------- | ------------------------------------------------------------- | --------- |---------- |---------- |---------- |
| `ORDER_LOC`     | table of coefficients from polynomial fits to order locations | FITS      |     `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT`         | master spectroscopic flat frame                           | FITS      |      `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER_LOC_RES` | visualisation of goodness of order edge fitting               | PDF       |      -     |     -      |     -      |
| `BKGROUND`      | Fitted intra-order image background                           | PDF       |       -    |     -      |     -      |


:::




:::{table} Output files for the `soxs_mflat` recipe for Xshooter UVB arm.
:name: soxs_mflat_output_uvvis
:widths: 15, 35, 10, 20, 10, 10

| Label               | Content                                                       | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------- | ------------------------------------------------------------- | --------- |---------- |---------- |---------- |
| `ORDER_LOC_DLAMP`     | table of coefficients from polynomial fits to order locations | FITS      |  `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT_DLAMP`         | UVB master spectroscopic flat frame (DLAMP)                   | FITS      |       `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER_LOC_QLAMP`     | table of coefficients from polynomial fits to order locations | FITS      |    `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT_QLAMP`         | UVB master spectroscopic flat frame (QLAMP)                   | FITS      |             `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER_LOC`           | table of coefficients from polynomial fits to order locations | FITS      |    `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT`               | UVB master spectroscopic flat frame                           | FITS      |           `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER_LOC_RES_DLAMP` | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |
| `BKGROUND_DLAMP`      | Fitted intra-order image background DLAMP                     | PDF       |    -       |    -       |     -      |
| `ORDER_LOC_RES_QLAMP` | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |
| `BKGROUND_QLAMP`      | Fitted intra-order image background QLAMP                     | PDF       |    -       |    -       |     -      |
| `ORDER_LOC_RES`       | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |


:::
