:::{table} Output files for the `soxs_mflat` recipe and their respective ESO PRO keywords (NIR and VIS).
:name: soxs_mflat_output
:widths: 15, 35, 10, 20, 10, 10

| Label         | Content                                                       | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------- | ------------------------------------------------------------- | --------- |---------- |---------- |---------- |
| `ORDER LOC`     | table of coefficients from polynomial fits to order locations | FITS      |     `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT`         | master spectroscopic flat frame                           | FITS      |      `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER LOC RES` | visualisation of goodness of order edge fitting               | PDF       |      -     |     -      |     -      |
| `BKGROUND`      | Fitted intra-order image background                           | PDF       |       -    |     -      |     -      |


:::




:::{table} Output files for the `soxs_mflat` recipe for Xshooter UVB arm.
:name: soxs_mflat_output_uvvis
:widths: 15, 35, 10, 20, 10, 10

| Label               | Content                                                       | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ------------------- | ------------------------------------------------------------- | --------- |---------- |---------- |---------- |
| `ORDER LOC DLAMP`     | table of coefficients from polynomial fits to order locations | FITS      |  `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT DLAMP`         | UVB master spectroscopic flat frame (DLAMP)                   | FITS      |       `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER LOC QLAMP`     | table of coefficients from polynomial fits to order locations | FITS      |    `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT QLAMP`         | UVB master spectroscopic flat frame (QLAMP)                   | FITS      |             `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER LOC`           | table of coefficients from polynomial fits to order locations | FITS      |    `ORDER_TAB_<ARM>`      |    `REDUCED`       |    `ECHELLE,SLIT`       |
| `MFLAT`               | UVB master spectroscopic flat frame                           | FITS      |           `MASTER_FLAT_<ARM>`     |    `REDUCED`       |     `IMAGE`     |
| `ORDER LOC RES DLAMP` | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |
| `BKGROUND DLAMP`      | Fitted intra-order image background DLAMP                     | PDF       |    -       |    -       |     -      |
| `ORDER LOC RES QLAMP` | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |
| `BKGROUND QLAMP`      | Fitted intra-order image background QLAMP                     | PDF       |    -       |    -       |     -      |
| `ORDER LOC RES`       | visualisation of goodness of order edge fitting               | PDF       |    -       |    -       |     -      |


:::
