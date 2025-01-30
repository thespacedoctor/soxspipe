:::{table} Output files for the `soxs_spatial_solution` recipe and their respective ESO PRO keywords.
:name: soxs_spatial_solution_output
:widths: 15, 35, 10, 20, 10, 10

| Label                  | Content                                                | Data Type | PRO CATG | PRO TYPE | PRO TECH |
| ---------------------- | ------------------------------------------------------ | --------- |---------- |---------- |---------- |
| `SPAT_SOL`               | Full dispersion-spatial solution                       | FITS      |    `DISP_TAB_<ARM>`       |    `REDUCED`       |     `ECHELLE,MULTI-PINHOLE`      |
| `2D_MAP`                 | 2D detector map of wavelength, slit position and order | FITS      |     `DISP_IMAGE_<ARM>`      |    `REDUCED`       |     `ECHELLE,MULTI-PINHOLE`      |
| `BKGROUND`               | Fitted intra-order image background                    | PDF       |    -       |    -       |     -      |
| `DISP_MAP_LINES`         | Dispersion solution fitted lines                       | FITS      |   -       |    -       |     -      |
| `DISP_MAP_LINES_MISSING` | Undetected arc lines                                   | FITS      |    -       |    -       |     -      |
| `DISP_MAP_RES`           | Dispersion solution QC plots                           | PDF       |      -       |    -       |     -      |

:::




