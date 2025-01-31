:::{table} Input files for the `soxs_mbias` recipe. The files are typically passed to the `soxs_mbias` recipe via a set-of-file (sof) file listing one file per line.
:name: table_soxs_mbias_input

| Data Type | Content | Related OB | Min. Frame Count |
|:----|:----|:---|:---|
| FITS images | raw bias frames (UV-VIS/AC exposures with exptime = 0) | `SOXS_img_cal_Bias`, `SOXS_gen_cal_VISBias` | 5 |


:::
