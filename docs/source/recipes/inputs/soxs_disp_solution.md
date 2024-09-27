:::{table} Input files for the `soxs_disp_solution` recipe. The files are typically passed to the `soxs_disp_solution` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_disp_solution_input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Arc Lamp through single pinhole mask | `SOXS_slt_cal_VISArcsPinhole`, `SOXS_slt_cal_NIRArcsPinhole`|
| FITS Image | Master Dark Frame (VIS only, optional) | - |
| FITS Image | Master Bias Frame (VIS only) | - |
| FITS Image | Dark frame (Lamp-Off) of equal exposure length as single pinhole frame (Lamp-On) (NIR only) | `SOXS_slt_cal_NIRArcsPinhole` |

:::



