:::{table} Input files for the `soxs_spatial_solution` recipe. The files are typically passed to the `soxs_spatial_solution` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_spatial_solution_input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Arc Lamp through multi-pinhole mask | `SOXS_slt_cal_VISArcsMultiplePinhole`, `SOXS_slt_cal_NIRArcsMultiplePinhole` |
| FITS Image | Master Dark Frame (VIS only, optional) | - |
| FITS Image | Master Bias Frame (VIS only) | - |
| FITS Image | Dark frame (Lamp-Off) of equal exposure length as multi-pinhole frame (Lamp-On) (NIR only) | `SOXS_slt_cal_NIRArcsMultiplePinhole` |
| FITS Table | First-guess [Dispersion Map](../files/dispersion_map.md) table |

:::


