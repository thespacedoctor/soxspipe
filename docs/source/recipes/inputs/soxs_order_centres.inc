:::{table} Input files for the `soxs_order_centres` recipe. The files are typically passed to the `soxs_order_centres` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_order_centres_input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Flat lamp through a single-pinhole mask | `SOXS_slt_cal_VISLampFlatPinhole`, `SOXS_slt_cal_NIRLampFlatPinhole` |
| FITS Image | Master Dark Frame (VIS only, optional) | - |
| FITS Image | Master Bias Frame (VIS only) | - |
| FITS Image | Dark frame (Lamp-Off) of equal exposure length as single-pinhole frame (Lamp-On) (NIR only) | `SOXS_slt_cal_NIRLampFlatPinhole` |
| FITS Table | First guess dispersion solution | - |


:::


