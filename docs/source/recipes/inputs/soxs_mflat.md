:::{table} Input files for the `soxs_mflat` recipe. The files are typically passed to the `soxs_mflat` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_mflat_input


| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS image |Raw flats frames (exposures with identical exposure time, slit-width and detectors readout parameters). UV-VIS requires separate sets D-Lamp and QTH-Lamp flats. | `SOXS_slt_cal_NIRLampFlat`, `SOXS_slt_cal_NIRLampFlatAtt`, `SOXS_slt_cal_VISLampFlat`, `SOXS_slt_cal_VISLampFlatAtt` |
| FITS Image |A master bias Frame (UV-VIS only) | - |
| FITS Image |A master dark frame or lamp-off frames with identical exposure times to the lamp-on flat frames (NIR only) | - |
| FITS Table |An order location table containing coefficients to the polynomial fit describing the order centre locations. UV-VIS requires separate tables for D-Lamp and QTH-Lamp. | - |

:::

