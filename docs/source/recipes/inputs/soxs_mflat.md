:::{table} Input files for the `soxs_mflat` recipe. The files are typically passed to the `soxs_mflat` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_mflat_input


| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS images | raw flats frames (exposures with identical exposure time and detectors readout parameters). UV-VIS requires separate sets D-Lamp and QTH-Lamp flats. | `SOXS_slt_cal_NIRLampFlat`, `SOXS_slt_cal_NIRLampFlatAtt`, `SOXS_slt_cal_VISLampFlat`, `SOXS_slt_cal_VISLampFlatAtt` |
| FITS Image | Master Bias Frame (UV-VIS only) | - |
| FITS Table | [order table](../files/order_table.md) containing coefficients to the polynomial fits describing the order centre locations. UV-VIS requires separate tables for D-Lamp and QTH-Lamp. | |

:::

