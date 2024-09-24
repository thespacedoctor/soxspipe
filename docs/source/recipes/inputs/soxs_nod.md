:::{table} Input files for the `soxs_nod` recipe. The files are typically passed to the `soxs_nod` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_nod_input

| Data Type | Content | Related OB |
|:----|:----|:---|
|FITS images|Raw science frames of targets observed in nodding mode |`SOXS_nod`|
| FITS Image | Master Bias (UVVIS only) | - |
| FITS Table | [order table](../files/order_table.md) containing coefficients to the polynomial fits describing the order centre locations. | |
| FITS Table | [Dispersion Map](../files/dispersion_map.md) table giving coefficients of polynomials describing 2D dispersion/spatial solution | |

:::


