:::{table} Input files for the `soxs_stare` recipe. The files are typically passed to the `soxs_stare` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_stare_input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Images| Raw science frames of targets observed in stare mode | `SOXS_slt_obs_StareSynchro`, `SOXS_slt_cal_SpecphotStdStare`, `SOXS_slt_cal_TelluricStdStare` |
| FITS Image | Master bias (UVVIS only) | - |
| FITS Image | Master dark frame (NIR only) | - |
| FITS Image | Master flat frame (optional) | - |
| FITS Table | Order location table containing coefficients to the polynomial fits describing the order locations. | - |
| FITS Table | Dispersion map table giving coefficients of polynomials describing 2D dispersion/spatial solution | - |
| FITS Image | Dispersion Map image with 3-extensions (wavelength, slit-position and echelle order number) | - |


:::


