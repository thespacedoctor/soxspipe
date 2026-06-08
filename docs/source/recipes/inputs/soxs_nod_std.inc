:::{table} Input files for the `soxs_nod_std` recipe. The files are typically passed to the `soxs_nod_std` recipe via a set-of-file (sof) file listing one file per line.
:name: soxs_nod_std_input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Raw spectroscopic standard star frames observed in nodding mode |`SOXS_slt_obs_AutoNodOnSlit`, `SOXS_slt_cal_TelluricStdNod`, `SOXS_slt_cal_SpecphotNod`|
| FITS Image | Master flat frame (optional) | - |
| FITS Table | order location table containing coefficients to the polynomial fits describing the order locations. | - |
| FITS Table | Dispersion map table giving coefficients of polynomials describing 2D dispersion/spatial solution | - |
| FITS Image | Dispersion map FITS image with 3-extensions (wavelength, slit-position and echelle order number) | - |


:::


