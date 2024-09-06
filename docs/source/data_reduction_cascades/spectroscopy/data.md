## Spectroscopic Data

### Daily Calibration Data

* Bias frames (UV-VIS only)
* Dark frames (NIR only)
* Through-slit flat field frames (UV-VIS/NIR)
* Single pinhole arc-lamp frames
* Single pinhole flat-field frame
* Multiple pinhole arc-lamp frames
* Standard star observations

### Static Calibration Data

* Reference bad pixel map 
* Line reference table to wavelength calibrate spectra
* Standard star flux table
* Atmospheric extinction table
* Sky Lines reference table

### Intermediate Data Products

* Master bias frame
* Master dark frame
* Master flat-frame
* Order location table
* Dispersion solution table
* Instrument response function

### Output Data Products

| Product                   | Description                                                                                                                                                                                                                                                                                                                                                                                                   |
| ------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| 1D Source Spectra         | 1D spectra in FITS binary table format, one for each arm.  Each FITS spectrum file will contain 4 extensions:  1. Wavelength- and flux-calibrated spectra with absolute flux correction via scaling to acquisition image source photometry, 2. an additional spectrum with correction for telluric absorption via MOLECFIT, 3. the variance array and 4. the sky-background spectra.      |
| 1D Merged Source Spectrum | 1D UV-VIS & NIR merged spectrum in FITS binary table format with PDF visualisation.  This spectrum will be rebinned to a common pixel scale foreach arm.  This spectrum file will also have the same 4 extensions described above.                                                                                                                                 |
| 2D Source Spectra         | A 2D FITS image for each spectral arm containing wavelength and flux calibrated spectra (no other corrections applied) allowing users to perform source extraction with their tool of choice. Note that rectification of the curved orders in the NIR introduces a source of correlated noise not present in extractions performed on the un-straightened orders as done by the pipeline. |
                                                                                                                                                     |
