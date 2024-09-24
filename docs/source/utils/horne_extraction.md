# horne_extraction

DAVE WAS HERE

The [`horne_extraction`](#soxspipe.commonutils.horne_extraction) utility performs an optimal extraction of a source's flux from an echelle order using the algorithms reported in {cite:t}`horne1986`.

## Input

| Frame.                   | Description                                   | 
| ------------------------ | --------------------------------------------- |
|  Detrended  2D spectrum | Frame containing the 2D spectrum (any observing mode)  |  
| Dispersion map | FITS table containing the computed dispersion map for the ARM|
| Order table | FITS table containing the spectral format computed for the ARM|

## Parameters

| Parameter                | Description                                   | Type  |
| ------------------------ | --------------------------------------------- | ----- |
|  horne-extraction-slit-length | size of the extraction box in the spatial direction  | int  |
| horne-extraction-profile-clipping-sigma | Sigma-clipping threshold for refining the fitting of polynomials used for determining the object profile  | int|
| horne-extraction-profile-clipping-iteration-count| Maximum number of iterations allowed to refine fitting of polynomials of the object profile |int |
| horne-extraction-profile-global-clipping-sigma| Sigma-clipping threshold used to exclude deviant pixels and cosmic rays during the extraction| int |

## Methods

The typical execution workflow is the following:

![](horne_extraction.png)

This utility, according to the prescription of Horne+86, follows those basic steps:

### 1 Determination of the spatial profile

Starting from the order_table computed by the `detect_continumm` utility, the horne_extraction utility runs along the spectral direction and takes, for each wavelength centred in the position measured by the detect_continuum utility, a window of 1x`horne-extraction-slit-length`. Then, the pixels are summed, and each pixel in the slice is normalized as 

$$
pixel_i = \frac{Flux_{pixel_i}}{\sum_{j=0}^{N}{Flux_{pixel_j}}}
$$

where N = `horne-extraction-slit-length`.


Then, a low-order polynomial is fitted for each pixel in the slice, modelling the value of the fractional flux received along the dispersion axis. The complete set of polynomials represents the object's profile along the spatial direction for each wavelength. For a slice length of N pixels, there will be N different polynomials, one per pixel.

Polynomials are fitted with subsequent iterations. Pixels for which residuals deviate more than the `horne-extraction-profile-clipping-sigma` standard deviation are removed before attempting a new fitting of the data. The procedure is repeated for a maximum of `horne-extraction-profile-clipping-iteration-count` times.

When the procedure above is completed, the actual extraction occurs as follows:

### 2 Extraction of the spectrum for each wavelength for each order.

Using the polynomials computed above, the extracted integrated flux for each wavelength in the order is calculated as follows:

$$
extractedFlux_{\lambda} = \frac{\sum_{i = 0}^{N}{P_{i}^{\lambda} (D_{i}^{\lambda} - S_{i}^{\lambda})}}{\sum_{i = 0}^{N}{(P_{i}^{\lambda})^2 (V_{i}^{\lambda})^{-1}}}
$$

Where $P_{i}^{\lambda}$ is the value of the polynomial for the pixel $i^{th}$ evaluated at wavelength $\lambda$, $D_{i}^{\lambda}$ is the raw pixel value on the 2D image, detrended and sky subtracted, of the science object and $S_{i}^{\lambda}$ is the estimated value, at the same pixel position, of the 2D image of the sky model computed in stare mode. When this utility is applied to nodding or offset mode images, this value is chosen to be zero. $V_{i}^{\lambda}$ is the so-called variance image and is computed as 

$$
V_{i}^{\lambda} = V_{0} + \frac{P_{i}^{\lambda} f_{\lambda} +  S_{i}^{\lambda}}{Q}
$$

Where $V_{0}$ is the square of the readout noise of the detector, $f_\{lambda}$ is the sum of the pixels within the extraction window at each wavelength, and $Q$ is the detector gain. 

Cosmic rays and uncorrected deviant pixels are removed from the extraction by comparing their value with a certain threshold. In particular, for each pixel $i$ in the extraction window the quantity

$$
\frac{(D_{i}^{\lambda} - P_{i}{\lambda} f_{\lambda} - S_{i}{\lambda})^2}{V_{i}^{\lambda}}
$$

is computed. The extraction does not include pixels for which this quantity exceeds the `horne-extraction-profile-global-clipping-sigma`.

### 3 Merging of spectral orders

The spectrum is extracted on each order separately. When all orders are extracted, they are merged together in a single spectrum. Since regions of different orders overlap in the wavelength space, they are first rectified on a common, equally spaced wavelength grid in order to merge together. Flux during resampling is conserved. 

## Output

| Data Type | Content |
| ------------------------ | --------------------------------------------- |
|FITS table |Extracted spectrum and merged orders|

 

:::{bibliography}
:filter: docname in docnames
:::



## Utility API



:::{autodoc2-object} soxspipe.commonutils.horne_extraction.horne_extraction
:::
