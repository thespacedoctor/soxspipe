# `horne_extraction` 

The purpose of the [`horne_extraction`](#soxspipe.commonutils.horne_extraction) utility is to perform optimal extraction on each spectral order applying the algorithms reported in Horne+86.

The typical execution workflow is the following:

![](horne_extraction.png)

This utility, according to the prescrption of Horne+86, follows those basic steps:

**1 Determination of the spatial profile

Starting from the order_table computed by the `detect_continumm` utility, the horne_extraction utility runs along the spectral direction and takes, for each wavelenght centered in the position measured by the detect_continuum utility, a window of 1x`horne-extraction-slit-length`. Then, the pixels are summed and each pixel in the slice is normalized as 

$$
pixel_i = \frac{Flux_{pixel_i}}{\sum_{j=0}^{N}{Flux_{pixel_j}}}
$$ 

where N = `horne-extraction-slit-length`.


Then, for each pixel in the slice along the dispersion a low order polynomial is fitted. This polynomial models the value of the fractional flux received in this pixel by the object along the dispersion. The full set of polynomials represent the profile of the object along the spatial direction for each wavelenght. For a slice length of N pixels there will be N different polynomials, one per each pixel.

Polynomials are fitted with subsequent iterations. Pixels for which they residuals deviate more than  `horne-extraction-profile-clipping-sigma` standard deviation are removed before attemping a new fitting of the data. The procedure is repeated for a maximum of `horne-extraction-profile-clipping-iteration-count` times.

When the procedure above is completed, the actual extraction takes place as follows:

**2 Extraction of the spectrum for each wavelength

Using the polynomials computed above, the extracted integrated flux for each wavelenght in the order is computed as

$$
extractedFlux_{\lambda} = \frac{\sum_{i = 0}^{N}{P_{i}^{\lambda} (D_{i}^{\lambda} - S_{i}^{\lambda})}}{\sum_{i = 0}^{N}{(P_{i}^{\lambda})^2 (V_{i}^{\lambda})^{-1}}}
$$

where $P_{i}^{\lambda}$ is the value of the polynomial for the pixel $i^{th}$ evaluated at wavelength $\lambda$, $D_{i}^{\lambda}$ is the raw pixel value on the 2D image, detrended and sky subtracted, of the science object and $S_{i}^{\lambda}$ is the estimated value, at the same pixel position, of the 2D image of the sky model computed in stare mode. When this utility is applied on nodding or offset mode images, this value is chosen to be zero. $V_{i}^{\lambda}$ is the so called variance image and it computed as 

$$
V_{i}^{\lambda} = V_{0} + \frac{|P_{i}{\lambda} f_{\lambda} +  S_{i}{\lambda}|{Q}
$$
where $V_{0}$ is the square of the readout noise of the detector, $f_\{lambda}$ is the sum of the pixels within the extraction window at each wavelength and Q is the detector gain. 

Cosmic rays and uncorrected deviant pixels are removed from the extraction by comparing their value with a certain threshold. In particular, for each pixel $i$ in the extraction window the quantity

$$
\frac{(D_{i}^{\lambda} - P_{i}{\lambda} f_{\lambda} - S_{i}{\lambda})^2}{V_{i}^{\lambda}}
$$

is computed. Pixels for which this quantity exceeds the `horne-extraction-profile-global-clipping-sigma` are not included in the extraction.
 