## `create_dispersion_map`

The [`create_dispersion_map`](#soxspipe.commonutils.create_dispersion_map) utility is used to search for arc-lines in the single/multi-pinhole arc-lamp frames and then iteratively fit a global polynomial dispersion solution (and spatial-solution in the case of multi-pinhole frame) with the observed line-positions. It is used by both the [`soxs_disp_solution`](../recipes/soxs_disp_solution.md)) and [`soxs_spatial_solution`](../recipes/soxs_spatial_solution.md)) solution recipes.

![](create_dispersion_map.png)

In the static calibration suite we have [Pinhole Maps](../files/pinhole_map.md) listing the wavelength $\lambda$, order number $n$ and slit position $s$ of the spectral lines alongside a first approximation of their ($X, Y$) pixel-positions on the detector.

If the input frame is a single-pinhole frame, we can filter the Pinhole Map to contain just the central pinhole positions. If however input is the multi-pinhole frame then we use the first guess Dispersion Map (created with [`soxs_disp_solution`](../recipes/soxs_disp_solution.md)) to calculate the shift between the predicted and the observed line positions for the central pinholes. We then update the Pinhole Map by applying the same shift to the other pinholes.

For each line in the Pinhole Map line-list:

* an image stamp centred on the predicted pixel-position ($X_o, Y_o$), of dimensions winX and winY, is generated from the pinhole calibration frame
* a sigma-clipped median pixel value is calculated and then subtracted from each stamp 
* DAOStarFinder is employed to search for the *observed* detector position ($X, Y$) of the arc-line via 2D Gaussian profile fitting on the stamp

We now have a list of arc-line wavelengths and their observed pixel-positions and the order they were detected in. These values are used to iteratively fit two polynomials that describe the global dispersion solution for the detector. In the case of the single-pinhole frames these are:

$$X = \sum\limits_{ij} c_{ij} \times n^i \times \lambda^j$$

$$Y = \sum\limits_{ij} c_{ij} \times n^i \times \lambda^j$$

where $\lambda$ is wavelength and $n$ is the echelle order number.

In the case of the multi-pinhole we also have the slit position $s$ and so adding a spatial solution to the dispersion solution:

$$X = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k$$

$$Y = \sum\limits_{ijk} c_{ijk} \times n^i \times \lambda^j \times s^k$$

Upon each iteration the residuals between the fits and the measured pixel-positions are calculated and sigma-clipping is employed to eliminate measurements that stray to far from the fit. Once the maximum number of iterations is reach, or all outlying lines have been clipped, the coefficients of the polynomials are written to a [Dispersion Map](../files/dispersion_map.md) file.

### 2D Image Map

[![](https://live.staticflickr.com/65535/51862169299_f6773a5b0f_b.jpg)](https://live.staticflickr.com/65535/51862169299_f6773a5b0f_b.jpg)


The [Dispersion Map](../files/dispersion_map.md) is used to generate a triple extension FITS file with each extension image exactly matching the dimensions of the detector. The first extension contains the wavelength value at the centre of each pixel location, the second the slit-position and the third the order number. The solutions for these images are iteratively converged on in a brute force manner (see workflow diagram below). These image maps are used in sky-background subtraction and object extraction utilities. 


![](create_dispersion_map_to_image.png)


:::{autodoc2-object} soxspipe.commonutils.create_dispersion_map.create_dispersion_map
:::


