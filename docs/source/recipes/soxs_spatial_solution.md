# soxs_spatial_solution


:::{include} ./descriptions/soxs_spatial_solution.inc
:::


## Input

:::{include} ./inputs/soxs_spatial_solution.inc
:::

:::{include} ./static_files/soxs_spatial_solution.inc
:::



## Parameters

:::{include} parameters/soxs_spatial_solution.inc
:::


## Method

The algorithm used in the `soxs_spatial_solution` recipe is shown in {numref}`soxs_spatial_solution_diagram`.

Having [prepared](../utils/prepare_frames.md) the multi-pinhole frame, the [bias and dark signatures are removed](../utils/detrend.md), and the frame is divided through by the master flat frame. The calibrated frame and the first-guess dispersion map are passed to the [`create_dispersion_map`](../utils/create_dispersion_map.md) utility, arc-lines from a static line-list are detected and fitted to produce a 2D dispersion solution covering both the spectral and spatial dimensions. Note that saturated arc-lines are either absent from the input line-list file or get clipped during the fitting process and so are not included in the final dispersion solution fitting. 

A three-extension FITS image with the same format as the detector is created using the complete dispersion solution. The first extension contains the wavelength values for each pixel, the second the slit positions, and the third the echelle order number. 


:::{figure-md} soxs_spatial_solution_diagram
![](soxs_spatial_solution.png){width=600px}

The `soxs_spatial_solution` recipe algorithm. At the top of the diagram, NIR input data is found on the right and VIS on the left. 
:::

## Output

:::{include} output/soxs_spatial_solution.inc
:::

## QC Metrics

:::{include} qcs/soxs_spatial_solution.inc
:::


:::{figure-md} soxs_spatial_solution_qc
![image-20250130134003533](../_images/image-20250130134003533.png)

A QC plot resulting from the `soxs_spatial_solution` recipe. The top-left panel shows an Xshooter VIS arc-lamp frame, taken with a multi-pinhole mask. The green circles represent arc lines detected in the image, and the blue circles and red crosses are lines that were detected but dropped as other pinholes of the same arc line failed to be detected or the lines were clipped during the polynomial fitting. The grey circles represent arc lines reported in the static calibration table that failed to be detected on the image. The bottom-left panel shows the same arc-lamp frame with the dispersion solution overlaid as a blue grid. Lines travelling along the dispersion axis (left to right) are lines of equal slit position, and lines travelling in the cross-dispersion direction (top to bottom) are lines of equal wavelength. The top-right panel shows the residuals of the dispersion solution fit, and the bottom-right panel shows the resolution measured for each line (as projected through the pinhole mask) with different colours for each echelle order and the mean order resolution in black.
:::




## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution
:::
