## soxs_disp_solution


:::{include} ../../../recipes/descriptions/soxs_disp_solution.inc
:::

### Usage

:::{include} ../../../recipes/cl_usage/soxs_disp_solution.inc
:::


### Parameters

:::{include} ../../../recipes/parameters/soxs_disp_solution.inc
:::

### Input

:::{include} ../../../recipes/inputs/soxs_disp_solution.inc
:::

### Output

:::{include} ../../../recipes/output/soxs_disp_solution.inc
:::

### QC Metrics

:::{include} ../../../recipes/qcs/soxs_disp_solution.inc
:::


The typical solution for the `soxs_disp_solution` recipe has sub-pixel residuals.


:::{figure-md} soxs_disp_solution_qc
![image-20250124164511814](../../../_images/image-20250124164511814.png)

A QC plot resulting from the `soxs_disp_solution` recipe as run on a SOXS NIR single pinhole arc lamp frame. A 'good' dispersion solution will have sub-pixel residuals (mean residuals $<$ 0.5 pixels). The top-left panel shows an SOXS NIR arc-lamp frame, taken with a single pinhole mask. The green circles represent arc lines detected in the image, and the blue circles and red crosses were detected but dropped due to poor DAOStarFinder fitting or clipped during the polynomial fitting, respectively. The grey circles represent arc lines reported in the static calibration table that failed to be detected on the image.  The bottom-left panel shows the same arc-lamp frame with the dispersion solution overlaid at the pixel locations modelled for the original lines in the line list. The top-right panels shows the residuals of the dispersion solution fit, and the final panels (bottom-right) the resolution measured for each line (as projected through the pinhole mask) with different colours for each echelle order and the mean order resolution in black.

:::
