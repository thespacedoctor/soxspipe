## soxs_spatial_solution

The [`soxs_spatial_solution`](../../../recipes/soxs_spatial_solution.md) recipe can be run with the following convention:

```bash
soxspipe [-Vx] spat_sol <inputFrames> [-o <outputDirectory> -s <pathToSettingsFile> --poly=<oowwss>]
```

To rerun a previously executed `soxs_spatial_solution` recipe, you can find the execution command at the end of the [recipe log file](../../logging.md) (found in the workspace `products/soxs_spatial_solution` directory). Use the `-x` flag to overwrite the product files if they already exist. For example, from the root of your workspace, you would run a command like:

```bash
soxspipe spat_sol sof/2021.09.05T15.26.34.694_VIS_1X1_FAST_SPAT_SOL_THAR_10.0S_XSHOOTER.sof -s ./sessions/base/soxspipe.yaml  -x
```

To adjust the default settings for the `soxs_spatial_solution` recipe, open the `soxspipe.yaml` file referenced in the command above in a text editor, navigate to the `soxs_spatial_solution` dictionary, save the file and rerun the recipe command. The settings' descriptions can be found in {numref}`soxs_spatial_solution_parameters`.

:::{include} ../../../recipes/parameters/soxs_spatial_solution.md
:::

Product files are written in the `products/soxs_spatial_solution`, and QC plots are in the `qc/soxs_spatial_solution` workspace directory. A report of the product files, QC plots and metrics is also printed to the terminal. The QC metrics calculated for `soxs_spatial_solution` are found in {numref}`soxs_spatial_solution_qc` and a typical QC plot in {numref}`soxs_spatial_solution_qc_fig`.

:::{include} ../../../recipes/qcs/soxs_spatial_solution.md
:::

:::{figure-md} soxs_spatial_solution_qc_fig
![image-20240924143842700](../../../_images/image-20240924143842700.png){width=600px}

A QC plot resulting from the `soxs_spatial_solution` recipe. The top panel shows an Xshooter VIS arc-lamp frame, taken with a multi-pinhole mask. The green circles represent arc lines detected in the image, and the blue circles and red crosses are lines that were detected but dropped as other pinholes of the same arc line failed to be detected or the lines were clipped during the polynomial fitting. The grey circles represent arc lines reported in the static calibration table that failed to be detected on the image. The second panel shows the same arc-lamp frame with the dispersion solution overlaid as a blue grid. Lines travelling along the dispersion axis (left to right) are lines of equal slit position, and lines travelling in the cross-dispersion direction (top to bottom) are lines of equal wavelength. The third panel shows the residuals of the dispersion solution fit, and the final panel shows the resolution measured for each line (as projected through the pinhole mask) with different colours for each echelle order and the mean order resolution in black.
:::
