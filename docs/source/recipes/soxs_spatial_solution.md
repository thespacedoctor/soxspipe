# soxs_spatial_solution

The ` soxs_spatial_solution` recipe enhances the wavelength solution achieved with [`soxs_disp_solution`](../recipes/soxs_disp_solution.md) by expanding the solution into the spatial dimension (along the slit). The resulting 2-dimensional solution accounts for any tilt in the spectral lines relative to the cross-dispersion axis. The recipe is similar in logic to [`soxs_disp_solution`](../recipes/soxs_disp_solution.md), but now samples arc lines along the slit in the cross-dispersion direction using a multi-pinhole slit mask.

## Input

:::{include} inputs/soxs_spatial_solution.md
:::


## Parameters

:::{include} parameters/soxs_spatial_solution.md
:::


## Method

The algorithm used in the `soxs_spatial_solution` recipe is shown in {numref}`soxs_spatial_solution_diagram`.

Having [prepared](../utils/prepare_frames.md) the multi-pinhole frame the [bias and dark signatures are removed](../utils/detrend.md) and the frame is divided through by the master flat frame. The calibrated frame and the first-guess dispersion map are passed to the [`create_dispersion_map`](../utils/create_dispersion_map.md) utility to produce a 2D dispersion solution covering both the spectral and spatial dimensions.

Using the complete dispersion solution, a three-extension FITS image with the same format as the detector is created. The first extension contains the wavelength values for each pixel, the second the slit positions, and the third the echelle order number.


:::{figure-md} soxs_spatial_solution_diagram
:target: soxs_spatial_solution.png
![](soxs_spatial_solution.png){width=600px}

The `soxs_spatial_solution` recipe algorithm.
:::

## Output

:::{include} output/soxs_spatial_solution.md
:::

## QC Metrics

[![](https://live.staticflickr.com/65535/51171156692_0588cc30d6_z.png)](https://live.staticflickr.com/65535/51171156692_0588cc30d6_o.png)

:::{include} qcs/soxs_spatial_solution.md
:::

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution
:::
