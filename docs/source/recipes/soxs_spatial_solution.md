# `soxs_spatial_solution` - PLANNED

The purpose of this recipe is to further enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit). This 2-dimensional solution will then account for any tilt in the spectral lines.[^20210305045943]

Each pinhole in the multi-pinhole mask is $\text{0.5"}$ in diameter and the 9 pinholes are evenly spaced along the $\text{11"}$ slit with a $\text{1.4"}$  gap between adjacent holes. This knowledge affords us the ability to now map the dispersion solution along the spatial direction.

## Input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Arc Lamp through multi-pinhole mask | `SOXS_slt_cal_VISArcsMultiplePinhole`, `SOXS_slt_cal_NIRArcsMultiplePinhole` |
| FITS Image | Master Dark Frame (VIS only) | - |
| FITS Image | Master Bias Frame (VIS only) | - |
| FITS Image | Dark frame (Lamp-Off) of equal exposure length as multi-pinhole frame (Lamp-On) (NIR only) | `SOXS_slt_cal_NIRArcsMultiplePinhole` |
| CSV File | First-guess [Dispersion Map](../files/dispersion_map.md) table |
| CSV File | [Pinhole Map](../files/pinhole_map.md) |

## Method

Having [prepared](../utils/prepare_frames.md) the multi-pinhole frame the [bias and dark signatures are removed](../utils/detrend.md) and the frame is divided through by the master flat frame. The calibrated frame and the first-guess dispersion map are passed to the [`create_dispersion_map`](../utils/create_dispersion_map.md) utility to produce a 2D dispersion solution covering both the spectral and spatial dimensions.

![](soxs_spatial_solution.png)

## Output

| Data Type | Content |
|:----|:----|
| CSV File (subject to change) | [Dispersion Map](../files/dispersion_map.md) table giving coefficients of polynomials describing 2D dispersion/spatial solution |

## QC Metrics

[![](https://live.staticflickr.com/65535/51171156692_0588cc30d6_z.png)](https://live.staticflickr.com/65535/51171156692_0588cc30d6_o.png)

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution
:::


[^20210305045943]: relative to the perpendicular of the dispersion direction
