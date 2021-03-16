## `soxs_spatial_solution` - PLANNED

The purpose of this recipe is to further enhance the wavelength solution achieved with `soxs_disp_solution` by expanding the solution into the spatial dimension (along the slit). This 2-dimensional solution will then account for any tilt in the spectral lines.[^20210305045943]

Each pinhole in the multi-pinhole mask is $\text{0.5"}$ in diameter and the 9 pinholes are evenly spaced along the $\text{11"}$ slit with a $\text{1.4"}$  gap between adjacent holes. This knowledge affords us the ability to now map the dispersion solution along the spatial direction.

As input the recipe receives the [Dispersion Map](../files/dispersion_map.md) table from `soxs_disp_solution`.

In the static calibration suite we have tables listing the wavelength $\lambda$, order number $n$ and slit position $s$ of the spectral lines alongside a first approximation of their ($X, Y$) pixel-positions on the detector (the same one used by `soxs_disp_solution`).

## Input

| Data Type | Content | Related OB |
|:----|:----|:---|
| Raw Calibration | Arc lamp frame with a multi-pinhole mask, Dark frame of equal exposure time as arc frame (NIR only) | `SOXS_slt_cal_VISArcsMultiplePinhole`, `SOXS_slt_cal_NIRArcsMultiplePinhole` |
| Intermediate Calibration | [Dispersion Map](../files/dispersion_map.md) table from `soxs_disp_solution`. Master bias and master dark frame (UV-VIS only). Master flat-frame. |  |
| Static Calibration | Multi-pinhole detector position map | |
| Raw Science | - | |

## Method

Having [prepared](../utils/prepare_frames.md) the multi-pinhole frame the [bias and dark signatures are removed](../utils/subtract_calibrations.md) and the frame is divided through by the master flat frame. The calibrated frame and the first-guess dispersion map are passed to the [`create_dispersion_map`](../utils/create_dispersion_map.md) utility to produce a 2D dispersion solution covering both the spectral and spatial dimensions.

![](soxs_spatial_solution.png)

## Output

* [Dispersion Map](../files/dispersion_map.md) table giving 2D dispersion solution.

## QC Metrics

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |


```eval_rst
.. autoclass:: soxspipe.recipes.soxs_spatial_solution
    :members:
```


[^20210305045943]: relative to the perpendicular of the dispersion direction
