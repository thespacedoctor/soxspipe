# soxs_order_centres

Starting with the first pass dispersion solution from the [`soxs_disp_solution`](#soxspipe.recipes.soxs_disp_solution) recipe, the [`soxs_order_centres`](#soxspipe.recipes.soxs_order_centres) recipe finds and fits a global polynomial model to the central trace of each echelle order.

## Input

:::{include} inputs/soxs_order_centres.md
:::


## Parameters


:::{include} parameters/soxs_order_centres.md
:::



## Method

The algorithm used in the `soxs_order_centres` recipe is shown in {numref}`soxs_order_centres_diagram`.



Once the single-pinhole flat-lamp frame has had the bias, dark and background subtracted it is passed to the [detect_continuum utility](../utils/detect_continuum.md) to fit the order centres.

:::{figure-md} soxs_order_centres_diagram
:target: soxs_order_centres.png
![](soxs_order_centres.png){width=600px}

The `soxs_order_centres` recipe algorithm.
:::

## Output


:::{include} output/soxs_order_centres.md
:::


## QC Metrics

Plots similar to the one below are generated after each execution of [`soxs_order_centres`](#soxspipe.recipes.soxs_order_centres). The residuals of a 'good' fit typically have a mean and standard-deviation <0.2px.


:::{figure-md} soxs_order_centres_qc
:target: ../_images/image-20240924101027298.png
![image-20240924101027298](../_images/image-20240924101027298.png){width=600px}

A QC plot resulting from the `soxs_order_centres` recipe as run on an SOXS NIR single pinhole QTH flat lamp frame. The top panel show the frame with green circles represent the locations on the cross-dispersion slices where a flux peak was detected. The red crosses show the centre of the slices where a peak failed to be detected. The second panel show the global polynomial fitted to the detected order-centre trace with the different colours representing individual echelle orders. The third row of panels show the fit residuals in the X and Y axes. The bottom panel shows the FWHM of the trace fits (in pixels) with respect to echelle order and wavelength.
:::



:::{include} qcs/soxs_order_centres.md
:::

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_order_centres.soxs_order_centres
:::
