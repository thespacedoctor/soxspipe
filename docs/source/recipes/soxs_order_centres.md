# soxs_order_centres

The purpose of the [`soxs_order_centres`](#soxspipe.recipes.soxs_order_centres) recipe is to find and fit the order centres with low-level polynomials.

## Input

:::{include} inputs/soxs_order_centres.md
:::


## Parameters


:::{include} parameters/soxs_order_centres.md
:::



## Method

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

Plots similar to the one below are generated after each execution of [`soxs_order_centres`](#soxspipe.recipes.soxs_order_centres).

[![](https://live.staticflickr.com/65535/50345130012_4e869a6a7f_b.png)](https://live.staticflickr.com/65535/50345130012_4e869a6a7f_o.png)

:::{include} qcs/soxs_order_centres.md
:::

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_order_centres.soxs_order_centres
:::
