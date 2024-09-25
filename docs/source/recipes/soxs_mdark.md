# soxs_mdark

Every raw CCD image contains counts resulting from a 'dark current', electrons released due to the thermal effects in the CCD material. For both the UVB-VIS (< 0.00012 $\rm{e^{-}/s/pixel}$) and NIR detectors (< 0.005  $\rm{e^{-}/s/pixel}$), the dark current is almost negligible. Not all pixels will have the same dark current, with some pixels having a higher-than-typical current. These are so-called 'hot-pixels', and these must be identified and recorded (using the [`create_noise_map`](../utils/create_noise_map.md) utility).

The purpose of the [`soxs_mdark`](#soxspipe.recipes.soxs_mdark) recipe is to generate a master-dark frame used to remove flux attributed to the dark current from other frames.

## Input


:::{include} inputs/soxs_mdark.md
:::

## Parameters

:::{include} parameters/soxs_mdark.md
:::


## Method

The raw dark frames are stacked together (using the [`clip_and_stack`](../utils/clip_and_stack.md) utility) into master-dark frames and, in the process, clipping rogue pixels from the individual raw frames and reducing the read-noise contribution. The algorithm used in the `soxs_mdark` recipe is shown in {numref}`soxs_mdark_diagram`.

:::{figure-md} soxs_mdark_diagram
![](../_images/soxs_mdark.png){width=600px}

The soxs-mdark recipe algorithm.
:::

## Output

:::{include} output/soxs_mdark.md
:::



## QC Metrics


:::{include} qcs/soxs_mdark.md
:::


## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_mdark.soxs_mdark
:::
