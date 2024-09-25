# soxs_stare

The `soxs_stare` recipe reduces object frames taken in stare mode. It models and removes the on-frame sky contribution to the flux. The object trace is then fitted and extracted using an optimal extraction routine.



## Input



:::{include} inputs/soxs_stare.md
:::


## Parameters

:::{include} parameters/soxs_stare.md
:::

## Method

The algorithm used in the `soxs_stare` recipe is shown in {numref}`soxs_stare_diagram`.


:::{figure-md} soxs_stare_diagram
:target: soxs_stare.png
![](soxs_stare.png){width=600px}

The `soxs_stare` recipe algorithm.
:::

If more than one stare mode frame is passed to the `soxs_stare` recipe, there is a call to [`clip_and_stack`](../utils/clip_and_stack.md) to combine the data into a single frame. The single stare-mode frame is detrended using the [`detrend`](../utils/detrend.md), optionally dividing by a master flat field and fitting and removing the background scattered light. The sky-flux is modelled and removed using the `subtract_sky` util, and finally, the object is optimally extracted using the [`horne_extraction`](../utils/horne_extraction.md) utility.

## Output

:::{include} output/soxs_stare.md
:::


## QC Metrics


:::{include} qcs/soxs_stare.md
:::


## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_stare.soxs_stare
:::
