# soxs_mbias

:::{include} ./descriptions/soxs_mbias.inc
:::


## Input

:::{include} ./inputs/soxs_mbias.inc
:::

:::{include} ./static_files/soxs_mbias.inc
:::


## Parameters

:::{include} ./parameters/soxs_mbias.inc
:::

## Method

The purpose of the [`soxs_mbias`](#soxspipe.recipes.soxs_mbias) recipe is to stack raw bias frames together (using the [`clip_and_stack`](../utils/clip_and_stack.md) utility) into master-bias frames, clipping rogue pixels from the individual raw frames and reducing the read-noise contribution. The algorithm used in the [`soxs_mbias`](#soxspipe.recipes.soxs_mbias) recipe is shown in {numref}`soxs_mbias_diagram`.

:::{figure-md} soxs_mbias_diagram
![](soxs_mbias.png){width=600px}

The `soxs-mbias` recipe algorithm.
:::

## Output

:::{include} ./output/soxs_mbias.inc
:::


## QC Metrics

:::{include} ./qcs/soxs_mbias.inc
:::


## Recipe API

:::{autodoc2-object} ./soxspipe.recipes.soxs_mbias.soxs_mbias
:::

