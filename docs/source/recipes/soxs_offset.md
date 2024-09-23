# soxs_offset

The purpose of the [`soxs_offset`](#soxspipe.recipes.soxs_offset) recipe is to reduced the science frames produced by the NTT and SOXS from a nodding mode observation block.

## Input

<!-- FIND OBs HERE : https://docs.google.com/spreadsheets/d/1-3VXkIWcydvpawwVl_C3pNTU3HgnElJaYFAKow65Fl8/edit#gid=0 -->

:::{include} inputs/soxs_offset.md
:::

## Parameters

## Parameters

:::{include} parameters/soxs_offset.md
:::


## Method

<!-- METHOD TEXT HERE, FOLLOWED BY WORKFLOW DIAGRAM -->

:::{figure-md} soxs_offset_diagram
:target: soxs_offset.png
![](soxs_offset.png){width=600px}

The `soxs_offset` recipe algorithm.
:::

## Output
 
:::{include} output/soxs_offset.md
:::


## QC Metrics

:::{include} qcs/soxs_offset.md
:::

## Recipe API

<!-- :::{autodoc2-object} soxspipe.recipes.soxs_nod.soxs_nod
:::
 -->
