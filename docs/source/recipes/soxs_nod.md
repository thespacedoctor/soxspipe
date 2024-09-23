# soxs_nod

The purpose of the [`soxs_nod`](#soxspipe.recipes.soxs_nod) recipe is to reduced the science frames produced by the NTT and SOXS from a nodding mode observation block.

## Input

<!-- FIND OBs HERE : https://docs.google.com/spreadsheets/d/1-3VXkIWcydvpawwVl_C3pNTU3HgnElJaYFAKow65Fl8/edit#gid=0 -->


:::{include} inputs/soxs_nod.md
:::

## Parameters

## Parameters

:::{include} parameters/soxs_nod.md
:::


## Method

<!-- METHOD TEXT HERE, FOLLOWED BY WORKFLOW DIAGRAM -->

:::{figure-md} soxs_nod_diagram
:target: soxs_nod.png
![](soxs_nod.png){width=600px}

The `soxs_nod` recipe algorithm.
:::

## Output
 
:::{include} output/soxs_nod.md
:::

## QC Metrics


:::{include} qcs/soxs_nod.md
:::

## Recipe API

<!-- :::{autodoc2-object} soxspipe.recipes.soxs_nod.soxs_nod
:::
 -->
