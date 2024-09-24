# soxs_stare

<!-- PURPOSE TEXT -->

## Input

<!-- FIND OBs HERE : https://docs.google.com/spreadsheets/d/1-3VXkIWcydvpawwVl_C3pNTU3HgnElJaYFAKow65Fl8/edit#gid=0 -->

:::{include} inputs/soxs_stare.md
:::


## Parameters

:::{include} parameters/soxs_stare.md
:::

## Method

<!-- METHOD TEXT HERE, FOLLOWED BY WORKFLOW DIAGRAM -->


:::{figure-md} soxs_stare_diagram
:target: soxs_stare.png
![](soxs_stare.png){width=600px}

The `soxs_stare` recipe algorithm.
:::

## Output

:::{include} output/soxs_stare.md
:::


## QC Metrics


:::{include} qcs/soxs_stare.md
:::


## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_stare.soxs_stare
:::
