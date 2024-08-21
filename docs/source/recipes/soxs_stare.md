# `soxs_stare` - PLANNED

<!-- PURPOSE TEXT -->

## Input

<!-- FIND OBs HERE : https://docs.google.com/spreadsheets/d/1-3VXkIWcydvpawwVl_C3pNTU3HgnElJaYFAKow65Fl8/edit#gid=0 -->

| Data Type | Content | Related OB |
|:----|:----|:---|
| | |

## Parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| stacked-clipping-sigma | number of σ deviations from the median *pixel* flux beyond which pixel is excluded from stack | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-iterations | number of σ-clipping iterations to perform before stacking | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
|   |   |   |   |

## Method

<!-- METHOD TEXT HERE, FOLLOWED BY WORKFLOW DIAGRAM -->

<!-- ![](soxs_stare.png) -->

## Output
 
| Data Type | Content |
|:----|:----|
| |

## QC Metrics

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_stare.soxs_stare
:::
