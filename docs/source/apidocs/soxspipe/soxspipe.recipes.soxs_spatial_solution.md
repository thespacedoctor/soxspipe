# {py:mod}`soxspipe.recipes.soxs_spatial_solution`

```{py:module} soxspipe.recipes.soxs_spatial_solution
```

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_spatial_solution <soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`parameterTuning <soxspipe.recipes.soxs_spatial_solution.parameterTuning>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.parameterTuning
    :summary:
    ```
````

### API

````{py:function} parameterTuning(p, log, recipeSettings, settings, multiPinholeFrame, disp_map_table, order_table, qc, products, sofName, lineDetectionTable)
:canonical: soxspipe.recipes.soxs_spatial_solution.parameterTuning

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.parameterTuning
```
````

`````{py:class} soxs_spatial_solution(log, settings=False, inputFrames=[], verbose=False, overwrite=False, create2DMap=True, polyOrders=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_spatial_solution.soxs_spatial_solution.xsh2soxs

````

`````
