# {py:mod}`soxspipe.recipes.soxs_disp_solution`

```{py:module} soxspipe.recipes.soxs_disp_solution
```

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_disp_solution <soxspipe.recipes.soxs_disp_solution.soxs_disp_solution>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.soxs_disp_solution
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`parameterTuning <soxspipe.recipes.soxs_disp_solution.parameterTuning>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.parameterTuning
    :summary:
    ```
````

### API

````{py:function} parameterTuning(p, log, recipeSettings, settings, pinholeFrame, qc, products, sofName, lineDetectionTable)
:canonical: soxspipe.recipes.soxs_disp_solution.parameterTuning

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.parameterTuning
```
````

`````{py:class} soxs_disp_solution(log, settings=False, inputFrames=[], verbose=False, overwrite=False, polyOrders=False, command=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.soxs_disp_solution
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_disp_solution.soxs_disp_solution.xsh2soxs

````

`````
