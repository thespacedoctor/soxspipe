# {py:mod}`soxspipe.recipes.soxs_order_centres`

```{py:module} soxspipe.recipes.soxs_order_centres
```

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_order_centres <soxspipe.recipes.soxs_order_centres.soxs_order_centres>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.soxs_order_centres
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`parameterTuning <soxspipe.recipes.soxs_order_centres.parameterTuning>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.parameterTuning
    :summary:
    ```
````

### API

````{py:function} parameterTuning(p, log, recipeSettings, settings, orderFrame, disp_map_table, orderPixelTable, qc, products, sofName, binx, biny)
:canonical: soxspipe.recipes.soxs_order_centres.parameterTuning

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.parameterTuning
```
````

`````{py:class} soxs_order_centres(log, settings=False, inputFrames=[], verbose=False, overwrite=False, polyOrders=False, command=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.soxs_order_centres
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.soxs_order_centres.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.soxs_order_centres.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_order_centres.soxs_order_centres.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_order_centres.soxs_order_centres.xsh2soxs

````

`````
