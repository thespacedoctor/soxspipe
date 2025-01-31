# {py:mod}`soxspipe.recipes.base_recipe`

```{py:module} soxspipe.recipes.base_recipe
```

```{autodoc2-docstring} soxspipe.recipes.base_recipe
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`base_recipe <soxspipe.recipes.base_recipe.base_recipe>`
  - ```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe
    :summary:
    ```
````

### API

`````{py:class} base_recipe(log, settings=False, inputFrames=False, verbose=False, overwrite=False, recipeName=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.base_recipe.base_recipe.clean_up

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.clean_up
```

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.base_recipe.base_recipe.clip_and_stack

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.clip_and_stack
```

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe.detrend

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.detrend
```

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.base_recipe.base_recipe.get_recipe_settings

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.get_recipe_settings
```

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe.prepare_frames

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.prepare_frames
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe.qc_median_flux_level

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.qc_median_flux_level
```

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe.qc_ron

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.qc_ron
```

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.base_recipe.base_recipe.report_output

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.report_output
```

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.base_recipe.base_recipe.subtract_mean_flux_level

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.subtract_mean_flux_level
```

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.base_recipe.base_recipe.update_fits_keywords

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.update_fits_keywords
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.base_recipe.base_recipe.xsh2soxs

```{autodoc2-docstring} soxspipe.recipes.base_recipe.base_recipe.xsh2soxs
```

````

`````
