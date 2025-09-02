# {py:mod}`soxspipe.recipes.soxs_mdark`

```{py:module} soxspipe.recipes.soxs_mdark
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mdark
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_mdark <soxspipe.recipes.soxs_mdark.soxs_mdark>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_mdark.soxs_mdark
    :summary:
    ```
````

### API

`````{py:class} soxs_mdark(log, settings=False, inputFrames=[], verbose=False, overwrite=False, command=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_mdark.soxs_mdark
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mdark.soxs_mdark.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_mdark.soxs_mdark.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_mdark.soxs_mdark.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_mdark.soxs_mdark.xsh2soxs

````

`````
