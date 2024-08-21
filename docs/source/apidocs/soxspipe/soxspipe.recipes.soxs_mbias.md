# {py:mod}`soxspipe.recipes.soxs_mbias`

```{py:module} soxspipe.recipes.soxs_mbias
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_mbias <soxspipe.recipes.soxs_mbias.soxs_mbias>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias
    :summary:
    ```
````

### API

`````{py:class} soxs_mbias(log, settings=False, inputFrames=[], verbose=False, overwrite=False)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias

Bases: {py:obj}`soxspipe.recipes._base_recipe_._base_recipe_`

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias.produce_product
```

````

````{py:method} qc_bias_structure(combined_bias_mean)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.qc_bias_structure

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias.qc_bias_structure
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.qc_median_flux_level

````

````{py:method} qc_periodic_pattern_noise(frames)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.qc_periodic_pattern_noise

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias.qc_periodic_pattern_noise
```

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_mbias.soxs_mbias.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_mbias.soxs_mbias.xsh2soxs

````

`````
