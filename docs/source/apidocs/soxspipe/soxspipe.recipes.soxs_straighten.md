# {py:mod}`soxspipe.recipes.soxs_straighten`

```{py:module} soxspipe.recipes.soxs_straighten
```

```{autodoc2-docstring} soxspipe.recipes.soxs_straighten
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_straighten <soxspipe.recipes.soxs_straighten.soxs_straighten>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_straighten.soxs_straighten
    :summary:
    ```
````

### API

`````{py:class} soxs_straighten(log, settings=False, inputFrames=[], verbose=False, overwrite=False)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten

Bases: {py:obj}`soxspipe.recipes._base_recipe_._base_recipe_`

```{autodoc2-docstring} soxspipe.recipes.soxs_straighten.soxs_straighten
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_straighten.soxs_straighten.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_straighten.soxs_straighten.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.report_output

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_straighten.soxs_straighten.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_straighten.soxs_straighten.xsh2soxs

````

`````
