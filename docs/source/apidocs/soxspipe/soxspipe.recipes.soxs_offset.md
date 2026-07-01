# {py:mod}`soxspipe.recipes.soxs_offset`

```{py:module} soxspipe.recipes.soxs_offset
```

```{autodoc2-docstring} soxspipe.recipes.soxs_offset
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_offset <soxspipe.recipes.soxs_offset.soxs_offset>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_offset.soxs_offset
    :summary:
    ```
````

### API

`````{py:class} soxs_offset(log, settings=False, inputFrames=[], verbose=False, overwrite=False, command=False, debug=False, turnOffMP=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset

Bases: {py:obj}`soxspipe.recipes.soxs_nod.soxs_nod`

```{autodoc2-docstring} soxspipe.recipes.soxs_offset.soxs_offset
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_offset.soxs_offset.__init__
```

````{py:method} clean_up(forceFail=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.detrend

````

````{py:method} flag_poor_data()
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.flag_poor_data

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.get_recipe_settings

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.prepare_frames

````

````{py:method} process_single_ab_nodding_cycle(aFrame, bFrame, locationSetIndex, orderTablePath, notFlattened=False, masterFlat=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.process_single_ab_nodding_cycle

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_offset.soxs_offset.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.report_output

````

````{py:method} stack_extractions(dataFrameList, notFlattened=False, orderJoins=None)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.stack_extractions

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.verify_input_frames

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_offset.soxs_offset.xsh2soxs

````

`````
