# {py:mod}`soxspipe.recipes.soxs_nod`

```{py:module} soxspipe.recipes.soxs_nod
```

```{autodoc2-docstring} soxspipe.recipes.soxs_nod
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_nod <soxspipe.recipes.soxs_nod.soxs_nod>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod
    :summary:
    ```
````

### API

`````{py:class} soxs_nod(log, settings=False, inputFrames=[], verbose=False, overwrite=False, command=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.__init__
```

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.detrend

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.get_recipe_settings

````

````{py:method} plot_stacked_spectrum_qc(stackedSpectrum)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.plot_stacked_spectrum_qc

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.plot_stacked_spectrum_qc
```

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.prepare_frames

````

````{py:method} process_single_ab_nodding_cycle(aFrame, bFrame, locationSetIndex, orderTablePath)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.process_single_ab_nodding_cycle

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.process_single_ab_nodding_cycle
```

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.report_output

````

````{py:method} stack_extractions(dataFrameList, postfix='')
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.stack_extractions

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.stack_extractions
```

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_nod.soxs_nod.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_nod.soxs_nod.xsh2soxs

````

`````
