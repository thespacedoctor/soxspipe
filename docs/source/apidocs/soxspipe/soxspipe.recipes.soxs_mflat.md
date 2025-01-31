# {py:mod}`soxspipe.recipes.soxs_mflat`

```{py:module} soxspipe.recipes.soxs_mflat
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`soxs_mflat <soxspipe.recipes.soxs_mflat.soxs_mflat>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`nearest_neighbour <soxspipe.recipes.soxs_mflat.nearest_neighbour>`
  - ```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.nearest_neighbour
    :summary:
    ```
````

### API

````{py:function} nearest_neighbour(singleValue, listOfValues)
:canonical: soxspipe.recipes.soxs_mflat.nearest_neighbour

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.nearest_neighbour
```
````

`````{py:class} soxs_mflat(log, settings=False, inputFrames=[], verbose=False, overwrite=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat

Bases: {py:obj}`soxspipe.recipes.base_recipe.base_recipe`

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.__init__
```

````{py:method} calibrate_frame_set()
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.calibrate_frame_set

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.calibrate_frame_set
```

````

````{py:method} clean_up()
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.clean_up

````

````{py:method} clip_and_stack(frames, recipe, ignore_input_masks=False, post_stack_clipping=True)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.clip_and_stack

````

````{py:method} detrend(inputFrame, master_bias=False, dark=False, master_flat=False, order_table=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.detrend

````

````{py:method} find_uvb_overlap_order_and_scale(dcalibratedFlats, qcalibratedFlats)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.find_uvb_overlap_order_and_scale

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.find_uvb_overlap_order_and_scale
```

````

````{py:method} get_recipe_settings()
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.get_recipe_settings

````

````{py:method} mask_low_sens_pixels(frame, orderTablePath, returnMedianOrderFlux=False, writeQC=True)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.mask_low_sens_pixels

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.mask_low_sens_pixels
```

````

````{py:method} normalise_flats(inputFlats, orderTablePath, firstPassMasterFlat=False, lamp='')
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.normalise_flats

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.normalise_flats
```

````

````{py:method} prepare_frames(save=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.prepare_frames

````

````{py:method} produce_product()
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.produce_product

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.produce_product
```

````

````{py:method} qc_median_flux_level(frame, frameType='MBIAS', frameName='master bias', medianFlux=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.qc_median_flux_level

````

````{py:method} qc_ron(frameType=False, frameName=False, masterFrame=False, rawRon=False, masterRon=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.qc_ron

````

````{py:method} report_output(rformat='stdout')
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.report_output

````

````{py:method} stitch_uv_mflats(medianOrderFluxDF, orderTablePath)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.stitch_uv_mflats

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.stitch_uv_mflats
```

````

````{py:method} subtract_mean_flux_level(rawFrame)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.subtract_mean_flux_level

````

````{py:method} update_fits_keywords(frame, rawFrames=False)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.update_fits_keywords

````

````{py:method} verify_input_frames()
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.verify_input_frames

```{autodoc2-docstring} soxspipe.recipes.soxs_mflat.soxs_mflat.verify_input_frames
```

````

````{py:method} xsh2soxs(frame)
:canonical: soxspipe.recipes.soxs_mflat.soxs_mflat.xsh2soxs

````

`````
