# {py:mod}`soxspipe.commonutils.subtract_sky`

```{py:module} soxspipe.commonutils.subtract_sky
```

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`subtract_sky <soxspipe.commonutils.subtract_sky.subtract_sky>`
  - ```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky
    :summary:
    ```
````

### API

`````{py:class} subtract_sky(log, settings, recipeSettings, objectFrame, twoDMap, qcTable, productsTable, dispMap=False, sofName=False, recipeName='soxs-stare')
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.__init__
```

````{py:method} add_data_to_placeholder_images(imageMapOrderDF, skymodelCCDData, skySubtractedCCDData, skySubtractedResidualsCCDData)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.add_data_to_placeholder_images

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.add_data_to_placeholder_images
```

````

````{py:method} adjust_tilt(orderDF, tck)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.adjust_tilt

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.adjust_tilt
```

````

````{py:method} calculate_residuals(skyPixelsDF, fluxcoeff, orderDeg, wavelengthDeg, slitDeg, writeQCs=False)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.calculate_residuals

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.calculate_residuals
```

````

````{py:method} clip_object_slit_positions(order_dataframes, aggressive=False)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.clip_object_slit_positions

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.clip_object_slit_positions
```

````

````{py:method} create_placeholder_images()
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.create_placeholder_images

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.create_placeholder_images
```

````

````{py:method} cross_dispersion_flux_normaliser(orderDF)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.cross_dispersion_flux_normaliser

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.cross_dispersion_flux_normaliser
```

````

````{py:method} determine_residual_floor(imageMapOrder, tck)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.determine_residual_floor

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.determine_residual_floor
```

````

````{py:method} fit_bspline_curve_to_sky(imageMapOrder)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.fit_bspline_curve_to_sky

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.fit_bspline_curve_to_sky
```

````

````{py:method} get_over_sampled_sky_from_order(imageMapOrder, clipBPs=True, clipSlitEdge=False)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.get_over_sampled_sky_from_order

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.get_over_sampled_sky_from_order
```

````

````{py:method} plot_image_comparison(objectFrame, skyModelFrame, skySubFrame)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.plot_image_comparison

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.plot_image_comparison
```

````

````{py:method} plot_sky_sampling(order, imageMapOrderDF, tck=False, knotLocations=False)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.plot_sky_sampling

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.plot_sky_sampling
```

````

````{py:method} rolling_window_clipping(imageMapOrderDF, windowSize, sigma_clip_limit=5, max_iterations=10)
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.rolling_window_clipping

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.rolling_window_clipping
```

````

````{py:method} subtract()
:canonical: soxspipe.commonutils.subtract_sky.subtract_sky.subtract

```{autodoc2-docstring} soxspipe.commonutils.subtract_sky.subtract_sky.subtract
```

````

`````
