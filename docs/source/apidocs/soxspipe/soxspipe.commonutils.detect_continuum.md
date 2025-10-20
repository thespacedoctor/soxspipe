# {py:mod}`soxspipe.commonutils.detect_continuum`

```{py:module} soxspipe.commonutils.detect_continuum
```

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`detect_continuum <soxspipe.commonutils.detect_continuum.detect_continuum>`
  - ```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum
    :summary:
    ```
````

### API

`````{py:class} detect_continuum(log, traceFrame, dispersion_map, settings=False, recipeSettings=False, recipeName=False, qcTable=False, productsTable=False, sofName=False, binx=1, biny=1, lampTag=False, locationSetIndex=False, orderPixelTable=False, startNightDate='')
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum

Bases: {py:obj}`soxspipe.commonutils.detect_continuum._base_detect`

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.__init__
```

````{py:method} calculate_residuals(orderPixelTable, coeff, axisACol, axisBCol, orderCol=False, writeQCs=False)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.calculate_residuals

````

````{py:method} create_pixel_arrays()
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.create_pixel_arrays

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.create_pixel_arrays
```

````

````{py:method} fit_1d_gaussian_to_slices(orderPixelTable, sliceLength, medianStddev=False)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.fit_1d_gaussian_to_slices

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.fit_1d_gaussian_to_slices
```

````

````{py:method} fit_global_polynomial(pixelList, axisACol='cont_x', axisBCol='cont_y', orderCol='order', exponentsIncluded=False, writeQCs=False)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.fit_global_polynomial

````

````{py:method} fit_order_polynomial(pixelList, order, axisBDeg, axisACol, axisBCol, exponentsIncluded=False)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.fit_order_polynomial

````

````{py:method} get()
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.get

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.get
```

````

````{py:method} plot_results(orderPixelTable, orderPolyTable, clippedData)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.plot_results

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.plot_results
```

````

````{py:method} sample_trace()
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.sample_trace

```{autodoc2-docstring} soxspipe.commonutils.detect_continuum.detect_continuum.sample_trace
```

````

````{py:method} write_order_table_to_file(frame, orderPolyTable, orderMetaTable)
:canonical: soxspipe.commonutils.detect_continuum.detect_continuum.write_order_table_to_file

````

`````
