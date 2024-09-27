# {py:mod}`soxspipe.commonutils.detect_order_edges`

```{py:module} soxspipe.commonutils.detect_order_edges
```

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`detect_order_edges <soxspipe.commonutils.detect_order_edges.detect_order_edges>`
  - ```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges
    :summary:
    ```
````

### API

`````{py:class} detect_order_edges(log, flatFrame, orderCentreTable, settings=False, recipeSettings=False, recipeName='soxs-mflat', verbose=False, qcTable=False, productsTable=False, tag='', sofName=False, binx=1, biny=1, extendToEdges=True, lampTag=False)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges

Bases: {py:obj}`soxspipe.commonutils._base_detect`

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges.__init__
```

````{py:method} calculate_residuals(orderPixelTable, coeff, axisACol, axisBCol, orderCol=False, writeQCs=False)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.calculate_residuals

````

````{py:method} determine_lower_upper_edge_pixel_positions(orderData)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.determine_lower_upper_edge_pixel_positions

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges.determine_lower_upper_edge_pixel_positions
```

````

````{py:method} determine_order_flux_threshold(orderData, orderPixelTable)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.determine_order_flux_threshold

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges.determine_order_flux_threshold
```

````

````{py:method} fit_global_polynomial(pixelList, axisACol='cont_x', axisBCol='cont_y', orderCol='order', exponentsIncluded=False, writeQCs=False)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.fit_global_polynomial

````

````{py:method} fit_order_polynomial(pixelList, order, axisBDeg, axisACol, axisBCol, exponentsIncluded=False)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.fit_order_polynomial

````

````{py:method} get()
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.get

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges.get
```

````

````{py:method} plot_results(orderPixelTableUpper, orderPixelTableLower, orderPolyTable, orderMetaTable, clippedDataUpper, clippedDataLower)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.plot_results

```{autodoc2-docstring} soxspipe.commonutils.detect_order_edges.detect_order_edges.plot_results
```

````

````{py:method} write_order_table_to_file(frame, orderPolyTable, orderMetaTable)
:canonical: soxspipe.commonutils.detect_order_edges.detect_order_edges.write_order_table_to_file

````

`````
