# {py:mod}`soxspipe.commonutils.horne_extraction`

```{py:module} soxspipe.commonutils.horne_extraction
```

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`horne_extraction <soxspipe.commonutils.horne_extraction.horne_extraction>`
  - ```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`create_cross_dispersion_slice <soxspipe.commonutils.horne_extraction.create_cross_dispersion_slice>`
  - ```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.create_cross_dispersion_slice
    :summary:
    ```
* - {py:obj}`extract_single_order <soxspipe.commonutils.horne_extraction.extract_single_order>`
  - ```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.extract_single_order
    :summary:
    ```
````

### API

````{py:function} create_cross_dispersion_slice(series)
:canonical: soxspipe.commonutils.horne_extraction.create_cross_dispersion_slice

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.create_cross_dispersion_slice
```
````

````{py:function} extract_single_order(crossDispersionSlices, log, ron, slitHalfLength, clippingSigma, clippingIterationLimit, globalClippingSigma, axisA, axisB)
:canonical: soxspipe.commonutils.horne_extraction.extract_single_order

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.extract_single_order
```
````

`````{py:class} horne_extraction(log, settings, recipeSettings, skyModelFrame, skySubtractedFrame, twoDMapPath, recipeName=False, qcTable=False, productsTable=False, dispersionMap=False, sofName=False, locationSetIndex=False)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.__init__
```

````{py:method} extract()
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.extract

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.extract
```

````

````{py:method} merge_extracted_orders(extractedOrdersDF)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.merge_extracted_orders

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.merge_extracted_orders
```

````

````{py:method} plot_extracted_spectrum_qc(uniqueOrders, extractions)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.plot_extracted_spectrum_qc

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.plot_extracted_spectrum_qc
```

````

````{py:method} plot_merged_spectrum_qc(merged_orders)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.plot_merged_spectrum_qc

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.plot_merged_spectrum_qc
```

````

````{py:method} residual_merge(group)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.residual_merge

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.residual_merge
```

````

````{py:method} weighted_average(group)
:canonical: soxspipe.commonutils.horne_extraction.horne_extraction.weighted_average

```{autodoc2-docstring} soxspipe.commonutils.horne_extraction.horne_extraction.weighted_average
```

````

`````
