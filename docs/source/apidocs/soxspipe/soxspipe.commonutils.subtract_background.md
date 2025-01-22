# {py:mod}`soxspipe.commonutils.subtract_background`

```{py:module} soxspipe.commonutils.subtract_background
```

```{autodoc2-docstring} soxspipe.commonutils.subtract_background
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`subtract_background <soxspipe.commonutils.subtract_background.subtract_background>`
  - ```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background
    :summary:
    ```
````

### API

`````{py:class} subtract_background(log, frame, orderTable, sofName=False, recipeName=False, settings=False, qcTable=False, productsTable=False, lamp='', startNightDate='')
:canonical: soxspipe.commonutils.subtract_background.subtract_background

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background.__init__
```

````{py:method} create_background_image(rowFitOrder, gaussianSigma)
:canonical: soxspipe.commonutils.subtract_background.subtract_background.create_background_image

```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background.create_background_image
```

````

````{py:method} mask_order_locations(orderPixelTable)
:canonical: soxspipe.commonutils.subtract_background.subtract_background.mask_order_locations

```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background.mask_order_locations
```

````

````{py:method} subtract()
:canonical: soxspipe.commonutils.subtract_background.subtract_background.subtract

```{autodoc2-docstring} soxspipe.commonutils.subtract_background.subtract_background.subtract
```

````

`````
