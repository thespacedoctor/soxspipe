# {py:mod}`soxspipe.commonutils.create_dispersion_map`

```{py:module} soxspipe.commonutils.create_dispersion_map
```

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`create_dispersion_map <soxspipe.commonutils.create_dispersion_map.create_dispersion_map>`
  - ```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`measure_line_position <soxspipe.commonutils.create_dispersion_map.measure_line_position>`
  - ```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.measure_line_position
    :summary:
    ```
* - {py:obj}`straighten_mph_sets <soxspipe.commonutils.create_dispersion_map.straighten_mph_sets>`
  - ```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.straighten_mph_sets
    :summary:
    ```
````

### API

`````{py:class} create_dispersion_map(log, settings, recipeSettings, pinholeFrame, firstGuessMap=False, orderTable=False, qcTable=False, productsTable=False, sofName=False, create2DMap=True, lineDetectionTable=False, startNightDate='', arcFrame=None)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.__init__
```

````{py:method} calculate_residuals(orderPixelTable, xcoeff, ycoeff, orderDeg, wavelengthDeg, slitDeg, writeQCs=False, pixelRange=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.calculate_residuals

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.calculate_residuals
```

````

````{py:method} convert_and_fit(order, bigWlArray, bigSlitArray, slitMap, wlMap, iteration, plots=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.convert_and_fit

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.convert_and_fit
```

````

````{py:method} create_new_static_line_list(dispersionMapPath)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.create_new_static_line_list

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.create_new_static_line_list
```

````

````{py:method} create_placeholder_images(order=False, plot=False, reverse=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.create_placeholder_images

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.create_placeholder_images
```

````

````{py:method} detect_pinhole_arc_lines(orderPixelTable, iraf=True, sigmaLimit=3, iteration=False, brightest=False, exclude_border=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.detect_pinhole_arc_lines

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.detect_pinhole_arc_lines
```

````

````{py:method} fit_polynomials(orderPixelTable, wavelengthDeg, orderDeg, slitDeg, missingLines=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.fit_polynomials

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.fit_polynomials
```

````

````{py:method} get()
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.get

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.get
```

````

````{py:method} get_predicted_line_list()
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.get_predicted_line_list

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.get_predicted_line_list
```

````

````{py:method} map_to_image(dispersionMapPath, orders=False)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.map_to_image

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.map_to_image
```

````

````{py:method} order_to_image(orderInfo)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.order_to_image

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.order_to_image
```

````

````{py:method} update_static_line_list_detector_positions(originalOrderPixelTable, dispersionMapPath)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.update_static_line_list_detector_positions

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.update_static_line_list_detector_positions
```

````

````{py:method} write_map_to_file(xcoeff, ycoeff, orderDeg, wavelengthDeg, slitDeg)
:canonical: soxspipe.commonutils.create_dispersion_map.create_dispersion_map.write_map_to_file

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.create_dispersion_map.write_map_to_file
```

````

`````

````{py:function} measure_line_position(stampInfo, log, windowHalf, iraf, sigmaLimit, brightest=False, exclude_border=False)
:canonical: soxspipe.commonutils.create_dispersion_map.measure_line_position

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.measure_line_position
```
````

````{py:function} straighten_mph_sets(group)
:canonical: soxspipe.commonutils.create_dispersion_map.straighten_mph_sets

```{autodoc2-docstring} soxspipe.commonutils.create_dispersion_map.straighten_mph_sets
```
````
