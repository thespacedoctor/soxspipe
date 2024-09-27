# {py:mod}`soxspipe.commonutils.polynomials`

```{py:module} soxspipe.commonutils.polynomials
```

```{autodoc2-docstring} soxspipe.commonutils.polynomials
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`chebyshev_order_wavelength_polynomials <soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials>`
  - ```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials
    :summary:
    ```
* - {py:obj}`chebyshev_order_xy_polynomials <soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials>`
  - ```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials
    :summary:
    ```
* - {py:obj}`chebyshev_xy_polynomial <soxspipe.commonutils.polynomials.chebyshev_xy_polynomial>`
  - ```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_xy_polynomial
    :summary:
    ```
````

### API

`````{py:class} chebyshev_order_wavelength_polynomials(log, orderDeg, wavelengthDeg, slitDeg, exponentsIncluded=False, axis=False)
:canonical: soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials.__init__
```

````{py:method} poly(orderPixelTable, *coeff)
:canonical: soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials.poly

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_wavelength_polynomials.poly
```

````

`````

`````{py:class} chebyshev_order_xy_polynomials(log, orderDeg, axisBDeg, axisB='y', axisBCol=False, orderCol=False, exponentsIncluded=False)
:canonical: soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials.__init__
```

````{py:method} poly(orderPixelTable, *coeff)
:canonical: soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials.poly

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_order_xy_polynomials.poly
```

````

`````

`````{py:class} chebyshev_xy_polynomial(log, y_deg, yCol=False, exponentsIncluded=False)
:canonical: soxspipe.commonutils.polynomials.chebyshev_xy_polynomial

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_xy_polynomial
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_xy_polynomial.__init__
```

````{py:method} poly(orderPixelTable, *coeff)
:canonical: soxspipe.commonutils.polynomials.chebyshev_xy_polynomial.poly

```{autodoc2-docstring} soxspipe.commonutils.polynomials.chebyshev_xy_polynomial.poly
```

````

`````
