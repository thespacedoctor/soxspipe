# {py:mod}`soxspipe.commonutils.response_function`

```{py:module} soxspipe.commonutils.response_function
```

```{autodoc2-docstring} soxspipe.commonutils.response_function
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`response_function <soxspipe.commonutils.response_function.response_function>`
  - ```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function
    :summary:
    ```
````

### API

`````{py:class} response_function(log, stdExtractionPath, recipeName, sofName, settings=False, qcTable=False, productsTable=False, startNightDate='', stdNotFlatExtractionPath='', orderJoins=None)
:canonical: soxspipe.commonutils.response_function.response_function

Bases: {py:obj}`object`

```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function.__init__
```

````{py:method} get()
:canonical: soxspipe.commonutils.response_function.response_function.get

```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function.get
```

````

````{py:method} plot_response_curve(stdExtWave, stdExtWaveNotFlat, stdExtFlux, binCentreWave, binCentreWaveOriginal, binIntegratedFlux, absToExtFluxRatio, responseFuncCoeffs, stdEfficiencyEstimate)
:canonical: soxspipe.commonutils.response_function.response_function.plot_response_curve

```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function.plot_response_curve
```

````

````{py:method} write_response_function_to_file(responseFuncCoeffs, polyOrder)
:canonical: soxspipe.commonutils.response_function.response_function.write_response_function_to_file

```{autodoc2-docstring} soxspipe.commonutils.response_function.response_function.write_response_function_to_file
```

````

`````
