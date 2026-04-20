# {py:mod}`soxspipe.commonutils.toolkit`

```{py:module} soxspipe.commonutils.toolkit
```

```{autodoc2-docstring} soxspipe.commonutils.toolkit
:allowtitles:
```

## Module Contents

### Classes

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`MaxFilter <soxspipe.commonutils.toolkit.MaxFilter>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.MaxFilter
    :summary:
    ```
````

### Functions

````{list-table}
:class: autosummary longtable
:align: left

* - {py:obj}`add_recipe_logger <soxspipe.commonutils.toolkit.add_recipe_logger>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.add_recipe_logger
    :summary:
    ```
* - {py:obj}`add_snr_efficiency_qcs <soxspipe.commonutils.toolkit.add_snr_efficiency_qcs>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.add_snr_efficiency_qcs
    :summary:
    ```
* - {py:obj}`calculate_rolling_snr <soxspipe.commonutils.toolkit.calculate_rolling_snr>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.calculate_rolling_snr
    :summary:
    ```
* - {py:obj}`create_dispersion_solution_grid_lines_for_plot <soxspipe.commonutils.toolkit.create_dispersion_solution_grid_lines_for_plot>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.create_dispersion_solution_grid_lines_for_plot
    :summary:
    ```
* - {py:obj}`cut_image_slice <soxspipe.commonutils.toolkit.cut_image_slice>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.cut_image_slice
    :summary:
    ```
* - {py:obj}`extinction_correction_factor <soxspipe.commonutils.toolkit.extinction_correction_factor>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.extinction_correction_factor
    :summary:
    ```
* - {py:obj}`frame_to_32 <soxspipe.commonutils.toolkit.frame_to_32>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.frame_to_32
    :summary:
    ```
* - {py:obj}`generic_quality_checks <soxspipe.commonutils.toolkit.generic_quality_checks>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.generic_quality_checks
    :summary:
    ```
* - {py:obj}`get_calibration_lamp <soxspipe.commonutils.toolkit.get_calibration_lamp>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.get_calibration_lamp
    :summary:
    ```
* - {py:obj}`get_calibrations_path <soxspipe.commonutils.toolkit.get_calibrations_path>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.get_calibrations_path
    :summary:
    ```
* - {py:obj}`plot_merged_spectrum_qc <soxspipe.commonutils.toolkit.plot_merged_spectrum_qc>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.plot_merged_spectrum_qc
    :summary:
    ```
* - {py:obj}`predict_product_path <soxspipe.commonutils.toolkit.predict_product_path>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.predict_product_path
    :summary:
    ```
* - {py:obj}`qc_settings_plot_tables <soxspipe.commonutils.toolkit.qc_settings_plot_tables>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.qc_settings_plot_tables
    :summary:
    ```
* - {py:obj}`quicklook_image <soxspipe.commonutils.toolkit.quicklook_image>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.quicklook_image
    :summary:
    ```
* - {py:obj}`read_spectral_format <soxspipe.commonutils.toolkit.read_spectral_format>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.read_spectral_format
    :summary:
    ```
* - {py:obj}`spectroscopic_image_quality_checks <soxspipe.commonutils.toolkit.spectroscopic_image_quality_checks>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.spectroscopic_image_quality_checks
    :summary:
    ```
* - {py:obj}`twoD_disp_map_image_to_dataframe <soxspipe.commonutils.toolkit.twoD_disp_map_image_to_dataframe>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.twoD_disp_map_image_to_dataframe
    :summary:
    ```
* - {py:obj}`unpack_order_table <soxspipe.commonutils.toolkit.unpack_order_table>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.unpack_order_table
    :summary:
    ```
* - {py:obj}`utility_setup <soxspipe.commonutils.toolkit.utility_setup>`
  - ```{autodoc2-docstring} soxspipe.commonutils.toolkit.utility_setup
    :summary:
    ```
````

### API

`````{py:class} MaxFilter(max_level)
:canonical: soxspipe.commonutils.toolkit.MaxFilter

```{autodoc2-docstring} soxspipe.commonutils.toolkit.MaxFilter
```

```{rubric} Initialization
```

```{autodoc2-docstring} soxspipe.commonutils.toolkit.MaxFilter.__init__
```

````{py:method} filter(record)
:canonical: soxspipe.commonutils.toolkit.MaxFilter.filter

```{autodoc2-docstring} soxspipe.commonutils.toolkit.MaxFilter.filter
```

````

`````

````{py:function} add_recipe_logger(log, productPath)
:canonical: soxspipe.commonutils.toolkit.add_recipe_logger

```{autodoc2-docstring} soxspipe.commonutils.toolkit.add_recipe_logger
```
````

````{py:function} add_snr_efficiency_qcs(log, spectrumDF, qcTable, orderJoins, recipeName, dateObs)
:canonical: soxspipe.commonutils.toolkit.add_snr_efficiency_qcs

```{autodoc2-docstring} soxspipe.commonutils.toolkit.add_snr_efficiency_qcs
```
````

````{py:function} calculate_rolling_snr(dataframe, flux_column, window_size)
:canonical: soxspipe.commonutils.toolkit.calculate_rolling_snr

```{autodoc2-docstring} soxspipe.commonutils.toolkit.calculate_rolling_snr
```
````

````{py:function} create_dispersion_solution_grid_lines_for_plot(log, dispMap, dispMapImage, associatedFrame, kw, skylines=False, slitPositions=False, slit_length=11)
:canonical: soxspipe.commonutils.toolkit.create_dispersion_solution_grid_lines_for_plot

```{autodoc2-docstring} soxspipe.commonutils.toolkit.create_dispersion_solution_grid_lines_for_plot
```
````

````{py:function} cut_image_slice(log, frame, width, length, x, y, sliceAxis='x', median=False, debug=False)
:canonical: soxspipe.commonutils.toolkit.cut_image_slice

```{autodoc2-docstring} soxspipe.commonutils.toolkit.cut_image_slice
```
````

````{py:function} extinction_correction_factor(wave, extinctionTablePath, airmass)
:canonical: soxspipe.commonutils.toolkit.extinction_correction_factor

```{autodoc2-docstring} soxspipe.commonutils.toolkit.extinction_correction_factor
```
````

````{py:function} frame_to_32(frame)
:canonical: soxspipe.commonutils.toolkit.frame_to_32

```{autodoc2-docstring} soxspipe.commonutils.toolkit.frame_to_32
```
````

````{py:function} generic_quality_checks(log, frame, settings, recipeName, qcTable)
:canonical: soxspipe.commonutils.toolkit.generic_quality_checks

```{autodoc2-docstring} soxspipe.commonutils.toolkit.generic_quality_checks
```
````

````{py:function} get_calibration_lamp(log, frame, kw)
:canonical: soxspipe.commonutils.toolkit.get_calibration_lamp

```{autodoc2-docstring} soxspipe.commonutils.toolkit.get_calibration_lamp
```
````

````{py:function} get_calibrations_path(log, settings)
:canonical: soxspipe.commonutils.toolkit.get_calibrations_path

```{autodoc2-docstring} soxspipe.commonutils.toolkit.get_calibrations_path
```
````

````{py:function} plot_merged_spectrum_qc(merged_orders, products, log, qcDir, filenameTemplate, noddingSequence, dateObs, arm, recipeName, orderJoins=False, debug=False, fluxCalibrated=False)
:canonical: soxspipe.commonutils.toolkit.plot_merged_spectrum_qc

```{autodoc2-docstring} soxspipe.commonutils.toolkit.plot_merged_spectrum_qc
```
````

````{py:function} predict_product_path(sofName, recipeName=False)
:canonical: soxspipe.commonutils.toolkit.predict_product_path

```{autodoc2-docstring} soxspipe.commonutils.toolkit.predict_product_path
```
````

````{py:function} qc_settings_plot_tables(log, qc, qcAx, settings, settingsAx)
:canonical: soxspipe.commonutils.toolkit.qc_settings_plot_tables

```{autodoc2-docstring} soxspipe.commonutils.toolkit.qc_settings_plot_tables
```
````

````{py:function} quicklook_image(log, CCDObject, show=True, ext='data', stdWindow=3, title=False, surfacePlot=False, dispMap=False, dispMapImage=False, inst=False, settings=False, skylines=False, saveToPath=False)
:canonical: soxspipe.commonutils.toolkit.quicklook_image

```{autodoc2-docstring} soxspipe.commonutils.toolkit.quicklook_image
```
````

````{py:function} read_spectral_format(log, settings, arm, dispersionMap=False, extended=True, binx=1, biny=1)
:canonical: soxspipe.commonutils.toolkit.read_spectral_format

```{autodoc2-docstring} soxspipe.commonutils.toolkit.read_spectral_format
```
````

````{py:function} spectroscopic_image_quality_checks(log, frame, orderTablePath, settings, recipeName, qcTable)
:canonical: soxspipe.commonutils.toolkit.spectroscopic_image_quality_checks

```{autodoc2-docstring} soxspipe.commonutils.toolkit.spectroscopic_image_quality_checks
```
````

````{py:function} twoD_disp_map_image_to_dataframe(log, slit_length, twoDMapPath, kw=False, associatedFrame=False, removeMaskedPixels=False, dispAxis='y')
:canonical: soxspipe.commonutils.toolkit.twoD_disp_map_image_to_dataframe

```{autodoc2-docstring} soxspipe.commonutils.toolkit.twoD_disp_map_image_to_dataframe
```
````

````{py:function} unpack_order_table(log, orderTablePath, extend=0.0, pixelDelta=1, binx=1, biny=1, prebinned=False, order=False, limitToDetectorFormat=False)
:canonical: soxspipe.commonutils.toolkit.unpack_order_table

```{autodoc2-docstring} soxspipe.commonutils.toolkit.unpack_order_table
```
````

````{py:function} utility_setup(log, settings, recipeName, startNightDate)
:canonical: soxspipe.commonutils.toolkit.utility_setup

```{autodoc2-docstring} soxspipe.commonutils.toolkit.utility_setup
```
````
