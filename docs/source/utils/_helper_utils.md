# Helper Utilities

## add_recipe_logger

Add a recipe-specific handler to the default logger that writes the recipe logs adjacent to the recipe project.

:::{autodoc2-object} soxspipe.commonutils.toolkit.add_recipe_logger
:::

## create_dispersion_solution_grid_lines_for_plot

Use a dispersion solution table and 2D image map to generate a dispersion solution grid to add to QC plots.

:::{autodoc2-object} soxspipe.commonutils.toolkit.create_dispersion_solution_grid_lines_for_plot
:::

## cut_image_slice

Cut and return an N-pixel wide and M-pixel long image slice, centred on a given coordinate from the input image frame.

:::{autodoc2-object} soxspipe.commonutils.toolkit.cut_image_slice
:::

## detector_lookup

Return a dictionary of detector characteristics and parameters.

:::{autodoc2-object} soxspipe.commonutils.detector_lookup.detector_lookup
:::

## dispersion_map_to_pixel_arrays

Use a dispersion solution to convert wavelength, slit-position and echelle order numbers to X, Y pixel positions.

:::{autodoc2-object} soxspipe.commonutils.dispersion_map_to_pixel_arrays.dispersion_map_to_pixel_arrays
:::

## filenamer

Given a FITS object (HDU list and header), use the SOXS file-naming scheme to return a filename to be used to save the FITS object to disk.

:::{autodoc2-object} soxspipe.commonutils.filenamer.filenamer
:::

## generic_quality_checks

Measure some very basic quality checks on a frame and return the QC table with results appended

:::{autodoc2-object} soxspipe.commonutils.toolkit.generic_quality_checks
:::

## get_calibration_lamp

Given a frame (CCDObject), determine which calibration lamp is used.

:::{autodoc2-object} soxspipe.commonutils.toolkit.get_calibration_lamp
:::

## get_calibrations_path

Return the root path to the static calibrations (ships alongside code).

:::{autodoc2-object} soxspipe.commonutils.toolkit.get_calibrations_path
:::

## keyword_lookup

Given a tag (internal to the pipeline), and an optional keyword index, return the FITS Header keyword for the selected instrument.

:::{autodoc2-object} soxspipe.commonutils.keyword_lookup.keyword_lookup
:::

## predict_product_path

Predict the path of the recipe product from a given SOF file name.

:::{autodoc2-object} soxspipe.commonutils.toolkit.predict_product_path
:::

## qc_settings_plot_tables

Generate the QC and settings tables at the bottom of the QC plots.

:::{autodoc2-object} soxspipe.commonutils.toolkit.qc_settings_plot_tables
:::

## quicklook_image

Generate a quick look image of a CCDObject - useful for development/debugging.

:::{autodoc2-object} soxspipe.commonutils.toolkit.quicklook_image
:::

## read_spectral_format

Return a spectral format table for your selected instrument containing key parameters about the detector.

:::{autodoc2-object} soxspipe.commonutils.toolkit.read_spectral_format
:::

## uncompress

Uncompress ESO `fits.Z` frames before processing them with the data-organiser.

:::{autodoc2-object} soxspipe.commonutils.uncompress.uncompress
:::

## spectroscopic_image_quality_checks

Perform some generic image quality checks and add to the QC output of the recipe.

:::{autodoc2-object} soxspipe.commonutils.toolkit.spectroscopic_image_quality_checks
:::

## twoD_disp_map_image_to_dataframe

Convert the 2D dispersion image map to a pandas dataframe.

:::{autodoc2-object} soxspipe.commonutils.toolkit.twoD_disp_map_image_to_dataframe
:::

## unpack_order_table

Unpack an order location table and return an `orderPolyTable` dataframe containing the polynomial coefficients for the order centres and edges, an `orderPixelTable` dataframe containing the pixel-coordinates for each order centre and edges, and finally, an `orderMetaTable` dataframe giving metadata about the frame binning and format

:::{autodoc2-object} soxspipe.commonutils.toolkit.unpack_order_table
:::
