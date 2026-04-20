:::{table} The `soxs_order_centres` recipe parameters.
:name: soxs_order_centres_parameters
:widths: 25, 30, 10, 15, 20
:class: longtable

| Parameter                            | Description                                                                                                  | Type  | Entry Point                   | Related Util                                             |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------ | ----- | ----------------------------- | -------------------------------------------------------- |
| `order_sample_count`                   | number of cross-order slices per order                                                                         | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| `slice_length`                         | length of each slice (pixels)                                                                                  | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| `slice_width`                          | width of each slice (pixels)                                                                                   | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| `peak_sigma_limit`                     | height gaussian peak must be above median flux to be "detected" by code (std via median absolute deviation). | float | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| `disp_axis_deg`                        | degree of y-component of global polynomial fit to order centres                                                | int   | settings file or command-line | [detect_continuum utility](../utils/detect_continuum.md) |
| `order_deg`                            | degree of echelle order number component of global polynomial fit to order centres                             | int   | settings file or command-line | [detect_continuum utility](../utils/detect_continuum.md) |
| `poly_fitting_residual_clipping_sigma` | sigma clipping limit when fitting global polynomial to order centres                                           | float | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| `poly_clipping_iteration_limit`        | maximum number of clipping iterations when fitting global polynomial to order centres                          | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |


:::

