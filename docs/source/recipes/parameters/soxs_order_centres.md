:::{table} The `soxs_order_centres` recipe parameters.
:name: soxs_order_centres_parameters

| Parameter                            | Description                                                                                                  | Type  | Entry Point                   | Related Util                                             |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------ | ----- | ----------------------------- | -------------------------------------------------------- |
| order-sample-count                   | number of cross-order slices per order                                                                       | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| slice-length                         | length of each slice (pixels)                                                                                | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| slice-width                          | width of each slice (pixels)                                                                                 | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| peak-sigma-limit                     | height gaussian peak must be above median flux to be "detected" by code (std via median absolute deviation). | float | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| disp-axis-deg                        | degree of y-component of global polynomial fit to order centres                                              | int   | settings file or command-line | [detect_continuum utility](../utils/detect_continuum.md) |
| order-deg                            | degree of echelle order number component of global polynomial fit to order centres                           | int   | settings file or command-line | [detect_continuum utility](../utils/detect_continuum.md) |
| poly-fitting-residual-clipping-sigma | sigma clipping limit when fitting global polynomial to order centres                                         | float | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |
| poly-clipping-iteration-limit        | maximum number of clipping iterations when fitting global polynomial to order centres                        | int   | settings file                 | [detect_continuum utility](../utils/detect_continuum.md) |


:::

