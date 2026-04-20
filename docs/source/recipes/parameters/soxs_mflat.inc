:::{table} The `soxs_mflat` recipe parameters.
:name: soxs_mflat_parameters
:widths: 25, 30, 10, 15, 20
:class: longtable

| Parameter                                   | Description                                                                                          | Type  | Entry Point                   | Related Util                                                       |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ----- | ----------------------------- | ------------------------------------------------------------------ |
| `subtract_background`                         | fit and subtract the intra-order background light                                                    | bool  | settings file                 | [`subtract_background`](../utils/subtract_background.md) |
| `stacked_clipping_sigma`                      | the sigma clipping limit used when stacking frames into a composite frame                            | float | settings file                 | [`clip_and_stack`](../utils/clip_and_stack.md)                     |
| `stacked_clipping_iterations`                 | the maximum sigma-clipping iterations used when stacking frames into a composite frame               | int   | settings file                 | [`clip_and_stack`](../utils/clip_and_stack.md)                     |
| `centre_order_window`                         | size of the window (in pixels) used to measure flux in the central band to determine median exposure | int   | settings file                 | -                                                                  |
| `slice_length_for_edge_detection`             | length of the cross_dispersion slices used to determine order edges                                  | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `slice_width_for_edge_detection`              | width of the cross_dispersion slices used to determine order edges                                   | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `min_percentage_threshold_for_edge_detection` | minimum value flux can drop to as percentage of central flux and be counted as an order edge         | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `max_percentage_threshold_for_edge_detection` | maximum value flux can climb to as percentage of central flux and be counted as an order edge        | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `disp_axis_deg`                               | degree of y-component of global polynomial fit to order edges                                        | int   | settings file or command-line | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `order_deg`                                   | degree of echelle order number component of global polynomial fit to order edges                     | int   | settings file or command-line | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `poly_fitting_residual_clipping_sigma`        | sigma clipping limit when fitting global polynomial to order edges                                   | float | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `poly_clipping_iteration_limit`               | maximum number of clipping iterations when fitting global polynomial to order edges                  | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| `low_sensitivity_clipping_sigma`              | pixels with a flux less than this many sigma (mad) below the median flux level added to the bp map   | float | settings file                 | -                                                                  |
| `scale_d2_to_qth`                             | scale d2 to qth lamp flats when stitching                                                            | bool  | settings file                 | -                                                                  |
| `background_subtraction: bspline_deg`         | degree of bsplines used to fit the inter-order background (if `subtract_background` == True)         | int  | settings file                  | [`subtract_background`](../utils/subtract_background.md) |
| `background_subtraction: gaussian_blur_sigma` | Standard deviation of Gaussian kernel used to smooth background image (if `subtract_background` == True)  | int  | settings file                 | [`subtract_background`](../utils/subtract_background.md) |

:::




