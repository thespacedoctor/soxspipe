:::{table} The `soxs_mflat` recipe parameters.
:name: soxs_mflat_parameters

| Parameter                                   | Description                                                                                          | Type  | Entry Point                   | Related Util                                                       |
| ------------------------------------------- | ---------------------------------------------------------------------------------------------------- | ----- | ----------------------------- | ------------------------------------------------------------------ |
| subtract_background                         | fit and subtract the intra-order background light                                                    | bool  | settings file                 | [`subtract_background`](../utils/subtract_background.md) |
| stacked-clipping-sigma                      | the sigma clipping limit used when stacking frames into a composite frame                            | float | settings file                 | [`clip_and_stack`](../utils/clip_and_stack.md)                     |
| stacked-clipping-iterations                 | the maximum sigma-clipping iterations used when stacking frames into a composite frame               | int   | settings file                 | [`clip_and_stack`](../utils/clip_and_stack.md)                     |
| centre-order-window                         | size of the window (in pixels) used to measure flux in the central band to determine median exposure | int   | settings file                 | -                                                                  |
| slice-length-for-edge-detection             | length of the cross_dispersion slices used to determine order edges                                  | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| slice-width-for-edge-detection              | width of the cross_dispersion slices used to determine order edges                                   | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| min-percentage-threshold-for-edge-detection | minimum value flux can drop to as percentage of central flux and be counted as an order edge         | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| max-percentage-threshold-for-edge-detection | maximum value flux can climb to as percentage of central flux and be counted as an order edge        | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| disp-axis-deg                               | degree of y-component of global polynomial fit to order edges                                        | int   | settings file or command-line | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| order-deg                                   | degree of echelle order number component of global polynomial fit to order edges                     | int   | settings file or command-line | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| poly-fitting-residual-clipping-sigma        | sigma clipping limit when fitting global polynomial to order edges                                   | float | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| poly-clipping-iteration-limit               | maximum number of clipping iterations when fitting global polynomial to order edges                  | int   | settings file                 | [`detect_order_edges`](../utils/detect_order_edges.md)             |
| low-sensitivity-clipping-sigma              | pixels with a flux less than this many sigma (mad) below the median flux level added to the bp map   | float | settings file                 | -                                                                  |
| scale-d2-to-qth                             | scale d2 to qth lamp flats when stitching                                                            | bool  | settings file                 | -                                                                  |



:::




