:::{table} The `soxs_mflat` recipe parameters.
:name: soxs_mflat_parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| frame-clipping-sigma     | number of σ from the median *frame* flux beyond which pixel is added to the bad-pixel mask    | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| clipping-iteration-count | number of sigma-clipping iterations to perform when added pixels to the bad-pixel mask | int   | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-sigma | number of σ deviations from the median *pixel* flux beyond which pixel is excluded from stack | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-iterations | number of σ-clipping iterations to perform before stacking | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| centre-order-window | the width of the slice to cut along the centre of each order when determining mean exposure level | int | settings file | -  |
| slice-length-for-edge-detection | length of image slice to take across orders when detecting edges | int | settings file |[`detect_order_edges`](../utils/detect_order_edges.md) |
| slice-width-for-edge-detection | width of image slice to take across orders when detecting edges | int | settings file |[`detect_order_edges`](../utils/detect_order_edges.md) |
| min-percentage-threshold-for-edge-detection | minimum value flux can drop to as percentage of central flux and be counted as an order edge | int |settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| max-percentage-threshold-for-edge-detection | maximum value flux can claim to as percentage of central flux and be counted as an order edge | int | settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| disp-axis-deg | degree of dispersion axis component of polynomial fit to order edges | int | settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| order-deg | degree of order component of polynomial fit to order edges | int | settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| poly-fitting-residual-clipping-sigma | number of σ deviations from the median fit residual beyond which individual data points are removed when iterating towards a fit of order edges | int | settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| poly-clipping-iteration-limit | number of sigma-clipping iterations to perform before settings on a polynomial fit for the order edges | int | settings file | [`detect_order_edges`](../utils/detect_order_edges.md) |
| low-sensitivity-clipping-sigma | number of σ deviations below the median flux of a master-flat frame beyond which a pixel is added to the bad-pixel mask | int | settings file | - |



:::
