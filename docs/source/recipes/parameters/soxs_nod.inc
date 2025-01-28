:::{table} The `soxs_nod` recipe parameters.
:name: soxs_nod_parameters
:widths: 25, 30, 10, 15, 20

| Parameter                                         | Description                                                                                                 | Type  | Entry Point   | Related Util                                       |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------- | ----- | ------------- | -------------------------------------------------- |
| use_flat                                          | divide image by master flat frame                                                                           | bool  | settings file | -                                                  |
| subtract_background                               | fit and subtract the intra-order background light                                                           | bool  | settings file | -                                                  |
| save_single_frame_extractions                     | If True, each A and B stacked frame for each offset is saved as 2D image                                    | bool  | settings file | -                                                  |
| stacked-clipping-sigma                            | the sigma clipping limit used when stacking frames into a composite frame                                   | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md)     |
| stacked-clipping-iterations                       | the maximum sigma-clipping iterations used when stacking frames into a composite frame                      | int   | settings file | [`clip_and_stack`](../utils/clip_and_stack.md)     |
| horne-extraction-slit-length                      | the length of the 'slit' used to collect object flux (in pixels)                                            | int   | settings file | [`horne_extraction`](../utils/horne_extraction.md) |
| horne-extraction-profile-global-clipping-sigma    | sigma clipping limit when fitting the object profile (global over the order)                                | float | settings file | [`horne_extraction`](../utils/horne_extraction.md) |
| horne-extraction-profile-clipping-sigma           | sigma clipping limit when fitting the dispersion-direction profiles of the object                           | float | settings file | [`horne_extraction`](../utils/horne_extraction.md) |
| horne-extraction-profile-clipping-iteration-count | maximum number of clipping iterations when fitting dispersion-direction profiles                            | int   | settings file | [`horne_extraction`](../utils/horne_extraction.md) |
| order-sample-count                                | number of cross-order slices per order                                                                      | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| slice-length                                      | length of each slice (pixels)                                                                               | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| slice-width                                       | width of each slice (pixels)                                                                                | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| peak-sigma-limit                                  | height gaussian peak must be above median flux to be "detected" by code (std via median absolute deviation) | float | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| order-deg                                         | degree of order-component of global polynomial fit to object trace                                          | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| disp-axis-deg                                     | degree of y-component of global polynomial fit to object trace                                              | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| poly-fitting-residual-clipping-sigma              | clipping limit (median and mad) when fitting global polynomial to object trace                              | float | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| poly-clipping-iteration-limit                     | maximum number of clipping iterations when fitting global polynomial to object trace                        | int   | settings file | [`detect_continuum`](../utils/detect_continuum.md) |
| background-subtraction: bspline-deg         | degree of bsplines used to fit the inter-order background (if `subtract_background` == True)         | int  | settings file                  | [`subtract_background`](../utils/subtract_background.md) |
| background-subtraction: gaussian-blur-sigma | Standard deviation of Gaussian kernel used to smooth background image (if `subtract_background` == True)  | int  | settings file                 | [`subtract_background`](../utils/subtract_background.md) |

:::
