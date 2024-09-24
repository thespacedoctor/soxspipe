:::{table} The `soxs_nod` recipe parameters.
:name: soxs_nod_parameters

| Parameter                   | Description                                                  | Type  | Entry Point   | Related Util                                   |
| --------------------------- | ------------------------------------------------------------ | ----- | ------------- | ---------------------------------------------- |
use_flat | If True, flat field correction if applied. If false, it will skipped | Boolean | Settings File | soxs_nod |
subtract_background | If True, background light subtraction is applied. If false, it will skipped | Boolean | Settings File | soxs_nod |
save_single_frame_extractions | If True, each A and B stacked frame for each offset iaresaved as 2D image | Boolean | Settings File | soxs_nod |
stacked-clipping-sigma | number of $\sigma$ deviations from the median *pixel* flux beyond which pixel is excluded from stack | Boolean | Setting File | [`clip_and_stack`](../utils/clip_and_stack.md)
stacked-clipping-iterations | number of $\sigma$-clipping iterations to perform before stacking | int | Setting File | [`clip_and_stack`](../utils/clip_and_stack.md)
horne-extraction-slit-length | the length of the 'slit' used to collect object flux (in pixels) | int | Setting File | [`horne_extraction`](../utils/horne_extraction.md)
horne-extraction-profile-global-clipping-sigma | the number of $\sigma$ deviations beyond the median pixel value for a pixel to be clipped (removing CRH and bad-pixels) during the fitting of the object profile| int | Setting File | [`horne_extraction`](../utils/horne_extraction.md)
horne-extraction-profile-clipping-sigma | the number of $\sigma$ deviations that residuals from the fitted, dispersion-direction profile need to be beyond for a pixel to be clipped during the estimation of the object profile | int |  Setting File |[`horne_extraction`](../utils/horne_extraction.md)
horne-extraction-profile-clipping-iteration-count | number of sigma-clipping iterations to perform during profile estimation of the object | int | Setting File | [`horne_extraction`](../utils/horne_extraction.md)
order-sample-count | Number of cross order slices per order | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
slice-length | Length of each slices used for the continuum detection | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
slice-width | Width of each slices used for the continuum detection | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
peak-sigma-limit | Height of the gaussian peak (in std dev) for continuum detection | float | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
order-deg | Degree of the order-component of the  global polynomial fit of the detected positions of the continuum | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
disp-axis-deg | Degree of the Y-component of the  global polynomial fit of the detected positions of the continuum | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
poly-fitting-residual-clipping-sigma | Clipping limit (MEDIAN AND MAD) when fitting global polynomial to the detected positions of the continuum | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)
poly-clipping-iteration-limit | Iteration limit on the polynomial fitting of the continuum | int | Setting File | [`detect_continuum`](../utils/detect_continuum.md)

:::
