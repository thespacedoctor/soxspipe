

:::{table} The `soxs_mdark` recipe parameters.
:name: soxs_mdark_parameters

| Parameter                   | Description                                                  | Type  | Entry Point   | Related Util                                   |
| --------------------------- | ------------------------------------------------------------ | ----- | ------------- | ---------------------------------------------- |
| frame-clipping-sigma        | number of σ from the median *frame* flux beyond which pixel is added to the bad-pixel mask | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| frame-clipping-iterations   | number of sigma-clipping iterations to perform when adding pixels to the bad-pixel mask | int   | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-sigma      | number of σ deviations from the median *pixel* flux beyond which pixel is excluded from stack | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-iterations | number of σ-clipping iterations to perform before stacking   | int   | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |

:::


