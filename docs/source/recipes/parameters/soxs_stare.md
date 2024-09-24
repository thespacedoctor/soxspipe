:::{table} The `soxs_stare` recipe parameters.
:name: soxs_stare_parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| stacked-clipping-sigma | number of σ deviations from the median *pixel* flux beyond which pixel is excluded from stack | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-iterations | number of σ-clipping iterations to perform before stacking | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
|   |   |   |   |

:::
