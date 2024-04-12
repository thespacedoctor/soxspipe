# `soxs_mbias`

A zero-second exposure will contain only read-noise and \~half of pixels within this Gaussian distribution centred around zero count will always contain negative flux. To avoid negative counts an offset *bias* voltage is applied at the amplifier stage so that even when no photons are detected the A/D converter will always register a positive value. This bias-voltage offset must be accounted for in the data reduction process. 

The purpose of the [`soxs_mbias`](../_api/soxspipe.recipes.soxs_mbias.html) recipe is to provide a [master-bias frame](../files/master_bias.md) that can be subtracted from science/calibration frames to remove the contribution of pixel counts resulting from the bias-voltage.

### Input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS images | raw bias frames (UV-VIS/AC exposures with exptime = 0) | `SOXS_img_cal_Bias`, `SOXS_gen_cal_VISBias` |

### Parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| frame-clipping-sigma     | number of σ from the median *frame* flux beyond which pixel is added to the bad-pixel mask    | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| clipping-iteration-count | number of sigma-clipping iterations to perform when added pixels to the bad-pixel mask | int   | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-sigma | number of σ deviations from the median *pixel* flux beyond which pixel is excluded from stack | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
| stacked-clipping-iterations | number of σ-clipping iterations to perform before stacking | float | settings file | [`clip_and_stack`](../utils/clip_and_stack.md) |
[Input parameters for the soxs-mbias recipe][tab:soxs-mbias-parameters]

### Method

The purpose of the [`soxs_mbias`](../_api/soxspipe.recipes.soxs_mbias.html) recipe is to stack raw bias-frames together (using the [`clip_and_stack`](../utils/clip_and_stack.md) utility) into [master-bias frames](../files/master_bias.md) and in the process clipping rogue pixels from the individual raw frames and reducing the read-noise contribution.

![](soxs_mbias.png)

### Output

| Data Type | Content |
|:----|:----|
| FITS image | Master bias frame (frame containing typical bias-voltage applied to the detector) | 

### QC Metrics

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |


### Recipe API

```eval_rst
.. autoclass:: soxspipe.recipes.soxs_mbias
    :members:
```


