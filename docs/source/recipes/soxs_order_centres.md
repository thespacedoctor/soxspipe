# `soxs_order_centres`

The purpose of the [`soxs_order_centres`](../_api/soxspipe.recipes.soxs_order_centres.html) recipe is to find and fit the order centres with low-level polynomials.

### Input

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS Image | Flat lamp through a single-pinhole mask | `SOXS_slt_cal_VISLampFlatPinhole`, `SOXS_slt_cal_NIRLampFlatPinhole` |
| FITS Image | Master Dark Frame (VIS only) | - |
| FITS Image | Master Bias Frame (VIS only) | - |
| FITS Image | Dark frame (Lamp-Off) of equal exposure length as single-pinhole frame (Lamp-On) (NIR only) | `SOXS_slt_cal_NIRLampFlatPinhole` |
| CSV File | First guess dispersion solution | - |

### Parameters

| Parameter                | Description                                   | Type  | Entry Point   | Related Util                                   |
| ------------------------ | --------------------------------------------- | ----- | ------------- | ---------------------------------------------- |
| order-sample-count  | number of times along the order in the dispersion direction to measure the order-centre trace  |  int | settings file |  [detect_continuum utility](../utils/detect_continuum.md) |
| peak-sigma-limit  |  minimum value a peak must be above the median value of pixel to be considered for order-trace fitting  | int | settings file  |  [detect_continuum utility](../utils/detect_continuum.md) |
| poly-deg  |  the order of polynomial used to fit the order-centre traces  | int | settings file  |  [detect_continuum utility](../utils/detect_continuum.md) |
| poly-fitting-residual-clipping-sigma  | sigma distance limit, where distance is the difference between the detected and polynomial fitted positions of the order-trace, outside of which to remove lines from the fit   | float   | settings file |  [detect_continuum utility](../utils/detect_continuum.md) | 
| clipping-iteration-limit  |  number of times to perform sigma-clipping (if number of clipped lines remains above zero) |  int   | settings file | [detect_continuum utility](../utils/detect_continuum.md) |

### Method

Once the single-pinhole flat-lamp frame has had the bias, dark and background subtracted it is passed to the [detect_continuum utility](../utils/detect_continuum.md) to fit the order centres.

![](soxs_order_centres.png)

### Output
 
| Data Type | Content |
|:----|:----|
| CSV File | [order table](../files/order_table.md) containing coefficients to the polynomial fits describing the order centre locations |

### QC Metrics

Plots similar to the one below are generated after each execution of [`soxs_order_centres`](../_api/soxspipe.recipes.soxs_order_centres.html).

[![](https://live.staticflickr.com/65535/50345130012_4e869a6a7f_b.png)](https://live.staticflickr.com/65535/50345130012_4e869a6a7f_o.png)

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |

### Recipe API

```eval_rst
.. autoclass:: soxspipe.recipes.soxs_order_centres
    :members:
```
