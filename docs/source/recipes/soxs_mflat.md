# `soxs_mflat`

The purpose of the [`soxs_mflat`](#soxspipe.recipes.soxs_mflat) recipe is to create a single normalised [master-flat frame](../files/master_flat.md) used to correct for non-uniformity in response to light across the detector plain.
 
Sources of this non-uniformity include varying pixel sensitivities, obstructions in the optical path (e.g. dust or pollen grains), vignetting at the edges of the detector. A flat-frame is ideally an image taken where the illumination is uniform across the light collecting pixels of the detector. This evenly exposed image can be used to identify irregularities in the response in the detector.

## Input

<!-- FIND OBs HERE : https://docs.google.com/spreadsheets/d/1-3VXkIWcydvpawwVl_C3pNTU3HgnElJaYFAKow65Fl8/edit#gid=0 -->

| Data Type | Content | Related OB |
|:----|:----|:---|
| FITS images | raw flats frames (exposures with identical exposure time and detectors readout parameters). UV-VIS requires separate sets D-Lamp and QTH-Lamp flats. | `SOXS_slt_cal_NIRLampFlat`, `SOXS_slt_cal_NIRLampFlatAtt`, `SOXS_slt_cal_VISLampFlat`, `SOXS_slt_cal_VISLampFlatAtt` |
| FITS Image | Master Bias Frame (UV-VIS only) | - |
| FITS Table | [order table](../files/order_table.md) containing coefficients to the polynomial fits describing the order centre locations. UV-VIS requires separate tables for D-Lamp and QTH-Lamp. | |

## Parameters

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


## Method

![](soxs_mflat.png)

The individual flat field frames need to have bias and dark signatures removed before they are combined. This is achieved with the [`detrend`](../utils/detrend.md) utility. Here is an example of one such calibrated flat frame:

[![](https://live.staticflickr.com/65535/51237891523_cc24b22ba7_b.jpg)](https://live.staticflickr.com/65535/51237891523_cc24b22ba7_b.jpg)

#### Normalising Exposure Levels in Individual Flat Frames

Once calibrated, exposure-levels in the individual flat frames need to be normalised as *total* illumination will vary from frame-to-frame. The individual frame exposure levels are calculated in two stages.

In the first stage the mean inner-order pixel-value across the frame is used as a *first approximation* of an individual frame's exposure level. To calculate this mean value, the order locations are used to identify a curved slice N-pixels wide centred on each of the order-centres and bad pixels are masked (see image below). The collected inner order pixel values are then sigma-clipped to excluded out-lying values and a mean value calculated.

[![](https://live.staticflickr.com/65535/51237856613_da50864d52_b.jpg)](https://live.staticflickr.com/65535/51237856613_da50864d52_b.jpg)

Individual frames are then divided through by their mean inner-order pixel value in this first attempt to normalise the exposure-levels of the frames. 

[![](https://live.staticflickr.com/65535/51238466119_6844cce08b_b.jpg)](https://live.staticflickr.com/65535/51238466119_6844cce08b_b.jpg)

The normalised flat-frames are then combined using the [`clip_and_stack`](../utils/clip_and_stack/.md) utility into a first-pass master-flat frame:

[![](https://live.staticflickr.com/65535/51237081792_6090fe56df_b.jpg)](https://live.staticflickr.com/65535/51237081792_6090fe56df_b.jpg)

The second stage then is to divide each original dark and bias subtracted flat frame by this first-pass master flat (see example below). This removes the *typical* cross-plane illumination and so now the mean inner-order pixel-value across the frame will give a much better estimate of each frame's intrinsic exposure level.

[![](https://live.staticflickr.com/65535/51237134772_bf081c62cc_b.jpg)](https://live.staticflickr.com/65535/51237134772_bf081c62cc_b.jpg)

The mean inner-order pixel-value is calculated again on this frame and the original dark and bias subtracted flat is re-normalised by divided through by this accurate measurement of its intrinsic exposure level.

[![](https://live.staticflickr.com/65535/51237233602_c3c96b3503_b.jpg)](https://live.staticflickr.com/65535/51237233602_c3c96b3503_b.jpg)

#### Building a Final Master-Flat Frame

These re-normalised flats are then combined for a second time into a master-flat frame.

[![](https://live.staticflickr.com/65535/51237949461_0bc4763a06_b.jpg)](https://live.staticflickr.com/65535/51237949461_0bc4763a06_b.jpg)

Finally order edges are located with the [`detect_order_edges `](../utils/detect_order_edges.md) utility and the inter-order area pixel value are set to 1. 

Low-sensitivity pixels are flagged and added to the bad-pixel map and a final master-flat frame written to file.

[![](https://live.staticflickr.com/65535/51239008475_b7c0aa33c7_b.jpg)](https://live.staticflickr.com/65535/51239008475_b7c0aa33c7_b.jpg)

#### UV Master Flat Frame Stitching

![](stitch_uv_mflats.png)

As the UV-VIS uses a combination of D-Lamp and QTH-Lamp flat sets, a further step is required to stitch the best orders from each of these master-flats together into a dual lamp master-flat.

For both the D-Lamp and QTH-Lamp master-flat frames, we have for each order the number of pixel positions that contributed to the final order-edge fit. We use these numbers to decide which orders to slice and stitch from the D-Lamp to the QTH-Lamp master flat frame. 

With a crossover order now selected, the median flux from a square window at the centre of this order in both D- and QTH frames is measured. Using the ratio of these fluxes the D-Lamp frame is scaled to the QTH-Lamp frame.

From the upper order-edge polynomial for the D-Lamp we define a curved, intra-order line 5 pixels above the upper edge of the crossover order selected previously. This line is used to slice and stitch the D-Lamp and QTH-Lamp orders together. This process is done on the flux images, error maps and bad-pixel maps. Typically the bluest orders from the D-Lamp will be selected with the remaining orders coming from the QTH-Lamp.

Finally, the combined normalised frames for both the D and QTH-Lamps are stacked to obtain a good level of flux in each order. This stacked frame is used to re-detect the order edges (the resulting order table is used going forward).


## Output
 
| Data Type | Content |
|:----|:----|
| master flat frame |  frame used correct for non-uniformity in response to light across the detector plain (including blaze) |

## QC Metrics

| Metric  | Description |
| :------------ | :----------- |
| TBC     | ...  |

## Recipe API

:::{autodoc2-object} soxspipe.recipes.soxs_mflat.soxs_mflat
:::
