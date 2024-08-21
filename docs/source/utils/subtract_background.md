# `subtract_background` 

The purpose of the [`subtract_background`](#soxspipe.commonutils.subtract_background) utility is to model the topology of the scattered background light within an image and then remove it.

Here's the workflow for subtracting a frame's background:

![](subtract_background.png)

Here's an example frame requiring the background scattered light to be fitted and removed:

[![](https://live.staticflickr.com/65535/51236661920_a034db6805_b.jpg)](https://live.staticflickr.com/65535/51236661920_a034db6805_b.jpg)

Having unpacked the order location table, a mask is created containing pixels lying within the order-locations. The mask is extended in either direction along the y-axis by a given fraction (recipe parameter) of its original length to make sure all inner-order flux is masked.

[![](https://live.staticflickr.com/65535/51235592981_e42aafbe63_b.jpg)](https://live.staticflickr.com/65535/51235592981_e42aafbe63_b.jpg)

The bad-pixel mask is merged with this inner-order mask. For each row in the masked frame, a bspline is fitted to the unmasked pixel fluxes to model the shape of the scattered background light along the row.

[![](https://live.staticflickr.com/65535/51234447259_1e10610df3_b.jpg)](https://live.staticflickr.com/65535/51234447259_1e10610df3_b.jpg)

For each row, flux-values are generated for all pixels in the row using the bspline fit and added to a first pass model background image. 

[![](https://live.staticflickr.com/65535/51236365319_8578f5f4f9_b.jpg)](https://live.staticflickr.com/65535/51236365319_8578f5f4f9_b.jpg)

It's essential that background fitting is accurate within the order locations, but not outside of these areas, so there is no cause for concern if the fit is poor at the edges and corners of the image.

A median-filter is applied to the image to remove the structure resulting from a bad row fits (light and dark lines along the x-axis).

[![](https://live.staticflickr.com/65535/51234884582_5fd181c063_b.jpg)](https://live.staticflickr.com/65535/51234884582_5fd181c063_b.jpg)

Finally, the modelled background image is subtracted from the original frame:

[![](https://live.staticflickr.com/65535/51236665480_2e269e7049_b.jpg)](https://live.staticflickr.com/65535/51236665480_2e269e7049_b.jpg)
