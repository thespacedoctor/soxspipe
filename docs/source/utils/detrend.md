# detrend


For the UVB-VIS arm, we will often need to scale the master-dark frame to match the exposure time of our science/calibration frame. If this is the case, the master-bias frame needs to be subtracted from the master-dark frame before scaling. However, if the master-dark frame has the same exposure time as your science/calibration frame, it can be subtracted directly from the frame, as this serves to remove the bias and dark-current contributions simultaneously. 

This logic is all housed within the [`detrend`](#soxspipe.recipes.base_recipe.base_recipe.detrend) method.


:::{figure-md} detrend
![](detrend.png){width=600px}

The algorithm used to detrend raw images.
:::


## Utility API

:::{autodoc2-object} soxspipe.recipes.base_recipe.base_recipe.detrend
:::
