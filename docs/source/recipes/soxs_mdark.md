# soxs_mdark - COMPLETED

Every raw CCD image contains counts resulting from a 'dark current', electrons released due to the thermal effects in the CCD material. For both the UVB-VIS (< 0.00012 e^-^/s/pixel) and NIR detectors (< 0.005 e^-^/s/pixel) the dark-current is almost negligible. Not all pixels will have the same dark-current, some will have a high than typical current. These are so-called 'hot-pixels' and it's important that these are identified and recorded (using the [`create_noise_map`](../utils/create_noise_map.md) utility).

The purpose of the [`soxs_mdark`](../_api/soxspipe.recipes.soxs_mdark.html) recipe is to stack raw dark-frames together (using the [`clip_and_stack`](../utils/clip_and_stack.md) utility) into [master-dark frames](../files/master_dark.md) and in the process clipping rogue pixels from the individual raw frames and reducing the read-noise contribution.

![](soxs_mdark.png)

```eval_rst
.. autoclass:: soxspipe.recipes.soxs_mdark
    :members:
```
