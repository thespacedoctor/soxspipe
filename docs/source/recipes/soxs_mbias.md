# soxs_mbias - COMPLETED

A zero-second exposure will contain only read-noise and \~half of pixels within this Gaussian distribution centred around zero count will always contain negative flux. To avoid negative counts an offset *bias* voltage is applied at the amplifier stage so that even when no photons are detected the A/D converter will always register a positive value. This bias-voltage offset must be accounted for in the data reduction process. 

The purpose of the [`soxs_mbias`](../_api/soxspipe.recipes.soxs_mbias.html) recipe is to stack raw bias-frames together (using the [`clip_and_stack`](../utils/clip_and_stack.md) utility) into [master-bias frames](../files/master_bias.md) and in the process clipping rogue pixels from the individual raw frames and reducing the read-noise contribution.

The master-bias frame can be subtracted from science/calibration frames to remove the contribution of pixel counts resulting from the bias-voltage.

![](soxs_mbias.png)

```eval_rst
.. autoclass:: soxspipe.recipes.soxs_mbias
    :members:
```


