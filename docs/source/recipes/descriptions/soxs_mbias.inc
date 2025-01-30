A zero-second exposure will contain only read-noise, and \~half of the pixels within this Gaussian distribution centred around zero counts will always contain negative flux. To avoid negative counts, an offset *bias* voltage is applied at the amplifier stage so that even when no photons are detected, the A/D converter will always register a positive value. This bias-voltage offset must be accounted for in the data reduction process. 

The purpose of the [`soxs_mbias`](#soxspipe.recipes.soxs_mbias) recipe is to provide a master-bias frame that can be subtracted from science/calibration frames to remove the contribution of pixel counts resulting from the bias voltage.


