#!/usr/bin/env python
"""
*Sub-class of CCDProc Combiner to fix error map combination*

Author
: David Young

Date Created
: October 27, 2022
"""

import os

from ccdproc import Combiner as OriginalCombiner

os.environ["TERM"] = "vt100"


class Combiner(OriginalCombiner):

    def average_combine(self):
        """ """
        import bottleneck as bn
        import numpy as np
        from astropy.nddata import CCDData

        data, masked_values, scale_func = self._combination_setup(
            None, bn.nanmean, None
        )

        mean = scale_func(data, axis=0)
        mask = masked_values == len(self.data_arr)

        # create the combined image with a dtype that matches the combiner
        combined_image = CCDData(
            np.asarray(mean, dtype=self.dtype), mask=mask, unit=self.unit
        )

        # update the meta data
        combined_image.meta["NCOMBINE"] = len(data)

        # return the combined image
        return combined_image
