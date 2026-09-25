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
        """*average-combine the stack, ignoring NaN pixels*

        **Return:**

        - ``combined_image`` -- a ``CCDData`` holding the NaN-ignoring mean of the
          stack. Pixels masked in every frame are masked in the result, and ``NCOMBINE``
          holds the number of frames.
        """
        import bottleneck as bn
        import numpy as np
        from astropy.nddata import CCDData

        data, masked_values, scale_func = self._combination_setup(
            None, bn.nanmean, None
        )

        mean = scale_func(data, axis=0)
        mask = masked_values == len(self.data_arr)

        # CREATE THE COMBINED IMAGE WITH A DTYPE THAT MATCHES THE COMBINER
        combined_image = CCDData(
            np.asarray(mean, dtype=self.dtype), mask=mask, unit=self.unit
        )

        # UPDATE THE META DATA
        combined_image.meta["NCOMBINE"] = len(data)

        # RETURN THE COMBINED IMAGE
        return combined_image
