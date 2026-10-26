"""The dispersion-map access of the faint-order trace.

``faint_order_trace`` needs only numpy and pandas, and takes a ``pixelsAt`` function wherever it needs the map. This
module builds that function from a dispersion map on disk.
"""

import numpy as np
import pandas as pd

from soxspipe.commonutils.dispersion_map_to_pixel_arrays import dispersion_map_to_pixel_arrays


def make_pixels_at(log, dispersionMap, mapBinning, axisA, axisB):
    """*build the function that turns an order, wavelengths and a slit position into frame pixels*

    The dispersion map gives pixels at the map's binning, so they are divided by ``mapBinning`` to land on the frame.
    Samples off the detector are kept, because a search window can reach past the edge of an order.

    **Key Arguments:**

    - ``log`` -- logger
    - ``dispersionMap`` -- path to the dispersion map
    - ``mapBinning`` -- the frame's binning over the map's binning, ``{"x": ..., "y": ...}``
    - ``axisA`` -- the cross-dispersion axis, ``"x"`` or ``"y"``
    - ``axisB`` -- the dispersion axis, ``"x"`` or ``"y"``

    **Return:**

    - ``pixelsAt`` -- the function ``pixelsAt(order, wavelengths, slitArcsec)``, giving the (axisA, axisB) pixels
    """

    def pixels_at(order, wavelengths, slitArcsec):
        wavelengths = np.asarray(wavelengths, dtype=float)
        table = pd.DataFrame(
            {
                "order": float(order),
                "wavelength": wavelengths,
                "slit_position": np.broadcast_to(np.asarray(slitArcsec, dtype=float), wavelengths.shape),
            }
        )
        table = dispersion_map_to_pixel_arrays(
            log=log, dispersionMapPath=dispersionMap, orderPixelTable=table, removeOffDetectorLocation=False
        )
        return (
            table[f"fit_{axisA}"].to_numpy(float) / mapBinning[axisA],
            table[f"fit_{axisB}"].to_numpy(float) / mapBinning[axisB],
        )

    return pixels_at
