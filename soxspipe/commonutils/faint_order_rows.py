"""How the global trace fit treats the synthetic rows of a transferred faint order.

The rows that ``faint_order_trace`` writes for a faint order lie exactly on the trace, so they fit well by construction.
They count toward the polynomial fit, but they must not set the clipping spread, be clipped, or pass for evidence of
a good fit in the residual statistics. A table without the ``faint_transfer`` column takes the unchanged paths.
"""

import numpy as np

TRANSFER_COLUMN = "faint_transfer"


def transferred_row_mask(table):
    """*mark the synthetic rows of a table*

    **Key Arguments:**

    - ``table`` -- a data-frame of trace samples, with or without a ``faint_transfer`` column

    **Return:**

    - ``isTransferred`` -- one boolean per row, or *None* when the table has no ``faint_transfer`` column

    A row is synthetic only when its flag is *True*. A missing flag (NaN, *None* or ``NA``) and any other value are
    measured rows.
    """
    if TRANSFER_COLUMN not in table.columns:
        return None
    return table[TRANSFER_COLUMN].eq(True).fillna(False).to_numpy(dtype=bool)


def clip_measured_residuals(residuals, isTransferred, sigmaLow, sigmaHigh):
    """*sigma-clip the residuals of the measured rows, and never clip the synthetic rows*

    The clipping centre and spread come from the measured rows only. Without a ``isTransferred`` mask the call is the
    plain sigma clip of all the residuals.

    **Key Arguments:**

    - ``residuals`` -- the fit residual of each row
    - ``isTransferred`` -- the mask of synthetic rows from ``transferred_row_mask``, or *None*
    - ``sigmaLow`` -- the lower clipping limit, in standard deviations
    - ``sigmaHigh`` -- the upper clipping limit, in standard deviations

    **Return:**

    - ``isClipped`` -- one boolean per row, *True* where the row is clipped
    """

    from astropy.stats import sigma_clip

    def clip(values):
        return sigma_clip(
            values, sigma_lower=sigmaLow, sigma_upper=sigmaHigh, maxiters=2, cenfunc="mean", stdfunc="std"
        )

    if isTransferred is None:
        return clip(residuals).mask
    isClipped = np.zeros(len(residuals), dtype=bool)
    isMeasured = ~isTransferred
    if isMeasured.any():
        isClipped[isMeasured] = clip(residuals[isMeasured]).mask
    return isClipped


def measured_residuals(residuals, isTransferred):
    """*the residuals the fit statistics are taken over: those of the measured rows*

    **Key Arguments:**

    - ``residuals`` -- the fit residual of each row
    - ``isTransferred`` -- the mask of synthetic rows from ``transferred_row_mask``, or *None*

    **Return:**

    - ``measured`` -- the residuals of the measured rows, or of all rows when the table has no synthetic rows or no
      measured rows
    """
    if isTransferred is None or isTransferred.all():
        return residuals
    return residuals[~isTransferred]
