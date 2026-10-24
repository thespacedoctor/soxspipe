"""Anchor the trace of a faint VIS order to its bright neighbour and refit it on binned profiles.

The slice-by-slice Gaussian fit of ``detect_continuum`` follows a bright order well and a faint order badly: at a
signal-to-noise ratio near 0.5 it fits noise, and a degree-5 polynomial through that noise wanders by arcseconds. The
functions here find the faint orders, take the slit position of the adjacent bright order as an anchor, collapse the
faint order into a few wavelength bins and fit one fixed-width Gaussian per bin inside a narrow window around the
anchor. The module needs numpy and pandas only, and runs no least-squares optimiser.
"""

from dataclasses import dataclass

# CODE FALLBACKS FOR SETTINGS FILES WRITTEN BEFORE THE FAINT-ORDER KEYS EXISTED. THE TRACE IS OFF BY DEFAULT, SO OLD
# SETTINGS FILES AND THE -STD RECIPES DO NOT CHANGE
FAINT_TRACE_DEFAULTS = {
    "faint-order-trace": False,
    "faint-order-min-detection-fraction": 0.5,
    "faint-order-max-scatter-px": 1.5,
    "faint-order-bin-count": 15,
    "faint-order-min-bin-snr": 3.0,
    "same-source-max-join-arcsec": 2.0,
}


@dataclass(frozen=True)
class faint_trace_settings:
    """The ``detect-continuum`` settings of the faint-order trace.

    ``isEnabled`` switches the whole feature. ``minDetectionFraction`` and ``maxScatterPx`` flag a faint order,
    ``binCount`` and ``minBinSnr`` shape the binned fit, and ``maxJoinArcsec`` is the same-source warning limit.
    """

    isEnabled: bool
    minDetectionFraction: float
    maxScatterPx: float
    binCount: int
    minBinSnr: float
    maxJoinArcsec: float

    @classmethod
    def from_block(cls, block):
        """*read the faint-order settings from a ``detect-continuum`` block, with code fallbacks*

        **Key Arguments:**

        - ``block`` -- the ``detect-continuum`` settings dictionary, or *None*

        **Return:**

        - ``settings`` -- the frozen settings, off when the block has no ``faint-order-trace`` key
        """
        values = {**FAINT_TRACE_DEFAULTS, **(block or {})}
        return cls(
            isEnabled=bool(values["faint-order-trace"]),
            minDetectionFraction=float(values["faint-order-min-detection-fraction"]),
            maxScatterPx=float(values["faint-order-max-scatter-px"]),
            binCount=int(values["faint-order-bin-count"]),
            minBinSnr=float(values["faint-order-min-bin-snr"]),
            maxJoinArcsec=float(values["same-source-max-join-arcsec"]),
        )
