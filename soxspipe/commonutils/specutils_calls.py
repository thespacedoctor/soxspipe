#!/usr/bin/env python
"""
*specutils calls with the warnings that come from specutils' own internals ignored for the call only*

Each helper ignores one known message, from one library, around one call. Nothing here changes
``warnings.filters`` once the call returns. A warning that is about the data stays visible.

Author
: David Young

Date Created
: November 19, 2026
"""

from typing import Any


def spectrum1d_quietly(**kwargs: Any) -> Any:
    """*build a* ``specutils.Spectrum1D`` *without gwcs' deprecation notice*

    **Key Arguments:**

    - ``kwargs`` -- the keyword arguments given to ``Spectrum1D``

    **Return:**

    - ``spectrum`` -- the ``Spectrum1D``

    **Usage:**

    ```python
    from soxspipe.commonutils.specutils_calls import spectrum1d_quietly
    spectrum = spectrum1d_quietly(flux=flux, spectral_axis=wavelength)
    ```
    """
    import warnings

    from astropy.utils.exceptions import AstropyDeprecationWarning
    from specutils import Spectrum1D

    with warnings.catch_warnings():
        # GWCS HOLDS THE SPECTRAL AXIS AND CALLS AN ASTROPY FUNCTION THAT ASTROPY DEPRECATED, SO SOXSPIPE CANNOT FIX IT
        warnings.filterwarnings(
            "ignore", message="The isiterable function is deprecated", category=AstropyDeprecationWarning
        )
        return Spectrum1D(**kwargs)


def fit_generic_continuum_quietly(spectrum: Any) -> Any:
    """*fit the generic continuum of a spectrum without astropy's notes about the fitter specutils chose*

    **Key Arguments:**

    - ``spectrum`` -- the ``Spectrum1D`` to fit

    **Return:**

    - ``continuum`` -- the fitted continuum model, as ``specutils.fitting.fit_generic_continuum`` returns it

    **Usage:**

    ```python
    from soxspipe.commonutils.specutils_calls import fit_generic_continuum_quietly
    continuum = fit_generic_continuum_quietly(spectrum)
    ```
    """
    import warnings

    from astropy.utils.exceptions import AstropyUserWarning
    from specutils.fitting import fit_generic_continuum

    from soxspipe.commonutils.toolkit import ignore_astropy_covariance_scaling_overflow

    with warnings.catch_warnings():
        # THE CONTINUUM IS A POLYNOMIAL, WHICH IS LINEAR IN ITS PARAMETERS, AND SPECUTILS PICKS THE NON-LINEAR FITTER.
        # THAT FITTER IS UNWEIGHTED, SO IT ALSO SCALES THE UNREAD COVARIANCE, WHICH CAN OVERFLOW
        warnings.filterwarnings("ignore", message="Model is linear in parameters", category=AstropyUserWarning)
        ignore_astropy_covariance_scaling_overflow()
        return fit_generic_continuum(spectrum)
