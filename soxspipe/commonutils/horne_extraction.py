#!/usr/bin/env python
# encoding: utf-8
"""
*perform optimal source extraction using the Horne method (Horne 1986)*

:Author:
    Marco Landoni & David Young

:Date Created:
    May 17, 2023
"""

# Please do import into the methods for speed reasons.

from fundamentals import tools
from builtins import object
import sys
import os
from astropy import units as u
from astropy.nddata import CCDData


os.environ['TERM'] = 'vt100'


class horne_extraction(object):
    """
    *The worker class for the horne_extraction module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary
        - ``skyModelFrame`` -- sky model frame
        - ``skySubFrame`` -- sky subtracted frame
        - ``recipeName`` -- name of the recipe as it appears in the settings dictionary
        - ``twoDMapPath`` -- 2D dispersion map image path
        - ``qcTable`` -- the data frame to collect measured QC metrics
        - ``productsTable`` -- the data frame to collect output products

    **Usage:**

    To setup your logger, settings and database connections, please use the ``fundamentals`` package (`see tutorial here <http://fundamentals.readthedocs.io/en/latest/#tutorial>`_). 

    To initiate a horne_extraction object, use the following:

    ```eval_rst
    .. todo::

        - add usage info
        - create a sublime snippet for usage
        - create cl-util for this class
        - add a tutorial about ``horne_extraction`` to documentation
        - create a blog post about what ``horne_extraction`` does
    ```

    ```python
    usage code 
    ```

    """

    def __init__(
            self,
            log,
            settings,
            skyModelFrame,
            skySubtractedFrame,
            twoDMapPath,
            recipeName=False,
            qcTable=False,
            productsTable=False,
            dispersionMap=False

    ):
        import numpy as np
        import pandas as pd

        self.log = log
        log.debug("instansiating a new 'horne_extraction' object")
        self.dispersionMap = dispersionMap
        self.settings = settings
        # xt-self-arg-tmpx

        self.skySubractedFrame = CCDData.read(skySubtractedFrame, hdu=0, unit=u.electron,
                                hdu_uncertainty='ERRS', hdu_mask='QUAL', hdu_flags='FLAGS',
                                key_uncertainty_type='UTYPE')
        from astropy.io import fits
        hdul = fits.open(twoDMapPath)

        # MAKE X, Y ARRAYS TO THEN ASSOCIATE WITH WL, SLIT AND ORDER
        xdim = hdul[0].data.shape[1]
        ydim = hdul[0].data.shape[0]

        xarray = np.tile(np.arange(0, xdim), ydim)
        yarray = np.repeat(np.arange(0, ydim), xdim)

        self.imageMap = pd.DataFrame.from_dict({
            "x": xarray,
            "y": yarray,
            "wavelength": hdul["WAVELENGTH"].data.flatten().byteswap().newbyteorder(),
            "slit_position": hdul["SLIT"].data.flatten().byteswap().newbyteorder(),
            "order": hdul["ORDER"].data.flatten().byteswap().newbyteorder(),
            "flux": self.skySubractedFrame.data.flatten().byteswap().newbyteorder()
        })

        self.imageMap.dropna(how="all", subset=["wavelength", "slit_position", "order"], inplace=True)



    def extract(self,order):
        """
        *get the horne_extraction object*

        **Return:**
            - ``horne_extraction``

        **Usage:**

        ```eval_rst
        .. todo::

            - add usage info
            - create a sublime snippet for usage
            - create cl-util for this method
            - update the package tutorial if needed
        ```

        ```python
        usage code 
        ```
        """
        self.log.debug('starting the ``extract`` method')

        from soxspipe.commonutils import detect_continuum
        import yaml
        import pandas as pd


        # DATAFRAMES TO COLLECT QCs AND PRODUCTS
        qc = pd.DataFrame({
            "soxspipe_recipe": [],
            "qc_name": [],
            "qc_value": [],
            "qc_unit": [],
            "qc_comment": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "to_header": []
        })
        products = pd.DataFrame({
            "soxspipe_recipe": [],
            "product_label": [],
            "file_name": [],
            "file_type": [],
            "obs_date_utc": [],
            "reduction_date_utc": [],
            "file_path": []
        })

        # HACK TO GET DETECT CONTINUUM TO RUN
        self.skySubractedFrame.header["ESO DPR TYPE"] = "LAMP,FLAT"
        self.skySubractedFrame.header["ESO DPR TECH"] = "IMAGE"
        self.settings["intermediate-data-root"] = "./"

        dct = yaml.safe_load('''
        soxs-order-centre:
            order-sample-count: 1000
            slice-length: 10
            peak-sigma-limit: 3
            disp-axis-deg: 5
            order-deg: 4
            poly-fitting-residual-clipping-sigma: 5.0
            poly-clipping-iteration-limit: 5
        ''')

        # MERGE DICTIONARIES (SECOND OVERRIDES FIRST)
        self.settings = {**self.settings, **dct}

        ## NOTE TO MARCO - I HAVE TRICKED THE CODE INTO THINKING THIS IS A LAMP-FLAT PINHOLE FRAME
        ## SO detect_continuum OUTPUTS 20190831T001327_NIR_MORDER_CENTRES_residuals.pdf AND 20190831T001327_NIR_ORDER_LOCATIONS.fits
        ## BUT THESE ARE REALLY THE OBJECT CONTINUUM NOT ORDER CENTRE TRACES

        detector = detect_continuum(
            log=self.log,
            pinholeFlat=self.skySubractedFrame,
            dispersion_map=self.dispersionMap,
            settings=self.settings,
            recipeName="soxs-order-centre",
            qcTable=qc,
            productsTable=products
        )


        productPath, qcTable, productsTable = detector.get()
        from soxspipe.commonutils.toolkit import unpack_order_table
        # UNPACK THE ORDER TABLE (CENTRE LOCATION ONLY AT THIS STAGE)
        orderPolyTable, orderPixelTable, orderMetaTable = unpack_order_table(
            log=log, orderTablePath=productPath)

        self.log.debug('completed the ``extract`` method')
        return None

    # xt-class-method
