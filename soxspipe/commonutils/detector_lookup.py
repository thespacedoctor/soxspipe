#!/usr/bin/env python
# encoding: utf-8
"""
*return a dictionary of detector characteristics and parameters*

Author
: David Young & Marco Landoni

Date Created
: August 13, 2020
"""
################# GLOBAL IMPORTS ####################
from fundamentals import tools
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'


class detector_lookup(object):
    """
    *return a dictionary of detector characteristics and parameters*

    **Key Arguments:**

    - ``log`` -- logger
    - ``settings`` -- the settings dictionary

    **Usage:**

    To initiate a detector_lookup object, use the following:

    ```python
    from soxspipe.commonutils import detector_lookup
    detector = detector_lookup(
        log=log,
        settings=settings
    ).get("NIR")
    print(detector["science-pixels"])
    ```

    """

    def __init__(
            self,
            log,
            settings=False,
    ):
        self.log = log
        log.debug("instantiating a new 'detector_lookup' object")
        self.settings = settings

        # SELECT THE INSTRUMENT AND READ THE KEYWORD DICTIONARY IN RESOURCES
        # FOLDER
        if "instrument" in settings:
            self.instrument = settings["instrument"]
        else:
            self.instrument = "soxs"
        self.dectDict = self._select_dictionary()

        return None

    def get(self,
            arm):
        """
        *return a dictionary of detector characteristics and parameters*

        **Key Arguments:**

        - ``arm`` -- the detector parameters to return
        """
        self.log.debug('starting the ``get`` method')

        arm = arm.upper()

        if arm not in self.dectDict:
            raise LookupError(f"the detector '{arm}' cannot be found in the detector parameters lookup file")

        self.log.debug('completed the ``get`` method')
        return self.dectDict[arm]

    def _select_dictionary(
            self):
        """*select the detector parameter dictionary based on the instrument passed via the settings*

        **Return:**

        - ``dectDict`` -- the python dictionary of detector parameters

        **Usage**

        ```python
        from soxspipe.commonutils import detector_lookup
        detector = detector_lookup(
            log=log,
            settings=settings
        )
        dectDict = detector._select_dictionary()
        ```
        """
        self.log.debug('starting the ``_select_dictionary`` method')

        # GENERATE PATH TO YAML DICTIONARY
        yamlFilePath = os.path.dirname(os.path.dirname(
            __file__)) + "/resources/" + self.instrument + "_detector_parameters.yaml"

        # YAML CONTENT TO DICTIONARY
        import yaml
        with open(yamlFilePath, 'r') as stream:
            dectDict = yaml.safe_load(stream)

        self.log.debug('completed the ``_select_dictionary`` method')
        return dectDict
