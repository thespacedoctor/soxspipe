#!/usr/bin/env python
# encoding: utf-8
"""
*Given a keyword token and instrument name return the exact FITS Header keyword*

:Author:
    David Young & Marco Landoni

:Date Created:
    February 26, 2020
"""
################# GLOBAL IMPORTS ####################
from builtins import object
import sys
import os
os.environ['TERM'] = 'vt100'
from fundamentals import tools


class keyword_lookup(object):
    """
    *The worker class for the keyword_lookup module*

    **Key Arguments:**
        - ``log`` -- logger
        - ``settings`` -- the settings dictionary

    **Usage**

    To initalise the keyword lookup object in your code add the following:

    ```python
    from soxspipe.commonutils import keyword_lookup
    kw = keyword_lookup(
        log=log,
        settings=settings
    ).get
    ```

    After this it's possible to either look up a single keyword using it's alias:

    ```python
    kw("DET_NDITSKIP")
    > "ESO DET NDITSKIP"
    ```

    or return a list of keywords:

    ```python
    kw(["PROV", "DET_NDITSKIP"])
    > ['PROV', 'ESO DET NDITSKIP']
    ```

    For those keywords that require an index it's possible to also pass the index to the `kw` function:

    ```python
    kw("PROV", 9)
    > 'PROV09'
    ```

    If a tag is not in the list of FITS Header keyword aliases in the configuration file a `LookupError` will be raised.
    """
    # Initialisation

    def __init__(
            self,
            log,
            settings=False,

    ):
        self.log = log
        log.debug("instansiating a new 'keyword_lookup' object")
        self.settings = settings
        # xt-self-arg-tmpx

        # SELECT THE INSTRUMENT AND READ THE KEYWORD DICTIONARY IN RESOURCES
        # FOLDER
        if "instrument" in settings:
            self.instrument = settings["instrument"]
        else:
            self.instrument = "soxs"
        self.kwDict = self._select_dictionary()

        return None

    def get(self,
            tag,
            index=False):
        """
        *given a tag, and optional keyword index, return the FITS Header keyword for the selected instrument*

        **Key Arguments:**
            - ``tag`` -- the keyword tag as set in the yaml keyword dictionary (e.g. 'SDP_KEYWORD_TMID' returns 'TMID'). Can be string or list of sttings.
            - ``index`` -- add an index to the keyword if not False (e.g. tag='PROV', index=3 returns 'PROV03') Default *False*

        **Return:**
            - ``keywords`` -- the FITS Header keywords. Can be string or list of sttings depending on format of tag argument

        **Usage**

        See docstring for the class
        """
        self.log.debug('starting the ``get`` method')

        # CONVERT STRING TO LIST OF ONE ITEM
        single = False
        if not isinstance(tag, list):
            single = True
            tag = [tag]

        # STRINGIFY INDEX
        if index:
            index = "%(index)0.2d" % locals()
        else:
            index = ""

        # LOOKUP KEYWORDS
        keywords = []
        for t in tag:
            if t not in self.kwDict:
                raise LookupError(
                    "%(tag)s is not in the list of known FITS Header keyword aliases" % locals())
            keywords.append(self.kwDict[t] + index)

        # RETURNING A SINGLE KEYWORD?
        if single:
            keywords = keywords[0]

        self.log.debug('completed the ``get`` method')
        return keywords

    def _select_dictionary(
            self):
        """*select the keyword dictionary based on the instrument passed via the settings*

        **Return:**
            - ``kwDict`` -- the python dictionary of keywords (key = tag, value = fits keyword)

        **Usage**

        ```python
        from soxspipe.commonutils import keyword_lookup
        this = keyword_lookup(
            log=log,
            settings=settings
        )
        kwDict = this._select_dictionary()
        ```
        """
        self.log.debug('starting the ``_select_dictionary`` method')

        # GENERATE PATH TO YAML DICTIONARY
        yamlFilePath = os.path.dirname(os.path.dirname(
            __file__)) + "/resources/" + self.instrument + "_keywords.yaml"

        # YAML CONTENT TO DICTIONARY
        import yaml
        with open(yamlFilePath, 'r') as stream:
            kwDict = yaml.load(stream)

        self.log.debug('completed the ``_select_dictionary`` method')
        return kwDict
