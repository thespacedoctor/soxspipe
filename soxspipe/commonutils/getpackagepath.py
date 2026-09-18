#!/usr/bin/env python
"""
*Get common file and folder paths for the host package*
"""

import os


def getpackagepath():
    """
    *Get the root path for this python package*

    Used in unit testing code

    **Return:**

    - ``packagePath`` -- the path to the package root, ending in ``/../`` (relative to this module's directory)
    """
    moduleDirectory = os.path.dirname(__file__)
    packagePath = os.path.dirname(__file__) + "/../"

    return packagePath
