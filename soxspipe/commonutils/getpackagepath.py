#!/usr/bin/env python
# encoding: utf-8
"""
*Get common file and folder paths for the host package*
"""
import os


def getpackagepath():
    """
    *Get the root path for this python package*

    Used in unit testing code
    """
    moduleDirectory = os.path.dirname(__file__)
    packagePath = os.path.dirname(__file__) + "/../"

    return packagePath
