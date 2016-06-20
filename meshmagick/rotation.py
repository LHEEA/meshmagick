#!/usr/bin/env python
#  -*- coding: utf-8 -*-
"""
This module defines many rotation representations.
"""


import numpy as np


class Rotation(object):
    """Base class for rotations"""
    def __init__(self):
        pass

    @property
    def matrix(self):
        return None
