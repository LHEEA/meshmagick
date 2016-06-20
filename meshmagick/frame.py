#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np


class Frame(object):

    def __init__(self, name, parent_frame=None, transform=None):

        self._name = str(name)

        if not parent_frame:
            self._parent = ground
        else:
            if not isinstance(parent_frame, Frame):
                raise ValueError, "A Frame parent frame must be a Frame object"
            self._parent = parent_frame

    def __str__(self):
        str_repr = "%s" % self._name
        return str_repr


class GroundFrame(Frame):
    def __init__(self):
        self._name = 'Ground_reference_frame'
        pass

    # TODO: overloader toutes les methodes afin de rendre impossible la modif des pptes




ground = GroundFrame()
