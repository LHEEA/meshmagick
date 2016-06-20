#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np


class Node(np.ndarray):

    def __new__(cls, name, coords):
        obj = np.asarray(coords, dtype=np.float).view(cls)
        obj._name = str(name)
        return obj

    def __str__(self):
        str_repr = "Node: '%s'\nCoords: %.3f, %.3f, %.3f" % (self._name, self[0], self[1], self[2])
        return str_repr

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, value):
        self._name = str(value)
        return

    @property
    def x(self):
        return self[0]

    @x.setter
    def x(self, value):
        self[0] = float(value)
        return
    @property
    def y(self):
        return self._coords[1]

    @y.setter
    def y(self, value):
        self._coords[1] = float(value)
        return

    @property
    def z(self):
        return self._coords[2]

    @z.setter
    def z(self, value):
        self._coords[2] = float(value)
        return
