#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *
import numpy as np

import sys




class Inertia(object):
    def __init__(self, mesh):
        self._mesh = mesh

        self.__internals__ = dict()

        self._compute_surface_integrals_vectorized()


    def reset(self):
        self.__internals__.clear()








if __name__ == '__main__':

    import mmio

    vertices, faces = mmio.load_VTP('SEAREV.vtp')
    mymesh = Mesh(vertices, faces)

    inertia = Inertia(mymesh)

