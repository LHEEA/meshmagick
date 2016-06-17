#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from mesh import *
from nodes import *
import numpy as np

# import sys


class InertiaTensor(np.ndarray):

    def __new__(cls, name=None):
        obj = np.zeros((6, 6), dtype=np.float).view(cls)
        obj._node = Node([0, 0, 0])
        obj._name = name
        return obj




    @property
    def rotational_inertia_matrix(self):
        return self[3:, 3:]

    # @property
    # def Ixx(self):
    #     return self.__internals__['inertia_matrix']







if __name__ == '__main__':

    import mmio

    vertices, faces = mmio.load_VTP('SEAREV.vtp')
    mymesh = Mesh(vertices, faces)

    # inertia = Inertia(mymesh)

