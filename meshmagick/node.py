#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy

# Node doit etre une fabrique

_ground_node_key_name = '___________'


class Node(np.ndarray):
    def __new__(cls, name, coords, parent='ground_node', fixed=True):
        if len(coords) != 3:
            raise ValueError, "Node coordinates must be an array like with 3 elements"
        obj = np.asarray(coords, dtype=np.float).view(cls)
        if not isinstance(parent, Node):
            if parent == 'ground_node':
                parent = GroundNode()
            else:
                if name == _ground_node_key_name:
                    return obj
                else:
                    raise TypeError, "A Node parent must be a Node object"
        obj._name = str(name)
        obj._parent = parent
        return obj

    def __str__(self):
        str_repr = "Node: '%s'\nCoords: %.3f, %.3f, %.3f" % (self._name, self[0], self[1], self[2])
        return str_repr

    __repr__ = __str__

    @property
    def parent(self):
        return self._parent

    @parent.setter
    def parent(self, _parent):
        if not isinstance(_parent, Node):
            raise TypeError, "Node parent must be a Node object."
        self._parent = _parent
        return

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
        return self[1]

    @y.setter
    def y(self, value):
        self[1] = float(value)
        return

    @property
    def z(self):
        return self[2]

    @z.setter
    def z(self, value):
        self[2] = float(value)
        return

    def copy(self):
        obj = deepcopy(self)
        obj._name = '_'.join((self._name, 'copy'))
        obj._parent = self._parent
        return obj

    def set_coords(self, coords):
        if len(coords) != 3:
            raise ValueError, "Node coordinates must be an array like with 3 elements"
        self[:] = np.asarray(coords, dtype=np.float)
        return

    def is_absolute(self):
        return isinstance(self._parent, GroundNode)

    def is_relative(self):
        return not self.is_absolute()

    def get_absolute_position(self):
        # CA NE FONCTIONNE PAS !!!!
        node = self
        position = self.asarray()
        while True:
            if isinstance(node._parent, GroundNode):
                break
            node = self._parent
            position



    def asarray(self):
        return np.asarray(self)

    def __add__(self, other):
        self[:] += other
        return self

    def __radd__(self, other):
        self[:] += other
        return self

    def __sub__(self, other):
        self[:] -= np.asarray(other)
        return self

    def __rsub__(self, other):
        self[:] *= -1
        self[:] += np.asarray(other)
        return self

    def __mul__(self, other):
        self[:] *= np.asarray(other)

class GroundNode(Node):
    """
    Reference Node for the world. Coordinates 0., 0., 0. and immutable

    It uses the Singleton design pattern to ensure a unique instance.
    """
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = Node.__new__(Node, _ground_node_key_name, [0., 0., 0.], None).view(cls)
            cls._instance._name = 'Ground_Node'
            cls._instance._parent = None
            cls._instance.setflags(write=False)
        return cls._instance
        # obj.set_option(write=False)

    @property
    def name(self):
        return 'Ground_Node'

    @property
    def parent(self):
        return None

    def is_absolute(self):
        return True

    def is_relative(self):
        return False

if __name__ == '__main__':
    mynode0 = Node('node', [0, 0, 0])

    # print mynode0 + [1, 2, 3]
    # print [1, 2, 3] + mynode0

    # mynode1 = Node('node1', np.random.rand(3), parent=mynode0)

    print [0, 1, 2] - mynode0
    # print mynode0 - [0, 1, 2]
    # a = mynode1 - [0, 1, 2]
    # print a
