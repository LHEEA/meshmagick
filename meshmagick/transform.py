#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
import base_classes as bc



class _Rotation(object):
    def __init__(self, rotation_matrix):
        self._matrix = np.array(rotation_matrix, dtype=np.float)





class _Rotation_axis_angle(object):
    def __init__(self, axis, angle):
        


class Rotation(object):
    __rotation_classes = {
        'axis_angle': _Rotation_axis_angle
    }

    @staticmethod
    def get_rotation(name, *args, **kwargs):
        rotation_class = Rotation.__rotation_classes.get(name.lower(), None)

        if rotation_class:
            return rotation_class(*args, **kwargs)
        raise NotImplemented("The requested rotation representation is not implemented")



class Transform(object):
    pass


if __name__ == '__main__':

    rot = Rotation.get_rotation('axis_angle', 1, 2)



