#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy
from itertools import count
import body

# _ground_frame_key_name = '___________'

class Frame(object):
    _ids = count(0)

    def __init__(self, name=None, parent_frame='ground'):

        self._id = self._ids.next()
        if not name:
            name = '_'.join(('Frame', str(self._id)))
        self._name = str(name)

        if not isinstance(parent_frame, Frame):
            if parent_frame == 'ground':
                parent_frame = ground
            else:
                raise TypeError, "parent_frame of a Frame must be a Frame object"
        self._parent_frame = parent_frame

    def __str__(self):
        str_repr = "Frame : '%s'" % self._name
        return str_repr

    __repr__ = __str__


class _Singleton(type):
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(_Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]

class GroundFrame(Frame):
    __metaclass__ = _Singleton

    def __init__(self):
        Frame.__init__(self, name='Ground_frame', parent_frame=self)
        self._parent_frame = None

    @property
    def name(self):
        return self._name


ground = GroundFrame()


class BodyFixedFrame(Frame):
    def __init__(self, parent_body, name=None, parent_frame='ground'):
        super(BodyFixedFrame, self).__init__(name=name, parent_frame=parent_frame)

        if not isinstance(parent_body, body.GenericBody):
            raise TypeError, "parent_body of BodyFixedFrame must be a Body instance"
        self._parent_body = parent_body

if __name__ == '__main__':

    mybody = body.GenericBody()
    frame = BodyFixedFrame(mybody)

    print frame



