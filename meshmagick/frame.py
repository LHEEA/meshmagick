#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import numpy as np
from copy import deepcopy
from itertools import count
from Singleton import _Singleton
# from body import * # FIXME: Voir pourquoi l'import selectif ne fonctionne pas!!
import body

class Frame(object):
    _ids = count(0)

    def __init__(self, name=None, parent_frame='ground'):

        self._id = self._ids.next()
        if not name:
            name = '_'.join(('Frame', str(self._id)))

        self._name = str(name)
        self.parent_frame = parent_frame

    def __str__(self):
        str_repr = "Frame : '%s'" % self._name
        return str_repr

    __repr__ = __str__

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, _name):
        self._name = str(_name)

    @property
    def parent_frame(self):
        return self._parent_frame

    @parent_frame.setter
    def parent_frame(self, _parent_frame):
        if not isinstance(_parent_frame, Frame):
            if _parent_frame == 'ground':
                _parent_frame = ground_frame
            elif _parent_frame is None:
                pass
            else:
                raise TypeError, "parent_frame of Frame must be a Frame object"
        self._parent_frame = _parent_frame
        return


class NullFrame(Frame):
    __metaclass__ = _Singleton

    def __init__(self):
        super(NullFrame, self).__init__(name='Null_frame', parent_frame=None)


null_frame = NullFrame()


class GroundFrame(Frame):
    __metaclass__ = _Singleton

    def __init__(self):
        super(GroundFrame, self).__init__(name='Ground_frame', parent_frame=self)
        self._parent_frame = null_frame # TODO: Definir une classe NullFrame afin de rester consistant ?

    @property
    def name(self):
        return self._name


ground_frame = GroundFrame()





class BodyFixedFrame(Frame):
    def __init__(self, name=None, parent_frame='ground', parent_body=None): # TODO: Remplacer le None par un
        super(BodyFixedFrame, self).__init__(name=name, parent_frame=parent_frame)
        if parent_body is None:
            parent_body = body.null_body
        self.parent_body = parent_body

    @property
    def parent_body(self):
        return self._parent_body

    @parent_body.setter
    def parent_body(self, _parent_body):
        if not isinstance(_parent_body, body.GenericBody):
            raise TypeError, "A parent body to attach to a Frame must be a GenericBody object"
        _parent_body.attach_frame(self)
        self._parent_body = _parent_body
        return





# if __name__ == '__main__':
#
#     mybody = body.GenericBody()
#
#     frame = BodyFixedFrame(parent_body=mybody)
#
#     print mybody.get_frames()



