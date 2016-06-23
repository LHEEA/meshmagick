#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from itertools import count
import frame
from Singleton import _Singleton

# TODO: mettre ground_frame, ground_body, null_frame et null_body dans un module dedie... Cela permettra de regler un
#  certain nombre de problemes de precedences d'import (peut-etre :/)
# Le module pourra s'appeler 'base' et comportera les types de base du framework dont GenericBody et Frame qu'on
# pourra renommer BaseFrame... 

class GenericBody(object):
    _ids = count(0)
    def __init__(self, name=None):
        self._id = self._ids.next()
        if not name:
            name = '_'.join(('Body', str(self._id)))
        self._name = str(name)

        self._attached_frames = list()
        pass

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, _name):
        self._name = str(_name)
        return


    def attach_frame(self, new_frame):
        if not isinstance(new_frame, frame.BodyFixedFrame):
            raise TypeError, "A frame attached to a GenericBody must be a Frame object"
        self._attached_frames.append(new_frame)
        return

    def get_frames(self):
        return self._attached_frames


class NullBody(GenericBody):
    """
    A body class representing an undefined body
    """
    __metaclass__ = _Singleton

    def __init__(self):
        super(NullBody, self).__init__(name='null_body')
        GenericBody._ids = count(self._id)
        self._id = -1

    @property
    def name(self):
        return self._name


null_body = NullBody()


class GroundBody(GenericBody):
    __metaclass__ = _Singleton

    def __init__(self):
        super(GroundBody, self).__init__(name='Ground_body')
        self._reference_frame = frame.ground_frame

