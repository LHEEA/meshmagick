#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import frame
import rotation as rot
import numpy as np


class Transform(object):
    """Base class for frame transformations"""
    def __init__(self, target_frame, rotation=None, translation=None, root_frame=None):
        if not isinstance(target_frame, frame.Frame):
            raise ValueError, "A Transform target frame must be a BaseFrame object"

        self._target_frame = target_frame

        if not root_frame:
            self._root_frame = frame.ground_frame
        else:
            if not isinstance(root_frame, frame.Frame):
                raise ValueError, "A Transform root frame must be a BaseFrame object"
            self._root_frame = root_frame

        if not rotation:
            rotation = rot.Rotation()
        else:
            # TODO : FINIR !!!
            if not isinstance(rotation, rot.Rotation):
                # TODO: construire un objet rotation ?
                pass
            else:
                pass

        if not translation:
            translation = np.zeros(3)

        self._rotation = rotation
        self._translation = translation



    # TODO: definir des operateurs !! *, +, /, \... plutot que directement des methodes !! --> plus succint
    # def motion_transform(self, motion_vector):
    #     pass
    #
    # def force_transform(self, force_vector):
    #     pass
    #
    # def inverse_motion_tranform(self, motion_vector):
    #     pass
    #
    # def inverse_force_transform(self, force_vector):
    #     pass

