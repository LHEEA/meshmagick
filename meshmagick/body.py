#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import base_classes as bc


class RigidBody(bc._BaseBody):
    def __init__(self, name=None):
        super(RigidBody, self).__init__(name=name)
        if name is None:
            self._name = '_'.join((self._name, '(rigid)'))

        self._inertial_frame = self.add_fixed_frame('inertial_frame')
        self._mesh_frame = bc.null_frame



