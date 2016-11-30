#!/usr/bin/env python
#  -*- coding: utf-8 -*-

from meshmagick.mmio import load_VTP
from meshmagick.mesh import Mesh

from os.path import curdir

vertices, faces = load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)

def test_bbox():
    assert (-5, 5, -5, 5, -10, 10) == cylinder.axis_aligned_bbox
    assert (-10, 10, -10, 10, -10, 10) == cylinder.squared_axis_aligned_bbox
