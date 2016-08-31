#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import mmio
import mesh

vertices, faces = mmio.load_VTP('tests/data/Cylinder.vtp')
cylinder = mesh.Mesh(vertices, faces)

def test_bbox():
    assert (-5, 5, -5, 5, -10, 10) == cylinder.axis_aligned_bbox
    assert (-10, 10, -10, 10, -10, 10) == cylinder.squared_axis_aligned_bbox
    
