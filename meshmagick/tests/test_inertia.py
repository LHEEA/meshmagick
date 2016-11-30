#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import meshmagick.mmio as mmio
from meshmagick.mesh import Mesh
from math import pi, fabs

vertices, faces = mmio.load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)

def test_cylinder_inertia():
    xmin, xmax, ymin, ymax, zmin, zmax = cylinder.axis_aligned_bbox
    R = (xmax-xmin)/2.
    h = zmax-zmin
    
    V = pi*R**2*h
    
    Ixx = Iyy = V*(R**2+h**2/3)/4
    Izz = V*R**2/2
    
    inertia = cylinder.eval_plain_mesh_inertias(rho_medium=1.)
    
    assert fabs(inertia.xx - Ixx) < 1000
    assert fabs(inertia.zz - Izz) < 1000
    
