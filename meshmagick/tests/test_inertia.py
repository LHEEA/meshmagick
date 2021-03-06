#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import meshmagick.mmio as mmio
from meshmagick.mesh import Mesh
from math import pi, fabs
import pytest

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
    
    assert pytest.approx(Ixx, rel=1e-1) == inertia.xx
    assert pytest.approx(Izz, rel=1e-1) == inertia.zz

def test_move_cog():
    xmin, xmax, ymin, ymax, zmin, zmax = cylinder.axis_aligned_bbox
    R = (xmax-xmin)/2.
    h = zmax-zmin
    
    V = pi*R**2*h
    
    Ixx = V*(R**2+h**2/3)/4
    Ixx_base = Ixx + V*(h/2)**2             # parallel axis theorem
    Ixx_base2 = 1/4*V*R**2 + 1/3*V*h**2

    inertia = cylinder.eval_plain_mesh_inertias(rho_medium=1.)
    inertia.set_cog((0,0,-h/2))
    inertia.shift_at_cog()

    assert pytest.approx(Ixx_base, rel=1e-1) == inertia.xx
