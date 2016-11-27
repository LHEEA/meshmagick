#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import mmio
from mesh import Mesh
from math import pi
from mesh_clipper import MeshClipper

vertices, faces = mmio.load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)

def test_cylinder_inertia():
    xmin, xmax, ymin, ymax, zmin, zmax = cylinder.axis_aligned_bbox
    R = (xmax-xmin)/2.
    h = zmax-zmin
    
    V = pi*R**2*h
    
    Ixx = Iyy = V*(R**2+h**2/3)/4
    Izz = V*R**2/2
    
    inertia_data = cylinder.eval_plain_mesh_inertias(rho_medium=1.)
    coeffs = inertia_data['coeffs']
    print coeffs['Ixx'] ,  Ixx
    print coeffs['Izz'] , Izz
    

if __name__ == '__main__':
    test_cylinder_inertia()
