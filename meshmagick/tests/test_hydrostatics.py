#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import mmio
import hydrostatics as hs
from mesh import Mesh
from math import pi, fabs


# Importing cylinder mesh
vertices, faces = mmio.load_VTP('tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)

vertices, faces = mmio.load_VTP('tests/data/SEAREV.vtp')
searev = Mesh(vertices, faces)

def test_cylinder_base_hydrostatics():

    output = hs.compute_hydrostatics(cylinder, 0)

    # Analytical results
    D = 10
    h = 10
    flottation_surface = pi*D*D/4
    immersed_volume = flottation_surface*h
    print
    wetted_surface = flottation_surface + pi*D*h

    assert fabs(output['Sf'] - flottation_surface) < 1
    assert fabs(output['Sw'] - wetted_surface) < 2
    assert fabs(output['Vw'] - immersed_volume) < 10


def test_searev_base_hydrostatics():
    output = hs.compute_hydrostatics(searev, 0)

    assert fabs(output['Sf'] - 300) < 1
    assert fabs(output['Sw'] - 550) < 1
    assert fabs(output['Vw'] - 1177) < 1


def test_displacement_equilibrium():
    # rho_water = 1023
    # D = 10
    # h = 12
    # disp = rho_water * pi*D*D*h/4 *1e-3
    
    disp = 2000
    CG = [-1, 0, -3]
    
    cylinder_eq, CGc = hs.compute_equilibrium(searev, disp, CG, rho_water=1023, anim=True)
    
    # print cylinder_eq.axis_aligned_bbox
    # output = hs.compute_hydrostatics(cylinder_eq, CG[-1], verbose=True)
