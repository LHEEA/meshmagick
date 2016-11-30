#!/usr/bin/env python
#  -*- coding: utf-8 -*-

import meshmagick.mmio as mmio
import meshmagick.hydrostatics as hs
from meshmagick.mesh import Mesh
from math import pi, fabs


# Importing cylinder mesh
vertices, faces = mmio.load_VTP('meshmagick/tests/data/Cylinder.vtp')
cylinder = Mesh(vertices, faces)

vertices, faces = mmio.load_VTP('meshmagick/tests/data/SEAREV.vtp')
searev = Mesh(vertices, faces)

def test_cylinder_base_hydrostatics():
    
    hydrostatics = hs.Hydrostatics(cylinder)
    hs_data = hydrostatics.hs_data
    
    # Analytical results
    D = 10
    h = 10
    flottation_surface = pi*D*D/4
    immersed_volume = flottation_surface*h
    wetted_surface = flottation_surface + pi*D*h

    assert fabs(hs_data['Aw'] - flottation_surface) < 1
    assert fabs(hs_data['Sw'] - wetted_surface) < 2
    assert fabs(hs_data['disp_volume'] - immersed_volume) < 10


def test_searev_base_hydrostatics():
    hydrostatics = hs.Hydrostatics(searev)
    hs_data = hydrostatics.hs_data

    assert fabs(hs_data['Aw'] - 300) < 1
    assert fabs(hs_data['Sw'] - 550) < 1
    assert fabs(hs_data['disp_volume'] - 1177) < 1


# def test_displacement_equilibrium():
#
#     disp = 2000
#     CG = [-1, 0, -3]
#
#     cylinder_eq, CGc = hs.compute_equilibrium(searev, disp, CG, rho_water=1023, anim=True)
    
